require 'childprocess'
require 'pathname'

# The main module, the two main utility methods offered are ::batch and
# ::experiment.
module BatchExperiment
  # The default callable object used by Comm2FnameConverter to convert
  # a command into a filename. Comm2FnameConverter don't create a sanitized
  # filename from the command string (it uses its first argument to do this,
  # whose default is FnameSanitizer).
  # Note that this is a pure function, so if the same command appears more than
  # one time, it will get the same name, it's Comm2FnameConverter that gives
  # multiple instances of the same command different names (by suffixing with
  # numbers).
  module FnameSanitizer
    def self.call(command)
      fname = command.strip
      fname.gsub!(/[^[:alnum:]]/, '_')
      fname.gsub!(/_+/, '_')
      fname.gsub!(/^_/, '')
      fname.gsub!(/_$/, '')
      fname
    end
  end

  # Converts a command to a filename using a given sanitizer, gives different
  # names to different calls with the same arguments. Example: if a call with
  # "sleep 1" yields "sleep_1", the second call with the same argument yields
  # "sleep_1.2", and so on. Note that this is done by remembering previous
  # calls, the object don't inspect the filesystem to check if that name
  # was or wasn't used.
  class Comm2FnameConverter
    # Creates a new Comm2FnameConverter, with no memory of any previous calls.
    #
    # @param sanitizer [#call] Callable object used to create a filename from
    #   the arguments passed to Comm2FnameConverter.call. This class expects
    #   that sanitizer has no internal state, so when an instance of this class
    #   is cloned, there's no problem with sharing the sanitizer between the
    #   clones. Default: BatchExperiment::FnameSanitizer.
    def initialize(sanitizer = FnameSanitizer)
      @num_times_seen = {}
      @sanitizer = sanitizer
    end

    # Takes a command, creates a fname for it, if this fname was already seen
    # before, returns the fname + ".N", where N is the number of times fname
    # was already seen.
    #
    # @param comm [String] A system command.
    # @return [String] The sanitized filename created from that command.
    # @note Note that different arguments can be reduced to the same
    #   sanitized filename and, if this happens, they will NOT overwrite
    #   each other. Example: 'echo "abc"' -> 'echo_abc'; 'echo abc' ->
    #   'echo_abc.2'.
    def call(comm)
      fname = @sanitizer.call(comm)
      if @num_times_seen.include? fname
        @num_times_seen[fname] += 1
        fname << ".#{@num_times_seen[fname]}"
      else
        @num_times_seen[fname] = 1
      end

      fname.clone
    end

    def initialize_clone(old)
      @num_times_seen = old.num_times_seen.clone
    end

    # To allow the initialize_clone implementation.
    protected
    attr_reader :num_times_seen
  end

  # INTERNAL USE ONLY.
  # Remove any finished commands from comms_running, insert the cpus
  # freed by the commands termination to the free_cpus, insert the
  # terminated commands on comms_executed.
  def self.update_finished(free_cpus, comms_running, comms_executed)
    comms_running.delete_if do | job |
      # Don't call '#exited?' twice, store value at variable. If you call
      # it twice it's possible to remove it from the list of running commands
      # without freeing a cpu, what will end locking all cpus forever.
      exited = job[:proc].exited?
      if exited
        free_cpus.push(job[:cpu])
        File.delete(job[:lockfname])
        comms_executed << job[:command]
      end
      exited # bool returned to delete_if
    end
  end

  # Takes a list of commands, execute them only on the designed core/cpus, and
  # kill them if the timeout expires, never lets a core/cpu rest for more than
  # a predetermined amount of seconds between a command and another. Partial
  # filenames are derived from the commands. Appending '.out' to one of the
  # partial filenames will give the filename were the command stdout was
  # redirected. The analogue is valid for '.err' and stderr. Right before a
  # command begans to run, a 'partial_filename.unfinished' file is created.
  # After the command ends its execution this file is removed. If the command
  # ends its execution by means of a timeout the file is also removed. The file
  # only remains if the batch procedure is interrupted (script was killed,
  # or system crashed). This '.unfinished' file will contain the process pid,
  # if the corresponding process started with success.
  #
  # @param commands [Array<String>] The shell commands.
  # @param conf [Hash] The configurations, as follows:
  #   -- cpus_available [Array<Fixnum>] CPU cores that can be used to run the
  #   commands. Required parameter. The cpu numbers begin at 0, despite what
  #   htop tells you.
  #   -- timeout [Number] Number of seconds before killing a command. Required
  #   parameter. Is the same for all the commands.
  #   -- time_fmt [String] A string in the time (external command) format. See
  #   http://linux.die.net/man/1/time. Default: 'ext_time: %e\next_mem: %M\n'.
  #   -- busy_loop_sleep [Number] How many seconds to wait before checking if
  #   a command ended execution. This time will be very close to the max time a
  #   cpu will remain vacant between two commands. Default: 0.1 (1/10 second).
  #   -- post_timeout [Number] A command isn't guaranteed to end after
  #   receiving a TERM signal. If the command hasn't stopped, waits
  #   post_timeout seconds before sending a KILL signal (give it a chance to
  #   end gracefully). Default: 5.
  #   -- converter [#call] The call method of this object should take a String
  #   and convert it (possibly losing information), to a valid filename. Used
  #   over the commands to define the output files of commands.
  #   Default: BatchExperiment::Comm2FnameConverter.new.
  #   -- skip_done_comms [FalseClass,TrueClass] Skip any command for what a
  #   corresponding '.out' file exists, except if both a '.out' and a
  #   '.unfinished' file exist, in the last case the command is always
  #   executed. If false, execute all commands and overwrite all ".out".
  #   Default: true.
  #   -- unfinished_ext [String] Extension to be used in place of
  #   '.unfinished'. Default: '.unfinished'.
  #   -- out_ext [String] Extension to be used in place of '.out'.
  #   Default: '.out'.
  #   -- err_ext [String] Extension to be used in place of '.err'.
  #   Default: '.err'.
  #   -- cwd [String] Command Working Directory. The path from where the
  #   commands will be executed. Default: './' (i.e. the same directory from
  #   where the ruby script was run).
  #   -- output_dir [String] The folder used to save the
  #   '.out'/'.err'/'.unfinished' files. Default: './' (i.e. the same directory
  #   from where the ruby script was run).
  #
  # @return [String] Which commands were executed. Can be different from
  #   the 'commands' argument if commands are skipped (see :skip_done_comms).
  #
  # @note If the same command is executed over the same file more than one
  #   time, then any run besides the first will have a numeric suffix.
  #   Example: "sleep 1" -> "sleep_1", "sleep 1" -> "sleep_1.2".
  #   For more info see the parameter conf\[:fname_sanitizer\], and its
  #   default value BatchExperiment::Comm2FnameConverter.new.
  # @note This procedure makes use of the following linux commands: time (not
  #   the bash internal one, but the package one, i.e.
  #   https://www.archlinux.org/packages/extra/x86_64/time/); timeout (from
  #   coreutils); taskset (from util-linux,
  #   https://www.archlinux.org/packages/core/x86_64/util-linux/); sh (the
  #   shell).
  # @note The command is executed inside a call to "sh -c command", so it has
  #   to be a valid sh command.
  # @note The output of the command "time -f conf\[:time_fmt\]" will be
  #   appended to the '.out' file of every command. If you set
  #   conf\[:time_fmt\] to a empty string only a newline will be appended.
  def self.batch(commands, conf)
    # Throw exceptions if required configurations aren't provided.
    if !conf[:cpus_available] then
      fail ArgumentError, 'conf[:cpus_available] not set'
    end
    fail ArgumentError, 'conf[:timeout] not set' unless conf[:timeout]

    # Initialize optional configurations with default values if they weren't
    # provided. Don't change the conf argument, only our version of conf.
    conf = conf.clone
    conf[:time_fmt]         ||= 'ext_time: %e\\next_mem: %M\\n'
    conf[:unfinished_ext]   ||= '.unfinished'
    conf[:out_ext]          ||= '.out'
    conf[:err_ext]          ||= '.err'
    conf[:busy_loop_sleep]  ||= 0.1
    conf[:post_timeout]     ||= 5
    conf[:converter]        ||= BatchExperiment::Comm2FnameConverter.new
    conf[:skip_done_comms]    = true if conf[:skip_done_comms].nil?
    conf[:cwd]              ||= './'
    conf[:output_dir]       ||= './'

    # Initialize main variables
    free_cpus = conf[:cpus_available].clone
    comms_running = []
    cpu = nil
    comms_executed = []

    commands.each do | command |
      commfname = conf[:converter].call(command)
      out_fname = conf[:output_dir] + commfname + conf[:out_ext]
      err_fname = conf[:output_dir] + commfname + conf[:err_ext]
      lockfname = conf[:output_dir] + commfname + conf[:unfinished_ext]

      if conf[:skip_done_comms] && File.exists?(out_fname)
        if File.exists?(lockfname)
          puts "Found file #{out_fname}, but a #{lockfname} also exists:"
          puts "Will execute command '#{command}' anyway."
        else
          puts "Found file #{commfname}, skipping command: #{command}"
          STDOUT.flush
          next
        end
      end

      puts "Waiting to execute command: #{command}"
      STDOUT.flush

      while free_cpus.empty? do
        sleep conf[:busy_loop_sleep]
        update_finished(free_cpus, comms_running, comms_executed)
      end

      cpu = free_cpus.pop

      cproc = ChildProcess.build(
        'taskset', '-c', cpu.to_s,
        'time', '-f', conf[:time_fmt], '--append', '-o', out_fname,
        'timeout', '--preserve-status', '-k', "#{conf[:post_timeout]}s",
          "#{conf[:timeout]}s",
        'sh', '-c', command
      )

      cproc.cwd = conf[:cwd]

      File.open(lockfname, 'w') {} # empty on purpose
      out = File.open(out_fname, 'w')
      err = File.open(err_fname, 'w')
      cproc.io.stdout = out
      cproc.io.stderr = err

      cproc.start

      comms_running << {
        proc: cproc,
        cpu: cpu,
        lockfname: lockfname,
        command: command
      }

      # The lock file now stores the process pid for debug reasons.
      File.open(lockfname, 'w') { | f | f.write cproc.pid }

      puts "command assigned to cpu#{cpu}"
      STDOUT.flush
    end

    until comms_running.empty? do
      sleep conf[:busy_loop_sleep]
      update_finished(free_cpus, comms_running, comms_executed)
    end

    comms_executed
  end

  # INTERNAL USE ONLY. gencommff: GENerate COMMands For Files.
  # Creates a hash with the generated commands as keys, and store (as the
  # respective value) the comm_info hash and the file (using a { comm_info: X,
  # filename: Y } structure).
  #
  # @param comm_info [Hash] A hash structure following the same format
  #   that the elements of the comms_info array parameter of #experiment.
  # @param files [Enumerable<String>] A list of strings that will replace
  #   comm_info[:pattern] at a copy of comm_info[:command].
  # @return [Hash<String, Hash>] A hash on the following format
  #   { expanded_command => { comm_info: comm_info, filename: f }, ...}
  def self.gencommff(comm_info, files)
    ret = {}
    comm = comm_info[:command]
    patt = comm_info[:pattern]
    files.each do | f |
      ret[comm.gsub(patt, f)] = { comm_info: comm_info, filename: f }
    end
    ret
  end

  # INTERNAL USE ONLY. Intercalate a variable number of variable sized arrays
  # in one array.
  #
  # @param [Array<Array<Object>>] xss An array of arrays.
  # @return [Array<Object>] An array of the same size as the sum of the size
  #   of all inner arrays. The values are the same (not copies) as the values
  #   of the array. Example: intercalate([[1, 4, 6, 7], [], [2, 5], [3]])
  #   returns [1, 2, 3, 4, 5, 6, 7].
  def self.intercalate(xss)
    ret = []
    xss = xss.map { | xs | xs.reverse }
    until xss.empty? do
      xss.delete_if do | xs |
        unless xs.empty?
          ret << xs.pop
        end
        xs.empty?
      end
    end
    ret
  end

  class ColumnSpecError < ArgumentError; end

  # INTERNAL USE ONLY. Check if the headers can be combined, if they can
  # return a shallow copy of the biggest header, otherwise throw an exception.
  #
  # @param headers [Array<Array<Comparable>>] An array of arrays of strings
  #   (or any object that implements '!=').
  # @return A shallow copy of the biggest inner array in headers. Only returns
  #   if for each position on the biggest inner array has the same value as
  #   that position on all the other arrays with at least that size.
  def self.merge_headers(headers)
    mer_size = headers.map { | h | h.size }.max
    merged_h = Array.new(mer_size)
    mer_size.times do | i |
      headers.each do | h |
        next if h.size < i
        if merged_h[i].nil?
          merged_h[i] = h[i]
        elsif merged_h[i] != h[i]
          raise ColumnSpecError, "Error: When using BatchExperiment::experiment"
            + " all the extractors have to agree on the columns they share." 
            + " In the specific case: the column nº #{i} was labeled as"
            + " '#{merged_h[i]}' on one extractor, and '#{h[i]}' on another,"
            + " this can be only a difference on notation ('time' vs 'Time'),"
            + " or can mean that in the same column two different kinds of data"
            + " are being presented. The program will be aborted. Check that."
        end
      end
    end
    merged_h
  end

  # Takes N shell commands and M files/parameters, execute each command of the
  # N commands over the M files, save the output of each command/file
  # combination, use objects provided with the command to extract relevant
  # information from the output file, and group those information in a CSV
  # file. Easier to understand seeing the sample_batch.rb example in action.
  #
  # @param comms_info [Array<Hash>] An array of hashs, each with the config
  #   needed to know how to deal with the command. Four required fields
  #   (all keys are symbols):
  #   command [String] A string with a sh shell command.
  #   pattern [String] A substring of command, will be replaced by the strings
  #   in the paramenter 'files'.
  #   extractor [#extract,#names] Object implementing the Extractor interface.
  #   prefix [String] A string that will be used on the 'algorithm' column
  #   to identify the used command.
  # @param batch_conf [Hash] Configuration used to call batch. See the
  #   explanation for parameter 'conf' on the documentation of the batch
  #   method. There are required fields for this hash parameter. Also, note
  #   that the batch_conf\[:converter\] should allow cloning without sharing
  #   mutable state. A converter clone is used by #experiment internally, it
  #   has to obtain the same results as the original copy (that is passed to
  #   BatchExperiment::batch).
  # @param conf [Hash] Lots of parameters. Here's a list:
  #   -- csvfname [String] The filename/filepath for the file that will contain
  #   the CSV data. Required field.
  #   separator [String] The separator used at the CSV file. Default: ';'.
  #   -- qt_runs [NilClass,Integer] If nil or one then each command is
  #   executed once. If is a number bigger than one, the command is executed
  #   that number of times. The batch_conf[:converter] will define the name
  #   that will be given to each run. Every file will appear qt_runs times on
  #   the filename column and, for the same file, the values on the run_number
  #   column will be the integer numbers between 1 and qt_runs (both
  #   inclusive). Default: nil.
  #   -- comms_order [:by_comm,:by_file,:random] The order the
  #   commands will be executed. Case by_comm: will execute the first command
  #   over all the files (using the files order), then will execute the
  #   second command over all files, and so on. Case by_file: will execute
  #   all the commands (using the comms_info order) over the first file,
  #   then will execute all the comands over the second file, and so on.
  #   Case random: will expand all the command/file combinations (replicating
  #   the same command qt_run times) and then will apply shuffle to this array,
  #   using the object passed to the rng parameter. This last option is the
  #   most adequate for statistical testing.
  #   -- rng [Nil,#rand] An object that implements the #rand method (behaves
  #   like an instance of the core Random class). If comms_order is random and
  #   rng is nil, will issue a warning remembering the default that was used.
  #   Default: Random.new(42).
  #   skip_commands [TrueClass, FalseClass] If true, will not execute the
  #   commands and assume that the outputs are already saved (on ".out" files).
  #   Will only execute the extractors over the already saved outputs, and
  #   create the CSV file from them. Default: false.
  #
  # @param files [Array<Strings>] The strings that will replace the :pattern
  #   on :command, for every element in comms_info. Can be a filename, or
  #   can be anything else (a numeric parameter, sh code, etc..), but we
  #   refer to them as files for simplicity and uniformity.
  #
  # @return [NilClass,Array<String>] The return of the internal #batch
  #   call. Returns nil if conf[:skip_commands] was set to true.
  #
  # @see BatchExperiment::batch
  # @note This command call ::batch internally.
  def self.experiment(comms_info, batch_conf, conf, files)
    # Throw exceptions if required configurations aren't provided.
    fail 'conf[:csvfname] is not defined' unless conf[:csvfname]

    # Initialize optional configurations with default values if they weren't
    # provided. Don't change the conf argument, only our version of conf.
    conf = conf.clone
    conf[:separator]    ||= ';'
    conf[:qt_runs]      ||= 1
    conf[:comms_order]  ||= :by_comm
    conf[:rng]          ||= Random.new(42)
    #conf[:skip_commands] defaults to false/nil

    # Get some of the batch config that we use inside here too.
    out_ext         = batch_conf[:out_ext]        || '.out'
    unfinished_ext  = batch_conf[:unfinished_ext] || '.unfinished'
    output_dir      = batch_conf[:output_dir]     || './'
    converter = batch_conf[:converter].clone unless batch_conf[:converter].nil?
    converter ||= BatchExperiment::Comm2FnameConverter.new

    # Expand all commands, combining command templates and files.
    comms_sets = []
    comms_info.each do | comm_info |
      comms_sets << gencommff(comm_info, files)
    end

    expanded_comms = comms_sets.map { | h | h.keys }
    # If each command should be run more than once...
    if conf[:qt_runs] > 1
      # ... we replace each single command by an array of qt_runs copies,
      # and then flatten the parent array.
      expanded_comms.map! do | a |
        a.map! { | c | Array.new(conf[:qt_runs], c) }.flatten!
      end
    end

    # At this moment the expanded_comms is an array of arrays, each internal
    # array has all the expanded commands of the one single command template
    # over all the files.
    # After the code block below, the expanded_comms will be an one-level array
    # of the expanded commands, in the order they will be executed.
    expanded_comms = case conf[:comms_order]
    when :by_comm # all runs of the first command template first
      expanded_comms.flatten!
    when :by_file # all runs over the first file first
      intercalate(expanded_comms)
    when :random  # a random order
      expanded_comms.flatten!.shuffle!(random: conf[:rng])
    end

    # Execute the commands (or not).
    ret = batch(expanded_comms, batch_conf) unless conf[:skip_commands]

    # Build header (first csv line, column names).
    header = ['algorithm', 'filename', 'run_number']
    header << merge_headers(comms_info.map { | c | c[:extractor].names })
    header = header.join(conf[:separator])

    # We need to merge the union of all comms_sets to query it.
    comm2origin = {}
    comms_sets.each do | h |
      comm2origin.merge!(h) do | k, v, v2 |
        puts "WARNING: The command expansion '#{k}' was generated more than once. The first time was by the template '#{v[:comm]}' and the file '#{v[:file]}', and this time by template '#{v2[:comm]}' and the file '#{v2[:file]}'. Will report on CSV as this command was generated by the template '#{v[:comm]}' and the file '#{v[:file]}'."
        v
      end
    end

    # Build body (inspect all output files and make csv lines).
    #
    # Body format: algorithm;filename;run_number;first extracted column; ...
    #
    # This means that the extractors have to agree on what is each column, two
    # different extractors have to extract the same kind of data at each column
    # (the first field returned by all extractors has to be, for example, cpu
    # time, the same applies for the remaining fields).
    # If one extractor extract more fields than the others this is not a
    # problem, if the second biggest extractor (in number of fields extract)
    # will extract, for example, 4 fields, and the biggest extract 6 fields,
    # the first 4 fields extracted by the biggest extractor have to be the same
    # as the ones on the second-biggest extractor. This way, all the lines will
    # have the kind of data on the first four columns (not counting the
    # algorithm, filename and run_number ones), and only lines provenient from
    # the biggest extractor will have data on the fifth and sixth columns.
    body = [header]
    times_found = {}
    expanded_comms.each do | exp_comm |
      run_info   = comm2origin[exp_comm]
      algorithm  = run_info[:comm_info][:prefix]
      filename   = run_info[:filename]

      times_found[exp_comm] ||= 0
      times_found[exp_comm]  += 1
      run_number = times_found[exp_comm]

      curr_line = [algorithm, filename, run_number]

      partial_fname = converter.call(exp_comm)
      out_fname = output_dir + partial_fname + out_ext
      lockfname = output_dir + partial_fname + unfinished_ext
      extractor = run_info[:comm_info][:extractor]

      if File.exists?(out_fname)
        if File.exists?(lockfname)
          puts "Ignored file '#{out_fname}' because there was a"
             + "  '#{lockfname}' file too."
        else
          f_content = File.open(out_fname, 'r') { | f | f.read }
          curr_line << extractor.extract(f_content)
        end
      end

      body << curr_line.join(conf[:separator])
    end
    body = body.join(conf[:separator] + "\n")

    # Write CSV data into a CSV file.
    File.open(conf[:csvfname], 'w') { | f | f.write(body) }

    return ret
  end
end

