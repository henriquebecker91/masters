require 'childprocess'
require 'pathname'

module Batch
  # The default callable class used by batch to convert a command into a
  # filename.
  class FilenameSanitizer
    def call(command)
      fname = command.strip
      fname.gsub!(/[^[:alnum:]]/, '_')
      fname.gsub!(/_+/, '_')
      fname
    end
  end

  # Internal use only. DO NOT DEPEND.
  # Remove the finished commands from commands_running, and insert the cpus
  # freed by the commands termination to the cpus_in_standby.
  def self.update_finished(cpus_in_standby, commands_running)
    commands_running.delete_if do | job |
      if job[:proc].exited?
        cpus_in_standby.push(job[:cpu])
        File.delete(job[:lockfname])
      end
      job[:proc].exited? # bool returned to delete_if
    end
  end

  # Takes a list of commands, execute them only on the designed core/cpus, and
  # kill them if the timeout expires, never lets a core/cpu rest for more than
  # conf[:busy_loop_sleep] seconds between a command and another. The return is
  # a list of partial filenames. Appending '.out' to one of the partial
  # filenames will give the filename were the command stdout was redirected.
  # The analogue is valid for '.err' and stderr. The first partial filename
  # corresponds to the first command in commands, and so on. Right before a
  # command begans to run, a "#{partial_filename}.#{conf[:unfinished_ext]}"
  # file is created. After the command ends its execution, this file is
  # removed. If the command ends its execution by means of a timeout the file
  # is also removed. The file only remains if the batch procedure is
  # interrupted (not a specific command).
  #
  # @param commands [Array<String>] The shell commands.
  # @param conf [Hash] The configurations, as follows:
  #   :cpus_available [Array<Fixnum>] Cpu cores that can be used to run the
  #   commands. Required parameter. The cpu numbers begin at 0, despite what
  #   htop tells you;
  #   :timeout [Number] Number of seconds before killing a command. Required
  #   parameter. Is the same for all the commands;
  #   :time_fmt [String] A string in the time (external command) format. See
  #   http://linux.die.net/man/1/time. Default: 'ext_time: %e\next_mem: %M\n'.
  #   :unfinished_ext [String] Extension to be used in place of '.unfinished'.
  #   Default: '.unfinished';
  #   :busy_loop_sleep [Number] How many seconds to wait before checking if a
  #   command ended execution. This is max time a cpu will be vacant between
  #   two commands. Default: 0.1;
  #   :post_timeout [Number] A command isn't guaranteed to end after receiving
  #   a TERM signal. If the command hasn't stopped, waits post_timeout seconds
  #   before sending a KILL signal (give it a chance to end gracefully).
  #   Default: 5;
  #   :fname_sanitizer [Callable Object] The call method of this object
  #   should take a String and convert it (possibly losing information), to a
  #   valid filename. Used over the commands to define the output files of
  #   commands.
  #   Default: Batch::FilenameSanitizer.new;
  #
  # @note This procedure was not designed to support equal commands (the last
  #   equal command executed will subscribe the '.out', '.err' and '.unfinished'
  #   files used by any previous equal command). But the parameter
  #   conf[:fname_sanitizer] can be used to circumvent the restriction over
  #   equal commands (if the object has state it can return a different
  #   filename for every time it's called with the same argument).
  # @note This procedure makes use of the following linux commands: time (not
  #   the bash internal one, but the package one, i.e.
  #   https://www.archlinux.org/packages/extra/x86_64/time/); timeout (from
  #   coreutils); taskset (from util-linux,
  #   https://www.archlinux.org/packages/core/x86_64/util-linux/); sh (the
  #   shell).
  # @note The command is executed inside a call to "sh -c command", so it has
  #   to be a valid sh command.
  # @note The output of the command "time -f #{conf[:time_fmt]}" will be
  #   appended to the '.out' file of every command. If you set conf[:time_fmt]
  #   to a empty string only a newline will be appended.
  #
  def self.batch(commands, conf)
    # Throw exceptions if required configurations aren't provided.
    fail "conf[:cpus_available] not set" unless conf[:cpus_available]
    fail "conf[:timeout] not set" unless conf[:timeout] 

    # Initialize optional configurations with default values if they weren't
    # provided. Don't change the conf argument, only our version of conf.
    conf = conf.clone
    conf[:time_fmt]           ||= 'ext_time: %e\\next_mem: %M\\n'
    conf[:unfinished_ext]     ||= '.unfinished'
    conf[:busy_loop_sleep]    ||= 0.1
    conf[:post_timeout]       ||= 5
    conf[:fname_sanitizer] ||= Batch::FilenameSanitizer.new

    # Initialize main variables
    cpus_in_standby = conf[:cpus_available].clone
    commands_running = []
    cpu = nil
    commfnames = []

    # For each command, put it to execute in a core (if there's a core
    # available), or wait until a core becomes available. The order in what the
    # commands are started is deterministic (is the order of commands). The
    # command ending order, and what core they execute, is not deterministic.
    commands.each do | command |
      puts "waiting to execute command: #{command}"
      STDOUT.flush

      while cpus_in_standby.empty? do
        sleep conf[:busy_loop_sleep]
        update_finished(cpus_in_standby, commands_running)
      end

      cpu = cpus_in_standby.pop 

      commfname = conf[:fname_sanitizer].call(command)

      cproc = ChildProcess.build(
        'taskset', '-c', cpu.to_s,
        'time', '-f', conf[:time_fmt], '--append', '-o', "#{commfname}.out",
        'timeout', '--preserve-status', '-k', "#{conf[:post_timeout]}s",
          "#{conf[:timeout]}s",
        'sh', '-c', command)

      lockfname = commfname + conf[:unfinished_ext]

      File.open(lockfname, 'w') {} # empty on purpose
      out = File.open(commfname + '.out', 'w')
      err = File.open(commfname + '.err', 'w')

      cproc.io.stdout = out
      cproc.io.stderr = err

      cproc.start

      commands_running << {
        proc: cproc,
        cpu: cpu,
        lockfname: lockfname
      }

      puts "command assigned to cpu#{cpu}"
      STDOUT.flush

      commfnames << commfname
    end

    until commands_running.empty? do
      sleep conf[:busy_loop_sleep]
      update_finished(cpus_in_standby, commands_running)
    end

    commfnames
  end

  # gencommff: GENerate COMMands For Files
  #
  # @param comm [String] A string with 'patt' as a substring.
  # @param patt [String] A string contained in 'comm'.
  # @param files [Enumerable<String>] A list of strings to substitute patt at
  #   comm.
  # @return [Array<String>] Example: gencommff('echo STR', 'STR', ['a', 'b',
  #   'c']) returns ['echo a', 'echo b', 'echo c'].
  def self.gencommff(comm, patt, files)
    ret = []
    files.each { | f | ret << comm.gsub(patt, f) }
    ret
  end

  # Intercalate a variable number of variable sized arrays in one array.
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

  # Takes N shell commands and M files/parameters, execute each command of the N commands over the M files, save the output of each command/file combination, use objects provided with the command to extract relevant information from the 
  def self.experiment(comms_info, batch_conf, conf, files)
    # We will make of it too, so it's useful for us to have it already defined
    conf[:fname_sanitizer] ||= Batch::FilenameSanitizer.new

    command_sets = []
    comms_info.each do | comm_info |
      command_sets << gencommff(comm_info[:command], comm_info[:pattern], files)
    end

    comm_list = if conf[:ic_comm]
      intercalate(command_sets)
    else
      command_sets.flatten
    end
    
    batch(comm_list, batch_conf)

    header = []
    comms_info.each do | comm_info |
      prefixed_names = comm_info[:extractor].names.map do | name |
        (comm_info[:prefix] + ' ') << name
      end
      header << prefixed_names
    end
    header = intercalate(header) if conf[:ic_columns]
    header = ['Filename'].concat(header).join(conf[:separator])

    body = [header]
    files.each_with_index do | inst_fname, j |
      line = []
      comms_info.each_with_index do | comm_info, i |
        output_fname = if conf[:ic_comm]
          comm_list[(j * comms_info.size) + i]
        else
          comm_list[(i * files.size) + j]
        end
        output_fname = conf[:fname_sanitizer].call(output_fname)
          
        f_content = File.open(output_fname + '.out') { | f | f.read }

        extractor = comm_info[:extractor]
        line << extractor.extract(f_content)
      end
      line = intercalate(line) if conf[:ic_columns]
      body << [inst_fname].concat(line).join(conf[:separator])
    end
    body = body.map! { | line | line << ';' }.join("\n")

    File.open(conf[:csvfname], 'w') { | f | f.write(body) }
  end
end

