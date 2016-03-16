require 'childprocess'
require 'pathname'

module Batch
  # I use this in the bash first to disable hyperthreading cores on my computer
  # for i in 4 5 6 7; do sudo sh -c "echo 0 > /sys/devices/system/cpu/cpu$i/online"; done

  def sanitize_filename(filename)
    name = filename.strip
    name.gsub!(/[^[:alnum:]]/, '_')
    name.gsub!(/_+/, '_')
    name
  end

  def clean_finished(cpus_in_standby, commands_running)
    commands_running.delete_if do | job |
      if job[:proc].exited?
        cpus_in_standby.push(job[:cpu])
        File.delete(job[:lockfname])
      end
      job[:proc].exited? # bool returned to delete_if
    end
  end

  # Takes a list of commands, execute them only on the designed core/cpus, and kill them if the timeout expires, never lets a core/cpu rest for more than conf[:busy_loop_sleep] seconds between a command and another. The return is a list of partial filenames. Appending '.out' to one of the partial filenames will give the filename were the command stdout was redirected. The analogue is valid for '.err' and stderr. The first partial filename corresponds to the first command in commands, and so on. Right before a command begans to run, a "#{partial_filename}.#{conf[:unfinished_ext]}" file is created. After the command ends its execution, this file is removed. If the command ends its execution by means of a timeout the file is also removed. The file only remains if the batch procedure is interrupted (not a specific command).
  #
  # @note This procedure was not designed to support equal commands (the last equal command executed will subscribe the '.out', '.err' and '.unfinished' files used by any previous equal command). 
  # @note This procedure makes use of the following linux commands: time (not the bash internal one, but the package one, i.e. https://www.archlinux.org/packages/extra/x86_64/time/); timeout (from coreutils); taskset (from util-linux, https://www.archlinux.org/packages/core/x86_64/util-linux/); sh (the shell).
  # @note The command is executed inside a call to "sh -c command", so it has to be a valid sh command.
  # @note 
  def batch(commands, conf)
    # Throw exceptions if required conigurations aren't provided.
    #conf[:timeout] 

    # Initialize optional configurations with default values if they aren't
    # provided. Don't change the conf argument, only our version of conf.
    conf = conf.clone
    conf[:time_fmt]         ||= 'ext_time: %e\\next_mem: %M\\n'
    conf[:unfinished_ext]   ||= '.unfinished'
    conf[:busy_loop_sleep]  ||= 0.1
    conf[:post_timeout]     ||= 5

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

      commfname = sanitize_filename(command)

      cproc = ChildProcess.build(
        'timeout', '-k', "#{conf[:post_timeout]}s", "#{conf[:timeout]}s",
        'taskset', '-c', cpu.to_s,
        'time', '-f', conf[:time_fmt], '--append', '-o', "#{commfname}.out",
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

    unless commands_running.empty? do
      sleep conf[:busy_loop_sleep]
      update_finished(cpus_available, commands_running)
    end

    commfnames
  end
end

