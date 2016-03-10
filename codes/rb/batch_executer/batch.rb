#!/usr/bin/ruby

require 'childprocess'
require 'pathname'

# Henrique Becker -- 01/09/2015
# Usage: batch.rb <file with one command by line> <CPUs to use, delimited by conmas, without spaces, ex:0,1,3 > <timeout to every command> <timout after timeout for kill -9> 
# Will create a file.out and a file.err with the stdout and stderr of every
# command (based in the command text), not designed to support equal commands
# (will subscribe the file with the output of the last equal command)

# I use this in the bash first to disable hyperthreading cores on my computer
# for i in 4 5 6 7; do sudo sh -c "echo 0 > /sys/devices/system/cpu/cpu$i/online"; done

def sanitize_filename(filename)
  name = filename.strip 
  name.gsub!(/[^[:alnum:]]/, '_')
  name.gsub!(/_+/, '_')
  name
end

def batch(commands, cpus_avaliable, timeout, post_timeout)
  cpus_in_standby = cpus_avaliable.clone
  commands_running = []
  cpu = nil

  commands.each do | command |
    puts "waiting to execute command: #{command}"
    STDOUT.flush
    while not cpu = cpus_in_standby.pop do
      commands_running.delete_if do | job |
        cpus_in_standby.push(job[:cpu]) if job[:proc].exited?
        job[:proc].exited?
      end
      sleep 1
    end
    
    cproc = ChildProcess.build('timeout', '-k',  post_timeout, timeout, 'taskset', '-c', cpu.to_s, 'sh', '-c', command)
    base_name = sanitize_filename(command)
    out = File.open(base_name << '.out', 'w')
    err = File.open(base_name << '.err', 'w')
    cproc.io.stdout = out
    cproc.io.stderr = err
    cproc.start

    commands_running << { proc: cproc, cpu: cpu }

    puts "command assigned to cpu#{cpu}"
    STDOUT.flush
  end
end

if ARGV.size != 4 then
  puts "usage: batch <name of the file with a bash command on each line> <conma-separated list of cpus to use> <timeout> <time after timeout to wait before killing the program with -9>"
else
  f = File.open(ARGV[0], 'r')
  $all_lines = f.read.lines
  f.close

  batch($all_lines, ARGV[1].split(',').map! { | x | x.to_i }, ARGV[2], ARGV[3])
end

