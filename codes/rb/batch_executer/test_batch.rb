#!/usr/bin/ruby

require 'require_relative'
require_relative './batch.rb'
require_relative './sample_extractors.rb'

# I run the three lines below in the console to disable hyperthreading cores on
# my computer before examining the cores with the top command.
# for i in 4 5 6 7; do
#   sudo sh -c "echo 0 > /sys/devices/system/cpu/cpu$i/online";
# done

comms_info = [{
  # String with command to be executed. Must have patt as substring.
  command: 'pyasukpt -src INST_FILE',
  # Substring present in 'command'. Often replaced by the instance filename.
  pattern: 'INST_FILE',
  # Extractor subclass object. Receives the output of the command and return
  # the most important fields.
  extractor: PyaExtractor.new,
  # String used to identify the command. Will be used to prefix the return of
  # extractor.names.
  prefix: 'PYAsUKP',
}, {
  command: 'run_ukp5.out INST_FILE',
  pattern: 'INST_FILE',
  extractor: UKP5Extractor.new,
  prefix: 'UKP5',
}]

execution_info = {
  # IDs of the CPU cores that can be used for executing tests.
  cpus_available: [1, 2, 3],
  # Maximum number of seconds that a command can run. After this a kill command
  # will be issued.
  timeout: 10,
  # Maximum number of seconds that a command can run after a kill command was
  # issued. After this a kill -9 command will be issued.
  post_timeout: 5,
}

conf = {
  # skip_commands # If true, will not execute the commands and assume that the
  # outputs are already saved. Will only execute the gatherer over the saved
  # outputs. NOT IMPLEMENTED YET.
  # The string that will be used as separator at the CSV.
  separator: ';',
  # csvfname: String (valid as a filename).
  # The name of the file where will be written the CSV file with all the
  # info gathered.
  csvfname: 'pya_site8.log',
  unfinished_ext: '.unfinished',
}

files = ['corepb.ukp', 'exnsd18.ukp', 'exnsd26.ukp', 'exnsdbis18.ukp', 'exnsd16.ukp', 'exnsd20.ukp', 'exnsdbis10.ukp', 'exnsds12.ukp']
path = ''
files.map! { | f | path + f }

experiment(comms_info, execution_info, conf, files)

