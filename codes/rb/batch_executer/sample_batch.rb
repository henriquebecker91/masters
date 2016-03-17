#!/bin/ruby

require 'require_relative'
require_relative './sample_extractors.rb'
require_relative './batch.rb'

comms_info = [{
  # String with command to be executed. Must have 'pattern' as substring.
  command: 'sleep 1 && echo X X',
  # Substring present in 'command'. Often replaced by the instance filename.
  pattern: 'X',
  # Extractor subclass object. Receives the output of the command and return
  # the most important fields.
  extractor: SampleExtractor.new,
  # String used to identify the command. Will be used to prefix the return of
  # extractor.names.
  prefix: 'alpha',
}, {
  command: 'sleep 3 && echo "banana X"',
  pattern: 'X',
  extractor: SampleExtractor.new, 
  prefix: 'beta',
}, {
  command: 'sleep 100 && echo "never gonna happen X"',
  pattern: 'X',
  extractor: SampleExtractor.new, 
  prefix: 'omega',
}]

execution_info = {
  # IDs of the CPU cores that can be used for executing tests.
  cpus_available: [1, 2, 3],
  # Maximum number of seconds that a command can run. After this a kill command
  # (TERM signal) will be issued.
  timeout: 5,
  # Maximum number of seconds that a command can run after a kill command was
  # issued. After this a kill -9 command (KILL signal) will be issued.
  post_timeout: 1,
}

conf = {
  #skip_commands # If true, will not execute the commands and assume that the
  # outputs are already saved. Will only execute the gatherer over the saved
  # outputs. NOT IMPLEMENTED YET.

  # The string that will be used as separator at the CSV.
  separator: ';',
  # csvfname: String (valid as a filename).
  # The name of the file where will be written the CSV file with all the
  # info gathered.
  csvfname: 'sample.log',
  unfinished_ext: '.unfinished',
  # Intercalate the data returned by the extractors. In other words, the csv
  # line for some file will not present all fields of the first command, then
  # all fields of the second command, ..., but instead will present the first
  # field of all commands, the second field of all commands, and so on.
  ic_columns: true,
  # Intercalate the commands execution. Instead of executing the first command
  # over all files first, execute all the commands over the first file first.
  # This was made to avoid confounding (statistical concept). If something
  # disrupts the processing power for some period of time, the effect will
  # probably be distributed between commands. The risk some algorithm seems
  # better or worse than it really is will be reduced. For example: you are
  # making tests at an notebook, the notebook becomes unplugged for a short
  # time. The cores will probably enter in energy saving mode and affect the
  # observed performance. If this happens when all tested commands are the
  # same, then will seem that that an command had a worse performance. If this
  # happens when the commands are intercalated, then maybe some instances will
  # seem harder than others (what is less problematic).
  ic_comm: true,
}

files = ['apple', 'orange'] # Applejack would be proud

Batch::experiment(comms_info, execution_info, conf, files)

