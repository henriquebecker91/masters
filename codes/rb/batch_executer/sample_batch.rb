#!/bin/ruby

require 'require_relative'
require_relative './sample_extractors.rb'
require_relative './batch.rb'

comms_info = [{
  # String with command to be executed. Must have 'pattern' as substring.
  command: 'sleep 1 && echo X X',
  # Substring present in 'command'. Often replaced by the instance filename.
  pattern: 'X',
  # Extractor object. Receives the output of the command and return
  # the most important fields.
  extractor: SampleExtractor.new,
  # String used to identify the command. Will be used to prefix the return of
  # extractor.names.
  prefix: 'doubled',
}, {
  command: 'sleep 3 && echo "banana X"',
  pattern: 'X',
  extractor: SampleExtractor.new,
  prefix: 'banana',
}, {
  command: 'sleep 100 && echo "never gonna happen X"',
  pattern: 'X',
  extractor: SampleExtractor.new,
  prefix: 'timeout',
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
  # The name of the file where will be written the CSV data.
  csvfname: 'sample.csv',
}

files = ['apple', 'orange'] # Applejack would be proud

Batch::experiment(comms_info, execution_info, conf, files)

