#!/bin/ruby

require 'batch_experiment'
require 'batch_experiment/sample_extractors'

comms_info = [{
  command: 'cat X',
  pattern: 'X',
  extractor: BatchExperiment::FirstLineExtractor,
  prefix: 'cat',
}, {
  command: 'echo y',
  pattern: 'y',
  extractor: BatchExperiment::FirstLineExtractor,
  prefix: 'echo',
}, {
  command: 'wc FILE',
  pattern: 'FILE',
  extractor: BatchExperiment::WcExtractor,
  prefix: 'wc',
}]

execution_info = {
  # IDs of the CPU cores that can be used for executing tests.
  cpus_available: [1, 2, 3],
  # Maximum number of seconds that a command can run. After this a kill command
  # (TERM signal) will be issued.
  timeout: 5,
}

conf = {
  # The name of the file where will be written the CSV data.
  csvfname: 'example.csv',
  # The columns will be ordered by command. All the columns of the first
  # command before the one from the second and so on.
  ic_columns: false,
}

files = ['bible.txt', 'taoteching.txt']

BatchExperiment::experiment(comms_info, execution_info, conf, files)

