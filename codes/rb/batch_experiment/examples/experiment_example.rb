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
}]

batch_conf = {
  # IDs of the CPU cores that can be used for executing tests.
  cpus_available: [1, 2, 3],
  # Maximum number of seconds that a command can run. After this a kill command
  # (TERM signal) will be issued.
  timeout: 5,
  # Object that gives the filename for storing the output of each run.
  converter: BatchExperiment::Comm2FnameConverter.new,
}

experiment_conf = {
  # The name of the file where will be written the CSV data.
  csvfname: 'example.csv',
  # Number of times the same command will be executed over the same file.
  qt_runs: 5,
  # Order of the commands execution
  comms_order: :random, #:by_comm, #:by_file,
  # Random seed (only used if comms_order is :random)
  rng: Random.new(0),
}

files = ['bible.txt', 'taoteching.txt']

BatchExperiment::experiment(comms_info, batch_conf, experiment_conf, files)

