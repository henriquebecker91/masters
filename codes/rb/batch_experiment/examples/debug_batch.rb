#!/bin/ruby

require_relative '../lib/batch_experiment.rb'

commands = []
10000.times { | n | commands << "sleep 1 && echo #{n}" }

conf = {
  cpus_available: [1, 2, 3],
  timeout: 5,
  post_timeout: 1,
}

BatchExperiment::batch(commands, conf)

