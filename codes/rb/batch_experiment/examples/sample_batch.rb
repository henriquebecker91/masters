#!/bin/ruby

require 'batch_experiment'

commands = [
  'sleep 2 && echo orange',
  'sleep 4 && echo banana',
  'sleep 100 && echo "never gonna happen"',
]

conf = {
  # IDs of the CPU cores that can be used for executing tests.
  cpus_available: [0, 1],
  # Maximum number of seconds that a command can run. After this a kill command
  # (TERM signal) will be issued.
  timeout: 5,
  # Maximum number of seconds that a command can run after a kill command was
  # issued. After this a kill -9 command (KILL signal) will be issued.
  post_timeout: 1,
}

BatchExperiment::batch(commands, conf)

