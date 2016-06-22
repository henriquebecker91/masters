#!/usr/bin/ruby

require 'batch_experiment'
require 'batch_experiment/sample_extractors'

# I run the three lines below in the console to disable hyperthreading cores on
# my computer before examining the cores with the top command.
# for i in 4 5 6 7; do
#   sudo sh -c "echo 0 > /sys/devices/system/cpu/cpu$i/online";
# done

comms_info = [{
  command: 'pyasukpt -src INST_FILE',
  pattern: 'INST_FILE',
  extractor: BatchExperiment::PyaExtractor,
  prefix: 'PYAsUKP',
}, {
  command: 'run_ukp5.out INST_FILE',
  pattern: 'INST_FILE',
  extractor: BatchExperiment::UKP5Extractor,
  prefix: 'UKP5',
}]

execution_info = {
  cpus_available: [1, 2, 3],
  timeout: 10,
  post_timeout: 2,
  cwd: '/home/henrique/AreaDeTrabalho/masters/data/ukp/',
  output_dir: '/home/henrique/AreaDeTrabalho/tmp/',
}

conf = {
  csvfname: 'pya_site8.csv',
  comms_order: :random,
  qt_runs: 10,
}

files = ['corepb.ukp', 'exnsd18.ukp', 'exnsd26.ukp', 'exnsdbis18.ukp', 'exnsd16.ukp', 'exnsd20.ukp', 'exnsdbis10.ukp', 'exnsds12.ukp']
# If you don't execute the script from the ukp files folder you need to put the
# folder relative or absolute path here (with trailing slash).

BatchExperiment::experiment(comms_info, execution_info, conf, files)

