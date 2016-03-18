#!/usr/bin/ruby

require_relative 'batch_experiment'
require_relative 'batch_experiment/sample_extractors'

# I run the three lines below in the console to disable hyperthreading cores on
# my computer before examining the cores with the top command.
# for i in 4 5 6 7; do
#   sudo sh -c "echo 0 > /sys/devices/system/cpu/cpu$i/online";
# done

comms_info = [{
  command: 'pyasukpt -src INST_FILE',
  pattern: 'INST_FILE',
  extractor: PyaExtractor.new,
  prefix: 'PYAsUKP',
}, {
  command: 'run_ukp5.out INST_FILE',
  pattern: 'INST_FILE',
  extractor: UKP5Extractor.new,
  prefix: 'UKP5',
}]

execution_info = {
  cpus_available: [1, 2, 3],
  timeout: 10,
  post_timeout: 5,
}

conf = { csvfname: 'pya_site8.csv' }

files = ['corepb.ukp', 'exnsd18.ukp', 'exnsd26.ukp', 'exnsdbis18.ukp', 'exnsd16.ukp', 'exnsd20.ukp', 'exnsdbis10.ukp', 'exnsds12.ukp']
# If you don't execute the script from the ukp files folder you need to put the
# folder relative or absolute path here (with trailing slash).
path = ''
files.map! { | f | path + f }

experiment(comms_info, execution_info, conf, files)

