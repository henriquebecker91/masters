comms_info = [{
  # String with command to be executed. Must have 'pattern' as substring.
  command: 'sleep 1 && echo X',
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
}

files = ['corepb.ukp', 'exnsd18.ukp', 'exnsd26.ukp', 'exnsdbis18.ukp', 'exnsd16.ukp', 'exnsd20.ukp', 'exnsdbis10.ukp', 'exnsds12.ukp']
path = ''
files.map! { | f | path + f }

experiment(comms_info, execution_info, conf, files)

