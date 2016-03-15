#!/usr/bin/ruby

require 'require_relative'
require_relative './info_extract.rb'
require_relative './batch_executer/batch.rb'
# Created at 10/03/2016 by Henrique Becker
# gencommff: GENerate COMMands For Files
# Receives N strings. The first one is intended to be a command (but can be anything really). The second one is a pattern that happens in the first string. All the other strings will replace the pattern in the first string generating an output of N-2 strings. 
# Ex.: 'echo STR' STR a b c -> 'echo a' 'echo b' 'echo c'

def gencommff(comm, patt, files)
  ret = []
  files.each { | f | ret << comm.gsub(patt, f) }
  ret
end

#def gencommsff(comms, patt, files, intercalate = true)
#  ret = []
#  if intercalate then
#    files.each do | f |
#      comms.each do | comm |
#        ret << comm.gsub(patt, f)
#      end
#    end
#  else
#    comms.each do | comm |
#      files.each do | file |
#        ret << comm.gsub(patt, f)
#      end
#    end
#  end
#  ret
#end

#def gather_from_lines(lines, gatherers, sep)
#  ret = []
#  gatherers.each | g | do 
#    ret << g(lines) << sep
#  end
#  ret.join
#end

#def gather_from_file(file, gatherers, sep)
#  file_lines = File.open(fname) { | f | f.read }.lines.to_a
#  gather_from_lines(file_lines, gatherers, sep)
#end

#def gather_from_files(files, gatherers, sep)
#  files.map do | file |
#    gather_from_file(file, gatherers)
#  end
#end

comm_info = {
  # String with command to be executed. Must have patt as substring.
  command: '',
  # Substring present in 'command'. Often replaced by the instance filename.
  pattern: '',
  # Extractor subclass object. Receives the output of the command and return
  # the most important fields.
  extractor: ,
  # String used to identify the command. Will be used to prefix the return of
  # extractor.names.
  prefix: '',
}

execution_info = {
  # IDs of the CPU cores that can be used for executing tests.
  cpus_available: [],
  # Maximum number of seconds that a command can run. After this a kill command
  # will be issued.
  timeout: ,
  # Maximum number of seconds that a command can run after a kill command was
  # issued. After this a kill -9 command will be issued.
  post_timeout: ,
}

conf = {
  # skip_commands # If true, will not execute the commands and assume that the outputs are already saved. Will only execute the gatherer over the saved outputs.
  # separator # The string that will be used as separator at the CSV.
  # csvfname: String (valid as a filename).
  # The name of the file where will be written the CSV file with all the
  # info gathered.
}

# Only works for arrays of the same size.
def intercalate(xss)
  group_size = xss.first.size
  num_groups = xss.size
  ret = Array.new(group_size * num_groups)
  xss.each_with_index do | xs, i |
    xs.each_with_index do | x, j |
      ret[(j * (group_size - 1)) + i] = x
    end
  end
  ret
end

# Only works for arrays of the same size.
#def regroup(xs, num_groups, group_size)
#  ret = Array.new(num_groups) { Array.new(group_size) }
#  group_size.times do | j |
#    num_groups.times do | i |
#      ret[i][j] = xs[(j * (group_size-1)) + i]
#    end
#  end
#  ret
#end

def experiment(comms_info, execution_info, conf, files)
  command_sets = []
  comms_info.each do | comm_info |
    command_sets << gencommsff(comm_info[:command], comm_info[:pattern], files)
  end
  intercalated_comms = intercalate(command_sets)
  
  output_files = batch(intercalated_comms, execution_info[:cpus_available], execution_info[:timeout], execution_info[:post_timeout])

  header = ['Filename']
  comms_info.each do | comm_info |
    prefixed_names = comm_info[:extractor].names.map do | name |
      (comm_info[:prefix] + ': ') << name
    end
    header.concat(prefixed_names)
  end
  header.join(conf[:separator])

  body = [header]
  files.each_with_index do | inst_fname, j |
    line = [inst_fname]
    comms_info.each_with_index do | comm_info, i |
      output_fname = output_files[j * (files.size-1)) + i]
      f_content = File.open(output_fname) { | f | f.read }

      extractor = comm_info[:extractor]
      line.concat(extractor.extract(f_content))
    end
    body << line.join(conf[:separator])
  end
  body.join("\n")

  File.new(conf[:csvfname], 'w') do { | f | f.write(body) } 
end

