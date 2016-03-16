require 'require_relative'
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

# Only works for arrays of the same size.
def intercalate(xss)
  group_size = xss.first.size
  num_groups = xss.size
  ret = Array.new(group_size * num_groups)
  xss.each_with_index do | xs, i |
    xs.each_with_index do | x, j |
      ret[(j * num_groups) + i] = x
    end
  end
  ret
end

def experiment(comms_info, execution_info, conf, files)
  command_sets = []
  comms_info.each do | comm_info |
    command_sets << gencommff(comm_info[:command], comm_info[:pattern], files)
  end
  intercalated_comms = intercalate(command_sets)

  output_files = batch(intercalated_comms, execution_info[:cpus_available], execution_info[:timeout], execution_info[:post_timeout])

  header = []
  comms_info.each do | comm_info |
    prefixed_names = comm_info[:extractor].names.map do | name |
      (comm_info[:prefix] + ' ') << name
    end
    header << prefixed_names
  end
  header = ['Filename'].concat(intercalate(header)).join(conf[:separator])

  body = [header]
  files.each_with_index do | inst_fname, j |
    line = []
    comms_info.each_with_index do | comm_info, i |
      output_fname = output_files[(j * comms_info.size) + i]
      f_content = File.open(output_fname + '.out') { | f | f.read }

      extractor = comm_info[:extractor]
      line << extractor.extract(f_content)
    end
    
    body << [inst_fname].concat(intercalate(line)).join(conf[:separator])
  end
  body = body.map! { | line | line << ';' }.join("\n")

  File.open(conf[:csvfname], 'w') { | f | f.write(body) }
end

