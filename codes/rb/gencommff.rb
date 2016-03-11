#!/usr/bin/ruby

# Created at 10/03/2016 by Henrique Becker
# gencommff: GENerate COMMands For Files
# Receives N strings. The first one is intended to be a command (but can be anything really). The second one is a pattern that happens in the first string. All the other strings will replace the pattern in the first string generating an output of N-2 strings. 
# Ex.: 'echo STR' STR a b c -> 'echo a' 'echo b' 'echo c'

def gencommff(comm, patt, files)
  ret = []
  files.each { | f | ret << comm.gsub(patt, f) }
  ret
end

def gencommsff(comms, patt, files, intercalate = true)
  ret = []
  if intercalate then
    files.each do | f |
      comms.each do | comm |
        ret << comm.gsub(patt, f)
      end
    end
  else
    comms.each do | comm |
      files.each do | file |
        ret << comm.gsub(patt, f)
      end
    end
  end
  ret
end

def experiment(comms, patt, repls, cpus, timeout, post_timeout)
end

log_lines = File.open(fname) { | f | f.read }.lines
# For when there's a field whose value is after '<field>: '.
def get_field(lines, field)
  lines.grep(/^#{field}: .*/).each { | l | l.match(/: (.*)/)[1] }
end

# For when there's a field whose value is in the next line.
def get_hfield(lines, field)
  lines[lines.find_index(field) + 1]
end

ukp5_labels = ['Seconds', 'ext_time', 'ext_mem', 'opt']
ukp5_gatherers = ukp5_labels.map { | l | lambda { | ls | get_field(ls, l) } }
pya_labels = ['Total Time ', 'ext_time', 'ext_mem']

def gather_field_info(lines, gatherers)
  ret = ''
  labels.each | l | do 
    ret << get_field(lines, l) << ';'
  end
  ret
end

