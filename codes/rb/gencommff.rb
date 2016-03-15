#!/usr/bin/ruby

require 'require_relative'
require_relative './info_extract.rb'
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

def gather_from_lines(lines, gatherers, sep)
  ret = []
  gatherers.each | g | do 
    ret << g(lines) << sep
  end
  ret.join
end

def gather_from_file(file, gatherers, sep)
  file_lines = File.open(fname) { | f | f.read }.lines.to_a
  gather_from_lines(file_lines, gatherers, sep)
end

def gather_from_files(files, gatherers, sep)
  files.map do | file |
    gather_from_file(file, gatherers)
  end
end

def Gatherer
  # Return the field names for each of the elements returned by
  # extract_field_values.
  def field_names
    fail 'This method should have been overwritten by a subclass.'
  end

  # Extract from the command output an array of values. This array has the same
  # size as the one returned by field_names.
  def extract_field_values(lines)
    fail 'This method should have been overwritten by a subclass.'
  def 
end

comms_info = {
  # Commands to be executed (array of strings without linebreaks).
  commands: [],
  # Substring present in all elements of 'commands'.
  # Will be replaced by the instance filename in every command.
  patt: '',
  # Array of objects with a call method.
  # All this objects will be called over the command output. Their returns will
  # be merged as a line in a CSV file.
  gatherers: [],
  # Array of strings without linebreaks.
  # Will be written at the top of the file that aggregates the info gathered
  # from all the commands outputs.
  header: [],
  # Separator to be used between the returns of the gatherers. And after the
  # last return.
  # The default is a semicolon.
  sep: '',
}

execution_info = {
  # IDs of the CPU cores that can be used for executing tests.
  cpus_available: [],
  # Maximum number of seconds that a command can run. After this a kill command will be issued.
  timeout: ,
  # Maximum number of seconds that a command can run after a kill command was issued. After this a kill -9 command will be issued.
  post_timeout: ,
}

  # logfname: String (valid as a filename).
  # The name of the file where will be written the CSV file with all the
  # info gathered.
def experiment(comms_info, execution_info, logfname, files)
  
end

