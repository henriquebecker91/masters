#!/usr/bin/ruby

# Created at 10/03/2016 by Henrique Becker
# gencommff: GENerate COMMands For Files
# Receives N strings. The first one is intended to be a command (but can be anything really). The second one is a pattern that happens in the first string. All the other strings will replace the pattern in the first string generating an output of N-2 strings. 
# Ex.: 'echo STR' STR a b c -> 'echo a' 'echo b' 'echo c'

def gencommff(s, patt, repls)
  ret = []
  repls.each { | repl | ret << s.gsub(patt, repl) }
  ret
end

def gencommsff(ss, patt, repls, intercalate = true)
  ret = []
  if intercalate then
    repls.each do | repl |
      ss.each do | s |
        ret << s.gsub(patt, repl)
      end
    end
  else
    ss.each do | s |
      repls.each do | repl |
        ret << s.gsub(patt, repl)
      end
    end
  end
  ret
end

# Put this on the scheduler, to not pollute the command name
def addtime(comm)
  "time -f 'Ext_time: %e\\nExt_mem: %M\\n' #{comm}"
end

