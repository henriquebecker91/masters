#!/usr/bin/ruby

# based on p. 221, Integer Programming, Robert S. Garfinkel
def ukp2(ukpi)
  n = ukpi[:n]
  c = ukpi[:c]
  items = ukpi[:items]

  # step one
  g = Array.new(c+1, 0)
  d = Array.new(c+1, n)
  items = items.sort_by { | i | i[:p]/i[:w] }
  y = 1
  j = 0

  # step two to six
  loop do
    print "\n#{g[y]} "
    loop do
      #step two
      if y - items[j][:w] >= 0 && j <= d[y - items[j][:w]] then
        #step three
        v = g[y - items[j][:w]] + items[j][:p]
        if v > g[y] then
          g[y] = v
          print "#{g[y]} "
          d[y] = j
        end
      end
      #step four
      break unless j < n-1 
      j += 1
    end
    #step five
    if y == c then
      break
    else
      j = 0
      y += 1
    end
  end
end

def ukp(ukpi)
  c = ukpi[:c]
  items = ukpi[:items]

  opt = Array.new(c+1, 0)
  items = items.sort_by { | i | i[:w] }

  (c+1).times do | y |
    items.each do | i |
      break if i[:w] > y
      x = i[:p] + opt[y - i[:w]]
      if opt[y] < x then
        opt[y] = x
      end
    end
  end

  opt
end

def read_ukp_instance(f)
  all_words = f.read.scan(/\w+/)
  ukpi = {
    n: all_words[0].to_i,
    c: all_words[1].to_i,
    items: []
  }
  fail "ukp instance format wrong - n: #{ukpi[:n]} - all_words.length: #{all_words.length} " unless (ukpi[:n]*2)+2 == all_words.length

  i = 2
  size = all_words.length
  while i <= size do
    item = { 
      p: all_words[i].to_i,
      w: all_words[i+1].to_i
    }
    ukpi[:items] << item
    i += 2
  end

  ukpi
end

def garfinkel_instance
  ukpi = {
    n: 4,
    c: 25,
    items: [
      { p: 11, w: 6},
      { p: 7, w: 4},
      { p: 5, w: 3},
      { p: 1, w: 1}
    ]
  }
end

def read_instance(filename)
  File.open(filename) do | f |
    ukpi = read_ukp_instance(f)
  end
end

#puts ARGV[0]
#puts teste2(ARGV[0])
puts ukp2(garfinkel_instance)

