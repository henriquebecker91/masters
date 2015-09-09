#!/usr/bin/ruby

# based on p. 221, Integer Programming, Robert S. Garfinkel
def ukp2(ukpi, ret_table = true)
  n = ukpi[:n]
  c = ukpi[:c]
  items = ukpi[:items]

  # step one
  g = Array.new(c+1, 0)
  d = Array.new(c+1, n-1)
  items = items.sort_by { | i | i[:p]/i[:w] }

  gy_col = nil
  dy_col = nil
  if ret_table then
    gy_col = Array.new(c+1) { [] }
    dy_col = Array.new(c+1) { [] }
  end

  # step two to six
  (c+1).times do | y |
    if ret_table then
      gy_col[y] << g[y]
      dy_col[y] << d[y]
    end
    items.each_with_index do | item, j |
      if y - item[:w] >= 0 && j <= d[y - item[:w]] then
        v = g[y - item[:w]] + item[:p]
        if v > g[y] then
          g[y] = v
          d[y] = j
          if ret_table then
            gy_col[y] << g[y]
            dy_col[y] << d[y]
          end
        end
      end
    end
  end

  table = "y\tg(y)\t\td(y)\n"
  if ret_table then
    (c+1).times do | y |
      table << "#{y}\t#{gy_col[y].join(' ')}\t\t#{dy_col[y].join(' ')}\n"
    end
  end
  
  if ret_table then
    table
  else
    g
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
  while i+1 <= size do
    item = { 
      w: all_words[i].to_i,
      p: all_words[i+1].to_i
    }
    ukpi[:items] << item
    i += 2
  end

  ukpi
end

def garfinkel_instance
  {
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

def my_instance
  {
    n: 3,
    c: 11,
    items: [
      { p: 21, w: 10},
      { p: 14, w: 7},
      { p: 8, w: 4},
    ]
  }
end

def read_instance(filename)
  File.open(filename) do | f |
    read_ukp_instance(f)
  end
end

puts ARGV[0]
puts ukp2(read_instance(ARGV[0])).last
#puts ukp2(garfinkel_instance)
#puts ukp2(my_instance, false)

