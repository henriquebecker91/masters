#!/usr/bin/ruby

def sort_items_by_profitability!(items)
  items.sort_by! { | i | Rational(i[:w],i[:p]) }
end

# Another crazy ideia I had. This one is like the recursive implementation of
# the garfinkel method, but without the recursion. Maybe it's what the
# pyasukp creators call sparsity. 
# The ideia is: Instead of doing the recursity and beggining from the y
# to the lesser values of y, we begin at 0 and generate lower bounds for
# many values of y by using the items. 
def ukp5(ukpi)
  n = ukpi[:n]
  c = ukpi[:c]
  items = ukpi[:items].clone

  max_w = items.max_by { | i | i[:w] }[:w]
  g = Array.new(c+max_w+1, 0)
  d = Array.new(c+max_w+1, n-1)
  sort_items_by_profitability!(items)

  items.each_with_index do | it, ix |
    wi = it[:w]
    pi = it[:p]
    if g[wi] < pi then
      g[wi] = pi
      d[wi] = ix
    end
  end

  (1..(c-1)).each do | y |
    next if g[y] == 0
    gy = g[y]
    dy = d[y]
    (0..dy).each do | ix |
      it = items[ix]
      pi = it[:p]
      wi = it[:w]
      ny = y + wi
      ogny = g[ny]
      ngny = gy + pi
      if ogny == 0 || ogny < ngny then
        g[ny] = ngny
        d[ny] = ix
      end
    end
  end

  g.slice!(c+1, max_w)
  g
end

def remove_simple_dominance(items)
#def remove_simple_dominance!(items, already_sorted = false)
#  items = sort_items_by_profitability!(items) unless already_sorted
  
  # adiciona um item arbitrário ao array de não dominados
  # a partir dai, para qualquer item do array original
  #   * verifica se não existe um elemento que domina ele já no array
  #     Se sim para aqui, senão continua
  #   * verifica se ele não domina algum elemento do array
  #     Se sim remove estes elementos do array
  #     Adiciona o elemento no array
  undominated = [items[0]]
  items.each do | i |
    dominated = false
    undominated.delete_if do | u |
      if i[:p] <= u[:p] && i[:w] >= u[:w] then
        dominated = true
        break
      elsif i[:p] >= u[:p] && i[:w] <= u[:w] then
        true
      else
        false
      end
    end
    undominated << i unless dominated
  end

  undominated
end

def ukp_rec2(ukpi)
  n = ukpi[:n]
  c = ukpi[:c]
  items = ukpi[:items].clone

  g = Array.new(c+1, nil)
  d = Array.new(c+1, n-1)
  sort_items_by_profitability!(items)

  ukp_aux2(items, c, g, d)

  g
end

def ukp_aux2(items, y, g, d)
  return g[y] if g[y]

  items.each_with_index do | item, j |
    wi = item[:w]
    pi = item[:p]
    if y - wi >= 0 && j <= d[y - wi] then
      v = (ukp_aux2(items, y - wi, g, d) or 0)
      v += pi
      if g[y].nil? || v > g[y] then
        g[y] = v
        d[y] = j
      end
    end
  end

  g[y]
end

def ukp_rec(ukpi)
  c = ukpi[:c]
  items = ukpi[:items]

  memo = Array.new(c+1, nil)
  items = items.sort_by { | i | i[:w] }

  ukp_aux(items, c, memo)

  memo
end

def ukp_aux(items, y, memo)
  return memo[y] if memo[y]

  opt = 0
  items.each do | i |
    break if i[:w] > y

    x = i[:p] + ukp_aux(items, y - i[:w], memo)
    opt = x if opt < x
  end

  memo[y] = opt
end

# Based on one crazy ideia I had. Probably is equivalent to the method of
# garfinkel. The ideia is: If you know the item with the best profitability,
# you know that every capacity y that is a perfect divisible by w_1 have
# the best filling using x_1 items where x_1 = y/w_1. That is indiscutible.
# After that, we know that every capacity y, that is perfect divisible by w_2,
# have the best filling using x_2 + x_1 where x_1 is the biggest possible
# value and x_2*a_2 + x_1*a_1 = y (not less or equal, but equal). And so on.
# NOTE: This don't solve the UKP for c, its solves the ukp for an arbitrary
# group of capacities that are easy to calculate (are perfect divisible for
# some item weight).
def ukp4(ukpi)
  n = ukpi[:n]
  c = ukpi[:c]
  items = ukpi[:items].clone

  g = Array.new(c+1, 0)
  sort_items_by_profitability!(items)

  items.each do | item |
    pi = item[:p]
    wi = item[:w]
    max_xi = Rational(c, wi).floor
    z = pi
    y = wi
    max_xi.times do
      if g[y] == 0 then
        g[y] = z
      else
        z = g[y]
      end
      z += pi
      y += wi
    end

    break if g[c] != 0
  end

  g
end

# based on p. 221, Integer Programming, Robert S. Garfinkel
# adding the use of y* introduced on page 223
def ukp3(ukpi, ret_table = true)
  n = ukpi[:n]
  c = ukpi[:c]
  items = ukpi[:items].clone

  # step one
  g = Array.new(c+1, 0)
  d = Array.new(c+1, n-1)
  sort_items_by_profitability!(items)

  gy_col = nil
  dy_col = nil
  if ret_table then
    gy_col = Array.new(c+1) { [] }
    dy_col = Array.new(c+1) { [] }
  end

  ys = y_star(items, true)
  min_necessary_steps = ys < c+1 ? ys : c+1
  # step two to six
  min_necessary_steps.times do | y |
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

  if min_necessary_steps != c+1 then
    extra_capacity = c - ys
    c1, a1 = items[0][:p], items[0][:w]
    qt_best_item_used = Rational(extra_capacity, a1).ceil
    space_used_by_best_item = qt_best_item_used*a1
    profit_generated_by_best_item = qt_best_item_used*c1
    g[c] = g[c-space_used_by_best_item] + profit_generated_by_best_item
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

# based on p. 221, Integer Programming, Robert S. Garfinkel
def ukp2(ukpi, ret_table = true)
  n = ukpi[:n]
  c = ukpi[:c]
  items = ukpi[:items].clone

  # step one
  g = Array.new(c+1, 0)
  d = Array.new(c+1, n-1)
  sort_items_by_profitability!(items)

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
      # the following are simple dominated by the above
      { p: 8, w: 4},
      { p: 7, w: 4},
      { p: 7, w: 5},
    ]
  }
end

def my_instance2
  {
    n: 2,
    c: 10,
    items: [
      { p: 9, w: 3 },
      { p: 20, w: 7 },
    ]
  }
end

def read_instance(filename)
  File.open(filename) do | f |
    read_ukp_instance(f)
  end
end

# based on p. 223, Integer Programming, Robert S. Garfinkel
def y_star(items, already_sorted = false)
  items = sort_items_by_profitability!(items.clone) unless already_sorted
  i1, i2 = items[0], items[1]
  c1, c2, a1, a2 = i1[:p], i2[:p], i1[:w], i2[:w]
  
  (Rational(c1)/(Rational(c1, a1) - Rational(c2, a2))).ceil
end

=begin
puts ARGV[0]
puts ukp2(read_instance(ARGV[0])).last
=end
=begin
instance = read_instance(ARGV[0])
puts y_star(instance[:items])
puts instance[:c]
puts Rational(instance[:c],y_star(instance[:items]))
=end

=begin
instance = read_instance(ARGV[0])
puts instance[:items].length
instance[:items] = remove_simple_dominance(instance[:items])
puts instance[:items].length
puts ukp3(instance, false).last
=end

#puts remove_simple_dominance(my_instance[:items])
#puts ukp2(my_instance, false).last
instance = read_instance(ARGV[0])
=begin
puts 'ukp'
t = Time.now
v = ukp(instance).last
puts Time.now - t
puts v
puts "---"
puts 'ukp2'
t = Time.now
v = ukp2(instance, false).last
puts Time.now - t
puts v
puts "---"
t = Time.now
puts 'ukp3'
v = ukp3(instance, false).last
puts Time.now - t
puts v
puts "---"
=end
=begin The recursive functions are overflowing the stack
t = Time.now
puts 'ukp_rec'
v = ukp_rec(instance).last
puts Time.now - t
puts v
puts "---"
t = Time.now
puts 'ukp_rec2'
v = ukp_rec2(instance).last
puts Time.now - t
puts v
puts "---"
=end
puts 'ukp5'
t = Time.now
v = ukp5(instance).last
puts Time.now - t
puts v

#puts sort_items_by_profitability!(garfinkel_instance[:items].clone)
#puts y_star(garfinkel_instance[:items])

