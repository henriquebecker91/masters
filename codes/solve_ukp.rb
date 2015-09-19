# Sort items on non-ascending order by profitability first and
# by non-descending weight after.
def sort_items_by_profitability!(items)
  items.sort_by! do | i |
    wi = i[:w]
    [Rational(wi,i[:p]), wi]
  end
end

# based on p. 223, Integer Programming, Robert S. Garfinkel
def y_star(items, already_sorted = false)
  items = sort_items_by_profitability!(items.clone) unless already_sorted
  i1, i2 = items[0], items[1]
  c1, c2, a1, a2 = i1[:p], i2[:p], i1[:w], i2[:w]
  
  (Rational(c1)/(Rational(c1, a1) - Rational(c2, a2))).ceil
end

# In truth, this removes simple and multiple dominance. It doesn't
# remove collective and threshold dominance.
def remove_simple_dominance(items, already_sorted = false)
  # the vector HAS TO BE sorted by non-increasing profitability first AND
  # non-decreasing weight after to this code work
  items = sort_items_by_profitability!(items.clone) unless already_sorted
  
  # If an item i dominate an item j, then pi/wi >= pj/wj. Note that NOT every
  # case we have an item i and item j and pi/wi >= pj/wj i dominates j
  # (if this was the case every UKP problem could be reduced to the best item).
  # This only means that, if i is the index of the item when the items are
  # ordered by pi/wi, then it can't be simple or multiple dominated by any
  # item of index j where j > i.
  undominated = [items[0]]
  # if this was C++ would be interesting to already reserve the n capacity
  # on the vector above
  items.each do | i |
    dominated = false
    wi = i[:w]
    pi = i[:p]
    undominated.each do | u |
      wu = u[:w]
      pu = u[:p]
      if Rational(wi,wu).floor * pu >= pi then
        dominated = true
        break
      end
    end
    undominated << i unless dominated
  end

  undominated
end

# Get the items used for the solution at capacity y when the g and d arrays
# that are used on the Garfinkel method are provided, in conjunction with the
# items sorted by the same order that was used by the method.
def get_used_items_for_y(y, items, g, d)
  used_items = {}
  dy = d[y]
  i = items[dy]
  wi = i[:w]
  while y - wi >= 0 do
    if used_items.has_key?(dy) then
      used_items[dy] += 1
    else
      used_items[dy] = 1
    end
    y -= wi
    dy = d[y]
    i = items[dy]
    wi = i[:w]
  end

  better_view = {}
  opt = 0
  used_items.each do | ix, qt |
    better_view[ix] = { qt: qt, p: items[ix][:p], w: items[ix][:w] }
    opt += items[ix][:p] * qt
  end

  better_view[:opt] = opt
  better_view
end

# Code to recover what capacity has the optimal solution when
# using the ukp5 method (the solution in g[y] of ukp5 can be a
# subotimal solution
def get_opt_y(c, items, g, d)
  w_min = items.min_by { | x | x[:w] }[:w]
  ix = c-w_min

  opt = g[ix]
  opt_y = ix
  ((ix+1)..c).each do | y |
    if g[y] > opt then
      opt = g[y]
      opt_y = y
    end
  end

  opt_y
end

# Another crazy ideia I had. This one is like the recursive implementation of
# the garfinkel method, but without the recursion. Maybe it's what the
# pyasukp creators call sparsity. 
# The ideia is: Instead of doing the recursity and beggining from the y
# to the lesser values of y, we begin at 0 and generate lower bounds for
# many values of y by using the items. 
def ukp5(ukpi, return_used_items = false)
  n = ukpi[:n]
  c = ukpi[:c]
  items = ukpi[:items].clone

  max_w = items.max_by { | i | i[:w] }[:w]
  g = Array.new(c+max_w+1, 0)
  d = Array.new(c+max_w+1, n-1)
  sort_items_by_profitability!(items)

  last_y_where_nonbest_item_was_used = 0
  items.each_with_index do | it, ix |
    wi = it[:w]
    pi = it[:p]
    if g[wi] < pi then
      g[wi] = pi
      d[wi] = ix
      last_y_where_nonbest_item_was_used = wi if wi > last_y_where_nonbest_item_was_used && ix != 0
    end
  end

  # This opt variable, and the "|| g[y] < opt" implements all the dominances,
  # the proof is simple:
  # (1) All the dominances (simple, multiple, collective and threshold) needs
  #     the existence of an item with better profitability than the dominated
  #     item.
  # (2) If opt > g[y] then there's a y' < y where g[y'] > g[y].
  # (3) The items are ordered by profitability.
  # (4) If g[y] is different from zero, then the d[y] is always the item of
  #     best profitability that was used on the solution of the knapsack with
  #     capacity y.
  # (5) All the dominances are special cases of a multiset of knapsack items
  #     being lighter and more profitable than another (i.e. better).
  #     The simple, multiple,
  #     and collective dominances are a multiset (respectively: of one, many of
  #     the same, and some differente items) being "better" than one single
  #     element. The threshold is the special case of some different items being
  #     better than a set of many of the same item.
  # (6) If we skip a y where g[y] != 0 then no capacity after y will
  #     use a multiset of items that include the multiset of items used
  #     to obtain g[y]. This happens because the solutions are build from
  #     inserting one more item on the current multiset. If we skip y
  #     no solution will be created from that multiset. Also, no multiset
  #     that is a subset of the multiset of y will follow a insertion route
  #     that allow it to become a superset of the multiset of y. Lets 
  #     suppose that the last item used to make the multiset of y was
  #     j, if instead of j we used i (to after use j) then the profitability
  #     of the whole multiset would be lower than the multiset of y, this is
  #     because of (3), ...
  # (X) When opt > g[y] there's a multiset of items that is better than
  #     another (2), what is a case that supersedes all dominances (5),
  #     by skipping that multiset we are skipping any solution that would
  #     have that multiset, 
  #     
  opt = 0
  (1..(c-1)).each do | y |
    next if g[y] == 0 || g[y] < opt
    break if last_y_where_nonbest_item_was_used < y

    opt = gy = g[y]
    dy = d[y]

    # this block is a copy-past of the loop bellow only for the best item
    bi = items[0]
    pb = bi[:p]
    wb = bi[:w]
    next_y = y + wb
    old_gny = g[next_y]
    new_gny = gy + pb
    if old_gny < new_gny then
      g[next_y] = new_gny
      d[next_y] = 0
    end

    (1..dy).each do | ix |
      it = items[ix]
      pi = it[:p]
      wi = it[:w]
      ny = y + wi
      ogny = g[ny]
      ngny = gy + pi
      if ogny < ngny then
        g[ny] = ngny
        d[ny] = ix
        last_y_where_nonbest_item_was_used = ny if ny > last_y_where_nonbest_item_was_used
      end
    end
  end

if last_y_where_nonbest_item_was_used < c-1 then
    y_ = last_y_where_nonbest_item_was_used
    while d[y_] != 0 do
      y_ += 1
    end
#    puts "Periodicity used - c: #{c}; last_y: #{y_}"

    extra_capacity = c - y_
    c1, a1 = items[0][:p], items[0][:w]
    qt_best_item_used = Rational(extra_capacity, a1).ceil
    space_used_by_best_item = qt_best_item_used*a1
    profit_generated_by_best_item = qt_best_item_used*c1

    opt_y = get_opt_y(c-space_used_by_best_item, items, g, d)
    g[c] = g[opt_y] + profit_generated_by_best_item
  end

  opt = g[c] if opt < g[c]

  if return_used_items then
    g.slice!(c+1, max_w)
    d.slice!(c+1, max_w)
    opt_y = get_opt_y(c, items, g, d)
    get_used_items_for_y(opt_y, items, g, d)
  else
    opt
  end
end

# Code based on p. 221, Integer Programming, Robert S. Garfinkel
# but using recursion
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

# Code for classical DP version using recursion
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

# Code based on p. 221, Integer Programming, Robert S. Garfinkel
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

# Code based on p. 221, Integer Programming, Robert S. Garfinkel
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

# Classical DP solution
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

# Code to read simplified UKP format (.sukp)
# Read the ukp instance from the file
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

# Read the ukp instance from the filename
def read_instance(filename)
  File.open(filename) do | f |
    read_ukp_instance(f)
  end
end

# Some toy instances
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
      { p: 14, w: 7},
      { p: 8, w: 4},
      # the following are simple dominated by the above
      { p: 8, w: 4},
      { p: 7, w: 4},
      { p: 7, w: 5},
      # the following item is multiple dominated by the item
      # on the second line
      { p: 15, w: 8},
      # the following item is collective dominated by the items
      # on the first and second line
      { p: 21, w: 11},
      # the following item is threshold dominated by the items on
      # the first and second line
      { p: 5, w: 3 }
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

