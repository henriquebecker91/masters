#!/usr/bin/ruby

# How to generate a worst case performance instance for B&B
#
# p1/w1 >= p2/w2 >= ... >= pn/wn
# (the index used to refer to the items follows the items efficiency order)
#
# We use [expression] to mean "the rounded down value returned by expression".
#
# opt sol: [c/wn] copies of item n
#
# Any instance where the optimal solution is comprised only of copies of the
# least efficient item is an instance with worst case-performance for B&B; as
# those methods will often begin with a solution comprised only of copies of
# the most efficient item and will gradually remove the most efficient items to
# check if the solution improves with copies of least efficient items (they use
# the space better).
#
# If pn*[c/wn] is the optimal solution value, then for any j (1..n-1),
# pn*[c/wn] >= pj*[c/wj]. Therefore, solutions comprised only of copies of a
# single item more efficient than n should waste so much space on c (i.e. c -
# [c/wj]*wj is sufficient big) to be dominated by the solution comprised only
# of copies of n. It's clear that c-[c/wj]*wj must be closer to zero than c -
# any-non-opt-solution-weight. The zero value can be achieved with wj == c, or
# the more general case were: (c mod wj) == 0. A simple way of testing this
# concept is having a set of items where all items have [c/wj] == 1 (in other
# words wj > c/2), and the weight of the item grows inversely to its efficiency
# (to a max of c) in a superlinear way.

# (j. 1..n) [c/wj] == 1
# p1 < p2 < ... < pn
# w1 < w2 < ... < wn
# p1/w1 > p2/w2 > ... pn/wn

n = 10000

#if 2*n + ((c/2) + 1) > c then
#  puts 'Isn\'t possible to use this formulae with those values of n and c. You should use a bigger c '
#end

c = n*n

puts "n: #{n}"
puts "c: #{c}"

#w = (c/2)+1
#p = w
#puts "begin data"
#n.times do
#  puts "#{w} #{p}"
#  w += 2
#  p += 1
#end
#puts "end data"

puts "begin data"
(1..n).each do | i |
  w = (c/i) + 1
  p = w - i
  puts "#{w} #{p}"
end
puts "end data"

