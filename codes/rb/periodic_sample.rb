#!/usr/bin/ruby


# Snippet of code used to understand better how periodicity works.
# Values for example: 11 6 & 7 4
puts "Enter best item profit:"
p1 = gets.to_i
puts "Enter best item weight:"
w1 = gets.to_i
puts "Enter second-best item profit:"
p2 = gets.to_i
puts "Enter second-best item weight:"
w2 = gets.to_i

y = [w1, w2].min
loop do
  # A valid value of solution using only the first element (lower bound)
  o1 = p1 * Rational(y, w1).floor
  # A not guaranteed to be valid upper bound without the best element
  o2 = Rational(p2*y, w2).ceil
  reached = o1 > o2
  op = reached ? '>' : '<='
  puts "y: #{y}\t\t#{p1}*[#{y}/#{w1}] #{op} #{p2}*#{y}/#{w2}\t\t#{o1} #{op} #{o2}"

  break if reached
  y += 1
end

