#!/usr/bin/ruby

require 'require_relative'
require_relative 'solve_ukp.rb'

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

=begin
puts instance[:items].length
u = remove_simple_dominance(instance[:items])
puts u.length
=end

#puts sort_items_by_profitability!(garfinkel_instance[:items].clone)
#puts y_star(garfinkel_instance[:items])
#puts ukp2(my_instance, false).last
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

