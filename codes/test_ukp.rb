#!/usr/bin/ruby

require 'require_relative'
require_relative 'solve_ukp.rb'

instances_data = [
#  { name: 'corepb', opt:  10077782 }, # too slow
  { name: 'exnsd16', opt: 1029680 },
  { name: 'exnsds12', opt: 3793952 },
  { name: 'teste2', opt: 225092 },
  { name: 'teste', opt: 135000 },
]

everything_ok = true
instances_data.each do | id |
  opt = id[:opt]
  path = "../data/sukp/#{id[:name]}.sukp"
  puts "Path: #{path}"
#  puts "Expected value: #{opt}"
  instance = read_instance(path)

  10.times do
    t = Time.now
    v = ukp5(instance)
  #  puts "Obtained value: #{v}"
  #  puts "y*: #{y_star(instance[:items])}"
    puts "Time used (in seconds): #{Time.now - t}"
  end
  puts
#  everything_ok = everything_ok && v == opt
end

if everything_ok then
  puts 'Everything is OK'
else
  puts 'The known optimal and the optimal obtained for some of the instances vary, something is wrong.'
end

