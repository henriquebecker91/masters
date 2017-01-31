#!/usr/bin/ruby

require 'require_relative'
require_relative 'realistic_random.rb'

# 80 number between 1 and 10^6 generated with random.org
seeds = [
  321446, 837176, 465040, 286473, 57635, 73266, 460116, 245741,
  90643, 36374, 56545, 770636, 776667, 458668, 104987, 529642,
  113608, 821582, 170676, 217121, 554787, 844842, 570795, 582092,
  266267, 522843, 225123, 803015, 421524, 619685, 256871, 594742,
  909073, 559133, 264409, 185965, 845341, 995156, 751140, 502089,
  926607, 446531, 903840, 742583, 512825, 574662, 732257, 639568,
  112983, 615173, 964192, 59271, 465129, 477565, 625841, 951699,
  651615, 895501, 650928, 906564, 12476, 575069, 776387, 746956,
  610976, 10474, 335984, 610281, 430955, 333973, 137962, 615086,
  750175, 608776, 157673, 593483, 18085, 927782, 961444, 31935,
]

folder = 'rr_new'
(0..7).to_a.each do | i |
  n    = (2**10) * (2**i)
  wmax = pmax = (2**10) * n
  wmin = pmin = wmax/(2**4)
  cmin = (2**1) * wmax
  cmax = cmin + wmin

  10.times do |j|
    seed = seeds[10*i + j]
    c = Random.new(seed).rand(cmin..cmax)
    content = realistic_random(n, wmin, wmax, pmin, pmax, c, seed)
    filename = "rr_n#{n}-#{j}-s#{seed}c#{c}.ukp"
    filepath = File.expand_path(filename, folder)
    File.open(filepath, 'w+') do | f |
      f.write(content)
    end
  end
end

