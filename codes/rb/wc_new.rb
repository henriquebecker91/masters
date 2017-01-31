#!/usr/bin/ruby

require 'require_relative'
require_relative 'weakly_correlated.rb'

# 80 number between 1 and 10^6 generated with random.org
seeds = [
  78961, 340121, 536975, 654779, 616759, 10073, 439442, 691754,
  255313, 92787, 214757, 840590, 660641, 183647, 640745, 896888,
  834596, 399531, 348379, 824748, 274868, 54056, 86112, 181876,
  365333, 343307, 54630, 974597, 350955, 304607, 105451, 202656,
  481858, 638403, 805035, 16755, 354971, 712381, 416795, 312207,
  199129, 981914, 727235, 720658, 307127, 788375, 863873, 17143,
  348889, 964122, 247818, 890585, 547121, 269047, 92151, 827396,
  744439, 60112, 670598, 744196, 496939, 861568, 743485, 537355,
  207758, 279743, 985715, 923897, 928260, 785591, 511491, 350718,
  964911, 639393, 611064, 791359, 135182, 168000, 585243, 59713,
]

folder = 'wc_new'
(0..7).to_a.each do | i |
  n     = (2**10) * (2**i)
  wmax  = (2**10) * n
  wmin  = wmax/(2**4)
  cmin  = (2**1) * wmax
  cmax  = cmin + wmin
  alpha = n/(2**3)

  10.times do |j|
    seed = seeds[10*i + j]
    c = Random.new(seed).rand(cmin..cmax)
    content = weakly_correlated(n, wmin, wmax, alpha, c, seed)
    filename = "wc_n#{n}-#{j}-s#{seed}c#{c}.ukp"
    filepath = File.expand_path(filename, folder)
    File.open(filepath, 'w+') do | f |
      f.write(content)
    end
  end
end

