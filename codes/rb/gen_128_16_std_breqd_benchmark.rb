#!/usr/bin/ruby

# Generates a bottom right ellipse quadrant shaped item distribution, where the
# lowest point is (1,1) and the highest is (wmax,pmax). 
# @note: The profit of an item is always computed based on its weight.
#   The parameter pmin isn't implemented, to avoid making the method more
#   complex. The best way to have an arbitrary pmin, while keeping the
#   distribution properties is to compute the items for pmax - pmin, and then
#   add pmin to the profit of each item.
# @note: If the profit values were floating point number with a large
#   (or infinite precision), there wouldn't be any simple/multiple/collective
#   dominated items, but as we use integers, small items can end up with the
#   same profit value because of rounding, allowing for some dominance.
# @note: For the programmers not familiar with ruby, the method below is using
#   infinite precision integers (ruby behaviour by default).
# @note: The seed is used to initialize a mersenne twister generator, that will
#   feed a sample method (not very reproducible, but avoid duplicated weights,
#   and avoid skewing the distribution when converting weights to the desired
#   range).
def gen_std_breqd_instance(n, c, wmax, pmax, seed)
  r = pmax/wmax # integer division
  r_2 = r**2
  pmax_2 = pmax**2
  rng = Random.new(seed)
  s = "# bottom right ellipse quadrant distribution-- n: #{n}; c: #{c}; wmax: #{wmax}; pmax: #{pmax}; seed: #{seed}\n"
  s << "n: #{n}\nc: #{c}\nbegin data\n"
  (1..wmax).to_a.sample(n, random: rng).each do | w |
    p = pmax - Math.sqrt(pmax_2 - (w**2 * r_2)).to_i
    s << "#{w} #{p}\n"
  end
  s << "end data\n"
  s
end

# The rationale of this method is avoid controlling for many variables, the
# other parameters for the instance generation are computed based on n.
# We do not allow for items bigger than c (i.e. wmax > c) as this would
# generate items that would not fit the knapsack.
# The parameter wmax will be equal to the capacity.
# The idea is an instance easy for B&B and hard for DP. The optimal solution
# probably will contain the best item, that will be, at the same time, the heaviest item, the most profitable item, and the most efficient item. This item will be found immediatly by B&B methods, that will begin testing solutions using it, and will need only to to fill a small gap left by it (excluding almost all items by the computation of bounds). DP methods will find this item only at end of the execution and, as the item efficiency increases with the weight, the DP methods can't exploit simple/multiple/collective dominance to speed-up computation (only threshold dominance, of the bigger items over solutions made of the smaller items).
# The capacity c has to be some orders of magnitude bigger than n. If c < n,
# there will be repeated items in an instance (wmax <= c < n, pidgeonhole). If c isn't many times greater than n, many items will be shared between instances generated with the same parameter n, also, we want DP algorithms to suffer because of the magnitude of c, that has little effect over B&B algorithms. So we opt for c = 128*n, this way less than 1% of the possible items for the interval wmin...wmax=c will be generated (values bigger than 128*wmax would be interesting, but could end up on capacity values that are too big for DP algorithms allocate sufficient RAM, what is not the idea). Finally, it's interesting to have the range pmin-pmax to be at least one order of magnitude greater than the range wmin-wmax. Otherwise, many items would dominate each other simply because the rounding of the exact value of p (that would be a floating point value) lose significant information. With pmax = 16*wmax, we have the first items with efficiency almost zero (smaller than one at least), and the efficiency will increases until almost 16.
def gen_128_16_std_breqd_instance(n, seed)
  c    = 128*n
  wmax = c
  pmax = 16*wmax
  gen_std_breqd_instance(n, c, wmax, pmax, seed)
end

# Generates one hundred instances, ten with n = 2^11, ten with n = 2^12, until
# ten with n = 2^20. Those instances follow the distribution defined by
# gen_bb_easy_dp_hard_std (i.e. c = 128*n, wmax = c, pmax = 16*wmax). The seeds
# vary between 0 and 9, for each ten instances with the same n. The 2*n
# instance (ex.: 4096) isn't guaranteed to have the 2048 weights present in the
# instance n (ex.: 2048), but it will probably share some values by simple
# probability.
def gen_128_16_std_breqd_benchmark(folder)
  10.times do | size |
    10.times do | seed |
      # 2048 * (2**0) = 2048 = 2**11 ~= 2*10^3
      # ...
      # 2048 * (2**9) = 1048576 = 2**20 ~= 10^6
      n = 2048 * (2**size)
      inst_str = gen_128_16_std_breqd_instance(n, seed)
      filepath = File.expand_path("128_16_std_breqd-n#{n}-s#{seed}.ukp", folder)
      File.open(filepath, 'w+') do | f |
        f.write(inst_str)
      end
    end
  end
end

# the original hundred instances were generated with
gen_128_16_std_breqd_benchmark('128_16_std_breqd_benchmark')
# on ruby 2.3.1p112 (2016-04-26 revision 54768) [x86_64-linux]

# Example of a case were the best item would not be in the solution:
#
## Utility method to get the profit of an item with weight w.
## Useful to do exploratory analysis.
#def p_of_w(w, wmax, pmax)
#  r = pmax/wmax # integer division
#  r_2 = r**2
#  pmax_2 = pmax**2
#
#  pmax - Math.sqrt(pmax_2 - (w**2 * r_2)).to_i
#end
#
#def print_p_of_w(w, wmax, pmax)
#  p = p_of_w(w, wmax, pmax)
#  puts "w: #{w} p: #{p} e: #{p.to_f / w.to_f}"
#end
#
# [384, 383, 129, 32].each { | x | print_p_of_w(x, 4*128, 4*128*16) }
# w: 384 p: 2774 e: 7.223958333333333
# w: 383 p: 2756 e: 7.195822454308094
# w: 129 p: 265 e: 2.054263565891473
# w: 32 p: 17 e: 0.53125
#
# With c = 512, the solution don't use the best item. If we use the best item
# we get:
# 384 + 4*32 = 512 (2774 + 4*17 = 2842) while
# 383 + 129 = 512 (2756 + 265 = 3021).

