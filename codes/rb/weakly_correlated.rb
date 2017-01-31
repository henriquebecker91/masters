#!/usr/bin/ruby

# Generic weakly correlated instance generator.
def weakly_correlated(n, wmin, wmax, alpha, c, seed)
  str = "# weakly correlated -- n: #{n}; wmin: #{wmin}; wmax: #{wmax}; alpha: #{alpha}; c: #{c}; seed: #{seed}\n"
  str << "n: #{n}\nc: #{c}\n"
  str << "begin data\n"

  rng = Random.new(seed)
  w = (wmin..wmax).to_a.sample(n, random: rng)
  p = w.map{ |wi| wi + rng.rand((-alpha)..alpha) }
  n.times { |i| str << "#{w[i]} #{p[i]}\n" }
  str << "end data\n"

  str
end


