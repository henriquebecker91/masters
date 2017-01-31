#!/usr/bin/ruby

# Realistic random instances. As described in page 21 of the "UNBOUNDED KNAPSACK PROBLEM: DYNAMIC PROGRAMMING REVISITED" (not the article, but the internal report N o PI-1152, with 26 pages and that is the first citation in the article).
def realistic_random(n, wmin, wmax, pmin, pmax, c, seed)
  str = "# realistic random -- n: #{n}; wmin: #{wmin}; wmax: #{wmax}; pmin: #{pmin}; pmax: #{pmax}; c: #{c}; seed: #{seed}\n"
  str << "n: #{n}\nc: #{c}\n"
  str << "begin data\n"

  rng = Random.new(seed)
  w = (wmin..wmax).to_a.sample(n, random: rng).sort!
  p = (pmin..pmax).to_a.sample(n, random: rng).sort!
  wp = []
  n.times { |i| wp << "#{w[i]} #{p[i]}" }
  str << wp.shuffle!(random: rng).join("\n") << "\n"
  str << "end data\n"

  str
end

