#!/usr/bin/ruby

(0..7).to_a.each do | i |
  n = (2**10) * (2**i)
  wmax = pmax = (2**10) * n
  wmin = pmin = wmax/(2**4)
  cmin = (2**1) * wmax
  cmax = cmin + wmin
  alpha = n/(2**3)
  puts <<-END
  i           = #{i}
  n           = (2**10) * (2**i) = #{n}
  wmax = pmax = (2**10) * n      = #{wmax}
  wmin = pmin = wmax/(2**4)      = #{wmin}
  cmin        = (2**1) * wmax    = #{cmin}
  cmax        = cmax + wmin      = #{cmax}
  alpha       = n/(2**3)         = #{alpha}

  END
end

