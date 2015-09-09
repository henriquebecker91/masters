#!/usr/bin/ruby

# putting the variables on the global scope
$n = nil
$c = nil
$p = ''
$w = ''

File.open(ARGV[0], 'r') do | f |
  $n = f.gets.to_i
  $c = f.gets.to_i

  # skip the comments (if any)
  i = 1
  while (i <= $n && line = f.gets)
    line_w = line.scan(/\w+/)
    $w << "#{i} #{line_w[0]}, "
    $p << "#{i} #{line_w[1]}, "
    i += 1
  end
end
# remove the last ", " (conma and space)
$p.chop!.chop!
$w.chop!.chop!

print <<END
/* UKP (Unbounded Knapsack Problem) */

param n, integer, > 0;
/* number of items */

param c, integer, > 0;
/* capacity of the knapsack */

param p{i in 1..n}, integer, > 0;
/* profit of putting the item on the knapsack */

param w{i in 1..n}, integer, > 0;
/* weight of the item  */

var x{i in 1..n}, integer, >= 0;
/* how many of each item to put on the knapsack */

maximize obj: sum{i in 1..n} p[i]*x[i];
/* objective is to maximize the profit of the items on the knapsack */

s.t. cap_con: sum{i in 1..n} w[i]*x[i] <= c;
/* we have to respect the knapsack limit */

solve;

for {i in 1..n} {
  printf 'Obj nº %d p %d w %d: %d\\n', i, p[i], w[i], x[i];
}
printf "Max Profit: %d\\n", sum {i in 1..n} p[i]*x[i];

data;

param n := #{$n};

param c := #{$c};

param p := #{$p};

param w := #{$w};

end;
END

