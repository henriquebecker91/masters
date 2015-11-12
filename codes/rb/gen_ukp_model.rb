#!/usr/bin/ruby

# putting the variables on the global scope
$n = nil
$c = nil
$p = ''
$w = ''

File.open(ARGV[0], 'r') do | f |
  bs = "[[:blank:]]*"
  bp = "[[:blank:]]+"
  nb = "([1-9][0-9]*)"
  comm_s = "[[:blank:]]*(#.*)?"
  mline_s = bs + "[mn]:" + bs + nb + comm_s
  cline_s = bs + "c:" + bs + nb + comm_s
  begin_data_s = bs + "begin" + bp + "data" + comm_s
  item_s = bs + nb + bp + nb + comm_s
  end_data_s = bs + "end" + bp + "data" + comm_s
  comm = /\A#{comm_s}\Z/
  mline = /\A#{mline_s}\Z/
  cline = /\A#{cline_s}\Z/
  begin_data = /\A#{begin_data_s}\Z/
  item = /\A#{item_s}\Z/
  end_data = /\A#{end_data_s}\Z/

  str = f.gets
  while comm.match(str) do str = f.gets end

  $n = mline.match(str)
  if $n
    $n = $n[1].to_i
  else
    fail "'n/m: <number>' wasn't found after comments"
  end

  str = f.gets
  while comm.match(str) do str = f.gets end

  $c = cline.match(str)
  if $c
    $c = $c[1].to_i
  else
    fail "'c: <number>' wasn't found after 'n/m: <number>'"
  end

  str = f.gets
  while comm.match(str) do str = f.gets end
  if !begin_data.match(str)
    fail "'begin data' wasn't found after 'c: <number>'"
  end

  i = 1
  while (i <= $n && str = f.gets)
    it = nil
    if it = item.match(str)
      $w << "#{i} #{it[1]}, "
      $p << "#{i} #{it[2]}, "
      i += 1
    else
      if end_data.match(str)
        warn "WARNING: 'end data' was found before #{$n} items were read, only #{i-1} items read"
        break
      elsif comm.match(str)
        warn  'WARNING: Comment line inside data section'
      else
        fail "ERROR: Strange line found inside data section: '#{str}'"
      end
    end
  end

  str = f.gets
  while comm.match(str) do str = f.gets end
  if !end_data.match(str) 
    warn "WARNING: first noncomment line after the last item read isn't 'end data':\n#{str}"
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
  # Ridiculous but necessary workaround: https://en.wikibooks.org/wiki/GLPK/GMPL_Workarounds
  for {{0}:  x[i] > 0}{ # IF condition THEN
    printf 'Obj nº %d p %d w %d: %d\\n', i, p[i], w[i], x[i];
  }                     # ENDIF
}

printf "Max Profit: %d\\n", sum {i in 1..n} p[i]*x[i];

data;

param n := #{$n};

param c := #{$c};

param p := #{$p};

param w := #{$w};

end;
END

