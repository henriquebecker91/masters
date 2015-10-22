def garfinkel_instance
  {
    n: 4,
    c: 25,
    items: [
      { p: 11, w: 6},
      { p: 7, w: 4},
      { p: 5, w: 3},
      { p: 1, w: 1}
    ]
  }
end

class String
  # Strip leading whitespace from each line that is the same as the 
  # amount of whitespace on the first line of the string.
  # Leaves _additional_ indentation on later lines intact.
  def unindent
    gsub /^#{self[/\A\s*/]}/, ''
  end
end

def sort_items_by_profitability!(items)
  items.sort_by! do | i |
    wi = i[:w]
    [Rational(wi,i[:p]), wi]
  end
end

def print_values(y, g, d, c, w_max)
    ys = gs = ds = ''
    lower = [(y/w_max)*w_max, c-2*w_max].min
    upper = [lower+2*w_max, c].min
    (lower..upper).each do | y2 |
      if y == y2 then
        ys += "\\textbf{#{y2}} &"
        gs += "\\textbf{#{g[y2]}} &"
        ds += "\\textbf{#{d[y2]+1}} & "
      else
        ys += "#{y2} &"
        gs += "#{g[y2]} &"
        ds += "#{d[y2]+1} & "
      end
    end

    puts <<-END.unindent
      \\begin{table}[h]
      \\centering
      \\caption{Memory at end of the outer loop (line \\ref{main_ext_loop_end}) where \(y = #{y}\)}
      \\label{mem_y_#{y}}
      \\begin{tabular}{l|#{'r'*((2*w_max)+1)}r}
      y: & #{ys}\\dots\\\\
      \\hline
      g: & #{gs}\\dots\\\\
      d: & #{ds}\\dots\\\\
      \\end{tabular}
      \\end{table}
    END
    puts
end

def ukp5(ukpi)
  n = ukpi[:n]
  c = ukpi[:c]
  items = ukpi[:items].clone
  sort_items_by_profitability!(items)
  w_max = items.max_by { | i | i[:w] }[:w]
  w_min = items.min_by { | i | i[:w] }[:w]

  g = Array.new(c+(w_max-w_min)+1, 0)
  d = Array.new(c+(w_max-w_min)+1, n-1)

  print_values(0, g, d, c, w_max)

  items.each_with_index do | it, ix |
    wi = it[:w]
    pi = it[:p]
    if g[wi] < pi then
      g[wi] = pi
      d[wi] = ix
    end
  end

  print_values(0, g, d, c, w_max)

  opt = 0
  (w_min..(c-w_min)).each do | y |
    next if g[y] <= opt

    opt = gy = g[y]
    dy = d[y]

    (0..dy).each do | ix |
      it = items[ix]
      pi = it[:p]
      wi = it[:w]
      ny = y + wi
      ogny = g[ny]
      ngny = gy + pi
      if ogny < ngny then
        g[ny] = ngny
        d[ny] = ix
      end
    end

    print_values(y, g, d, c, w_max)
  end

  y_opt = nil
  ((c-w_min+1)..c).each do | y |
    if g[y] > opt then
      opt = g[y]
      y_opt = y
    end
  end
  
  puts "opt: #{opt}"
  puts "y_opt: #{y_opt}"
end

ukp5(garfinkel_instance)

