#include <algorithm>
#include <boost/rational.hpp>

#include "ukp.hpp"

using namespace std;
using namespace boost;

pair<size_t,size_t> minmax_item_weight(vector<item_t> &items) {
  size_t min, max;
  min = max = items[0].w;
  for (auto it = items.begin()+1; it != items.end(); ++it) {
    size_t x = (*it).w;
    if (x < min) min = x;
    else if (x > max) max = x;
  }
  return make_pair(min,max);
}

/* Profitability nonascending first weight nondescending after */
bool item_order(const item_t& i, const item_t& j) {
  rational<size_t> ri(i.w, i.p);
  rational<size_t> rj(j.w, j.p);

  return ri == rj ? i.w < j.w : ri < rj;
}

void sort_items_by_profitability(vector<item_t> &items) {
  sort(items.begin(), items.end(), item_order);

  return;
}

size_t get_opt_y(size_t c, const vector<item_t> &items, const vector<size_t> &g, const vector<size_t> &d, size_t w_min) {
  size_t ix = c - w_min;

  size_t opt = g[ix];
  size_t opt_y = ix;

  for (size_t y = ix+1; y <= c; ++y) {
    if (g[y] > opt) {
      opt = g[y];
      opt_y = y;
    }
  }

  return opt_y;
}

/* This function reorders the ukpi.items vector, if you don't want this pass a copy of
 * the instance or pass it already ordered and true for the parameter already_sorted.
 */
void ukp5(ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted/* = false*/) {
  size_t n = ukpi.items.size();
  size_t c = ukpi.c;
  vector<item_t> &items(ukpi.items);
  if (!already_sorted) sort_items_by_profitability(ukpi.items);

  auto minmax_w = minmax_item_weight(items);
  size_t min_w = minmax_w.first, max_w = minmax_w.second;

  vector<size_t> &g = sol.g;
  vector<size_t> &d = sol.d;
  size_t &opt = sol.opt;
  opt = 0;

  g.assign(c+max_w+1, 0);
  d.assign(c+max_w+1, n-1);
  
  size_t last_y_where_nonbest_item_was_used = 0;

  for (size_t i = 0; i < n; ++i) {
    size_t pi = items[i].p;
    size_t wi = items[i].w;
    if (g[wi] < pi) {
      g[wi] = pi;
      d[wi] = i;
      if (wi > last_y_where_nonbest_item_was_used && i != 0) {
        last_y_where_nonbest_item_was_used = wi;
      }
    }
  }

  for (size_t y = 0; y < c; ++y) {
    if (g[y] == 0 || g[y] < opt) continue;
    if (last_y_where_nonbest_item_was_used < y) break;

    size_t gy, dy;
    opt = gy = g[y];
    dy = d[y];

    /* this block is a copy-past of the loop bellow only for the best item */
    item_t bi = items[0];
    size_t pb = bi.p;
    size_t wb = bi.w;
    size_t next_y = y + wb;
    size_t old_gny = g[next_y];
    size_t new_gny = gy + pb;
    if (old_gny < new_gny) {
      g[next_y] = new_gny;
      d[next_y] = 0;
    }

    for (size_t ix = 0; ix <= dy; ++ix) {
      item_t it = items[ix];
      size_t pi = it.p;
      size_t wi = it.w;
      size_t ny = y + wi;
      size_t ogny = g[ny];
      size_t ngny = gy + pi;
      if (ogny < ngny) {
        g[ny] = ngny;
        d[ny] = ix;
        if (ny > last_y_where_nonbest_item_was_used) last_y_where_nonbest_item_was_used = ny;
      }
    } 
  }

  if (last_y_where_nonbest_item_was_used < c-1) {
    size_t y_ = last_y_where_nonbest_item_was_used;
    while (d[y_] != 0) y_ += 1;

    size_t extra_capacity = c - y_;
    size_t c1, a1;
    c1 = items[0].p;
    a1 = items[0].w;

    size_t qt_best_item_used = extra_capacity / a1;

    size_t profit_generated_by_best_item = qt_best_item_used*c1;
    size_t space_used_by_best_item = qt_best_item_used*a1;

    size_t opt_y = get_opt_y(c-space_used_by_best_item, items, g, d, min_w);
    g[c] = g[opt_y] + profit_generated_by_best_item;
  }

  if (opt < g[c]) opt = g[c];

  return;
}

void read_sukp_instance(istream &in, ukp_instance_t &ukpi) {
  size_t n;
  in >> n;
  in >> ukpi.c;
  ukpi.items.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    item_t tmp;
    in >> tmp.w;
    in >> tmp.p;
    ukpi.items.push_back(tmp);
  }

  return;
}

void write_sukp_instance(ostream &out, ukp_instance_t &ukpi) {
  size_t n = ukpi.items.size();
  out << n << endl;
  out << ukpi.c << endl;

  for (size_t i = 0; i < n; ++i) {
    item_t tmp = ukpi.items[i];
                out << tmp.w << "\t" << tmp.p << endl;
  }

  return;
}
