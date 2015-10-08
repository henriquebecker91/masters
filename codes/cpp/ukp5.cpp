#include "ukp5.hpp"

#if defined(CHECK_PERIODICITY) && (defined(INT_EFF) || defined(FP_EFF))
  #error CHECK PERIODICITY ONLY MAKE SENSE IF THE ORDER IS NOT AN APPROXIMATION
#endif

using namespace std;

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

/* This function reorders the ukpi.items vector, if you don't want this pass a
 * copy of the instance or pass it already ordered by non-decreasing efficiency
 * and true for the parameter already_sorted.
 */
void ukp5(ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted/* = false*/) {
  size_t n = ukpi.items.size();
  size_t c = ukpi.c;
  vector<item_t> &items(ukpi.items);
  if (!already_sorted) sort_by_efficiency(ukpi.items);

  auto minmax_w = minmax_item_weight(items);
  size_t min_w = minmax_w.first, max_w = minmax_w.second;

  vector<size_t> &g = sol.g;
  vector<size_t> &d = sol.d;
  size_t &opt = sol.opt;
  opt = 0;

  /* After the block bellow we can safely assume that there are at least two
   * items, and one of them is smaller than c*/
  switch (n) {
    case 0: return;
    case 1: opt = (c % items[0].w)* items[0].p; return;
  }
  if (c < min_w) return;

  g.assign(c+1+(max_w-min_w), 0);
  d.assign(c+1+(max_w-min_w), n-1);
  
  #ifdef CHECK_PERIODICITY
  size_t last_y_where_nonbest_item_was_used = 0;
  #endif

  /* this block is a copy-past of the loop bellow only for the best item */
  size_t wb = items[0].w;
  g[wb] = items[0].p;;
  d[wb] = 0;

  for (size_t i = 0; i < n; ++i) {
    size_t pi = items[i].p;
    size_t wi = items[i].w;
    if (g[wi] < pi) {
      g[wi] = pi;
      d[wi] = i;
      #ifdef CHECK_PERIODICITY
      if (wi > last_y_where_nonbest_item_was_used) {
        last_y_where_nonbest_item_was_used = wi;
      }
      #endif
    }
  }

  opt = 0;
  for (size_t y = min_w; y <= c-min_w; ++y) {
    if (g[y] <= opt) continue;
    #ifdef CHECK_PERIODICITY
    if (last_y_where_nonbest_item_was_used < y) break;
    #endif

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

    for (size_t ix = 1; ix <= dy; ++ix) {
      item_t it = items[ix];
      size_t pi = it.p;
      size_t wi = it.w;
      size_t ny = y + wi;
      size_t ogny = g[ny];
      size_t ngny = gy + pi;
      if (ogny < ngny) {
        g[ny] = ngny;
        d[ny] = ix;
        #ifdef CHECK_PERIODICITY
        if (ny > last_y_where_nonbest_item_was_used) last_y_where_nonbest_item_was_used = ny;
        #endif
      }
    } 
  }

  #ifdef CHECK_PERIODICITY
  if (last_y_where_nonbest_item_was_used < c-min_w) {
    size_t y_ = last_y_where_nonbest_item_was_used;
    while (d[y_] != 0) ++y_;

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
  #endif

  size_t y_opt;
  for (size_t y = c-min_w+1; y <= c; ++y) {
    if (opt < g[y]) {
      opt = g[y];
      y_opt = y;
    }
  }

  return;
}

