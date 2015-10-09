#include "ukp5.hpp"

#if defined(CHECK_PERIODICITY) && (defined(INT_EFF) || defined(FP_EFF))
  #error CHECK PERIODICITY ONLY MAKE SENSE IF THE ORDER IS NOT AN APPROXIMATION
#endif

using namespace std;

struct shared_data_t {
  vector<size_t> d;
};

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

pair<size_t, size_t> get_opts(size_t c, const vector<item_t> &items, const vector<size_t> &g, const vector<size_t> &d, size_t w_min) {
  size_t opt = 0;
  size_t y_opt;

  for (size_t y = (c-w_min)+1; y <= c; ++y) {
    if (opt < g[y]) {
      opt = g[y];
      y_opt = y;
    }
  }

  return make_pair(opt, y_opt);
}

void ukp5_phase2(const vector<item_t> &items, const shared_data_t &sd, ukp_solution_t &sol) {
  const vector<size_t> &d = sd.d;
  size_t n = items.size();
  vector<size_t> qts_its(n, 0);

  size_t y_opt = sol.y_opt;
  size_t dy_opt;
  while (y_opt != 0) {
    dy_opt = d[y_opt];
    y_opt -= items[dy_opt].w;
    ++qts_its[dy_opt];
  }

  for (size_t i = 0; i < n; ++i) {
    if (qts_its[i] > 0) {
      sol.used_items.emplace_back(items[i], qts_its[i]);
    }
  }
  sol.used_items.shrink_to_fit();

  return;
}

void ukp5_phase1(ukp_instance_t &ukpi, shared_data_t &sd, ukp_solution_t &sol, bool already_sorted/* = false*/) {
  size_t n = ukpi.items.size();
  size_t c = ukpi.c;
  vector<item_t> &items = ukpi.items;
  if (!already_sorted) sort_by_efficiency(ukpi.items);

  auto minmax_w = minmax_item_weight(items);
  size_t min_w = minmax_w.first, max_w = minmax_w.second;

  vector<size_t> g(c+1+(max_w-min_w), 0);
  vector<size_t> &d = sd.d;
  d.assign(c+1+(max_w-min_w), n-1);
  size_t &y_opt = sol.y_opt;
  size_t &opt = sol.opt;
  opt = 0;
  
  #ifdef CHECK_PERIODICITY
  size_t last_y_where_nonbest_item_was_used = 0;
  /* We could pre-compute the maximum weight at each point
   * to accelerate the periodicity check, but this would make
   * the looser, and the periodicity check doesn't seem to
   * consume much of the algorithm time
  vector<size_t> max_ws;
  max_ws.reserve(n);
  max_ws.push_back(items[0].w);
  for (size_t i = 1; i < n; ++i) {
    max_ws.push_back(max(max_ws[i-1], items[i].w));
  }*/
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
    y_opt = y;

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

    size_t c1, a1;
    c1 = items[0].p;
    a1 = items[0].w;
    size_t extra_capacity = c - y_;

    size_t qt_best_item_used = extra_capacity / a1;
    size_t profit_generated_by_best_item = qt_best_item_used*c1;
    size_t space_used_by_best_item = qt_best_item_used*a1;

    auto opts = get_opts(c-space_used_by_best_item, items, g, d, max_w);
    opt = opts.first + profit_generated_by_best_item;
    y_opt = opts.second + space_used_by_best_item;
  } else {
    auto opts = get_opts(c, items, g, d, min_w);
    opt = opts.first;
    y_opt = opts.second;
  }
  #else
  auto opts = get_opts(c, items, g, d, min_w);
  opt = opts.first;
  y_opt = opts.second;
  #endif

  return;
}

/* This function reorders the ukpi.items vector, if you don't want this pass a
 * copy of the instance or pass it already ordered by non-decreasing efficiency
 * and true for the parameter already_sorted.
 */
void ukp5(ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted/* = false*/) {
  shared_data_t sd;
  ukp5_phase1(ukpi, sd, sol, already_sorted);
  ukp5_phase2(ukpi.items, sd, sol);
}

