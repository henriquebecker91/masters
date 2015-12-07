#ifndef HBM_UKP5_HPP
#define HBM_UKP5_HPP

#include "ukp_common.hpp"

#ifdef HBM_PROFILE
  #include <chrono>
#endif

namespace hbm {
  namespace hbm_ukp5_impl {
    using namespace std;

    #ifdef HBM_PROFILE
    using namespace std::chrono;

    template<typename W, typename P, typename I>
    void ukp5_gen_stats(W c, I n, W w_min, W w_max, const vector<P> &g, const vector<I> &d, solution_t<W, P, I> &sol) {
      sol.g = g;
      sol.d = d;
      sol.c = c;
      sol.n = n;
      sol.w_min = w_min;
      sol.w_max = w_max;

      sol.non_skipped_d.assign(c-w_min + 1, n-1);

      sol.qt_i_in_dy.assign(n, 0);
      sol.qt_i_in_dy[n-1] = w_min;

      sol.qt_gy_zeros = w_min;
      P opt = 0;
      sol.qt_non_skipped_ys = 0;
      sol.qt_inner_loop_executions = 0;

      for (W y = w_min; y <= c-w_min; ++y) {
        ++sol.qt_i_in_dy[d[y]];
        if (g[y] > opt) {
          ++sol.qt_non_skipped_ys;
          sol.qt_inner_loop_executions += d[y];
          sol.non_skipped_d[y] = d[y];
          opt = g[y];
        }
        if (g[y] == 0) ++sol.qt_gy_zeros;
        if (d[y] != 0 && d[y] != (n-1)) sol.last_dy_non_zero_non_n = y;
      }
      return;
    }
    #endif //HBM_PROFILE

    template<typename W, typename P>
    pair<W, W> minmax_item_weight(vector< item_t<W, P> > &items) {
      W min, max;
      min = max = items[0].w;
      for (auto it = items.begin()+1; it != items.end(); ++it) {
        W x = it->w;
        if (x < min) min = x;
        else if (x > max) max = x;
      }
      return make_pair(min, max);
    }

    template<typename W, typename P>
    pair<P, W> get_opts(W c, const vector<P> &g, W w_min) {
      P opt = 0;
      /* Dont need to be initialized, but initializing to stop compiler
       * warning messages */
      W y_opt = 0;

      for (W y = (c-w_min)+1; y <= c; ++y) {
        if (opt < g[y]) {
          opt = g[y];
          y_opt = y;
        }
      }

      return make_pair(opt, y_opt);
    }

    template<typename W, typename P, typename I>
    void ukp5_phase2(const vector< item_t<W, P> > &items, const vector<I> &d, solution_t<W, P, I> &sol) {
      I n = items.size();
      vector<I> qts_its(n, 0);

      W y_opt = sol.y_opt;
      I dy_opt;
      while (y_opt != 0) {
        dy_opt = d[y_opt];
        y_opt -= items[dy_opt].w;
        ++qts_its[dy_opt];
      }

      for (I i = 0; i < n; ++i) {
        if (qts_its[i] > 0) {
          sol.used_items.emplace_back(items[i], qts_its[i], i);
        }
      }
      sol.used_items.shrink_to_fit();

      return;
    }

    template<typename W, typename P, typename I>
    void ukp5_phase1(const instance_t<W, P> &ukpi, vector<P> &g, vector<I> &d, solution_t<W, P, I> &sol, W w_min, W w_max) {
      const W &c = ukpi.c;
      const I &n = ukpi.items.size();
      const vector< item_t<W, P> > &items = ukpi.items;

      W &y_opt = sol.y_opt;
      P &opt = sol.opt;

      opt = 0;

      #ifdef HBM_CHECK_PERIODICITY
      W last_y_where_nonbest_item_was_used = 0;
      #endif

      /* this block is a copy-past of the loop bellow only for the best item
       * its utility is to simplify the code when HBM_CHECK_PERIODICITY is defined
       */
      W wb = items[0].w;
      g[wb] = items[0].p;;
      d[wb] = 0;

      for (W i = 0; i < n; ++i) {
        P pi = items[i].p;
        W wi = items[i].w;
        if (g[wi] < pi) {
          g[wi] = pi;
          d[wi] = i;
          #ifdef HBM_CHECK_PERIODICITY
          if (wi > last_y_where_nonbest_item_was_used) {
            last_y_where_nonbest_item_was_used = wi;
          }
          #endif
        }
      }

      opt = 0;
      for (W y = w_min; y <= c-w_min; ++y) {
        if (g[y] <= opt) continue;
        #ifdef HBM_CHECK_PERIODICITY
        if (last_y_where_nonbest_item_was_used < y) break;
        #endif

        P gy;
        I dy;
        opt = gy = g[y];
        dy = d[y];

        /* this block is a copy-past of the loop bellow only for the best item
         * its utility is to simplify the code when HBM_CHECK_PERIODICITY is defined
         */
        item_t<W, P> bi = items[0];
        W wb = bi.w;
        P pb = bi.p;
        W next_y = y + wb;
        P old_gny = g[next_y];
        P new_gny = gy + pb;
        if (old_gny < new_gny) {
          g[next_y] = new_gny;
          d[next_y] = 0;
        }

        for (I ix = 1; ix <= dy; ++ix) {
          item_t<W, P> it = items[ix];
          I wi = it.w;
          P pi = it.p;
          W ny = y + wi;
          P ogny = g[ny];
          P ngny = gy + pi;
          if (ogny < ngny) {
            g[ny] = ngny;
            d[ny] = ix;
            #ifdef HBM_CHECK_PERIODICITY
            if (ny > last_y_where_nonbest_item_was_used) last_y_where_nonbest_item_was_used = ny;
            #endif
          }
        } 
      }

      #ifdef HBM_CHECK_PERIODICITY
      if (last_y_where_nonbest_item_was_used < c-w_min) {
        W y_ = last_y_where_nonbest_item_was_used;
        while (d[y_] != 0) ++y_;

        W a1 = items[0].w;
        W extra_capacity = c - y_;

        W space_used_by_best_item = extra_capacity - (extra_capacity % a1);

        auto opts = get_opts(c-space_used_by_best_item, g, w_max);
        opt = opts.first;
        y_opt = opts.second;
        #ifdef HBM_PROFILE
        sol.last_y_value_outer_loop = last_y_where_nonbest_item_was_used+1;
        #endif
      } else {
        auto opts = get_opts(c, g, w_min);
        opt = opts.first;
        y_opt = opts.second;
        #ifdef HBM_PROFILE
        sol.last_y_value_outer_loop = c-w_min;
        #endif
      }
      #else
      auto opts = get_opts(c, g, w_min);
      opt = opts.first;
      y_opt = opts.second;
      #endif //HBM_CHECK_PERIODICITY

      return;
    }

    template<typename W, typename P, typename I>
    void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      #ifdef HBM_PROFILE
      steady_clock::time_point all_ukp5_begin = steady_clock::now();
      steady_clock::time_point begin = steady_clock::now();
      #endif
      if (!already_sorted) sort_by_eff(ukpi.items);
      #ifdef HBM_PROFILE
      sol.sort_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();
      #endif

      #ifdef HBM_PROFILE
      begin = steady_clock::now();
      #endif
      W c = ukpi.c;
      I n = ukpi.items.size();
      auto minw_max = minmax_item_weight(ukpi.items);
      W w_min = minw_max.first, w_max = minw_max.second;
      #ifdef HBM_PROFILE
      sol.linear_comp_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();
      #endif

      #ifdef HBM_PROFILE
      begin = steady_clock::now();
      #endif

      #ifdef HBM_PROFILE
      /* Use the solution fields instead of local variables to propagate
       * the array values. The arrays will be dumped to files, making possible
       * study them with R or other tool. */
      vector<P> &g = sol.g;
      vector<I> &d = sol.d;
      g.assign(c+1+(w_max-w_min), 0);
      d.assign(c+1+(w_max-w_min), n-1);
      #else
      vector<P> g(c+1+(w_max-w_min), 0);
      vector<I> d(c+1+(w_max-w_min), n-1);
      #endif

      #ifdef HBM_PROFILE
      sol.vector_alloc_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();
      #endif

      #ifdef HBM_PROFILE
      begin = steady_clock::now();
      #endif
      ukp5_phase1(ukpi, g, d, sol, w_min, w_max);
      #ifdef HBM_PROFILE
      sol.phase1_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();
      begin = steady_clock::now();
      #endif
      ukp5_phase2(ukpi.items, d, sol);
      #ifdef HBM_PROFILE
      sol.phase2_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();
      sol.total_time = duration_cast<duration<double>>(steady_clock::now() - all_ukp5_begin).count();
      #endif

      #ifdef HBM_PROFILE
      ukp5_gen_stats(c, n, w_min, w_max, g, d, sol);
      #endif

      #ifdef HBM_CHECK_PERIODICITY
      /* If we use our periodicity check the sol.used_items constructed by
       * ukp5_phase2 doesn't include the copies of the best item used to
       * fill all the extra_capacity. This solves it. */
      W qt_best_item_inserted_by_per = (c - sol.y_opt)/ukpi.items[0].w;
      sol.opt += static_cast<P>(qt_best_item_inserted_by_per)*ukpi.items[0].p;
      sol.y_opt+= qt_best_item_inserted_by_per*ukpi.items[0].w;
      /* Our periodicity check will always get a partial solution that have at
       * least one of the best item. And ukp5_phase2 will populate the
       * sol.used_items in the same order the sol.items are, this way we
       * have the guarantee that the first element of sol.used_items
       * exist and it's the best item. */
      sol.used_items[0].qt += qt_best_item_inserted_by_per;
      #endif

      return;
    }
  }

  /* This function reorders the ukpi.items vector, if you don't want this pass a
   * copy of the instance or pass it already ordered by non-decreasing eff
   * and true for the parameter already_sorted.
   */
  template<typename W, typename P, typename I>
  void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted = false) {
    hbm_ukp5_impl::ukp5(ukpi, sol, already_sorted);
  }
}

#endif //HBM_UKP5_HPP
