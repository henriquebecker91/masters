#ifndef HBM_UKP5_HPP
#define HBM_UKP5_HPP

#include <sstream>  // For stringstream
#include <iostream>
#include <iomanip>  // For precision related routines

#include "ukp_common.hpp"

#ifdef HBM_PROFILE
  #include <chrono>
#endif

#if !defined(HBM_PROFILE) && defined(HBM_DUMP)
  #error The HBM_DUMP flag can only be used with the HBM_PROFILE flag
#endif

#if defined(HBM_PROFILE) && defined(HBM_DUMP)
  #include <boost/filesystem.hpp>
#endif

#ifdef HBM_PROFILE
  #ifndef HBM_PROFILE_PRECISION
    #define HBM_PROFILE_PRECISION 5
  #endif
#endif

namespace hbm {
  namespace hbm_ukp5_impl {
    using namespace std;

    template <typename W, typename P, typename I>
    struct ukp5_extra_info_t : extra_info_t {
      #ifdef HBM_CHECK_PERIODICITY
      /// @attention The last_y_value_outer_loop field exists only if
      /// HBM_CHECK_PERIODICITY is defined.

      /// The last capacity computed before detecting periodicity and stoping.
      W last_y_value_outer_loop;
      #endif // HBM_CHECK_PERIODICITY

      #ifdef HBM_PROFILE
      // Time of each phase
      double sort_time;        ///< Time used sorting items.
      double vector_alloc_time;///< Time used allocating vectors for DP.
      double linear_comp_time; ///< Time used by linear time preprocessing.
      double phase1_time;      ///< Time used by ukp5 phase1 (find optimal).
      double phase2_time;      ///< Time used by ukp5 phase2 (assemble solution).
      double total_time;       ///< Time used by all previous steps.
      // Some data about instance
      I n;   ///< Instance number of items.
      W c,   ///< Instance capacity.
      w_min, ///< Instance smallest item weight.
      w_max; ///< Instance biggest item weight.
      // Some data about structures manipulated by ukp5
      std::vector<P> g; ///< The vector of size c+w_max+1 and profit values.
      std::vector<I> d; ///< The vector of size c+w_max+1 and item index values.
      // Some statistics
      /// Number of items sized vector, with the quantity of each value i in dy.
      std::vector<W> qt_i_in_dy;
      /// Same as d, but without the positions skipped.
      std::vector<I> non_skipped_d;
      /// Last position of the d vector that wasn't zero or n (number of items).
      W last_dy_non_zero_non_n;
      /// Quantity of g positions that weren't skipped by ukp5.
      W qt_non_skipped_ys;
      /// Quantity of zeros in g.
      W qt_gy_zeros;
      /// How many times ukp5 phase 1 inner loop executed. Sum of non_skipped_d.
      W qt_inner_loop_executions;
      #endif // HBM_PROFILE

      virtual string gen_info(void) {
        /* TEMPORARY DISABLE UNTIL FIND WAY TO USE THIS AGAIN
        #if defined(HBM_PROFILE) && defined(HBM_DUMP)
        using namespace boost::filesystem;
        path my_path(spath);
        path filename = my_path.filename();
        const string  gd_path = "./g_dump_" + filename.native() + ".dat",
                      dd_path = "./d_dump_" + filename.native() + ".dat",
                      nsd_path = "./nsd_dump_" + filename.native() + ".dat",
                      dqt_path = "./dqt_dump_" + filename.native() + ".dat";
        dump(gd_path, "y\tgy", g);
        dump(dd_path, "y\tdy", d);
        //dump(nsd_path, "y\tdy", non_skipped_d);
        //dump(dqt_path, "i\tqt_in_d", qt_i_in_dy);
        #endif*/
        stringstream cout("");
        #ifdef HBM_CHECK_PERIODICITY
        cout << "last_y_value_outer_loop: " << last_y_value_outer_loop << endl;
        #endif
        #ifdef HBM_PROFILE
        cout << "c: " << c << endl;
        cout << "n: " << n << endl;
        cout << "w_min: " << w_min << endl;
        cout << "w_max: " << w_max << endl;
        cout << "last_dy_non_zero_non_n: " << last_dy_non_zero_non_n << endl;
        cout << "qt_non_skipped_ys: " << qt_non_skipped_ys << endl;
        cout << "qt_gy_zeros: " << qt_gy_zeros << endl;
        cout << "qt_inner_loop_executions: " << qt_inner_loop_executions << endl;
        cout << "qt_inner_loop_executions/qt_non_skipped_ys: " << ((long double) qt_inner_loop_executions)/((long double) qt_non_skipped_ys) << endl;
        cout << "qt_inner_loop_executions/c: " << ((long double) qt_inner_loop_executions)/((long double) c) << endl;
        cout << "(qt_inner_loop_executions/qt_non_skipped_ys)/n: " << ((long double) qt_inner_loop_executions)/((long double) qt_non_skipped_ys)/((long double)n)  << endl;
        cout << "(qt_inner_loop_executions/c)/n: " << ((long double) qt_inner_loop_executions)/((long double) c)/((long double)n) << endl;

        const double &stime = sort_time, &vtime = vector_alloc_time,
          &lctime = linear_comp_time, &p1time = phase1_time,
          &p2time = phase2_time, &ttime = total_time;

        streamsize old_precision = cout.precision(HBM_PROFILE_PRECISION);
        const int two_first_digits_and_period = 3;
        const int percent_size = two_first_digits_and_period+HBM_PROFILE_PRECISION;
        ios_base::fmtflags old_flags = cout.setf(std::ios::fixed, std:: ios::floatfield);
        char old_fill = cout.fill(' ');

        cout << "Sort time: " << stime << "s (";
        cout << setw(percent_size) << (stime/ttime)*100;
        cout << "%)" << endl;
        cout << "Vect time: " << vtime << "s (";
        cout << setw(percent_size) << (vtime/ttime)*100;
        cout << "%)" << endl;
        cout << "O(n) time: " << lctime << "s (";
        cout << setw(percent_size) << (lctime/ttime)*100;
        cout << "%)" << endl;
        cout << "pha1 time: " << p1time << "s (";
        cout << setw(percent_size) << (p1time/ttime)*100;
        cout << "%)" << endl;
        cout << "pha2 time: " << p2time << "s (";
        cout << setw(percent_size) << (p2time/ttime)*100;
        cout << "%)" << endl;
        cout << "Sum times: " << ttime << "s" << endl;

        cout.fill(old_fill);
        cout.setf(old_flags);
        cout.precision(old_precision);
        #endif

        return cout.str();
      }
    };

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
      } else {
        auto opts = get_opts(c, g, w_min);
        opt = opts.first;
        y_opt = opts.second;
      }
      #ifdef HBM_PROFILE
      shared_ptr< ukp5_extra_info_t<W, P, I> > info =
        dynamic_pointer_cast<ukp5_extra_info_t<W, P, I> >(sol.extra_info);
      info->last_y_value_outer_loop = last_y_where_nonbest_item_was_used+1;
      #endif

      #else //HBM_CHECK_PERIODICITY
      auto opts = get_opts(c, g, w_min);
      opt = opts.first;
      y_opt = opts.second;
      #endif //HBM_CHECK_PERIODICITY

      return;
    }

    #ifdef HBM_PROFILE
    using namespace std::chrono;

    template<typename W, typename P, typename I>
    void ukp5_gen_stats(W c, I n, W w_min, W w_max, const vector<P> &g, const vector<I> &d, ukp5_extra_info_t<W, P, I>* &info) {
      info->g = g;
      info->d = d;
      info->c = c;
      info->n = n;
      info->w_min = w_min;
      info->w_max = w_max;

      info->non_skipped_d.assign(c-w_min + 1, n-1);

      info->qt_i_in_dy.assign(n, 0);
      info->qt_i_in_dy[n-1] = w_min;

      info->qt_gy_zeros = w_min;
      P opt = 0;
      info->qt_non_skipped_ys = 0;
      info->qt_inner_loop_executions = 0;

      for (W y = w_min; y <= c-w_min; ++y) {
        ++(info->qt_i_in_dy[d[y]]);
        if (g[y] > opt) {
          ++(info->qt_non_skipped_ys);
          info->qt_inner_loop_executions += d[y];
          info->non_skipped_d[y] = d[y];
          opt = g[y];
        }
        if (g[y] == 0) ++(info->qt_gy_zeros);
        if (d[y] != 0 && d[y] != (n-1)) info->last_dy_non_zero_non_n = y;
      }
      return;
    }
    #endif //HBM_PROFILE
    
    template<typename W, typename P, typename I>
    void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, I sort_k_most_eff) {
      // This pointer stores all the extra info that the UKP5
      // generate when HBM_PROFILE is defined
      ukp5_extra_info_t<W, P, I>* ptr = new ukp5_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(ptr);
      sol.extra_info = shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      steady_clock::time_point all_ukp5_begin = steady_clock::now();
      steady_clock::time_point begin = steady_clock::now();
      #endif
      sort_by_eff(ukpi.items, sort_k_most_eff);
      #ifdef HBM_PROFILE
      ptr->sort_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();
      #endif

      #ifdef HBM_PROFILE
      begin = steady_clock::now();
      #endif
      W c = ukpi.c;
      I n = ukpi.items.size();
      auto minw_max = minmax_item_weight(ukpi.items);
      W w_min = minw_max.first, w_max = minw_max.second;
      #ifdef HBM_PROFILE
      ptr->linear_comp_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();
      #endif

      #ifdef HBM_PROFILE
      begin = steady_clock::now();
      #endif

      #ifdef HBM_PROFILE
      /* Use the solution fields instead of local variables to propagate
       * the array values. The arrays will be dumped to files, making possible
       * study them with R or other tool. */
      vector<P> &g = ptr->g;
      vector<I> &d = ptr->d;
      g.assign(c+1+(w_max-w_min), 0);
      d.assign(c+1+(w_max-w_min), n-1);
      #else
      vector<P> g(c+1+(w_max-w_min), 0);
      vector<I> d(c+1+(w_max-w_min), n-1);
      #endif

      #ifdef HBM_PROFILE
      ptr->vector_alloc_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();
      #endif

      #ifdef HBM_PROFILE
      begin = steady_clock::now();
      #endif
      ukp5_phase1(ukpi, g, d, sol, w_min, w_max);
      #ifdef HBM_PROFILE
      ptr->phase1_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();
      begin = steady_clock::now();
      #endif
      ukp5_phase2(ukpi.items, d, sol);
      #ifdef HBM_PROFILE
      ptr->phase2_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();
      ptr->total_time = duration_cast<duration<double>>(steady_clock::now() - all_ukp5_begin).count();
      #endif

      #ifdef HBM_PROFILE
      ukp5_gen_stats(c, n, w_min, w_max, g, d, ptr);
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

    template<typename W, typename P, typename I>
    void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, double sort_k_most_eff) {
      /* This procedure doesn't call itself. It calls the "I sort_k_most_eff"
       * overloaded version. */
      if (sort_k_most_eff < 0.0 || sort_k_most_eff > 1.0) {
        cerr << "WARNING: ukp5 (double sort_k_most_eff overloaded version): "
                "the value of sort_k_most_eff must be between 0 and 1. "
                "The sort_k_most_eff value was: " << sort_k_most_eff << ". "
                "Will execute ukp5 with sort_k_most_eff = 1."
        << endl;
        hbm::hbm_ukp5_impl::ukp5(ukpi, sol, static_cast<I>(ukpi.items.size()));
      } else {
        double dsize = static_cast<double>(ukpi.items.size());
        I k = static_cast<I>(dsize*sort_k_most_eff);
        hbm::hbm_ukp5_impl::ukp5(ukpi, sol, k);
      }
    }

    template<typename W, typename P, typename I>
    void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, char** argv) {
      /* This procedure doesn't call itself. It calls the "I sort_k_most_eff"
       * or the "double sort_k_most_eff" overloaded versions. */
      if (argc == 0) {
        hbm::hbm_ukp5_impl::ukp5(ukpi, sol, static_cast<I>(ukpi.items.size()));
      } else if (argc == 1) {
        double percent;
        from_string(argv[0], percent);
        hbm::hbm_ukp5_impl::ukp5(ukpi, sol, percent);
      } else if (argc > 1) {
        cout << "WARNING: more than one extra parameter to ukp5." << endl;
        cout << "usage: a.out data.ukp [k]" << endl;
        cout << "       k: A real number between 0 and 1. The items will be"
                " reordered before ukp5 is called, as std::partial_sort"
                " was called with the 'middle' argument being position k*"
                "<number of items>. If k is not provided we assume 1 (the"
                " entire vector will be ordered by efficiency). Use 0 if the"
                " items are already ordered by efficiency in the instance"
                " file, or you don't want to sort the items by efficiency."
                << endl;
        hbm::hbm_ukp5_impl::ukp5(ukpi, sol, static_cast<I>(ukpi.items.size()));
      }
    }
  }

  template<typename W, typename P, typename I>
  void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, I sort_k_most_eff) {
    hbm_ukp5_impl::ukp5(ukpi, sol, sort_k_most_eff);
  }

  template<typename W, typename P, typename I>
  void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, double sort_k_most_eff) {
    hbm_ukp5_impl::ukp5(ukpi, sol, sort_k_most_eff);
  }

  template<typename W, typename P, typename I>
  void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, char** argv) {
    hbm_ukp5_impl::ukp5(ukpi, sol, argc, argv);
  }
}

#endif //HBM_UKP5_HPP
