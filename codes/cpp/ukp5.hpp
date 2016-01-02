#ifndef HBM_UKP5_HPP
#define HBM_UKP5_HPP

#include <sstream>  // For stringstream
#include <fstream>  // ofstream for dump
#include <iostream>
#include <iomanip>  // For precision related routines
#include <regex>    // For command-line argument handmade verification
#include <chrono>
#include <boost/filesystem.hpp>

// Includes for command-line arguments parsing
// Not used yet. Will use in the future.
//#include <boost/program_options/option.hpp>
//#include <boost/program_options/options_description.hpp>
//#include <boost/program_options/parsers.hpp>

#include "ukp_common.hpp"
#include "dominance.hpp"

#ifndef HBM_PROFILE_PRECISION
  #define HBM_PROFILE_PRECISION 5
#endif

namespace hbm {
  /// Class with the UKP5 extra parameters.
  /// It's initialized with the default values. You only need
  /// to touch it if you want to tweak the UKP5 inner workings.
  template <typename I>
  struct ukp5_conf_t {
    bool apply_smdom{false};
    bool apply_smdom_before_sort{true};
    bool use_percent{true};
    union {
      double sort_percent{1.0};
      I sort_k_most_eff;
    };
    bool create_dumps{false};
    /// Dump paths, if empty they are set by ukp5 (argc/argv version).
    std::string gd_path, dd_path/*, nsd_path, dqt_path*/;
  };

  namespace hbm_ukp5_impl {
    using namespace std;
    using namespace std::regex_constants;
    using namespace std::chrono;
    //namespace po = boost::program_options;

    template <typename W, typename P, typename I>
    struct ukp5_extra_info_t : extra_info_t {
      #ifdef HBM_CHECK_PERIODICITY
      /// @attention The last_y_value_outer_loop field exists only if
      /// HBM_CHECK_PERIODICITY is defined.

      /// The last capacity computed before detecting periodicity and stoping.
      W last_y_value_outer_loop;
      #endif // HBM_CHECK_PERIODICITY

      // Time of each phase
      double sort_time;        ///< Time used sorting items.
      double smdom_time;       ///< Time used removing simple/multiple dominance.
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

      /// If true, will try to create dump files on gd_path and dd_path.
      bool create_dumps;
      string gd_path, ///< Path where the g vector will be dumped.
                  dd_path; ///< Path where the d vector will be dumped.

      template <typename W_OR_P>
      static void dump(const string &path,
                       const string &header,
                       const vector<W_OR_P> &v) {
        ofstream f(path, ofstream::out|ofstream::trunc);
        if (f.is_open())
        {
          f << header << endl;
          for (size_t y = 0; y < v.size(); ++y) {
            f << y << "\t" << v[y] << endl;
          }
        } else {
          cerr << "Couldn't open file: " << path << endl;
        }
      }

      virtual string gen_info(void) {
        ukp5_gen_stats();

        stringstream cout("");
        #ifdef HBM_CHECK_PERIODICITY
        cout << "last_y_value_outer_loop: " << last_y_value_outer_loop << endl;
        #endif
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

        const double &stime = sort_time, &dtime = smdom_time,
          &vtime = vector_alloc_time, &lctime = linear_comp_time,
          &p1time = phase1_time, &p2time = phase2_time, &ttime = total_time;

        streamsize old_precision = cout.precision(HBM_PROFILE_PRECISION);
        const int two_first_digits_and_period = 3;
        const int percent_size = two_first_digits_and_period+HBM_PROFILE_PRECISION;
        ios_base::fmtflags old_flags = cout.setf(std::ios::fixed, std:: ios::floatfield);
        char old_fill = cout.fill(' ');

        cout << "Sort time: " << stime << "s (";
        cout << setw(percent_size) << (stime/ttime)*100.0;
        cout << "%)" << endl;
        cout << "Dom  time: " << dtime << "s (";
        cout << setw(percent_size) << (dtime/ttime)*100.0;
        cout << "%)" << endl;
        cout << "Vect time: " << vtime << "s (";
        cout << setw(percent_size) << (vtime/ttime)*100.0;
        cout << "%)" << endl;
        cout << "O(n) time: " << lctime << "s (";
        cout << setw(percent_size) << (lctime/ttime)*100.0;
        cout << "%)" << endl;
        cout << "pha1 time: " << p1time << "s (";
        cout << setw(percent_size) << (p1time/ttime)*100.0;
        cout << "%)" << endl;
        cout << "pha2 time: " << p2time << "s (";
        cout << setw(percent_size) << (p2time/ttime)*100.0;
        cout << "%)" << endl;
        cout << "Sum times: " << ttime << "s" << endl;

        cout.fill(old_fill);
        cout.setf(old_flags);
        cout.precision(old_precision);

        if (create_dumps) {
          dump(gd_path, "y\tgy", g);
          dump(dd_path, "y\tdy", d);
          //dump(nsd_path, "y\tdy", non_skipped_d);
          //dump(dqt_path, "i\tqt_in_d", qt_i_in_dy);
        }

        return cout.str();
      }

      /// This method only works if the following variables are
      /// set first: n, c, w_min, w_max, g and d.
      void ukp5_gen_stats(void) {
        non_skipped_d.assign(c-w_min + 1, n-1);

        qt_i_in_dy.assign(n, 0);
        qt_i_in_dy[n-1] = w_min;

        qt_gy_zeros = w_min;
        P opt = 0;
        qt_non_skipped_ys = 0;
        qt_inner_loop_executions = 0;

        for (W y = w_min; y <= c-w_min; ++y) {
          ++(qt_i_in_dy[d[y]]);
          if (g[y] > opt) {
            ++(qt_non_skipped_ys);
            qt_inner_loop_executions += d[y];
            non_skipped_d[y] = d[y];
            opt = g[y];
          }
          if (g[y] == 0) ++(qt_gy_zeros);
          if (d[y] != 0 && d[y] != (n-1)) last_dy_non_zero_non_n = y;
        }
        return;
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

      shared_ptr< ukp5_extra_info_t<W, P, I> > info =
        dynamic_pointer_cast<ukp5_extra_info_t<W, P, I> >(sol.extra_info);
      info->last_y_value_outer_loop = last_y_where_nonbest_item_was_used+1;

      #else //HBM_CHECK_PERIODICITY
      auto opts = get_opts(c, g, w_min);
      opt = opts.first;
      y_opt = opts.second;
      #endif //HBM_CHECK_PERIODICITY

      return;
    }
    
    template<typename W, typename P, typename I>
    void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, const ukp5_conf_t<I> &conf) {
      steady_clock::time_point all_ukp5_begin = steady_clock::now();
      steady_clock::time_point begin; // Only declare

      // This pointer stores all the extra info that the UKP5 will print after.
      ukp5_extra_info_t<W, P, I>* ptr = new ukp5_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(ptr);
      sol.extra_info = shared_ptr<extra_info_t>(upcast_ptr);

      if (conf.create_dumps) {
        ptr->create_dumps = conf.create_dumps;
        ptr->gd_path = conf.gd_path;
        ptr->dd_path = conf.dd_path;
      }

      if (conf.apply_smdom && conf.apply_smdom_before_sort) { 
        begin = steady_clock::now();
        vector< item_t<W, P> > undominated;
        smdom(ukpi.items, undominated, 0, true);
        if (undominated.size() < ukpi.items.size()) {
          ukpi.items = undominated;
        }
        ptr->smdom_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();
      }

      // How many elements we will partially sort
      I k;
      if (conf.use_percent) {
        // If the conf says we sort a fraction of the list, we compute
        // how many elements make that fraction.
        if (conf.sort_percent < 0.0 || conf.sort_percent > 1.0) {
          cerr << "WARNING: " << __func__ << ": "
                  "the value of conf.sort_percent must be between 0 and 1. "
                  "The conf.sort_percent value was: " << conf.sort_percent <<
                  ". Will execute ukp5 with conf.sort_percent = 1."
          << endl;
          k =  static_cast<I>(ukpi.items.size());
        } else {
          double dsize = static_cast<double>(ukpi.items.size());
          k = static_cast<I>(dsize*conf.sort_percent);
        }
      } else {
        k = conf.sort_k_most_eff;
      }

      begin = steady_clock::now();
      sort_by_eff(ukpi.items, k);
      ptr->sort_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();

      if (conf.apply_smdom && !conf.apply_smdom_before_sort) { 
        begin = steady_clock::now();
        vector< item_t<W, P> > undominated;
        smdom(ukpi.items, undominated, k, true);
        if (undominated.size() < ukpi.items.size()) {
          ukpi.items = undominated;
        }
        ptr->smdom_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();
      }

      begin = steady_clock::now();
      W c = ukpi.c;
      I n = ukpi.items.size();
      auto min_max_w = minmax_item_weight(ukpi.items);
      W w_min = min_max_w.first, w_max = min_max_w.second;
      ptr->linear_comp_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();
      ptr->c = c;
      ptr->n = n;
      ptr->w_min = w_min;
      ptr->w_max = w_max;

      begin = steady_clock::now();
      /* Use the solution fields instead of local variables to propagate
       * the array values. The arrays will be dumped to files, making possible
       * study them with R or other tool. */
      vector<P> &g = ptr->g;
      vector<I> &d = ptr->d;
      g.assign(c+1+(w_max-w_min), 0);
      d.assign(c+1+(w_max-w_min), n-1);
      ptr->vector_alloc_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();

      begin = steady_clock::now();
      ukp5_phase1(ukpi, g, d, sol, w_min, w_max);
      ptr->phase1_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();

      begin = steady_clock::now();
      ukp5_phase2(ukpi.items, d, sol);
      ptr->phase2_time = duration_cast<duration<double>>(steady_clock::now() - begin).count();

      #ifdef HBM_CHECK_PERIODICITY
      // If we use our periodicity check the sol.used_items constructed by
      // ukp5_phase2 doesn't include the copies of the best item used to
      // fill all the extra_capacity. This solves it.
      // NOTE: this works even if the array isn't sorted
      W qt_best_item_inserted_by_per = (c - sol.y_opt)/ukpi.items[0].w;
      sol.opt += static_cast<P>(qt_best_item_inserted_by_per)*ukpi.items[0].p;
      sol.y_opt+= qt_best_item_inserted_by_per*ukpi.items[0].w;
      // Our periodicity check will always get a partial solution that have at
      // least one of the best item. And ukp5_phase2 will populate the
      // sol.used_items in the same order the sol.items are, this way we
      // have the guarantee that the first element of sol.used_items
      // exist and it's the best item.
      sol.used_items[0].qt += qt_best_item_inserted_by_per;
      #endif

      ptr->total_time = duration_cast<duration<double>>(steady_clock::now() - all_ukp5_begin).count();

      return;
    }

    void usage(const char *func_name) {
      cout << "Usage info from: " << func_name << "." << endl;
      cout << "usage: ./a.out data.ukp [k apply_dom apply_dom_before_sort"
              "create_dumps]"
              << endl;
      cout << "       k: A real number between 0 and 1. The items will be"
              " reordered before ukp5 is called, as std::partial_sort"
              " was called with the 'middle' argument being position k*"
              "<number of items>. If k is not provided we assume 1 (the"
              " entire vector will be ordered by efficiency). Use 0 if the"
              " items are already ordered by efficiency in the instance"
              " file, or you don't want to sort the items by efficiency."
              << endl;
      cout << "       apply_dom: A boolean value (0 or 1)."
              " If it's one, simple/multiple dominance is applied to"
              " the items before executing ukp5. The default value is"
              " zero (UKP5 already removes dominance on its own way,"
              " if the items are sorted by non-increasing efficiency)."
              << endl;
      cout << "       apply_dom_before_sort: A boolean value (0 or 1)."
              " Apply dominance before or after sorting. Dominance after"
              " sort is useful only if you will sort all the items ("
              " k = 1), because in this case a faster version of dominance"
              " can be used. But: if there was few dominated items, would"
              " be better to not use dominance, and if there's a lot, the"
              " sorting cost could be reduced by applying dominance first."
              << endl;
      cout << "       create_dumps: A boolean value (0 or 1). If true will"
              " set the solution_t.extra_info field to save the value of"
              " the g and d vectors on a format that can be easily read"
              " by R. Can create very big files."
              << endl;
      cout << "NOTE: the arguments are taken by position, if you want to"
              " specify the second argument value, you need to specify the"
              " first one too." << endl;
    }

    template<typename W, typename P, typename I>
    void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
      // The empty constructor already has default values for the fields
      // below we only change what was explicited by command-line
      ukp5_conf_t<I> conf;

      const regex float_0_1("(0|1|0\\.[0-9]*[1-9])", nosubs),
                  bool_0_1("[01]", nosubs);

      switch (argc) {
        case 6:
        if (!regex_match(argv[5], bool_0_1)) {
          cerr << "Parameter create_dumps isn't valid. "
                  "It's '" << argv[5] << "' when the valid values "
                  "are 0 or 1."
          << endl;
          usage(__func__);
        } else {
          from_string(argv[5], conf.create_dumps);
        }
        case 5:
        if (!regex_match(argv[4], bool_0_1)) {
          cerr << "Parameter apply_dom_before_sort isn't valid. "
                  "It's '" << argv[4] << "' when the valid values "
                  "are 0 or 1."
          << endl;
          usage(__func__);
        } else {
          from_string(argv[4], conf.apply_smdom_before_sort);
        }
        case 4:
        if (!regex_match(argv[3], bool_0_1)) {
          cerr << "Parameter apply_dom isn't valid. "
                  "It's '" << argv[3] << "' when the valid values "
                  "are 0 or 1."
          << endl;
          usage(__func__);
        } else {
          from_string(argv[3], conf.apply_smdom);
        }
        case 3:
        if (!regex_match(argv[2], float_0_1)) {
          cerr << "Parameter k isn't valid. "
                  "It's '" << argv[2] << "' when the valid values "
                  "are real numbers between 0.0 and 1.0."
          << endl;
          usage(__func__);
        } else {
          from_string(argv[2], conf.sort_percent);
        }
        case 2:
        // exists only to avoid giving waning if the method was called
        // without parameters
        break;

        case 1:
        case 0:
        cerr << "WARNING: " << __func__ << ": argc < 2?! Executing method"
                " with default parameter values."
        << endl;
        break;

        default:
        cout << "WARNING: " << __func__ << ": more than " << argc << " "
              "parameters to ukp5. Will execute with the default values."
        << endl;
        usage(__func__);
        break;
      }

      if (conf.create_dumps) {
        using namespace boost::filesystem;
        path my_path = path(string(argv[1]));
        path fname = my_path.filename();
        if (conf.gd_path.empty()) {
          conf.gd_path = "./g_dump_" + fname.native() + ".dat";
        }
        if (conf.dd_path.empty()) {
          conf.dd_path = "./d_dump_" + fname.native() + ".dat";
        }
        //conf.nsd_path = "./nsd_dump_" + filename.native() + ".dat";
        //conf.dqt_path = "./dqt_dump_" + filename.native() + ".dat";
      }

      hbm_ukp5_impl::ukp5(ukpi, sol, conf);
    }
  }

  template<typename W, typename P, typename I>
  void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, const ukp5_conf_t<I> &conf) {
    hbm_ukp5_impl::ukp5(ukpi, sol, conf);
  }

  template<typename W, typename P, typename I>
  void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_ukp5_impl::ukp5(ukpi, sol, argc, argv);
  }
}

#endif //HBM_UKP5_HPP
