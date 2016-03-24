#ifndef HBM_UKP5_HPP
#define HBM_UKP5_HPP

#include <sstream>  // For stringstream
#include <fstream>  // ofstream for dump
#include <iostream>
#include <iomanip>  // For precision related routines
#include <regex>  // For command-line argument handmade verification
#include <algorithm>  // For find
#include <chrono>     
#include <boost/filesystem.hpp>

// Includes for command-line arguments parsing
// Not used yet. Will use in the future.
//#include <boost/program_options/option.hpp>
//#include <boost/program_options/options_description.hpp>
//#include <boost/program_options/parsers.hpp>

#include "ukp_common.hpp"
#include "dominance.hpp"
#include "periodicity.hpp"

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
    /// If true, will try to create dump files on the paths.
    /// The ukp5 (argc/argv version) function set this paths,
    /// if they are empty.
    bool create_dumps{false};
    
    std::string
      gd_path, ///< Path where the g vector will be dumped.
      dd_path, ///< Path where the d vector will be dumped.
      nsd_path, ///< Path where non_skipped_d will be dumped.
      dqt_path; ///< Path where qt_i_in_dy will be dumped.

    bool use_y_star_per{true};

    /// Write human-readable object representation to a stream.
    ///
    /// We take an stream as a parameter to allow people changing
    /// the precision the numbers should be shown.
    ///
    /// @param out An ostream where the object representation will be
    ///   outputed.
    void print(std::ostream &out = std::cout) const {
      // TODO: create a to_string overloaded funtion for this object
      // and use this function with an stringstream as implementation.
      #define HBM_PRINT_VAR(var) out << #var ": " << var << std::endl
      HBM_PRINT_VAR(apply_smdom);
      if (apply_smdom) HBM_PRINT_VAR(apply_smdom_before_sort);

      HBM_PRINT_VAR(use_percent);
      if (use_percent) {
        HBM_PRINT_VAR(sort_percent);
      } else {
        HBM_PRINT_VAR(sort_k_most_eff);
      }

      HBM_PRINT_VAR(create_dumps);
      if (create_dumps) {
        HBM_PRINT_VAR(gd_path);
        HBM_PRINT_VAR(dd_path);
        HBM_PRINT_VAR(nsd_path);
        HBM_PRINT_VAR(dqt_path);
      }

      HBM_PRINT_VAR(use_y_star_per);
    }
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

      /// The UKP5 used config.
      ukp5_conf_t<I> conf;

      // Time of each phase
      double sort_time;        ///< Time used sorting items.
      double smdom_time;       ///< Time used removing simple/multiple dominance.
      double vector_alloc_time;///< Time used allocating vectors for DP.
      double linear_comp_time; ///< Time used by linear time preprocessing.
      double phase1_time;      ///< Time used by ukp5 phase1 (find optimal).
      double phase2_time;      ///< Time used by ukp5 phase2 (assemble solution).
      double total_time;       ///< Time used by all previous steps.
      // Some data about instance
      I n;      ///< Instance number of items.
      W c,      ///< Instance capacity.
      y_bound,  ///< Capacity used if use_y_star_per is true.
      w_min,    ///< Instance smallest item weight.
      w_max;    ///< Instance biggest item weight.
      // Some data about structures manipulated by ukp5
      std::vector<P> g; ///< The vector with profit values.
      std::vector<I> d; ///< The vector with item index values.
      // Some statistics
      /// Number of items sized vector, with the quantity of each value i in dy.
      std::vector<W> qt_i_in_dy;
      /// Same as d, but without the positions skipped.
      std::vector<I> non_skipped_d;
      /// Last position of the d vector that wasn't zero or n (item indexes).
      W last_dy_non_zero_non_n;
      /// Quantity of g positions that weren't skipped by ukp5.
      W qt_non_skipped_ys;
      /// Quantity of zeros in g.
      W qt_gy_zeros;
      /// How many times ukp5 phase 1 inner loop executed. Sum of non_skipped_d.
      W qt_inner_loop_executions;

      /// If an file existed in path, its content is dicarded.
      /// The value of header is written in the file at path, then
      /// two columns are written in the file. The first column is
      /// composed from the vector indexes (will vary between 0 and
      /// v.size()-1). The second column has the vector values.
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
          cerr << __func__ << ": Couldn't open file: " << path << endl;
        }
      }

      /// Generates a good number of stats about the UKP5 execution.
      /// Probably there's no way of adptating this method if the
      /// UKP5 arrays are changed to slices.
      virtual string gen_info(void) {
        ukp5_gen_stats();

        stringstream out("");
        out << "ukp5 used conf follows:" << endl;
        conf.print(out);
        out << "end of ukp5 conf" << endl;

        #ifdef HBM_CHECK_PERIODICITY
        out << "last_y_value_outer_loop: " << last_y_value_outer_loop << endl;
        #endif
        out << "y_bound: " << y_bound << endl;
        out << "c: " << c << endl;
        out << "n: " << n << endl;
        out << "w_min: " << w_min << endl;
        out << "w_max: " << w_max << endl;
        out << "last_dy_non_zero_non_n: " << last_dy_non_zero_non_n << endl;
        out << "qt_non_skipped_ys: " << qt_non_skipped_ys << endl;
        out << "qt_gy_zeros: " << qt_gy_zeros << endl;
        out << "qt_inner_loop_executions: " << qt_inner_loop_executions << endl;
        out << "qt_inner_loop_executions/qt_non_skipped_ys: " << ((long double) qt_inner_loop_executions)/((long double) qt_non_skipped_ys) << endl;
        out << "qt_inner_loop_executions/c: " << ((long double) qt_inner_loop_executions)/((long double) c) << endl;
        out << "(qt_inner_loop_executions/qt_non_skipped_ys)/n: " << ((long double) qt_inner_loop_executions)/((long double) qt_non_skipped_ys)/((long double)n)  << endl;
        out << "(qt_inner_loop_executions/c)/n: " << ((long double) qt_inner_loop_executions)/((long double) c)/((long double)n) << endl;

        const double &stime = sort_time, &dtime = smdom_time,
          &vtime = vector_alloc_time, &lctime = linear_comp_time,
          &p1time = phase1_time, &p2time = phase2_time, &ttime = total_time;
        const double sum_time =
          stime + dtime + vtime + lctime + p1time + p2time;

        streamsize old_precision = out.precision(HBM_PROFILE_PRECISION);
        const int two_first_digits_and_period = 3;
        const int percent_size = two_first_digits_and_period+HBM_PROFILE_PRECISION;
        ios_base::fmtflags old_flags = out.setf(std::ios::fixed, std:: ios::floatfield);
        char old_fill = out.fill(' ');

        out << "Sort time: " << stime << "s (";
        out << setw(percent_size) << (stime/ttime)*100.0;
        out << "%)" << endl;
        out << "Dom  time: " << dtime << "s (";
        out << setw(percent_size) << (dtime/ttime)*100.0;
        out << "%)" << endl;
        out << "Vect time: " << vtime << "s (";
        out << setw(percent_size) << (vtime/ttime)*100.0;
        out << "%)" << endl;
        out << "O(n) time: " << lctime << "s (";
        out << setw(percent_size) << (lctime/ttime)*100.0;
        out << "%)" << endl;
        out << "pha1 time: " << p1time << "s (";
        out << setw(percent_size) << (p1time/ttime)*100.0;
        out << "%)" << endl;
        out << "pha2 time: " << p2time << "s (";
        out << setw(percent_size) << (p2time/ttime)*100.0;
        out << "%)" << endl;
        out << "Sum  time: " << sum_time << "s (";
        out << setw(percent_size) << (sum_time/ttime)*100.0;
        out << "%)" << endl;
        out << "All  time: " << ttime << "s" << endl;

        out.fill(old_fill);
        out.setf(old_flags);
        out.precision(old_precision);

        if (conf.create_dumps) {
          dump(conf.gd_path, "y\tgy", g);
          dump(conf.dd_path, "y\tdy", d);
          dump(conf.nsd_path, "y\tdy", non_skipped_d);
          dump(conf.dqt_path, "i\tqt_in_d", qt_i_in_dy);
        }

        return out.str();
      }

      /// This method set the majority of the values of this stats class. All
      /// the values set by it need extra computation, and because of this they
      /// are computed here, after the UKP5 executed, and without affecting its
      /// time measurement. This method only works if the following variables
      /// are set first: n, c, y_bound, w_min, w_max, g and d. It can't receive
      /// this values by constructor as g and d are initialized empty and
      /// changed by UKP5 (they aren't a copy of the vectors used by ukp5, they
      /// are the vectors used by ukp5).
      void ukp5_gen_stats(void) {
        non_skipped_d.assign(y_bound-w_min + 1, n-1);

        qt_i_in_dy.assign(n, 0);
        qt_i_in_dy[n-1] = w_min;

        qt_gy_zeros = w_min;
        P opt = 0;
        qt_non_skipped_ys = 0;
        qt_inner_loop_executions = 0;

        for (W y = w_min; y <= y_bound-w_min; ++y) {
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
  
    /// Gets the minimal and maximal item weights and return them.
    template<typename W, typename P>
    pair<W, W> minmax_item_weight(const vector< item_t<W, P> > &items) {
      W min, max;
      min = max = items[0].w;
      for (auto it = items.begin()+1; it != items.end(); ++it) {
        W x = it->w;
        if (x < min) min = x;
        else if (x > max) max = x;
      }
      return make_pair(min, max);
    }

    /// Internal use helper struct, DO NOT DEPEND. Used to gather
    /// all O(n) information about the instance on only one pass.
    template<typename W, typename P, typename I>
    struct bag_t {
      W w_max, w_min;
      item_t<W, P> b1, b2;
      I b1_ix;
    };

    /// Gather the w_min, w_max, the two most efficient items
    /// (b1 and b2, respectively) and the index of b1 in one 
    /// pass by the items vector. Return all this info in the
    /// bag struct.
    template<typename W, typename P, typename I>
    bag_t<W, P, I> minmax_w_and_two_best_items(const vector< item_t<W, P> > &items) {
      bag_t<W, P, I> b;
      b.w_min = b.w_max = items[0].w;
      b.b1 = b.b2 = items[0];
      b.b1_ix = 0;
      for (auto it = items.begin()+1; it != items.end(); ++it) {
        // Two best items and index of the best
        item_t<W, P> i = *it;
        if (i < b.b1) {
          b.b1_ix = it - items.begin();
          b.b2 = b.b1;
          b.b1 = i;
        } else if (i < b.b2) {
          b.b2 = i;
        }
        // Weight max/min
        W x = i.w;
        if (x < b.w_min) b.w_min = x;
        else if (x > b.w_max) b.w_max = x;
      }

      return b;
    }

    /// @brief Examine the range of the g vector where the
    /// optimal solution value is guaranteed to be, and retrieve
    /// it. Can only be used after ukp5_phase1 executed over g.
    ///
    /// @param c The knapsack capacity.
    /// @param g The g vector used by UKP5, after ukp5_phase1
    ///   was executed over it (i.e. after being populated with
    ///   solutions profit values).
    /// @param w_min The smallest item weight on the knapsack.
    /// @return The optimal solution value and its weight,
    ///   respectively.
    template<typename W, typename P>
    pair<P, W> get_opts(W c, const vector<P> &g, W w_min) {
      P opt = 0;
      // Dont need to be initialized, but initializing to stop
      // compiler warning messages
      W y_opt = 0;

      for (W y = (c-w_min)+1; y <= c; ++y) {
        if (opt < g[y]) {
          opt = g[y];
          y_opt = y;
        }
      }

      return make_pair(opt, y_opt);
    }

    /// Retrieves the optimal solution (as an item multiset). 
    /// In other words, it sets sol.used_items.
    /// 
    /// Follows an explanation of this procedure. The ukp5_phase1
    /// populates the "d" array with the index of items used in
    /// solutions (including the optimal solutions). At d[y] rests
    /// the index of an item present in a solution with weight equal
    /// to y. The ukp5_phase1 also sets sol.y_opt to the weight of the
    /// lightest optimal solution. Therefore, d[sol.y_opt] has the index
    /// of an item used on the lightest optimal solution (between all
    /// solutions with the same weight, d stores an item index used by the
    /// best of them). After subtracting the weight of the item with index
    /// d[sol.y_opt] we have the weight of the remaining part of the
    /// solution. We repeat this process until we had found every item that
    /// composes the lightest optimal solution.
    /// 
    /// @param items The same items, in the same order, that were passed
    ///   to ukp5_phase1.
    /// @param d The d vector set by ukp5_phase1.
    /// @param sol The struct with the field used_items, that is set by this
    ///   procedure.
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

    /// Retrieves the optimal solution (as a profit and weight value).
    /// In other words, it sets sol.opt and sol.y_opt. Also the g and d vectors
    /// populated by it are used to get the optimal solution item multiset
    /// (ukp5_phase2) and the ukp5 stats (ukp5_extra_info_t.gen_info) THIS IS
    /// THE CORE OF THE UKP5 ALGORITHM. Is in this phase that the optimal
    /// solution value is obtained, and is in this phase that the ukp5 method
    /// spends 99% of the processing time.
    /// 
    /// Follows an explanation of this procedure. The ideia is to store at g/d
    /// all the solutions composed by one unique item (the weight is always
    /// stored implicity as the g/d index). After this we iterate g and stop at
    /// the already stored solutions to create new solutions combining "all"
    /// the items with the current solution. We don't combine all items, we
    /// combine only the ones until d[y], because this prunes symmetry (all
    /// solutions are yet generated but we avoid the same solution generated in
    /// many different orders, i.e. instead of {1,2,3}, {2,1,3}, ..., and
    /// {3,2,1} we have only {3,2,1}). Also we skip some inferior/dominated
    /// solutions (if g[y] < opt then the solution in g[y] is less profitable
    /// than a previous solution, and therefore can be discarded without loss
    /// to the optimality). The periodicity check simply checks if between the
    /// current position (y) and the last (c-w_min+w_max) we will add an item
    /// that is not the best item (g[y] > 0 && d[y] > 0). If not, then we know
    /// that the remaining capacity will be filled with copies of the first
    /// item (index zero). This is the gist of it.
    /// 
    /// @param ukpi The UKP instance.
    /// @param g This vector will be used to store the profit value of
    ///   knapsack solutions. It can't have non-zero values, or the wrong
    ///   size (less than c+w_max-w_min).
    /// @param d This vector will be used to store the index of the last
    ///   item used to compose a solution whose profit value is stored at g.
    ///   The data of this vector makes reference to the order of ukpi.items.
    ///   It can't have non-zero values, or the wrong size (less than
    ///   c+w_max-w_min).
    /// @param sol This procedure sets the fields opt and y_opt of the
    ///   sol parameter. If HBM_CHECK_PERIODICITY is defined, the field
    ///   extra_info is expect to be set (already point to a valid
    ///   ukp5_extra_info_t<W, P, I> object), as the field
    ///   sol.extra_info->last_y_value_outer_loop will be set by the procedure.
    /// @param w_min The minimal item weight between the items in ukpi.items.
    /// @param w_max The maximal item weight between the items in ukpi.items.
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

      // This block is a copy-past of the loop below, but only for the best
      // item. Its utility is to simplify the code when HBM_CHECK_PERIODICITY
      // is defined.
      W wb = items[0].w;
      g[wb] = items[0].p;;
      d[wb] = 0;

      for (I i = 1; i < n; ++i) {
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

        // This block is a copy-past of the loop below, but only for the best
        // item. Its utility is to simplify the code when HBM_CHECK_PERIODICITY
        // is defined.
        item_t<W, P> bi = items[0];
        W wb = bi.w;
        P pb = bi.p;
        W next_y = y + wb;
        P old_gny = g[next_y];
        P new_gny = gy + pb;
        if (old_gny <= new_gny) {
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
          } else if (ogny == ngny && d[ny] > ix) {
            g[ny] = ngny;
            d[ny] = ix;
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

    double inline difftime_between_now_and(const steady_clock::time_point &begin) {
      return duration_cast< duration<double> >(steady_clock::now() - begin)
             .count();
    }

    template<typename W, typename P, typename I>
    void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, const ukp5_conf_t<I> &conf) {
      steady_clock::time_point all_ukp5_begin = steady_clock::now();
      steady_clock::time_point begin; // Only declare.
      // Saves the original c and n, before any periodicity test
      // or dominance item removal is done.
      const W c = ukpi.c;
      const I n = ukpi.items.size();

      // This pointer stores all the extra info that the UKP5 will print or
      // dump after.
      ukp5_extra_info_t<W, P, I>* ptr = new ukp5_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(ptr);
      sol.extra_info = shared_ptr<extra_info_t>(upcast_ptr);

      // We print the configuration used later.
      ptr->conf = conf;

      // Apply or not the simple/multiple dominance removal.
      if (conf.apply_smdom && conf.apply_smdom_before_sort) {
        begin = steady_clock::now();
        vector< item_t<W, P> > undominated;
        smdom(ukpi.items, undominated, 0, true);
        if (undominated.size() < ukpi.items.size()) {
          ukpi.items = undominated;
        }
        ptr->smdom_time = difftime_between_now_and(begin);
      }

      // k: how many of the most efficient items will be moved to front
      // (this 'if' doesn't sort anything, only computes the value
      // of k)
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
      ptr->sort_time = difftime_between_now_and(begin);

      if (conf.apply_smdom && !conf.apply_smdom_before_sort) {
        begin = steady_clock::now();
        vector< item_t<W, P> > undominated;
        smdom(ukpi.items, undominated, k, true);
        if (undominated.size() < ukpi.items.size()) {
          ukpi.items = undominated;
        }
        ptr->smdom_time = difftime_between_now_and(begin);
      }

      begin = steady_clock::now();
      W w_min, w_max; // Always initialized below.
      item_t<W, P> b1, b2; // Initialized if conf.use_y_star_per is true.
      // Changed if conf.use_y_star_per is true, and the bound computed
      // is smaller than c. Initialized here to remove compiler warnings.
      W y_bound = c;
      I b1_ix = 0; // Initialized to remove compiler warnings.

      if (conf.use_y_star_per) {
        if (k < 2) {
          auto bag = minmax_w_and_two_best_items<W, P, I>(ukpi.items);
          w_min = bag.w_min;
          w_max = bag.w_max;
          b1_ix = bag.b1_ix;
          b1 = bag.b1;
          b2 = bag.b2;
        } else {
          auto min_max_w = minmax_item_weight(ukpi.items);
          w_min = min_max_w.first;
          w_max = min_max_w.second;
          b1_ix = 0;
          b1 = ukpi.items[0];
          b2 = ukpi.items[1];
        }

        W y_ = y_star(b1, b2);
        y_bound = refine_y_star(y_, c, b1.w);
      } else {
        auto min_max_w = minmax_item_weight(ukpi.items);
        w_min = min_max_w.first;
        w_max = min_max_w.second;
      }
      ptr->linear_comp_time = difftime_between_now_and(begin);

      // Set some variables for post printing
      ptr->n = n;
      ptr->c = c;
      ptr->y_bound = y_bound;
      ptr->w_min = w_min;
      ptr->w_max = w_max;

      begin = steady_clock::now();
      // Use solution_t.extra_info fields instead of local variables to
      // propagate the array values. The arrays will be dumped to files,
      // making possible study them with R or other tool.
      vector<P> &g = ptr->g;
      vector<I> &d = ptr->d;
      g.assign(y_bound+1+(w_max-w_min), 0);
      d.assign(y_bound+1+(w_max-w_min), n-1);
      ptr->vector_alloc_time = difftime_between_now_and(begin);

      begin = steady_clock::now();
      if (conf.use_y_star_per) ukpi.c = y_bound;
      ukp5_phase1(ukpi, g, d, sol, w_min, w_max);
      if (conf.use_y_star_per) ukpi.c = c;
      ptr->phase1_time = difftime_between_now_and(begin);

      begin = steady_clock::now();
      ukp5_phase2(ukpi.items, d, sol);
      ptr->phase2_time = difftime_between_now_and(begin);

      // This can be a little weird at first but we have two periodicity
      // checks. The one right below is inside ukp5_phase1, and only works
      // if the first element of the vector is the best one. If the first
      // item isn't the best one it simply do nothing (the code work
      // normally). The one after this one uses the y* (y_star) bound,
      // as this ones changes ukpi.c to the value of y_bound before passing
      // it to ukp5_phase1, we need to do the first one based on the y_bound
      // value and then compute the second one based on the original c.
      #ifdef HBM_CHECK_PERIODICITY
      // If we use our periodicity check the sol.used_items constructed by
      // ukp5_phase2 doesn't include the copies of the best item used to
      // fill all the extra_capacity. This solves it.
      W qt_best_item_inserted_by_per = (y_bound - sol.y_opt)/ukpi.items[0].w;
      sol.opt += static_cast<P>(qt_best_item_inserted_by_per)*ukpi.items[0].p;
      sol.y_opt += qt_best_item_inserted_by_per*ukpi.items[0].w;
      // Our periodicity check will always stop with a partial solution
      // that have at least one of the best item. And ukp5_phase2 will
      // populate the sol.used_items in the same order the sol.items are,
      // this way we have the guarantee that the first element of
      // sol.used_items exist and it's the best item.
      // NOTE: this works even if the array isn't sorted, because if the
      // best item isn't the first, then our periodicity check never
      // stops the computation before the end. If our periodicity check
      // stopped the computation (the statement c - sol.y_opt > items[0].w
      // will be true in this case), we can safely assume that the
      // best item is the first on the solution.
      sol.used_items[0].qt += qt_best_item_inserted_by_per;
      #endif

      // This is the second periodicity check we talked before.
      if (conf.use_y_star_per && y_bound < c) {
        W qt_b = (c - sol.y_opt)/b1.w;
        sol.opt += static_cast<P>(qt_b) * b1.p;
        sol.y_opt += qt_b * b1.w;

        auto eq_b1 = [&](const itemqt_t<W, P, I> &i) { return i.it == b1; };
        auto bp = find_if(sol.used_items.begin(), sol.used_items.end(), eq_b1);
        if (bp == sol.used_items.end()) {
          // The best item wasn't used on the reduced capacity optimal
          // solution, we need to add it.
          auto bqt = itemqt_t<W, P, I>(b1, qt_b, b1_ix);
          sol.used_items.push_back(bqt);
        } else {
          bp->qt += qt_b;
        }
      }

      ptr->total_time = difftime_between_now_and(all_ukp5_begin);

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

      const regex float_0_1("(0|1|0\\.[0-9]+|1\\.0+)", nosubs),
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
        // Exists only to avoid giving warning if the method was called
        // without any extra options.
        break;

        case 1:
        case 0:
        cerr << "WARNING: " << __func__ << ": argc < 2?! Executing method"
                " with default parameter values."
        << endl;
        usage(__func__);
        break;

        default:
        cout << "WARNING: " << __func__ << ": more than " << argc << " "
              "parameters to ukp5. Will execute with the default values,"
              " not with the values provided by command line."
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
        if (conf.nsd_path.empty()) {
          conf.nsd_path = "./nsd_dump_" + fname.native() + ".dat";
        }
        if (conf.dqt_path.empty()) {
          conf.dqt_path = "./dqt_dump_" + fname.native() + ".dat";
        }
      }

      hbm_ukp5_impl::ukp5(ukpi, sol, conf);
    }
  }

  /// Solves an UKP instance by the UKP5 algorithm, and stores the results (and
  /// stats) at sol.
  ///
  /// @param ukpi The UKP instance to be solved.
  /// @param sol The object where the results will be written.
  /// @param conf An object with the UKP5 options.
  ///
  /// @see ukp5_conf_t For the explanation of the UKP5 options.
  template<typename W, typename P, typename I>
  void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, const ukp5_conf_t<I> &conf) {
    hbm_ukp5_impl::ukp5(ukpi, sol, conf);
  }

  /// For information on argv valid arguments consult the procedure "usage",
  /// or call this function with argv = null_ptr and argc = 0.
  template<typename W, typename P, typename I>
  void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_ukp5_impl::ukp5(ukpi, sol, argc, argv);
  }
}

#endif //HBM_UKP5_HPP
