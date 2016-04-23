#ifndef HBM_GREENDP_HPP
#define HBM_GREENDP_HPP

#include "ukp_common.hpp"
#include "periodicity.hpp"
#include "wrapper.hpp"
//#include "prettyprint.hpp"

#include <vector> // For vector
#include <forward_list>
#include <algorithm>  // For reverse
#include <boost/math/common_factor_rt.hpp>  // For boost::math::gcd
#include <chrono> // For steady_clock::
#include <boost/rational.hpp> // For rational<T> used at greendp2

#ifndef HBM_PROFILE_PRECISION
  #define HBM_PROFILE_PRECISION 5
#endif

#ifndef HBM_PRINT_VAR
  #define HBM_PRINT_VAR(var) out << #var ": " << var << std::endl
#endif

namespace hbm {
  template <typename W, typename P, typename I>
  struct greendp_extra_info_t : extra_info_t {
    /// The last capacity computed before detecting periodicity and stoping.
    W last_y_value_outer_loop{0};
    #ifdef HBM_PROFILE
    double sort_time{0};        ///< Time used sorting items.
    double vector_alloc_time{0};///< Time used allocating vectors for DP.
    double linear_comp_time{0}; ///< Time used by linear time preprocessing.
    double dp_time{0};     ///< Time used creating partial solutions.
    double sol_time{0};    ///< Time used to assemble solution.
    double bound_time{0};  ///< Time used computing bounds.
    double total_time{0};  ///< Time used by all the algorithm.
    #endif //HBM_PROFILE
    I n{0};      ///< Instance number of items.
    W c{0};      ///< Instance capacity.

    virtual std::string gen_info(void) {
      std::stringstream out("");

      HBM_PRINT_VAR(last_y_value_outer_loop);
      HBM_PRINT_VAR(c);
      HBM_PRINT_VAR(n);

      #ifdef HBM_PROFILE
      const double sum_time = sort_time + vector_alloc_time +
        linear_comp_time + dp_time + sol_time + bound_time;

      std::streamsize old_precision = out.precision(HBM_PROFILE_PRECISION);
      const int two_first_digits_and_period = 3;
      const int percent_size = two_first_digits_and_period+HBM_PROFILE_PRECISION;
      std::ios_base::fmtflags old_flags = out.setf(std::ios::fixed, std:: ios::floatfield);
      char old_fill = out.fill(' ');

      #ifndef HBM_PRINT_TIME
        #define HBM_PRINT_TIME(name, var)\
          out << name << " time: " << var << "s (";\
          out << std::setw(percent_size) << (var/total_time)*100.0;\
          out << "%)" << std::endl
      #endif

      HBM_PRINT_TIME("Sort", sort_time);
      HBM_PRINT_TIME("Vect", vector_alloc_time);
      HBM_PRINT_TIME("O(n)", linear_comp_time);
      HBM_PRINT_TIME("dp  ", dp_time);
      HBM_PRINT_TIME("bd  ", bound_time);
      HBM_PRINT_TIME("sol ", sol_time);
      HBM_PRINT_TIME("Sum ", sum_time);
      HBM_PRINT_TIME("All ", total_time);

      out.fill(old_fill);
      out.setf(old_flags);
      out.precision(old_precision);
      #endif //HBM_PROFILE

      return out.str();
    }
  };

  /// Same as greendp_extra_info_t, created only as a contingence if, in the
  /// future, we have different stats between greendp and mgreendp.
  template <typename W, typename P, typename I>
  struct mgreendp_extra_info_t : greendp_extra_info_t<W, P, I> {};

  /// Same as greendp_extra_info_t, created only as a contingence if, in the
  /// future, we have different stats between greendp and greendp1.
  template <typename W, typename P, typename I>
  struct greendp1_extra_info_t : greendp_extra_info_t<W, P, I> {
    W t; ///< Value of constant t (it's constant for each instance).
    W m; ///< Final value of variable m (number of generated solutions).
    /// GCD of all profit values. If this is bigger than one, then all the
    /// computations used profits divided by gcd_n, and at the end of the
    /// algorithm, z and the optimal solution items were returned to they
    /// original value to be displayed correctly.
    P gcd_c;

    std::string gen_info(void) {
      std::string s = greendp_extra_info_t<W, P, I>::gen_info();
      std::stringstream out("");

      HBM_PRINT_VAR(t);
      HBM_PRINT_VAR(m);
      HBM_PRINT_VAR(gcd_c);

      return s + out.str();
    }
  };

  template <typename W, typename P, typename I>
  struct mgreendp1_extra_info_t : greendp1_extra_info_t<W, P, I> {};

  template <typename W, typename P, typename I>
  struct greendp2_extra_info_t : extra_info_t {
    #ifdef HBM_PROFILE
    double sort_time{0};        ///< Time used sorting items.
    double vector_alloc_time{0};///< Time used allocating vectors for DP.
    double linear_comp_time{0}; ///< Time used by linear time preprocessing.
    double dom_time{0};    ///< Time used dominance tests.
    double dp_time{0};     ///< Time used creating partial solutions.
    double sol_time{0};    ///< Time used to assemble solution.
    double total_time{0};  ///< Time used by all the algorithm.
    #endif //HBM_PROFILE
    I n{0}; ///< Instance number of items.
    I n2{0};///< Number of items not excluded before starting.
    W c{0}; ///< Instance capacity.
    W m; ///< Final value of variable m (number of generated solutions).
    /// GCD of all profit values. If this is bigger than one, then all the
    /// computations used profits divided by gcd_n, and at the end of the
    /// algorithm, z and the optimal solution items were returned to they
    /// original value to be displayed correctly.
    P gcd_c;

    virtual std::string gen_info(void) {
      std::stringstream out("");

      HBM_PRINT_VAR(c);
      HBM_PRINT_VAR(n);
      HBM_PRINT_VAR(n2);
      HBM_PRINT_VAR(m);

      #ifdef HBM_PROFILE
      const double sum_time = sort_time + vector_alloc_time +
        linear_comp_time + dp_time + dom_time + sol_time;

      std::streamsize old_precision = out.precision(HBM_PROFILE_PRECISION);
      const int two_first_digits_and_period = 3;
      const int percent_size = two_first_digits_and_period+HBM_PROFILE_PRECISION;
      std::ios_base::fmtflags old_flags = out.setf(std::ios::fixed, std:: ios::floatfield);
      char old_fill = out.fill(' ');

      #ifndef HBM_PRINT_TIME
        #define HBM_PRINT_TIME(name, var)\
          out << name << " time: " << var << "s (";\
          out << std::setw(percent_size) << (var/total_time)*100.0;\
          out << "%)" << std::endl
      #endif

      HBM_PRINT_TIME("Sort", sort_time);
      HBM_PRINT_TIME("Vect", vector_alloc_time);
      HBM_PRINT_TIME("O(n)", linear_comp_time);
      HBM_PRINT_TIME("Dom ", linear_comp_time);
      HBM_PRINT_TIME("dp  ", dp_time);
      HBM_PRINT_TIME("sol ", sol_time);
      HBM_PRINT_TIME("Sum ", sum_time);
      HBM_PRINT_TIME("All ", total_time);

      out.fill(old_fill);
      out.setf(old_flags);
      out.precision(old_precision);
      #endif //HBM_PROFILE

      return out.str();
    }
  };

  template <typename W, typename P, typename I>
  struct mgreendp2_extra_info_t : greendp2_extra_info_t<W, P, I> {};

  namespace hbm_greendp_impl {
    using namespace std;
    using namespace std::chrono;
    using namespace boost;
    using namespace boost::math;

    /// Compute the gcd for a variable quantity of numbers.
    /// @note This functions was designed with positive integers on mind.
    template<typename ForwardIterator, typename T>
    inline T gcd_n(ForwardIterator begin, ForwardIterator end) {
      if (begin == end) {
        return 1;
      }

      if (begin + 1 == end) {
        return *begin;
      }

      T acc = boost::math::gcd(*begin, *(begin+1));
      ++begin;

      for (; begin != end && acc > 1; ++begin) {
        acc = boost::math::gcd(acc, *begin);
      }

      return acc;
    }

    /// Compute k_max, also, if returns false the greendp procedure (caller)
    /// is to be stopped.
    template<typename W, typename P, typename I>
    inline bool compute_k_max(
      const I &m,
      const vector<P> &b_,
      const vector<W> &a,
      const vector<P> &c,
      const W &am,
      const P &cm,
      const W &b,
      const P &z,
      const P &d,
      const W &k,
      const W &lambda,
      W &k_max
    ) {
      // STEP 1 (Routine k_max)
      I r = m - 1;
      while (r > 0 && b_[r] > cm*b - am*(z + d)) --r;

      if (r == 0) return false;

      // STEP 2 (Routine k_max)
      const P delta = z - cm*(b/am) + d;
      const P l_side = (c[r]*lambda - a[r]*delta)/b_[r];
      const P r_side = b/am;
      k_max = l_side < r_side ? l_side : r_side;

      if (k > k_max) return false;

      return true;
    }

    /// Compute bi_qt_lb, also, if returns false the greendp procedure (caller)
    /// is to be stopped.
    template<typename W, typename P>
    inline bool compute_bi_qt_lb(const vector<P> &b,
                          const vector< item_t<W, P> > &items,
                          const W &c,
                          const P &opt,
                          const P &d,
                          const W &bi_qt,
                          W &bi_qt_lb) {
      // STEP 1 (Routine k_max)
      const item_t<W, P> &bi = items.front();
      const size_t n = items.size();

      size_t r = 1;
      const P eff_diff_bi_and_opt_plus_one = bi.p*c - (opt + d)*bi.w;
      while (r < n && b[r] > eff_diff_bi_and_opt_plus_one) ++r;
      // Let's examine the condition above:
      // b[r] > eff_diff_bi_and_opt_plus_one | original
      // bi.p*it.w - it.p*bi.w > bi.p*c - (opt + d)*bi.w | expansion
      // (bi.p*it.w)/bi.w - it.p > (bi.p*c)/bi.w - (opt + d) | div by bi.w
      // (bi.p/bi.w * it.w) - it.p > (bi.p/bi.w * c) - (opt + d) | rearrange
      // So, basically what's happening here is that the we are multiplying the
      // efficiency of the best item (bi.p/bi.w) with two weights (the weight
      // of an item on the left side, and the knapsack max weight at the right
      // side) and subtracting the 'natural' profit of that weight (the profit
      // of an item on the left side, and the best solution found until now at
      // the right side). For the same weight (i.e it.w == bi.w and bi.w == c)
      // we have bi.p - it.p and bi.p - (opt + d). This way we can easily
      // verify that for a more efficient 'it' or opt we have a smaller value
      // (always above zero, as bi is the most efficient item and (bi.p/bi.w)*c
      // is the relaxed solution). If b[r] is greater than
      // eff_diff_bi_and_opt_plus_one this means 'it' is LESS efficient than
      // the solution, and therefore can't improve the solution. Remember that
      // the current solution is f[c-bi_qt*bi.w]+bi_qt*bi.p, or in other words,
      // the guaranteed optimal solution value for c-bi_qt*bi.w more copies of
      // the most efficient item. The current solution can be more efficient
      // than the second most efficient item (items[1]) because opt is
      // reasonably efficient and some copies of bi are sufficient to make the
      // current solution more efficient than the second most efficient item.
      // Also, as f[c-bi_qt*bi.w] is a guaranteed optimal value, the only way
      // to insert a copy of item 'it' to the solution is removing one or more
      // copies of the best item.
      if (r == n) return false;
      // if r == n, then no item can improve the solution

      // lambda is the remaining space if we fill the capacity c
      // with the most efficient item.
      const W lambda = c - (c / bi.w) * bi.w;
      const W max_bi_qt = c/bi.w;

      // delta is the profit value of lambda in the current solution, if we
      // assume that after lambda the remaining capacity is filled with copies
      // of the best item.
      const P delta = (opt + d) - bi.p*max_bi_qt;
      const P l_side = (items[r].p*lambda - items[r].w*delta)/b[r];
      bi_qt_lb = max_bi_qt - std::min(l_side, max_bi_qt);

      if (bi_qt < bi_qt_lb) return false;

      return true;
    }

    /// Compute next_z, also, if returns false the greendp procedure (caller)
    /// is to be stopped.
    template<typename W, typename P, typename I>
    bool next_z(const I m,
                const vector<P> &b_,
                const vector<W> &a,
                const vector<P> &c,
                const W am,
                const P cm,
                const W b,
                P &z,
                const P d,
                W &k,
                const W lambda,
                W &k_max,
                const vector<P> &f,
                const W y,
                W &upper_l) {
      // STEP 1 (Routine next z)
      if (f[y] + cm*(b/am) - cm*k > z) {
        z = f[y] + cm*(b/am) - cm*k;
        if (!compute_k_max(m, b_, a, c, am, cm, b, z, d, k, lambda, k_max)) {
          return false;
        }
      }

      // STEP 2 (Routine next z)
      k = k + 1;
      if (k > k_max) return false;

      upper_l = upper_l + am;
      return true;
    }

    // This function is used to initialize b as const.
    // If no item share the same efficiency as the best item
    // then all b[i] > 0 for i=1..(n-1) (b[0] = 0 always)
    template<typename W, typename P>
    vector<P> compute_b(const vector< item_t<W, P> > &items) {
      const size_t n = items.size();
      const item_t<W, P> &bi = items.front();

      vector<P> b;
      b.reserve(items.size());
      b.push_back(0);

      for (size_t i = 1; i < n; ++i) {
        const item_t<W, P> &it = items[i];
        b.push_back(bi.p*it.w - it.p*bi.w);
      }

      return b;
    }

    // This function is used to initialize b_ as const.
    template<typename W, typename P, typename I>
    vector<W> compute_b_(const I m, const vector<W> &a, const vector<P> &c) {
      vector<P> b_(m); // it goes until m-1
      b_[0] = 0; // does not exist, notation begins at 1

      for (I i = 1; i <= m - 1; ++i) {
        b_[i] = c[m]*a[i] - c[i]*a[m];
      }

      return b_;
    }

    /// Modified greendp, or Modern greendp: essentially the same algorithm,
    /// but using loops (instead of gotos), and without trying to be a carbon
    /// copy of the article notation and step organization.
    template<typename W, typename P, typename I>
    void mgreendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      mgreendp_extra_info_t<W, P, I>* eip = new mgreendp_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_mgreendp_begin = steady_clock::now();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      HBM_START_TIMER();
      vector< item_t<W, P> > &items = ukpi.items;
      if (!already_sorted) sort_by_eff(items);
      HBM_STOP_TIMER(eip->sort_time);

      HBM_START_TIMER();
      const I n = eip->n = static_cast<I>(items.size());
      const W c = eip->c = ukpi.c;
      // bi stands for 'best item'
      const item_t<W, P> bi = items[0];

      vector<P> p;
      p.reserve(n);
      for (I i = 0; i < n; ++i) p.push_back(items[i].p);
      const P d = gcd_n<typename vector<P>::iterator, P>(p.begin(), p.end());
      HBM_STOP_TIMER(eip->linear_comp_time);

      HBM_START_TIMER();
      vector<P> f(c + 1, 0);
      myvector<I> i;
      i.resize(c + 1); // this do not initialize the vector
      HBM_STOP_TIMER(eip->vector_alloc_time);

      HBM_START_TIMER();
      // z is another name for sol.opt
      P &z = sol.opt;
      // We init z with the profit value of the solution composed by the
      // maximum number of best items that the knapsack can hold.
      W bi_qt = (c / bi.w);
      z = bi.p * bi_qt;

      // STEP 1
      const vector<P> b = compute_b(items);
      HBM_STOP_TIMER(eip->linear_comp_time);

      HBM_START_TIMER();
      W bi_qt_lb;

      if (!compute_bi_qt_lb(b, items, c, z, d, bi_qt, bi_qt_lb)) {
        return;
      }
      HBM_STOP_TIMER(eip->bound_time);

      HBM_START_TIMER();
      for (I j = 1; j < n; ++j) {
        const W first_y_reserved_for_best_items = (c - bi_qt_lb*bi.w) + 1;
        const item_t<W, P> it = items[j];
        if (it.w < first_y_reserved_for_best_items && it.p > f[it.w]) {
          f[it.w] = it.p;
          i[it.w] = j;
        }
      }
      HBM_STOP_TIMER(eip->dp_time);

      W &y = eip->last_y_value_outer_loop = 1;
      W y_opt_no_bi = 0;
      W bi_qt_in_z = 0;
      P max_previous_fy = 0;
      for (; bi_qt >= bi_qt_lb; --bi_qt) {
        HBM_START_TIMER();
        for (; y < c - bi_qt*bi.w; ++y) {
          if (f[y] <= max_previous_fy) continue;

          max_previous_fy = f[y];

          // STEP 2a, 2b, 2c, 2d
          // This is very similar to the loop over items of UKP5.
          // The difference is the j=1 and the first_y_reserved_for_best_items
          const W first_y_reserved_for_best_items = (c - bi_qt_lb*bi.w) + 1;
          for (I j = 1; j <= i[y]; ++j) {
            const item_t<W, P> &it = items[j];
            const W new_y = y + it.w;
            if (new_y < first_y_reserved_for_best_items
                && it.p + f[y] > f[new_y]) {
              f[new_y] = it.p + f[y];
              i[new_y] = j;
            }
          }
        }
        HBM_STOP_TIMER(eip->dp_time);

        // STEP 3d
        // STEP 1 (Routine next z)
        // Note that now y == c - bi_qt*bi.w,
        // and max_previous_fy == max(f[0..y-1])
        HBM_START_TIMER();
        const P z_ = std::max(max_previous_fy, f[y]) + bi.p*bi_qt;
        if (z_ > z) {
          z = z_;
          y_opt_no_bi = y;
          bi_qt_in_z = bi_qt;

          if (!compute_bi_qt_lb(b, items, c, z, d, bi_qt, bi_qt_lb)) {
            break;
          }
        }
        HBM_STOP_TIMER(eip->bound_time);
        // The W type can be unsigned, so we have to stop before it underflows
        if (bi_qt == 0) break;
      }

      HBM_START_TIMER();
      vector<I> qts_its(n, 0);

      W y_opt = y_opt_no_bi;
      while (y_opt > 0 && f[y_opt] < z - bi_qt_in_z*bi.p) --y_opt;
      sol.y_opt = y_opt + bi_qt_in_z*bi.w;

      I dy_opt;
      while (y_opt != 0) {
        dy_opt = i[y_opt];
        y_opt -= items[dy_opt].w;
        ++qts_its[dy_opt];
      }

      for (I i = 0; i < n; ++i) {
        if (qts_its[i] > 0) {
          sol.used_items.emplace_back(items[i], qts_its[i], i);
        }
      }

      if (bi_qt_in_z > 0) {
        auto bqt = itemqt_t<W, P, I>(bi, bi_qt_in_z, 0);
        sol.used_items.push_back(bqt);
      }

      sol.used_items.shrink_to_fit();
      HBM_STOP_TIMER(eip->sol_time);
      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_mgreendp_begin);
      #endif
    }

    /// The modernized version (without goto) of the second algorithm presented
    /// at "On Equivalent Knapsack Problems", (H.  Greenberg).
    template<typename W, typename P, typename I>
    void mgreendp2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      mgreendp2_extra_info_t<W, P, I>* eip =
        new mgreendp2_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_greendp1_begin = steady_clock::now();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      // BEFORE STEP 1 (notation already introduced by article)
      HBM_START_TIMER();
      auto &items = ukpi.items;
      // In this paper the items are sorted by non-decreasing efficiency
      // and the ties are solved by ordering by non-increasing weight
      // (or profit).
      if (!already_sorted) sort_by_eff(items);
      HBM_STOP_TIMER(eip->sort_time);

      HBM_START_TIMER();
      // In this algorithm 'n' can change, as we can remove items of the
      // item list. We will use eip->n2 to denote the new n.
      const I n = eip->n = static_cast<I>(items.size());
      const W b = eip->c = ukpi.c;

      vector<W> a(n + 1);
      a[0] = 0; // Does not exist, notation begins at 1
      vector<P> c(n + 1);
      c[0] = 0; // Does not exist, notation begins at 1

      for (I i = 0; i < n; ++i) {
        item_t<W, P> it = items[i];
        a[i+1] = it.w;
        c[i+1] = it.p;
      }
      const P gcd_c = gcd_n<typename vector<P>::iterator, P>
                           (c.begin() + 1, c.end());
      if (gcd_c > 1) for (P &cj: c) cj /= gcd_c;
      eip->gcd_c = gcd_c;
      const W a1 = a[1];
      const P c1 = c[1];

      // PARAGRAPH ALGORITHM 2, before step 1
      myvector< rational<P> > d_;
      d_.resize(n + 1);
      const rational<P> frac_a1_c1 = rational<P>(a1, c1);
      for (I j = 2; j <= n; ++j) {
        // On the article this line of the algorithm includes divisions that
        // aren't rounded down. This would call for floating pointer
        // arithmetic. However, those values are compared for equality later.
        // On the article's example, the numbers are shown as fractions, what
        // implies that absolute precision can be necessary. With this in mind,
        // I chose to represent the values of d_ as fractions.
        // UPDATE: When executing this algorithm over the example provided
        // on the paper, the optimal solution is reached after comparing
        // two equal fractions (absolute precision is necessary).
        d_[j] = a[j] - (frac_a1_c1*c[j]);
        // this value is always non-negative (proof follows):
        // previous knowledge: c1/a1 >= cj/aj
        // we can derive:
        // c1 >= (cj/aj)*a1 (multiply by a1)
        // c1*aj >= cj*a1   (multiply by aj, basic efficiency comparation)
        // c1*aj/cj >= a1   (divide by cj)
        // aj/cj >= a1/c1   (divide by c1, inverse of the first statement)
        // aj >= (a1/c1)*cj (aj - (a1/c1)*cj is always non-negative)
        // Meaning of d_: d_ seems to be the difference between how
        // much the item would weight if it had the same efficiency as the best
        // item, and the actual weight of the item. If it had the same
        // efficiency as the best item, this would mean its efficiency would be
        // greater, and for the same profit the weight would be smaller (so the
        // difference with its original weight is always non-negative).
      }
      HBM_STOP_TIMER(eip->linear_comp_time);

      HBM_START_TIMER();
      // I reordered the tests to make them clear.
      vector<bool> items_to_remove(n + 1, false);
      I qt_items_to_remove = 0;
      for (I j = 2; j <= n; ++j) {
        if (c[j] % c1 == 0) {
          items_to_remove[j] = true;
          ++qt_items_to_remove;
          continue;
        }
        // The step one of the article don't define the range of 'i'
        // (an item index variable). I do believe that i=1 is excluded
        // because otherwise d_[1] would be referenced (and it's undefined).
        // I also believe that 'i' == 'j' can be skipped, as the conditions
        // for exclusion only evaluate to true if the items are different
        // (i.e. skipping i == j makes no difference).
        for (I i = 2; i <= n; ++i) {
          if (i == j) continue;
          if (c[i] % c1 == c[j] % c1) {
            if (d_[i] < d_[j]) {
              items_to_remove[j] = true;
              ++qt_items_to_remove;
              break;
            } else if (d_[i] == d_[j] && c[i] <= c[j]) {
              items_to_remove[j] = true;
              ++qt_items_to_remove;
              break;
            }
          }
        }
      }

      // The number of items that weren't excluded by step one.
      // And after, the number os solutions generated.
      I m;
      // Here we remove the items maked for removal.
      // Scope to guarantee vector deallocation.
      {
        m = eip->n2 = n - qt_items_to_remove;

        vector<W> a_;
        a_.reserve(m + 1);
        a_.push_back(0);

        vector<P> c_;
        c_.reserve(m + 1);
        c_.push_back(0);

        myvector< rational<P> > tmp_d_;
        tmp_d_.reserve(m + 1);
        tmp_d_.push_back(0);

        for (I j = 1; j <= n; ++j) {
          if (!items_to_remove[j]) {
            a_.push_back(a[j]);
            c_.push_back(c[j]);
            tmp_d_.push_back(d_[j]);
          }
        }

        a_.swap(a);
        c_.swap(c);
        tmp_d_.swap(d_);
      }
      HBM_STOP_TIMER(eip->dom_time);

      HBM_START_TIMER();
      P x, y;
      rational<P> d;
      //step_1: // unused step
      // We will call c' as c_prime. The c_prime saves the value of
      // c[j] % c1 on position j. It's updated for j=2..n (the initial
      // items) and after this isn't changed anymore.
      myvector<P> c_prime;
      c_prime.resize(m + 1);
      // This isn't explicity said on the algorithm, but the size of all
      // vectors (except by upper_e) is c1. This is because x (the variable
      // used to index all vectors except by upper_e) is always smaller than c1
      // (as it's modulo c1).
      // As there's an array 'y' and a variable 'y', we will call
      // the array 'y_'.
      // y_ saves the profit value of solutions (z's), the key to access the
      // value is the saved value (z) % c1.
      vector<P> y_(c1, 0);
      // Every uppercase letter will be upper_x.
      // Each of the three vectors below store info related to a solution.
      // If the solution profit value is z, the info is stored at position
      // z % c1 of the array. The upper_g stores a fraction (meaning what?),
      // the upper_c stores the index of the solution on upper_e, and upper_d
      // stores the last item used to make the solution.
      vector<rational<P> > upper_g(c1, 0);
      vector<I> upper_c(c1, 0);
      vector<I> upper_d(c1, 0);

      // For upper_e we need a structure between a vector and a list. We want
      // fast access to one point inside the container and to its end (could be
      // done by list or vector). We want to clean from memory the elements
      // already iterated by our iterator in a O(1) way (list). We want to make
      // better use of cache and do not allocate memory at every new item
      // inserted at the container's end (vector). So the ideia is a vector of
      // vectors, the inner vectors have a max size, and we can ask to clean an
      // entire inner vector (i.e. chunk) easily.
      // If the P value is 64 bits, a 1000000 chunk_size is about 7MiB.
      // TODO: make this a configurable parameter
      // chunk-based upper_e code below
      const size_t chunk_size = 1000000;
      size_t h = 0, l = 0;
      vector< vector<P> > upper_e;
      upper_e.emplace_back();
      upper_e.front().reserve(chunk_size);
      upper_e.front().emplace_back(0);
      HBM_STOP_TIMER(eip->vector_alloc_time);

      HBM_START_TIMER();
      const rational<P> infinity(numeric_limits<P>::max(), 1);
      for (P k = 1; k < c1; ++k) {
        // The positions where k == c_prime[j] for j=2..m will be overwritten
        upper_g[k] = infinity;
      }

      for (I j = 2; j <= m; ++j) {
        c_prime[j] = c[j] % c1;
        upper_g[c_prime[j]] = d_[j];
        y_[c_prime[j]] = c[j];
        upper_c[c_prime[j]] = j;
        upper_d[c_prime[j]] = j;

        if (upper_e.back().size() == chunk_size) {
          upper_e.emplace_back();
          upper_e.back().reserve(chunk_size);
        }
        upper_e.back().emplace_back(c_prime[j]);
      }

      // Iterate through the generated solutions, and combine them with items
      // to generate new solutions. The initial solutions are single items
      // solutions (already present at the upper_* arrays). The 'm' variable
      // stores the number of generated solutions (begins at n, the number
      // of single item solutions).
      for (I j = 2; j <= m; ++j) {
        if (++l == chunk_size) {
          l = 0;
          // Clean chunk with index h, that will not be referenced anymore
          vector<I>().swap(upper_e[h]);
          ++h;
        }
        P i = upper_e[h][l];

        // If this solution 'i' isn't the last already generated solution to
        // have the same 'profit value % c1' then we skip it.
        if (upper_c[i] != j) continue;

        //step_2c:
        // Iterate the items, generating new solutions with symmetry pruning.
        // Combine a solution 'i' with the item 'k' (where k is always equal to
        // or lower than the lowest index of an item on solution 'i').
        for (I k = 2; k <= upper_d[i]; ++k) {
          d = upper_g[i] + d_[k];
          y = y_[i] + c[k];
          x = y % c1;

          //step_2d:
          // x != 0: if x == 0 then the solution profit is perfect divisible
          //  by c1, and you will not use that solution (because using
          //  multiple copies of the best item will be better).
          // d <= upper_g[x]: d is 'the weight of the solution -(minus)
          //  the weight of a relaxed solution with the same profit value
          //  (i.e. multiple copies of the best item, with the last copy
          //  possibly fractionary). 'd' will be smaller for lighter solutions
          //  or more efficient solutions. This way we are trying to obtain
          //  the most efficient solution for the periodic capacity value x. If
          //  the most efficient solution for y % c1 = x is upper_g[x] then we
          //  can backtrack that solution and combine with copies of c1 (for
          //  the remaining profit).
          // (d != upper_g[x] || y < y_[x]): if d is different from upper_g[x]
          //  then we can be safe upper_g[x] is bigger than d; if d is equal to
          //  upper_g[x] (they are equally efficient?), we will chose the one
          //  with the smaller profit value (and consequently smaller weight),
          //  that can be combined with more copies of the best item (or
          //  combined with more solutions).
          if (x != 0 && d <= upper_g[x] && (d != upper_g[x] || y < y_[x])) {
            y_[x] = y;
            upper_g[x] = d;
            upper_d[x] = k;

            // the code below is equivalent to: upper_e[m] = x
            if (upper_e.back().size() == chunk_size) {
              upper_e.emplace_back();
              upper_e.back().reserve(chunk_size);
            }
            upper_e.back().emplace_back(x);

            upper_c[x] = ++m;
          }
        }
      }
      // m isn't modified anymore, we can save it
      eip->m = m;
      HBM_STOP_TIMER(eip->dp_time);

      //step_3:
      HBM_START_TIMER();
      vector<W> x_(m + 1, 0);
      // z is initialized with the best possible value (an upper bound,
      // the relaxed solution of the problem rounded down).
      sol.opt = (c1*b)/a1;
      P &z = sol.opt;

      // As c is more used as the array of item profit values than a simple
      // variable, we use c to denote the array, and c_ to denote the variable
      // used to backtrack the found optimal solution.
      P c_ = 0;
      do {
        x = z % c1;
        // If y_[x] is bigger than z, then y_[x] = z + w*c1, with w > 0; (as
        // the two are in the same modulo c1 residue, and y_[x] is greater).
        // As z begins at the upper bound and it's decremented one unit at
        // time, the only way to missing the z + w*c1 value (that's in the same
        // modulo c1 residue), is that z have begun smaller than y_[x] and,
        // therefore, y_[x] stores a solution over the upper bound (invalid).
        if (y_[x] > z) {
          // The article says: "If x1 < 0, G(x) does not produce F(z)."
          // This means that the periodic solution to the position
          // c1 has a profit bigger than the upper bound. This is possible
          // because the loop that generates new solutions don't care about
          // the solution weight, so invalid solutions (bigger than the
          // capacity) can be generated. This is done by design, the ideia
          // of this algorithm is to compute all optimal solutions (without
          // the best item) to the modulo residue c1. These solutions can be
          // smaller than c1 or multiple times bigger, but only the most
          // efficient solution with the same 'solution profit % c1' (i.e.
          // module residue) will be saved. After we have all the optimal
          // solutions for the c1 module residue, we can guess profit values
          // (beggining with the upper bound and descending). We get the
          // module residue corresponding to the guess profit value, and
          // combine it with copies of the best item until filling
          // the capacity. The problem is that the optimal solution for an
          // specific module residue 'x' can be only reached for capacities
          // bigger than 'b', and therefore we can't solve the problem
          // combining copies of the best item with a periodic solution
          // (optimal solution for a specific module residue value).
          // To solve an instance that ends on this branch of the algorithm
          // we would need to use another algorithm that don't assumes
          // that the c1 periodicity is reached within the 'b' capacity value.
          cout << "negative x1: G(x) does not produce F(z)" << endl;
          cout << "can't find a solution" << endl;
          cout << "see comment on this 'if' for a detailed explanation" << endl;
          HBM_STOP_TIMER(eip->sol_time);
          #ifdef HBM_PROFILE
          eip->total_time = difftime_between_now_and(all_greendp1_begin);
          #endif
          z = 0;
          return;
        }
        // We fill any remaining profit with copies of the best item.
        x_[1] = (z - y_[x])/c1;

        // frac_a1_c1*z == a1*(z/c1) == how many would weight this solution
        // if it was made entirely of copies of the best item (and allowing
        // to fraction the best item). upper_g[x] == the difference between
        // the real weight of the solution, and the weight if it was made of
        // copies of the best item (and allowing to fraction it).
        // frac_a1_c1*z + upper_g == the true weight of the solution
        // Here we are simply cheking if the weight of the current solution
        // is smaller than the capacity (i.e. if the current solution is
        // valid).
        if (frac_a1_c1*z + upper_g[x] <= b) {
          // If the solution is valid, begin the backtrack procedure to
          // assemble the solution.
          c_ = z;
          break;
        } else {
          // If the solution isn't valid, inspect the next one (the one with
          // profit one unity smaller).
          z = z - 1;
          // If no valid solution with bigger profit than a solution composed
          // entirely of copies the best item (not allowing fractions) was
          // found, then the optimal solution is composed entirely of copies
          // of the best item.
          if (z <= c1*(b/a1)) {
            z = c1*(b/a1);
            sol.used_items.emplace_back(items[0], b/a1, 1);
            return;
          }
        }
      } while (true);

      while (c_ != 0) {
        // Backtrack procedure to assemble the subset of the found optimal
        // solution that isn't composed of copies of the best item.
        // Quantity of the best item already included before.
        I k = upper_d[x];
        ++x_[k];
        c_ = y_[x] - c[k];
        x = c_ % c1;
      }

      if (gcd_c > 1) for (P &cj: c) cj *= gcd_c;
      z *= gcd_c;
      //stop:
      // Put the optimal solution on our format.
      for (I k = 1; k <= eip->n2; ++k) {
        if (x_[k] > 0) {
          sol.used_items.emplace_back(item_t<W, P>(a[k], c[k]), x_[k], k);
          sol.y_opt += x_[k]*a[k];
        }
      }
      HBM_STOP_TIMER(eip->sol_time);
      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_greendp1_begin);
      #endif
    }

    /// The second algorithm presented at "On Equivalent Knapsack Problems",
    /// (H.  Greenberg).
    template<typename W, typename P, typename I>
    void greendp2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      greendp2_extra_info_t<W, P, I>* eip =
        new greendp2_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_greendp1_begin = steady_clock::now();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      // I tried to implement the algorithm in a way that anyone can verify
      // this is the same algorithm described at the paper. The goto construct
      // is used because of this. The only ommited goto's are the ones that
      // jump to the next step (that are plainly and clearly unecessary here).
      // BEFORE STEP 1 (notation already introduced by article)
      HBM_START_TIMER();
      auto &items = ukpi.items;
      // In this paper the items are sorted by non-decreasing efficiency
      // and the ties are solved by ordering by non-increasing weight
      // (or profit).
      if (!already_sorted) sort_by_eff(items);
      HBM_STOP_TIMER(eip->sort_time);

      HBM_START_TIMER();
      // CHANGING DATA STRUCTURES TO HAVE A NOTATION SIMILAR TO THE ARTICLE
      // In this algorithm 'n' can change, as we can remove items of the
      // item list.
      const I n = eip->n = static_cast<I>(items.size());
      const W b = eip->c = ukpi.c;

      vector<W> a(n + 1);
      a[0] = 0; // Does not exist, notation begins at 1
      vector<P> c(n + 1);
      c[0] = 0; // Does not exist, notation begins at 1

      for (I i = 0; i < n; ++i) {
        item_t<W, P> it = items[i];
        a[i+1] = it.w;
        c[i+1] = it.p;
      }
      const P gcd_c = gcd_n<typename vector<P>::iterator, P>
                           (c.begin() + 1, c.end());
      if (gcd_c > 1) for (P &cj: c) cj /= gcd_c;
      eip->gcd_c = gcd_c;
      const W a1 = a[1];
      const P c1 = c[1];
      HBM_STOP_TIMER(eip->linear_comp_time);

      HBM_START_TIMER();
      // PARAGRAPH ALGORITHM 2, before step 1
      myvector< rational<P> > d_;
      d_.resize(n + 1);
      const rational<P> frac_a1_c1 = rational<P>(a1, c1);
      for (I j = 2; j <= n; ++j) {
        // On the article this line of the algorithm includes divisions that
        // aren't rounded down. This would call for floating pointer
        // arithmetic. However, those values are compared for equality later.
        // On the article's example, the numbers are shown as fractions, what
        // implies that absolute precision can be necessary. With this in mind,
        // I choose to represent the values of d_ as fractions.
        d_[j] = a[j] - (frac_a1_c1*c[j]);
        // this value is always non-negative (proof follows):
        // previous knowledge: c1/a1 >= cj/aj
        // we can derive:
        // c1 >= (cj/aj)*a1 (multiply by a1)
        // c1*aj >= cj*a1   (multiply by aj, basic efficiency comparation)
        // c1*aj/cj >= a1   (divide by cj)
        // aj/cj >= a1/c1   (divide by c1, inverse of the first statement)
        // aj >= (a1/c1)*cj (aj - (a1/c1)*cj is always non-negative)
      }

      // I reordered the tests to make them clear.
      vector<bool> items_to_remove(n + 1, false);
      I qt_items_to_remove = 0;
      for (I j = 2; j <= n; ++j) {
        if (c[j] % c1 == 0) {
          items_to_remove[j] = true;
          ++qt_items_to_remove;
          continue;
        }
        // The step one of the article don't define the range of 'i'
        // (an item index variable). I do believe that i=1 is excluded
        // because otherwise d_[1] would be referenced (and it's undefined).
        // I also believes that 'i' == 'j' can be skipped, as the conditions
        // for exclusion only evaluate to true if the items are different.
        for (I i = 2; i <= n; ++i) {
          if (i == j) continue;
          if (c[i] % c1 == c[j] % c1) {
            if (d_[i] < d_[j]) {
              items_to_remove[j] = true;
              ++qt_items_to_remove;
              break;
            } else if (d_[i] == d_[j] && c[i] <= c[j]) {
              items_to_remove[j] = true;
              ++qt_items_to_remove;
              break;
            }
          }
        }
      }

      eip->n2 = n - qt_items_to_remove;
      // The number of items that weren't excluded by step one.
      I m;
      // Scope to guarantee vector deallocation.
      {
        m = n - qt_items_to_remove;

        vector<W> a_;
        a_.reserve(m + 1);
        a_.push_back(0);

        vector<P> c_;
        c_.reserve(m + 1);
        c_.push_back(0);

        myvector< rational<P> > tmp_d_;
        tmp_d_.reserve(m + 1);
        tmp_d_.push_back(0);

        for (I j = 1; j <= n; ++j) {
          if (!items_to_remove[j]) {
            a_.push_back(a[j]);
            c_.push_back(c[j]);
            tmp_d_.push_back(d_[j]);
          }
        }

        a_.swap(a);
        c_.swap(c);
        tmp_d_.swap(d_);
      }
      HBM_STOP_TIMER(eip->dom_time);

      HBM_START_TIMER();
      // We need to declare the variables before we start jumping
      P c_, i, x, y;
      I j, k;
      rational<P> d;
      //step_1: // unused step
      HBM_START_TIMER();
      myvector<P> c_prime; // We will call c' as c_prime
      c_prime.resize(m + 1);
      // Every uppercase letter will be upper_x
      //map<P, rational<P> > upper_g;
      //map<P, I> upper_c;
      //map<P, I> upper_d;
      //map<I, P> upper_e;

      // The x variable used to index all arrays but upper_e have (c1 - 1) as
      // an upper bound. So the arrays upper bound size is c1.

      // As there's an array 'y' and a variable 'y', we will call
      // the array 'y_'
      vector<P> y_(c1, 0);
      // Every uppercase letter will be upper_x
      vector<rational<P> > upper_g(c1, 0);
      vector<I> upper_c(c1, 0);
      vector<I> upper_d(c1, 0);
      map<I, P> upper_e;
      //vector<P> upper_e;
      //upper_e.push_back(0); // positiom 0
      //upper_e.push_back(0); // position 1
      HBM_STOP_TIMER(eip->vector_alloc_time);

      HBM_START_TIMER();
      const rational<P> infinity(numeric_limits<P>::max(), 1);
      for (P k = 1; k < c1; ++k) {
        // The positions where k == c_prime[j] for j=2..m will be overwritten
        upper_g[k] = infinity;
      }

      for (I j = 2; j <= m; ++j) {
        c_prime[j] = c[j] - (c[j]/c1)*c1;
        upper_g[c_prime[j]] = d_[j];
        y_[c_prime[j]] = c[j];
        upper_c[c_prime[j]] = j;
        upper_d[c_prime[j]] = j;
        upper_e[j] = c_prime[j];
      }

      j = 1;

      step_2a:
      j = j + 1;
      if (j > m) goto step_3;

      //step_2b: // unused label
      i = upper_e[j];
      if (upper_c[i] != j) goto step_2a;
      k = 1;

      step_2c:
      k = k + 1;
      if (k > upper_d[i]) goto step_2a;
      d = upper_g[i] + d_[k];
      y = y_[i] + c[k];
      x = y - (y/c1)*c1;

      //step_2d: //unused label
      if (x == 0 || d > upper_g[x] || (d == upper_g[x] && y >= y_[x])) {
        //goto step_2c;
      } else {
        m = m + 1;
        upper_c[x] = m;
        //upper_e.push_back(x);
        upper_e[m] = x;
        upper_g[x] = d;
        upper_d[x] = k;
        y_[x] = y;
        //goto step_2c;
      }

      goto step_2c;
      // Because the goto above we only get on step_3 jumping directly to it.

      step_3:
      HBM_STOP_TIMER(eip->dp_time);
      HBM_START_TIMER();
      // We initialize the vector here, as x[1] is already used on the next
      // step.
      vector<W> x_(m + 1, 0);
      sol.opt = (c1*b)/a1;
      P &z = sol.opt;
      // m isn't modified anymore, we can save it
      eip->m = m;

      step_3a:
      x = z - (z/c1)*c1;
      if (z < y_[x]) {
        // The article says: "If x1 < 0, G(x) does not produce F(z)."
        // This means that the method failed and can't possibly
        // find a solution?
        cout << "negative x1: G(x) does not produce F(z)" << endl;
        cout << "can't find a solution" << endl;
        goto stop;
      }
      x_[1] = (z - y_[x])/c1;

      //step_3b: //unused label
      // As c is more used as the array of item profit values than a simple
      // variable, we use c to denote the array, and c_ to denote the variable.
      if (b < frac_a1_c1*z + upper_g[x]) {
        goto step_3c;
      } else {
        c_ = z;
        // The article asks to set x[j] = 0, j=2..n, on this step, we already
        // did this before.
        goto step_4a;
      }

      step_3c:
      z = z - 1;
      if (z > c1*(b/a1)) goto step_3a;
      else {
        z = c1*(b/a1);
        x_[1] = (b/a1);
        cout << "stopped as filled only by the best item" << endl;
        // The algorithm don't give the next step here.
        goto stop;
      }

      step_4a:
      k = upper_d[x];
      x_[k] = x_[k] + 1;
      c_ = y_[x] - c[k];

      //step_4b: // unused label
      if (c_ > 0) {
        x = c_ - (c_/c1)*c1;
        goto step_4a;
      }

      stop:
      if (gcd_c > 1) for (P &cj: c) cj *= gcd_c;
      z *= gcd_c;
      // Put the optimal solution on our format.
      for (I k = 1; k <= m; ++k) {
        if (x_[k] > 0) {
          sol.used_items.emplace_back(item_t<W, P>(a[k], c[k]), x_[k], k);
          sol.y_opt += x_[k]*a[k];
        }
      }

      HBM_STOP_TIMER(eip->sol_time);
      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_greendp1_begin);
      #endif
    }

    /// The modernized version of the first algorithm presented at "On
    /// Equivalent Knapsack Problems", (H.  Greenberg).
    template<typename W, typename P, typename I>
    void mgreendp1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      mgreendp1_extra_info_t<W, P, I>* eip =
        new mgreendp1_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_greendp1_begin = steady_clock::now();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      // BEFORE STEP 1 (notation already introduced by article)
      HBM_START_TIMER();
      auto &items = ukpi.items;
      // In this paper the items are sorted by non-decreasing efficiency
      // and the ties are solved by ordering by non-increasing weight
      // (or profit).
      if (!already_sorted) sort_by_eff(items);
      HBM_STOP_TIMER(eip->sort_time);

      HBM_START_TIMER();
      const I n = eip->n = static_cast<I>(items.size());
      // We use b to represent the capacity, as in the article, because
      // there's already three distinct 'c' variables on the algorithm
      // (as described in the article).
      const W b = eip->c = ukpi.c;

      // Items weight array.
      vector<W> a(n + 1);
      a[0] = 0; // Does not exist, notation begins at 1
      // Items profit array.
      vector<P> c(n + 1);
      c[0] = 0; // Does not exist, notation begins at 1

      for (I i = 0; i < n; ++i) {
        item_t<W, P> it = items[i];
        a[i+1] = it.w;
        c[i+1] = it.p;
      }
      const P gcd_c = gcd_n<typename vector<P>::iterator, P>
                           (c.begin() + 1, c.end());
      if (gcd_c > 1) for (P &cj: c) cj /= gcd_c;
      eip->gcd_c = gcd_c;
      const W a1 = a[1];
      const P c1 = c[1];

      // The constant 't' is a value strictly greater than the optimal solution
      // value. It's used to multiply the items weight, this way we can codify
      // both profit and weight in a single number (p + t*w = x, where x is a
      // combined number). To extract them back use: p = x % t; w = x / t. We
      // will refer to these numbers (that represents both a weight and a
      // profit) as combined numbers. The items are represented this way, as
      // the solutions (that are simply the sum of many items).
      const P t = (c1*b)/a1 + 1;
      // the real upper bound t value is dependent on gdc_c
      eip->t = (c1*gcd_c*b)/a1 + 1;

      // z: A solution profit value. The algorithm will end with the
      //    optimal solution value inside this variable.
      P &z = sol.opt;
      HBM_STOP_TIMER(eip->linear_comp_time);

      HBM_START_TIMER();
      // STEP 1
      // The upper_f maps a profit value to a combined number, initially this
      // means that giving an item profit will return the combined number of
      // the item. After, this means that giving z will return you the
      // associated solution (sum of the combined numbers that make that
      // solution).
      vector<P> upper_f(t, 0);
      // The algorithm have an lowercase d, an uppercase D, and a d with
      // subscript, we use "d_" to denote the d with subscript.
      // The d_ array maps the index of an item to its combined number
      // representation.
      myvector<P> d_;
      d_.resize(n + 1);
      d_[0] = 0; // Does not exist, notation begins at 1

      // The upper_c and upper_e are inverted indexes of each other.
      // upper_c maps profit values to indexes, and upper_e maps
      // indexes to profit values. The problem is: the indexes don't
      // are item indexes, they don't end at n, they end at 'm'.
      // The 'm' variable is the number of solutions generated by the
      // algorithm. So the indexes aren't only item indexes, but are
      // solution indexes.
      vector<I> upper_c(t, 0);

      // For upper_e we need a structure between a vector and a list. We want
      // fast access to one point inside the container and to its end (could be
      // done by list or vector). We want to clean from memory the elements
      // already iterated by our iterator in a O(1) way (list). We want to make
      // better use of cache and do not allocate memory at every new item
      // inserted at the container's end (vector). So the ideia is a vector of
      // vectors, the inner vectors have a max size, and we can ask to clean an
      // entire inner vector (i.e. chunk) easily. I let the vector and
      // forward_list versions commented to allow alternating between them. The
      // forward_list version has the smallest memory footprint but don't have
      // a very good performance. The vector version has slight better time
      // than the chunk-based version, but it will consume many and many GiB
      // easily, where the chunk based memory footprint can be tweaked by
      // chunk_size and have a good performance with a very small memory
      // footprint. The vector version performance also depends on the
      // re-allocation strategy, and the initial reserved amount.
      // If the P value is 64 bits, a 1000000 chunk_size is about 7MiB.
      // TODO: make this a configurable parameter
      // chunk-based upper_e code below
      const size_t chunk_size = 1000000;
      size_t h = 0, l = 0;
      vector< vector<P> > upper_e;
      upper_e.emplace_back();
      upper_e.front().reserve(chunk_size);
      upper_e.front().emplace_back(0);

      // forward_list upper_e code below
      //forward_list<P> upper_e;
      //auto last_z_value = upper_e.before_begin();
      //last_z_value = upper_e.insert_after(last_z_value, 0);

      // vector upper_e code below
      //vector<P> upper_e;
      //upper_e.push_back(0);

      // This array (upper_d) is used to backtrack the optimal solution.
      // It maps solution profits (z) to the index of the last item used in
      // that solution.
      vector<I> upper_d(t, 0);
      HBM_STOP_TIMER(eip->vector_alloc_time);

      HBM_START_TIMER();
      for (I j = 1; j <= n; ++j) {
        // inverted indexes initialization
        upper_c[c[j]] = j;
        if (upper_e.back().size() == chunk_size) {
          upper_e.emplace_back();
          upper_e.back().reserve(chunk_size);
        }
        upper_e.back().emplace_back(c[j]);

        // forward_list upper_e code below
        //last_z_value = upper_e.insert_after(last_z_value, c[j]);

        // vector code below
        //upper_e.push_back(c[j]);

        // Combined numbers for items (d_), and upper_f is initialized with
        // them (in upper_f they are single item solutions).
        d_[j] = c[j] + t*a[j];
        upper_f[c[j]] = d_[j];

        // This array (upper_d) is initialized like upper_c, but after upper_c
        // will receive m values (index of solutions), and upper_d receives k
        // values (index of items).
        upper_d[c[j]] = j;
      }
      // The index of the current solution.
      W j = 0;
      // The number of solutions (optimal or not) generated by the algorithm.
      W m = n;

      // The j is the index of the current solution (initially we will consider
      // every single item as an single item solution). The m is the total
      // number of solutions generated. The m value will increase while the
      // algorithm executes, as we will iterate over the existing solutions
      // with k and generate new solutions combining the current solution k
      // with the items. When j reaches m we have generated enough solutions
      // relevant to find one optimal solution between them, and we can stop.
      while (++j <= m) {
        // chunk-based upper_e code below
        if (++l == chunk_size) {
          l = 0;
          // Clean chunk with index h, that will not be referenced anymore
          vector<I>().swap(upper_e[h]);
          ++h;
        }
        P z_curr_sol = upper_e[h][l];

        // forward_list upper_e code below
        //if (upper_e.empty()) break;
        //P z_curr_sol = upper_e.front();
        //upper_e.pop_front();

        // vector upper_e code below
        //P z_curr_sol = upper_e[j];

        // Skips a solution if there's already another solution with the same
        // profit value and a greater position index (j). Don't make sense to
        // create new solutions from two solutions with the same profit value,
        // we will process only one solution for each solution tied with the
        // same profit.
        if (upper_c[z_curr_sol] != j) continue;

        // This loop combines existing solutions with new items generating more
        // solutions.
        // As upper_d is like 'd' in ukp5, and stores the index of the last
        // item used in a solution, iterating until upper_d is the same as
        // pruning symmetric solutions? Seems so.
        P k = 0;
        while (++k <= upper_d[z_curr_sol]) {
          // Combines the solution with z == i with the item k, and
          // extracts z back from the new solution.
          // d is the combined number for the newly generated solution
          P d = upper_f[z_curr_sol] + d_[k];
          z = d - (d/t)*t;

          // We will only save the newly generated solution d, if its profit
          // (z) is bigger than zero and there isn't a solution with that
          // profit saved at upper_f, or the solution with the same profit
          // already saved at upper_f has more weight than ours (it can be
          // invalid, or simply a solution with more weight for the same
          // profit).
          if (z != 0 && (upper_f[z] == 0 || upper_f[z] > d)) {
            // increment the number of solutions by one
            m = m + 1;
            // chunk-based upper_e code below
            if (upper_e.back().size() == chunk_size) {
              upper_e.emplace_back();
              upper_e.back().reserve(chunk_size);
            }
            upper_e.back().emplace_back(z);

            // forward_list upper_e code below
            //last_z_value = upper_e.insert_after(last_z_value, z);

            // vector upper_e code below
            //upper_e.push_back(z);

            upper_c[z] = m;
            upper_f[z] = d;
            // THIS LINE WAS MISSING FROM THE ORIGINAL ARTICLE
            upper_d[z] = k;
          }
        }
      }
      // m isn't modified anymore, we can save it
      eip->m = m;
      HBM_STOP_TIMER(eip->dp_time);

      // We start with the greatest possible value for the optimal
      // solution, and then we decrement until finding a valid solution
      // with that profit value.
      HBM_START_TIMER();
      z = t - 1;
      // The values for the decision variables. If x[k] == y, then
      // there's y copies of the item k on an optimal solution.
      // The x[n+1] stores how much empty space there's in an optimal
      // solution (the difference between the weight of the optimal
      // solution we found and the capacity of the instance knapsack).
      vector<W> x(n + 2, 0); // It's indexed until x + 1

      bool only_best_item_used = false;
      // If upper_f[z] is greater than z + t*b, then the combined number
      // returned by upper_f[z] isn't a valid solution. The reason for that
      // is simple: the profit value of upper_f[z] is equal to z
      // (z + t*b < z + t*b' == upper_f[z]), so if upper_f[z] is bigger
      // than z + t*b, is because b' is bigger than b, and then the solution
      // is heavier than the knapsack capacity (i.e. an invalid solution).
      while (upper_f[z] > z + t*b) {
        z = z - 1;
        // If we end up with z being perfectly divisible by the best item
        // profit, then it's clear that filling the knapsack with copies
        // of the best item is an optimal solution.
        if (z == c1*(b/a1)) {
          x[1] = b/a1;
          x[n + 1] = b - a1*(b/a1);
          only_best_item_used = true;
          break;
        }
      }

      if (!only_best_item_used) {
        // Store at x[n + 1] how much of the the capacity was left unfilled
        // by the optimal solution.
        x[n + 1] = (z + t*b - upper_f[z])/t;

        // As c is more used as the array of item profit values than a simple
        // variable, we use c to denote the array, and c_ to denote the
        // variable.
        P c_ = z;

        // c_ begins with the optimal solution profit and then it's decremented
        // for every item we put on the x array. It's a typical backtrack for
        // assembling the solution (only that this one is by the profit, not
        // by the weight).
        do {
          I k = upper_d[c_];
          x[k] = x[k] + 1;
          // This part is a little different from greendp1 because we are
          // avoiding underflow (the P type used by c_ can be unsigned).
          // I don't know in what circunstances c_ - c[k] is expected to
          // be negative, but the original algorithm considered this a
          // possibility.
          if (c_ < c[k]) {
            c_ = c_ + t - c[k];
          } else {
            c_ = c_ - c[k];
          }
        } while (c_ != 0);
      }

      z *= gcd_c;
      // Put the optimal solution on our format.
      for (I k = 1; k <= n; ++k) {
        if (x[k] > 0) {
          sol.used_items.emplace_back(items[k-1], x[k], k);
          sol.y_opt += x[k]*items[k-1].w;
        }
      }
      HBM_STOP_TIMER(eip->sol_time);

      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_greendp1_begin);
      #endif
      return;
    }

    /// The first algorithm presented at "On Equivalent Knapsack Problems",
    /// (H. Greenberg). See mgreendp1 for a explained version.
    template<typename W, typename P, typename I>
    void greendp1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      greendp1_extra_info_t<W, P, I>* eip =
        new greendp1_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_greendp1_begin = steady_clock::now();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      // I tried to implement the algorithm in a way that anyone can verify
      // this is the same algorithm described at the paper. The goto construct
      // is used because of this. The only ommited goto's are the ones that
      // jump to the next step (that are plainly and clearly unecessary here).
      // BEFORE STEP 1 (notation already introduced by article)
      HBM_START_TIMER();
      auto &items = ukpi.items;
      // In this paper the items are sorted by non-decreasing efficiency
      // and the ties are solved by ordering by non-increasing weight
      // (or profit).
      if (!already_sorted) sort_by_eff(items);
      HBM_STOP_TIMER(eip->sort_time);

      HBM_START_TIMER();
      // CHANGING DATA STRUCTURES TO HAVE A NOTATION SIMILAR TO THE ARTICLE
      const I n = eip->n = static_cast<I>(items.size());
      const W b = eip->c = ukpi.c;

      vector<W> a(n + 1);
      a[0] = 0; // Does not exist, notation begins at 1
      vector<P> c(n + 1);
      c[0] = 0; // Does not exist, notation begins at 1

      for (I i = 0; i < n; ++i) {
        item_t<W, P> it = items[i];
        a[i+1] = it.w;
        c[i+1] = it.p;
      }
      const P gcd_c = gcd_n<typename vector<P>::iterator, P>
                           (c.begin() + 1, c.end());
      if (gcd_c > 1) for (P &cj: c) cj /= gcd_c;
      eip->gcd_c = gcd_c;
      const W a1 = a[1];
      const P c1 = c[1];

      P t = (c1*b)/a1 + 1;
      // the real upper bound t value is dependent on gdc_c
      eip->t = (c1*gcd_c*b)/a1 + 1;

      P i, k, d, &z = sol.opt;
      HBM_STOP_TIMER(eip->linear_comp_time);

      HBM_START_TIMER();
      // STEP 1
      vector<P> upper_f(t, 0);
      vector<I> upper_c(t, 0);
      // The algorithms have an lowercase d, an uppercase D, and a d with
      // subscript, we use "d_" to denote the d with subscript.
      myvector<P> d_;
      d_.resize(n + 1);
      d_[0] = 0; // Does not exist, notation begins at 1
      vector<I> upper_d(t, 0);
      // Here we use a map for upper_e, this is terrible inneficient,
      // but is the structure with the closest notation to the one used
      // on the article.
      map<P, I> upper_e;
      for (I j = 1; j <= n; ++j) {
        d_[j] = c[j] + t*a[j];
        upper_c[c[j]] = j;
        upper_d[c[j]] = j;
        upper_e[j] = c[j];
        upper_f[c[j]] = d_[j];
      }
      vector<W> x(n + 2, 0); // It's indexed until x + 1
      W j = 0;
      W m = n;
      HBM_STOP_TIMER(eip->vector_alloc_time);

      HBM_START_TIMER();
      step_2a:
      j = j + 1;
      if (j > m) goto step_3;

      //step_2b: // unused label
      //i = mfwd(upper_e, j, static_cast<size_t>(0));
      i = upper_e[j];
      if (upper_c[i] != j) {
        goto step_2a;
      } else {
        k = 0;
      }

      step_2c:
      k = k + 1;
      if (k > upper_d[i]) {
        goto step_2a;
      } else {
        if (k > n) cout << "k: " << endl;
        d = upper_f[i] + d_[k];
        z = d - (d/t)*t;
      }

      //step_2d: // unused label
      if (z == 0 || (0 < upper_f[z] && upper_f[z] <= d)) {
        //goto step_2c;
      } else {
        m = m + 1;
        //if (upper_c.count(z)) cout << "old m: " << upper_c[z] << endl << "new m: " << m << endl;
        upper_c[z] = m;
        upper_e[m] = z;
        upper_f[z] = d;
        // THIS LINE WAS MISSING FROM THE ARTICLE
        upper_d[z] = k;
        //goto step_2c;
      }
      goto step_2c;

      // Because the goto above, the code below can only be accessed by jumping
      // directly into step_3.
      step_3:
      // As m isn't modified anymore, we save it here
      eip->m = m;
      HBM_STOP_TIMER(eip->dp_time);
      HBM_START_TIMER();
      z = t - 1;
      // As c is more used as the array of item profit values than a simple
      // variable, we use c to denote the array, and c_ to denote the variable.
      P c_;
      step_3a:
      if (z + t*b < upper_f[z]) {
        goto step_3b;
      } else {
        x[n + 1] = (z + t*b - upper_f[z])/t;
        c_ = z;
        goto step_3c;
      }

      step_3b:
      z = z - 1;
      if (z == c1*(b/a1)) {
        goto step_3e;
      } else {
        goto step_3a;
      }

      step_3c:
      k = upper_d[c_];
      x[k] = x[k] + 1;

      //step_3d: // unused label
      c_ = c_ - c[k];
      if (c_ < 0) {
        c_ = c_ + t;
        goto step_3c;
      } else if (c_ > 0) {
        goto step_3c;
      } else {
        goto stop;
      }

      step_3e:
      x[1] = b/a1;
      x[n + 1] = b - a1*(b/a1);

      stop:
      z *= gcd_c;
      // Put the optimal solution on our format.
      for (I l = 1; l <= n; ++l) {
        if (x[l] > 0) {
          sol.used_items.emplace_back(items[l-1], x[l], l);
          sol.y_opt += x[l]*items[l-1].w;
        }
      }
      HBM_STOP_TIMER(eip->sol_time);

      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_greendp1_begin);
      #endif
      return;
    }

    // The dynamic programming algorithm presented at "A Better Step-off
    // Algorithm for the Knapsack Problem" (H. Greenberg).
    template<typename W, typename P, typename I>
    void greendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      greendp_extra_info_t<W, P, I>* eip = new greendp_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_greendp_begin = steady_clock::now();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      // I tried to implement the algorithm in a way that anyone can verify
      // this is the same algorithm described at the paper. The goto construct
      // is used because of this. The only ommited goto's are the ones that
      // jump to the next step (that are plainly and clearly unecessary here).
      // BEFORE STEP 1 (on the bold 'Algorithm.' block)
      HBM_START_TIMER();
      auto &items = ukpi.items;
      if (!already_sorted) {
        sort_by_eff(items);
        reverse(items.begin(), items.end());
      }
      HBM_STOP_TIMER(eip->sort_time);

      HBM_START_TIMER();
      // CHANGING DATA STRUCTURES TO HAVE A NOTATION SIMILAR TO THE ARTICLE
      const I m = eip->n = static_cast<I>(items.size());
      const W b = eip->c = ukpi.c;

      vector<W> a(m + 1);
      a[0] = 0; // Does not exist, notation begins at 1
      vector<P> c(m + 1);
      c[0] = 0; // Does not exist, notation begins at 1
      // As b is the capacity, b_ will be used where b was subscribed
      vector<P> b_(m + 1);
      b_[0] = 0; // Does not exist, notation begins at 1

      const W am = items[m - 1].w;
      const P cm = items[m - 1].p;
      for (I i = 0; i < m; ++i) {
        item_t<W, P> it = items[i];
        a[i+1] = it.w;
        c[i+1] = it.p;
      }

      // Starts at index 1, this is the cause of the +1
      const W d = gcd_n<typename vector<P>::iterator, P>(c.begin() + 1, c.end());

      // lambda is the remaining space if we fill the capacity b
      // with the most efficient item.
      const W lambda = b - (b / a[m]) * a[m];
      W upper_l = lambda;
      HBM_STOP_TIMER(eip->linear_comp_time);

      HBM_START_TIMER();
      vector<P> f(b + 1, 0);
      vector<I> i(b + 1);
      i[0] = 1;
      HBM_STOP_TIMER(eip->vector_alloc_time);

      HBM_START_TIMER();
      P &z = sol.opt; // z is another name for sol.opt
      z = cm * (b / am);
      W k = 0;

      // STEP 1
      for (I i = 1; i < m; ++i) {
        b_[i] = cm*a[i] - c[i]*am;
      }
      HBM_STOP_TIMER(eip->linear_comp_time);

      HBM_START_TIMER();
      W &y = eip->last_y_value_outer_loop = 0;
      W k_max;
      if (!compute_k_max(m, b_, a, c, am, cm, b, z, d, k, lambda, k_max)) {
        HBM_STOP_TIMER(eip->bound_time);
        goto stop;
      }
      HBM_STOP_TIMER(eip->bound_time);

      // STEP 2a
      I j, v;
      step_2a:
      HBM_START_TIMER();
      // if DEBUG
      j = i[y];

      // STEP 2b
      step_2b:
      if (y + a[j] <= lambda + am*k_max)  {
        v = c[j] + f[y];
        goto step_2c;
      } else {
        goto step_2d;
      }

      // STEP 2c
      step_2c:
      if (v >= f[y + a[j]]) {
        f[y + a[j]] = v;
        i[y + a[j]] = j;
      }

      // STEP 2d
      step_2d:
      if (j < m - 1) {
        j = j + 1;
        goto step_2b;
      }
      HBM_STOP_TIMER(eip->dp_time);

      // STEP 3a
      step_3a:
      HBM_START_TIMER();
      y = y + 1;

      // STEP 3b
      if (f[y] > f[y - 1]) {
        HBM_STOP_TIMER(eip->dp_time);
        goto step_3c;
      } else {
        f[y] = f[y - 1];
        i[y] = m + 1;
        HBM_STOP_TIMER(eip->dp_time);
        goto step_3d;
      }

      // STEP 3c
      step_3c:
      HBM_START_TIMER();
      if (y == upper_l) {
        if (!next_z(m, b_, a, c, am, cm, b, z, d, k, lambda, k_max, f, y, upper_l)) {
          HBM_STOP_TIMER(eip->bound_time);
          goto stop;
        }
      }
      HBM_STOP_TIMER(eip->bound_time);
      goto step_2a;
      // Because the goto above, the code below can only be accessed by jumping
      // directly into step_3d

      // STEP 3d
      step_3d:
      HBM_START_TIMER();
      if (y == upper_l) {
        if (!next_z(m, b_, a, c, am, cm, b, z, d, k, lambda, k_max, f, y, upper_l)) {
          HBM_STOP_TIMER(eip->bound_time);
          goto stop;
        }
      }
      HBM_STOP_TIMER(eip->bound_time);
      goto step_3a;
      // Because the goto above, the code below can only be accessed by jumping
      // directly into stop
      stop:
      sol.y_opt = y;

      HBM_START_TIMER();
      vector<I> qts_its(m, 0);

      W y_opt = lambda;
      W bi_qt = b/am;
      for (; y_opt < b && f[y_opt] != z - bi_qt*cm; y_opt += am, --bi_qt);
      while (y_opt > 0 && f[y_opt-1] == f[y_opt]) --y_opt;
      sol.y_opt = y_opt + bi_qt*am;

      I dy_opt;
      while (y_opt != 0) {
        dy_opt = i[y_opt];
        y_opt -= a[dy_opt];
        ++qts_its[dy_opt];
      }

      for (I x = 1; x < m; ++x) {
        if (qts_its[x] > 0) {
          sol.used_items.emplace_back(items[x-1], qts_its[x], x);
        }
      }

      if (bi_qt > 0) {
        auto bqt = itemqt_t<W, P, I>(items[m-1], bi_qt, m);
        sol.used_items.push_back(bqt);
      }

      sol.used_items.shrink_to_fit();
      HBM_STOP_TIMER(eip->sol_time);
      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_greendp_begin);
      #endif
    }

    template<typename W, typename P, typename I>
    struct mgreendp1_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_greendp_impl::mgreendp1(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "mgreendp1";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void mgreendp1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol,
                int argc, argv_t argv) {
      simple_wrapper(mgreendp1_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }

    template<typename W, typename P, typename I>
    struct mgreendp_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_greendp_impl::mgreendp(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "mgreendp";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void mgreendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol,
                int argc, argv_t argv) {
      simple_wrapper(mgreendp_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }

    template<typename W, typename P, typename I>
    struct mgreendp2_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_greendp_impl::mgreendp2(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "mgreendp2";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void mgreendp2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol,
                  int argc, argv_t argv) {
      simple_wrapper(mgreendp2_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }

    template<typename W, typename P, typename I>
    struct greendp2_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_greendp_impl::greendp2(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "greendp2";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void greendp2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol,
                  int argc, argv_t argv) {
      simple_wrapper(greendp2_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }

    template<typename W, typename P, typename I>
    struct greendp1_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_greendp_impl::greendp1(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "greendp1";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void greendp1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol,
                  int argc, argv_t argv) {
      simple_wrapper(greendp1_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }

    template<typename W, typename P, typename I>
    struct greendp_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_greendp_impl::greendp(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "greendp";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void greendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol,
                int argc, argv_t argv) {
      simple_wrapper(greendp_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }
  }

  /// A rewritten version of the second algorithm presented at "On Equivalent
  /// Knapsack Problems" (H. Greenberg), an attemp to adapt the algorithm to
  /// structured programming.
  ///
  /// @see greendp2
  template<typename W, typename P, typename I>
  void mgreendp2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_greendp_impl::mgreendp2(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items isn't sorted by non-decreasing weight. If it's
  /// ommited the ukpi.items is sorted by non-decreasing weight.
  ///
  /// @see main_take_path
  /// @see mgreendp2(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void mgreendp2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_greendp_impl::mgreendp2(ukpi, sol, argc, argv);
  }

  /// Solves an UKP instance by one of the dynamic programming algorithm
  /// presented at "On Equivalent Knapsack Problems" from H. Greenberg (this is
  /// the one presented as "Algorithm 2"), and stores the results at sol.
  ///
  /// @note IMPORTANT: only works with integers, as it relies on the rounding
  ///   behaviour of integers on divisions.
  /// @param ukpi The UKP instance to be solved.
  /// @param sol The object where the results will be written.
  /// @param already_sorted If the ukpi.items vector needs to be sorted by
  ///   non-decreasing efficiency. An if there is more than one most efficient
  ///   item, the one with the smallest profit (or weight) should be the first.
  template<typename W, typename P, typename I>
  void greendp2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_greendp_impl::greendp2(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items isn't sorted by non-decreasing weight. If it's
  /// ommited the ukpi.items is sorted by non-decreasing weight.
  ///
  /// @see main_take_path
  /// @see greendp2(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void greendp2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_greendp_impl::greendp2(ukpi, sol, argc, argv);
  }

  /// Solves an UKP instance by one of the dynamic programming algorithm
  /// presented at "On Equivalent Knapsack Problems" from H. Greenberg (this is
  /// the one presented as "Algorithm 1"), and stores the results at sol.
  ///
  /// @note IMPORTANT: only works with integers, as it relies on the rounding
  ///   behaviour of integers on divisions.
  /// @param ukpi The UKP instance to be solved.
  /// @param sol The object where the results will be written.
  /// @param already_sorted If the ukpi.items vector needs to be sorted by
  ///   non-decreasing efficiency. An if there is more than one most efficient
  ///   item, the one with the smallest profit (or weight) should be the first.
  template<typename W, typename P, typename I>
  void greendp1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_greendp_impl::greendp1(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items isn't sorted by non-decreasing weight. If it's
  /// ommited the ukpi.items is sorted by non-decreasing weight.
  ///
  /// @see main_take_path
  /// @see greendp1(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void greendp1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_greendp_impl::greendp1(ukpi, sol, argc, argv);
  }

  /// A rewritten version of the first algorithm presented at "On Equivalent
  /// Knapsack Problems" (H. Greenberg), an attemp to adapt the algorithm to
  /// structured programming.
  ///
  /// @see greendp1
  template<typename W, typename P, typename I>
  void mgreendp1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_greendp_impl::mgreendp1(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items isn't sorted by non-decreasing weight. If it's
  /// ommited the ukpi.items is sorted by non-decreasing weight.
  ///
  /// @see main_take_path
  /// @see mgreendp1(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void mgreendp1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_greendp_impl::mgreendp1(ukpi, sol, argc, argv);
  }

  /// A rewritten version of the algorithm presented at "A Better Step-off
  /// Algorithm for the Knapsack Problem", an attemp to adapt the algorithm
  /// to structured programming.
  ///
  /// @see greendp
  template<typename W, typename P, typename I>
  void mgreendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_greendp_impl::mgreendp(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items isn't sorted by non-decreasing weight. If it's
  /// ommited the ukpi.items is sorted by non-decreasing weight.
  ///
  /// @see main_take_path
  /// @see mgreendp(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void mgreendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_greendp_impl::mgreendp(ukpi, sol, argc, argv);
  }

  /// Solves an UKP instance by the dynamic programming algorithm presented at
  /// "A Better Step-off Algorithm for the Knapsack Problem", and stores the
  /// results at sol.
  ///
  /// @note IMPORTANT: only works with integers, as it relies on the rounding
  ///   behaviour of integers on divisions.
  /// @param ukpi The UKP instance to be solved.
  /// @param sol The object where the results will be written.
  /// @param already_sorted If the ukpi.items vector needs to be sorted by
  ///   non-decreasing efficiency.
  template<typename W, typename P, typename I>
  void greendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_greendp_impl::greendp(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items isn't sorted by non-decreasing weight. If it's
  /// ommited the ukpi.items is sorted by non-decreasing weight.
  ///
  /// @see main_take_path
  /// @see greendp(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void greendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_greendp_impl::greendp(ukpi, sol, argc, argv);
  }
}

#endif //HBM_GREENDP_HPP

