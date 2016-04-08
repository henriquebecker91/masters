#ifndef HBM_GREENDP_HPP
#define HBM_GREENDP_HPP

#include "ukp_common.hpp"
#include "periodicity.hpp"
#include "wrapper.hpp"

#include <vector> // For vector
#include <algorithm>  // For reverse
#include <boost/math/common_factor_rt.hpp>  // For boost::math::gcd

#define PRINT_VAR(x) cout << #x << ": " << x << endl

namespace hbm {
  namespace hbm_greendp_impl {
    using namespace std;
    using namespace boost::math;

    /// Compute the gcd for a variable quantity of numbers.
    /// @note This functions was designed with positive integers on mind.
    template<typename ForwardIterator, typename T>
    T gcd_n(ForwardIterator begin, ForwardIterator end) {
      if (begin == end) {
        return 1;
      }

      if (begin + 1 == end) {
        return *begin;
      }

      T acc = gcd(*begin, *(begin+1));
      ++begin;

      for (; begin != end && acc > 1; ++begin) {
        acc = gcd(acc, *begin);
      }
      
      return acc;
    }

    /// Compute k_max, also, if returns false the greendp procedure (caller)
    /// is to be stopped.
    template<typename W, typename P, typename I>
    bool compute_k_max( const I m,
                        const vector<P> &b_,
                        const vector<W> &a,
                        const vector<P> &c,
                        const W am,
                        const P cm,
                        const W b,
                        const P z,
                        const P d,
                        const W k,
                        const W lambda,
                        W &k_max) {
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
                          const W c,
                          const P opt,
                          const P d,
                          const W bi_qt,
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

    /// Compute k_max, also, if returns false the greendp procedure (caller)
    /// is to be stopped.
    /*template<typename W, typename P>
    bool compute_k_max( const vector<P> &b,
                        const vector< item_t<W, P> > &items,
                        const W c,
                        const P opt,
                        const P d,
                        const W k,
                        const W lambda,
                        W &k_max) {
      // STEP 1 (Routine k_max)
      const item_t<W, P> &bi = items[0];

      size_t r = 1;
      while (r < items.size() && b[r] > bi.p*c - bi.w*(opt + d)) ++r;

      if (r == items.size()) return false;

      // STEP 2 (Routine k_max)
      const P delta = opt - bi.p*(c/bi.w) + d; 
      const P l_side = (items[r].p*lambda - items[r].w*delta)/b[r];
      const P max_bi_qt = c/bi.w;
      k_max = l_side < max_bi_qt ? l_side : max_bi_qt;

      if (k > k_max) return false;

      return true;
    }*/

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

      //const P bi_g = gcd(bi.w, bi.p);
      //const item_t<W, P> bi_n(bi.w / bi_g, bi.p / bi_g);
      for (size_t i = 1; i < n; ++i) {
        const item_t<W, P> &it = items[i];
        b.push_back(bi.p*it.w - it.p*bi.w);
        //const P it_g = gcd(it.w, it.p);
        //const item_t<W, P> it_n(it.w / it_g, it.p / it_g);
        //cout << "it_g: " << it_g << endl;
        //cout << "b_n[" << i << "]: " << bi_n.p*it_n.w - it_n.p*bi_n.w << endl;
      }
      //cout << "bi_g: " << bi_g << endl;

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
      vector< item_t<W, P> > &items = ukpi.items;
      if (!already_sorted) sort_by_eff(items);

      const I n = static_cast<I>(items.size());
      const W c = ukpi.c;
      // bi stands for 'best item'
      const item_t<W, P> bi = items[0];

      vector<P> p;
      p.reserve(n);
      for (I i = 0; i < n; ++i) p.push_back(items[i].p);
      const P d = gcd_n<typename vector<P>::iterator, P>(p.begin(), p.end());

      vector<P> f(c + 1, 0);
      myvector<I> i;
      i.resize(c + 1);

      // z is another name for sol.opt
      P &z = sol.opt;
      // We init z with the profit value of the solution composed by the 
      // maximum number of best items that the knapsack can hold.
      W bi_qt = (c / bi.w);
      z = bi.p * bi_qt;

      // STEP 1
      const vector<P> b = compute_b(items);

      W bi_qt_lb;

      if (!compute_bi_qt_lb(b, items, c, z, d, bi_qt, bi_qt_lb)) {
        return;
      }

      for (I j = 1; j < n; ++j) {
        const W first_y_reserved_for_best_items = (c - bi_qt_lb*bi.w) + 1;
        const item_t<W, P> it = items[j];
        if (it.w < first_y_reserved_for_best_items && it.p > f[it.w]) {
          f[it.w] = it.p;
          i[it.w] = j;
        }
      }

      // TODO create variable at extra_info to show the last y value
      W y = 1;
      W y_opt_no_bi = 0;
      W bi_qt_in_z = 0;
      P max_previous_fy = 0;
      for (; bi_qt >= bi_qt_lb; --bi_qt) {
        for (; y < c - bi_qt*bi.w; ++y) {
          if (f[y] <= max_previous_fy) continue;

          max_previous_fy = f[y];

          // STEP 2a, 2b, 2c, 2d
          // This is very similar to the loop over items of UKP5.
          // The difference is the j=1 and the first_y_reserved_for_best_items
          for (I j = 1; j <= i[y]; ++j) {
            const item_t<W, P> it = items[j];
            const W new_y = y + it.w;
            const W first_y_reserved_for_best_items = (c - bi_qt_lb*bi.w) + 1;
            if (new_y < first_y_reserved_for_best_items
                && it.p + f[y] > f[new_y]) {
              f[new_y] = it.p + f[y];
              i[new_y] = j;
            }
          }
        }

        // STEP 3d
        // STEP 1 (Routine next z)
        // Note that now y == c - bi_qt*bi.w,
        // and max_previous_fy == max(f[0..y-1])
        const P z_ = std::max(max_previous_fy, f[y]) + bi.p*bi_qt;
        if (z_ > z) {
          z = z_;
          y_opt_no_bi = y;
          bi_qt_in_z = bi_qt;
          
          if (!compute_bi_qt_lb(b, items, c, z, d, bi_qt, bi_qt_lb)) {
            break;
          }
        }
        // The W type can be unsigned, so we have to stop before it underflows
        if (bi_qt == 0) break;
      }

      vector<I> qts_its(n, 0);

      W y_opt = y_opt_no_bi;
      while (y_opt > 0 && f[y_opt] < z - bi_qt_in_z*bi.p) --y_opt;

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
    }

    template<typename W, typename P, typename I>
    void greendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // I tried to implement the algorithm in a way that anyone can verify
      // this is the same algorithm described at the paper. The goto construct
      // is used because of this. The only ommited goto's are the ones that
      // jump to the next step (that are plainly and clearly unecessary here).
      // BEFORE STEP 1 (on the bold 'Algorithm.' block)
      auto &items = ukpi.items;
      if (!already_sorted) {
        sort_by_eff(items);
        reverse(items.begin(), items.end());
      }

      // CHANGING DATA STRUCTURES TO HAVE A NOTATION SIMILAR TO THE ARTICLE
      const I m = static_cast<I>(items.size());
      const W b = ukpi.c;

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
      vector<P> f(b + 1, 0);
      vector<I> i(b + 1);
      i[0] = 1;

      P &z = sol.opt; // z is another name for sol.opt
      z = cm * (b / am);
      W k = 0;

      // STEP 1
      for (I i = 1; i < m; ++i) {
        b_[i] = cm*a[i] - c[i]*am;
      }
      W y = 0;
      W k_max;
      if (!compute_k_max(m, b_, a, c, am, cm, b, z, d, k, lambda, k_max)) {
        goto stop;
      }

      // STEP 2a
      I j, v;
      step_2a:
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
      
      // STEP 3a
      step_3a:
      y = y + 1;

      // STEP 3b
      if (f[y] > f[y - 1]) {
        goto step_3c;
      } else {
        f[y] = f[y - 1];
        i[y] = m + 1;
        goto step_3d;
      }

      // STEP 3c
      step_3c:
      if (y == upper_l) {
        if (!next_z(m, b_, a, c, am, cm, b, z, d, k, lambda, k_max, f, y, upper_l)) {
          goto stop;
        }
      }
      goto step_2a;
      // Because the goto above, the code below can only be accessed by jumping
      // directly into step_3d

      // STEP 3d
      step_3d:
      if (y == upper_l) {
        if (!next_z(m, b_, a, c, am, cm, b, z, d, k, lambda, k_max, f, y, upper_l)) {
          goto stop;
        }
      }
      goto step_3a;
      // Because the goto above, the code below can only be accessed by jumping
      // directly into stop
      stop:
      sol.y_opt = y;

      vector<I> qts_its(m, 0);

      W y_opt = lambda;
      W bi_qt = b/am;
      for (; y_opt < b && f[y_opt] != z - bi_qt*cm; y_opt += am, --bi_qt);
      while (y_opt > 0 && f[y_opt-1] == f[y_opt]) --y_opt;

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

