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
    void gcd_n(ForwardIterator begin, ForwardIterator end, T &acc) {
      if (begin == end) {
        acc = 1;
        return;
      }
      if (begin + 1 == end) {
        acc = *begin;
        return;
      }

      acc = gcd(*begin, *(begin+1));
      ++begin;

      for (; begin != end && acc > 1; ++begin) {
        acc = gcd(acc, *begin);
      }
      
      return;
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

    /// Compute k_max, also, if returns false the greendp procedure (caller)
    /// is to be stopped.
    template<typename W, typename P>
    bool compute_k_max( const vector<P> &b,
                        const vector< item_t<W, P> > &items,
                        const W c,
                        const P opt,
                        const P d,
                        const W k,
                        const W lambda,
                        W &k_max) {
      // STEP 1 (Routine k_max)
      const item_t<W, P> &bi = items.back();

      size_t r = items.size() - 1;
      while (r > 0 && b[r] > bi.p*c - bi.w*(opt + d)) --r;


      if (r == 0) return false;

      // STEP 2 (Routine k_max)
      const P delta = opt - bi.p*(c/bi.w) + d; 
      const P l_side = (items[r].p*lambda - items[r].w*delta)/b[r];
      const P max_bi_qt = c/bi.w;
      k_max = l_side < max_bi_qt ? l_side : max_bi_qt;

      if (k > k_max) return false;

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
    template<typename W, typename P>
    vector<P> compute_b(const vector< item_t<W, P> > &items) {
      const size_t n = items.size();
      const item_t<W, P> &bi = items.back();

      vector<P> b;
      b.reserve(items.size() - 1);
      b.push_back(0);

      for (size_t i = 1; i < n; ++i) {
        const item_t<W, P> &it = items[i-1];
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
      vector< item_t<W, P> > &items = ukpi.items;
      // Greenberg used the items on non-decreasing efficiency order
      if (!already_sorted) {
        sort_by_eff(items);
        reverse(items.begin(), items.end());
      }

      // CHANGING DATA STRUCTURES TO HAVE A NOTATION SIMILAR TO THE ARTICLE
      const I n = static_cast<I>(items.size());
      const W c = ukpi.c;

      /*vector<W> a(n + 1);
      a[0] = 0; // Does not exist, notation begins at 1*/
      vector<P> p;
      p.reserve(n + 1);
      p.push_back(0); // Does not exist, notation begins at 1

      //const W am = ukpi.items[n - 1].w;
      //const P cm = ukpi.items[n - 1].p;
      const item_t<W, P> bi = items[n - 1];
      for (I i = 0; i < n; ++i) {
        p.push_back(items[i].p);
        //item_t<W, P> it = ukpi.items[i];
        //a[i+1] = it.w;
        //p[i+1] = it.p;
      }

      // Starts at index 1, this is the cause of the +1
      W tmp;
      gcd_n(p.begin() + 1, p.end(), tmp);
      const W d = tmp;

      // lambda is the remaining space if we fill the capacity b
      // with the most efficient item.
      const W lambda = c - (c / bi.w) * bi.w;
      W upper_l = lambda;
      vector<P> f(c + 1, 0);
      myvector<I> i;
      i.resize(c + 1);


      // z is another name for sol.opt
      P &z = sol.opt;
      // We init z with the profit value of the solution composed by the 
      // maximum number of best items that the knapsack can hold.
      z = bi.p * (c / bi.w);
      W k = 0;

      // STEP 1
      // As b is the capacity, b will be used where b was subscribed
      const vector<P> b = compute_b(items);

      W k_max;
      if (!compute_k_max(b, items, c, z, d, k, lambda, k_max)) {
        return;
      }

      for (I j = 1; j <= n - 1; ++j) {
        const W first_y_reserved_for_best_items = lambda + bi.w*k_max + 1;
        const item_t<W, P> it = items[j];
        if (it.w < first_y_reserved_for_best_items && it.p >= f[it.w]) {
          f[it.w] = it.p;
          i[it.w] = j;
        }
      }

      bool enumerate_solutions;
      for (W y = 1; k <= k_max; ++y) {
        enumerate_solutions = true;
        // STEP 3b
        if (f[y] <= f[y - 1]) {
          f[y] = f[y - 1];
          i[y] = n + 1;
          enumerate_solutions = false;
        }

        // STEP 3d
        if (y == upper_l) {
          // STEP 1 (Routine next z)
          const P z_ = f[y] + bi.p*(c/bi.w) - bi.p*k;
          if (z_ > z) {
            z = z_;
            if (!compute_k_max(b, items, c, z, d, k, lambda, k_max)) {
              return;
            }
          }

          // STEP 2 (Routine next z)
          k = k + 1;
          if (k > k_max) return;

          upper_l = upper_l + bi.w;
        }

        // STEP 2a, 2b, 2c, 2d
        // This is very similar to the loop over items of UKP5.
        // The difference is the "m-1" and the "lambda + am*k_max".
        if (enumerate_solutions) {
          for (I j = i[y]; j <= n - 1; ++j) {
            const item_t<W, P> it = items[j];
            const W new_y = y + it.w;
            const W first_y_reserved_for_best_items = lambda + bi.w*k_max + 1;
            const P new_p = it.p + f[y];
            const P old_p = f[new_y];
            if (new_y < first_y_reserved_for_best_items && new_p >= old_p) {
              f[new_y] = new_p;
              i[new_y] = j;
            }
          }
        }
      }
    }

    template<typename W, typename P, typename I>
    void greendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // I tried to implement the algorithm in a way that anyone can verify
      // this is the same algorithm described at the paper. The goto construct
      // is used because of this. The only ommited goto's are the ones that
      // jump to the next step (that are plainly and clearly unecessary here).
      // BEFORE STEP 1 (on the bold 'Algorithm.' block)
      if (!already_sorted) {
        sort_by_eff(ukpi.items);
        reverse(ukpi.items.begin(), ukpi.items.end());
      }

      // CHANGING DATA STRUCTURES TO HAVE A NOTATION SIMILAR TO THE ARTICLE
      const I m = static_cast<I>(ukpi.items.size());
      const W b = ukpi.c;

      vector<W> a(m + 1);
      a[0] = 0; // Does not exist, notation begins at 1
      vector<P> c(m + 1);
      c[0] = 0; // Does not exist, notation begins at 1
      // As b is the capacity, b_ will be used where b was subscribed
      vector<P> b_(m + 1);
      b_[0] = 0; // Does not exist, notation begins at 1

      const W am = ukpi.items[m - 1].w;
      const P cm = ukpi.items[m - 1].p;
      for (I i = 0; i < m; ++i) {
        item_t<W, P> it = ukpi.items[i];
        a[i+1] = it.w;
        c[i+1] = it.p;
      }

      // Starts at index 1, this is the cause of the +1
      W tmp;
      gcd_n(c.begin() + 1, c.end(), tmp);
      const W d = tmp;

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
//      for (size_t x = 0; x <= y+1; ++x) {
//        cout << "f[" << x << "]: " << f[x] << "\ti[" << x << "]: " << i[x] << endl;
//      }
      // TODO: ASSEMBLE SOLUTION
      /*I n = ukpi.items.size();
      vector<I> qts_its(n, 0);

      W y_opt = y;
//      W y_opt = 53312;
      while (y_opt != 0 && i[y_opt] == m + 1) --y_opt;
      W old_y_opt = y_opt;
      I dy_opt;
      while (y_opt != 0) {
        dy_opt = i[y_opt] - 1;
        y_opt -= ukpi.items[dy_opt].w;
        ++qts_its[dy_opt];
      }

      for (I i = 0; i < n; ++i) {
        if (qts_its[i] > 0) {
          sol.used_items.emplace_back(ukpi.items[i], qts_its[i], i + 1);
        }
      }
      sol.used_items.shrink_to_fit();

      W qtd_best_item = (b - old_y_opt)/am;
      if (qtd_best_item > 0) {
        sol.used_items.emplace_back(ukpi.items[n-1], qtd_best_item, m);
      }*/
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

