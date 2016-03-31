#ifndef HBM_GREENDP_HPP
#define HBM_GREENDP_HPP

#include "ukp_common.hpp"
#include "periodicity.hpp"
#include "wrapper.hpp"

#include <vector> // For vector
#include <algorithm>  // For reverse
#include <boost/math/common_factor_rt.hpp>  // For boost::math::gcd

namespace hbm {
  namespace hbm_greendp_impl {
    using namespace std;
    using namespace boost::math;

    /// Compute the gcd for a variable quantity of numbers.
    /// @note This functions was designed with positive integers on mind.
    template<typename T>
    T gcd_n(typename vector<T>::iterator begin,
            typename vector<T>::iterator end) {
      if (begin == end) return 1;
      if (begin + 1 == end) return *begin;

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
      while (r > 0 && b_[r] <= cm*b - am*(z + d)) --r;

      if (r == 0) return false;

      // STEP 2 (Routine k_max)
      const P delta = z - cm*(b/am) + d;
      const P l_side = (c[r]*lambda - a[r]*delta)/b_[r];
      const P r_side = b/am;
      k_max = l_side > r_side ? l_side : r_side;

      if (k > k_max) return false;

      return true;
    }

    /// Compute next_z, also, if returns false the greendp procedure (caller)
    /// is to be stopped.
    template<typename W, typename P, typename I>
    bool compute_k_max( const I m,
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
                        W upper_l) {
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

    template<typename W, typename P, typename I>
    void greendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // BEFORE STEP 1 (on the bold 'Algorithm.' block)
      if (!already_sorted) {
        sort_by_eff(ukpi.items);
        reverse(ukpi.items.begin(), ukpi.items.end());
      }

      // CHANGING NOTATION TO BE SIMILAR TO THE ARTICLE
      const I m = static_cast<I>(ukpi.items.size());
      const W b = ukpi.c;

      vector<W> a(m + 1);
      a[0] = 0; // Does not exist, notation begins at 1
      vector<P> c(m + 1);
      c[0] = 0; // Does not exist, notation begins at 1
      vector<P> b_(m + 1);
      b_[0] = 0; // Does not exist, notation begins at 1

      const W am = ukpi.items[m - 1].p;
      const P cm = ukpi.items[m - 1].p;
      for (I i = 0; i < m; ++i) {
        item_t<W, P> it = ukpi.items[i];
        a[i+1] = it.w;
        c[i+1] = it.p;
        // Setting b_[i] is already STEP 1, but it's here to
        // reuse the loop
        b_[i+1] = cm*a[i] - c[i]*am;
      }

      // Starts at index 1, this is the cause of the +1
      const W d = gcd_n(c.begin() + 1, c.end());

      // STEP 1
      W y = 0;
      W lambda = b - (b / a[m]) * a[m];
      W upper_l = lambda;
      vector<P> f(c + 1, 0);
      vector<I> i(c + 1);
      i[0] = 1;

      P z = c[m] * (b / a[m]);
      W k = 0;

      W k_max;
      if (!compute_k_max(m, b_, a, c, am, cm, b, z, d, k, lambda, k_max)) {
        return;
      }

      // STEP 2a
      I j, v;
      step_2a:
      j = i[y];

      // STEP 2b
      if (y + a[j] < lambda + am*k_max)  {
        v = c[j] + f[y];
        goto step_2c;
      } else {
        goto step_2d;
      }

      // STEP 2c
      step_2c:
      if (v > f[y + a[j]]) {
        f[y + a[j]] = v;
        i[y + a[j]] = j;
        goto step_2d;
      } else {
        
      }

      // STEP 2d
      //step_2d:
      
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

