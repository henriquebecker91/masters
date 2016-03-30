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

    /// Compute the gcd for a variable quantity of numbers.
    /// @note This functions was designed with positive integers on mind.
    template<typename T>
    T gcd_n(vector<T>::iterator begin, vector<T>::iterator end) {
    //size_t gcd_n(vector<size_t>::iterator begin, vector<size_t>::iterator end) {
      /*if (begin == end) return 1;
      if (begin + 1 == end) return *begin;

      T acc = gcd(*begin, *(begin+1));
      ++begin;

      for (; begin != end && acc > 1; ++begin) {
        acc = gcd(acc, *begin);
      }
      
      return acc;*/
      return 1;
    }

    /// Compute k_max, also, if returns false the greendp procedure (caller)
    /// is to be stopped.
    template<typename W, typename P, typename I>
    void compute_k_max(const vector<I> &i, W &k_max) {

    }

    template<typename W, typename P, typename I>
    void greendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      if (!already_sorted) {
        sort_by_eff(ukpi.items);
        reverse(ukpi.items.begin(), ukpi.items.end());
      }

      // Change to notation used at the article
      const I m = static_cast<I>(ukpi.items.size());
      const W b = ukpi.c;
      vector<W> a(m + 1);
      a[0] = 0; // Does not exist, notation begins at 1
      vector<P> c(m + 1);
      c[0] = 0; // Does not exist, notation begins at 1
      for (I j = 0; j < m; ++j) {
        item_t<W, P> it = ukpi.items[j];
        a[j+1] = it.w;
        c[j+1] = it.p;
      }

      W y = 0;
      W lambda = b - (b / a[m]) * a[m];
      W upper_l = lambda;
      vector<P> f(c + 1, 0);
      vector<I> i(c + 1);
      i[0] = 1;
      P z = c[m] * (b / a[m]);
      W k = 0;

      W k_max;
      if (!compute_k_max(i, k_max)) goto stop;

      
      stop:
      return;
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

