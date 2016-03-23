#ifndef HBM_GARDP_HPP
#define HBM_GARDP_HPP

#include "ukp_common.hpp"
#include "wrapper.hpp"

namespace hbm {
  namespace hbm_gardp_impl {
    using namespace std;

    template<typename W, typename P, typename I>
    void gardp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      const W c = ukpi.c;
      const I n = ukpi.items.size();

      if (!already_sorted) sort_by_eff(ukpi.items);

      std::vector<P> g(c + 1, 0);
      std::vector<I> d(c + 1, n-1);

      for (W y = 0; y <= c; ++y) {
        for (I j = 0; j < n; ++j) {
          item_t<W, P> it = ukpi.items[j];
          if (y >= it.w && j <= d[y - it.w]) {
            P x = g[y - it.w] + it.p;
            if (x > g[y]) {
              g[y] = x;
              d[y] = j;
            }
          }
        }
      }

      // Found the optimal solution value
      sol.opt = g[c];

      // Find  the weight of the smallest optimal solution
      sol.y_opt = c;
      while (sol.y_opt > 0 && g[sol.y_opt-1] == g[sol.y_opt]) --sol.y_opt;

      // Assemble the smallest optimal solution
      vector<I> qts_its(n, 0);

      W y_opt = sol.y_opt;
      I dy_opt;
      while (y_opt != 0) {
        dy_opt = d[y_opt];
        y_opt -= ukpi.items[dy_opt].w;
        ++qts_its[dy_opt];
      }

      for (I i = 0; i < n; ++i) {
        if (qts_its[i] > 0) {
          sol.used_items.emplace_back(ukpi.items[i], qts_its[i], i);
        }
      }
      sol.used_items.shrink_to_fit();
    }

    template<typename W, typename P, typename I>
    struct gardp_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_gardp_impl::gardp(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "gardp";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void gardp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol,
                int argc, argv_t argv) {
      simple_wrapper(gardp_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }
  }

  /// Solves an UKP instance by the dynamic programming algorithm presented at
  /// p. 221, Integer Programming, Robert S. Garfinkel, and stores the results
  /// at sol. 
  ///
  /// @param ukpi The UKP instance to be solved.
  /// @param sol The object where the results will be written.
  /// @param already_sorted If the ukpi.items vector needs to be sorted by
  ///   non-decreasing efficiency.
  template<typename W, typename P, typename I>
  void gardp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_gardp_impl::gardp(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items isn't sorted by non-decreasing weight. If it's
  /// ommited the ukpi.items is sorted by non-decreasing weight.
  ///
  /// @see main_take_path
  /// @see gardp(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void gardp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_gardp_impl::gardp(ukpi, sol, argc, argv);
  }
}

#endif //HBM_GARDP_HPP
