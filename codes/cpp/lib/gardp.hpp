#ifndef HBM_GARDP_HPP
#define HBM_GARDP_HPP

#include "ukp_common.hpp"
#include "periodicity.hpp"
#include "wrapper.hpp"
#include "type_name.hpp"

namespace hbm {
  template <typename W, typename P, typename I>
  struct gardp_extra_info_t : extra_info_t {
    gardp_extra_info_t(void) { }

    virtual std::string gen_info(void) {
      return std::string("algorithm_name: gardp\n")
           + "git_head_at_compilation: " + HBM_GIT_HEAD_AT_COMPILATION + "\n"
           + "type_W: " + hbm::type_name<W>::get() + "\n"
           + "type_P: " + hbm::type_name<P>::get() + "\n"
           + "type_I: " + hbm::type_name<I>::get() + "\n";
    }
  };

  namespace hbm_gardp_impl {
    using namespace std;

    template<typename W, typename P, typename I>
    void gardp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      gardp_extra_info_t<W, P, I>* eip = new gardp_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      const W original_c = ukpi.c;
      W c = ukpi.c;
      const I n = ukpi.items.size();

      if (!already_sorted) sort_by_eff(ukpi.items);

      I b1_ix = 0;
      item_t<W, P> b1 = ukpi.items[0];
      item_t<W, P> b2 = ukpi.items[1];

      W y_ = y_star(b1, b2);
      W y_bound = refine_y_star(y_, c, b1.w);
      if (y_bound < c) c = y_bound;

      std::vector<P> g(c + 1, 0);
      std::vector<I> d(c + 1, n-1);

      for (W y = 0; y <= c; ++y) {
        for (I j = 0; j < n; ++j) {
          item_t<W, P> it = ukpi.items[j];
          if (y >= it.w && j <= d[y - it.w]) {
            const P x = g[y - it.w] + it.p;
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

      if (y_bound < original_c) {
        W qt_b = (original_c - sol.y_opt)/b1.w;
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
  ///   non-increasing efficiency.
  template<typename W, typename P, typename I>
  void gardp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_gardp_impl::gardp(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items isn't sorted by non-increasing efficiency. If it's
  /// ommited the ukpi.items is sorted by non-increasing efficiency.
  ///
  /// @see main_take_path
  /// @see gardp(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void gardp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_gardp_impl::gardp(ukpi, sol, argc, argv);
  }
}

#endif //HBM_GARDP_HPP
