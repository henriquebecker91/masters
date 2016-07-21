#ifndef HBM_PERIODICITY_HPP
#define HBM_PERIODICITY_HPP

#include <type_traits>  // For is_integral, is_same
#include <limits> // for numeric_limits
#include <boost/multiprecision/cpp_int.hpp> // For cpp_rational (infinite
                                            // precision boost rational)
                                            // used at y_star
#include <boost/math/common_factor_rt.hpp>  // For boost::math::lcm used
                                            // at huangtang
#include "ukp_common.hpp"
#include "wrapper.hpp"

namespace hbm {
  template <typename W>
  struct per_extra_info_t : extra_info_t {
    std::string info;

    per_extra_info_t(W original_cap, W y_cap) {
      info = std::string("algorithm_name: y_star_periodicity_bound\n")
             + "Original capacity: " + std::to_string(original_cap) + "\n"
             + "y* capacity: " + std::to_string(y_cap) + "\n";
    }

    virtual std::string gen_info(void) {
      return info;
    }
  };

  namespace hbm_periodicity_impl {
    using namespace std;
    using namespace boost;
    using namespace boost::multiprecision;

    // From "A constructive periodicity bound for the unbounded
    //  knapsack problem"
    // Its complexity is bigger than O(n^2), and it seems worse than
    // y*. The overall impression is that it is terrible.
    template <typename W, typename P>
    W huangtang(instance_t<W, P> &ukpi, bool already_sorted = false) {
      auto &items = ukpi.items;
      if (!already_sorted) sort_by_eff(items);

      size_t n = items.size();
      W h0 = 0;
      for (size_t j = 1; j < n; ++j) {
        W min = boost::math::lcm(items[0].w, items[j].w) - items[j].w;
        for (size_t i = 1; i < j; ++i) {
          W x = boost::math::lcm(items[i].w, items[j].w) - items[j].w;
          if (x < min) min = x;
        }
        h0 += min;
      }

      return h0 + 1;
    }

    template <typename W, typename P>
    W y_star(const item_t<W, P> &b, const item_t<W, P> &b2) {
      assert(b < b2);
      if (std::is_integral<W>::value && std::is_same<W, P>::value) {
        // Without the castings, the compiler gives conversion warnings
        // that don't will ever happen. If P is a floating point number
        // this 'if' never executes. The castings have literally no
        // effect since inside this if W and P are the same type.
        W w1 = b.w, w2 = b2.w;
        W p1 = static_cast<W>(b.p), p2 = static_cast<W>(b2.p);

        cpp_rational r_p1 = static_cast<W>(b.p);
        cpp_rational r1(p1, w1);
        cpp_rational r2(p2, w2);

        // If the two numbers are equal we would divide by
        // zero later. The closest value we have to infinity
        // is the best return.
        if (r1 == r2) return numeric_limits<W>::max();

        // Always positive: the r1 efficiency is bigger than the r2 efficiency
        cpp_rational d = r1 - r2;
        // Final value, the "plus one" is to avoid getting 1 less
        // than the real value when we cast back to W
        cpp_rational y = (r_p1 / d) + 1;

        return y <= numeric_limits<W>::max() ? static_cast<W>(y) : numeric_limits<W>::max();
      } else if (std::is_floating_point<P>::value) {
        W w1 = b.w, w2 = b2.w;
        P p1 = b.p, p2 = b2.p;

        P e1 = p1 / static_cast<P>(w1);
        P e2 = p2 / static_cast<P>(w2);

        if (e1 - e2 < numeric_limits<P>::epsilon())
          return numeric_limits<W>::max();

        return static_cast<W>(p1/(e1 - e2)) + 1;
      } else {
        cerr << __func__ << ": W and P aren't valid types. " << endl;
        exit(EXIT_FAILURE);
      }
    }

    template <typename W>
    W refine_y_star(W y_, W c, W w_b) {
      if (y_ > c) return c;
      W qt_b = ((c - y_) / w_b) + 1;
      return c - qt_b*w_b;
    }

    template <typename W, typename P>
    W y_star(vector< item_t<W, P> > &items, bool already_sorted = false) {
      if (!already_sorted) sort_by_eff(items, 2u);

      return hbm_periodicity_impl::y_star(items[0], items[1]);
    }

    template <typename W, typename P>
    W y_star(instance_t<W, P> &ukpi, bool already_sorted = false) {
      return y_star(ukpi.items, already_sorted);
    }

//    template <typename W, typename P, typename I>
//    W run_with_y_star(void(*ukp_solver)(instance_t<W, P> &, solution_t<W, P, I> &, void*),
//      instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, void* ukp_solver_extra_params) {
//      W y_ = y_star(ukpi, false);
//
//      if (y_ >= ukpi.c) {
//        (*ukp_solver)(ukpi, sol, ukp_solver_extra_params);
//        return y_;
//      }
//
//      vector< item_t<W, P> > &items(ukpi.items);
//
//      W old_c = ukpi.c;
//
//      W w1, p1;
//      w1 = items[0].w;
//      p1 = items[0].p;
//
//      W qt_best_item_used = (old_c - y_)/w1;
//      P profit_generated_by_best_item = static_cast<P>(qt_best_item_used)*p1;
//      W space_used_by_best_item = qt_best_item_used*w1;
//
//      ukpi.c = old_c - space_used_by_best_item;
//
//      (*ukp_solver)(ukpi, sol, true);
//
//      sol.opt += profit_generated_by_best_item;
//
//      return y_;
//    }

    template <typename W, typename P, typename I>
    void y_star_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted = false) {
      sol.show_only_extra_info = true;
      I y_star_cap = hbm_periodicity_impl::y_star(ukpi, already_sorted);

      per_extra_info_t<W>* ptr =
        new per_extra_info_t<W>(ukpi.c, y_star_cap);

      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(ptr);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      return;
    }

    template<typename W, typename P, typename I>
    struct y_star_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_periodicity_impl::y_star_wrapper(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "y_star";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void y_star_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
      simple_wrapper(y_star_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }
  }

  /// Computes the y* periodicity bound. It's guaranteed that any optimal
  /// solution for a capacity bigger than this bound have at least one copy of
  /// the best item. Taken from Garfinkel and Nemhauser at "Integer
  /// Programming", p. 223.
  ///
  /// @param b The best item, i.e. the most efficient one.
  /// @param b2 The second best item, i.e. the second
  ///   most efficient one.
  ///
  /// @return The first capacity value that have guarantee that any optimal
  ///   solution will contain a copy of the best item.
  template <typename W, typename P>
  W y_star(const item_t<W, P> &b, const item_t<W, P> &b2) {
    return hbm_periodicity_impl::y_star(b, b2);
  }

  /// Based on the y* bound, gives you the capacity for what you
  /// should compute the UKP solution, to after fill with remaining
  /// space with copies of the best item.
  ///
  /// Don't use the y* bound as capacity. Use the value
  /// returned by this function. This value is guaranteed to be
  /// equal to or smaller than y*, and the difference between it
  /// and c is always a multiple of w_b (if y_ > c it will be c,
  /// and the difference will be zero, that is a multiple of any
  /// number). This way, to know how many copies of the best item
  /// you should add to the solution you only need to compute this
  /// value divided by w_b.
  ///
  /// @param y_ The value obtained by y_star.
  /// @param c The original capacity value of the instance.
  /// @param w_b The weight of the best item, i.e the most
  ///   efficient one.
  ///
  /// @return A safe capacity to compute the result, and then fill the
  ///   remaining space with exactly (c - <this return value>)/best_item.w
  ///   copies of the best item.
  template <typename W>
  W refine_y_star(W y_, W c, W w_b) {
    return hbm_periodicity_impl::refine_y_star(y_, c, w_b);
  }

  /// A terrible bound only implemented to access that it is indeed
  /// terrible. Kept only for comparison.
  template <typename W, typename P>
  W huangtang(instance_t<W, P> &ukpi, bool already_sorted = false) {
    return hbm_periodicity_impl::huangtang(ukpi, already_sorted);
  }

  /// Executes y_star over the two best items of ukpi.items. Assumes that ukpi
  /// has at least two items, and that if already_sorted is true, the items are
  /// ordered by non-increasing efficiency.
  ///
  /// @see y_star(const item_t<W, P> &b, const item_t<W, P> &b2)
  template <typename W, typename P>
  W y_star(instance_t<W, P> &ukpi, bool already_sorted = false) {
    return hbm_periodicity_impl::y_star(ukpi, already_sorted);
  }

//  template <typename W, typename P, typename I = size_t>
//  W run_with_y_star(void(*ukp_solver)(instance_t<W, P> &, solution_t<W, P, I> &, bool),
//    instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted = false) {
//    return hbm_periodicity_impl::run_with_y_star(ukp_solver, ukpi, sol, already_sorted);
//  }

  /// Convenience overload, executes y_star over ukpi.items and saves the
  /// "solution" (the result of y_star, not an optimal solution) to sol. A hack
  /// used to allow y_star to make use of main_take_path and benchmark_pyasukp
  /// procedures.
  ///
  /// @see main_take_path
  /// @see benchmark_pyasukp
  /// @see per_extra_info_t
  template <typename W, typename P, typename I = size_t>
  void y_star_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted = false) {
    hbm_periodicity_impl::y_star_wrapper(ukpi, sol, already_sorted);
  }

  /// Other convenience overload, executes y_star over ukpi.items and saves the
  /// "solution" (the result of y_star, not an optimal solution) to sol. A hack
  /// used to allow y_star to make use of main_take_path and benchmark_pyasukp
  /// procedures.
  ///
  /// @see main_take_path
  /// @see benchmark_pyasukp
  /// @see per_extra_info_t
  template<typename W, typename P, typename I = size_t>
  void y_star_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_periodicity_impl::y_star_wrapper(ukpi, sol, argc, argv);
  }
}

#endif //HBM_PERIODICITY_HPP
