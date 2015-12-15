#ifndef HBM_PERIODICITY_HPP
#define HBM_PERIODICITY_HPP

#include <type_traits>
#include <boost/rational.hpp>
#include "ukp_common.hpp"

namespace hbm {
  namespace hbm_periodicity_impl {
    using namespace std;
    using namespace boost;

    template <typename W, typename P>
    W y_star(instance_t<W, P> &ukpi, bool already_sorted/* = false*/) {
      static_assert(std::is_integral<W>::value &&
                    std::is_integral<P>::value &&
                    std::is_same<W, P>::value,
                    "For now, y_star<W, P> only compiles if the W and P"
                    " types are the same integral type.");
      vector< item_t<W, P> > &items(ukpi.items);
      if (!already_sorted) sort_by_eff(ukpi.items);

      item_t<W, P> i1 = items[0], i2 = items[1];
      W w1 = i1.w, w2 = i2.w, p1 = i1.p, p2 = i2.p;
      
      rational<W> r_p1(i1.p);
      rational<W> r1(p1, w1);
      rational<W> r2(p2, w2);

      return rational_cast<W>(r_p1/(r1 - r2)) + 1;
    }

    template <typename W, typename P, typename I>
    W run_with_y_star(void(*ukp_solver)(instance_t<W, P> &, solution_t<W, P, I> &, bool),
      instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted/* = false*/) {
      vector< item_t<W, P> > &items(ukpi.items);
      if (!already_sorted) sort_by_eff(ukpi.items);

      W y_ = y_star(ukpi, already_sorted);

      if (y_ >= ukpi.c) {
        (*ukp_solver)(ukpi, sol, true);
        return y_;
      }

      W old_c = ukpi.c;

      W w1, p1;
      w1 = items[0].w;
      p1 = items[0].p;

      W qt_best_item_used = (old_c - y_)/w1;
      P profit_generated_by_best_item = static_cast<P>(qt_best_item_used)*p1;
      W space_used_by_best_item = qt_best_item_used*w1;

      ukpi.c = old_c - space_used_by_best_item;

      (*ukp_solver)(ukpi, sol, true);

      sol.opt += profit_generated_by_best_item;

      return y_;
    }

    template <typename W, typename P, typename I>
    void y_star_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted = false) {
      //(void) run_with_y_star(&ukp5, ukpi, sol, false);
      sol.opt = hbm_periodicity_impl::y_star(ukpi, already_sorted);

      return;
    }

    template<typename W, typename P, typename I = size_t>
    void y_star_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, char** argv) {
      // This function don't call itself. It call its overloaded variant
      // where the third parameter is a bool.
      if (argc == 0) {
        y_star_wrapper(ukpi, sol);
      } else if (argc == 1) {
        if ("--already-sorted" == argv[0]) {
          y_star_wrapper(ukpi, sol, true)
        } else {
          cerr << #__func__" (argc/argv overload): parameter error:"
                  " The only allowed flag is --already-sorted."
                  " The flag received was \"" << argv[0] << 
                  "\". Executing the algorithm as no"
                  " flags were given. " << endl;
          y_star_wrapper(ukpi, sol)
        }
      } else {
          cerr << #__func__"(argc/argv overload): parameter error: Only one"
                  " flag is allowed. The allowed flag is "
                  "--already-sorted. The first flag received was \""
                  << argv[0] << "\". Executing the algorithm as no"
                  " flags were given. " << endl;
          y_star_wrapper(ukpi, sol)
      }
    }
  }
  }

  /* Assumes that ukpi has at least two items */
  template <typename W, typename P>
  W y_star(instance_t<W, P> &ukpi, bool already_sorted = false) {
    hbm_periodicity_impl::y_star(ukpi, already_sorted);
  }
  template <typename W, typename P, typename I>
  W run_with_y_star(void(*ukp_solver)(instance_t<W, P> &, solution_t<W, P, I> &, bool),
    instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted = false) {
    hbm_periodicity_impl::run_with_y_star(ukp_solver, ukpi, sol, already_sorted);
  }
  template <typename W, typename P, typename I>
  void y_star_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted = false) {
    hbm_periodicity_impl::y_star_wrapper(ukpi, sol, already_sorted);
  }

  template<typename W, typename P, typename I = size_t>
  void y_star_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, char** argv) {
    hbm_periodicity_impl::y_star_wrapper(ukpi, argc, argc, argv);
  }
}

#endif //HBM_PERIODICITY_HPP
