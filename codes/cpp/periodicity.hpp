#ifndef HBM_PERIODICITY_HPP
#define HBM_PERIODICITY_HPP

#include <type_traits>
#include <boost/rational.hpp>
#include "ukp_common.hpp"
#include "wrapper.hpp"

namespace hbm {
  template <typename W>
  struct per_extra_info_t : extra_info_t {
    std::string info;

    per_extra_info_t(W original_cap, W y_cap) {
      info = "Original capacity: " + std::to_string(original_cap) + "\n"
             + "y* capacity: " + std::to_string(y_cap) + "\n";
    }

    virtual std::string gen_info(void) {
      return info;
    }
  };

  namespace hbm_periodicity_impl {
    using namespace std;
    using namespace boost;

    template <typename W, typename P>
    W y_star(item_t<W, P> b, item_t<W, P> b2) {
      if (std::is_integral<W>::value && std::is_integral<P>::value &&
          std::is_same<W, P>::value) {
        // Without the pragma, the compiler gives conversion warnings
        // that don't will ever happen. If P is a floating point number
        // this 'if' never executes.
        #pragma GCC diagnostic ignored "-Wfloat-conversion"
        W w1 = b.w, w2 = b2.w, p1 = b.p, p2 = b2.p;
        
        rational<W> r_p1(b.p);
        rational<W> r1(p1, w1);
        rational<W> r2(p2, w2);

        return rational_cast<W>(r_p1/(r1 - r2)) + 1;
        // Return diagnostic to normal
        #pragma GCC diagnostic warning "-Wfloat-conversion"
      } else if (std::is_floating_point<P>::value) {
        W w1 = b.w, w2 = b2.w;
        P p1 = b.p, p2 = b2.p;

        P e1 = p1 / static_cast<P>(w1);
        P e2 = p2 / static_cast<P>(w2);

        return static_cast<W>(p1/(e1 - e2)) + 1;
      } else {
        cerr << __func__ << ": W and P aren't valid types. " << endl;
        exit(EXIT_FAILURE);
      }
    }

    template <typename W, typename P>
    W y_star(vector< item_t<W, P> > &items, bool already_sorted = false) {
      if (!already_sorted) sort_by_eff(items, 2u);

      return y_star(items[0], items[1]);
    }

    template <typename W, typename P>
    W y_star(instance_t<W, P> &ukpi, bool already_sorted = false) {
      return y_star(ukpi.items, already_sorted);
    }

    template <typename W, typename P, typename I>
    W run_with_y_star(void(*ukp_solver)(instance_t<W, P> &, solution_t<W, P, I> &, void*),
      instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, void* ukp_solver_extra_params) {
      W y_ = y_star(ukpi, false);

      if (y_ >= ukpi.c) {
        (*ukp_solver)(ukpi, sol, ukp_solver_extra_params);
        return y_;
      }

      vector< item_t<W, P> > &items(ukpi.items);

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
      sol.show_only_extra_info = true;
      I y_star_cap = hbm_periodicity_impl::y_star(ukpi, already_sorted);

      per_extra_info_t<W>* ptr =
        new per_extra_info_t<W>(ukpi.c, y_star_cap);

      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(ptr);
      sol.extra_info = shared_ptr<extra_info_t>(upcast_ptr);

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

  /* Assumes that ukpi has at least two items */
  template <typename W, typename P>
  W y_star(instance_t<W, P> &ukpi, bool already_sorted = false) {
    hbm_periodicity_impl::y_star(ukpi, already_sorted);
  }
  template <typename W, typename P, typename I = size_t>
  W run_with_y_star(void(*ukp_solver)(instance_t<W, P> &, solution_t<W, P, I> &, bool),
    instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted = false) {
    hbm_periodicity_impl::run_with_y_star(ukp_solver, ukpi, sol, already_sorted);
  }
  template <typename W, typename P, typename I = size_t>
  void y_star_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted = false) {
    hbm_periodicity_impl::y_star_wrapper(ukpi, sol, already_sorted);
  }

  template<typename W, typename P, typename I = size_t>
  void y_star_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_periodicity_impl::y_star_wrapper(ukpi, sol, argc, argv);
  }
}

#endif //HBM_PERIODICITY_HPP
