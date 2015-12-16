#ifndef HBM_DOMINANCE_HPP
#define HBM_DOMINANCE_HPP

#include "ukp_common.hpp"
#include "wrapper.hpp"

namespace hbm {
  namespace hbm_dominance_impl {
    using namespace std;

    template < typename P, typename W, typename I = size_t >
    void sel_not_sm_dom(vector < item_t < W, P > >&items,
                        vector < item_t < W, P > >&undominated,
                        bool already_sorted = false) {
      // the vector HAS TO BE sorted by non-increasing profitability first
      // AND non-decreasing weight after for this code work
      if (!already_sorted) sort_by_eff(items);

      undominated.reserve(items.size() + undominated.size());
      undominated.push_back(items.front());

      // If an item i dominate an item j, then pi/wi >= pj/wj.
      // (Note that the reverse isn't always true, if pi/wi >= pj/wj then
      // i can dominate j, or not. If the reverse was true then every UKP
      // problem could be reduced to the best item.)
      // The first statement only means that: if i is the index of an item
      // and the items are ordered by efficiency, then it can't be simple
      // or multiple dominated by any item of index j where j > i.
      // This ordering guarantees that the most efficient item is never
      // dominated and that we never need to remove items from our
      // undominated items list.
      for (auto i : items) {
        bool i_is_dominated = false;
        for (auto u : undominated) {
          // The line below only works if we assume that W is an integer
          // number and that when an integer number divides another the
          // result is truncated.
          if (static_cast<P>(i.w/u.w) * u.p >= i.p) {
            i_is_dominated = true;
            break;
          }
        }
        if (!i_is_dominated) undominated.push_back(i);
      }
    }

    template<typename W, typename P, typename I>
    void sel_not_sm_dom_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted = false) {
      vector< item_t<W, P> > undominated;
      hbm_dominance_impl::sel_not_sm_dom(ukpi.items, undominated, already_sorted);
      sol.opt = undominated.size();
      sol.y_opt = 0;
      sol.last_y_value_outer_loop = 0;
    }

    template<typename W, typename P, typename I>
    struct sel_not_sm_dom_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_dominance_impl::sel_not_sm_dom_wrapper(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "sel_not_sm_dom";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void sel_not_sm_dom_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, char** argv) {
      simple_wrapper(sel_not_sm_dom_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }
  }
  
  /// @brief Copy the items that aren't simple or multiple dominated
  ///   from a vector to another.
  ///
  /// Maybe you wat to use shrink_to_fit over undominated after this
  /// procedure. We reserve memory as there wasn't simple or multiple
  /// dominated items in 'items'.
  ///
  /// @param items Original set of items. Isn't const because can
  ///   be sorted by efficiency.
  /// @param undominated The non-dominated items will be added at
  ///   the end of this vector. It doesn't need to be empty, but if
  ///   it isn't, the items already there will be used for the
  ///   dominance computation. To have the guarantee that the result
  ///   is correct all the items in undominated must be at least as 
  ///   efficient as the most efficient item in the 'items' argument,
  ///   and items already in undominated can't dominate other items in
  ///   undominated.
  /// @param already_sorted If is false, the 'items' argument is
  ///   modified (will be sorted by efficiency).
  template < typename P, typename W, typename I = size_t >
  void sel_not_sm_dom(std::vector < item_t < W, P > >&items,
                      std::vector < item_t < W, P > >&undominated,
                      bool already_sorted = false) {
    hbm_dominance_impl::sel_not_sm_dom(items, undominated, already_sorted);
  }

  template<typename W, typename P, typename I>
  void sel_not_sm_dom_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted = false) {
    hbm_dominance_impl::sel_not_sm_dom_wrapper(ukpi, sol, already_sorted);
  }

  template<typename W, typename P, typename I = size_t>
  void sel_not_sm_dom_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, char** argv) {
    hbm_dominance_impl::sel_not_sm_dom_wrapper(ukpi, sol, argc, argv);
  }
}

#endif //HBM_DOMINANCE_HPP
