#ifndef HBM_CHUNG_HPP
#define HBM_CHUNG_HPP

#include "ukp_common.hpp"

namespace hbm {
  namespace hbm_chung_impl {
    using namespace std;

    template<typename W, typename P, typename I>
    void gen_chung_hard_inst(I n, W a, P z, W b, instance_t<W, P> &ukpi) {
      ukpi.c = b;
      ukpi.items.resize(n);
      for (I j = 0; j < n; ++j) {
        ukpi.items[j].w = a + j;
        ukpi.items[j].p = z + static_cast<P>(a + j);
      }
    }
  }

  // -------------------- EXTERNAL FUNCTIONS --------------------

  /// Create instances by the method described on "A Hard Knapsack Problem",
  /// by Chung et al. The description of the method wasn't from "A Hard
  /// Knapsack Problem" but from "A New Knapsack Solution Approach by Integer
  /// Equivalent Aggregation and Consistency Determination", as the former
  /// is behind a paywall. All items follow the formulae:
  /// a_j = a + j - 1
  /// c_j = a_j + z
  /// j = 1, 2, ..., n
  /// 
  /// @param n The number of items.
  /// @param a Value added to an item index to get its weight.
  /// @param z Value added to an item weight to get its profit value.
  /// @param b The knapsack capacity.
  /// @param ukpi Instance object. Will be overwritten with generated instance.
  /// @note Instance returned will be ordered by ascending weight, maybe by
  ///   efficiency too, depending on the arguments.
  template<typename W, typename P, typename I>
  void gen_chung_hard_inst(I n, W a, P z, W b, instance_t<W, P> &ukpi) {
    hbm_chung_impl::gen_chung_hard_inst<W, P, I>(n, a, z, b, ukpi);
  }
}

#endif //HBM_CHUNG_HPP

