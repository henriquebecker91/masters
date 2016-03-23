#ifndef HBM_UKP5_HPP
#define HBM_UKP5_HPP

#include "ukp_common.hpp"

namespace hbm {

  namespace hbm_dpgar_impl {
    using namespace std;

  }

  /// Solves an UKP instance by the dynamic programming algorithm , and stores the results (and
  /// stats) at sol.
  ///
  /// @param ukpi The UKP instance to be solved.
  /// @param sol The object where the results will be written.
  /// @param conf An object with the UKP5 options.
  ///
  /// @see ukp5_conf_t For the explanation of the UKP5 options.
  template<typename W, typename P, typename I>
  void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, const ukp5_conf_t<I> &conf) {
    hbm_ukp5_impl::ukp5(ukpi, sol, conf);
  }

  /// For information on argv valid arguments consult the procedure "usage",
  /// or call this function with argv = null_ptr and argc = 0.
  template<typename W, typename P, typename I>
  void ukp5(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_ukp5_impl::ukp5(ukpi, sol, argc, argv);
  }
}

