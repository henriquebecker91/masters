#ifndef HBM_PERIODICITY_HPP
#define HBM_PERIODICITY_HPP
#include "ukp_common.hpp"

namespace hbm {
  size_t y_star(instance &ukpi, bool already_sorted/* = false*/);
  void y_star_wrapper(instance &ukpi, solution &sol, bool already_sorted/* = false*/);
  size_t run_with_y_star(void(*ukp_solver)(instance &, solution &, bool), instance &ukpi, solution &sol, bool already_sorted/* = false*/);
}

#endif //HBM_PERIODICITY_HPP
