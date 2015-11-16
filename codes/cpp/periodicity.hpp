#ifndef HBM_PERIODICITY_HPP
#define HBM_PERIODICITY_HPP
#include "ukp_common.hpp"

namespace hbm {
  size_t y_star(ukp_instance_t &ukpi, bool already_sorted/* = false*/);
  void y_star_wrapper(ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted/* = false*/);
  size_t run_with_y_star(void(*ukp_solver)(ukp_instance_t &, ukp_solution_t &, bool), ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted/* = false*/);
}

#endif //HBM_PERIODICITY_HPP
