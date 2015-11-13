#ifndef HBM_PERIODICITY_HPP
#define HBM_PERIODICITY_HPP
#include "ukp_common.hpp"

size_t y_star(ukp_instance_t &ukpi, bool already_sorted/* = false*/);
void y_star_wrapper(ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted/* = false*/);

#endif //HBM_PERIODICITY_HPP
