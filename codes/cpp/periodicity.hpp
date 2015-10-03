#ifndef __PERIODICITY_HPP_
#define __PERIODICITY_HPP_

#include "ukp_common.hpp"

size_t y_star(ukp_instance_t &ukpi, bool already_sorted/* = false*/);
void y_star_wrapper(ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted/* = false*/);

#endif //__PERIODICITY_HPP_ 

