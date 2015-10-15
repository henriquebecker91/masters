#ifndef __TEST_COMMON_HPP_
#define __TEST_COMMON_HPP_

#include "ukp_common.hpp"
#include <chrono>

struct run_t {
  ukp_solution_t result;
  std::chrono::duration<double> time;
};

struct instance_data_t {
  std::string name;
  size_t expected_opt;
};

int run_ukp(void(*ukp_solver)(ukp_instance_t &, ukp_solution_t &, bool), const std::string& path, run_t &run);
int benchmark_pyasukp(void(*ukp_solver)(ukp_instance_t &, ukp_solution_t &, bool));
int main_take_path(void(*ukp_solver)(ukp_instance_t &, ukp_solution_t &, bool), int argc, char** argv);

#endif //__TEST_COMMON_HPP_ 

