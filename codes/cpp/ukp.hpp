#ifndef __UKP_HPP_
#define __UKP_HPP_

#include <vector>
#include <cstdlib>
#include <istream>

struct item_t {
  size_t p;
  size_t w;
};

struct ukp_instance_t {
  size_t c;
  std::vector<item_t> items;
};

struct ukp_solution_t {
  std::vector<size_t> g;
  std::vector<size_t> d;
  size_t opt;
};

void ukp5(ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted = false);

void read_sukp_instance(std::istream &in, ukp_instance_t &ukpi);

void write_sukp_instance(std::ostream &out, ukp_instance_t &ukpi);

#endif //__UKP_HPP_ 

