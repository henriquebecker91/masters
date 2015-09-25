#ifndef __UKP_COMMON_HPP_
#define __UKP_COMMON_HPP_

#include <vector>
#include <cstdlib>
#include <istream>

struct item_t {
  size_t w;
  size_t p;

  bool operator==(const item_t &o) const {
    return p == o.p && w == o.w;
  }
};

struct ukp_instance_t {
  size_t c;
  std::vector<item_t> items;
};

struct ukp_solution_t {
  std::vector<size_t> g;
  std::vector<size_t> d;
  std::vector<item_t> res;
  size_t opt;
};

void read_sukp_instance(std::istream &in, ukp_instance_t &ukpi);

void write_sukp_instance(std::ostream &out, ukp_instance_t &ukpi);

#endif //__UKP_COMMON_HPP_ 

