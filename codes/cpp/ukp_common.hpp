#ifndef __UKP_COMMON_HPP_
#define __UKP_COMMON_HPP_

#include <vector>
#include <cstdlib>
#include <istream>
#include <boost/rational.hpp>

#define fp_type float

struct item_t {
  size_t w;
  size_t p;
  //boost::rational<size_t> efficiency;
  //fp_type efficiency;
  //size_t efficiency;

  inline item_t(void) {}
  inline item_t(const size_t &w, const size_t &p) : w(w), p(p) {}
  //inline item_t(const size_t &w, const size_t &p) : w(w), p(p), efficiency(p, w) {}
  /*inline item_t(size_t w, size_t p) : w(w), p(p) {
    efficiency = ((fp_type)p) / ((fp_type)w);
  }*/
  /*item_t(size_t w, size_t p) : w(w), p(p) {
    efficiency = (p << 32) / w;
  }*/

  inline bool operator==(const item_t &o) const {
    return p == o.p && w == o.w;
  }

  /*inline bool operator<(const item_t &o) const {
    return efficiency == o.efficiency ? w < o.w : efficiency > o.efficiency;
  }*/

  /*inline bool operator<(const item_t &o) const {
    return efficiency > o.efficiency;
  }*/

  inline bool operator<(const item_t &o) const {
    size_t a = o.p * w, b = p * o.w;
    return a == b ? w < o.w : a < b;
  }

  /*inline bool operator<(const item_t &o) const {
    size_t a = o.p * w, b = p * o.w;
    return a < b;
  }*/

  /*inline bool operator>>(int s) const {
    return efficiency >> s;
  }*/
};

struct ukp_instance_t {
  size_t c;
  std::vector<item_t> items;
  /*item_t best_item;
  item_t second_best_item;*/
};

struct ukp_solution_t {
  std::vector<size_t> g;
  std::vector<size_t> d;
  std::vector<item_t> res;
  size_t opt;
};

void read_sukp_instance(std::istream &in, ukp_instance_t &ukpi);

void write_sukp_instance(std::ostream &out, ukp_instance_t &ukpi);

bool efficiency_order(const item_t& i, const item_t& j);

void sort_by_efficiency(std::vector<item_t> &items);
#endif //__UKP_COMMON_HPP_ 

