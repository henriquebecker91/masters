#ifdef INT_EFF
#include <boost/sort/spreadsort/spreadsort.hpp>
#else
#include <algorithm>
#endif

#include "ukp_common.hpp"

using namespace std;

void read_sukp_instance(istream &in, ukp_instance_t &ukpi) {
  size_t n;
  in >> n;
  in >> ukpi.c;
  ukpi.items.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    size_t w, p;
    in >> w;
    in >> p;
    ukpi.items.emplace_back(w, p);
  }

  return;
}

void write_sukp_instance(ostream &out, ukp_instance_t &ukpi) {
  size_t n = ukpi.items.size();
  out << n << endl;
  out << ukpi.c << endl;

  for (size_t i = 0; i < n; ++i) {
    item_t tmp = ukpi.items[i];
    out << tmp.w << "\t" << tmp.p << endl;
  }

  return;
}

void sort_by_efficiency(vector<item_t> &items) {
  #ifdef INT_EFF
  boost::sort::spreadsort::integer_sort(items.begin(), items.end());
  #else
  std::sort(items.begin(), items.end());
  #endif
  return;
}

