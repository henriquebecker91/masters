#include "ukp_common.hpp"

using namespace std;

void read_sukp_instance(istream &in, ukp_instance_t &ukpi) {
  size_t n;
  in >> n;
  in >> ukpi.c;
  ukpi.items.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    item_t tmp;
    in >> tmp.w;
    in >> tmp.p;
    ukpi.items.push_back(tmp);
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

