#include <algorithm>
#include <boost/rational.hpp>

#include "ukp_common.hpp"

using namespace std;
using namespace boost;

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

/* Efficiency nonascending first weight nondescending after */
/*bool efficiency_order(const item_t& i, const item_t& j) {
  rational<size_t> ri(i.w, i.p);
  rational<size_t> rj(j.w, j.p);

  return ri == rj ? i.w < j.w : ri < rj;
}*/

//#include <boost/sort/spreadsort/spreadsort.hpp>

void sort_by_efficiency(vector<item_t> &items) {
  sort(items.begin(), items.end());
  //sort::spreadsort::spreadsort(items.begin(), items.end());
  return;
}

