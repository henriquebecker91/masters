#include <boost/rational.hpp>
#include "bb.hpp"

using namespace std;
using namespace boost;

void inner_bb(vector<item_t>::iterator begin, vector<item_t>::iterator end, const size_t w, const size_t p, const size_t c, const rational<size_t> bi_eff, size_t &bp) {
  // There are no items
  if (begin == end) return;

  // If the solution is invalid or is impossible to get a solution better
  // than the bound from this partial solution, then stop
  if (w > c || p + (c-w)*bi_eff < bp) return;

  // If the solution is valid and better than the bound, update the bound
  if (p > bp) bp = p;

  // Solutions that does not use one more of the item pointed by begin
  inner_bb(begin + 1, end, w, p, c, bi_eff, bp);
  // Solutions that use one more of the item pointed by begin
  inner_bb(begin, end, w + begin->w, p + begin->p, c, bi_eff, bp);
}

void bb(ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted/* = false*/) {
  if (!already_sorted) sort_by_efficiency(ukpi.items);
  const size_t wb = ukpi.items[0].w, pb = ukpi.items[0].p;
  const rational<size_t> bi_eff(pb, wb);

  sol.opt = (ukpi.c/wb)*pb;
  inner_bb(ukpi.items.begin(), ukpi.items.end(), 0, 0, ukpi.c, bi_eff, sol.opt);
}

