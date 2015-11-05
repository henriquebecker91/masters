#include <boost/rational.hpp>
#include <iostream>
#include "bb.hpp"

using namespace std;
using namespace boost;

/* TODO: Instead of starting at zero elements of some item i,
 * start with the greatest possible quantity of elements of that
 * type, and then decreases the number (if it's greater than zero)
 * This helps mainly because we examine first the solutions with the
 * least of capacity left (most probability of being optimal).
 * TODO: The bound computed is very stupid. The items are already
 * ordered by efficiency, then if the first item can't be used anymore,
 * why use it in the bound? The upper bound over the remaining capacity
 * must be done by using the current (or next) item.
 * TODO: The type of pruning done by MTU/MT2 is a periodicity prunning,
 * if you use n-1 of the least efficient item used and cover the 
 * remaining capacity at its fullest (relaxed covering) with the most
 * efficient of the remaining items (that are all less efficient than the
 * current), you can known if the there are possible solution that take
 * less than n items of the best item, this is the same that verifying
 * if for some capacity is more interesting add one more item of i
 * or use the remaining items.
 */

struct bound_info_t {
  size_t bw, bp, bw2, bp2;
  rational<size_t> be3;
};

/* The items must be sorted by nondecreasing efficiency and
 * if tied, by nonincreasing weight, there must be at least 3
 * different items */
void gen_bound_info(const vector<item_t> &items, bound_info_t &bi) {
  bi.bw = items[0].w;
  bi.bp = items[0].p;
  auto it = items.begin() + 1, end = items.end();
  bool f2 = false;
  for (; it < end; ++it) {
    if (it->w < bi.bw) continue;
    if (!f2) {
      f2 = true;
      bi.bw2 = it->w;
      bi.bp2 = it->p;
    } else {
      bi.be3 = rational<size_t>(it->p, it->w);
      break;
    }
  }
}

rational<size_t> lower_bound(size_t c, const bound_info_t &bi) {
  rational<size_t> lb((c/bi.bw)*bi.bp);
  c = c % bi.bw;
  lb += (c/bi.bw2)*bi.bp2;
  c = c % bi.bw2;
  lb += c*bi.be3;

  return lb;
}

item_t heuristic_upper_bound(size_t c, const bound_info_t &bi) {
  item_t i(0, 0);
  i.p += (c/bi.bw)*bi.bp;
  i.w += (c/bi.bw)*bi.bw;
  c = c % bi.bw;
  i.p += (c/bi.bw2)*bi.bp2;
  i.w += (c/bi.bw2)*bi.bw2;
  c = c % bi.bw2;
  i.p += (c/bi.be3.denominator())*bi.be3.numerator();
  i.w += (c/bi.be3.denominator())*bi.be3.denominator();

  return i;
}

void inner_bb(vector<item_t>::iterator begin, vector<item_t>::iterator end, const size_t w_min, const size_t w, const size_t p, const size_t c, const bound_info_t &bi, size_t &bp) {
  // There are no items
  /*cout << "w: " << w << endl;
  cout << "p: " << p << endl;
  cout << "bp: " << bp << endl;
  if (w < c) cout << "pp: " << rational_cast<size_t, size_t>((c-w)*bi_eff) + 1 << endl;
  cout << endl;*/
  if (begin == end) return;

  // If the solution is invalid or is impossible to get a solution better
  // than the bound from this partial solution, then stop
  if (w > c || p + lower_bound(c-w, bi) < bp) return;

  // If the solution is valid and better than the bound, update the bound
  if (p > bp) bp = p;

  // Avoid unnecessary branching
  if (w + w_min > c) return;

  // Solutions that use one more of the item pointed by begin
  inner_bb(begin, end, w_min, w + begin->w, p + begin->p, c, bi, bp);
  // Solutions that does not use one more of the item pointed by begin
  inner_bb(begin + 1, end, w_min, w, p, c, bi, bp);
}

pair<size_t,size_t> minmax_item_weight(vector<item_t> &items) {
  size_t min, max;
  min = max = items[0].w;
  for (auto it = items.begin()+1; it != items.end(); ++it) {
    size_t x = (*it).w;
    if (x < min) min = x;
    else if (x > max) max = x;
  }
  return make_pair(min,max);
}

void bb(ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted/* = false*/) {
  if (!already_sorted) sort_by_efficiency(ukpi.items);
  const auto minw_max = minmax_item_weight(ukpi.items);
  const size_t w_min = minw_max.first, w_max = minw_max.second;

  bound_info_t bi;
  gen_bound_info(ukpi.items, bi);

  item_t ub = heuristic_upper_bound(ukpi.c, bi);

  sol.opt = ub.p;
  inner_bb(ukpi.items.begin(), ukpi.items.end(), w_min, ub.w, ub.p, ukpi.c, bi, sol.opt);
}

