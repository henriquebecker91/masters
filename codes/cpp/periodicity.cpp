#include "periodicity.hpp"

#include <boost/rational.hpp>

using namespace std;
using namespace boost;
using namespace hbm;

size_t hbm::y_star(ukp_instance_t &ukpi, bool already_sorted = false) {
  vector<item_t> &items(ukpi.items);
  if (!already_sorted) sort_by_eff(ukpi.items);

  item_t i1 = items[0], i2 = items[1];
  //item_t i1 = ukpi.best_item, i2 = ukpi.second_best_item;
  size_t w1 = i1.w, w2 = i2.w, p1 = i1.p, p2 = i2.p;
  
  rational<size_t> r_p1(i1.p);
  rational<size_t> r1(p1, w1);
  rational<size_t> r2(p2, w2);

  return rational_cast<size_t>(r_p1/(r1 - r2)) + 1;
}

size_t hbm::run_with_y_star(void(*ukp_solver)(ukp_instance_t &, ukp_solution_t &, bool), ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted = false) {
  vector<item_t> &items(ukpi.items);
  if (!already_sorted) sort_by_eff(ukpi.items);

  size_t y_ = y_star(ukpi, already_sorted);

  if (y_ >= ukpi.c) {
    (*ukp_solver)(ukpi, sol, true);
    return y_;
  }

  size_t old_c = ukpi.c;

  size_t w1, p1;
  w1 = items[0].w;
  p1 = items[0].p;

  size_t qt_best_item_used = (old_c - y_)/w1;
  size_t profit_generated_by_best_item = qt_best_item_used*p1;
  size_t space_used_by_best_item = qt_best_item_used*w1;

  ukpi.c = old_c - space_used_by_best_item;

  (*ukp_solver)(ukpi, sol, true);

  sol.opt += profit_generated_by_best_item;

  return y_;
}

void hbm::y_star_wrapper(ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted = false) {
  //(void) run_with_y_star(&ukp5, ukpi, sol, false);
  sol.opt = y_star(ukpi, already_sorted);

  return;
}

