#include <boost/rational.hpp>
#include <iostream>
#include <cmath>
#include "bb.hpp"

using namespace std;
using namespace boost;
using namespace hbm;

/* TODO: Instead of starting at zero elements of some item i,
 * start with the greatest possible quantity of elements of that
 * type, and then decreases the number (if it's greater than zero)
 * This helps mainly because we examine first the solutions with the
 * least of capacity left (most probability of being optimal).
 * TODO: The bound computed is very stupid. The items are already
 * ordered by eff, then if the first item can't be used anymore,
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

/* WARNING: THIS CODE DON'T WORK IF ANY OF THE HBM_*_EFF MACROS IS DEFINED */
/* When any of these macros is defined the efficiency tpe is already defined.
 * The code has yet to be adapted to use the efficiency type defined by
 * ukp_common.hpp */
typedef long double efficiency;

#ifdef HBM_NOT_DEFINED
template <class T>
struct solution {
  const vector<item_t> * const items;
  vector<T> s;
  weight w
  profit p;
  const weight c;

  solution(const vector<item_t> * const items, const weight c) : items(items), s(items->size(), 0), w(0), p(0), c(c) {};
  inline T get_qt_item_i(const size_t i) { return s[i]; }
  inline void increment(const size_t i, const T new_qt = 1) {
    T old_qt = s[i], diff;
    diff = new_qt - old_qt;
    weight new_w = w + diff*items->at(i).w;
    if (new_w > c) return;
    w = new_w;
    p += diff*items->at(i).p;
  }
  inline void decrement(const size_t i, const T new_qt = 1) {
    T old_qt = s[i], diff;
    diff = old_qt - new_qt;
    w -= diff*items->at(i).w;
    p -= diff*items->at(i).p;
  }
  inline void set_qt_item_i(const size_t i, const T new_qt) {
    T old_qt = s[i];
    if (old_qt > new_qt) {
      decrement(i, new_qt);
    } else {
      increment(i, new_qt);
    }
  }
  inline void unsafe_increment(const size_t i, const T new_qt = 1) {
    T old_qt = s[i], diff;
    diff = (new_qt - old_qt);
    w += diff*items->at(i).w;
    p += diff*items->at(i).p;
  }
  inline void unsafe_set_qt_item_i(const size_t i, const T new_qt) {
    T old_qt = s[i];
    if (old_qt > new_qt) {
      decrement(i, new_qt);
    } else {
      unsafe_increment(i, new_qt);
    }
  }
  inline bool is_valid(void) { return w <= c; }
  inline bool is_invalid(void) { return w > c; }

  /* returns the number of items added */
  inline size_t fill_with_item_i(const size_t i) {
    weight diff = (c-w)/items->at(i).w;
    unsafe_set_qt_item_i(i, diff);
    return diff;
  }
  /* returns the item index of the last item used */
  inline size_t fill_in_order(const size_t from = 0) {
    weight i = from, size = items->size(), last_item_type_used = 0;
    for (; i < size; ++i) {
      if (fill_with_item_i(i) > 0) last_item_type_used = i;
    }
    return last_item_type_used;
  }
};
typedef solution<size_t> solution_t;
#endif /* HBM_NOT_DEFINED */

struct stack_t {
  size_t i;
  weight ws;
  profit ps;

  stack_t(void) {}
  stack_t(size_t i, weight ws, profit ps) : i(i), ws(ws), ps(ps) {}

  inline void set(size_t i_, weight ws_, profit ps_) {
    i = i_;
    ws = ws_;
    ps = ps_;
  }
};

inline profit upper_bound(const vector<efficiency> &effs, const vector<item_t> &items, size_t i, weight c) {
  profit u = 0;
  for (; i < items.size(); ++i) {
    weight wi = items[i].w;
    profit pi = items[i].p;
    weight qt = c/wi;
    u += static_cast<profit>(qt)*pi;
    c -= qt*wi;
    if (qt == 0) {
      //for (; i < items.size(); ++i) if (c >= items[i].w) break;
      u += static_cast<profit>(effs[i+1]*static_cast<efficiency>(c));
      break;
    }
  }
 
  return u;
}

void inner_bb(const vector<item_t> &items, const weight c, const vector<weight> &w_mins, const vector<efficiency> &effs, profit &bp) {
  stack_t s(0, 0, 0);
  vector<stack_t> stack;
  stack.push_back(s);

  size_t nodes_popped = 0;
  while (stack.size() > 0) {
    ++nodes_popped;
    s = stack.back();
    stack.pop_back();
    size_t &i = s.i;
    weight &ws = s.ws;
    profit &ps = s.ps;
    /*cout << "i: " << i << endl;
    cout << "ws: " << ws << endl;
    cout << "ps: " << ps << endl;
    cout << endl;*/

    if (i >= items.size()) continue;

    // If is impossible to get a solution better than the bound
    // from this partial solution, then stop
    //size_t j = i + 1;
    //for (; j < items.size(); ++j) if ((c-ws) >= items[j].w) break;
    //if (ps + floor(((double)(c-ws))*effs[j]) <= bp) continue;
    //if (ps + ceil(((double)(c-ws))*effs[i+1]) <= bp) continue;
    if (ps + upper_bound(effs, items, i + 1, c-ws) <= bp) continue;

    // If the solution is better than the bound, update the bound
    if (ps > bp) bp = ps;
     
    // Branch only if it is possible to add items to the solution
    //if (c - ws < w_mins[i+1]) continue;

    // Test solutions with every possible quantity of i
    weight wi = items[i].w;
    profit pi = items[i].p;
    weight max_qt = (c-ws)/wi;
    for (weight qt = 0; qt <= max_qt ; ++qt) {
      stack.emplace_back(i+1, ws + qt*wi, ps + static_cast<profit>(qt)*pi);
    }
  }
  cout << "Nodes popped: " << nodes_popped << endl;
}

/*
void inner_bb(const vector<item_t> &items, const weight c, const vector<weight> &w_mins, const vector<efficiency> &effs, profit &bp, size_t i, weight ws, profit ps) {
  if (i >= items.size()) return;

  // Fill the capacity with item i
  weight wi = items[i].w;
  profit pi = items[i].p;
  weight qt = (c-ws)/wi;
  ws += qt*wi;
  ps += qt*pi;

  // If is impossible to get a solution better than the bound
  // from this partial solution, then stop
  if (ps + static_cast<profit>(static_cast<efficiency>(c-ws)*effs[i+1]) <= bp) return;

  // If the solution is better than the bound, update the bound
  if (ps > bp) bp = ps;
   
  // Branch only if it is possible to add items to the solution
  if (c - ws < w_mins[i+1]) return;

  // Test solutions with every possible quantity of i
  for (; qt > 0; --qt) {
    inner_bb(items, c, w_mins, effs, bp, i+1, ws, ps);
    // Decrement the weight and the profit of the item removed
    ws -= wi;
    ps -= pi;
  }
  inner_bb(items, c, w_mins, effs, bp, i+1, ws, ps);
}*/

/*void inner_bb(const weight w_min, const vector<efficiency> &effs, profit &bp, profit u, size_t i, solution_t s) {
  // If is impossible to get a solution better than the bound
  // from this partial solution, then stop
  if (s.p + static_cast<profit>(static_cast<efficiency>(s.c-s.w)*effs[i+1]) <= bp) return;

  // If the solution is better than the bound, update the bound
  if (s.p > bp) bp = s.p;

  while (s.get_qt_item_i(i) == 0 && i <= u) --i;
  while (s.get_qt_item_i(i) > 0) {
    s.decrement(i);
    size_t j = i;
    solution_t s_ = s;
    j = s_.fill_in_order(j+1);
    inner_bb(w_min, effs, bp, j, s_);
  }
}*/

vector<efficiency> efficiencies(vector<item_t> &items) {
  vector<efficiency> effs;
  effs.reserve(items.size());
  for (auto it = items.begin()+1; it <items.end(); ++it) {
    //effs.emplace_back(it->p, it->w);
    effs.emplace_back(static_cast<efficiency>(it->p)/static_cast<efficiency>(it->w));
  }
  return effs;
}

vector<weight> weight_mins(vector<item_t> &items) {
  vector<weight> w_mins(items.size(), 0);
  w_mins.back() = items.back().w;
  for (quantity i = items.size() - 1; i > 0; --i) {
    w_mins[i] = min(items[i].w, w_mins[i+1]);
  }
  w_mins[0] = min(items[0].w, w_mins[1]);
  return w_mins;
}

void hbm::bb(ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted/* = false*/) {
  if (!already_sorted) sort_by_eff(ukpi.items);

  vector<weight> w_mins = weight_mins(ukpi.items);
  w_mins.push_back(0);
  vector<efficiency> effs = efficiencies(ukpi.items);
  effs.emplace_back(1);

  sol.opt = 0;
  inner_bb(ukpi.items, ukpi.c, w_mins, effs, sol.opt);
}

