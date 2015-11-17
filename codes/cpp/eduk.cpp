#include "eduk.hpp"
#include <algorithm>
#include <memory>
#include <forward_list>

#define HBM_MAX_MEMORY_WASTED_BY_S 10

using namespace std;
using namespace hbm;

bool weight_order(const item_t& i, const item_t& j) {
  return i.w < j.w;
}

void sort_by_weigth(vector<item_t> &items) {
  sort(items.begin(), items.end(), weight_order);

  return;
}

struct LazyList {
    virtual bool get_next(item_t& x) = 0;
};

void consume(LazyList * const l, vector<item_t> &res) {
  item_t i;
  bool has_head = l->get_next(i);
  while (has_head) {
    res.push_back(i);
    has_head = l->get_next(i);
  }
}

/* Makes reference to the v argument, v can't be freed before 
 * this struct */
struct Lazyfy : LazyList {
  vector<item_t>::const_iterator begin;
  vector<item_t>::const_iterator end;
  
  Lazyfy(void) {
  }

  Lazyfy(const vector<item_t> &v) {
    begin = v.cbegin();
    end = v.cend();
  }

  bool get_next(item_t& i) {
    if (begin != end) {
      i = *begin;
      ++begin;
      return true;
    } else {
      return false;
    }
  }
};

struct AddTest : LazyList {
  item_t g;
  LazyList * l;
  weight c;

  AddTest(void) {
    #ifdef HBM_TWO_MULT_COMP
    g = {0, 0};
    #elif defined(HBM_INT_EFF) || defined(HBM_FP_EFF) || defined(HBM_RATIONAL_EFF)
    g = {0, 0, 0} /* define efficiency manually as zero */
    #endif

    l = nullptr;
    c = 0;
  }

  AddTest(item_t g, LazyList * l, weight c) : g(g), l(l), c(c) { }

  bool get_next(item_t& i) {
    bool has_next = l->get_next(i);
    i.w += g.w;
    i.p += g.p;
    return has_next && i.w <= c;
  }
};

struct Filter : LazyList {
  LazyList * l;
  profit limit;

  Filter(void) {
    limit = 0;
    l = nullptr;
  }

  Filter(LazyList * l, profit limit) : l(l), limit(limit) { }

  bool get_next(item_t& i) {
    bool has_next = l->get_next(i);
    while (has_next && i.p <= limit) {
      has_next = l->get_next(i);
    } 
    if (has_next) {
      limit = i.p;
    }

    return has_next;
  }
};

struct Merge : LazyList {
  LazyList * l1, * l2;
  item_t h1, h2;
  bool has_h1, has_h2, has_old_h1, has_old_h2;

  Merge(void) {
    has_h1 = has_h2 = has_old_h1 = has_old_h2 = false;
    l1 = l2 = nullptr;
    #ifdef HBM_TWO_MULT_COMP
    h1 = h2 = {0, 0};
    #elif defined(HBM_INT_EFF) || defined(HBM_FP_EFF) || defined(HBM_RATIONAL_EFF)
    h1 = h2 = {0, 0, 0} /* define efficiency manually as zero */
    #endif
  }

  Merge(LazyList * l1, LazyList * l2) : l1(l1), l2(l2) {
    has_old_h1 = has_old_h2 = false;
  }

  bool get_next(item_t& i) {
    if (!has_old_h1) has_h1 = l1->get_next(h1);
    if (!has_old_h2) has_h2 = l2->get_next(h2);

    if ((has_h1 || has_old_h1) && (has_h2 || has_old_h2)) {
      if (h1.w < h2.w) {
        i = h1;
        has_old_h1 = false;
        has_old_h2 = true;
      } else if (h1.w > h2.w) {
        i = h2;
        has_old_h1 = true;
        has_old_h2 = false;
      } else {
        i.w = h1.w;
        i.p = max(h1.p, h2.p);
        has_old_h1 = false;
        has_old_h2 = false;
      }
    } else if (has_h1 || has_old_h1) {
      i = h1;
      has_old_h1 = false;
    } else if (has_h2 || has_old_h2) {
      i = h2;
      has_old_h2 = false;
    } else { /* The only possibility is: all the four flags are false */
      return false;
    }

    return true;
  }
};

struct AddHead : LazyList {
  item_t original_head;
  LazyList * l;

  bool has_original_head;

  AddHead(void) {
    l = nullptr;
    has_original_head = false;
    original_head = {0, 0};
  }

  AddHead(const item_t &head, LazyList * l) : original_head(head), l(l) {
    has_original_head = true;
  }

  bool get_next(item_t& i) {
    if (has_original_head) {
      i = original_head;
      has_original_head = false;
      return true;
    } else {
      return l->get_next(i);
    }
  }
};

struct S;
struct S {
  unique_ptr<S> s_pred_k;
  Lazyfy l_items;
  AddHead addhead;
  AddTest addtest;
  Merge merge;
  Filter filter;

  bool empty;

  vector<item_t> computed;

  struct iterator : LazyList {
    S * s;
    size_t ix;

    iterator(S * s) : s(s) {
      ix = 0;
    }

    bool get_next(item_t &i) {
      // Using the handmade garbage_collect was not paying off
      //s->garbage_collect();
      if (ix < s->computed.size() || s->has_next()) {
        i = s->computed[ix];
        ++ix;
        return true;
      } else {
        return false;
      }
    }
  };

  /* Only to avoid allocating dynamic memory by hand
   * HAS TO BE A LIST, LOST TWO HOURS DEBUGGING BECAUSE
   * USING VECTOR IS A BAD IDEA (WHEN SIZE BECOMES BIGGER THAN CAPACITY
   * THE INTERNAL DATA IS REALLOCATED AND THE OLD POINTERS TO INTERNAL
   * OBJECTS LOSE THEIR MEANING)
   * */
  forward_list<S::iterator> its;

  S(quantity k, weight c, const vector<item_t> &items) {
    if (k > 0) { 
      empty = false;

      s_pred_k = unique_ptr<S>(new S(k-1, c, items));
      l_items = Lazyfy(items);
      item_t zero(0, 0);
      addhead = AddHead(zero, this->begin());
      addtest = AddTest(items[k-1], &addhead, c);
      merge = Merge(s_pred_k->begin(), &addtest);
      filter = Filter(&merge, 0);
    } else {
      empty = true;
    }
  }

  bool has_next(void) {
    if (empty) return false;

    item_t next;
    empty = !filter.get_next(next);
    if (!empty) computed.push_back(next);
    else { // Is this correct?
      computed.clear();
      computed.shrink_to_fit();
    }

    return !empty;
  }

  iterator * begin(void) {
    iterator tmp(this);
    its.push_front(tmp);
    return &its.front();
  }

  void garbage_collect(void) {
    auto it = its.cbegin();
    size_t m = it->ix;
    ++it;

    for (; it != its.cend(); ++it) {
      m = min(m, it->ix);
    }

    /* Another possibility is to verify if (m > computed.size()/MAGIC_CONST) */
    if (m > HBM_MAX_MEMORY_WASTED_BY_S) {
      vector<item_t> new_v(computed.cbegin()+m, computed.cend());
      computed.swap(new_v);

      for (auto it = its.begin(); it != its.end(); ++it) {
        it->ix = it->ix - m;
      }
    }
  }
};

void hbm::eduk(ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted/* = false*/) {
  quantity n = ukpi.items.size();
  weight c = ukpi.c;
  vector<item_t> &items(ukpi.items);
  if (!already_sorted) sort_by_weigth(ukpi.items);

  S s = S(n, c, items);
  auto it = s.begin();

  vector<item_t> res;
  consume(it, res);

  sol.opt = res.back().p;
}

