#include "eduk.hpp"
#include <algorithm>

using namespace std;

bool weight_order(const item_t& i, const item_t& j) {
  return i.w < j.w;
}

void sort_by_weigth(vector<item_t> &items) {
  sort(items.begin(), items.end(), weight_order);

  return;
}

struct LazyList {
    virtual bool get_next(item_t& x) {
      return false;
    };
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
  LazyList * const l;
  size_t c;

  bool has_head;
  item_t head;

  AddTest(item_t g, LazyList * const l, size_t c) : g(g), l(l), c(c) {
    has_head = l->get_next(head);
    head.w += g.w;
    head.p += g.p;
    if (!has_head || head.w > c) {
      has_head = false;
    } 
  }

  bool get_next(item_t& i) {
    if (has_head) i = head;
    else return false;

    has_head = l->get_next(head);
    head.w += g.w;
    head.p += g.p;
    if (!has_head || head.w > c) {
      has_head = false;
    }
    return true;
  }
};

struct Filter : LazyList {
  size_t limit;
  LazyList * const l;

  bool has_head;
  item_t head;

  Filter(LazyList * const l) : l(l) {
    has_head = l->get_next(head);
    if (has_head) limit = head.p;
  }

  bool get_next(item_t& i) {
    if (has_head) i = head;
    else return false;

    has_head = l->get_next(head);
    while (has_head && head.p <= limit) {
      has_head = l->get_next(head);
    } 
    if (has_head) {
      limit = head.p;
    }

    return true;
  }
};

struct Merge : LazyList {
  LazyList * const l1, * const l2;
  item_t h1, h2;
  bool has_h1, has_h2;

  Merge(LazyList * const l1, LazyList * const l2) : l1(l1), l2(l2) {
    has_h1 = l1->get_next(h1);
    has_h2 = l2->get_next(h2);
  }

  bool get_next(item_t& i) {
    if (has_h1 && has_h2) {
      if (h1.w > h2.w) {
        i = h1;
        has_h1 = l1->get_next(h1);
      } else if (h1.w < h2.w) {
        i = h2;
        has_h2 = l2->get_next(h2);
      } else {
        i.w = h1.w;
        i.p = max(h1.p, h2.p);
        has_h1 = l1->get_next(h1);
        has_h2 = l2->get_next(h2);
      }
    } else if (has_h1) {
      i = h1;
      has_h1 = l1->get_next(h1);
    } else if (has_h2) {
      i = h2;
      has_h2 = l2->get_next(h2);
    } else { /* The only possibility is: has_h1 == false && has_h2 == false */
      return false;
    }
    return true;
  }
};

struct AddHead : LazyList {
  const item_t original_head;
  LazyList * const l;

  bool has_original_head;

  AddHead(const item_t &head, LazyList * const l) : original_head(head), l(l) {
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

struct S : LazyList {
  LazyList s_pred_k, l_items, addhead, addtest, merge, filter;

  bool has_head;
  item_t head;

  S(size_t k, size_t c, const vector<item_t> &items) {
    if (k > 0) { 
      s_pred_k = S(k-1, c, items);
      l_items = Lazyfy(items);
      item_t zero;
      zero.w = 0;
      zero.p = 0;
      addhead = AddHead(zero, this);
      addtest = AddTest(items[k-1], &addhead, c);
      merge = Merge(&s_pred_k, &addtest);
      filter = Filter(&merge);

      has_head = filter.get_next(head);
    } else {
      has_head = false;
    }
  }

  bool get_next(item_t &i) {
    if (has_head) i = head;
    else return false;

    has_head = filter.get_next(head);

    return true;
  }
};

void eduk(ukp_instance_t &ukpi, ukp_solution_t &sol, bool already_sorted) {
    vector<item_t> &res = sol.res;
    size_t n = ukpi.items.size();
    size_t c = ukpi.c;
    vector<item_t> &items(ukpi.items);
    if (!already_sorted) sort_by_weigth(ukpi.items);

    S s = S(n, c, items);
    consume(&s, res);

    //sol.opt = res[res.size() - 1].p;
}

