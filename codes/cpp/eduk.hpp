#ifndef HBM_EDUK_HPP
#define HBM_EDUK_HPP

#include <algorithm>
#include <memory>
#include <forward_list>

#include "ukp_common.hpp"

#ifndef HBM_MAX_MEMORY_WASTED_BY_S 
  #define HBM_MAX_MEMORY_WASTED_BY_S 10
#endif

namespace hbm {
  namespace hbm_eduk_impl {
    using namespace std;
    using namespace hbm;

    template<typename W, typename P>
    bool weight_order(const item_t<W, P>& i, const item_t<W, P>& j) {
      return i.w < j.w;
    }

    template<typename W, typename P>
    void sort_by_weigth(vector< item_t<W, P> > &items) {
      sort(items.begin(), items.end(), weight_order<W, P>);

      return;
    }

    template<typename W, typename P>
    struct LazyList {
        virtual bool get_next(item_t<W, P>& x) = 0;
    };

    template<typename W, typename P>
    void consume(LazyList<W, P> * const l, vector< item_t<W, P> > &res) {
      item_t<W, P> i;
      bool has_head = l->get_next(i);
      while (has_head) {
        res.push_back(i);
        has_head = l->get_next(i);
      }
    }

    /* Makes reference to the v argument, v can't be freed before 
     * this struct */
    template<typename W, typename P>
    struct Lazyfy : LazyList<W, P> {
      typename vector< item_t<W, P> >::const_iterator begin;
      typename vector< item_t<W, P> >::const_iterator end;
      
      Lazyfy(void) {
      }

      Lazyfy(const vector< item_t<W, P> > &v) {
        begin = v.cbegin();
        end = v.cend();
      }

      bool get_next(item_t<W, P>& i) {
        if (begin != end) {
          i = *begin;
          ++begin;
          return true;
        } else {
          return false;
        }
      }
    };

    template<typename W, typename P>
    struct AddTest : LazyList<W, P> {
      item_t<W, P> g;
      LazyList<W, P> * l;
      W c;

      AddTest(void) {
        #ifdef HBM_TWO_MULT_COMP
        g = {0, 0};
        #elif defined(HBM_INT_EFF) || defined(HBM_FP_EFF) || defined(HBM_RATIONAL_EFF)
        g = {0, 0, 0} /* define efficiency manually as zero */
        #endif

        l = nullptr;
        c = 0;
      }

      AddTest(item_t<W, P> g, LazyList<W, P> * l, W c) : g(g), l(l), c(c) { }

      bool get_next(item_t<W, P>& i) {
        bool has_next = l->get_next(i);
        i.w += g.w;
        i.p += g.p;
        return has_next && i.w <= c;
      }
    };

    template<typename W, typename P>
    struct Filter : LazyList<W, P> {
      LazyList<W, P> * l;
      P limit;

      Filter(void) {
        limit = 0;
        l = nullptr;
      }

      Filter(LazyList<W, P> * l, P limit) : l(l), limit(limit) { }

      bool get_next(item_t<W, P>& i) {
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

    template<typename W, typename P>
    struct Merge : LazyList<W, P> {
      LazyList<W, P> * l1, * l2;
      item_t<W, P> h1, h2;
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

      Merge(LazyList<W, P> * l1, LazyList<W, P> * l2) : l1(l1), l2(l2) {
        has_old_h1 = has_old_h2 = false;
      }

      bool get_next(item_t<W, P>& i) {
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

    template<typename W, typename P>
    struct AddHead : LazyList<W, P> {
      item_t<W, P> original_head;
      LazyList<W, P> * l;

      bool has_original_head;

      AddHead(void) {
        l = nullptr;
        has_original_head = false;
        #ifdef HBM_TWO_MULT_COMP
        original_head = {0, 0};
        #elif defined(HBM_INT_EFF) || defined(HBM_FP_EFF) || defined(HBM_RATIONAL_EFF)
        original_head = {0, 0, 0} /* define efficiency manually as zero */
        #endif
      }

      AddHead(const item_t<W, P> &head, LazyList<W, P> * l) : original_head(head), l(l) {
        has_original_head = true;
      }

      bool get_next(item_t<W, P>& i) {
        if (has_original_head) {
          i = original_head;
          has_original_head = false;
          return true;
        } else {
          return l->get_next(i);
        }
      }
    };

    template<typename W, typename P, typename I = size_t>
    struct S;
    template<typename W, typename P, typename I>
    struct S {
      unique_ptr< S<W, P, I> > s_pred_k;
      Lazyfy<W, P>  l_items;
      AddHead<W, P> addhead;
      AddTest<W, P> addtest;
      Merge<W, P>   merge;
      Filter<W, P>  filter;

      bool empty;

      vector< item_t<W, P> > computed;

      struct iterator : LazyList<W, P> {
        S<W, P, I> * s;
        size_t ix;

        iterator(S<W, P, I> * s) : s(s) {
          ix = 0;
        }

        bool get_next(item_t<W, P> &i) {
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
      forward_list< S<W, P, I>::iterator > its;

      S(I k, W c, const vector< item_t<W, P> > &items) {
        if (k > 0) { 
          empty = false;

          s_pred_k = unique_ptr< S<W, P, I> >(new S<W, P, I>(k-1, c, items));
          l_items = Lazyfy<W, P>(items);
          #ifdef HBM_TWO_MULT_COMP
          item_t<W, P> zero = {0, 0};
          #elif defined(HBM_INT_EFF) || defined(HBM_FP_EFF) || defined(HBM_RATIONAL_EFF)
          item_t<W, P> zero = {0, 0, 0}; /* define efficiency manually as zero */
          #endif
          addhead = AddHead<W, P>(zero, this->begin());
          addtest = AddTest<W, P>(items[k-1], &addhead, c);
          merge = Merge<W, P>(s_pred_k->begin(), &addtest);
          filter = Filter<W, P>(&merge, 0);
        } else {
          empty = true;
        }
      }

      bool has_next(void) {
        if (empty) return false;

        item_t<W, P> next;
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
          vector< item_t<W, P> > new_v(computed.cbegin()+m, computed.cend());
          computed.swap(new_v);

          for (auto it = its.begin(); it != its.end(); ++it) {
            it->ix = it->ix - m;
          }
        }
      }
    };

    template<typename W, typename P, typename I>
    void eduk(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted/* = false*/) {
      I n = ukpi.items.size();
      W c = ukpi.c;
      vector< item_t<W, P> > &items(ukpi.items);
      if (!already_sorted) sort_by_weigth(ukpi.items);

      S<W, P, I> s = S<W, P, I>(n, c, items);
      auto it = s.begin();

      vector< item_t<W, P> > res;
      consume(it, res);

      sol.opt = res.back().p;
    }
  }

  template<typename W, typename P, typename I = size_t>
  void eduk(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted = false) {
    hbm_eduk_impl::eduk(ukpi, sol, already_sorted);
  }
}

#endif //HBM_EDUK_HPP
