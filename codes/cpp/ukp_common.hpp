#ifndef HBM_UKP_COMMON_HPP
#define HBM_UKP_COMMON_HPP

/* includes that don't depend on type definitions */
#include <vector>
#include <iostream>
#include <stdexcept> /* for runtime_error */

/* type definitions */
#include "ukp_types.hpp"

/* includes that depend on type definitions */
#if defined(HBM_TWO_MULT_COMP) || defined(HBM_INT_EFF)
#include <utility> /* to specialize swap */
#endif

namespace hbm {
  /* The following code doesn't allow any two of the following macros to be
   * defined at the same time: HBM_TWO_MULT_COMP, HBM_INT_EFF, HBM_FP_EFF
   * and HBM_RATIONAL_EFF. Also it defines the efficiency type based on 
   * the macro value*/
  #ifndef HBM_RATIONAL_EFF
    #ifndef HBM_FP_EFF
      #ifndef HBM_INT_EFF
        #ifndef HBM_TWO_MULT_COMP
          #error ONE OF: HBM_RATIONAL_EFF, HBM_FP_EFF, HBM_INT_EFF OR HBM_TWO_MULT_COMP MUST BE \
                 DEFINED, THE RECOMMENDED IS HBM_TWO_MULT_COMP
        #endif
      #else /*HBM_INT_EFF*/
        #ifdef HBM_TWO_MULT_COMP
          #error HBM_INT_EFF AND HBM_TWO_MULT_COMP CANNOT BE DEFINED AT THE SAME TIME
        #endif
        typedef int_eff efficiency;
      #endif /*HBM_INT_EFF*/
    #else /*HBM_FP_EFF*/
      #ifdef HBM_TWO_MULT_COMP
        #error HBM_FP_EFF AND HBM_TWO_MULT_COMP CANNOT BE DEFINED AT THE SAME TIME
      #endif
      #ifdef HBM_INT_EFF
        #error HBM_FP_EFF AND HBM_INT_EFF CANNOT BE DEFINED AT THE SAME TIME
      #endif
      typedef fp_eff efficiency;
    #endif /*HBM_FP_EFF*/
  #else /*HBM_RATIONAL_EFF*/
    #ifdef HBM_TWO_MULT_COMP
      #error HBM_RATIONAL_EFF AND HBM_TWO_MULT_COMP CANNOT BE DEFINED AT THE SAME TIME
    #endif
    #ifdef HBM_INT_EFF
      #error HBM_RATIONAL_EFF AND HBM_INT_EFF CANNOT BE DEFINED AT THE SAME TIME
    #endif
    #ifdef HBM_FP_EFF
      #error HBM_RATIONAL_EFF AND HBM_FP_EFF CANNOT BE DEFINED AT THE SAME TIME
    #endif
    typedef rational_eff efficiency;
  #endif /*HBM_RATIONAL_EFF*/


  struct item_t {
    weight w;
    profit p;
    #if defined(HBM_RATIONAL_EFF) || defined(HBM_INT_EFF) || defined(HBM_FP_EFF)
    efficiency eff;
    #endif

    inline item_t(void) {}
    #ifdef HBM_TWO_MULT_COMP
    inline item_t(const weight &w, const profit &p) : w(w), p(p) {}
    #elif defined(HBM_RATIONAL_EFF)
    inline item_t(const weight &w, const profit &p) : w(w), p(p), eff(p, w) {}
    #elif defined(HBM_INT_EFF)
    inline item_t(const weight &w, const profit &p) : w(w), p(p) {
      eff = (p << 32) / w;
    }
    #elif defined(HBM_FP_EFF)
    inline item_t(const weigth &w, const profit &p) : w(w), p(p) {
      eff = static_cast<efficiency>(p) / static_cast<efficiency>(w);
    }
    #endif

    inline bool operator==(const item_t &o) const {
      return p == o.p && w == o.w;
    }

    /* Sort by non-increasing eff, if the efficiences are equal
     * sort by non-decreasing weight
     */
    #if defined(HBM_TWO_MULT_COMP)
    inline bool operator<(const item_t &o) const {
      profit a = p * static_cast<profit>(o.w),
             b = o.p * static_cast<profit>(w);
      return a > b || (a == b && w < o.w); 
    }
    #elif defined(HBM_RATIONAL_EFF) || defined(HBM_FP_EFF) || defined(HBM_INT_EFF)
    inline bool operator<(const item_t &o) const {
      return eff > o.eff || (eff == o.eff && w < o.w);
    }
    #endif

    #ifdef HBM_INT_EFF
    inline size_t operator>>(const int s) const {
      /* NOTE: this operator is needed for boost::sort::spreadsort::integer_sort
       * that is the sort we use when HBM_INT_EFF is defined.
       * As we sort by non-increasing eff, we cannot return
       * (eff >> s), because the integer_sort will compare the value
       * trying to sort them in non-decreasing order, as (a < b) iff (~a >= ~b)
       * this solves the problem. */
      return (~eff) >> s;
    }
    #endif
  };

  struct ukp_instance_t {
    weight c;
    std::vector<item_t> items;
  };

  struct ukp_itemqt_t {
    item_t it; /* the item */
    weight qt; /* its quantity in the result */
    itemix ix; /* its index in the ordered items list */

    inline ukp_itemqt_t(const item_t &it, const size_t qt, const size_t ix) : it(it), qt(qt), ix(ix) {}

    void print(std::ostream &cout = std::cout) const {
      cout << "ix: " << ix << " qt: " << qt << " w: " << it.w << " p: " << it.p << std::endl;
    }
  };

  struct ukp_solution_t {
    profit opt;
    weight y_opt;
    std::vector<ukp_itemqt_t> used_items;
    #ifdef HBM_PROFILE
    /* Time of each phase */
    double sort_time;
    double vector_alloc_time;
    double linear_comp_time;
    double phase1_time;
    double phase2_time;
    double total_time;
    /* Some data about instance */
    itemix n;
    weight c, w_min, w_max;
    /* Some data about structures manipulates by ukp5 */
    weight last_dy_non_zero_non_n;
    weight qt_non_skipped_ys;
    weight qt_gy_zeros;
    weight qt_inner_loop_executions;
    std::vector<itemix> qt_i_in_dy;
    std::vector<profit> g;
    std::vector<itemix> d;
    std::vector<itemix> non_skipped_d;
    #endif /* HBM_PROFILE */
    #if defined(HBM_CHECK_PERIODICITY) || defined(HBM_CHECK_PERIODICITY_FAST)
    weight last_y_value_outer_loop;
    #endif /* PERIODICITY */
  };

  struct ukp_read_error : std::runtime_error {
    explicit ukp_read_error (const std::string &s) noexcept : std::runtime_error(s) {};
    explicit ukp_read_error (const char* s) noexcept : runtime_error(s) {};
  };

  void read_ukp_instance(std::istream &in, ukp_instance_t &ukpi);
  void read_sukp_instance(std::istream &in, ukp_instance_t &ukpi);

  void write_sukp_instance(std::ostream &out, ukp_instance_t &ukpi);

  void sort_by_eff(std::vector<item_t> &items);
}

/* Use an optimized swap for the item class (improves sorting time),
 * but only if this is possible (all members have the operator ^= defined) */
#if (defined(HBM_TWO_MULT_COMP) || defined(HBM_INT_EFF)) && !defined(HBM_NO_XOR_SWAP)
#define HBM_XORSWAP(a, b) ((a)^=(b),(b)^=(a),(a)^=(b))
namespace std {
  template <>
  inline void swap(hbm::item_t& a, hbm::item_t& b) noexcept
  {
    HBM_XORSWAP(a.w, b.w);
    HBM_XORSWAP(a.p, b.p);
    #ifdef HBM_INT_EFF
    HBM_XORSWAP(a.eff, b.eff);
    #endif
  }
}
#endif //XOR_SWAP SPECIALIZATION

#endif //HBM_UKP_COMMON_HPP
