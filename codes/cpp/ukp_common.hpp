#ifndef __UKP_COMMON_HPP_
#define __UKP_COMMON_HPP_

#define FP_EFF_TYPE float
#define INT_EFF_TYPE size_t
#define RATIONAL_EFF_TYPE boost::rational<size_t>

#ifndef RATIONAL_EFF
  #ifndef FP_EFF
    #ifndef INT_EFF
      #ifndef TWO_MULT_COMP
        #error ONE OF: RATIONAL_EFF, FP_EFF, INT_EFF OR TWO_MULT_COMP MUST BE \
               DEFINED, THE RECOMMENDED IS TWO_MULT_COMP
      #endif
    #else /*INT_EFF*/
      #ifdef TWO_MULT_COMP
        #error INT_EFF AND TWO_MULT_COMP CANNOT BE DEFINED AT THE SAME TIME
      #endif
      #define EFF_TYPE INT_EFF_TYPE
    #endif /*INT_EFF*/
  #else /*FP_EFF*/
    #ifdef TWO_MULT_COMP
      #error FP_EFF AND TWO_MULT_COMP CANNOT BE DEFINED AT THE SAME TIME
    #endif
    #ifdef INT_EFF
      #error FP_EFF AND INT_EFF CANNOT BE DEFINED AT THE SAME TIME
    #endif
    #define EFF_TYPE FP_EFF_TYPE
  #endif /*FP_EFF*/
#else /*RATIONAL_EFF*/
  #ifdef TWO_MULT_COMP
    #error RATIONAL_EFF AND TWO_MULT_COMP CANNOT BE DEFINED AT THE SAME TIME
  #endif
  #ifdef INT_EFF
    #error RATIONAL_EFF AND INT_EFF CANNOT BE DEFINED AT THE SAME TIME
  #endif
  #ifdef FP_EFF
    #error RATIONAL_EFF AND FP_EFF CANNOT BE DEFINED AT THE SAME TIME
  #endif
  #define EFF_TYPE RATIONAL_EFF_TYPE
#endif /*RATIONAL_EFF*/

#include <vector>
#include <istream>
#include <stdexcept> /* for runtime_error */

#if defined(TWO_MULT_COMP) || defined(INT_EFF)
#include <utility> /* to specialize swap */
#endif

#ifdef RATIONAL_EFF
#include <boost/rational.hpp>
#endif

struct item_t {
  size_t w;
  size_t p;
  #if defined(RATIONAL_EFF) || defined(INT_EFF) || defined(FP_EFF)
  EFF_TYPE efficiency;
  #endif

  inline item_t(void) {}
  #ifdef TWO_MULT_COMP
  inline item_t(const size_t &w, const size_t &p) : w(w), p(p) {}
  #elif defined(RATIONAL_EFF)
  inline item_t(const size_t &w, const size_t &p) : w(w), p(p), efficiency(p, w) {}
  #elif defined(INT_EFF)
  inline item_t(size_t w, size_t p) : w(w), p(p) {
    efficiency = (p << 32) / w;
  }
  #elif defined(FP_EFF)
  inline item_t(size_t w, size_t p) : w(w), p(p) {
    efficiency = ((EFF_TYPE)p) / ((EFF_TYPE)w);
  }
  #endif

  inline bool operator==(const item_t &o) const {
    return p == o.p && w == o.w;
  }

  /* Sort by non-increasing efficiency, if the efficiences are equal
   * sort by non-decreasing weight
   */
  #if defined(TWO_MULT_COMP)
  inline bool operator<(const item_t &o) const {
    size_t a = p * o.w, b = o.p * w;
    return  a > b || (a == b && w < o.w); 
  }
  #elif defined(RATIONAL_EFF) || defined(FP_EFF) || defined(INT_EFF)
  inline bool operator<(const item_t &o) const {
    return efficiency > o.efficiency || (efficiency == o.efficiency && w < o.w);
  }
  #endif

  #ifdef INT_EFF
  inline size_t operator>>(const int s) const {
    /* As we sort by non-increasing efficiency, we cannot return
     * (efficiency >> s), because the integer_sort will compare the value
     * trying to sort them in non-decreasing order, as (a < b) iff (~a >= ~b)
     * this solves the problem. */
    return (~efficiency) >> s;
  }
  #endif
};

#if (defined(TWO_MULT_COMP) || defined(INT_EFF)) && !NO_XOR_SWAP
#define XORSWAP(a, b) ((a)^=(b),(b)^=(a),(a)^=(b))
namespace std {
  template <>
  inline void swap(item_t& a, item_t& b) noexcept
  {
    XORSWAP(a.w, b.w);
    XORSWAP(a.p, b.p);
    #ifdef INT_EFF
    XORSWAP(a.efficiency, b.efficiency);
    #endif
  }
}
#endif

struct ukp_instance_t {
  size_t c;
  std::vector<item_t> items;
};

struct ukp_itemqt_t {
  item_t it;
  size_t qt;

  inline ukp_itemqt_t(void) {}
  inline ukp_itemqt_t(const item_t &it, const size_t &qt) : it(it), qt(qt) {}
};

struct ukp_solution_t {
  size_t opt;
  size_t y_opt;
  std::vector<ukp_itemqt_t> used_items;
  #ifdef PROFILE
  /* Time of each phase */
  double sort_time;
  double vector_alloc_time;
  double linear_comp_time;
  double phase1_time;
  double phase2_time;
  double total_time;
  /* Some data about instance */
  size_t c, n, w_min, w_max;
  /* Some data about structures manipulates by ukp5 */
  size_t last_dy_non_zero_non_n;
  size_t qt_non_skipped_ys;
  size_t qt_gy_zeros;
  size_t qt_inner_loop_executions;
  std::vector<size_t> qt_i_in_dy;
  std::vector<size_t> g;
  std::vector<size_t> d;
  std::vector<size_t> non_skipped_d;
  #endif /* PROFILE */
  #if defined(CHECK_PERIODICITY) || defined(CHECK_PERIODICITY_FAST)
  size_t last_y_value_outer_loop;
  #endif /* PERIODICITY */
};

struct ukp_read_error : std::runtime_error {
  explicit ukp_read_error (const std::string &s) noexcept : std::runtime_error(s) {};
  explicit ukp_read_error (const char* s) noexcept : runtime_error(s) {};
};

void read_ukp_instance(std::istream &in, ukp_instance_t &ukpi);
void read_sukp_instance(std::istream &in, ukp_instance_t &ukpi);

void write_sukp_instance(std::ostream &out, ukp_instance_t &ukpi);

void sort_by_efficiency(std::vector<item_t> &items);
#endif //__UKP_COMMON_HPP_ 

