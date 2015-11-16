#ifndef HBM_UKP_TYPES_HPP
#define HBM_UKP_TYPES_HPP

/* See the README for the explanation of every macro and basic type */

/* Efficiency related definitions
 * Only one of the bellow can be defined at the same time */
#define HBM_TWO_MULT_COMP
//#define HBM_INT_EFF
//#define HBM_FP_EFF
//#define HBM_RATIONAL_EFF

/* Sort related definitions */
//#define HBM_NO_XOR_SWAP

/* Other definitions */
//#define HBM_PROFILE
//#define HBM_DUMP
//#define HBM_CHECK_PERIODICITY

/* Only avoiding unnecessary include */
#ifdef HBM_RATIONAL_EFF
  #include <boost/rational.hpp>
#endif

/* Basic types */
namespace hbm {
  typedef size_t profit;
  typedef size_t weight;
  
  #ifdef HBM_INT_EFF
  typedef weight int_eff;
  #endif
  #ifdef HBM_FP_EFF
  typedef float fp_eff;
  #endif
  #ifdef HBM_RATIONAL_EFF
  typedef boost::rational<weight> rational_eff;
  #endif
}

#endif // HBM_UKP_TYPES_HPP

