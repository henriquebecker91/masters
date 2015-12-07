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

/* Avoids unnecessary include */
#ifdef HBM_RATIONAL_EFF
  #include <boost/rational.hpp>
#endif

/* Basic types */
namespace hbm {
  /* The code makes some assumptions about the profit and weight types,
   * if you break these the code probably won't compile or work as
   * intended.
   * 0) The profit and weight types behave as primitive numbers in c++
   *    would behave (define the same operators, you can output them
   *    with an ostream, etc...).
   * 1) For UKP5 (dynamic programming) the weight must always be an
   *    integral number (it is used to index arrays).
   * 2) The int_eff should be an integral type, fp_eff should be a
   *    floating point type.
   * 3) If profit or weight isn't integral, the macro HBM_NO_XOR_SWAP
   *    should be defined (it avoids applying ^= over those datatypes).
   * 3) If profit and weight aren't the same type:
   * 3.a) You probably shouldn't use HBM_RATIONAL_EFF.
   * 3.b) You should understand that any operation between a value of
   *      type profit and one of type weight (like multiplication, or
   *      division) will convert the weight type to the profit type
   *      before performing the operation.
   */
  typedef size_t profit;
  typedef size_t weight;

  /* This type is used to loop between the items and to store the 
   * index of items. The only assumption made is that it's integral
   * and sufficient large to contain the 'n' value (number of items).
   * The quantity of items in a solution is of type weight, as it is
   * the result of a capacity value divided by a item weight (both capacity
   * and item weight are of the weight type). And if the item weight is 1
   * (one) the the quantity is equal to the capacity. */
  typedef size_t itemix;
  
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

