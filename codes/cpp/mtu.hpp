#ifndef HBM_MTU_HPP
#define HBM_MTU_HPP

#include "ukp_common.hpp"
#include "wrapper.hpp"

#include <vector>
#include <chrono>
#include <boost/rational.hpp>

#ifndef HBM_PROFILE_PRECISION
  #define HBM_PROFILE_PRECISION 5
#endif

#ifndef HBM_PRINT_VAR
  #define HBM_PRINT_VAR(var) out << #var ": " << var << std::endl
#endif

namespace hbm {
  template <typename W, typename P, typename I>
  struct mtu1_extra_info_t : extra_info_t {
    #ifdef HBM_PROFILE
    double sort_time{0};        ///< Time used sorting items.
    double vector_alloc_time{0};///< Time used allocating vectors for DP.
    double linear_comp_time{0}; ///< Time used by linear time preprocessing.
    double dp_time{0};     ///< Time used creating partial solutions.
    double sol_time{0};    ///< Time used to assemble solution.
    double total_time{0};  ///< Time used by all the algorithm.
    #endif //HBM_PROFILE
    I n{0};      ///< Instance number of items.
    W c{0};      ///< Instance capacity.

    virtual std::string gen_info(void) {
      std::stringstream out("");

      HBM_PRINT_VAR(c);
      HBM_PRINT_VAR(n);

      #ifdef HBM_PROFILE
      const double sum_time = sort_time + vector_alloc_time +
        linear_comp_time + dp_time + sol_time;

      std::streamsize old_precision = out.precision(HBM_PROFILE_PRECISION);
      const int two_first_digits_and_period = 3;
      const int percent_size = two_first_digits_and_period+HBM_PROFILE_PRECISION;
      std::ios_base::fmtflags old_flags = out.setf(std::ios::fixed, std:: ios::floatfield);
      char old_fill = out.fill(' ');

      #ifndef HBM_PRINT_TIME
        #define HBM_PRINT_TIME(name, var)\
          out << name << " time: " << var << "s (";\
          out << std::setw(percent_size) << (var/total_time)*100.0;\
          out << "%)" << std::endl
      #endif

      HBM_PRINT_TIME("Sort", sort_time);
      HBM_PRINT_TIME("Vect", vector_alloc_time);
      HBM_PRINT_TIME("O(n)", linear_comp_time);
      HBM_PRINT_TIME("DP  ", dp_time);
      HBM_PRINT_TIME("sol ", sol_time);
      HBM_PRINT_TIME("Sum ", sum_time);
      HBM_PRINT_TIME("All ", total_time);

      out.fill(old_fill);
      out.setf(old_flags);
      out.precision(old_precision);
      #endif //HBM_PROFILE

      return out.str();
    }
  };

  template <typename W, typename P, typename I>
  struct mtu1_extra_info_t : mtu2_extra_info_t<W, P, I> {};

  namespace hbm_mtu_impl {
    using namespace std;
    using namespace std::chrono;
    using namespace boost::math;

    #ifdef HBM_CHECK_OVERFLOW
    // Overflow check for multiplying two numbers, built-in style. We are
    // considering that the A/B/C types can be signed or unsigned, but the
    // values of both variables 'a' and 'b' will be non-negative. Also, this
    // is slow, and only used if someone want to compile with overflow check
    // enabled and isn't on g++/icc (for those two compiler we use a built-in
    // compiler function).
    inline bool mul_overflow_check<A, B, C>(const A &a, const B &b, C &c) {
      c = a*b;
      C max = numeric_limits<C>::max();
      if (b != 0 && a > max / b) {
        return true;
      }
      return false;
    }
      // Stolen from: https://sourceforge.net/p/predef/wiki/Compilers/
      #if defined(__GNUC__)
        #if defined(__GNUC_PATCHLEVEL__)
          #define __GNUC_VERSION__ (__GNUC__ * 10000 \
                                    + __GNUC_MINOR__ * 100 \
                                    + __GNUC_PATCHLEVEL__)
        #else
          #define __GNUC_VERSION__ (__GNUC__ * 10000 \
                                    + __GNUC_MINOR__ * 100)
        #endif
      #endif

      #if __GNUC_VERSION__ >= 50000
        // Will save the multiplication of 'a' and 'b' in 'c', and return
        // non-zero on the case of an overflow (and zero otherwise).
        #undef HBM_CHECK_OVERFLOW
        #define HBM_CHECK_OVERFLOW(a, b, c) __builtin_mull_overflow(a, b, &c)
      #elif
        // Same working as described the alternative version described above.
        #define HBM_CHECK_OVERFLOW(a, b, c) mul_overflow_check(a, b, c)
      #endif
    #endif //HBM_CHECK_OVERFLOW

    template<typename W, typename P>
    inline P u0(const item_t<W, P> &bi) {
      /* Overflow will not be tested here. We will test it before MTU1/MTU2
       * begins (one time per call) and then decide what to do. Theorically
       * if the multiplication between the best item profit by the capacity
       * don't overflows then no overflow will occur on MTU1/2.
       * #ifdef HBM_CHECK_OVERFLOW
      P ret;
      if (HBM_CHECK_OVERFLOW(bi.w, bi.p, ret)) {
      } else {
        cout << "The computation of the u0 bound has overflowed." << endl;
        cout << "Overflow values: bi.w = " << bi.w << ", bi.p = " << bi.p
             << endl;
        
      }
      #elif*/
      //return (c*p)/;
    }

    inline P u1_(const vector<W> &w, const vector<P> &p) {
      //return ;
    }

    inline P u3(const vector<W> &w, const vector<P> &p) {
      return max(u1_, u0);
    }

    // TODO: quick comment on the method
    template<typename W, typename P, typename I>
    void mtu1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      mtu1_extra_info_t<W, P, I>* eip = new mtu1_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_mtu2_begin = steady_clock::now();

      const W c = eip->c = ukpi.c;
      const I n = eip->n = ukpi.items.size();

      P z = 0, z_ = 0;
      W c_ = c;
      vector<W> w;
      w.reserve(n + 2);
      vector<P> p;
      p.reserve(n + 2);
      for (auto &it : items) {
        w.push_back(it.w);
        p.push_back(it.p);
      }
      w.push_back(numeric_limits<W>::max());
      p.push_back(0);

      vector<W> x(n+1, 0);

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      // METHOD CODE

      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_mtu2_begin);
      #endif
    }

    // TODO: quick comment on the method
    template<typename W, typename P, typename I>
    void mtu2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      mtu2_extra_info_t<W, P, I>* eip = new mtu2_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_mtu2_begin = steady_clock::now();

      const W c = ukpi.c;
      const I n = ukpi.items.size();
      I k;
      // v is the size the core problem begins and grows at each iteration.
      // TODO: MAKE THIS A PARAMETER
      I v = 100;
      do {
        k += min(k + v, n);

      } while ();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      // METHOD CODE

      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_mtu2_begin);
      #endif
    }

    // -------------------- WRAPPERS --------------------
    template<typename W, typename P, typename I>
    struct mtu1_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_mtu1_impl::mtu1(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "mtu1";
        return name;
      }
    };

    template<typename W, typename P, typename I>
    struct mtu2_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_mtu_impl::mtu2(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "mtu2";
        return name;
      }
    };
  }
  // -------------------- EXTERNAL FUNCTIONS --------------------

  /// Solves an UKP instance by the 
  /// TODO: ALGORITH_NAME presented at PAPER_NAME
  /// , and stores the results at sol.
  ///
  /// @note TODO: it only work with integers? there's some other caveat?
  /// @param ukpi The UKP instance to be solved.
  /// @param sol The object where the results will be written.
  /// @param already_sorted If the ukpi.items vector needs to be sorted by
  ///   TODO: what kind of ordering?
  template<typename W, typename P, typename I>
  void mtu1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_mtu_impl::mtu1(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items will NOT be sorted by non-increasing efficiency. If
  /// it's ommited the ukpi.items will be sorted by non-increasing efficiency.
  /// If more than one item have the same efficiency, then those items will be
  /// sorted by weight (in relation to each other).
  ///
  /// @see main_take_path
  /// @see mtu2(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void mtu1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_mtu_impl::mtu1(ukpi, sol, argc, argv);
  }

  /// Solves an UKP instance by the 
  /// TODO: ALGORITH_NAME presented at PAPER_NAME
  /// , and stores the results at sol.
  ///
  /// @note TODO: it only work with integers? there's some other caveat?
  /// @param ukpi The UKP instance to be solved.
  /// @param sol The object where the results will be written.
  /// @param already_sorted If the ukpi.items vector needs to be sorted by
  ///   TODO: what kind of ordering?
  template<typename W, typename P, typename I>
  void mtu2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_mtu_impl::mtu2(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items will NOT be sorted by non-increasing efficiency. If
  /// it's ommited the ukpi.items will be sorted by non-increasing efficiency.
  /// If more than one item have the same efficiency, then those items will be
  /// sorted by weight (in relation to each other).
  ///
  /// @see main_take_path
  /// @see mtu2(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void mtu2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_mtu_impl::mtu2(ukpi, sol, argc, argv);
  }
}

#endif //HBM_MTU_HPP

