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
    inline P u3(
      const item_t<W, P> &bi,
      const item_t<W, P> &bi2,
      const item_t<W, P> &bi3,
      const W &c
    ) {
      // The numbers at the right are the equation reference number on
      // "Knapsack Problems" from S. Martello and P. Toth (section 3.6).
      // Auxiliar values
      W c_ = c % bi.w;                              // 3.18
      W c_prime = c_ % bi2.w;                       // 3.22
      P z_prime = (c/bi.w)*bi.p + (c_/bi2.w)*bi2.p; // 3.21

      // u0 computation
      P u0 = z_prime + (c_prime*bi3.p)/bi3.w; // 3.23

      // u1_ computation (all bellow is 3.25)
      // We use the remainder of bi.w (rem) to simulate the
      // rounding-up of a result used in the article.
      W w2_less_c_prime = bi2.w - c_prime;
      W rem = w2_less_c_prime % bi1.w;
      W quo = w2_less_c_prime / bi1.w;
      W res = rem ? quo + 1 : quo;

      // On the equation immediatly below, p2/w2 is shown as a fraction
      // multiplying the rest of the equation. So there's no round down in the
      // equation but there's one round down here (assuming that W and P are
      // integers). This is safe to do because this value is subtracted by an
      // already rounded value (integer) and then the result is rounded down
      // (would have thrown away the truncated value anyway).
      P l = ((c_prime + res*bi.w1)*bi2.p)/bi2.w;
      P u1_ = z_prime + (l - res*bi.p);

      // u3 computation
      return max(u1_, u0); // 3.26
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
      const auto &items = ukpi.items;

      if (!already_sorted) {
        sort_by_eff(ukpi.items);
      }

      // Put items on article's notation
      vector<W> w;
      w.reserve(n + 2);
      vector<P> p;
      p.reserve(n + 2);
      for (auto &it : items) {
        w.push_back(it.w);
        p.push_back(it.p);
      }

      // We need to initialize some variables before we start jumping.
      W y;
      P u;

      // 1. Initialize
      step_1:
      P z = 0, z_ = 0;
      W c_ = c;

      // TODO: At the end verify if this can be done without adding this
      // dummy item.
      w.push_back(numeric_limits<W>::max());
      p.push_back(0);

      vector<W> x_(n+1, 0);
      P upper_u = u3(items[0], items[1], items[2], c);

      myvector<W> m;
      m.resize(n+1); // This resize don't initialize the m contents with zero

      m[n] = w[n + 1]; // Changing the dummy item would change there
      for (I k = n - 1; k > 0; --k) m[k] = min(m[k + 1], w[k + 1]);

      I j = 1;

      // 2. build a new current solution
      step_2:
      while (w[j] > c_)
        if (z >= z_ + (c_*p[j+1])/w[j+1]) goto step_5; else ++j;

      y = c_/w[j];
      u = ((c_ - y*w[j])*p[j+1])/w[j+1];

      if (z >= z_ + y*p[j] + u) goto step_5;
      if (u == 0) goto step_4;

      // 3. save the current solution
      step_3:
      c_ = c_ - y*w[j];
      z_ = z_ + y*p[j];
      x_[j] = y;
      j = j + 1;
      if (c_ >= m[j-1]) goto 2;
      if (z >= z_) goto 5;
      y = 0;

      // 4. update the best solution so far
      step_4:
      z = z_ + y*p[j];
      for (I k = 1; k < j; ++k) x[k] = x_[k];
      // STOPPED WORK HERE

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

