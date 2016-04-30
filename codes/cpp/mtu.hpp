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
  struct mtu2_extra_info_t : mtu1_extra_info_t<W, P, I> {};

  namespace hbm_mtu_impl {
    using namespace std;
    using namespace std::chrono;
    using namespace boost;

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

    // "Knapsack Problems" p. 93, equation 3.19 
    template<typename W, typename P>
    inline P u1(
      const item_t<W, P> &bi,
      const item_t<W, P> &bi2,
      const W &c
    ) {
      W c_ = c % bi.w;
      return (c/bi.w)*bi.p + (c_*bi2.p)/bi2.w;
    }

    template<typename W, typename P>
    inline P u3(
      const item_t<W, P> &bi,
      const item_t<W, P> &bi2,
      const item_t<W, P> &bi3,
      const W &c
    ) {
      // The numbers at the right are the equations reference number on
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
      W rem = w2_less_c_prime % bi.w;
      W quo = w2_less_c_prime / bi.w;
      W res = rem ? quo + 1 : quo;

      // On the equation immediatly below, at the book, p2/w2 is shown as a
      // fraction multiplying the rest of the equation. So there's no round
      // down in the equation but there's one round down here (assuming that
      // W and P are integers). This is safe to do because this value is
      // subtracted by an already rounded value (integer) and then the result
      // is rounded down (would have thrown away the smaller than on value
      // anyway).
      P l = ((c_prime + res*bi.w)*bi2.p)/bi2.w;
      P u1_ = z_prime + (l - res*bi.p);

      // u3 computation
      return max(u1_, u0); // 3.26
    }

    // TODO: write short note
    // Note that this method has index notation beggining on 1.
    // It's expected that the size of w and p exceeds n by at least one
    // position (i.e. w and p have at least size n+1).
    template<typename W, typename P, typename I>
    void inner_mtu1(
      myvector<W> &w,
      myvector<P> &p,
      const I &n,
      const W &c,
      P &z,
      myvector<W> &x
    ) {
      const item_t<W, P> bi(w[1], p[1]), bi2(w[2], p[2]), bi3(w[3], p[3]);

      W y, i;
      P u;
      I h;

      // 1. Initialize
      //step_1: // unused step
      z = 0;
      P z_ = 0;
      W c_ = c;

      // TODO: At the end verify if this can be done without adding this
      // dummy item.
      p[n+1] = 0;
      w[n+1] = numeric_limits<W>::max();

      vector<W> x_(n+1, 0);
      P upper_u = u3(bi, bi2, bi3, c);

      myvector<W> m;
      m.resize(n+1); // This resize don't initialize the m contents with zero

      m[n] = w[n + 1]; // Changing the dummy item would change there
      for (I k = n - 1; k > 0; --k) m[k] = min(m[k + 1], w[k + 1]);

      I j = 1;

      // 2. build a new current solution
      step_2:
      while (w[j] > c_)
        if (z >= z_ + (c_*p[j+1])/w[j+1]) goto step_5; else j = j + 1;

      y = c_/w[j];
      u = ((c_ - y*w[j])*p[j+1])/w[j+1];

      if (z >= z_ + y*p[j] + u) goto step_5;
      if (u == 0) goto step_4;

      // 3. save the current solution
      //step_3: // unused step
      c_ = c_ - y*w[j];
      z_ = z_ + y*p[j];
      x_[j] = y;
      j = j + 1;
      if (c_ >= m[j-1]) goto step_2;
      if (z >= z_) goto step_5;
      y = 0;

      // 4. update the best solution so far
      step_4:
      z = z_ + y*p[j];
      for (I k = 1; k < j; ++k) x[k] = x_[k];
      x[j] = y;
      for (I k = j + 1; k <= n; ++k) x[k] = 0;
      if (z == upper_u) goto stop;

      // 5. backtrack
      step_5:
      i = 0; // We use zero as an sentinel value
      for (I k = j - 1; k > 0; --k) if (x_[k] > 0) { i = k; break; }
      if (i == 0) goto stop; // "if no such i then return"
      c_ = c_ + w[i];
      z_ = z_ - p[i];
      x_[i] = x_[i] - 1;
      if (z >= z_ + (c_*p[i+1])/w[i+1]) {
        c_ = c_ + w[i]*x_[i];
        z_ = z_ - p[i]*x_[i];
        x_[i] = 0;
        j = i;
        goto step_5;
      }
      j = i + 1;
      if (c_ - w[i] >= m[i]) goto step_2;
      h = i;

      // 6. try to replace one item of type i with items of type h
      step_6:
      h = h + 1;
      if (z >= z_ + (c_*p[h])/w[h]) goto step_5;
      if (w[h] == w[i]) goto step_6;
      if (w[h] > w[i] ) {
        if (w[h] > c_ || z >= z_ + p[h]) goto step_6;
        z = z_ + p[h];
        for (I k = 1; k <= n; ++k) x[k] = x_[k];
        x[h] = 1;
        if (z == upper_u) goto stop;
        i = h;
        goto step_6;
      } else {
        if (c_ - w[h] < m[h-1]) goto step_6;
        j = h;
        goto step_2;
      }

      stop:;
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
      steady_clock::time_point all_mtu1_begin = steady_clock::now();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      const W c = eip->c = ukpi.c;
      const I n = eip->n = ukpi.items.size();
      const auto &items = ukpi.items;

      if (!already_sorted) {
        sort_by_eff(ukpi.items);
      }

      // Put items on article's notation
      myvector<W> w;
      w.resize(n + 2);
      w[0] = 0; // position does not exist, notation begins at 1
      myvector<P> p;
      p.resize(n + 2);
      p[0] = 0; // position does not exist, notation begins at 1
      item_t<W, P> it;
      for (I i = 0, j = 1; i < n; ++i, ++j) {
        it = items[i];
        w[j] = it.w;
        p[j] = it.p;
      }

      // We need to initialize some variables before we start jumping.
      myvector<W> x;
      x.resize(n+1); // This resize don't initialize the m contents with zero

      P z;
      inner_mtu1(w, p, n, c, z, x);
      /*W y, i;
      P u;
      I h;

      // 1. Initialize
      //step_1: // unused step
      P z = 0, z_ = 0;
      W c_ = c;

      // TODO: At the end verify if this can be done without adding this
      // dummy item.
      p.push_back(0);
      w.push_back(numeric_limits<W>::max());

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
        if (z >= z_ + (c_*p[j+1])/w[j+1]) goto step_5; else j = j + 1;

      y = c_/w[j];
      u = ((c_ - y*w[j])*p[j+1])/w[j+1];

      if (z >= z_ + y*p[j] + u) goto step_5;
      if (u == 0) goto step_4;

      // 3. save the current solution
      //step_3: // unused step
      c_ = c_ - y*w[j];
      z_ = z_ + y*p[j];
      x_[j] = y;
      j = j + 1;
      if (c_ >= m[j-1]) goto step_2;
      if (z >= z_) goto step_5;
      y = 0;

      // 4. update the best solution so far
      step_4:
      z = z_ + y*p[j];
      for (I k = 1; k < j; ++k) x[k] = x_[k];
      x[j] = y;
      for (I k = j + 1; k <= n; ++k) x[k] = 0;
      if (z == upper_u) goto stop;

      // 5. backtrack
      step_5:
      i = 0; // We use zero as an sentinel value
      for (I k = j - 1; k > 0; --k) if (x_[k] > 0) { i = k; break; }
      if (i == 0) goto stop; // "if no such i then return"
      c_ = c_ + w[i];
      z_ = z_ - p[i];
      x_[i] = x_[i] - 1;
      if (z >= z_ + (c_*p[i+1])/w[i+1]) {
        c_ = c_ + w[i]*x_[i];
        z_ = z_ - p[i]*x_[i];
        x_[i] = 0;
        j = i;
        goto step_5;
      }
      j = i + 1;
      if (c_ - w[i] >= m[i]) goto step_2;
      h = i;

      // 6. try to replace one item of type i with items of type h
      step_6:
      h = h + 1;
      if (z >= z_ + (c_*p[h])/w[h]) goto step_5;
      if (w[h] == w[i]) goto step_6;
      if (w[h] > w[i] ) {
        if (w[h] > c_ || z >= z_ + p[h]) goto step_6;
        z = z_ + p[h];
        for (I k = 1; k <= n; ++k) x[k] = x_[k];
        x[h] = 1;
        if (z == upper_u) goto stop;
        i = h;
        goto step_6;
      } else {
        if (c_ - w[h] < m[h-1]) goto step_6;
        j = h;
        goto step_2;
      }*/

      stop:
      sol.opt = 0;
      sol.y_opt = 0;
      for (I k = 1; k <= n; ++k) {
        if (x[k] > 0) {
          sol.used_items.emplace_back(item_t<W, P>(w[k], p[k]), x[k], k);
          sol.opt += x[k]*p[k];
          sol.y_opt += x[k]*w[k];
        }
      }
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

      const W c = eip->c = ukpi.c;
      const I n = eip->n = ukpi.items.size();
      auto &items = ukpi.items;

      P z;
      // v is the size the core problem begins and grows at each iteration.
      // TODO: MAKE THIS A PARAMETER
      if (n < 100) {
        mtu1(ukpi, sol, already_sorted);
        return;
      }
      I v = max(100, n/100);

      if (!already_sorted) {
        // The G, E, and overlined E sets are unecessary. The code does the
        // almost the same as partial ordering the first 'v' elements.
        sort_k_most_eff(items.begin(), items.end(), v);
      }

      // The core will be kept a profit and a weight array, to avoid having
      // to convert for each call to inner_mtu1. The arrays are of size 'n + 2'
      // because: 1º) notation starts at 1; 2º) inner_mtu1 adds a dummy
      // item at the end.
      myvector<W> core_p;
      core_p.resize(n + 2); // this resize don't initialize any values
      core_p[0] = 0; // position does not exist, notation begins at 1
      myvector<P> core_w;
      core_w.resize(n + 2); // this resize don't initialize any values
      core_w[0] = 0; // position does not exist, notation begins at 1
      item_t<W, P> aux_item;
      for (I i = 0, j = 1; i < v; ++i, ++j) {
        aux_item = items[i];
        core_w[j] = aux_item.w;
        core_p[j] = aux_item.p;
      }

      // MTU2 describe some steps on a high level. We tried to implemet those
      // on the most efficient fashion possible. One of the high-level steps
      // that are hard to manage in an efficient fashion is keeping non_core
      // update, it have elements removed frequently (from front and any
      // position in the middle of it).
      myvector< item_t<W, P> > non_core;
      non_core.resize(n - v);
      const size_t qt_bytes_to_copy = non_core.size()*sizeof(item_t<W, P>);
      memcpy(non_core.data(), items.data() + v, qt_bytes_to_copy);
      // The vector<bool> isn't ideal here, as it trades memory consumption for
      // time. And we are optimizing this algorithm for time.
      myvector<char> items_to_remove(non_core.size(), 0);
      size_t qts_item_to_remove = 0;
      myvector< item_t<W, P> > aux_non_core;

      // The solution vector.
      myvector<W> x;
      x.resize(n + 1);
      memset(x.data(), 0, (v+1)*sizeof(W));

      inner_mtu1(core_w, core_p, v, c, z, x);

      P upper_u3 = u3(items[0], items[1], items[2], c);

      I non_core_start = 0;
      size_t non_core_size = non_core.size();

      I k = v;
      while (z != upper_u3 && non_core.size() > 0) {
        P u, pj;
        W wj;
        for (size_t j = non_core_start; j < non_core_size; ++j) {
          pj = non_core[j].p;
          wj = non_core[j].w;
          u = pj + u1(items[0], items[1], c - wj);
          if (u > z) u = pj + u3(items[0], items[1], items[2], c - wj);
          // the 'if' below can't be an 'else' of the 'if' above, don't try
          // the variable 'u' changes value on the 'if' above
          if (u <= z) {
            items_to_remove[j] = true;
            ++qt_items_to_remove;
          }
        }

        // Here we remove the dominated items from non_core (together with
        // the v first items that were added to core, and that are not part
        // of non_core anymore), and update the relevant variables.
        {
          aux_non_core.resize(non_core_size - qt_items_to_remove);

          for (size_t j = non_core_start; j < non_core_size; ++j) {
            if (!items_to_remove[j]) { 
              aux_non_core[j - non_core_start] = non_core[j];
            }
          }

          non_core.swap(aux_non_core);
          non_core_start = 0;
          non_core_size = non_core.size();
          memset(items_to_remove.data(), 0, non_core_size*sizeof(char));
          qt_items_to_remove = 0;
        }

        // Now we get the v most efficient items on non_core, and add them to
        // the core. We do not remove those items of non_core, but we set
        // non_core_start that will remove them for us on the next loop.
        if (!already_sorted) {
          sort_k_most_eff(non_core.begin(), non_core.end(), v);
        }
        const I next_slice_size = min(v, non_core_size);
        for (I i = 0; i < next_slice_size; ++i, ++k) {
          aux_item = items[i];
          core_w[k] = aux_item.w;
          core_p[k] = aux_item.p;
        }
        non_core_start = v;

        // We do not need to clean the x solution array, inner_mtu1 already
        // does this.
        //memset(x.data(), 0, k*sizeof(W));
        z = 0;

        inner_mtu1(core_w, core_p, k, c, z, x);
      } 

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
        hbm_mtu_impl::mtu1(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "mtu1";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void mtu1(
      instance_t<W, P> &ukpi,
      solution_t<W, P, I> &sol,
      int argc,
      argv_t argv
    ) {
      simple_wrapper(mtu1_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }

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

    template<typename W, typename P, typename I = size_t>
    void mtu2(
      instance_t<W, P> &ukpi,
      solution_t<W, P, I> &sol,
      int argc,
      argv_t argv
    ) {
      simple_wrapper(mtu2_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }
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

