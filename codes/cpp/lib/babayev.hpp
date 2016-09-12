// CHECKLIST:
// 3º) Implement your method (replace METHOD CODE with your code).
// 4º) Search every TODO, follow it, and delete it.
// 5º) Make a copy and change the relevant bits on a run_*.cpp and maybe the
//    test_*.cpp (basically execute steps 1 and 2 on them, including the
//    filename).
// 6º) Change the Makefile to compile your run_*.cpp and test_*.cpp too.
// 7º) Delete this comment.
// Notes: if you want to pass commandline flags you can remove the wrapper_t
//  and make your own function that receives argv/argc.
#ifndef HBM_BABAYEV_HPP
#define HBM_BABAYEV_HPP

#include "ukp_common.hpp"
#include "wrapper.hpp"
#include "type_name.hpp"
#include "mtu.hpp" // babayev make use of the u3 upper bound

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
  struct babayev_extra_info_t : extra_info_t {
    // TODO: check if this struct is adequated or not. If you want to have
    // those profiling times you will need to put HBM_START_TIMER and
    // HBM_STOP_TIMER on your code (at METHOD CODE). Check the other codes
    // to see how it works. Also, the variables n and c need to be set.
    std::string algorithm_name{"babayev"};
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

      HBM_PRINT_VAR(algorithm_name);
      out <<  "git_head_at_compilation: "
          << HBM_GIT_HEAD_AT_COMPILATION << std::endl;
      out << "type_W: " << hbm::type_name<W>::get() << std::endl;
      out << "type_P: " << hbm::type_name<P>::get() << std::endl;
      out << "type_I: " << hbm::type_name<I>::get() << std::endl;
      HBM_PRINT_VAR(c);
      HBM_PRINT_VAR(n);

      #ifdef HBM_PROFILE
      out << "HBM_PROFILE: defined" << std::endl;
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
      #else //HBM_PROFILE
      out << "HBM_PROFILE: undefined" << std::endl;
      #endif //HBM_PROFILE

      return out.str();
    }
  };

  namespace hbm_babayev_impl {
    using namespace std;
    using namespace std::chrono;
    using namespace boost;

    // This methods isn't described exactly in the paper. By trial and error,
    // we determined it as a lower bound version of u3. Yet, at cost of an
    // scan over the items (O(n)) we can get a better lower bound in many
    // cases.
    template<typename W, typename P>
    P greedy_lower_bound(instance_t<W, P> &ukpi, bool already_sorted) {
      /*P z = 0;
      W c_ = ukpi.c;
      if (!already_sorted) sort_by_eff(ukpi.items);
      for (auto& i : ukpi.items) {
        auto qt = c_ / i.w;
        z += qt*i.p;
        c_ -= qt*i.w;
      }

      return z;*/
      auto &v = ukpi.items;
      W c = ukpi.c;
      W c1 = c % v[0].w;
      W c2 = c1 % v[1].w;
      P z = (c/v[0].w)*v[0].p + (c1/v[1].w)*v[1].p + (c2/v[2].w)*v[2].p;

      return z;
    }

    // TODO: quick comment on the method
    template<typename W, typename P, typename I>
    void babayev(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      babayev_extra_info_t<W, P, I>* eip = new babayev_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_babayev_begin = steady_clock::now();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      // NOTE: the paper's notation is used here, so 'c' isn't used to denote
      // the knapsack capacity (as in other algorithms of this repository) but
      // to denote the optimal profit value (or its lower/upper bound, or a
      // tentative value for it). The variable 'b' is used to denote the
      // knapsack capacity.
      const W b = ukpi.c;
      const I n = ukpi.items.size();
      // STEP 1: Determine upper_c -- an upper bound for c_star,
      //   upper_c >= c_star;

      // Get the three most efficient items. TODO: do this in a better way?
      P lower_c = 0, c = 0, upper_c = 0;
      if (!already_sorted) sort_by_eff(ukpi.items); // already sort for the greedy algorithm
      sol.y_opt = 0;
      sol.opt = upper_c = u3(ukpi.items[0], ukpi.items[1], ukpi.items[2], ukpi.c);
      lower_c = greedy_lower_bound(ukpi, true);
      // STEP 2: c := upper_c
      c = upper_c;

      // STEP 3: Compute beta and test consistency of Eq. 7;
      // TODO: On page four, it's said that M1 is more efficient on a
      // certain condition than M2, and vice-versa. Check the condition and
      // implement both methods after.
      rational<P> R(ukpi.items[0].p, ukpi.items[0].w); // defined at page 3

      // Computing w and v with M1, see first paragraph of "Procedure M1" (p.3)
      // for reference.

      // Eq. 18
      P k = 0;
      if (R.denominator() == 1) {
        k = R.numerator() + 1;
      } else {
        k = rational_cast<P>(R) + 1;
      }
      // Eq. 16 (multiplied by -1) and adding one for rounding purposes
      P v1 = - (((upper_c - 1) / k) + 1);
      // Eq. 17 (multiplied by -1) and adding one for rounding purposes
      P v2 = - (rational_cast<P>((b*R - lower_c - 1)/(k - R)) + 1);
      // We want the "smallest absolute value" between two negative numbers.
      P v = max(v1, v2);

      // Eq. 15
      W w = -v*k + 1;

      // Important: don't forget we need to check if they are coprime with gcd.
      // END OF PROCEDURE M1 (w and v defined)
      auto &out = cout;
      HBM_PRINT_VAR(v);
      HBM_PRINT_VAR(w);
      HBM_PRINT_VAR(k);
      HBM_PRINT_VAR(R);

      myvector<P> agg_eq;
      agg_eq.resize(n + 1);
      agg_eq[n] = w; // slack variable / dummy item
      P alpha1 = w;
      for (I i = 0; i < n; ++i) {
        auto it = ukpi.items[i];
        agg_eq[i] = it.p*v + it.w*w;
        if (alpha1 > agg_eq[i]) alpha1 = agg_eq[i];
        cout << "agg_eq[" << i << "]: " << agg_eq[i] << endl;
      }
      cout << "agg_eq[n]: " << w << endl;
      cout << "agg_eq[n+1]: " << w*b << v << "*c" << endl;
      HBM_PRINT_VAR(alpha1);

      // IMPORTANT: the lower bound lower_c on their computational results
      // is given by solving the reduced problem (1% most efficient) exactly
      // (page 6).
      // IMPORTANT: two approachs for solving the consistency problem are
      // presented. On page 6, it's said that approach two was used to
      // generate the results presented.

      vector<P> r(alpha1, 0);

      // The algorithm isn't fully specified on the paper "A New Knapsack
      // Approach by Integer Equivalent Aggregation and Consistency
      // Determination", it references the paper "Integer Programming over a
      // Finite Additive Group" for core of the algorithm (the consistency
      // testing). We know the algorithm is right until there because the
      // values match the ones given at the the paper ("A New [...]"). Using
      // GLPK and a simple mathprog model we can verify that a consistency
      // check as described on the paper really would yield the results
      // presented. I abandoned this implementation task because there wasn't
      // enough time before ending master's.
      //
      // param beta, integer, > 0;
      //
      // param n, integer, > 0;
      // /* number of items */
      //
      // set I := 1..n;
      // /* set of items */
      //
      // param a{i in 1..n}, integer, >= 0;
      // /* the aggregate value between weight and profit */
      //
      // var x{i in I}, integer, >= 0;
      // /* x[i] = m, means m copies of item i are in the knapsack */
      //
      // minimize obj: sum{i in I} a[i]*x[i];
      // /* objective is to maximize the profit, but for this we minimize the aggregate function value */
      //
      // s.t. the_only_constraint: sum{i in I} a[i] * x[i] = beta;
      //
      // solve;
      //
      // for {i in I} { printf 'Item nº %d: %d\n', i, x[i]; }
      // printf "Result: %d\n", sum {i in I} a[i]*x[i];
      //
      // data;
      //
      // param beta := 229; /* will fail (unfeasible) */
      // /* param beta := 267;*/ /* will find one solution */
      // param n := 6;
      // param a := 1 48, 2 43, 3 127, 4 41, 5 121, 6 115;
      //
      // end;

      // STEP 4: If Eq. 7 is inconsistent, then set c := c - 1 and return to
      //  step 3.
      // STEP 5: If Eq. 7 is consistent, then c_star = c and the corresponding
      //   equation of Eq. 7 is a solution to the original knapsack problem
      //   (1)--(3). End.

      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_babayev_begin);
      #endif
    }

    // -------------------- WRAPPERS --------------------
    template<typename W, typename P, typename I>
    struct babayev_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_babayev_impl::babayev(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "babayev";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void babayev(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol,
                int argc, argv_t argv) {
      simple_wrapper(babayev_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }
  }
  // -------------------- EXTERNAL FUNCTIONS --------------------

  /// Solves an UKP instance by the unnamed algorithm presented at
  /// "A New Knapsack Approach by Integer Equivalent Aggregation and
  /// Consistency Determination", and stores the results at sol. The
  /// name 'babayev' is the surname of the paper's first author.
  ///
  /// @note Unfinished implementation. Currently waiting a response from
  ///   Babayev and Jennifer Ryan/Sanchez.
  /// @param ukpi The UKP instance to be solved.
  /// @param sol The object where the results will be written.
  /// @param already_sorted If the ukpi.items vector needs to be sorted by
  ///   efficiency.
  template<typename W, typename P, typename I>
  void babayev(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_babayev_impl::babayev(ukpi, sol, already_sorted);
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
  /// @see babayev(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void babayev(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_babayev_impl::babayev(ukpi, sol, argc, argv);
  }
}

#endif //HBM_BABAYEV_HPP

