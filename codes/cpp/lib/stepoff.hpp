#ifndef HBM_STEPOFF_HPP
#define HBM_STEPOFF_HPP

#include "ukp_common.hpp"
#include "wrapper.hpp"
#include "type_name.hpp"

#include <chrono>
#include <iomanip>

#ifndef HBM_PROFILE_PRECISION
  #define HBM_PROFILE_PRECISION 5
#endif

#ifndef HBM_PRINT_VAR
  #define HBM_PRINT_VAR(var) out << #var ": " << var << std::endl
#endif

namespace hbm {
  template <typename W, typename P, typename I>
  struct ordered_step_off_extra_info_t : extra_info_t {
    std::string algorithm_name{"ordered_step_off"};
    #ifdef HBM_PROFILE
    double sort_time{0};        ///< Time used sorting items.
    double vector_alloc_time{0};///< Time used allocating vectors for DP.
    double linear_comp_time{0}; ///< Time used by linear time preprocessing.
    double dp_time{0};     ///< Time used creating partial solutions.
    double sol_time{0};    ///< Time used to assemble solution.
    double total_time{0};  ///< Time used by all the algorithm.
    #endif //HBM_PROFILE
    I n{0};  ///< Instance number of items.
    W c{0};  ///< Instance capacity.

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

  template <typename W, typename P, typename I>
  struct terminating_step_off_extra_info_t : extra_info_t {
    std::string algorithm_name{"terminating_step_off"};
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
    W x0{0}; ///< Last capacity computed (the algorithm can terminate early). 

    virtual std::string gen_info(void) {
      std::stringstream out("");

      HBM_PRINT_VAR(algorithm_name);
      out <<  "git_head_at_compilation: "
          << HBM_GIT_HEAD_AT_COMPILATION << std::endl;
      out << "type_W: " << hbm::type_name<W>::get() << std::endl;
      out << "type_P: " << hbm::type_name<P>::get() << std::endl;
      out << "type_I: " << hbm::type_name<I>::get() << std::endl;
      HBM_PRINT_VAR(n);
      HBM_PRINT_VAR(c);
      HBM_PRINT_VAR(x0);

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

  namespace hbm_stepoff_impl {
    using namespace std;
    using namespace std::chrono;

    // Code based in the algorithm presented in page 1059 (15) of
    // "The theory and computation of knapsack functions" by Gilmore
    // and Gomory.
    template<typename W, typename P, typename I>
    void ordered_step_off(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      ordered_step_off_extra_info_t<W, P, I>* eip = new ordered_step_off_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_ordered_step_off_begin = steady_clock::now();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      // n == m, in the original article
      const I n = eip->n = static_cast<I>(ukpi.items.size());
      // c == L, in the original article
      const W c = eip->c = ukpi.c;
      myvector<W> w; // w == l, in the original article
      myvector<P> p; // p == pi (?), in the original article

      HBM_START_TIMER();
      //I: // not used
      if (!already_sorted) {
        // this sorts by non-increasing efficiency
        sort_by_eff(ukpi.items);
        // we need non-decreasing, so we reverse the array
        reverse(ukpi.items.begin(), ukpi.items.end());
      }
      HBM_STOP_TIMER(eip->sort_time);

      HBM_START_TIMER();
      // break the base-0 items array in two arrays base-1
      w.resize(n+1);
      p.resize(n+1);
      for (I i = 0, i_ = 1; i < n; ++i, ++i_) {
        auto it = ukpi.items[i];
        w[i_] = it.w;
        p[i_] = it.p;
      }
      HBM_STOP_TIMER(eip->linear_comp_time);

      HBM_START_TIMER();
      vector<P> F_star(c + 1, static_cast<P>(0));
      myvector<I> l_star;
      l_star.resize(c+1); // myvector does not initialize the memory in resize
      l_star[0] = static_cast<I>(1);
      W y = static_cast<W>(0); // y == x2, in the article
      I j;
      P V;
      HBM_STOP_TIMER(eip->vector_alloc_time);

      HBM_START_TIMER();
      //II:
      II_1:
      j = l_star[y];
      II_2:
      if (y + w[j] <= c) {
        V = p[j] + F_star[y];
        goto II_3;
      } else {
        goto II_4;
      }
      II_3:
      // The original algorithm has 'V >= ...' and no 'else if' in this step.
      // This small change makes the algorithm to be 20~40 times faster in
      // some instances. For example, the original version takes 24.0287
      // seconds or 17238430636 inner loop operations to solve the
      // instance sc_a-5n10000wmin110000-17-c9492344.ukp; the modified
      // version (including the 'else if') takes 0.869493 seconds ot
      // 608781965 inner loop operations to solve the same instance.
      if (V > F_star[y + w[j]]) {
        F_star[y + w[j]] = V;
        l_star[y + w[j]] = j;
      } else if (V == F_star[y + w[j]] && l_star[y + w[j]] < j) {
        l_star[y + w[j]] = j;
      }
      II_4:
      if (j < n) {
        j = j + 1;
        goto II_2;
      }

      //III:
      III_1:
      if (y < c) {
        y = y + 1;
        goto III_2;
      } else {
        goto stop;
      }
      III_2:
      if (F_star[y] > F_star[y - 1]) {
        goto II_1;
      } else {
        F_star[y] = F_star[y - 1];
        l_star[y] = n + 1;
        goto III_1;
      }

      stop:
      HBM_STOP_TIMER(eip->dp_time);
      HBM_START_TIMER();
      // store the solution in our format
      vector<I> qts_its(n, 0);
      I ly;
      while (y != 0 && l_star[y] == n + 1) --y;
      sol.y_opt = y;
      while (y != 0) {
        ly = l_star[y];
        y -= w[ly];
        ++qts_its[ly - 1];
      }

      sol.opt = 0;
      for (I i = 0; i < n; ++i) {
        if (qts_its[i] > 0) {
          auto it = ukpi.items[i];
          sol.used_items.emplace_back(it, qts_its[i], i);
          sol.opt += it.p * static_cast<P>(qts_its[i]);
        }
      }
      sol.used_items.shrink_to_fit();
      HBM_STOP_TIMER(eip->sol_time);

      // debug dump
      //auto &out = std::cout;
      //size_t num_ops = 0;
      //for (size_t y = 0; y <= c; ++y) {
        //cout << y << ";" << (l_star[y] == n + 1 ? 0 : n - l_star[y]) << endl;
        //num_ops += (n + 1) - l_star[y];
      //}
      //HBM_PRINT_VAR(num_ops);

      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_ordered_step_off_begin);
      #endif
    }

    // Code based in the algorithm presented in page 1059 (15) of
    // "The theory and computation of knapsack functions" (by Gilmore
    // and Gomory) and improved by the changes presented in 1061 (17).
    template<typename W, typename P, typename I>
    void terminating_step_off(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      terminating_step_off_extra_info_t<W, P, I>* eip = new terminating_step_off_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_terminating_step_off_begin = steady_clock::now();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      // n == m, in the original article
      const I n = eip->n = static_cast<I>(ukpi.items.size());
      // c == L, in the original article
      const W c = eip->c = ukpi.c;
      myvector<W> w; // w == l, in the original article
      myvector<P> p; // p == pi (?), in the original article
      myvector<W> lambda_; // max{l_i, i>= l*(x2)}, in the original article

      HBM_START_TIMER();
      //I: // not used
      if (!already_sorted) {
        // this sorts by non-increasing efficiency
        sort_by_eff(ukpi.items);
        // we need non-decreasing, so we reverse the array
        reverse(ukpi.items.begin(), ukpi.items.end());
      }
      HBM_STOP_TIMER(eip->sort_time);

      HBM_START_TIMER();
      // break the base-0 items array in two arrays base-1
      w.resize(n+1);
      p.resize(n+1);
      lambda_.resize(n+1);
      for (I i = 0, i_ = 1; i < n; ++i, ++i_) {
        auto it = ukpi.items[i];
        w[i_] = it.w;
        p[i_] = it.p;
      }
      lambda_[n] = w[n];
      for (I i = n-1; i > 0; --i) {
        lambda_[i] = max(lambda_[i+1], w[i]);
      }
      HBM_STOP_TIMER(eip->linear_comp_time);

      HBM_START_TIMER();
      vector<P> F_star(c + 1, static_cast<P>(0));
      myvector<I> l_star;
      l_star.resize(c+1); // myvector does not initialize the memory in resize
      l_star[0] = static_cast<I>(1);
      W y = static_cast<W>(0); // y == x2, in the article
      W x_star = 0;
      W lambda = lambda_[1];
      I j;
      P V;
      HBM_STOP_TIMER(eip->vector_alloc_time);

      HBM_START_TIMER();
      //II:
      II_1:
      j = l_star[y];
      II_2:
      if (y + w[j] <= c) {
        V = p[j] + F_star[y];
        goto II_3;
      } else {
        goto II_4;
      }
      II_3:
      // The original algorithm has 'V >= ...' and no 'else if' in this step.
      // This small change makes the algorithm to be 20~40 times faster in
      // some instances. For example, the original version takes 24.0287
      // seconds or 17238430636 inner loop operations to solve the
      // instance sc_a-5n10000wmin110000-17-c9492344.ukp; the modified
      // version (including the 'else if') takes 0.869493 seconds ot
      // 608781965 inner loop operations to solve the same instance.
      if (V > F_star[y + w[j]]) {
        F_star[y + w[j]] = V;
        l_star[y + w[j]] = j;
      } else if (V == F_star[y + w[j]] && l_star[y + w[j]] < j) {
        l_star[y + w[j]] = j;
      }
      II_4:
      if (j < n) {
        j = j + 1;
        goto II_2;
      }

      //III:
      III_1:
      if (y < x_star + lambda && y < c) {
        y = y + 1;
        goto III_2;
      } else {
        goto stop;
      }
      III_2:
      //III_3:
      if (F_star[y] > F_star[y - 1]) {
        // The 'if' below could be outside of the 'if' above, but this would
        // make necessary to initialize all of l_star with n+1 before starting
        // the algorithm. In the original paper this 'if' below was step III.2,
        // probably G&G initialized the l_star positions with n+1.
        if (l_star[y] < n) {
          x_star = y;
          lambda = lambda_[l_star[y]];
        }
        goto II_1;
      } else {
        F_star[y] = F_star[y - 1];
        l_star[y] = n + 1;
        goto III_1;
      }

      stop:
      HBM_STOP_TIMER(eip->dp_time);
      HBM_START_TIMER();
      // store the solution in our format
      vector<I> qts_its(n, 0);
      eip->x0 = x_star + lambda;
      if (y == eip->x0) { // if the algorithm terminated early
        qts_its[n-1] = ((c - y)/w[n]) + 1;
        y = c - qts_its[n-1]*w[n];
      }
      while (y != 0 && l_star[y] == n + 1) --y;
      sol.y_opt = y + qts_its[n-1]*w[n];
      while (y != 0) {
        I ly = l_star[y];
        y -= w[ly];
        ++qts_its[ly - 1]; // return solution to base zero
      }

      sol.opt = 0;
      for (I i = 0; i < n; ++i) {
        if (qts_its[i] > 0) {
          auto it = ukpi.items[i];
          sol.used_items.emplace_back(it, qts_its[i], i);
          sol.opt += it.p * static_cast<P>(qts_its[i]);
        }
      }
      sol.used_items.shrink_to_fit();
      HBM_STOP_TIMER(eip->sol_time);

      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_terminating_step_off_begin);
      #endif
    }

    // -------------------- WRAPPERS --------------------
    // ORDERED STEP-OFF
    template<typename W, typename P, typename I>
    struct ordered_step_off_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_stepoff_impl::ordered_step_off(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "ordered_step_off";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void ordered_step_off(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol,
                int argc, argv_t argv) {
      simple_wrapper(ordered_step_off_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }

    // TERMINATING STEP-OFF
    template<typename W, typename P, typename I>
    struct terminating_step_off_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_stepoff_impl::terminating_step_off(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "terminating_step_off";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void terminating_step_off(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol,
                int argc, argv_t argv) {
      simple_wrapper(terminating_step_off_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }
  }
  // -------------------- EXTERNAL FUNCTIONS --------------------

  /// Solves an UKP instance by the `ordered step-off method' presented in
  /// "The theory and computation of knapsack functions" by Gilmore and Gomory,
  /// and stores the results at sol.
  ///
  /// @note: The method only works with integer weights.
  /// @param ukpi The UKP instance to be solved.
  /// @param sol The object where the results will be written.
  /// @param already_sorted If the ukpi.items vector needs to be sorted by
  ///   non-decreasing efficiency. An if there is more than one most efficient
  ///   item, the one with the smallest profit (or weight) should be the first.
  template<typename W, typename P, typename I>
  void ordered_step_off(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_stepoff_impl::ordered_step_off(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items will NOT be sorted by non-decreasing efficiency. If
  /// it's ommited the ukpi.items will be sorted by non-increasing efficiency.
  /// If more than one item have the same efficiency, then those items will be
  /// sorted by weight (in relation to each other).
  ///
  /// @see main_take_path
  /// @see ordered_step_off(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void ordered_step_off(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_stepoff_impl::ordered_step_off(ukpi, sol, argc, argv);
  }

  /// Solves an UKP instance by the `ordered step-off method' presented in
  /// "The theory and computation of knapsack functions" by Gilmore and Gomory,
  /// and stores the results at sol.
  ///
  /// @note: The method only works with integer weights.
  /// @param ukpi The UKP instance to be solved.
  /// @param sol The object where the results will be written.
  /// @param already_sorted If the ukpi.items vector needs to be sorted by
  ///   non-decreasing efficiency. An if there is more than one most efficient
  ///   item, the one with the smallest profit (or weight) should be the first.
  template<typename W, typename P, typename I>
  void terminating_step_off(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_stepoff_impl::terminating_step_off(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items will NOT be sorted by non-decreasing efficiency. If
  /// it's ommited the ukpi.items will be sorted by non-increasing efficiency.
  /// If more than one item have the same efficiency, then those items will be
  /// sorted by weight (in relation to each other).
  ///
  /// @see main_take_path
  /// @see terminating_step_off(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void terminating_step_off(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_stepoff_impl::terminating_step_off(ukpi, sol, argc, argv);
  }
}

#endif //HBM_STEPOFF_HPP


