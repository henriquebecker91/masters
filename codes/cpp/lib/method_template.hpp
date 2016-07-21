// CHECKLIST:
// 1º) Use your editor search/replace to replace greendp with the method name.
// 2º) Use your editor search/replace to replace GREENDP with the method name
//    (all in caps).
// 3º) Implement your method (replace METHOD CODE with your code).
// 4º) Search every TODO, follow it, and delete it.
// 5º) Make a copy and change the relevant bits on a run_*.cpp and maybe the
//    test_*.cpp (basically execute steps 1 and 2 on them, including the
//    filename).
// 6º) Change the Makefile to compile your run_*.cpp and test_*.cpp too.
// 7º) Delete this comment.
// Notes: if you want to pass commandline flags you can remove the wrapper_t
//  and make your own function that receives argv/argc.
#ifndef HBM_GREENDP_HPP
#define HBM_GREENDP_HPP

#include "ukp_common.hpp"
#include "wrapper.hpp"
#include "type_name.hpp"

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
  struct greendp_extra_info_t : extra_info_t {
    // TODO: check if this struct is adequated or not. If you want to have
    // those profiling times you will need to put HBM_START_TIMER and
    // HBM_STOP_TIMER on your code (at METHOD CODE). Check the other codes
    // to see how it works. Also, the variables n and c need to be set.
    std::string algorithm_name{"greendp"};
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
      out << "type_W: " << hbm::type_name<W>::get() << std::endl;
      out << "type_P: " << hbm::type_name<P>::get() << std::endl;
      out << "type_I: " << hbm::type_name<I>::get() << std::endl;
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

  namespace hbm_greendp_impl {
    using namespace std;
    using namespace std::chrono;
    using namespace boost::math;

    // TODO: quick comment on the method
    template<typename W, typename P, typename I>
    void greendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      greendp_extra_info_t<W, P, I>* eip = new greendp_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_greendp_begin = steady_clock::now();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      // METHOD CODE

      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_greendp_begin);
      #endif
    }

    // -------------------- WRAPPERS --------------------
    template<typename W, typename P, typename I>
    struct greendp_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_greendp_impl::greendp(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "greendp";
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
  void greendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_greendp_impl::greendp(ukpi, sol, already_sorted);
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
  /// @see greendp(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void greendp(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_greendp_impl::greendp(ukpi, sol, argc, argv);
  }
}

#endif //HBM_GREENDP_HPP


