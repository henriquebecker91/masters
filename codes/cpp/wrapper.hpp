#ifndef HBM_WRAPPER_HPP
#define HBM_WRAPPER_HPP

#include "ukp_common.hpp"

namespace hbm {
  /// @brief Subclass this class to use simple_wrapper.
  template <typename W, typename P, typename I>
  struct wrapper_t {
    /// What your method do.
    virtual void operator()(instance_t<W, P> &, solution_t<W, P, I> &sol, bool) const {
      // To make easier to see that the operator was not redefined we set
      // some fields as 42. If this appear the programmmer will search
      // for 42 in the sources and find this line. We don't use a macro
      // for this magic number as we want the programmer to find this line.
      sol.opt = 42;
      sol.y_opt = 42;
    }
    /// Its name.
    virtual const std::string& name(void) const {
      static const std::string name = "<This should be a function name, "
        "but a programmer forgot to redefine the " + std::string(__func__)
        + " virtual method.>";
      return name;
    }
  };

  namespace hbm_wrapper_impl {
    using namespace std;

    template<typename W, typename P, typename I = size_t>
    void simple_wrapper(const wrapper_t<W, P, I> &wp, instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, char** argv) {
      static const string ALREADY_SORTED = "--already-sorted";
      if (argc == 0) {
        wp(ukpi, sol, false);
      } else if (argc >= 1) {
        if (argc > 1) {
          cerr << wp.name() << " (simple_wrapper): parameter warning: "
                  "Only one flag is allowed, but argc is " << argc
                  << ". " << endl;
        }

        if (ALREADY_SORTED == argv[0]) {
          wp(ukpi, sol, true);
        } else {
          cerr << wp.name() << " (simple_wrapper): parameter warning:"
                  " The only allowed flag is '" << ALREADY_SORTED <<
                  "'. The first flag received was \"" << argv[0] <<
                  "\". Executing the algorithm as no"
                  " flags were given. " << endl;
          wp(ukpi, sol, false);
        }
      } else {
        cerr << "simple_wrapper (error): Negative argc?! " <<
                wp.name() << " wasn't called. " << endl;
      }
    }
  }

  /// @brief Takes a wrapper_t and calls it to solve the ukpi and write
  ///   the results on sol. If the first parameter (argv[0]) is
  ///   "--already-sorted", then wp is called with true as the last
  ///   argument, otherwise the last argument will be false.
  ///
  /// It is necessary to stress that the procedure expect "--already-sorted"
  /// to be the argv[0], therefore the program name and previous parameters
  /// should be skipped first. The procedure writes warning and error
  /// messages on cerr.
  ///
  /// @param wp A subclass of wrapper_t with operator() and name redefined.
  /// @param ukpi An UKP instance.
  /// @param sol Where the results will be written.
  /// @param argc The number of arguments, expected to be 0 or 1.
  /// @param argv The arguments, expected to be ["--already-sorted"]
  ///   if argc is 1.
  template<typename W, typename P, typename I = size_t>
  void simple_wrapper(const wrapper_t<W, P, I> &wp,
                      instance_t<W, P> &ukpi,
                      solution_t<W, P, I> &sol,
                      int argc, char** argv) {
    hbm_wrapper_impl::simple_wrapper(wp, ukpi, sol, argc, argv);
  }
}

#endif //HBM_WRAPPER_HPP
