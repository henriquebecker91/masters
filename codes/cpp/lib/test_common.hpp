#ifndef HBM_TEST_COMMON_HPP
#define HBM_TEST_COMMON_HPP

#include "ukp_common.hpp"

#include <chrono>
#include <fstream>  // ifstream for reading UKP files
#include <array>
#include <iostream>

namespace hbm {
  /// @brief A type that contains a solution, and how much time was needed
  ///   to obtain it.
  template <typename W, typename P, typename I>
  struct run_t {
    /// An UKP solution.
    solution_t<W, P, I> result;
    /// The time needed to obtain the solution.
    std::chrono::duration<double> time; 
  };

  /// Inner test_common implementations. Do not depend.
  namespace hbm_test_common_impl {
    using namespace std;
    using namespace std::chrono;

    template <typename W, typename P, typename I>
    int run_ukp(ukp_solver_t<W, P, I> ukp_solver, int argc, argv_t argv, run_t<W, P, I> &run) {
      ifstream f(argv[1]);

      if (f.is_open())
      {
        instance_t<W, P> ukpi;
        solution_t<W, P, I> &ukps = run.result;

        hbm::hbm_ukp_common_impl::read_ukp_instance(f, ukpi);

        steady_clock::time_point t1 = steady_clock::now();
        (*ukp_solver)(ukpi, ukps, argc, argv);
        steady_clock::time_point t2 = steady_clock::now();
        run.time = duration_cast<duration<double>>(t2 - t1);
      } else {
        cout << "Couldn't open file: " << argv[1] << endl;
        return EXIT_FAILURE;
      }

      return EXIT_SUCCESS;
    }

    /// @brief Internal datatype used to store an instance name and its
    ///   optimal solution value.
    template <typename P>
    struct instance_data_t {
      std::string name;
      P expected_opt;
    };

    template <typename W, typename P, typename I>
    int benchmark_pyasukp(ukp_solver_t<W, P, I> ukp_solver, int argc, argv_t argv) {
      array<instance_data_t<P>, 8> instances_data = {{
        { "exnsd16", 1029680 },
        { "exnsd18", 1112131 },
        { "exnsd20",  1026086 },
        { "exnsd26", 1027564 },
        { "exnsdbis10", 1028035 },
        { "exnsdbis18",  1037156 },
        { "exnsds12", 3793952 },
        { "corepb", 10077782 }
      }};

      bool everything_ok = true;
      for (size_t i = 0; i < instances_data.size(); ++i) {
        const string path = "../../data/ukp/" + instances_data[i].name + ".ukp";
        cout << path << endl;

        // Insert the path as the second argument. The second argument is
        // expected to be the instance name by all the ukp_solver_t methods.
        const char ** hijacked_argv = new const char *[argc+1];
        hijacked_argv[0] = argv[0];
        hijacked_argv[1] = path.c_str();
        copy_n(&argv[1], argc - 1, &hijacked_argv[2]);
        int hijacked_argc = argc + 1;

        run_t<W, P, I> run;
        int status = hbm_test_common_impl::run_ukp(ukp_solver, hijacked_argc, hijacked_argv, run);

        delete[] hijacked_argv;
        if (status == EXIT_SUCCESS) {
          P expected = instances_data[i].expected_opt;
          P obtained = run.result.opt;
          cout << "Expected result: " << expected << endl;
          cout << "Obtained result: " << obtained << endl;
          everything_ok = everything_ok && expected == obtained;

          cout << "Seconds: " << run.time.count() << endl;

          cout << endl;
        } else {
          return EXIT_FAILURE;
        }
      }

      if (everything_ok) {
        cout << "Everything is ok." << endl;
      } else {
        cout << "The expected and obtained optimal values of some instances differ, check." << endl;
      }

      return EXIT_SUCCESS;
    }

    template <typename W, typename P, typename I>
    int main_take_path(ukp_solver_t<W, P, I> ukp_solver, int argc, argv_t argv) {
      run_t<W, P, I> run;

      if (argc < 2) {
        cout << "main_take_path: called with argc < 2, the expected usage is: "
             << "./<binary_name> path/file.ukp <other method specific params>."
             << endl << "main_take_path: will call method anyway to give it a "
             << "chance to show its specific usage."
        << endl;
        instance_t<W, P> ukpi;
        ukp_solver(ukpi, run.result, argc, argv);
        return EXIT_FAILURE;
      }

      const string spath(argv[1]);
      cout << "instance_filepath: " << spath << endl;

      // Call the function without giving the program name or instance filename
      int status = hbm::hbm_test_common_impl::run_ukp(ukp_solver, argc, argv, run);

      if (status == EXIT_FAILURE) {
        cout << "There was some problem with this instance. " << endl;
        return EXIT_FAILURE;
      }

      const auto &res = run.result;

      if (!res.show_only_extra_info) {
        // The printing of all the fields except extra_info should
        // be inside this condition
        cout << "opt:    " << res.opt << endl;
        cout << "y_opt:  " << res.y_opt << endl;
        for (auto it = res.used_items.cbegin(); it != res.used_items.cend(); ++it) {
          it->print();
        }
      }

      // The extra info and the timing are always shown
      string extra_info = res.extra_info->gen_info();
      if (!extra_info.empty()) cout << extra_info << endl;
      cout << "Seconds: " << run.time.count() << endl;
      cout << endl;

      return EXIT_SUCCESS;
    }
  }

  /// Run a custom procedure over the ukp instance file at argv[1].
  ///
  /// The time for reading the instance isn't measured, only the time taken by
  /// the ukp_solver is measured.
  ///
  /// @param ukp_solver A procedure that takes an instance_t, some extra
  ///   parameters coded as argc/argv and writes the result on a solution_t.
  /// @param argc Number of arguments stored on argv.
  /// @param argv Arguments passed by command-line. It's expected that the
  ///   argv[0] will be a program name, and the argv[1] will be the instance
  ///   filename. The instance data will be read from argv[1], and the argc
  ///   and argv will be passed to ukp_solver without change.
  /// @param run Where the solution_t written by ukp_solver and the time
  ///   used by the the ukp_solver will be stored.
  ///
  /// @return EXIT_SUCCESS or EXIT_FAILURE. It's a failure if the file
  ///   at argv[1] can't be opened.
  ///
  /// @exception ukp_read_error If the instance format is wrong.
  /// @see read_ukp_instance For the ukp format.
  template <typename W, typename P, typename I>
  int run_ukp(ukp_solver_t<W, P, I> ukp_solver, int argc, hbm::argv_t argv, run_t<W, P, I> &run) {
    return hbm_test_common_impl::run_ukp(ukp_solver, argc, argv, run);
  }

  /// Procedure used to test if a solver procedure is working.
  ///
  /// This procedure is very limited. It simply try to read the
  /// eight benchmark instances of PYAsUKP from ../../data/ukp/
  /// and execute ukp_solver over each one. Everything is
  /// hardcoded. This procedure aim is only to do a quick
  /// empirical test to check if the ukp_solver solver function
  /// remains working after a change. Used in the development
  /// of the solver functions.
  /// 28/12/2015: Now updated to receive and repass command-line
  /// arguments.
  ///
  /// @param ukp_solver A procedure that takes an instance_t, some
  ///   extra parameters coded as argc/argv and writes the result
  ///   on a solution_t.
  /// @param argc The number of arguments stored on argv. The argc passed
  ///   to ukp_solver is argc+1. See reasononig below (in argv parameter
  ///   explanation).
  /// @param argv The command-line arguments. It's expected that
  ///   argv[0] is the program name. The argv passed to the
  ///   ukp_solver parameter isn't this argv. It's a copy of this argv
  ///   with the instance filename inserted as argv[1] (old
  ///   parameter argv[1] will become argv[2] and so on).
  ///
  /// @return EXIT_SUCCESS or EXIT_FAILURE. It's a failure if the
  ///   PYAsUKP benchmark files aren't found at the hardcoded path.
  /// @exception ukp_read_error If the instance format is wrong.
  template <typename W, typename P, typename I>
  int benchmark_pyasukp(ukp_solver_t<W, P, I> ukp_solver, int argc, hbm::argv_t argv) {
    return hbm_test_common_impl::benchmark_pyasukp(ukp_solver, argc, argv);
  }

  /// @brief Takes the first argument (argv[1]) as the instance file and
  ///   gives the remaining arguments to the custom procedure.
  ///
  /// The position argv[0] is expected to have the executable name.
  /// The first argument is argv[1] (it's expected to be the instance
  /// filename). The custom procedure receive argc and argv as they are
  /// passed (the executable filename and the instance filename are NOT
  /// skipped). The instance_t passed to the ukp_solver is read from
  /// argv[1]. It's NOT expected from the ukp_solver to try reading the
  /// file itself (ukp_solver should probably ignore argv[1]).
  ///
  /// @note The main difference between this function and run_ukp
  ///   is that this function will output extra info to the stdout (cout).
  ///   The run_ukp only outputs a message if an error occur.
  ///   This function calls run_ukp internally.
  ///
  /// @param ukp_solver A procedure that takes an instance_t, some
  ///   parameters coded as argc/argv and writes the result on a
  ///   solution_t.
  /// @param argc The command-line arguments number. This is expected
  ///   to be greater than two (at least one argument has to exist,
  ///   the instance filename).
  /// @param argv The command-line arguments. The first argument
  ///   (argv[1] is expected to be a filename of an ukp format instance).
  ///
  /// @return EXIT_SUCCESS or EXIT_FAILURE. It's a failure if the instance
  ///   filename expected at argv[1] (first command line parameter) can't
  ///   be opened, or if argc < 2 (there's no argv[1]).
  ///
  /// @exception ukp_read_error If the instance format is wrong.
  /// @see read_ukp_instance For the ukp format.
  template <typename W, typename P, typename I>
  int main_take_path(ukp_solver_t<W, P, I> ukp_solver, int argc, hbm::argv_t argv) {
    return hbm_test_common_impl::main_take_path(ukp_solver, argc, argv);
  }
}

#endif //HBM_TEST_COMMON_HPP
