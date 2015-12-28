#ifndef HBM_TEST_COMMON_HPP
#define HBM_TEST_COMMON_HPP

#include <chrono>
#include <fstream>
#include <array>
#include <iostream>

#include "ukp_common.hpp"

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
    int run_ukp(ukp_solver_t<W, P, I> ukp_solver, int argc, char ** argv, const string& path, run_t<W, P, I> &run) {
      ifstream f(path);

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
        cout << "Couldn't open file: " << path << endl;
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
    int benchmark_pyasukp(ukp_solver_t<W, P, I> ukp_solver, int argc, char **argv) {
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

        run_t<W, P, I> run;
        // Call the function without any flags (argc = 0, argv = NULL)
        int status = hbm::hbm_test_common_impl::run_ukp(ukp_solver, argc-1, argv+1, path, run);

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

    #if defined(HBM_PROFILE) && defined(HBM_DUMP)
    template <typename W_OR_P>
    void dump(const string &path, const string &header, const vector<W_OR_P> &v) {
      ofstream f(path, ofstream::out|ofstream::trunc);
      if (f.is_open())
      {
        f << header << endl;
        for (size_t y = 0; y < v.size(); ++y) {
          f << y << "\t" << v[y] << endl;
        }
      } else {
        cerr << "Couldn't open file: " << path << endl;
      }
    }
    #endif

    template <typename W, typename P, typename I>
    int main_take_path(ukp_solver_t<W, P, I> ukp_solver, int argc, char** argv) {
      if (argc < 2) {
        cout << "usage: a.out data.ukp [extra specific method params]" << endl;
        return EXIT_FAILURE;
      }

      string spath(argv[1]);
      cout << spath << endl;

      run_t<W, P, I> run;

      // Call the function without giving the program name or instance filename
      int status = hbm::hbm_test_common_impl::run_ukp(ukp_solver, argc-2, argv+2, spath, run);

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

  /// Run a custom procedure over the ukp instance file at path.
  ///
  /// The time for reading the instance isn't measured, only the time taken
  /// by the ukp_solver is measured.
  ///
  /// @param ukp_solver A procedure that takes an instance_t, some
  ///   extra parameters coded as argc/argv and writes the result
  ///   on a solution_t.
  /// @param argc Number of extra parameters of ukp_solver.
  /// @param argv Extra parameters given without change to ukp_solver.
  /// @param path A string with the path to an instance on the UPK format.
  /// @param run Where the solution_t written by ukp_solver and the time
  ///   used by the the ukp_solver will be stored.
  ///
  /// @return EXIT_SUCCESS or EXIT_FAILURE. It's a failure if the file
  ///   at path can't be opened.
  ///
  /// @exception ukp_read_error If the instance format is wrong.
  /// @see read_ukp_instance For the ukp format.
  template <typename W, typename P, typename I>
  int run_ukp(ukp_solver_t<W, P, I> ukp_solver, int argc, char **argv, const std::string& path, run_t<W, P, I> &run) {
    return hbm_test_common_impl::run_ukp(ukp_solver, argc, argv, path, run);
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
  ///
  /// @param ukp_solver A procedure that takes an instance_t, some
  ///   extra parameters coded as argc/argv and writes the result
  ///   on a solution_t.
  ///
  /// @return EXIT_SUCCESS or EXIT_FAILURE. It's a failure if the
  ///   PYAsUKP benchmark files aren't found at the hardcoded path.
  /// @exception ukp_read_error If the instance format is wrong.
  template <typename W, typename P, typename I>
  int benchmark_pyasukp(ukp_solver_t<W, P, I> ukp_solver, int argc, char** argv) {
    return hbm_test_common_impl::benchmark_pyasukp(ukp_solver, argc, argv);
  }

  /// @brief Takes the first argument (argv[1]) as the instance file and
  ///   gives the remaining arguments to the custom procedure.
  ///
  /// The custom procedure receive argc as argc-2 and argv as argv+2
  ///   (the executable filename and the instance filename are skipped).
  /// 
  /// @param ukp_solver A procedure that takes an instance_t, some
  ///   extra parameters coded as argc/argv and writes the result
  ///   on a solution_t. Remember that the argc and argv given to
  ///   the ukp_solver are a subset of the ones given to main_take_path.
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
  int main_take_path(ukp_solver_t<W, P, I> ukp_solver, int argc, char** argv) {
    return hbm_test_common_impl::main_take_path(ukp_solver, argc, argv);
  }
}

#endif //HBM_TEST_COMMON_HPP
