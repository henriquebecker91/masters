#include <iostream>
#include <fstream>
#include <array>
#include <iomanip>

#include "test_common.hpp"

#ifdef HBM_PROFILE
  #include <boost/filesystem.hpp>
  using namespace boost::filesystem;
#endif
#ifndef HBM_PROFILE_PRECISION
  #define HBM_PROFILE_PRECISION 5
#endif

using namespace std;
using namespace std::chrono;

int run_ukp(void(*ukp_solver)(ukp_instance_t &, ukp_solution_t &, bool), const string& path, run_t &run) {
  ifstream f(path);

  if (f.is_open())
  {
    ukp_instance_t ukpi;
    ukp_solution_t &ukps = run.result;

    read_ukp_instance(f, ukpi);

    steady_clock::time_point t1 = steady_clock::now();
    (*ukp_solver)(ukpi, ukps, false);
    steady_clock::time_point t2 = steady_clock::now();
    run.time = duration_cast<duration<double>>(t2 - t1);
  } else {
    cout << "Couldn't open file: " << path << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int benchmark_pyasukp(void(*ukp_solver)(ukp_instance_t &, ukp_solution_t &, bool)) {
  array<instance_data_t, /*sizeof(instance_data_t)**/8> instances_data = {{
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
    string path = "../../data/ukp/" + instances_data[i].name + ".ukp";
    cout << path << endl;

    run_t run;
    int status = run_ukp(ukp_solver, path, run);

    if (status == EXIT_SUCCESS) {
      size_t expected = instances_data[i].expected_opt;
      size_t obtained = run.result.opt;
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

#if !defined(HBM_PROFILE) && defined(HBM_DUMP)
  #error The HBM_DUMP flag can only be used with the HBM_PROFILE flag
#endif

#if defined(HBM_PROFILE) && defined(HBM_DUMP)
void dump(const string &path, const string &header, const vector<size_t> &v) {
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

int main_take_path(void(*ukp_solver)(ukp_instance_t &, ukp_solution_t &, bool), int argc, char** argv) {
  if (argc != 2) {
    cout << "usage: a.out data.sukp" << endl;
    return EXIT_FAILURE;
  }

  string spath(argv[1]);
  cout << spath << endl;

  run_t run;

  int status = run_ukp(ukp_solver, spath, run);

  if (status == EXIT_SUCCESS) {
    const auto &res = run.result;

    #if defined(HBM_PROFILE) && defined(HBM_DUMP)
    path my_path(spath);
    path filename = my_path.filename();
    const string  gd_path = "./g_dump_" + filename.native() + ".dat",
                  dd_path = "./d_dump_" + filename.native() + ".dat",
                  nsd_path = "./nsd_dump_" + filename.native() + ".dat",
                  dqt_path = "./dqt_dump_" + filename.native() + ".dat";
    dump(gd_path, "y\tgy", res.g);
    dump(dd_path, "y\tdy", res.d);
    //dump(nsd_path, "y\tdy", res.non_skipped_d);
    //dump(dqt_path, "i\tqt_in_d", res.qt_i_in_dy);
    #endif

    cout << "opt:    " << res.opt << endl;
    cout << "y_opt:  " << res.y_opt << endl;
    #if defined(HBM_CHECK_PERIODICITY) || defined(HBM_CHECK_PERIODICITY_FAST)
    cout << "last_y_value_outer_loop: " << res.last_y_value_outer_loop << endl;
    #endif
    for (auto it = res.used_items.cbegin(); it != res.used_items.cend(); ++it) {
      it->print();
    }
    #ifdef HBM_PROFILE
    cout << "c: " << res.c << endl;
    cout << "n: " << res.n << endl;
    cout << "w_min: " << res.w_min << endl;
    cout << "w_max: " << res.w_max << endl;
    cout << "last_dy_non_zero_non_n: " << res.last_dy_non_zero_non_n << endl;
    cout << "qt_non_skipped_ys: " << res.qt_non_skipped_ys << endl;
    cout << "qt_gy_zeros: " << res.qt_gy_zeros << endl;
    cout << "qt_inner_loop_executions: " << res.qt_inner_loop_executions << endl;
    cout << "qt_inner_loop_executions/qt_non_skipped_ys: " << ((long double) res.qt_inner_loop_executions)/((long double) res.qt_non_skipped_ys) << endl;
    cout << "qt_inner_loop_executions/c: " << ((long double) res.qt_inner_loop_executions)/((long double) res.c) << endl;
    cout << "(qt_inner_loop_executions/qt_non_skipped_ys)/n: " << ((long double) res.qt_inner_loop_executions)/((long double) res.qt_non_skipped_ys)/((long double)res.n)  << endl;
    cout << "(qt_inner_loop_executions/c)/n: " << ((long double) res.qt_inner_loop_executions)/((long double) res.c)/((long double)res.n) << endl;

    const double &stime = res.sort_time, &vtime = res.vector_alloc_time,
      &lctime = res.linear_comp_time, &p1time = res.phase1_time,
      &p2time = res.phase2_time, &ttime = res.total_time;

    streamsize old_precision = cout.precision(HBM_PROFILE_PRECISION);
    int percent_size = 3+HBM_PROFILE_PRECISION;
    ios_base::fmtflags old_flags = cout.setf(std::ios::fixed, std:: ios::floatfield);
    char old_fill = cout.fill(' ');

    cout << "Sort time: " << stime << "s (";
    cout << setw(percent_size) << (stime/ttime)*100;
    cout << "%)" << endl;
    cout << "Vect time: " << vtime << "s (";
    cout << setw(percent_size) << (vtime/ttime)*100;
    cout << "%)" << endl;
    cout << "O(n) time: " << lctime << "s (";
    cout << setw(percent_size) << (lctime/ttime)*100;
    cout << "%)" << endl;
    cout << "pha1 time: " << p1time << "s (";
    cout << setw(percent_size) << (p1time/ttime)*100;
    cout << "%)" << endl;
    cout << "pha2 time: " << p2time << "s (";
    cout << setw(percent_size) << (p2time/ttime)*100;
    cout << "%)" << endl;
    cout << "Sum times: " << ttime << "s" << endl;

    cout.fill(old_fill);
    cout.setf(old_flags);
    cout.precision(old_precision);
    #endif
    cout << "Seconds: " << run.time.count() << endl;
    cout << endl;
    return EXIT_SUCCESS;
  } else {
    cout << "There was some problem with this instance" << endl;
    return EXIT_FAILURE;
  }
}

