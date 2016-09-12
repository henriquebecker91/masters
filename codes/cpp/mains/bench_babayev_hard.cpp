#include "chung.hpp"
#include "ukp5.hpp"
#include "mtu.hpp"

#include <cstdlib>

using namespace std;
using namespace hbm;

// This is based on the five hard instances presented at "A New Knapsack
// Solution Approach by Integer Equivalent Aggregation and Consistency
// Determination". The five instances are created following the Chung
// distribution (from "A Hard Knapsack Problem").
int main(void) {
  // int is used because one instance has z == -4 and the values
  // magnitude is very low (int can handle them without problem).
  instance_t<int, int> h[5];

  gen_chung_hard_inst(  7, 100, 100, 9999, h[0]);
  gen_chung_hard_inst( 20, 100,   3, 9999, h[1]);
  gen_chung_hard_inst(100, 100, 100, 6060, h[2]);
  gen_chung_hard_inst(100,  60, 100, 3599, h[3]);
  gen_chung_hard_inst( 20,  82,  -4, 9700, h[4]);

  vector<ukp_solver_t<int, int, size_t>> algs;
  vector<string> alg_names;
  algs.push_back(&ukp5<int, int, size_t>);
  alg_names.push_back("ukp5");
  algs.push_back(&mtu1<int, int, size_t>);
  alg_names.push_back("mtu1");
  algs.push_back(&mtu2<int, int, size_t>);
  alg_names.push_back("mtu2");

  const char * fake_args[] = {"./fake_binary", "fake_filename.ukp"};
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < algs.size(); ++j) {
      solution_t<int, int, size_t> sol;
      (*algs[j])(h[i], sol, 2, fake_args);
      cout << alg_names[j] << endl;
      sol.print();
    }
  }

  return EXIT_SUCCESS;
}

