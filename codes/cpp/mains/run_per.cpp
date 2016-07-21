#include "periodicity.hpp"
#include "test_common.hpp"

// Computes the y* periodicity bound described by Garfinkel and Nemhauser
// at "Integer Programming", p. 223. Don't solve the UKP instance, only
// shows its periodicity bound (a knapsack capacity value). The integer
// and fractional versions are differently implemented, so we test both.
int main(int argc, char** argv) {
  int exit_code;

  std::cout << "hbm::main_take_path<size_t, size_t, size_t>(&hbm::y_star_wrapper, argc, argv);" << std::endl;
  exit_code = hbm::main_take_path<size_t, size_t, size_t>(&hbm::y_star_wrapper, argc, argv);
  if (exit_code != EXIT_SUCCESS) return exit_code;

  std::cout << "hbm::main_take_path<size_t, double, size_t>(&hbm:y_star_wrapper:, argc, argv);" << std::endl;
  return hbm::main_take_path<size_t, long double, size_t>(&hbm::y_star_wrapper, argc, argv);
}

