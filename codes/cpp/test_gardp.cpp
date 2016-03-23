#include "gardp.hpp"
#include "test_common.hpp"

// Execute the dynamic programming algorithm presented at p. 221, Integer
// Programming, Robert S. Garfinkel, over a hardcoded set of instances. Used to
// check if it is still working after a change. All the instances have only
// integers, but we execute the gardp version where the profit values are
// stored as doubles to test it.
int main(int argc, char** argv) {
  int exit_code;

  std::cout << "hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::gardp, argc, argv)" << std::endl;
  exit_code = hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::gardp, argc, argv);
  if (exit_code != EXIT_SUCCESS) return exit_code;

  std::cout << "hbm::benchmark_pyasukp<size_t, double, size_t>(&hbm::gardp)" << std::endl;
  return hbm::benchmark_pyasukp<size_t, double, size_t>(&hbm::gardp, argc, argv);
}

