#include "stepoff.hpp"
#include "test_common.hpp"

// Execute the `ordered step-off' algorithm presented in "The theory and
// computation of knapsack functions" (by Gilmore and Gomory) over a hardcoded
// set of instances. Used to check if it is still working after a change.
int main(int argc, char** argv) {
  std::cout << "hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::ordered_step_off, argc, argv)" << std::endl;
  return hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::ordered_step_off, argc, argv);
}

