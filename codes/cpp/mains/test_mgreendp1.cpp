#include "greendp.hpp"
#include "test_common.hpp"

// Execute the modern version of the first dynamic programming algorithm
// presented at "On Equivalent Knapsack Problems" (H. Greenberg) over a
// hardcoded set of instances. Used to check if it is still working after a
// change.
int main(int argc, char** argv) {
  std::cout << "hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::greendp1, argc, argv)" << std::endl;
  return hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::mgreendp1, argc, argv);
}

