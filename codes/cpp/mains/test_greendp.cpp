#include "greendp.hpp"
#include "test_common.hpp"

// Execute the dynamic programming algorithm presented at 
// "A Better Step-off Algorithm for the Knapsack Problem" over a hardcoded set
// of instances. Used to check if it is still working after a change.
int main(int argc, char** argv) {
  std::cout << "hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::greendp, argc, argv)" << std::endl;
  return hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::greendp, argc, argv);
}

