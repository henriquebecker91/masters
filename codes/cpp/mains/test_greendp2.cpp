#include "greendp.hpp"
#include "test_common.hpp"

// Execute the second dynamic programming algorithm presented at "On Equivalent
// Knapsack Problems" (H. Greenberg) over a hardcoded set of instances. Used
// to check if it is still working after a change.
// NOTE: this algorithm fails for corepb.ukp this is expected. This algorithm
// don't work for every instance.
int main(int argc, char** argv) {
  std::cout << "hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::greendp2, argc, argv)" << std::endl;
  return hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::greendp2, argc, argv);
}

