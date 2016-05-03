#include "greendp.hpp"
#include "test_common.hpp"

// Execute the modernized version of the dynamic programming algorithm
// presented at "A Better Step-off Algorithm for the Knapsack Problem" over an
// UKP instance. Takes the name of a file in the ".ukp" format. Other options
// should be consulted at greendp.hpp.
int main(int argc, char** argv) {
  return hbm::main_take_path<size_t, size_t, size_t>(&hbm::mgreendp, argc, argv);
}


