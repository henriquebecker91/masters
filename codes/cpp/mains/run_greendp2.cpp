#include "greendp.hpp"
#include "test_common.hpp"

// Solves an UKP instance by the greendp2 algorithm (the second algorithm
// presented at "On Equivalent Knapsack Problems", H. Greenberg). Takes the
// name of a file in the ".ukp" format. Other options should be consulted at
// greendp.hpp.
int main(int argc, char** argv) {
  return hbm::main_take_path<size_t, size_t, size_t>(&hbm::greendp2, argc, argv);
}

