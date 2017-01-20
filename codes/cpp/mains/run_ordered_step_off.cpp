#include "stepoff.hpp"
#include "test_common.hpp"

// Execute the `ordered step-off' algorithm presented in "The theory and
// computation of knapsack functions" by Gilmore and Gomory, over an UKP
// instance. Takes the name of a file in the ".ukp" format. Other options
// should be consulted at stepoff.hpp.
int main(int argc, char** argv) {
  return hbm::main_take_path<size_t, size_t, size_t>(&hbm::ordered_step_off, argc, argv);
}

