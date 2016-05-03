#include "mtu.hpp"
#include "test_common.hpp"

// Solves an UKP instance by the mtu1 algorithm (created by Martello and Toth
// and described in the book "Knapsack Problems" p. 96). Takes the name of a
// file in the ".ukp" format. Other options should be consulted at mtu.hpp.
int main(int argc, char** argv) {
  return hbm::main_take_path<size_t, size_t, size_t>(&hbm::mtu1, argc, argv);
}

