#include "ukp5.hpp"
#include "test_common.hpp"

// Solves an UKP instance by the UKP5 algorithm. Probably print lots
// of extra information in addition to the optimal solution.
// Takes the name of a file in the ".ukp" format. Other options
// should be consulted at ukp5.hpp.
int main(int argc, char** argv) {
    return hbm::main_take_path<size_t, size_t, size_t>(&hbm::ukp5, argc, argv);
}

