#include "cplex_ukp_model.hpp"
#include "test_common.hpp"

// Solves an UKP instance with the cplex solver. Probably print lots
// of extra information in addition to the optimal solution.
// Takes the name of a file in the ".ukp" format. Other options
// should be consulted at cplex_ukp_model.hpp.
int main(int argc, char** argv) {
    return hbm::main_take_path<size_t, size_t, size_t>(&hbm::cplex_solve_ukp, argc, argv);
}

