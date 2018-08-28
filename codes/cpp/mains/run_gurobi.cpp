#include "gurobi_ukp_model.hpp"
#include "test_common.hpp"

// Solves an UKP instance withe the gurobi solver. Probably print lots
// of extra information in addition to the optimal solution.
// Takes the name of a file in the ".ukp" format. Other options
// should be consulted at gurobi_ukp_model.hpp.
int main(int argc, char** argv) {
    return hbm::main_take_path<size_t, size_t, size_t>(&hbm::gurobi_solve_ukp, argc, argv);
}

