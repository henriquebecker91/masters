#include "test_common.hpp"
#include "eduk.hpp"

// Solves an UKP instance by using an algorithm adapted from the paper: "Sparse
// knapsack algo-tech-cuit and its synthesis". The algorithm described at the
// paper was functional (ML). This is the result from a try to make it
// imperative (C++). The result is terribly slow, and consumed too much memory
// on big instances to be viable. Takes an instance filename in the ".ukp"
// format as first argument. Other options see eduk.hpp.
int main(int argc, hbm::argv_t argv) {
    return hbm::main_take_path<size_t, size_t, size_t>(&hbm::eduk, argc, argv);
}

