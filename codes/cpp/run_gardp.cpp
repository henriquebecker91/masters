#include "gardp.hpp"
#include "test_common.hpp"

// Solves an UKP instance by the dynamic programming algorithm presented at p.
// 221, Integer Programming, Robert S. Garfinkel. Takes the name of a file in
// the ".ukp" format. Other options see gardp.hpp.
int main(int argc, char** argv) {
    return hbm::main_take_path<size_t, size_t, size_t>(&hbm::gardp, argc, argv);
}

