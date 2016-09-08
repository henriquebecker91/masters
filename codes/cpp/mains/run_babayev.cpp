#include "babayev.hpp"
#include "test_common.hpp"

// NOTE: work in progress, for now is returning the u3 bound of Martello and
// Toth.
// NOTE: using long long ints because negative numbers are extensively used
// by the algorithm, also operations between profits and weights are common
// (so it's good to have the same type for both).
int main(int argc, char** argv) {
    return hbm::main_take_path<long long int, long long int, size_t>(&hbm::babayev, argc, argv);
}

