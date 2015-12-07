#include "test_common.hpp"
#include "eduk.hpp"

int main(int argc, char** argv) {
    return hbm::main_take_path<size_t, size_t, size_t>(&hbm::eduk, argc, argv);
}

