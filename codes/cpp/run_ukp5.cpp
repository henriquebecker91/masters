#include "ukp5.hpp"
#include "test_common.hpp"

int main(int argc, char** argv) {
    return hbm::main_take_path<size_t, size_t, size_t>(&hbm::ukp5, argc, argv);
}

