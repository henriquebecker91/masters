#include "mtu.hpp"
#include "test_common.hpp"

int main(int argc, char** argv) {
  return hbm::main_take_path<size_t, size_t, size_t>(&hbm::mtu1, argc, argv);
}

