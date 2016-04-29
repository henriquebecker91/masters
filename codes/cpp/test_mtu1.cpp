#include "mtu.hpp"
#include "test_common.hpp"

int main(int argc, char** argv) {
  std::cout << "hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::mtu1, argc, argv)" << std::endl;
  return hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::mtu1, argc, argv);
}

