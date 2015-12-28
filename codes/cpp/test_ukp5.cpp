#include "ukp5.hpp"
#include "test_common.hpp"

int main(int argc, char** argv) {
  int exit_code;

  std::cout << "hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::ukp5, argc, argv)" << std::endl;
  exit_code = hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::ukp5, argc, argv);
  if (exit_code != EXIT_SUCCESS) return exit_code;

  std::cout << "hbm::benchmark_pyasukp<size_t, double, size_t>(&hbm::ukp5)" << std::endl;
  return hbm::benchmark_pyasukp<size_t, double, size_t>(&hbm::ukp5, argc, argv);
}

