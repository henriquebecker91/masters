#include "babayev.hpp"
#include "test_common.hpp"

// NOTE: work in progress, for now is returning the u3 bound of Martello and
// Toth.
int main(int argc, char** argv) {
  std::cout << "hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::babayev, argc, argv)" << std::endl;
  return hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::babayev, argc, argv);
}

