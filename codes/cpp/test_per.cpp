#include "periodicity.hpp"
#include "test_common.hpp"

int main(int argc, char** argv) {
  return hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::y_star_wrapper);
}

