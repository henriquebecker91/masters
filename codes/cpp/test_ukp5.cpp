#include "ukp5.hpp"
#include "test_common.hpp"

int main(int argc, char** argv) {
  return hbm::benchmark_pyasukp<size_t, size_t, size_t>(&hbm::ukp5);
}

