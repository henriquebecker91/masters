#include "eduk.hpp"
#include "test_common.hpp"

int main(int argc, char** argv) {
  return hbm::benchmark_pyasukp(&hbm::eduk<hbm::weight, hbm::profit>);
}

