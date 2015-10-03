#include "periodicity.hpp"
#include "test_common.hpp"

int main(int argc, char** argv) {
  return benchmark_pyasukp(&y_star_wrapper);
}

