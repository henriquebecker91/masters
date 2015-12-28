#include "periodicity.hpp"
#include "test_common.hpp"

int main(int argc, char** argv) {
  int exit_code;

  std::cout << "hbm::<size_t, size_t, size_t>(&hbm::y_star_wrapper, argc, argv);" << std::endl;
  exit_code = hbm::main_take_path<size_t, size_t, size_t>(&hbm::y_star_wrapper, argc, argv);
  if (exit_code != EXIT_SUCCESS) return exit_code;

  std::cout << "hbm::<size_t, double, size_t>(&hbm:y_star_wrapper:, argc, argv);" << std::endl;
  return hbm::main_take_path<size_t, long double, size_t>(&hbm::y_star_wrapper, argc, argv);
}

