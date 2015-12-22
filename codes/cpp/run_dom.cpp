#include "dominance.hpp"
#include "test_common.hpp"

using namespace std;
using namespace hbm;

int main(int argc, char** argv) {
  int exit_code;

  std::cout << "hbm::main_take_path<size_t, size_t, size_t>"
               "(&hbm::smdom_sort_wrapper, argc, argv)"
               << std::endl;
  exit_code = hbm::main_take_path<size_t, size_t, size_t>
    (&hbm::smdom_sort_wrapper, argc, argv);

  if (exit_code != EXIT_SUCCESS) return exit_code;

  std::cout << "hbm::main_take_path<size_t, size_t, size_t>"
               "(&hbm::smdom_no_sort_wrapper, argc, argv)"
               << std::endl;
  exit_code = hbm::main_take_path<size_t, size_t, size_t>
    (&hbm::smdom_no_sort_wrapper, argc, argv);

  return exit_code;
}

