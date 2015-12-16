#include "dominance.hpp"
#include "test_common.hpp"

using namespace std;
using namespace hbm;

int main(int argc, char** argv) {
    return hbm::main_take_path<size_t, size_t, size_t>(&hbm::sel_not_sm_dom_wrapper, argc, argv);
}

