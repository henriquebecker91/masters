#include <cstdlib>
#include <iostream>
#include "ukp_common.hpp"

using namespace std;

// Convert a instance from sukp format to ukp format, by standard input and
// output.
int main(int argc, char** argv) {
  if (argc != 1) {
    cout << "usage: " << argv[0] << " < file.sukp > file.ukp" << endl;
    cout << "This program reads an UKP instance on the 'sukp' format from the"
         << " standard input, and writes the same instance on the 'ukp'"
         << " format to the standard output." << endl;
    return EXIT_FAILURE;
  }

  hbm::sukp2ukp(cin, cout);
}

