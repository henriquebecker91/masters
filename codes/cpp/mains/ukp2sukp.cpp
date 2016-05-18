#include <cstdlib>
#include <iostream>
#include "ukp_common.hpp"

using namespace std;

// Convert a instance from ukp format to sukp format, by standard input and
// output.
int main(int argc, char** argv) {
  if (argc != 1) {
    cout << "usage: " << argv[0] << " < file.ukp > file.sukp" << endl;
    cout << "This program reads an UKP instance on the 'ukp' format from the"
         << " standard input, and writes the same instance on the 'sukp'"
         << " format to the standard output." << endl;
    return EXIT_FAILURE;
  }

  hbm::ukp2sukp(cin, cout);
}

