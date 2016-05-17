#include <cstdlib>
#include <iostream>

using namespace std;

// TODO
int main(int argc, char** argv) {
  if (argc != 1) {
    cout << "usage: " << argv[0] << " < file.sukp > file.ukp" << endl;
    // TODO: Fix line length
    cout << "This program reads an UKP instance on the 'sukp'"
         << " format (from the standard input), and writes the same instance on the"
         << " 'ukp' format to the standard output.";
    return EXIT_FAILURE;
  }

  sukp2ukp(cin, cout);
}

