#include <cstdlib>
#include <iostream>

#ifndef HBM_GIT_HEAD_AT_COMPILATION
  #define HBM_GIT_HEAD_AT_COMPILATION "compiler call forgot to set \
HBM_GIT_HEAD_AT_COMPILATION"
#endif

using namespace std;

int main (void)
{
  cout << HBM_GIT_HEAD_AT_COMPILATION << endl;

	return EXIT_SUCCESS;
}

