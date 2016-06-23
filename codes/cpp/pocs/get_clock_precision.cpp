#include <cstdlib>
#include <iostream>
#include <chrono>
#include <ratio>

using namespace std;
using namespace std::chrono;

// NOTE: The compiled binary will only display the result for the
// machine/compiler that it was compiled. So, after compiled, the binary
// will show always the same value independently from the machine it is
// being run. To measure for your compiler/machine you need to compile
// this source again.
int main (void)
{
  cout  << "The steady_clock period is: "
        << steady_clock::period::num << "/" << steady_clock::period::den
        << " seconds."<< endl;
  if (ratio_greater<steady_clock::period, hours::period>::value) {
    cout << "It don't measure hours with precision." << endl;
  } else if (ratio_equal<hours::period, steady_clock::period>::value) {
    cout << "The precision is exactly the needed for measuring hours." << endl;
  } else if (ratio_greater<steady_clock::period, minutes::period>::value) {
    cout << "It don't measure minutes with precision." << endl;
  } else if (ratio_equal<minutes::period, steady_clock::period>::value) {
    cout << "The precision is exactly the needed for measuring minutes." << endl;
  } else if (ratio_greater<steady_clock::period, seconds::period>::value) {
    cout << "It don't measure seconds with precision." << endl;
  } else if (ratio_equal<seconds::period, steady_clock::period>::value) {
    cout << "The precision is exactly the needed for measuring seconds." << endl;
  } else if (ratio_greater<steady_clock::period, milliseconds::period>::value) {
    cout << "It don't measure milliseconds with precision." << endl;
  } else if (ratio_equal<milliseconds::period, steady_clock::period>::value) {
    cout << "The precision is exactly the needed for measuring milliseconds." << endl;
  } else if (ratio_greater<steady_clock::period, microseconds::period>::value) {
    cout << "It don't measure microseconds with precision." << endl;
  } else if (ratio_equal<microseconds::period, steady_clock::period>::value) {
    cout << "The precision is exactly the needed for measuring microseconds." << endl;
  } else if (ratio_greater<steady_clock::period, nanoseconds::period>::value) {
    cout << "It don't measure nanoseconds with precision." << endl;
  } else if (ratio_equal<nanoseconds::period, steady_clock::period>::value) {
    cout << "The precision is exactly the needed for measuring nanoseconds." << endl;
  } else {
    cout << "The precision is greater than the one needed to measure nanoseconds." << endl;
  }

	return EXIT_SUCCESS;
}

