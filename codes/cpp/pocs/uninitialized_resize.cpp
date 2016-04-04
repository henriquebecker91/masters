#include <cstdlib>
#include <iostream>
#include <vector>
#include <chrono>     

using namespace std;
using namespace std::chrono;

double inline difftime_between_now_and(const steady_clock::time_point &begin) {
  return duration_cast< duration<double> >(steady_clock::now() - begin)
         .count();
}

// Stolen from: http://stackoverflow.com/questions/21028299
// Allocator adaptor that interposes construct() calls to
// convert value initialization into default initialization.
template <typename T, typename A=std::allocator<T>>
class default_init_allocator : public A {
  typedef std::allocator_traits<A> a_t;
public:
  template <typename U> struct rebind {
    using other =
      default_init_allocator<
        U, typename a_t::template rebind_alloc<U>
      >;
  };

  using A::A;

  template <typename U>
  void construct(U* ptr)
    noexcept(std::is_nothrow_default_constructible<U>::value) {
    ::new(static_cast<void*>(ptr)) U;
  }
  template <typename U, typename...Args>
  void construct(U* ptr, Args&&... args) {
    a_t::construct(static_cast<A&>(*this),
                   ptr, std::forward<Args>(args)...);
  }
};

int main(void)
{
  vector<size_t> v; cout << "vector<size_t> v;" << endl;
  v.reserve(10);    cout << "v.reserve(10);" << endl;
  cout << "After setting with data()[i]: " << endl;
  for (size_t i = 0; i < 10; ++i) {
    cout << "v.data()[" << i << "]: " << (v.data()[i] = i) << endl;
  }
  v.assign(5, 42);  cout << "v.assign(5, 42);" << endl;
  v.resize(10);     cout << "v.resize(10);" << endl;
  for (size_t i = 0; i < 10; ++i) {
    cout << "v.data()[" << i << "]: " << v.data()[i] << endl;
  }
  
  vector<size_t, default_init_allocator<size_t>> v2;
  cout << "vector<size_t, default_init_allocator<size_t>> v2;" << endl;
  v2.reserve(10); cout << "v.reserve(10);" << endl;
  cout << "After setting with data()[i]: " << endl;
  for (size_t i = 0; i < 10; ++i) {
    cout << "v2.data()[" << i << "]: " << (v2.data()[i] = i) << endl;
  }
  v2.assign(5, 42);  cout << "v2.assign(5, 42);" << endl;
  v2.resize(10);     cout << "v2.resize(10);" << endl;
  for (size_t i = 0; i < 10; ++i) {
    cout << "v2.data()[" << i << "]: " << v2.data()[i] << endl;
  }

  cout << "---------- Benchmark ----------" << endl;

  steady_clock::time_point begin = steady_clock::now();
  v.resize(10000000);
  cout << "v.resize(10^7): " << difftime_between_now_and(begin) << endl;

  begin = steady_clock::now();
  v2.resize(10000000);
  cout << "v2.resize(10^7): " << difftime_between_now_and(begin) << endl;

  /* On an ASUS R552JK-CN159H (Intel Core i7-4700HQ, 6M Cache, 3.40 GHz):
   * g++ -std=c++11 -O0 uninitialized_resize.cpp && ./a.out | tail -n 3
   * ---------- Benchmark ----------
   * v.resize(10^7):  0.027077
   * v2.resize(10^7): 0.145314
   *
   * g++ -std=c++11 -O3 uninitialized_resize.cpp && ./a.out | tail -n 3
   * ---------- Benchmark ----------
   * v.resize(10^7):  0.0151531
   * v2.resize(10^7): 8.313e-06
   */
	return EXIT_SUCCESS;
}

