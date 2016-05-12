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

int main(void)
{
  const size_t size = 10000000;
  vector<int> v(size, 0);

  cout << "---------- Benchmark ----------" << endl;

  steady_clock::time_point begin = steady_clock::now();
  for (int& i : v) { ++i; }
  cout << "for_each c++11:     " << difftime_between_now_and(begin) << endl;

  begin = steady_clock::now();
  for (auto it = v.data(); it != v.data() + size; ++it) { ++(*it); }
  cout << "ptr nocache:        " << difftime_between_now_and(begin) << endl;

  begin = steady_clock::now();
  const auto v_end_ptr = v.data() + size;
  for (auto it = v.data(); it != v_end_ptr; ++it) { ++(*it); }
  cout << "ptr cache:          " << difftime_between_now_and(begin) << endl;

  begin = steady_clock::now();
  for (auto it = v.begin(); it != v.end(); ++it) { ++(*it); }
  cout << "iterator nocache:   " << difftime_between_now_and(begin) << endl;

  begin = steady_clock::now();
  const auto v_end = v.end();
  for (auto it = v.begin(); it != v_end; ++it) { ++(*it); }
  cout << "iterator cache:     " << difftime_between_now_and(begin) << endl;

  begin = steady_clock::now();
  for (size_t i = 0; i < v.size(); ++i) { ++(v[i]); }
  cout << "operator[] nocache: " << difftime_between_now_and(begin) << endl;

  begin = steady_clock::now();
  for (size_t i = 0; i < size; ++i) { ++(v[i]); }
  cout << "operator[] cache:   " << difftime_between_now_and(begin) << endl;

  begin = steady_clock::now();
  for (size_t i = 0; i < v.size(); ++i) { ++(v.at(i)); }
  cout << "at nocache:         " << difftime_between_now_and(begin) << endl;

  begin = steady_clock::now();
  for (size_t i = 0; i < size; ++i) { ++(v.at(i)); }
  cout << "at cache:           " << difftime_between_now_and(begin) << endl;

  /* On an ASUS R552JK-CN159H (Intel Core i7-4700HQ, 6M Cache, 3.40 GHz),
   * with gcc version 6.1.1 20160501 (x86_64-pc-linux-gnu):
   *
   * g++ -O0 vector_index_vs_iterator.cpp && ./a.out 
   * ---------- Benchmark ----------
   * for_each c++11:     0.0855757
   * ptr nocache:        0.0493906
   * ptr cache:          0.019197
   * iterator nocache:   0.126014
   * iterator cache:     0.082235
   * operator[] nocache: 0.046754
   * operator[] cache:   0.0316647
   * at nocache:         0.130557
   * at cache:           0.105731
   * g++ -O3 vector_index_vs_iterator.cpp && ./a.out 
   * ---------- Benchmark ----------
   * for_each c++11:     0.00472208
   * ptr nocache:        0.00462304
   * ptr cache:          0.00498309
   * iterator nocache:   0.00514497
   * iterator cache:     0.00495087
   * operator[] nocache: 0.00447077
   * operator[] cache:   0.00446857
   * at nocache:         0.00472267
   * at cache:           0.00589979
   *
   * The author of this code has no ideia on how to explain those results.
   */
	return EXIT_SUCCESS;
}

