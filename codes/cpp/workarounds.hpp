#ifndef HBM_WORKAROUNDS_HPP
#define HBM_WORKAROUNDS_HPP

#include <sstream>

namespace hbm {
  /// @brief Convert the string to type T.
  ///
  /// This works as we were reading the string from the standard input with
  /// the >> operator. So this can be used with any type that defines an >> 
  /// operator.
  ///
  /// @tparam T A type that defines the >> operator.
  /// @param s The string representation of a value of type T.
  /// @param t The variable that will receive the value represented in
  ///   the string.
  template <class T>
  inline void from_string(const std::string& s, T& t)
  {
    std::stringstream ss(s);
    ss >> t;
  }
}

#endif // HBM_WORKAROUNDS_HPP

