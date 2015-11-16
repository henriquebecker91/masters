#ifndef HBM_WORKAROUNDS_HPP
#define HBM_WORKAROUNDS_HPP

#include <sstream>

namespace hbm {
  /* It's utter ridiculous, but seems that is already 2015 and this function
   * isn't some sort of standard library utility. */
  template <class T>
  inline void from_string(const std::string& s, T& t)
  {
    std::stringstream ss(s);
    ss >> t;
  }
}

#endif // HBM_WORKAROUNDS_HPP

