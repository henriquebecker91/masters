#ifndef HBM_TYPE_NAME_HPP
#define HBM_TYPE_NAME_HPP

#include <string>
#include <type_traits>
#include <iostream>
#include <cstdint>

namespace hbm {
  // This code was shamelessly stolen from the post of bames53 (on
  // http://arstechnica.com/civis/viewtopic.php?f=20&t=1192725)
  // by Henrique Becker (on 20/07/2016).
  template<typename T> struct type_name {
      static std::string get() {
          if (std::is_enum<T>::value) return "enum";
          if (std::is_class<T>::value) return "class";
          if (std::is_union<T>::value) return "union";
          return "unknown";
      }
  };

  #define HBM_TYPE_NAME(T) template<> struct type_name<T> { static std::string get() { return #T; } }
  HBM_TYPE_NAME(char);
  HBM_TYPE_NAME(unsigned char);
  HBM_TYPE_NAME(signed char);
  HBM_TYPE_NAME(char16_t);
  HBM_TYPE_NAME(char32_t);
  HBM_TYPE_NAME(bool);
  HBM_TYPE_NAME(unsigned int);
  HBM_TYPE_NAME(int);
  HBM_TYPE_NAME(unsigned short);
  HBM_TYPE_NAME(unsigned long);
  HBM_TYPE_NAME(unsigned long long);
  HBM_TYPE_NAME(long);
  HBM_TYPE_NAME(long long);
  HBM_TYPE_NAME(short);
  HBM_TYPE_NAME(wchar_t);
  HBM_TYPE_NAME(float);
  HBM_TYPE_NAME(double);
  HBM_TYPE_NAME(long double);
  HBM_TYPE_NAME(void);
  #undef HBM_TYPE_NAME

  template<typename T> struct type_name<T*> {
      static std::string get() { return "pointer to " + type_name<T>::get(); }
  };

  template<typename T> struct type_name<T&> {
      static std::string get() { return "reference to " + type_name<T>::get(); }
  };

  template<typename T> struct type_name<T&&> {
      static std::string get() { return "rvalue reference to " + type_name<T>::get(); }
  };

  template<typename T> struct type_name<T const> {
      static std::string get() { return "const " + type_name<T>::get(); }
  };

  template<typename T> struct type_name<T volatile> {
      static std::string get() { return "volatile " + type_name<T>::get(); }
  };

  template<typename T> struct type_name<T const volatile> {
      static std::string get() { return "const volatile " + type_name<T>::get(); }
  };

  template<typename T, typename C> struct type_name<T C::*> {
      static std::string get() { return "pointer to " + type_name<C>::get() + " member " + type_name<T>::get(); }
  };

  template<typename T, std::intmax_t N> struct type_name<T[N]> {
      static std::string get() { return "array of " + std::to_string(N) + " " + type_name<T>::get(); }
  };

  // take care of annoying rule that array of cv-qualified type is the
  // same as cv-qualified array
  template<typename T, std::intmax_t N> struct type_name<T const [N]> {
      static std::string get() { return "array of " + std::to_string(N) + " const " + type_name<T>::get(); }
  };

  template<typename T, std::intmax_t N> struct type_name<T volatile [N]> {
      static std::string get() { return "array of " + std::to_string(N) + " volatile " + type_name<T>::get(); }
  };

  template<typename T, std::intmax_t N> struct type_name<T const volatile [N]> {
      static std::string get() { return "array of " + std::to_string(N) + " const volatile " + type_name<T>::get(); }
  };

  template<typename... Args> struct type_name_list { static std::string get() { return ""; } };
  template<typename Arg> struct type_name_list<Arg> { static std::string get() { return type_name<Arg>::get(); } };
  template<typename Arg, typename... Args> struct type_name_list<Arg,Args...> { static std::string get() { return type_name<Arg>::get() + ", " + type_name_list<Args...>::get(); } };

  template<typename RetT, typename... Args> struct type_name<RetT(Args...)> {
      static std::string get() { return "function taking (" + type_name_list<Args...>::get() + ") and returning " + type_name<RetT>::get(); }
  };
}

#endif // HBM_TYPE_NAME_HPP

