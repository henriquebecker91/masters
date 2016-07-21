#include <string>
#include <type_traits>
#include <iostream>
#include <cstdint>

// This code was shamelessly stolen from the post of bames53 (on
// http://arstechnica.com/civis/viewtopic.php?f=20&t=1192725)
// by Henrique Becker on 20/07/2016.

template<typename T> struct type_name {
    static std::string get() {
        if (std::is_enum<T>::value) return "enum";
        if (std::is_class<T>::value) return "class";
        if (std::is_union<T>::value) return "union";
        return "unknown";
    }
};

#define TYPE_NAME(T) template<> struct type_name<T> { static std::string get() { return #T; } }
TYPE_NAME(char);
TYPE_NAME(unsigned char);
TYPE_NAME(signed char);
TYPE_NAME(char16_t);
TYPE_NAME(char32_t);
TYPE_NAME(bool);
TYPE_NAME(unsigned int);
TYPE_NAME(int);
TYPE_NAME(unsigned short);
TYPE_NAME(unsigned long);
TYPE_NAME(unsigned long long);
TYPE_NAME(long);
TYPE_NAME(long long);
TYPE_NAME(short);
TYPE_NAME(wchar_t);
TYPE_NAME(float);
TYPE_NAME(double);
TYPE_NAME(long double);
TYPE_NAME(void);
#undef TYPE_NAME

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

// take care of annoying rule that array of cv-qualified type is the same as cv-qualified array
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

int main() {
    std::cout << type_name<int*>::get() << '\n';
    std::cout << type_name<int[10]>::get() << '\n';
    std::cout << type_name<int const volatile()>::get() << '\n';
    struct S {};
    std::cout << type_name<int volatile const (* volatile (&(S::*)(int)))[10]>::get() << '\n';
}
