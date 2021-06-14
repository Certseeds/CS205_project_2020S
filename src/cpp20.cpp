//
// Created by nanos on 2021/6/14.
//
#include <concepts>
#include <complex>
#include <iostream>

using namespace std;
namespace cpp20 {
    // `T` is not limited by any constraints.
    template<typename T>
    concept always_satisfied = true;
// Limit `T` to integrals.
    template<typename T>
    concept integral = std::is_integral_v<T>;
// Limit `T` to both the `integral` constraint and signedness.
    template<typename T>
    concept signed_integral = integral<T> && std::is_signed_v<T>;
// Limit `T` to both the `integral` constraint and the negation of the `signed_integral` constraint.
    template<typename T>
    concept unsigned_integral = integral<T> && !signed_integral<T>;

    template<typename T>
    concept is_complex = same_as<complex<typename T::value_type>, T>;

    template<typename T>
    concept realimag = requires(T f) {
        f.real();
        f.imag();
        f /= 1.0f;
        constructible_from<T, typename T::value_type, typename T::value_type>;
        same_as<complex<typename T::value_type>, T>;
    };

    template<typename T>
    struct p {
        T t1{};
        T t2{};

        p() {
            std::cout << "list zero " << std::endl;
        }

        p(T t1, T t2) {
            std::cout << "list double T " << std::endl;
        }

        p(const initializer_list<T> &list) {
            std::cout << "list const one& " << std::endl;
        }

        p(const initializer_list<initializer_list<T>> &list) {
            std::cout << "list two " << std::endl;
        }

        template<typename T1=T>
        requires (!realimag<T>)
        void link() {
            std::cout << "not " << std::endl;
        }

        template<typename T1=T>
        requires realimag<T>
        void link() {
            std::cout << "yes " << std::endl;
        }
    };

    int main() {
        std::complex<int32_t> cint{0, 0};
        int yesint = 114514;
        p<complex<int>> p1{cint, cint};
        p<int> p2{yesint, yesint};
        p1.link();
        p2.link();
        p<double> p3{1, 2, 3, 4};
        p<double> p4{{1, 2, 3, 4}};
        return 0;
    }
}

int main() {
    return cpp20::main();
}