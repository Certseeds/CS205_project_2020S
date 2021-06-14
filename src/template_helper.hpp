/*  CS205_C_CPP
    Copyright (C) 2020  nanoseeds

    CS205_C_CPP is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    CS205_C_CPP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    */
/**
 * @Github: https://github.com/Certseeds/CS205_C_CPP
 * @Organization: SUSTech
 * @Author: nanoseeds
 * @Date: 2020-05-24 17:56:20
 * @LastEditors  : nanoseeds
 */
#ifndef CS205_C_CPP_CS205_PROJECT_2020S_SRC_TEMPLATE_HELPER_H
#define CS205_C_CPP_CS205_PROJECT_2020S_SRC_TEMPLATE_HELPER_H

#include <type_traits>

// #define Multiply_Result_t_Macro decltype(std::declval<T1>() * std::declval<T2>())

template<typename T1, typename T2>
constexpr bool is_same() { return std::is_same<T1, T2>::value; }

template<typename T>
concept is_complex = requires(T f) {
    f.real();
    f.imag();
    f /= 1.0f;
    std::constructible_from<T, typename T::value_type, typename T::value_type>;
    std::same_as<std::complex<typename T::value_type>, T>;
};
template<typename T>
concept not_complex = requires(T f) {
    !is_complex<T>;
};
template<typename T>
struct complex_inside_type : std::false_type {
    //static_assert(!is_complex<T>(), "complex_inside_type: std::false_type");
    using Type = T;
};
template<typename T>
struct complex_inside_type<std::complex<T>> : std::true_type {
    //static_assert(is_complex<T>(), "complex_inside_type<std::complex<T>>: std::true_type");
    using Type = T;
};

template<typename T>
using complex_inside_type_t = typename complex_inside_type<T>::Type;

template<typename T1, typename T2>
struct Minus_Result {
    using Type = decltype(std::declval<T1>() - std::declval<T2>());
};

template<typename T1, typename T2>
using Minus_Result_t = typename Minus_Result<T1, T2>::Type;


template<typename T1, typename T2>
struct Add_Result {
    using Type = decltype(std::declval<T1>() + std::declval<T2>());
};

template<typename T1, typename T2>
using Add_Result_t = typename Add_Result<T1, T2>::Type;


template<typename T1, typename T2>
struct Multiply_Result {
    using Type = decltype(std::declval<T1>() * std::declval<T2>());
};

template<typename T1, typename T2>
using Multiply_Result_t = typename Multiply_Result<T1, T2>::Type;

template<typename T1, typename T2>
struct Divide_Result {
    using Type = decltype(std::declval<T1>() / std::declval<T2>());
};
template<typename T1, typename T2>
using Divide_Result_t = typename Divide_Result<T1, T2>::Type;

template<typename T1>
using divide_t = decltype(std::declval<T1>() / std::declval<double>());


template<typename T>
concept has_conj = requires(T f) {
    f.conj();
};

template<typename T>
T from_char_array(unsigned char const *buffer) {
    T will_return;
    auto *dp = reinterpret_cast<unsigned char *>(&will_return);
    std::copy(buffer, buffer + sizeof(T), dp);
    return will_return;
}

#define MY_IF0(...) typename std::enable_if<(bool)(__VA_ARGS__), int >::type
#define MY_IF(...) MY_IF0(__VA_ARGS__) = 0
#endif //CS205_C_CPP_CS205_PROJECT_2020S_SRC_TEMPLATE_HELPER_H
