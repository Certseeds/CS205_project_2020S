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

#define Multiply_Result_t_Macro decltype(std::declval<T1>() * std::declval<T2>())

template<typename T1, typename T2>
constexpr bool is_same() { return std::is_same<T1, T2>::value; }

template<typename T>
struct is_complex_imp : std::false_type {
};
template<typename T>
struct is_complex_imp<std::complex<T>> : std::true_type {
};

template<typename T>
constexpr bool is_complex() {
    return is_complex_imp<T>();
}

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

template<typename T, typename = void>
struct has_conj_imp : std::false_type {
};
template<typename T>
struct has_conj_imp<T, std::void_t<decltype(std::declval<T>().conj())>> : std::true_type {
};

template<typename T>
constexpr bool has_conj() {
    return has_conj_imp<T>();
}

#define MY_IF0(...) typename std::enable_if<(bool)(__VA_ARGS__), int >::type
#define MY_IF(...) MY_IF0(__VA_ARGS__) = 0
#endif //CS205_C_CPP_CS205_PROJECT_2020S_SRC_TEMPLATE_HELPER_H
