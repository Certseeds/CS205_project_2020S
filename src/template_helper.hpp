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
#ifndef CS205_C_CPP_CS205_PROJECT_2020S_SRC_TEMPLATE_HELPER_HPP
#define CS205_C_CPP_CS205_PROJECT_2020S_SRC_TEMPLATE_HELPER_HPP

#include <type_traits>
#include <concepts>

using std::type_identity;
using std::is_same;
using std::is_same_v;

template<bool B, typename T, typename F>
struct conditional : type_identity<T> {};

template<typename T, typename F>
struct conditional<false, T, F> : type_identity<F> {};

template<bool B, typename T, typename F>
using conditional_t = typename conditional<B, T, F>::type;

template<typename T>
concept is_complex = requires(T f) {
	f.real();
	f.imag();
	std::conj(f);
	f /= 1.0f;
	std::constructible_from<T, typename T::value_type, typename T::value_type>;
	std::same_as<std::complex<typename T::value_type>, T>;
};
template<typename T>
struct complex_inside : type_identity<T> {};
template<typename T>
struct complex_inside<std::complex<T>> : type_identity<T> {};

template<typename T>
using complex_inside_t = typename complex_inside<T>::type;

template<typename T1, typename T2>
struct Minus_Result : type_identity<decltype(std::declval<T1>() - std::declval<T2>())> {
};

template<typename T1, typename T2>
using Minus_Result_t = typename Minus_Result<T1, T2>::type;


template<typename T1, typename T2>
struct Add_Result : type_identity<decltype(std::declval<T1>() + std::declval<T2>())> {
};

template<typename T1, typename T2>
using Add_Result_t = typename Add_Result<T1, T2>::type;


template<typename T1, typename T2>
struct Multiply_Result : type_identity<decltype(std::declval<T1>()*std::declval<T2>())> {
};

template<typename T1, typename T2>
using Multiply_Result_t = typename Multiply_Result<T1, T2>::type;

template<typename T1>
using divide_t = decltype(std::declval<T1>()/std::declval<double>());


template<typename T>
concept has_conj = requires(T f) { f.conj(); };

template<typename T>
T from_char_array(unsigned char const *buffer) {
	T will_return;
	auto *dp = reinterpret_cast<unsigned char *>(&will_return);
	std::copy(buffer, buffer + sizeof(T), dp);
	return will_return;
}

template<typename T, typename U, typename... Rest>
struct is_one_of : conditional_t<
	is_same_v<T, U>, std::true_type, is_one_of<T, Rest...>> {
};

template<typename T, typename U>
struct is_one_of<T, U> : conditional_t<
	is_same_v<T, U>, std::true_type, std::false_type> {
};
template<typename T>
concept opencv_type = requires(T f) {
	is_one_of<T, uint8_t, int8_t, uint16_t, int16_t, int32_t, float, double>::value;
};
#endif //CS205_C_CPP_CS205_PROJECT_2020S_SRC_TEMPLATE_HELPER_HPP
