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
 * @Date: 2020-05-07 11:30:27 
 * @LastEditors  : nanoseeds
 */
#include "./../../catch.hpp"
#include "./Matrix.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <complex>

using std::cout;
using std::endl;
using std::string;
using Catch::Matchers::Equals;


TEST_CASE("test 0", "[test1]") {
    auto function = [](const std::vector<int32_t> &t1,
                       const std::vector<int32_t> &t2) -> bool {
        assert(t1.size() == t2.size());
        for (uint32_t i = 0; i < t1.size(); i++) {
            if (t1[i] > t2[i]) {
                return true;
            } else if (t1[i] < t2[i]) {
                return false;
            }
        }
        return true;
    };
    REQUIRE(function({CATCH_VERSION_MAJOR, CATCH_VERSION_MINOR, CATCH_VERSION_PATCH}, {2, 12, 1}));
    REQUIRE(function({2, 13, 0}, {2, 12, 1}));
    REQUIRE(function({3, 0, 12}, {2, 12, 1}));
}

TEST_CASE("default nonparametric Constructor", "[test 1]") {
    Matrix<int32_t> m1;
    cout << m1;
    Matrix<int32_t> *m2_p = new Matrix<int32_t>();
    cout << *m2_p;
    delete m2_p;
}

TEST_CASE("default parameterized Constructor", "[test 1]") {
    Matrix<int64_t> m1(2, 2);
    cout << m1;
    Matrix<uint32_t> m2(5, 4);
    cout << m2;
    Matrix<std::complex<int32_t>> mat_comp(3, 4);
    cout << mat_comp;
}

TEST_CASE("Copy Constructor", "[test 1]") {
    Matrix<int64_t> m1(2, 2);
    Matrix<uint32_t> m2(5, 4);
    Matrix<std::complex<int32_t>> mat_comp(3, 4);
    Matrix<int64_t> m3(m1);
    cout << m3;
    Matrix<uint32_t> m4 = m2;
    cout << m4;
    CHECK(!m1.isEmpty());
    CHECK(!m2.isEmpty());
    CHECK(!m3.isEmpty());
    CHECK(!m4.isEmpty());
}

TEST_CASE("Move Constructor", "[test 1]") {
    Matrix<uint32_t> m1(5, 4);
    cout << m1;
    CHECK(!m1.isEmpty());
    Matrix<uint32_t> m2 = std::move(m1);
    cout << m2;
    CHECK(m1.isEmpty());
    CHECK(!m2.isEmpty());
    Matrix<int32_t> m3 = std::move(Matrix<int32_t>(2, 3));
    cout << m3;
    CHECK(!m3.isEmpty());
}

TEST_CASE("Copy Assignment operator", "[test 1]") {
    Matrix<uint64_t> m1;
    Matrix<uint64_t> m2;
    m2 = m1;
    CHECK(!m1.isEmpty());
    CHECK(!m2.isEmpty());
}

TEST_CASE("Move Assignment operator", "[test 2]") {
    Matrix<uint64_t> m1(3, 7);
    Matrix<uint64_t> m2(7,3);
    CHECK(!m1.isEmpty());
    CHECK(!m2.isEmpty());
    m2 = std::move(m1);
    CHECK(m1.isEmpty());
    CHECK(!m2.isEmpty());
}