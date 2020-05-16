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
#include <random>

using std::cout;
using std::endl;
using std::string;
using Catch::Matchers::Equals;


TEST_CASE("test 0", "[test 0]") {
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
    CHECK(!m1.is_empty());
    CHECK(!m2.is_empty());
    CHECK(!m3.is_empty());
    CHECK(!m4.is_empty());
}

TEST_CASE("Move Constructor", "[test 1]") {
    Matrix<uint32_t> m1(5, 4);
    cout << m1;
    CHECK(!m1.is_empty());
    Matrix<uint32_t> m2 = std::move(m1);
    cout << m2;
    CHECK(m1.is_empty());
    CHECK(!m2.is_empty());
    Matrix<int32_t> m3 = std::move(Matrix<int32_t>(2, 3));
    cout << m3;
    CHECK(!m3.is_empty());
}

TEST_CASE("test vector<vector<T>> constructor") {
    vector<vector<int32_t>> v1 = {{1, 2},
                                  {3, 4},
                                  {5, 6},
                                  {7, 8},
                                  {9, 10}};
    vector<vector<int32_t>> v2 = {{1, 1, 4},
                                  {5, 1, 4}};
    vector<vector<int32_t>> v3 = {{1, 1, 9},
                                  {1, 2, 0},
                                  {1, 1, 4}};
    vector<vector<int32_t>> v4 = {{1, 10},
                                  {3, 4},
                                  {5, 6},
                                  {7, 8},
                                  {9, 10}};
    auto vci = v4[1].begin();
    Matrix<int32_t> m1(v1);
    Matrix<int32_t> m2 = Matrix<int32_t>(v2);
    auto m3 = Matrix<int32_t>(v3);
    auto m4((Matrix<int32_t>(v4)));
    cout << m4;
}

TEST_CASE("Copy Assignment operator", "[test 1]") {
    Matrix<uint64_t> m1(3, 4);
    Matrix<uint64_t> m2;
    CHECK(!m1.is_empty());
    CHECK(m2.is_empty());
    m2 = m1;
    CHECK(!m1.is_empty());
    CHECK(!m2.is_empty());
}

TEST_CASE("Move Assignment operator", "[test 1]") {
    Matrix<uint64_t> m1(3, 7);
    Matrix<uint64_t> m2(7, 3);
    CHECK(!m1.is_empty());
    CHECK(!m2.is_empty());
    m2 = std::move(m1);
    CHECK(m1.is_empty());
    CHECK(!m2.is_empty());
}

TEST_CASE("initializer_list", "[test 1]") {
    Matrix<int64_t> m1{{1, 2, 3, 4, 54},
                       {5, 2, 3, 3, 41}};
    CHECK(!m1.is_empty());
    cout << m1 << endl;
    Matrix<int64_t> m2 = {{1, 1, 4, 5, 1, 4},
                          {1, 9, 1, 9, 8, 1, 0}};
    CHECK(m2.is_empty());
    Matrix<int64_t> m3{};
    CHECK(m3.is_empty());
    cout << m3;
    Matrix<int64_t> m4 = {1, 2, 3, 4, 5};
    CHECK(!m4.is_empty());
    cout << m4;
    Matrix<int64_t> m5 = {1};
    CHECK(!m5.is_empty());
    cout << m5;
    // can not Matrix<int64_t> m6 = {};
}

TEST_CASE("no vectors negative", "[test 1]") {
    Matrix<std::complex<int32_t>> m1(5, 2);
    Matrix<int32_t> m2(4, 6);
    Matrix<int_fast16_t> m3(1, 4);
    Matrix<double> m4(0, -1);
    Matrix<float> m5(-1, -1);
    Matrix<int_least32_t> m6(-2, 9);
    CHECK(m1.rows() >= 0);
    CHECK(m1.cols() >= 0);
    CHECK(m2.rows() >= 0);
    CHECK(m2.cols() >= 0);
    CHECK(m3.rows() >= 0);
    CHECK(m3.cols() >= 0);
    CHECK(m4.rows() >= 0);
    CHECK(m4.cols() >= 0);
    CHECK(m5.rows() >= 0);
    CHECK(m5.cols() >= 0);
    CHECK(m6.rows() >= 0);
    CHECK(m6.cols() >= 0);
}

TEST_CASE("zeros and ones", "[test 2]") {
    Matrix<int16_t> m1 = Matrix<int16_t>::zeros(3, 3);
    cout << m1;
    Matrix<int32_t> m2 = Matrix<int32_t>::ones(4, 4);
    cout << m2;
    Matrix<int64_t> m3 = Matrix<int64_t>::values(5, 5, 114);
    cout << m3;
}

TEST_CASE("eye and eye_value", "[test 2]") {
    Matrix<int16_t> m1 = Matrix<int16_t>::eye(4);
    cout << m1;
    Matrix<int32_t> m2 = Matrix<int32_t>::eye_value(4, 4);
    cout << m2;
}

TEST_CASE("test is empty", "[test 2]") {
    CHECK(!Matrix<int16_t>::eye(4).is_empty());
    CHECK(!Matrix<int32_t>::eye_value(4, 4).is_empty());
    CHECK(Matrix<int16_t>::zeros(0, 0).is_empty());
    CHECK(Matrix<int32_t>::ones(0, 0).is_empty());
}

// UNTODO test the equal function.
TEST_CASE("size_equal function", "[test 2]") {
    auto m1 = Matrix<std::complex<int32_t>>::eye_value(5, 2);
    auto m2 = Matrix<int32_t>::zeros(5, 5);
    CHECK(size_equal(m1, m2));
    CHECK(size_equal(m1, m1));
    CHECK(size_equal(m2, m2));
    auto m3 = Matrix<int_fast16_t>::ones(3, 4);
    auto empty = Matrix<double>(0, 0);
    CHECK(!size_equal(m3, empty));
    CHECK(size_equal(m3, m3));
    CHECK(!size_equal(empty, empty));
    auto empty_2 = Matrix<float>(0, 0);
    CHECK(!size_equal(empty, empty_2));
    auto temp = std::move(empty_2);
    CHECK(!size_equal(temp, empty));
    CHECK(!size_equal(temp, empty_2));
    CHECK(!size_equal(m1, empty));
    CHECK(!size_equal(m2, empty));
}

TEST_CASE("inside_equal function", "[test 2]") {
    Matrix<int32_t> m1 = {{1, 2},
                          {3, 4},
                          {5, 6},
                          {7, 8},
                          {9, 10}};
    Matrix<int32_t> m2 = {{1, 1, 4},
                          {5, 1, 4}};
    Matrix<int_fast16_t> m3 = {{1, 1, 9},
                               {1, 2, 0},
                               {1, 1, 4}};
    Matrix<int32_t> m4 = {{1, 10},
                          {3, 4},
                          {5, 6},
                          {7, 8},
                          {9, 10}};
    Matrix<int32_t> m5 = {{1, 3, 4},
                          {5, 1, 4}};
    Matrix<int_fast16_t> m6 = {{1, 1, 9},
                               {1, 2, 0},
                               {4, 1, 4}};
    CHECK(!Matrix<int32_t>::inside_equal(m1, m4));
    CHECK(!Matrix<int32_t>::inside_equal(m2, m5));
    CHECK(!Matrix<int_fast16_t>::inside_equal(m3, m6));
}
// UNTODO and rows functiont
TEST_CASE("rows function", "[test 2]") {
    auto m1 = Matrix<std::complex<int32_t>>::eye_value(5, 2);
    auto m2 = Matrix<int32_t>::zeros(5, 5);
    auto m3 = Matrix<int_fast16_t>::ones(3, 4);
    auto empty = Matrix<double>(0, 0);
    auto empty_2 = Matrix<float>(0, 0);
    auto m4 = Matrix<int_least32_t>(3, 6);
    CHECK(m1.rows() == 5);
    CHECK(m2.rows() == 5);
    CHECK(m3.rows() == 3);
    CHECK(empty.rows() == 0);
    CHECK(empty_2.rows() == 0);
    CHECK(m4.rows() == 3);
    CHECK(Matrix<u_int64_t>::ones(3, 4).rows() == 3);
    CHECK(Matrix<int16_t>::eye(4).rows() == 4);
    CHECK(Matrix<int32_t>::eye_value(4, 4).rows() == 4);
    CHECK(!Matrix<int16_t>::zeros(0, 0).rows());
    CHECK(!Matrix<int32_t>::ones(0, 0).rows());
    CHECK(!Matrix<int32_t>::ones(0, 5).rows());
    CHECK(Matrix<int32_t>::ones(5, 0).rows() == 5);
}
// UNTODO test for cols funciton
TEST_CASE("cols function", "[test 2]") {
    auto m1 = Matrix<std::complex<int32_t>>::eye_value(5, 2);
    auto m2 = Matrix<int32_t>::zeros(5, 5);
    auto m3 = Matrix<int_fast16_t>::ones(3, 4);
    auto empty = Matrix<double>(0, 0);
    auto empty_2 = Matrix<float>(0, 0);
    auto m4 = Matrix<int_least32_t>(3, 6);
    CHECK(m1.cols() == 5);
    CHECK(m2.cols() == 5);
    CHECK(m3.cols() == 4);
    CHECK(empty.cols() == 0);
    CHECK(empty_2.cols() == 0);
    CHECK(m4.cols() == 6);
    CHECK(Matrix<u_int64_t>::ones(3, 4).cols() == 4);
    CHECK(Matrix<int16_t>::eye(4).cols() == 4);
    CHECK(Matrix<int32_t>::eye_value(4, 4).cols() == 4);
    CHECK(!Matrix<int16_t>::zeros(0, 0).cols());
    CHECK(!Matrix<int32_t>::ones(0, 0).cols());
    CHECK(!Matrix<int32_t>::ones(0, 5).cols());
    CHECK(!Matrix<int32_t>::ones(5, 0).cols());
}

// UNTODO test transpose
TEST_CASE("transpose", "[test 3]") {
    Matrix<std::complex<int32_t>> m1(5, 2);
    Matrix<int32_t> m2(4, 6);
    Matrix<int_fast16_t> m3(1, 4);
    Matrix<double> m4(0, -1);
    Matrix<float> m5(-1, -1);
    Matrix<int_least32_t> m6(-2, 9);
    CHECK(Matrix<u_int64_t>::ones(3, 4).transpose().rows() == 4);
    CHECK(Matrix<int16_t>::eye(4).transpose().rows() == 4);
    CHECK(Matrix<int32_t>::eye_value(4, 4).transpose().rows() == 4);
    CHECK(!Matrix<int16_t>::zeros(0, 0).transpose().rows());
    CHECK(!Matrix<int32_t>::ones(0, 0).transpose().rows());
    CHECK(!Matrix<int32_t>::ones(0, 5).transpose().rows());
    CHECK(!Matrix<int32_t>::ones(5, 0).transpose().rows());
    CHECK(Matrix<u_int64_t>::ones(3, 4).transpose().cols() == 3);
    CHECK(Matrix<int16_t>::eye(4).transpose().cols() == 4);
    CHECK(Matrix<int32_t>::eye_value(4, 4).transpose().cols() == 4);
    CHECK(!Matrix<int16_t>::zeros(0, 0).transpose().cols());
    CHECK(!Matrix<int32_t>::ones(0, 0).transpose().cols());
    CHECK(!Matrix<int32_t>::ones(0, 5).transpose().cols());
    CHECK(!Matrix<int32_t>::ones(5, 0).transpose().cols());
    CHECK(m1.transpose().cols() == 5);
    CHECK(m1.transpose().rows() == 2);
    CHECK(m2.transpose().cols() == 4);
    CHECK(m2.transpose().rows() == 6);
    CHECK(m3.transpose().cols() == 1);
    CHECK(m3.transpose().rows() == 4);
    CHECK(m4.transpose().cols() == 0);
    CHECK(m4.transpose().rows() == 0);
    CHECK(m5.transpose().cols() == 0);
    CHECK(m5.transpose().rows() == 0);
    CHECK(m6.transpose().cols() == 0);
    CHECK(m6.transpose().rows() == 0);
}

TEST_CASE("conj test", "[test 3]") {
    Matrix<std::complex<int32_t>> m1 = {{std::complex(2, 3),  std::complex(3, 5)},
                                        {std::complex(-2, 3), std::complex(-2, -3)},
                                        {std::complex(2, -3), std::complex(2, 3)}};
    Matrix<std::complex<double>> m2 = {{std::complex(2.5f, -3.0f), std::complex(3.0f, 5.0f)},
                                       {std::complex(-2.3f, 3.2f), std::complex(-2.0f, -3.0f)},
                                       {std::complex(2.0f, -3.0f), std::complex(2.5f, 3.0f)}};
    Matrix<std::complex<int32_t>> m1_conj = {{std::complex(2, -3),  std::complex(3, -5)},
                                             {std::complex(-2, -3), std::complex(-2, 3)},
                                             {std::complex(2, 3),   std::complex(2, -3)}};
    Matrix<std::complex<double>> m2_conj = {{std::complex(2.5f, 3.0f),   std::complex(3.0f, -5.0f)},
                                            {std::complex(-2.3f, -3.2f), std::complex(-2.0f, 3.0f)},
                                            {std::complex(2.0f, 3.0f),   std::complex(2.5f, -3.0f)}};
    Matrix<int32_t> m3 = {{2,  12},
                          {6,  8},
                          {10, 12},
                          {14, 16},
                          {18, 20}};
    cout << conj(m1);
    cout << conj(m2);
    cout << conj(m3);
    CHECK(Matrix<std::complex<int32_t>>::inside_equal(conj(m1), m1_conj));
    CHECK(Matrix<std::complex<double>>::inside_equal(conj(m2), m2_conj));
    CHECK(Matrix<int32_t>::inside_equal(conj(m3), m3));
}
// UNTODO for operator+
TEST_CASE("operator plus", "[test 3]") {
    Matrix<int32_t> m1 = {{1, 2},
                          {3, 4},
                          {5, 6},
                          {7, 8},
                          {9, 10}};
    Matrix<int32_t> m4 = {{1, 10},
                          {3, 4},
                          {5, 6},
                          {7, 8},
                          {9, 10}};
    Matrix<int32_t> m7 = {{2,  12},
                          {6,  8},
                          {10, 12},
                          {14, 16},
                          {18, 20}};
    Matrix<int32_t> m2 = {{1, 1, 4},
                          {5, 1, 4}};
    Matrix<int32_t> m5 = {{1, 3, 4},
                          {5, 1, 4}};
    Matrix<int32_t> m8 = {{2,  4, 8},
                          {10, 2, 8}};
    Matrix<int_fast16_t> m3 = {{1, 1, 9},
                               {1, 2, 0},
                               {1, 1, 4}};
    Matrix<int_fast16_t> m6 = {{1, 1, 9},
                               {1, 1, 0},
                               {2, 1, 6}};
    Matrix<int_fast16_t> m9 = {{2, 2, 18},
                               {2, 3, 0},
                               {3, 2, 10}};
    CHECK(Matrix<int32_t>::inside_equal(m1 + m4, m7));
    CHECK(Matrix<int32_t>::inside_equal(m2 + m5, m8));
    CHECK(Matrix<int_fast16_t>::inside_equal(m3 + m6, m9));
}
// UNTODO for operator-
TEST_CASE("operator minus", "[test 3]") {
    Matrix<int32_t> m1 = {{1, 2},
                          {3, 4},
                          {5, 6},
                          {7, 8},
                          {9, 10}};
    Matrix<int32_t> m4 = {{1, 10},
                          {3, 4},
                          {5, 6},
                          {7, 8},
                          {9, 10}};
    Matrix<int32_t> m7 = {{0, -8},
                          {0, 0},
                          {0, 0},
                          {0, 0},
                          {0, 0}};
    Matrix<int32_t> m2 = {{1, 1, 4},
                          {5, 1, 4}};
    Matrix<int32_t> m5 = {{1, 3, 4},
                          {5, 1, 4}};
    Matrix<int32_t> m8 = {{0, -2, 0},
                          {0, 0,  0}};
    Matrix<int_fast16_t> m3 = {{1, 1, 9},
                               {1, 2, 0},
                               {1, 1, 4}};
    Matrix<int_fast16_t> m6 = {{1, 1, 9},
                               {1, 1, 0},
                               {2, 1, 6}};
    Matrix<int_fast16_t> m9 = {{0,  0, 0},
                               {0,  1, 0},
                               {-1, 0, -2}};
    CHECK(Matrix<int32_t>::inside_equal(m1 - m4, m7));
    CHECK(Matrix<int32_t>::inside_equal(m2 - m5, m8));
    CHECK(Matrix<int_fast16_t>::inside_equal(m3 - m6, m9));
}
// UNTODO for operator* (Matrix & number)
TEST_CASE("operator multiply_matrix_number", "[test 3]") {
    Matrix<int32_t> m1 = {{1, 2},
                          {3, 4},
                          {5, 6},
                          {7, 8},
                          {9, 10}};
    Matrix<int32_t> m4 =
            {{2,  4},
             {6,  8},
             {10, 12},
             {14, 16},
             {18, 20}};
    Matrix<int32_t> m2 = {{1, 1, 4},
                          {5, 1, 4}};
    Matrix<int32_t> m5 = {{2,  2, 8},
                          {10, 2, 8}};
    Matrix<int_fast16_t> m3 = {{1, 1, 9},
                               {1, 2, 0},
                               {1, 1, 4}};
    Matrix<int_fast16_t> m6 = {{0, 0, 0},
                               {0, 0, 0},
                               {0, 0, 0}};
    CHECK(Matrix<int32_t>::inside_equal(m1 * 2, m4));
    CHECK(Matrix<int32_t>::inside_equal(2 * m1, m4));
    CHECK(Matrix<int32_t>::inside_equal(m2 * 2, m5));
    CHECK(Matrix<int32_t>::inside_equal(2 * m2, m5));
    CHECK(Matrix<int_fast16_t>::inside_equal(m3 * 0, m6));
    CHECK(Matrix<int_fast16_t>::inside_equal(0 * m3, m6));
    Matrix<std::complex<int32_t>> mc1 = Matrix<std::complex<int32_t>>::ones(5, 4);
    Matrix<int32_t> m10 = Matrix<int32_t>::ones(5, 4);
    mc1 = mc1 * std::complex<int32_t>(3, 4);
    auto m11 = m10 * std::complex<int32_t>(3, 4);
    auto m12 = std::complex<int32_t>(3, 4) * m10;
    cout << mc1;
    cout << m11;
    cout << m12;
}
// TODO for operator* (matrix * matrix)
TEST_CASE("operator multiply_matrix_matrix_same", "[test 3]") {
    std::uniform_int_distribution<int32_t> range1(3, 6);
    std::random_device r;
    std::default_random_engine e1(r());
    int32_t t1 = range1(e1);
    int32_t t2 = range1(e1);
    int32_t t3 = range1(e1);
    vector<vector<int32_t>> vec1(t1, vector<int32_t>(t2));
    vector<vector<int32_t>> vec2(t2, vector<int32_t>(t3));
    for (auto &item : vec1) {
        for (auto &item2 : item) {
            item2 = range1(e1);
        }
    }
    for (auto &item : vec2) {
        for (auto &item2 : item) {
            item2 = range1(e1);
        }
    }
    cout << Matrix<int32_t>(vec1) * Matrix<int32_t>(vec2);
    vector<int32_t> vec3 = {1, 2, 3, 4};
    vector<int32_t> vec4 = {1, 2, 3, 4};
    cout << std::inner_product(vec1[0].cbegin(), vec1[0].cend(), vec1[0].cbegin(), 0);
    vector<vector<int32_t>> vec5 = {{1, 2, 3},
                                    {4, 5, 6}};
    vector<vector<std::complex<int32_t>>> vec_c6 = {{std::complex<int32_t>(9, 10),  std::complex<int32_t>(11, 12)},
                                                    {std::complex<int32_t>(13, 14), std::complex<int32_t>(15, 61)}};
    //cout << Matrix<int32_t>(vec5) * Matrix<std::complex<int32_t>>(vec_c6);

}
// UNTODO for operator mul
TEST_CASE("operator mul", "[test 3]") {
    Matrix<int32_t> m1 = {{1, 2},
                          {3, 4},
                          {5, 6},
                          {7, 8},
                          {9, 10}};
    Matrix<int32_t> m4 = {{1, 10},
                          {3, 4},
                          {5, 6},
                          {7, 8},
                          {9, 10}};
    Matrix<int32_t> m7 = {{1,  20},
                          {9,  16},
                          {25, 36},
                          {49, 64},
                          {81, 100}};
    Matrix<int32_t> m2 = {{1, 1, 4},
                          {5, 1, 4}};
    Matrix<int32_t> m5 = {{1, 3, 4},
                          {5, 1, 4}};
    Matrix<int32_t> m8 = {{1,  3, 16},
                          {25, 1, 16}};
    Matrix<int_fast16_t> m3 = {{1, 1, 9},
                               {1, 2, 0},
                               {1, 1, 4}};
    Matrix<int_fast16_t> m6 = {{1, 1, 9},
                               {1, 1, 0},
                               {2, 1, 6}};
    Matrix<int_fast16_t> m9 = {{1, 1, 81},
                               {1, 2, 0},
                               {2, 1, 24}};
    CHECK(Matrix<int32_t>::inside_equal(m1.mul(m4), m7));
    CHECK(Matrix<int32_t>::inside_equal(m2.mul(m5), m8));
    CHECK(Matrix<int_fast16_t>::inside_equal(m3.mul(m6), m9));
}
// UNTODO for operator dot
TEST_CASE("operator dot", "[test 3]") {
    Matrix<int32_t> m1 = {{1, 2, 3},
                          {4, 5, 6}};
    Matrix<int32_t> m2 = {{1, 2, 3},
                          {4, 5, 6}};
    cout << m1.dot(m2);
    CHECK(m1.dot(m2) == 91);
}

TEST_CASE("Kronecker product", "[test 3]") {
    Matrix<int32_t> m1 = {{1, 2},
                          {3, 1}};
    Matrix<int32_t> m2 = {{0, 3},
                          {2, 1}};
    Matrix<int32_t> result = {{0, 3, 0, 6},
                              {2, 1, 4, 2},
                              {0, 9, 0, 3},
                              {6, 3, 2, 1}};
    auto m3 = kron(m1, m2);
    CHECK(Matrix<int32_t>::inside_equal(m3, result));
    cout << m3;
}

TEST_CASE("vector * matrix", "[test 3]") {
    Matrix<int32_t> m1 = {{0, 3, 0},
                          {2, 1, 4}};
    vector<int32_t> vec1 = {4, 5, 6};
    vector<int32_t> vec2 = {4, 5};
    vector<std::complex<int32_t>> vec3 = {std::complex<int64_t>(3, 4), std::complex<int64_t>(5, 6)};
    vector<std::complex<int32_t>> vec4 = {std::complex<int64_t>(3, 4), std::complex<int64_t>(5, 6),
                                          std::complex<int64_t>(7, 98)};
    cout << m1 * vec1;
    cout << vec2 * m1;
    cout << m1 * vec4;
    // decltype(std::declval<std::complex<int64_t>>() * std::declval<int64_t>()) temp2 = {1, 2};
    cout << vec3 * m1;
}

TEST_CASE("vector cross vector", "[test 3]") {
    vector<int32_t> vec1 = {4, 5, 6};
    vector<std::complex<int32_t>> vec2 = {std::complex<int64_t>(3, 4), std::complex<int64_t>(5, 6),
                                          std::complex<int64_t>(7, 98)};
    auto vec3 = cross(vec1, vec2);
}

TEST_CASE("test for trace", "[test 5]") {
    Matrix<int32_t> m1 = {{1, 2},
                          {3, 1}};
    Matrix<int32_t> m2 = {{0, 3},
                          {2, 1}};
    Matrix<int32_t> result = {{0, 3, 0, 6},
                              {2, 1, 4, 2},
                              {0, 9, 0, 3},
                              {6, 3, 2, 1}};
    cout << m1.trace() << "\n";
    cout << m1.determinant() << "\n";
    cout << m2.trace() << "\n";
    cout << m2.determinant() << "\n";
    cout << result.trace() << "\n";
    cout << result.determinant() << "\n";
}