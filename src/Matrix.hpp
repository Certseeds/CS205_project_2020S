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
 * @Date: 2020-05-07 15:46:11 
 * @LastEditors  : nanoseeds
 */
#ifndef CS205_C_CPP_CS205_PROJECT_2020S_SRC_MATRIX_HPP
#define CS205_C_CPP_CS205_PROJECT_2020S_SRC_MATRIX_HPP

#include <vector>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <functional>
#include <complex>
#include <utility>
#include <type_traits>
#include <iostream>
#include <cstring>

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

#define MY_IF0(...) typename std::enable_if<(bool)(__VA_ARGS__), int >::type
#define MY_IF(...) MY_IF0(__VA_ARGS__) = 0

using std::vector;
using std::endl;

/** @related: get from https://blog.csdn.net/qq_31175231/article/details/77479692
 */

template<typename T>
class Matrix {
private:
    vector<vector<T>> vec;
public:

//explicit Matrix() { this->vec = vector<vector<T>>(1, vector<T>(1)); };
    explicit Matrix() : Matrix(0, 0) {};

    explicit Matrix(int32_t rows, int32_t cols);

// 拷贝构造函数 Copy Constructor
    Matrix(const Matrix<T> &mat);

// 移动构造函数
    Matrix(Matrix<T> &&mat) noexcept;

// 拷贝赋值运算符 Copy Assignment operator
    Matrix<T> &operator=(const Matrix<T> &mat);

// 移动赋值运算符 Move Assignment operator
    Matrix<T> &operator=(Matrix<T> &&mat)
    noexcept;

    Matrix<T>(const vector<vector<T>> &vec);

    Matrix<T>(vector<vector<T>> &&vec);

    Matrix<T>(const std::initializer_list<std::initializer_list<T>> &list);

// only for
    Matrix<T>(const std::initializer_list<T> &list);

    inline T get_inside(int32_t rows, int32_t cols) const;

    static Matrix<T> zeros(int32_t rows, int32_t cols);

    static Matrix<T> ones(int32_t rows, int32_t cols);

    static Matrix<T> values(int32_t rows = 0, int32_t cols = 0, T = static_cast<T>(0));

    static Matrix<T> eye(int32_t s);

    static Matrix<T> eye_value(int32_t s, T t);

    static bool inside_equal(const Matrix<T> &mat1, const Matrix<T> &mat2);

    friend std::ostream &operator<<(std::ostream &output, const Matrix<T> &mat) {
        for (const auto &i: mat.vec) {
            for (const auto &j:i) {
                output << j << " ";
            }
            output << endl;
        }
        return output;
    }

    inline int32_t rows() const;

    inline int32_t cols() const;

    Matrix<T> transpose() const;

    ~Matrix() = default;

    inline bool is_empty() const;

    inline bool is_square() const { return this->rows() == this->cols(); };

    inline static Matrix<T>
    operator_table(const Matrix<T> &mat1, const Matrix<T> &mat2,
                   const std::function<T(const T &t1, const T &t2)> &func);

    inline static Matrix<T>
    operator_table(const Matrix<T> &mat1, const std::function<T(const T &t1)> &func);

    typename vector<T>::const_iterator get_row_iter_begin(int32_t rows) const;

    typename vector<T>::const_iterator get_row_iter_end(int32_t rows) const;

    /**@related: https://stackoverflow.com/questions/30736951/templated-class-check-if-complex-at-compile-time
     * */
    template<typename T1 = T, MY_IF(is_complex<T1>())>
    Matrix<T1> conj() const {
        vector<vector<T1>> vvc(this->rows(), vector<T1>(this->cols()));
        for (int i = 0; i < this->rows(); ++i) {
            for (int j = 0; j < this->cols(); ++j) {
                vvc[i][j] = std::conj(this->vec[i][j]);
            }
        }
        return Matrix<T1>(std::move(vvc));
    }

    T max() const;

    T min() const;

    T sum() const;

    template<typename T1 = T, MY_IF(!is_complex<T1>())>
    auto avg() -> decltype(std::declval<T1>() / std::declval<double>()) const {
        T sums = this->sum();
        return sums / static_cast<double_t>(this->rows() * this->cols());
    }

    //template< class T1 = T ,typename std::enable_if<is_complex<T1>(),int>::type = 0>
    template<typename T1 = T, MY_IF(is_complex<T1>())>
    T1 avg() const {
        T1 sums = this->sum();
        double_t size = static_cast<double_t>(this->rows() * this->cols());
        return T1 (sums.real() / size, sums.imag() / size);
    }

    T row_max(int32_t row) const;

    T col_max(int32_t col) const;

    T row_min(int32_t row) const;

    T col_min(int32_t col) const;

    T row_sum(int32_t row) const;

    T col_sum(int32_t col) const;

    template<typename T1 = T, MY_IF(!is_complex<T1>())>
    auto row_avg(int32_t row) -> decltype(std::declval<T1>() / std::declval<double>()) const {
        if (row > this->rows()) {
            return -1;
        }
        return this->row_sum(row) / static_cast<double_t>(this->cols());
//        return 1;
    }

    template<typename T1 = T, MY_IF(is_complex<T1>())>
    T1 row_avg(int32_t row) const {
        if (row > this->rows()) {
            return static_cast<T1>(-1);
        }
        T1 sums = this->row_sum(row);
        double_t size = static_cast<double_t>(this->cols());
        return T1 (sums.real() / size, sums.imag() / size);
        // return std::complex(static_cast<double>(2), static_cast<double >(3));
    }

    template<typename T1 = T, MY_IF(!is_complex<T1>())>
    auto col_avg(int32_t col) -> decltype(std::declval<T1>() / std::declval<double>()) const {
        if (col <= 0 || col > vec[0].size()) {
            return -1;
        }
        return this->col_sum(col) / static_cast<double>(this->rows());
    }

    template<typename T1 = T, MY_IF(is_complex<T1>())>
    T1 col_avg(int32_t col) const {
        if (col <= 0 || col > vec[0].size()) {
            return -1;
        }
        T1 sums = this->col_sum(col);
        double_t size = static_cast<double_t>(this->rows());
        return T1 ( sums.real() / size, sums.imag() / size);
    }

// Matrix_n_m, Matrix_n_m, result is Matrix_N_M
    Matrix<T> mul(const Matrix<T> &mat2);

    T dot(const Matrix<T> &mat2);

    T trace() const;

    T determinant() const;

    void set_inside(int32_t row, int32_t col, T ele);

    void set_element(T ele);

    vector<vector<double_t> > transform() const;

    Matrix<double_t > Householder(int32_t col, int32_t ele) const;

    Matrix<double_t > Hessenberg() const;

    Matrix<double_t > Givens(int32_t col, int32_t begin, int32_t end) const;

    Matrix<double_t > QR_iteration() const;

    double_t* eigenvalue() const;

    Matrix<double_t > LU_decomposition() const;

    Matrix<double_t > Gauss() const;

    Matrix<double_t > row_elimination(int32_t row, double_t ele) const;

    Matrix<double_t > row_elimination(int32_t row, int32_t col, int32_t remove_row) const;

    Matrix<T> row_exchange(int32_t row1, int32_t row2) const;

    Matrix<double_t > Elimination() const;

    Matrix<double_t > Nullspace()const;

    Matrix<double_t > eigenvector() const;

    Matrix<T> joint(Matrix<T> matrix) const;

};


template<typename T>
Matrix<T>::Matrix(const Matrix &mat) {
    this->vec = vector<vector<T >>(mat.vec);
}

template<typename T>
Matrix<T>::Matrix(int32_t rows, int32_t cols) {
    rows = rows > 0 ? rows : 0;
    cols = cols > 0 ? cols : 0;
    this->vec = vector<vector<T >>(rows, vector<T>(cols));
}

template<typename T>
Matrix<T>::Matrix(Matrix &&mat)
noexcept {
    this->
            vec = std::move(mat.vec);
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &mat) {
    if (this != &mat) {
        this->vec = vector<vector<T >>(mat.vec);
    }
    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(Matrix &&mat)
noexcept {
    if (this != &mat) {
        this->
                vec = std::move(mat.vec);
    }
    return *this;
}

template<typename T>
Matrix<T>::Matrix(const vector<vector<T>> &vec) {
    this->vec = vector<vector<T >>(vec);
}

template<typename T>
Matrix<T>::Matrix(vector<vector<T >> &&vec) {
    this->vec = std::move(vec);
}

template<typename T>
Matrix<T>::Matrix(const std::initializer_list<T> &list) {
    this->vec = vector<vector<T >>(1);
    this->vec[0] = list;
}

template<typename T>
Matrix<T>::Matrix(const std::initializer_list<std::initializer_list<T>> &list) {
    this->vec = vector<vector<T >>(0);
    for (auto i = list.begin(); i < list.end() - 1; i++) {
        if ((*i).size() != (*(i + 1)).size()) {
            return;
        }
    }
    this->vec = vector<vector<T >>(list.size());
    uint32_t row = 0;
    for (const auto &i:list) {
        this->vec[row++] = i;
    }
}

template<typename T>
T Matrix<T>::get_inside(int32_t rows, int32_t cols) const {
    return this->vec[rows][cols];
}

template<typename T>
typename vector<T>::const_iterator Matrix<T>::get_row_iter_begin(int32_t rows) const {
    if (rows >= this->rows()) {
        return this->vec.front().cbegin();
    }
    return this->vec[rows].cbegin();
}

template<typename T>
typename vector<T>::const_iterator Matrix<T>::get_row_iter_end(int32_t rows) const {
    if (rows >= this->rows()) {
        return this->vec.back().end();
    }
    return this->vec[rows].end();
}

template<typename T>
inline bool Matrix<T>::is_empty() const {
    return vec.empty() || vec.front().empty();
}

template<typename T>
Matrix<T> Matrix<T>::zeros(int32_t rows, int32_t cols) {
    return Matrix<T>::values(rows, cols, static_cast<T>(0));
}

template<typename T>
Matrix<T> Matrix<T>::ones(int32_t rows, int32_t cols) {
    return Matrix<T>::values(rows, cols, static_cast<T>(1));
}

template<typename T>
Matrix<T> Matrix<T>::values(int32_t rows, int32_t cols, T t) {
    Matrix<T> will_return(rows, cols);
    for (auto &i:will_return.vec) {
        i = vector<T>(cols, t);
    }
    return will_return;
}

template<typename T>
Matrix<T> Matrix<T>::eye_value(int32_t s, T t) {
    Matrix<T> will_return(s, s);
    for (int32_t i = 0; i < s; ++i) {
        will_return.vec[i][i] = t;
    }
    return will_return;
}

template<typename T>
Matrix<T> Matrix<T>::eye(int32_t s) {
    return Matrix<T>::eye_value(s, static_cast<T>(1));
}

template<typename T>
inline Matrix<T>
Matrix<T>::operator_table(const Matrix<T> &mat1, const Matrix<T> &mat2,
                          const std::function<T(const T &t1, const T &t2)> &func) {
    if (!size_equal(mat1, mat2)) {
        exit(0);
    }
    Matrix<T> will_return(mat1.rows(), mat1.cols());
    for (int32_t i = 0; i < mat1.rows(); ++i) {
        std::transform(mat1.vec[i].begin(), mat1.vec[i].end(),
                       mat2.vec[i].begin(), will_return.vec[i].begin(), func);
    }
    return will_return;
}

template<typename T>
Matrix<T> Matrix<T>::operator_table(const Matrix<T> &mat1, const std::function<T(const T &t1)> &func) {
    Matrix<T> will_return(mat1.rows(), mat1.cols());
    for (int32_t i = 0; i < mat1.rows(); ++i) {
        std::transform(mat1.vec[i].begin(), mat1.vec[i].end(), will_return.vec[i].begin(), func);
    }
    return will_return;
}

template<typename T>
inline int32_t Matrix<T>::rows() const {
    return this->vec.size();
}

/**
 * if the it's empty,then is_empty is true,
 * */
template<typename T>
inline int32_t Matrix<T>::cols() const {
    if (this->rows() == 0) {
        return 0;
    }
    return this->vec.front().size();
}


template<typename T>
Matrix<T> Matrix<T>::transpose() const {
    if (this->is_empty()) {
        return Matrix<T>();
    }
    Matrix<T> will_return(this->cols(), this->rows());
    for (int32_t i = 0; i < this->rows(); ++i) {
        for (int32_t j = 0; j < this->cols(); ++j) {
            will_return.vec[j][i] = this->vec[i][j];
        }
    }
    return will_return;
}

/**
 * matrix + matrix, must equal.
 *  * input mat1,mat2 and will_return's type is same.
 * */
template<typename T>
Matrix<T> operator+(const Matrix<T> &mat1, const Matrix<T> &mat2) {
    return Matrix<T>::operator_table(mat1, mat2, std::plus<>());
    // [](const T &t1, const T &t2) -> T { return t1 + t2; }
}

/**
 * matrix - matrix, must equal.
 * input mat1,mat2 and will_return's type is same.
 * */
template<typename T>
Matrix<T> operator-(const Matrix<T> &mat1, const Matrix<T> &mat2) {
    return Matrix<T>::operator_table(mat1, mat2, std::minus<>());
}

/**
 * Matrix / number, return a Matrix.
 * same type.
 * */
template<typename T>
Matrix<T> operator/(const Matrix<T> &mat1, const T &t2) {
    return Matrix<T>::operator_table(mat1, [t2](const T &t1) { return t1 / t2; });
}

/**
 * matrix * number, need T1 can * T2
 * */
template<typename T1, typename T2>
auto operator*(const Matrix<T1> &mat1, const T2 &t2) -> Matrix<Multiply_Result_t_Macro> {
    vector<vector<Multiply_Result_t<T1, T2>>> temp(mat1.rows(), vector<Multiply_Result_t<T1, T2>>
            (mat1.cols()));
    for (uint32_t i = 0; i < temp.size(); ++i) {
        for (uint32_t j = 0; j < temp[i].size(); ++j) {
            temp[i][j] = mat1.get_inside(i, j) * t2;
        }
    }
    //  Matrix<decltype(mat1.get_type() * t2)> will_return(mat1.rows(), mat1.cols());
    return Matrix<Multiply_Result_t<T1, T2>>(std::move(temp));
    //for (int32_t i = 0; i < mat1.rows(); ++i) {
    //    std::transform(mat1.vec[i].begin(), mat1.vec[i].end(), will_return.vec[i].begin(), func);
    //}
    //return Matrix<decltype(mat1.get_type() * t2)>::operator_table(mat1, [&t2](const T1 &t1) { return t1 * t2; });
    // return Matrix<T1>::operator_table(mat1, [&t2](const T1 &t1) { return t1 * t2; });
}

/**
 * matrix * vector
 * @param1: Matrix<T1> m_n
 * @Param2; vector<T2> length is n
 * @return: Matrix<decltype(T1*T2)> m_1,(rows is m,cols is 1)
 * */
template<typename T1, typename T2>
auto operator*(const Matrix<T1> &mat1, const vector<T2> &t2) -> Matrix<Multiply_Result_t_Macro> {
    if (mat1.cols() != t2.size()) {
        // TODO
    }
    vector<vector<Multiply_Result_t<T1, T2>>> temp(1, vector<Multiply_Result_t<T1, T2>>
            (mat1.rows()));
    for (uint32_t i = 0; i < temp.size(); ++i) {
        temp[0][i] = std::inner_product(mat1.get_row_iter_begin(i), mat1.get_row_iter_end(i), t2.cbegin(),
                                        static_cast<Multiply_Result_t<T1, T2>>(0));
    }
    return Matrix<Multiply_Result_t<T1, T2>>(std::move(temp));
}

/**
 * matrix * vector
 * @Param1; vector<T1> length is m
 * @param2: Matrix<T2> m_n
 * @return: Matrix<decltype(T1*T2)> 1_n,(rows is 1,cols is n)
 * */
template<typename T1, typename T2>
auto
operator*(const vector<T1> &t1, const Matrix<T2> &mat2) -> Matrix<Multiply_Result_t_Macro> {
    if (t1.size() != mat2.rows()) {
        // TODO
    }
    vector<vector<Multiply_Result_t<T1, T2>>> temp(mat2.cols(), vector<Multiply_Result_t<T1, T2>>
            (1));
    auto transfor = mat2.transpose();
    for (int32_t i = 0; i < transfor.rows(); ++i) {
        temp[i][0] = std::inner_product(transfor.get_row_iter_begin(i), transfor.get_row_iter_end(i), t1.cbegin(),
                                        static_cast<Multiply_Result_t<T1, T2>>(0));
    }
    return Matrix<Multiply_Result_t<T1, T2>>(std::move(temp));
}

/**
 *  number * matrix , need T1 can * T2
 * */
template<typename T1, typename T2>
auto operator*(const T1 &t1, const Matrix<T2> &mat2) -> Matrix<Multiply_Result_t_Macro> {
    return mat2 * t1;
}

template<typename T>
Matrix<T> operator*(const Matrix<T> &mat1, const Matrix<T> &mat2) {
    if (mat1.cols() != mat2.rows()) {
        // TODO
    }
    Matrix<T> temp = mat2.transpose();
    vector<vector<T >> will_return(mat1.rows(), vector<T>(mat2.cols()));
    for (int32_t i = 0; i < mat1.rows(); ++i) {
        for (int32_t j = 0; j < mat2.cols(); ++j) {
            will_return[i][j] = std::inner_product(mat1.get_row_iter_begin(i), mat1.get_row_iter_end(i),
                                                   temp.get_row_iter_begin(j), static_cast<T>(0));
            // transform_reduce can not use.
        }
    }
    return Matrix<T>(std::move(will_return));
}

// Matrix_n_m, Matrix_n_m, result is Matrix_N_M
template<typename T>
Matrix<T> Matrix<T>::mul(const Matrix<T> &mat2) {
    return Matrix<T>::operator_table(*this, mat2, std::multiplies<>());
}

// TODO T1 and T2
template<typename T>
T Matrix<T>::dot(const Matrix<T> &mat2) {
    if (!size_equal(*this, mat2)) {
        // TODO
    }
    T will_return(0);
    for (int32_t i = 0; i < mat2.rows(); ++i) {
        will_return = std::inner_product(this->get_row_iter_begin(i), this->get_row_iter_end(i),
                                         mat2.get_row_iter_begin(i), will_return);
    }
    return will_return;
}

/**
 * kron
 * @param1: Matrix_m_n T1
 * @param2: Matrix_k_p T2
 * @return: Matrix_(mk)_(np) decltype(T1*T2)
 * */
template<typename T1, typename T2>
auto kron(const Matrix<T1> &mat1, const Matrix<T2> &mat2) -> Matrix<Multiply_Result_t_Macro> {
    vector<vector<Multiply_Result_t<T1, T2>>> will_return(
            mat1.rows() * mat2.rows(), vector<Multiply_Result_t<T1, T2>>
                    (mat1.cols() * mat2.cols()));
    for (int32_t i = 0; i < mat1.rows(); ++i) {
        for (int32_t j = 0; j < mat1.cols(); ++j) {
            for (int32_t k = 0; k < mat2.rows(); ++k) {
                for (int32_t m = 0; m < mat2.cols(); ++m) {
                    will_return[i * mat1.rows() + k][j * mat1.cols() + m] =
                            mat1.get_inside(i, j) * mat2.get_inside(k, m);
                }
            }
        }
    }
    return Matrix<Multiply_Result_t<T1, T2>>(std::move(will_return));
}

/**
 * kron
 * @param1: vector_1_3 T1
 * @param2: vector_1_3 T2
 * @return: vector_1_3 decltype(T1*T2)
 * */
template<typename T1, typename T2>
auto cross(const vector<T1> &vec1, const vector<T2> &vec2) -> vector<Multiply_Result_t_Macro> {
    vector<Multiply_Result_t<T1, T2>> will_return(3);
    if (vec1.size() == 3 && vec2.size() == 3) {
        will_return[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
        will_return[1] = vec1[0] * vec2[2] - vec1[2] * vec2[0];
        will_return[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    }
    return std::move(will_return);
}

/**
 * Equal must use to compare two un_empty matrix,
 * so, if one of the matrix is empty , it will be false.
 * only un_empty, rows and columns both equal can be equal.
 */
template<typename T1, typename T2>
bool size_equal(const Matrix<T1> &mat1, const Matrix<T2> &mat2) {
    return !mat1.is_empty() && !mat2.is_empty()  \
 && mat1.rows() == mat2.rows() \
 && mat1.cols() == mat2.cols();
}

/** judge is equal
 * @param1: Matrix_m_n T
 * @param2: Matrix_m_n T
 * @return: bool
 * */
template<typename T>
bool Matrix<T>::inside_equal(const Matrix<T> &mat1, const Matrix<T> &mat2) {
    if (!size_equal(mat1, mat2)) {
        return false;
    }
    for (int32_t i = 0; i < mat1.rows(); ++i) {
        if (!std::equal(mat1.vec[i].begin(), mat1.vec[i].end(), mat2.vec[i].begin(), mat2.vec[i].end())) {
            return false;
        }
    }
    return true;
}

template<typename T>
T Matrix<T>::trace() const {
    T will_return(0);
    if (this->is_square()) {
        for (int32_t i = 0; i < this->rows(); ++i) {
            will_return += this->vec[i][i];
        }
    }
    return will_return;
}

template<typename T>
T determinant_in(const vector<vector<T>> &matrix) {
    if (matrix.size() == 1) {
        return matrix.front().front();
    }
    uint32_t size_m = matrix.size();
    vector<vector<T >> submatrix(size_m - 1, vector<T>(size_m - 1, static_cast<T>(0)));
    T will_return(0);
    for (uint32_t i = 0; i < size_m; ++i) {
        for (uint32_t j = 0; j < size_m - 1; ++j) {
            for (uint32_t k = 0; k < size_m - 1; ++k) {
                submatrix[j][k] = matrix[(((i > j) ? 0 : 1) + j)][k + 1];
            }
        }
        will_return += ((i % 2) ? -1 : 1) * matrix[i].front() * determinant_in(submatrix);
    }
    return will_return;
}

template<typename T>
T Matrix<T>::determinant() const {
    T will_return(0);
    if (!this->is_square()) {
        // TODO
        return will_return;
    }
    return determinant_in(this->vec);
}

template<typename T>
T Matrix<T>::max() const {
    if (!this->is_empty()) {
        T will_return = this->vec.front().front();
        for (const auto &iter:this->vec) {
            auto max_v = *std::max_element(std::begin(iter), std::end(iter));
            will_return = std::max(max_v, will_return);
        }
        return will_return;
    }
    // TODO False;
    return static_cast<T>(0);
}

template<typename T>
T Matrix<T>::min() const {
    if (!this->is_empty()) {
        T will_return = this->vec.front().front();
        for (const auto &iter:this->vec) {
            auto min_v = *std::min_element(std::begin(iter), std::end(iter));
            will_return = std::min(min_v, will_return);
        }
        return will_return;
    }
    // TODO False;
    return static_cast<T>(0);
}

template<typename T>
T Matrix<T>::sum() const {
    T will_return = static_cast<T>(0);
    for (const auto &item : this->vec) {
        will_return += std::accumulate(std::begin(item), std::end(item), static_cast<T>(0));
    }
    return will_return;
}

template<typename T>
T Matrix<T>::row_max(int32_t row) const {
    if (row <= 0 || row > this->rows()) {
        // TODO;
        return -1;
    }
    return (*std::max_element(this->get_row_iter_begin(row - 1),
                              this->get_row_iter_end(row - 1)));
}

template<typename T>
T Matrix<T>::row_min(int32_t row) const {
    if (row <= 0 || row > this->rows()) {
        // TODO;
        return -1;
    }
    return (*std::min_element(this->get_row_iter_begin(row - 1),
                              this->get_row_iter_end(row - 1)));
}

template<typename T>
T Matrix<T>::row_sum(int32_t row) const {
    if (row <= 0 || row > vec.size()) {
        return -1;
    }
    return std::accumulate(this->get_row_iter_begin(row - 1), this->get_row_iter_end(row - 1), static_cast<T>(0));
}

template<typename T>
T Matrix<T>::col_max(int32_t col) const {
    if (col <= 0 || col > vec[0].size()) {
        return -1;
    }
    T max_v = vec.front()[col - 1];
    for (int i = 0; i < this->rows(); ++i) {
        max_v = std::max(max_v, vec[i][col - 1]);
    }
    return max_v;
}

template<typename T>
T Matrix<T>::col_min(int32_t col) const {
    if (col <= 0 || col > vec[0].size()) {
        return -1;
    }
    T min_v = vec.front()[col - 1];
    for (int i = 0; i < this->rows(); ++i) {
        min_v = std::min(min_v, vec[i][col - 1]);
    }
    return min_v;
}

template<typename T>
T Matrix<T>::col_sum(int32_t col) const {
    if (col <= 0 || col > vec[0].size()) {
        return -1;
    }
    T sum(0);
    for (int32_t i = 0; i < this->rows(); ++i) {
        sum += vec[i][col - 1];
    }
    return sum;
}

template<typename T>
Matrix<double_t > Matrix<T>::Householder(int32_t col, int32_t ele) const {
    double_t square = 0;
    for (int i = ele - 1; i < vec.size(); ++i) {
        square += pow(vec[i][col - 1],2);
    }
    double_t mod = vec[ele - 1][col - 1] > 0 ? pow(square,0.5) : -pow(square,0.5);

    double_t modulus = mod*(mod + vec[ele - 1][col - 1]);
    Matrix<double_t > U (vec.size(),1);
    for (int j = 0; j < vec.size(); ++j) {
        if(j > ele - 1){
            U.set_inside(j,0,vec[j][col - 1]);
        }
        else{
            U.set_inside(j,0,0);
        }
    }
    U.set_inside(ele - 1,0,vec[ele - 1][col - 1] + mod);

    Matrix<double_t > H = Matrix<double_t >::eye(vec.size()) - U*U.transpose()/modulus;

    return H;
}

template<typename T>
Matrix<double_t> Matrix<T>::Hessenberg() const{
    if(vec.size() != vec[0].size()){
        return Matrix<double_t >::eye_value(vec.size(),0);
    }

    Matrix<double_t > left_H = Matrix<double_t >::eye(vec.size());
    Matrix<double_t > right_H = Matrix<double_t >::eye(vec.size());
    vector<vector<double_t> > v = this->transform();

    Matrix<double_t > H = left_H * Matrix<double_t >(v) * right_H;

    for (int i = 1; i < vec.size() - 1; ++i) {
        if (H.get_inside(i + 1, i - 1) == 0) {
            continue;
        }
        left_H = left_H * H.Householder(i, i + 1);
        right_H = H.Householder(i, i + 1) * right_H;
        H = left_H * Matrix<double_t>(v) * right_H;
    }
    return H;
}

template<typename T>
void Matrix<T>::set_inside(int32_t row, int32_t col, T ele) {
    vec[row][col] = ele;
}

template<typename T>
void Matrix<T>::set_element(T ele) {
    for (int i = 0; i < vec.size(); ++i) {
        for (int j = 0; j < vec[0].size(); ++j) {
            vec[i][j] = ele;
        }
    }
}

template<typename T>
Matrix<double_t> Matrix<T>::Givens(int32_t col, int32_t begin, int32_t end) const {
    Matrix<double_t > R = Matrix<double_t >::eye(vec.size());
    double_t r = pow(pow(vec[begin - 1][col - 1],2) + pow(vec[end - 1][col - 1],2),0.5);
    double_t c = 1;
    double_t s = 0;
    if(r != 0){
        c = vec[begin - 1][col - 1]/r;
        s = vec[end - 1][col - 1]/r;
    }

    R.set_inside(begin - 1, begin - 1, c);
    R.set_inside(begin - 1, end - 1, s);
    R.set_inside(end - 1, begin - 1, -s);
    R.set_inside(end - 1, end - 1, c);
    return R;
}

template<typename T>
Matrix<double_t> Matrix<T>::QR_iteration() const {
    Matrix<double_t > R = this->Hessenberg();
    Matrix<double_t > Q = Matrix<double_t >::eye(vec.size());
    for (int i = 1; i < vec.size(); ++i) {
        Matrix<double_t > temp_R = R.Givens(i,i,i + 1);
        R = temp_R * R;
        Q = Q * temp_R.transpose();
    }
    Matrix<double_t > H = R * Q;
    return H;
}

template<typename T>
double_t* Matrix<T>::eigenvalue() const{
    if(vec.size() != vec[0].size()){
        return nullptr;
    }
    auto * eigenvalues = new double_t[vec.size()];
    int32_t iter_times = 150;
    Matrix<double_t > H = this->QR_iteration();
    for (int i = 0; i < iter_times; ++i) {
        H = H.QR_iteration();
    }
    for (int j = 0; j < vec.size(); ++j) {
        eigenvalues[j] = H.get_inside(j,j);
    }
    return eigenvalues;
}

template<typename T>
Matrix<double_t> Matrix<T>::eigenvector()const {
    double_t* eigenvalues = this->eigenvalue();
    vector<vector<double_t> > v = this->transform();
    std::sort(eigenvalues,eigenvalues + vec.size());

    int32_t begin = 0;
    int32_t end = 1;
    int32_t repeat = 0;
    while(end < vec.size()){
        if(round(eigenvalues[begin]*10)/10 == round(eigenvalues[end]*10)/10){
            end++;
            repeat++;
        }
        else{
            eigenvalues[begin + 1] = eigenvalues[end];
            begin++;
            end++;
        }
    }

    Matrix<double_t > vectors(vec.size(),vec.size());
    int32_t size = vec.size() - repeat;
    for (int i = 0; i < size; ++i) {
        Matrix<double_t > eigenmatrix(v);
        double_t eigenvalue = eigenvalues[i];
        for (int j = 0; j < vec.size(); ++j) {
            eigenmatrix.set_inside(j, j, round((eigenmatrix.get_inside(j, j) - eigenvalue)*10)/10);
        }
        if(i == 0){
            vectors = eigenmatrix.Nullspace();
        }
        else{
            vectors = vectors.joint(eigenmatrix.Nullspace());
        }
    }
    delete [] eigenvalues;
    return vectors;
}

template<typename T>
Matrix<T> Matrix<T>::joint(Matrix<T> matrix) const {
    if(vec.size() != matrix.rows() && vec.size() != 0){
        return *this;
    }
    int32_t row_size = matrix.rows();
    int32_t col_size = vec[0].size() + matrix.cols();
    Matrix<T> joint(row_size,col_size);

    for (int i = 0; i < row_size; ++i) {
        for (int j = 0; j < col_size; ++j) {
            if(j < vec[0].size()){
                joint.set_inside(i,j,vec[i][j]);
            }
            else{
                joint.set_inside(i,j,matrix.get_inside(i,j - vec[0].size()));
            }
        }
    }
    return joint;
}

template<typename T>
Matrix<double_t> Matrix<T>::Gauss() const {
    vector<vector<double_t> > v = this->transform();
    Matrix<double_t > REF(v);
    int32_t row = 0;
    int32_t find = 0;

    while(row < vec.size()){
        double_t element = 1;
        bool flag = false;
        for (int i = find; i < vec.size(); ++i) {
            if(REF.get_inside(i,row) != 0){
                element = REF.get_inside(i,row);
                REF = REF.row_exchange(i,find);
                flag = true;
                break;
            }
        }
        if(flag){
            REF = REF.row_elimination(find,element);
            for (int i = find + 1; i < vec.size(); ++i) {
                REF = REF.row_elimination(i,row,find);
            }
            find++;
        }
        row++;
    }
    return REF;
}
/*
template<typename T>
Matrix<double_t> Matrix<T>::LU_decomposition() const{
    Matrix<double_t > L = Matrix<double_t >::eye(vec.size());
    Matrix<double_t > U = Matrix<double_t >::eye_value(vec.size(),0);
    if(vec.size() != vec[0].size()){
        return U;
    }

    if(vec[0][0] != 0){
        for (int i = 0; i < vec.size(); ++i) {
            L.set_inside(i,0,(double_t)vec[i][0]/vec[0][0]);
        }
    }
    for (int j = 0; j < vec[0].size(); ++j) {
        U.set_inside(0,j,vec[0][j]);
    }

    for (int k = 1; k < vec.size(); ++k) {
        for (int j = k; j < vec.size(); ++j) {
            double_t temp = 0;
            for (int t = 0; t < k; ++t) {
                temp += L.get_inside(k,t) * U.get_inside(t,j);
            }
            U.set_inside(k,j,vec[k][j] - temp);
        }
        for (int i = k + 1; i < vec.size(); ++i) {
            double_t temp = 0;
            for (int t = 0; t < k; ++t) {
                temp += L.get_inside(i,t) * U.get_inside(t,k);
            }
            if(U.get_inside(k,k) != 0){
            L.set_inside(i,k,(vec[i][k]-temp)/U.get_inside(k,k));
            }
        }
    }
    return U;
}*/

template<typename T>
Matrix<double_t > Matrix<T>::Elimination() const {
    Matrix<double_t > R = this->Gauss();
    for (int i = 0; i < vec.size(); ++i) {
        for (int j = 0; j < vec[0].size(); ++j) {
            if(R.get_inside(i,j) != 0){
                for (int k = i - 1; k >= 0; --k) {
                    R = R.row_elimination(k,j,i);
                }
                break;
            }
        }
    }
    return R;
}

template<typename T>
Matrix<double_t> Matrix<T>::Nullspace()const {
    Matrix<double_t > matrix = this->Elimination();
    int32_t pivots[vec[0].size()];
    memset(pivots, 0, sizeof(pivots));
    int32_t null_num = vec[0].size();

    for (int i = 0; i < vec.size(); ++i) {
        for (int j = 0; j < vec[0].size(); ++j) {
            if(matrix.get_inside(i,j) != 0){
                pivots[j] = 1;
                null_num--;
                break;
            }
        }
    }

    Matrix<double_t > vectors(vec[0].size(),null_num);
    vectors.set_element(0);
    int32_t col = 0;
    for (int i = 0; i < vec[0].size(); ++i) {
        if(pivots[i] == 0){
            for (int j = 0; j < null_num; ++j) {
                if(j == col){
                    vectors.set_inside(i,j,1);
                }
            }
            col++;
        }
        else{
            for (int k = 0; k < vec[0].size(); ++k) {
                if(pivots[k] == 0){
                    for (int j = 0; j < null_num; ++j) {
                        vectors.set_inside(i,j,0 - matrix.get_inside(i,k));
                    }
                }
            }
        }
    }
    return vectors;
}

template<typename T>
vector<vector<double_t> > Matrix<T>::transform() const{
    vector<vector<double_t >> v(vec.size());
    for (int i = 0; i < vec.size(); ++i) {
        v[i] = vector<double_t >(vec[0].size());
        for (int j = 0; j < vec[0].size(); ++j) {
            v[i][j] = vec[i][j];
        }
    }
    return v;
}

template<typename T>
Matrix<T> Matrix<T>::row_exchange(int32_t row1, int32_t row2)const {
    Matrix<T> matrix(vec);
    for (int i = 0; i < vec[0].size(); ++i) {
        matrix.set_inside(row1,i,vec[row2][i]);
        matrix.set_inside(row2,i,vec[row1][i]);
    }
    return matrix;
}

template<typename T>
Matrix<double_t> Matrix<T>::row_elimination(int32_t row, double_t ele) const{
    Matrix<double_t > matrix(this->transform());
    for (int i = 0; i < vec[0].size(); ++i) {
        if(ele <= pow(10,-5) && ele >= 0 - pow(10,-5)){
            matrix.set_inside(row,i,0);
        }
        else{
            matrix.set_inside(row,i,(double_t)this->get_inside(row,i)/ele);
        }
    }
    return matrix;
}

template<typename T>
Matrix<double_t> Matrix<T>::row_elimination(int32_t row, int32_t col, int32_t remove_row)const {
    Matrix<double_t > matrix(this->transform());
    double_t temp = vec[row][col];
    for (int i = 0; i < vec[0].size(); ++i) {
        matrix.set_inside(row,i,this->get_inside(row,i) - temp*vec[remove_row][i]);
    }
    return matrix;
}


#endif //CS205_C_CPP_CS205_PROJECT_2020S_SRC_MATRIX_HPP
