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

#define Multiply_Result_t_Macro decltype(std::declval<T1>() * std::declval<T2>())

template<typename T1, typename T2>
struct Multiply_Result {
    using Type = decltype(std::declval<T1>() * std::declval<T2>());
};

template<typename T1, typename T2>
using Multiply_Result_t = typename Multiply_Result<T1, T2>::Type;

using std::vector;
using std::endl;

/** @related: get from https://blog.csdn.net/qq_31175231/article/details/77479692
 */
//template<typename T>
//struct has_member_conj {
//private:
//    template<typename U>
//    static auto Check(int) -> decltype(std::declval<U>().conj(), std::true_type());
//
//    template<typename U>
//    static std::false_type Check(...);
//
//public:
//    enum {
//        value = std::is_same<decltype(Check<T>(0)), std::true_type>::value
//    };
//};

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
    Matrix<T> &operator=(Matrix<T> &&mat) noexcept;

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

    inline bool isEmpty() const;

    inline static Matrix<T>
    operator_table(const Matrix<T> &mat1, const Matrix<T> &mat2,
                   const std::function<T(const T &t1, const T &t2)> &func);

    inline static Matrix<T>
    operator_table(const Matrix<T> &mat1, const std::function<T(const T &t1)> &func);

    typename vector<T>::const_iterator get_row_iter_begin(int32_t rows) const;

    typename vector<T>::const_iterator get_row_iter_end(int32_t rows) const;


    // T get_type() const;

    // Matrix_n_m, Matrix_n_m, result is Matrix_N_M
    Matrix<T> mul(const Matrix<T> &mat2);

    T dot(const Matrix<T> &mat2);
};


template<typename T>
Matrix<T>::Matrix(const Matrix &mat) {
    this->vec = vector<vector<T>>(mat.vec);
}

template<typename T>
Matrix<T>::Matrix(int32_t rows, int32_t cols) {
    rows = rows > 0 ? rows : 0;
    cols = cols > 0 ? cols : 0;
    this->vec = vector<vector<T>>(rows, vector<T>(cols));
}

template<typename T>
Matrix<T>::Matrix(Matrix &&mat) noexcept {
    this->vec = std::move(mat.vec);
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &mat) {
    if (this != &mat) {
        this->vec = vector<vector<T>>(mat.vec);
    }
    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(Matrix &&mat) noexcept {
    if (this != &mat) {
        this->vec = std::move(mat.vec);
    }
    return *this;
}

template<typename T>
Matrix<T>::Matrix(const vector<vector<T>> &vec) {
    this->vec = vector<vector<T>>(vec);
}

template<typename T>
Matrix<T>::Matrix(vector<vector<T>> &&vec) {
    this->vec = std::move(vec);
}

template<typename T>
Matrix<T>::Matrix(const std::initializer_list<T> &list) {
    this->vec = vector<vector<T>>(1);
    this->vec[0] = list;
}

template<typename T>
Matrix<T>::Matrix(const std::initializer_list<std::initializer_list<T>> &list) {
    this->vec = vector<vector<T>>(0);
    for (auto i = list.begin(); i < list.end() - 1; i++) {
        if ((*i).size() != (*(i + 1)).size()) {
            return;
        }
    }
    this->vec = vector<vector<T>>(list.size());
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
inline bool Matrix<T>::isEmpty() const {
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
 * if the it's empty,then isEmpty is true,
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
    if (this->isEmpty()) {
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

template<typename T>
Matrix<T> conj(const Matrix<T> &mat_cp) {
    return mat_cp;
}

template<typename T>
Matrix<std::complex<T>> conj(const Matrix<std::complex<T>> &mat_cp) {
    vector<vector<std::complex<T>>> vvc(mat_cp.rows(), vector<std::complex<T>>(mat_cp.cols()));
    for (int i = 0; i < mat_cp.rows(); ++i) {
        for (int j = 0; j < mat_cp.cols(); ++j) {
            vvc[i][j] = std::conj(mat_cp.get_inside(i, j));
        }
    }
    return Matrix<std::complex<T>>(std::move(vvc));
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
    vector<vector<Multiply_Result_t<T1, T2>>> temp(mat1.rows(), vector<Multiply_Result_t<T1, T2>>(mat1.cols()));
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
auto
operator*(const Matrix<T1> &mat1, const vector<T2> &t2) -> Matrix<Multiply_Result_t_Macro> {
    if (mat1.cols() != t2.size()) {
        // TODO
    }
    vector<vector<Multiply_Result_t<T1, T2>>> temp(1, vector<Multiply_Result_t<T1, T2>>(mat1.rows()));
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
    vector<vector<Multiply_Result_t<T1, T2>>> temp(mat2.cols(), vector<Multiply_Result_t<T1, T2>>(1));
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
    vector<vector<T>> will_return(mat1.rows(), vector<T>(mat1.cols()));
    for (int32_t i = 0; i < mat1.rows(); ++i) {
        for (int32_t j = 0; j < mat1.cols(); ++j) {
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
            mat1.rows() * mat2.rows(), vector<Multiply_Result_t<T1, T2>>(mat1.cols() * mat2.cols()));
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
    return !mat1.isEmpty() && !mat2.isEmpty()  \
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

//template<typename T>
//T Matrix<T>::get_type() const {
//    if (this->isEmpty()) {
//        return static_cast<T>(0);
//    }
//    return this->vec.front().front();
//}
//

#endif //CS205_C_CPP_CS205_PROJECT_2020S_SRC_MATRIX_HPP
