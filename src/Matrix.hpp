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

#include <algorithm>
#include <complex>
#include <cstdint>
#include <functional>
#include <numeric>
#include <utility>
#include <vector>
#include <valarray>

#include "./template_helper.h"

#include <opencv2/opencv.hpp>

namespace Mat_pro {
    using std::endl;
    using std::vector;
    using std::valarray;
    using cv::Mat;

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
        Matrix<T> &operator=(Matrix<T> &&mat) noexcept;

        explicit Matrix<T>(const vector<vector<T>> &vec);

        explicit Matrix<T>(vector<vector<T>> &&vec);

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
            for (const auto &i : mat.vec) {
                for (const auto &j : i) {
                    output << j << " ";
                }
                output << endl;
            }
            return output;
        }


        [[nodiscard]] inline int32_t rows() const;

        [[nodiscard]] inline int32_t cols() const;

        Matrix<T> transpose() const;

        ~Matrix() = default;

        [[nodiscard]] inline bool is_empty() const;

        [[nodiscard]] inline bool is_square() const { return this->rows() == this->cols(); };

        inline static Matrix<T>
        operator_table(const Matrix<T> &mat1, const Matrix<T> &mat2,
                       const std::function<T(const T &t1, const T &t2)> &func);

        inline static Matrix<T>
        operator_table(const Matrix<T> &mat1, const std::function<T(const T &t1)> &func);

        typename vector<T>::const_iterator get_row_iter_begin(int32_t rows) const;

        typename vector<T>::const_iterator get_row_iter_end(int32_t rows) const;

        /**@related: https://stackoverflow.com/questions/30736951/templated-class-check-if-complex-at-compile-time
         * */
        Matrix<T> conj() const;

        T max() const;

        T min() const;

        T sum() const;

        template<typename T1 = T, MY_IF(!is_complex<T1>())>
        auto avg() -> Divide_Result_t<T1, double_t> const {
            T sums = this->sum();
            return sums / static_cast<double_t>(this->rows() * this->cols());
        }

        //template< class T1 = T ,typename std::enable_if<is_complex<T1>(),int>::type = 0>
        template<typename T1 = T, MY_IF(is_complex<T1>())>
        T1 avg() const {
            T1 sums = this->sum();
            double_t size = static_cast<double_t>(this->rows() * this->cols());
            return T1{sums.real() / size, sums.imag() / size};
        }

        T row_max(int32_t row) const;

        T col_max(int32_t col) const;

        T row_min(int32_t row) const;

        T col_min(int32_t col) const;

        T row_sum(int32_t row) const;

        T col_sum(int32_t col) const;

        template<typename T1 = T, MY_IF(!is_complex<T1>())>
        auto row_avg(int32_t row) -> Divide_Result_t<T1, double_t> const {
            if (row > this->rows()) {
                return -1;
            }
            return this->row_sum(row) / static_cast<double_t>(this->cols());
            //        return 1;
        }

        template<typename T1 = T, MY_IF(is_complex<T1>())>
        T1 row_avg(int32_t row) const {
            if (row > this->rows() || this->is_empty()) {
                return static_cast<T1>(-1);
            }
            T1 sums = this->row_sum(row);
            double_t col = static_cast<double_t>(this->cols());
            return T1{sums.real() / col, sums.imag() / col};
            // return std::complex(static_cast<double>(2), static_cast<double >(3));
        }

        template<typename T1 = T, MY_IF(!is_complex<T1>())>
        auto col_avg(int32_t col) -> Divide_Result_t<T1, double_t> const {
            if (col <= 0 || col > this->cols()) {
                return -1;
            }
            return this->col_sum(col) / static_cast<double>(this->rows());
        }

        template<typename T1 = T, MY_IF(is_complex<T1>())>
        T1 col_avg(int32_t col) const {
            if (col <= 0 || col > this->cols()) {
                return -1;
            }
            T1 sums = this->col_sum(col);
            double_t size = static_cast<double_t>(this->rows());
            return T1{sums.real() / size, sums.imag() / size};
        }

        // Matrix_n_m, Matrix_n_m, result is Matrix_N_M
        Matrix<T> mul(const Matrix<T> &mat2);

        T dot(const Matrix<T> &mat2);

        T trace() const;

        T determinant() const;

        Matrix<T> convolution(const Matrix<T> &kernel, int32_t padding = 0, int32_t stride = 1) const;

        Matrix<T> reshape(int32_t row, int32_t col) const;  //重整
        Matrix<T> slice(int32_t row1, int32_t row2) const;  //切片
        Matrix<T> slice(int32_t row1, int32_t row2, int32_t col1, int32_t col2) const;  //切片
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
        for (auto i = list.begin(); i != list.end() - 1; i++) {
            if ((*i).size() != (*(i + 1)).size()) {
                return;
            }
        }
        this->vec = vector<vector<T>>(list.size());
        for (auto i = list.begin(); i != list.end(); i++) {
            this->vec[i - list.begin()] = *i;
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
        for (auto &i : will_return.vec) {
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
    auto operator*(const Matrix<T1> &mat1, const T2 &t2) {
        vector<vector<Multiply_Result_t<T1, T2>>> temp(mat1.rows(), vector<Multiply_Result_t<T1, T2>>(mat1.cols()));
        for (uint32_t i = 0; i < temp.size(); ++i) {
            for (uint32_t j = 0; j < temp[i].size(); ++j) {
                temp[i][j] = mat1.get_inside(i, j) * t2;
            }
        }
        return Matrix<Multiply_Result_t<T1, T2>>(std::move(temp));
    }

/**
 * matrix * vector
 * @param1: Matrix<T1> m_n
 * @Param2; vector<T2> length is n
 * @return: Matrix<decltype(T1*T2)> m_1,(rows is m,cols is 1)
 * */
    template<typename T1, typename T2>
    auto operator*(const Matrix<T1> &mat1, const vector<T2> &t2) -> Matrix<Multiply_Result_t<T1, T2>> {
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
    auto operator*(const vector<T1> &t1, const Matrix<T2> &mat2) -> Matrix<Multiply_Result_t<T1, T2>> {
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
    // used to  -> Matrix<Multiply_Result_t_Macro> ;
    template<typename T1, typename T2>
    auto operator*(const T1 &t1, const Matrix<T2> &mat2) {
        return mat2 * t1;
    }

    template<typename T>
    Matrix<T> operator*(const Matrix<T> &mat1, const Matrix<T> &mat2) {
        if (mat1.cols() != mat2.rows()) {
            // TODO
        }
        Matrix<T> temp = mat2.transpose();
        vector<vector<T>> will_return(mat1.rows(), vector<T>(mat2.cols()));
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
    auto kron(const Matrix<T1> &mat1, const Matrix<T2> &mat2) -> Matrix<Multiply_Result_t<T1, T2>> {
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
    auto cross(const vector<T1> &vec1, const vector<T2> &vec2) -> vector<Multiply_Result_t<T1, T2>> {
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
        return !mat1.is_empty() && !mat2.is_empty() && mat1.rows() == mat2.rows() && mat1.cols() == mat2.cols();
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
        vector<vector<T>> submatrix(size_m - 1, vector<T>(size_m - 1, static_cast<T>(0)));
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
    Matrix<T> Matrix<T>::conj() const {
        vector<vector<T>> vvc(this->rows(), vector<T>(this->cols()));
        if constexpr (is_complex<T>()) {
            for (int32_t i = 0; i < this->rows(); ++i) {
                for (int32_t j = 0; j < this->cols(); ++j) {
                    vvc[i][j] = std::conj(this->vec[i][j]);
                }
            }
        } else if constexpr (has_conj<T>()) {
            for (int32_t i = 0; i < this->rows(); ++i) {
                for (int32_t j = 0; j < this->cols(); ++j) {
                    vvc[i][j] = this->vec[i][j].conj();
                }
            }
        } else if constexpr (!is_complex<T>() && !has_conj<T>()) {
            vvc = this->vec;
        }
        return Matrix<T>(std::move(vvc));
    }

    template<typename T>
    T Matrix<T>::max() const {
        if (!this->is_empty()) {
            T will_return = this->vec.front().front();
            for (const auto &iter : this->vec) {
                will_return = std::max(*std::max_element(std::begin(iter), std::end(iter)), will_return);
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
            for (const auto &iter : this->vec) {
                will_return = std::min(*std::min_element(std::begin(iter), std::end(iter)), will_return);
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
        return *std::max_element(this->get_row_iter_begin(row - 1),
                                 this->get_row_iter_end(row - 1));
    }

    template<typename T>
    T Matrix<T>::row_min(int32_t row) const {
        if (row <= 0 || row > this->rows()) {
            // TODO;
            return -1;
        }
        return *std::min_element(this->get_row_iter_begin(row - 1),
                                 this->get_row_iter_end(row - 1));
    }

    template<typename T>
    T Matrix<T>::row_sum(int32_t row) const {
        if (row <= 0 || row > this->rows()) {
            return -1;
        }
        return std::accumulate(this->get_row_iter_begin(row - 1), this->get_row_iter_end(row - 1), static_cast<T>(0));
    }

    template<typename T>
    T Matrix<T>::col_max(int32_t col) const {
        if (col <= 0 || col > this->cols()) {
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
        if (col <= 0 || col > this->cols()) {
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
        if (col <= 0 || col > this->cols()) {
            return -1;
        }
        T sum(0);
        for (int32_t i = 0; i < this->rows(); ++i) {
            sum += vec[i][col - 1];
        }
        return sum;
    }

/**get convolution of Matrix and kernel
 * @param1: this: Matrix: m_n T
 * @param2: mat2: Matrix: f1_f2 T, addvise f1=f2 and they are odd.
 * @param3: padding: int32_t the size of padding.<br> default = 0
 * @param4: stride: int32_t the stride of each step.<br> default = 1
 * @return: Matrix: x1_x2 T</br>
 * <br>x1 = ⌊(m+2p-f_1)/s⌋+1</br>
 * <br>x2 = ⌊(n+2p-f_2)/s⌋+1</br>
 * */
    template<typename T>
    Matrix<T> Matrix<T>::convolution(const Matrix<T> &kernel, int32_t padding, int32_t stride) const {
        if (padding <= 0 || stride <= 0 || this->rows() + 2 * padding > kernel.rows() ||
            this->cols() + 2 * padding > kernel.cols()) {
            // TODO
        }
        //padding = std::max(padding, std::max(kernel.rows(), kernel.cols()));
        int32_t new_row = floor((this->rows() + 2 * padding - kernel.rows()) / stride) + 1;
        int32_t new_col = floor((this->cols() + 2 * padding - kernel.cols()) / stride) + 1;
        vector<vector<T>> will_return(new_row, vector<T>(new_col, static_cast<T>(0)));
        vector<vector<T>> big_vec(this->rows() + padding * 2, vector<T>(this->cols() + padding * 2, static_cast<T>(0)));
        for (int i = padding; i < this->rows() + padding; ++i) {
            for (int j = padding; j < this->cols() + padding; ++j) {
                big_vec[i][j] = this->vec[i - padding][j - padding];
            }
        }
        for (int k = 0; k < new_row; k++) {
            for (int i = 0; i < new_col; i++) {
                for (int j = 0; j < kernel.rows(); ++j) {
                    will_return[k][i] += std::inner_product(std::begin(kernel.vec[j]), std::end(kernel.vec[j]),
                                                            std::begin(big_vec[k * stride + j]) + i * stride,
                                                            static_cast<T>(0));
                }
            }
        }
        return Matrix<T>(std::move(will_return));
    }

    template<typename T>
    Matrix<T> Matrix<T>::reshape(int32_t row, int32_t col) const {
        int32_t col_num = this->cols();
        int32_t num = this->rows() * col_num;
        if (row * col != num || num <= 0) {
            // TODO, there should be error
            return Matrix<T>(*this);
        }
        vector<vector<T>> res(row, vector<T>(col, static_cast<T>(0)));
        for (int i = 0; i < num; i++) {
            res[i / col][i % col] = this->get_inside(i / this->cols(), i % col_num);
        }
        return Matrix<T>(std::move(res));
    }

    template<typename T>
    Matrix<T> Matrix<T>::slice(int32_t row1, int32_t row2) const {
        if (row1 < 0 || row2 >= this->rows()) {
            return Matrix<T>(*this);
        }
        vector<vector<T>> will_return(row2 - row1 + 1, vector<T>(this->cols()));
        for (int32_t i = row1; i < row2; i++) {
            will_return[i - row1] = this->vec[i];
        }
        return Matrix<T>(std::move(will_return));
    }

    template<typename T>
    Matrix<T> Matrix<T>::slice(int32_t row1, int32_t row2, int32_t col1, int32_t col2) const {
        if (row1 < 0 || row2 >= this->rows()
            || col1 < 0 || col2 >= this->cols()
            || row1 > row2 || col1 > col2) {
            return Matrix<T>(*this);
        }
        vector<vector<T>> will_return(row2 - row1 + 1, vector<T>(col2 - col1 + 1));
        for (int32_t i = row1; i < row2; i++) {
            for (int32_t j = col1; j < col2; j++) {
                will_return[i - row1][j - col1] = vec[i][j];
            }
        }
        return Matrix<T>(std::move(will_return));
    }

    template<typename T>
    Matrix<T> cv_to_mat(const Mat &m) {
        vector<vector<T>> will_return(m.rows, vector<T>(m.cols * m.channels()));
        for (int i = 0; i < m.rows; ++i) {
            auto temp_head = m.ptr(i);
            for (int j = 0; j < m.cols * m.channels(); ++j) {
                switch (m.type() % 8) {
                    case 5: {
                        will_return[i][j] = *(float *) (temp_head + j * m.elemSize1());
                        break;
                    }
                    case 6: {
                        will_return[i][j] = *(double *) (temp_head + j * m.elemSize1());
                        break;
                    }
                    default: {
                        will_return[i][j] = temp_head[j * m.elemSize1()];
                        break;
                    }
                }
            }
        }
        return Matrix<T>(std::move(will_return));
    }

    template<typename T, MY_IF(
            is_same<T, uint8_t>()
            || is_same<T, int8_t>()
            || is_same<T, uint16_t>()
            || is_same<T, int16_t>()
            || is_same<T, int32_t>()
            || is_same<T, float>()
            || is_same<T, double>())>
    Mat mat_to_cv(const Matrix<T> &matrix, int32_t demen) {
        int type = 7;
        if constexpr(is_same<T, uint8_t>()) {
            type = 0;
        } else if constexpr(is_same<T, int8_t>()) {
            type = 1;
        } else if constexpr(is_same<T, uint16_t>()) {
            type = 2;
        } else if constexpr(is_same<T, int16_t>()) {
            type = 3;
        } else if constexpr(is_same<T, int32_t>()) {
            type = 4;
        } else if constexpr(is_same<T, float>()) {
            type = 5;
        } else if constexpr(is_same<T, double>()) {
            type = 6;
        }
        if (matrix.cols() % demen != 0) {
            // TODO
        }

        Mat will_return(matrix.rows(), matrix.cols() / demen, type + (demen - 1) * 8);
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                will_return.at<T>(i, j) = matrix.get_inside(i, j);
            }
        }
        return will_return;
    }
}
#endif  //CS205_C_CPP_CS205_PROJECT_2020S_SRC_MATRIX_HPP
