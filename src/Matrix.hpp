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

using std::vector;
using std::endl;

template<class T>
class Matrix {
private:
    vector<vector<T>> vec;
public:
    //explicit Matrix() { this->vec = vector<vector<T>>(1, vector<T>(1)); };
    explicit Matrix() : Matrix(0, 0) {};

    explicit Matrix(size_t rows, size_t cols);

    // 拷贝构造函数 Copy Constructor
    Matrix(const Matrix<T> &mat);

    // 移动构造函数
    Matrix(Matrix<T> &&mat) noexcept;

    // 拷贝赋值运算符 Copy Assignment operator
    Matrix<T> &operator=(const Matrix<T> &mat);

    // 移动赋值运算符 Move Assignment operator
    Matrix<T> &operator=(Matrix<T> &&mat) noexcept;

    Matrix<T>(const std::initializer_list<std::initializer_list<T>> &list);

    // only for
    Matrix<T>(const std::initializer_list<T> &list);

    static Matrix<T> zeros(size_t rows, size_t cols);

    static Matrix<T> ones(size_t rows, size_t cols);

    static Matrix<T> values(size_t rows = 0, size_t cols = 0, T = static_cast<T>(0));

    static Matrix<T> eye(size_t s);

    static Matrix<T> eye_value(size_t s, T t);

    friend std::ostream &operator<<(std::ostream &output, const Matrix<T> &mat) {
        for (const auto &i: mat.vec) {
            for (const auto &j:i) {
                output << j << " ";
            }
            output << endl;
        }
        return output;
    }

    inline size_t rows() const;

    inline size_t cols() const;

    Matrix<T> operator+(const Matrix<T> &mat1) const;

    Matrix<T> operator-(const Matrix<T> &mat1) const;

    Matrix<T> operator*(T &mat1) const;

    Matrix<T> operator*(const Matrix<T> &mat1) const;

    Matrix<T> mul(const Matrix<T> &mat1) const;

    Matrix<T> operator/(const Matrix<T> &mat1) const;

    inline Matrix<T> operator_table(const Matrix<T> &mat1, const std::function<void()> &) const;

    inline bool Equal(const Matrix<T> &mat1) const;

    ~Matrix() = default;

    inline bool isEmpty();
};


template<class T>
Matrix<T>::Matrix(const Matrix &mat) {
    this->vec = vector<vector<T>>(mat.vec);
}

template<class T>
Matrix<T>::Matrix(size_t rows, size_t cols) {
    rows = rows > 0 ? rows : 0;
    cols = cols > 0 ? cols : 0;
    this->vec = vector<vector<T>>(rows, vector<T>(cols));
}

template<class T>
Matrix<T>::Matrix(Matrix &&mat) noexcept {
    this->vec = std::move(mat.vec);
}

template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &mat) {
    if (this != &mat) {
        this->vec = vector<vector<T>>(mat.vec);
    }
    return *this;
}

template<class T>
Matrix<T> &Matrix<T>::operator=(Matrix &&mat) noexcept {
    if (this != &mat) {
        this->vec = std::move(mat.vec);
    }
    return *this;
}

template<class T>
Matrix<T>::Matrix(const std::initializer_list<T> &list) {
    this->vec = vector<vector<T>>(1);
    this->vec[0] = list;
}

template<class T>
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


template<class T>
inline bool Matrix<T>::isEmpty() {
    return vec.empty() || vec[0].empty();
}

template<class T>
Matrix<T> Matrix<T>::zeros(size_t rows, size_t cols) {
    return Matrix<T>::values(rows, cols, static_cast<T>(0));
}

template<class T>
Matrix<T> Matrix<T>::ones(size_t rows, size_t cols) {
    return Matrix<T>::values(rows, cols, static_cast<T>(1));
}

template<class T>
Matrix<T> Matrix<T>::values(size_t rows, size_t cols, T t) {
    Matrix<T> will_return(rows, cols);
    for (auto &i:will_return.vec) {
        for (auto &j:i) {
            j = t;
        }
    }
    return will_return;
}

template<class T>
Matrix<T> Matrix<T>::eye_value(size_t s, T t) {
    Matrix<T> will_return(s, s);
    for (size_t i = 0; i < s; ++i) {
        will_return.vec[i][i] = t;
    }
    return will_return;
}

template<class T>
Matrix<T> Matrix<T>::eye(size_t s) {
    return Matrix<T>::eye_value(s, static_cast<T>(1));
}

template<class T>
inline Matrix<T> Matrix<T>::operator_table(const Matrix<T> &mat1, const std::function<void()> &func) const {
    if (!this->Equal(mat1)) {
        exit(0);
    }
    Matrix<T> will_return(*this);
    for (size_t i = 0; i < this->vec.size(); ++i) {
        std::transform(will_return.vec.begin() + i, will_return.vec.end(),
                       mat1.begin() + i, mat1.vec.end(),
                       func());
    }
    return will_return;
}

template<class T>
inline size_t Matrix<T>::rows() const {
    return this->vec.size();
}

template<class T>
inline size_t Matrix<T>::cols() const {
    return !this->isEmpty() || this->vec[0].size();
}

// TODO
template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &mat1) const {
    return operator_table(std::plus<T>());
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &mat1) const {
    return operator_table(std::minus<T>());
}

// Matrix_n_m, Matrix_n_m, result is Matrix_N_M
template<class T>
Matrix<T> Matrix<T>::mul(const Matrix<T> &mat1) const {
    return operator_table(std::multiplies<T>());
}

/**
 * Equal must use to compare two un_empty matrix,
 * so, if one of the matrix is empty , it will be false.
 * only rows and columns both equal can be equal.
 */
template<class T>
inline bool Matrix<T>::Equal(const Matrix<T> &mat1) const {
    return !this->isEmpty() && !mat1.isEmpty()  \
 && this->rows() == mat1.vec.rows() \
 && this->cols() == mat1.cols();
}

// TODO, choose equal from above and below.
template<typename T1, typename T2>
bool Equal(const Matrix<T1> &mat1, const Matrix<T2> &mat2) {
    return !mat1.isEmpty() && !mat2.isEmpty()  \
 && mat1.rows() == mat2.rows() \
 && mat1.cols() == mat2.cols();
}

#endif //CS205_C_CPP_CS205_PROJECT_2020S_SRC_MATRIX_HPP
