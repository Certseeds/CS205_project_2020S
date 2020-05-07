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

    explicit Matrix(int64_t rows, int64_t col);

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

    friend std::ostream &operator<<(std::ostream &output, const Matrix<T> &mat) {
        for (const auto &i: mat.vec) {
            for (const auto &j:i) {
                output << j << " ";
            }
            output << endl;
        }
        return output;
    }

    ~Matrix() = default;

    bool isEmpty();
};

template<class T>
Matrix<T>::Matrix(int64_t rows, int64_t cols) {
    rows = rows > 0 ? rows : 0;
    cols = cols > 0 ? cols : 0;
    this->vec = vector<vector<T>>(rows, vector<T>(cols));
}


template<class T>
Matrix<T>::Matrix(const Matrix &mat) {
    this->vec = vector<vector<T>>(mat.vec);
}

template<class T>
Matrix<T>::Matrix(Matrix &&mat) noexcept {
    this->vec = vector<vector<T>>(std::move(mat.vec));
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
bool Matrix<T>::isEmpty() {
    return (vec.empty() || vec[0].empty());
}


#endif //CS205_C_CPP_CS205_PROJECT_2020S_SRC_MATRIX_HPP
