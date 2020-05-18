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
#include <functional>
#include <complex>

using std::vector;
using std::endl;

template<class T>
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

    T get_type() const;

    T max() const;

    T min() const;

    T sum() const;

    double_t avg() const;

    std::complex<double_t> complex_avg() const;

    T row_max(int32_t row) const;

    T col_max(int32_t col) const;

    T row_min(int32_t row) const;

    T col_min(int32_t col) const;

    T row_sum(int32_t row) const;

    T col_sum(int32_t col) const;

    double_t row_avg(int32_t row) const;

    std::complex<double_t> complex_row_avg(int32_t row) const;

    std::complex<double_t> complex_col_avg(int32_t col) const;

    double_t col_avg(int32_t col) const;

    // Matrix_n_m, Matrix_n_m, result is Matrix_N_M
    Matrix<T> mul(const Matrix<T> &mat2);
};


template<class T>
Matrix<T>::Matrix(const Matrix &mat) {
    this->vec = vector<vector<T>>(mat.vec);
}

template<class T>
Matrix<T>::Matrix(int32_t rows, int32_t cols) {
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
Matrix<T>::Matrix(const vector<vector<T>> &vec) {
    this->vec = vector<vector<T>>(vec);
}

template<class T>
Matrix<T>::Matrix(vector<vector<T>> &&vec) {
    this->vec = std::move(vec);
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
T Matrix<T>::get_inside(int32_t rows, int32_t cols) const {
    return this->vec[rows][cols];
}

template<class T>
inline bool Matrix<T>::isEmpty() const {
    return vec.empty() || vec.front().empty();
}

template<class T>
Matrix<T> Matrix<T>::zeros(int32_t rows, int32_t cols) {
    return Matrix<T>::values(rows, cols, static_cast<T>(0));
}

template<class T>
Matrix<T> Matrix<T>::ones(int32_t rows, int32_t cols) {
    return Matrix<T>::values(rows, cols, static_cast<T>(1));
}

template<class T>
Matrix<T> Matrix<T>::values(int32_t rows, int32_t cols, T t) {
    Matrix<T> will_return(rows, cols);
    for (auto &i:will_return.vec) {
        i = vector<T>(cols, t);
    }
    return will_return;
}

template<class T>
Matrix<T> Matrix<T>::eye_value(int32_t s, T t) {
    Matrix<T> will_return(s, s);
    for (int32_t i = 0; i < s; ++i) {
        will_return.vec[i][i] = t;
    }
    return will_return;
}

template<class T>
Matrix<T> Matrix<T>::eye(int32_t s) {
    return Matrix<T>::eye_value(s, static_cast<T>(1));
}

template<class T>
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

template<class T>
Matrix<T> Matrix<T>::operator_table(const Matrix<T> &mat1, const std::function<T(const T &t1)> &func) {
    Matrix<T> will_return(mat1.rows(), mat1.cols());
    for (int32_t i = 0; i < mat1.rows(); ++i) {
        std::transform(mat1.vec[i].begin(), mat1.vec[i].end(), will_return.vec[i].begin(), func);
    }
    return will_return;
}

template<class T>
inline int32_t Matrix<T>::rows() const {
    return this->vec.size();
}

/**
 * if the it's empty,then isEmpty is true,
 * */
template<class T>
inline int32_t Matrix<T>::cols() const {
    if (this->rows() == 0) {
        return 0;
    }
    return this->vec.front().size();
}


template<class T>
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


/**
 * matrix + matrix, must equal.
 * */
template<class T>
Matrix<T> operator+(const Matrix<T> &mat1, const Matrix<T> &mat2) {
    return Matrix<T>::operator_table(mat1, mat2, std::plus<>());
    // [](const T &t1, const T &t2) -> T { return t1 + t2; }
}

/**
 * matrix - matrix, must equal.
 * */
template<class T>
Matrix<T> operator-(const Matrix<T> &mat1, const Matrix<T> &mat2) {
    return Matrix<T>::operator_table(mat1, mat2, std::minus<>());
}

/**
 * matrix * number, need T1 can * T2
 * */
template<typename T1, typename T2>
auto operator*(const Matrix<T1> &mat1, const T2 &t2) -> Matrix<decltype(mat1.get_type() * t2)> {
    vector<vector<decltype(mat1.get_type() * t2)>> temp(mat1.rows(),
                                                        vector<decltype(mat1.get_type() * t2)>(mat1.cols()));
    for (uint32_t i = 0; i < temp.size(); ++i) {
        for (uint32_t j = 0; j < temp[i].size(); ++j) {
            temp[i][j] = mat1.get_inside(i, j) * t2;
        }
    }
    //  Matrix<decltype(mat1.get_type() * t2)> will_return(mat1.rows(), mat1.cols());
    return Matrix<decltype(mat1.get_type() * t2)>(std::move(temp));
    //for (int32_t i = 0; i < mat1.rows(); ++i) {
    //    std::transform(mat1.vec[i].begin(), mat1.vec[i].end(), will_return.vec[i].begin(), func);
    //}
    //return Matrix<decltype(mat1.get_type() * t2)>::operator_table(mat1, [&t2](const T1 &t1) { return t1 * t2; });
    // return Matrix<T1>::operator_table(mat1, [&t2](const T1 &t1) { return t1 * t2; });
}

/**
 *  number * matrix , need T1 can * T2
 * */
template<typename T1, typename T2>
auto operator*(const T1 &t1, const Matrix<T2> &mat2) -> Matrix<decltype(mat2.get_type() * t1)> {
    return mat2 * t1;
}

template<class T>
Matrix<T> operator*(const Matrix<T> &mat1, const Matrix<T> &mat2) {
    if (mat1.cols() != mat2.rows()) {
        // TODO
    }
    Matrix<T> temp = mat2.transpose();
    Matrix<T> will_return(mat1.rows(), mat2.cols());
    for (int32_t i = 0; i < will_return.rows(); ++i) {
        for (int32_t j = 0; j < will_return.cols(); ++j) {
            will_return.vec[i][j] = std::inner_product(std::begin(will_return.vec) + i, std::end(will_return.vec),
                                                       std::begin(temp.vec) + j, std::end(temp.vec), 0);
        }
    }
    return will_return;
}

// Matrix_n_m, Matrix_n_m, result is Matrix_N_M
template<class T>
Matrix<T> Matrix<T>::mul(const Matrix<T> &mat2) {
    return Matrix<T>::operator_table(*this, mat2, std::multiplies<>());
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

template<typename T>
bool Matrix<T>::inside_equal(const Matrix<T> &mat1, const Matrix<T> &mat2) {
    if (!size_equal(mat1, mat2)) {
        return false;
    }
    for (int32_t i = 0; i < mat1.rows(); ++i) {
        if (!std::equal(mat1.vec[i].begin(), mat1.vec[i].end(),
                        mat2.vec[i].begin(), mat2.vec[i].end())) {
            return false;
        }
    }
    return true;
}

template<class T>
T Matrix<T>::get_type() const {
    if (this->isEmpty()) {
        return static_cast<T>(0);
    }
    return this->vec.front().front();
}

template<class T>
T Matrix<T>::max() const {
    T Max = vec[0][0];
    for (int i = 0; i < vec.size(); ++i) {
        auto max = *std::max_element(std::begin(vec[i]),std::end(vec[i]));
        if(max > Max)
            Max = max;
    }
    return Max;
}

template<class T>
T Matrix<T>::min() const {
    T Min = vec[0][0];
    for (int i = 0; i < vec.size(); ++i) {
        auto min = *std::min_element(std::begin(vec[i]),std::end(vec[i]));
        if(min < Min)
            Min = min;
    }
    return Min;
}

template<class T>
T Matrix<T>::sum() const {
    T sum = 0;
    for (int i = 0; i < vec.size(); ++i) {
        for (int j = 0; j < vec[i].size(); ++j) {
            sum += vec[i][j];
        }
    }
    return sum;
}

template<class T>
double_t Matrix<T>::avg() const {
    int size = vec.size()*vec[0].size();
    double_t avg = (double_t)this->sum()/size;
    return avg;
}

template<class T>
std::complex<double_t> Matrix<T>::complex_avg() const {
    int size = vec.size()*vec[0].size();
    double_t img = (double_t)this->sum().imag()/size;
    double_t real = (double_t)this->sum().real()/size;
    std::complex<double_t> avg;
    avg.imag(img);
    avg.real(real);
    return avg;
}

template<class T>
T Matrix<T>::row_max(int32_t row) const {
    if(row <= 0 || row > vec.size()){
        return -1;
    }
    T max = vec[row-1][0];
    for (int i = 0; i < vec[row-1].size(); ++i) {
        if(max < vec[row-1][i])
            max = vec[row-1][i];
    }
    return max;
}

template<class T>
T Matrix<T>::row_min(int32_t row) const {
    if(row <= 0 || row > vec.size()){
        return -1;
    }
    T min = vec[row-1][0];
    for (int i = 0; i < vec[row-1].size(); ++i) {
        if(min > vec[row-1][i])
            min = vec[row-1][i];
    }
    return min;
}

template<class T>
T Matrix<T>::row_sum(int32_t row) const {
    if(row <= 0 || row > vec.size()){
        return -1;
    }
    T sum = 0;
    for (int i = 0; i < vec[row-1].size(); ++i) {
        sum += vec[row-1][i];
    }
    return sum;
}

template<class T>
double_t Matrix<T>::row_avg(int32_t row) const {
    if(row <= 0 || row > vec.size()){
        return -1;
    }
    double_t avg = (double_t)this->row_sum(row)/vec[row-1].size();
    return avg;
}

template<class T>
std::complex<double_t> Matrix<T>::complex_row_avg(int32_t row) const {
    if(row <= 0 || row > vec.size()){
        return -1;
    }
    int size = vec[row-1].size();
    double_t img = (double_t)this->row_sum(row).imag()/size;
    double_t real = (double_t)this->row_sum(row).real()/size;
    std::complex<double_t> avg;
    avg.imag(img);
    avg.real(real);
    return avg;
}

template<class T>
T Matrix<T>::col_max(int32_t col) const {
    if(col <= 0 || col > vec[0].size()){
        return -1;
    }
    T max = vec[0][col-1];
    for (int i = 0; i < vec.size(); ++i) {
        if(max < vec[i][col-1])
            max = vec[i][col-1];
    }
    return max;
}

template<class T>
T Matrix<T>::col_min(int32_t col) const {
    if(col <= 0 || col > vec[0].size()){
        return -1;
    }
    T min = vec[0][col-1];
    for (int i = 0; i < vec.size(); ++i) {
        if(min > vec[i][col-1])
            min = vec[i][col-1];
    }
    return min;
}

template<class T>
T Matrix<T>::col_sum(int32_t col) const {
    if(col <= 0 || col > vec[0].size()){
        return -1;
    }
    T sum = 0;
    for (int i = 0; i < vec.size(); ++i) {
        sum += vec[i][col-1];
    }
    return sum;
}

template<class T>
double_t Matrix<T>::col_avg(int32_t col) const {
    if(col <= 0 || col > vec[0].size()){
        return -1;
    }
    double_t avg = (double_t)this->col_sum(col)/vec.size();
    return avg;
}

template<class T>
std::complex<double_t> Matrix<T>::complex_col_avg(int32_t col) const {
    if(col <= 0 || col > vec[0].size()){
        return -1;
    }
    int size = vec.size();
    double_t img = (double_t)this->col_sum(col).imag()/size;
    double_t real = (double_t)this->col_sum(col).real()/size;
    std::complex<double_t> avg;
    avg.imag(img);
    avg.real(real);
    return avg;
}

#endif //CS205_C_CPP_CS205_PROJECT_2020S_SRC_MATRIX_HPP
