/*
 * @Github: https://github.com/Certseeds/CS205_C_CPP
 * @Organization: SUSTech
 * @Author: nanoseeds
 * @Date: 2020-08-21 16:41:03
 * @LastEditors: nanoseeds
 * @LastEditTime: 2020-08-21 17:24:01
 */
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
#include <complex>
#include <iostream>

#include "./../../include/catch_main.hpp"
#include "./Matrix.hpp"

using namespace Mat_pro;
using Catch::Matchers::Equals;
using std::cout;
using std::endl;
using std::string;

TEST_CASE("test 0", "[end test]") {
    Mat img = cv::imread("./3c011ba75965bcfc.jpg");
    if (img.empty()) {
        std::cout << "can not load image " << endl;
    }
    cout << img.type() << "\n";
    cout << img.rows << "\n";
    cout << img.cols << "\n";
    auto matrix_pic1 = cv_to_mat<uint16_t>(img);
    cout << matrix_pic1.rows() << "\n";
    cout << matrix_pic1.cols() << "\n";
    auto m2 = Matrix<uint16_t>::values(5, 5, 1);
    Matrix<uint16_t> m3 = {1};
    auto result = matrix_pic1.convolution_mul(m2, 4, 1, 3);
    result = result / static_cast<uint16_t>(25);
    cout << result.rows() << " " << result.cols() << " \n";
    cv::imwrite("store.png", img);
    cv::imwrite("store_origin.jpg", mat_to_cv(matrix_pic1, 3));
    cv::imwrite("store_mask.jpg", mat_to_cv(result, 3));
    cv::waitKey(0);

    Mat green = cv::imread("./green.jpg");
    Mat red2 = cv::imread("./red2.jpg");
    if (green.empty() || red2.empty()) {
        std::cout << "can not load image " << endl;
    }
    auto green_matrix = cv_to_mat<uint16_t>(green);
    auto red_matrix = cv_to_mat<uint16_t>(red2);
    auto result2 = green_matrix + red_matrix - 255;
    cv::imwrite("store_add_result.jpg", mat_to_cv(result2, 3));
}
//
//string getCwd() {
//    //获取当前工作目录
//    string path;
//    path = getcwd(NULL, 0);
//    return path;
//}