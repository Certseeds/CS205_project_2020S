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
 * @Date: 2020-05-07 17:48:55 
 * @LastEditors  : nanoseeds
 */
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace std;
using namespace cv;

int main() {
    complex temp(1,2);
    cout << temp.real() << endl;
    cout << temp.imag() << endl;
    Mat mat1 = Mat(3, 4, CV_8U);
    Mat mat2 = Mat(3, 4, CV_8U);
    Mat mat3 = Mat().zeros(3, 4, CV_8U);
    Mat mat4 = Mat().ones(3, 4, CV_8U);
    Mat mat5 = Mat().ones(3, 5, CV_8UC(2));
    cout << mat1 << endl;
    cout << mat2 << endl;
    cout << mat3 << endl;
    cout << mat4 << endl;
    cout << mat4.rows;
    cout << mat4.cols;
    cout << mat4 + mat4 << endl;
    cout << mat4+mat5 <<endl;
    cout << mat5.col(0).row(0) << endl;
    cout << mat5.cols << " " << mat5.rows;
    return 0;
}