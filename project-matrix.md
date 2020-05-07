<!--
 * @Github: https://github.com/Certseeds/CS205_C_CPP
 * @Organization: SUSTech
 * @Author: nanoseeds
 * @Date: 2020-05-07 10:09:50
 * @LastEditors: nanoseeds
 * @LastEditTime: 2020-05-07 10:10:52
 * @License: CC-BY-NC-SA_V4_0 or any later version 
 -->
## Building a library for matrix computation

1. Matrix is an important concept introduced in linear algebra. Matrix calculation is widely used in many practical applications, such as image processing and machine learning. Programmers can indeed use many different existing libraries, and in certain cases, programmers are required to design their own matrix calculation libraries for specific implementations. This project will build a new library (do not attempt to directly copy codes from other existing library) that can perform the following operations on the matrix:
    1.  It supports all matrix sizes, from small fixed-size matrices to arbitrarily large dense matrices, and even sparse matrices (Add: try to use efficient ways to store the sparse matrices).
    2.  It supports all standard numeric types, including std::complex, integers, and is easily extensible to custom numeric types.
    3.  It supports matrix and vector arithmetic, including addition, subtraction, scalar multiplication, scalar division, transposition, conjugation, element-wise multiplication, matrix-matrix multiplication, matrix-vector multiplication, dot product and cross product.
    4.  It supports basic arithmetic reduction operations, including finding the maximum value, finding the minimum value, summing all items, calculating the average value (all supporting axis-specific and all items).
    5.  It supports computing eigenvalues and eigenvectors, calculating traces, computing inverse and computing determinant.
    6.  It supports the operations of reshape and slicing.
    7.  It supports convolutional operations of two matrices.
    8.  It supports to transfer the matrix from OpenCV to the matrix of this library and vice versa.
    9.  It should process likely exceptions as much as possible.
