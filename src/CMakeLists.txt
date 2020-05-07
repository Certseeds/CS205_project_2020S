cmake_minimum_required(VERSION 3.10.2)
project(CS205_project_2020S LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable(${PROJECT_NAME}
        ./../../catch_main.cpp
        ./test_matrix_1.cpp
        )
