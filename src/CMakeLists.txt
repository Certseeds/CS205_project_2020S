cmake_minimum_required(VERSION 3.10.2)
project(CS205_project_2020S LANGUAGES CXX)

set(OpenCV_INCLUDE_DIRS /usr/local/include/opencv2)
FIND_PACKAGE(OpenCV REQUIRED)
message(STATUS "OpenCV library status:")
message(STATUS "    version: ${OpenCV_VERSION}")
message(STATUS "    libraries: ${OpenCV_LIBS}")
message(STATUS "    include path: ${OpenCV_INCLUDE_DIRS}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include_directories(
        ${PROJECT_SOURCE_DIR}
        /usr/local/include/opencv2
)

add_executable(${PROJECT_NAME}_
        ../../include/catch_main.hpp
        #./test_matrix_2.cpp
        ./test_matrix_1.cpp
        ./Matrix.hpp
        )
# test_matrix_1 for unit test
# test_matrix_2 for demo.
target_link_libraries(${PROJECT_NAME}_ ${OpenCV_LIBS})

add_executable(${PROJECT_NAME}_test_opencv ./opencv_matrix.cpp)
target_link_libraries(${PROJECT_NAME}_test_opencv ${OpenCV_LIBS})
#include_directories(${OpenCV_INCLUDE_DIRS})