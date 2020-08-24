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
include_directories(./include)
include_directories(./third_party)

add_executable(${PROJECT_NAME}_unit_test
        ./test_matrix_1.cpp)
target_link_libraries(${PROJECT_NAME}_unit_test ${OpenCV_LIBS})
# test_matrix_1 for unit test

add_executable(${PROJECT_NAME}_demo_test
        ./test_matrix_2.cpp)
# test_matrix_2 for demo.
target_link_libraries(${PROJECT_NAME}_demo_test ${OpenCV_LIBS})

add_executable(${PROJECT_NAME}_test_opencv ./opencv_matrix.cpp)
target_link_libraries(${PROJECT_NAME}_test_opencv ${OpenCV_LIBS})
#include_directories(${OpenCV_INCLUDE_DIRS})