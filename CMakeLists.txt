cmake_minimum_required(VERSION 3.10.2)
project(CS205_project_2020S LANGUAGES CXX)

include(${CMAKE_SOURCE_DIR}/cmake/OpenCVPath.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/CMakeHelper.cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

Output_variable(CMAKE_BUILD_TYPE)
#set(OpenCV_LIBS "${OpenCV_DIR}/x64/vc16/lib")
set(OpenCV_BINS "${OpenCV_DIR}/x64/vc16/bin") # work for windows10-2004

Output_variable(OpenCV_FOUND)
Output_Paths()
find_package(OpenCV REQUIRED)
Output_Paths()
Output_variable(OpenCV_FOUND)

if (OpenCV_FOUND)
    add_definitions(-D_HAVE_OPENCV_)
endif ()

include_directories(third_party)
include_directories(include)
add_executable(${PROJECT_NAME}_unit_test ${CMAKE_SOURCE_DIR}/src/test_matrix_1.cpp)
target_link_libraries(${PROJECT_NAME}_unit_test ${OpenCV_LIBS})
enable_testing()
add_test(${PROJECT_NAME}_CTest ${PROJECT_NAME}_unit_test)

#target_link_libraries(${PROJECT_NAME}_unit_test ${OpenCV_LIBS})
# test_matrix_1 for unit test

#add_executable(${PROJECT_NAME}_demo_test ./test_matrix_2.cpp)
# test_matrix_2 for demo.
#target_link_libraries(${PROJECT_NAME}_demo_test ${OpenCV_LIBS})

# add_executable(${PROJECT_NAME}_test_opencv ./opencv_matrix.cpp)
# target_link_libraries(${PROJECT_NAME}_test_opencv ${OpenCV_LIBS})
#include_directories(${OpenCV_INCLUDE_DIRS})