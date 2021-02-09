cmake_minimum_required(VERSION 3.16.6)
project(CS205_project_2020S LANGUAGES CXX)

include(${CMAKE_SOURCE_DIR}/cmake/CMakeHelper.cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

Output_variable(CMAKE_BUILD_TYPE)
#set(OpenCV_LIBS "${OpenCV_DIR}/x64/vc16/lib")

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    include(${CMAKE_SOURCE_DIR}/cmake/windows_OpenCV.cmake)
    MESSAGE("will use opencv")
else ()
endif ()

if (OpenCV_FOUND)
    add_definitions(-D_HAVE_OPENCV_)
endif ()

include_directories(third_party)
add_subdirectory(include)
add_subdirectory(src)
