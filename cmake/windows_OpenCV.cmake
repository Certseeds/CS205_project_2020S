include(${CMAKE_SOURCE_DIR}/cmake/OpenCVPath.cmake)
set(OpenCV_BINS "${OpenCV_DIR}/x64/vc16/bin") # work for windows10-2004

Output_variable(OpenCV_FOUND)
Output_Paths()
find_package(OpenCV)
Output_Paths()
Output_variable(OpenCV_FOUND)
