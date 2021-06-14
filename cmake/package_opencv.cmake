include(${CMAKE_SOURCE_DIR}/cmake/OpenCVPath.cmake)
IF (WIN32)
    set(OpenCV_BINS "${OpenCV_DIR}/x64/vc16/bin")
    # work for windows10-2004,vs2019
ELSEIF (APPLE)
    message(FATAL_ERROR "THIS CMAKE SCRIPT HAVEND TEST IN APPLE")
ELSEIF (UNIX)
    message(FATAL_ERROR "This cmake script have not test opencv in unix-like system")
ENDIF ()


