add_executable(${PROJECT_NAME}_unit_test ${CMAKE_CURRENT_SOURCE_DIR}/test_matrix_1.cpp)
target_link_libraries(${PROJECT_NAME}_unit_test
        PUBLIC
        Catch
        )
IF (OpenCV_FOUND)
    target_link_libraries(${PROJECT_NAME}_unit_test
            PUBLIC
            ${OpenCV_LIBS}
            )
    MESSAGE("TARGET LINK OPENCV_LIBS")
endif ()
add_test(${PROJECT_NAME}_CTest1 ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_unit_test
        COMMAND ${PROJECT_NAME}_unit_test)
add_test(${PROJECT_NAME}_CTest2 ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_unit_test
        COMMAND ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_unit_test
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(${PROJECT_NAME}_CTest3 ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_unit_test
        COMMAND ${PROJECT_NAME}_unit_test
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
# PS, just show how to use cmake
MESSAGE(${CMAKE_BINARY_DIR}) # use this
MESSAGE(${CMAKE_CURRENT_BINARY_DIR}) # not this