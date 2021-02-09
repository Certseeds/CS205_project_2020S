add_executable(${PROJECT_NAME}_unit_test ${CMAKE_CURRENT_SOURCE_DIR}/test_matrix_1.cpp)

target_link_libraries(${PROJECT_NAME}_unit_test
        PUBLIC
        ${OpenCV_LIBS}
        Catch
        )

enable_testing()
add_test(${PROJECT_NAME}_CTest ${PROJECT_NAME}_unit_test)
