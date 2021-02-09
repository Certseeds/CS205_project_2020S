set(target_name Catch)

add_library(${target_name} INTERFACE)
target_include_directories(${target_name}
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR})
unset(${target_name})
