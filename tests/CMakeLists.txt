MESSAGE("** SUPPORT INCLUDE => ${CMAKE_BINARY_DIR}")
configure_file("${CMAKE_SOURCE_DIR}/tests/support.hpp.in" "${CMAKE_BINARY_DIR}/tests/support.hpp")