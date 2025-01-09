# for google-tests
project(Google_tests)

macro(addToIt name)
    #    if(EXISTS "${TEST_DIR}/${name}.cpp.inp")
    #        configure_file(${name}.cpp.in "${CMAKE_CURRENT_BINARY_DIR}/${name}.cpp")
    #    endif()
    if(EXISTS "${TEST_DIR}/${name}.cpp")
        list(APPEND CPP_FILES "${TEST_DIR}/${name}.cpp")
    elseif(EXISTS "${TEST_DIR}/${name}.h")
        list(APPEND CPP_FILES "${TEST_DIR}/${name}.h")
    endif()

endmacro(addToIt)


set(TEST_DIR "${CMAKE_SOURCE_DIR}/tests")
MESSAGE("** TEST DIR ${TEST_DIR}")
include("${TEST_DIR}/CMakeLists.txt") # this will run configure to make the support.hpp file

MESSAGE("MODULES PATH ${CMAKE_MODULE_PATH}")

add_subdirectory(google-tests)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

MESSAGE("CMAKE_BINARY_DIR : ${CMAKE_BINARY_DIR}/tests") # this is required for running support.hpp that is made in the tests directory
include_directories("${CMAKE_BINARY_DIR}/tests")

set(CPP_FILES)
addToIt(SphericalHarmonics)
addToIt(SphericalBessels)
addToIt(AtomisticModel)
addToIt(Fit)
addToIt(Waters)
addToIt(functions_WETSAXS)
addToIt(SBESJY)
addToIt(Ensemble)

include_directories(../src/base)
include_directories(../src/fsrc)

MESSAGE("** SOURCES => ${SOURCES}")
MESSAGE("** MAIN DIR => ${CMAKE_SOURCE_DIR}")

set(TEST_OBJ ${CMAKE_PROJECT_NAME}_test)

add_executable(${TEST_OBJ} ${CPP_FILES})

target_link_libraries(${TEST_OBJ} PRIVATE ${BOOST_FILESYSTEM} ${BOOST_PO} ${BOOST_RANDOM} ${BOOST_REGEX} ${BOOST_SYSTEM} ${SASTools_LIB} ${CMAKE_PROJECT_NAME}_lib nnls_lib gtest gtest_main)
# target_link_libraries(${TEST_OBJ} gtest gtest_main)