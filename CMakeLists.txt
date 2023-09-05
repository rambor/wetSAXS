cmake_minimum_required(VERSION 3.21)
project(wetSAXS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_Release "${CMAKE_CXX_FLAGS_RELEASE} -Wall -g1 -DNDEBUG -O3 ")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native -pthread -ffast-math -fopenmp")
MESSAGE("** SYSTEM : ${CMAKE_SYSTEM_NAME}")

SET (SASTOOLS_INCLUDE_DIRS $ENV{HOME}/usr/local/include)
SET(USER_LOCAL_DIR $ENV{HOME}/usr/local)
include_directories(${SASTOOLS_INCLUDE_DIRS})

find_library(SASTools_LIB libSASTools.a PATHS ~/usr/local/lib ~/lib)
#
# BOOST SECTION
#
# set Boost directory
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)

# Specify the -DBOOST_LIBRARYDIR=/Users/xos81802/src_libs/boost/lib in shell
# Boost root variable set as environment variable in CLION Preferences
find_package(Boost)
MESSAGE("** BOOST ROOT => ${BOOSTROOT}")
MESSAGE("** BOOST ROOT => ${BOOST_ROOT}")
MESSAGE("** BOOST INCLUDE DIRECTORY: ${Boost_INCLUDE_DIR}")
MESSAGE("** BOOST LIB DIRECTORY: ${Boost_LIBRARY_DIR}")
MESSAGE("** BOOST LIB DIRECTORY: ${BOOST_LIBRARYDIR}")
include_directories(${Boost_INCLUDE_DIR})

set(BOOST_FILESYSTEM ${BOOST_LIBRARYDIR}/libboost_filesystem.a )
set(BOOST_PO ${BOOST_LIBRARYDIR}/libboost_program_options.a )
set(BOOST_RANDOM ${BOOST_LIBRARYDIR}/libboost_random.a )
set(BOOST_REGEX ${BOOST_LIBRARYDIR}/libboost_regex.a )
set(BOOST_SYSTEM ${BOOST_LIBRARYDIR}/libboost_system.a )

MESSAGE("** BOOST SPECIAL => ${BOOST_FILESYSTEM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_PO}")
MESSAGE("** BOOST SPECIAL => ${BOOST_RANDOM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_REGEX}")
MESSAGE("** BOOST SPECIAL => ${BOOST_SYSTEM}")
MESSAGE("** BOOST INCLUDE DIRECTORY: ${Boost_INCLUDE_DIR}")

SET(SOURCES)
SET(HEADERS)

add_subdirectory(src/base)

include_directories(src/base)
add_executable(wetSAXS main.cpp)

target_link_libraries(wetSAXS
        PUBLIC
        ${BOOST_FILESYSTEM} ${BOOST_PO} ${BOOST_RANDOM} ${BOOST_REGEX} ${BOOST_SYSTEM}
        ${SASTools_LIB} ${CMAKE_PROJECT_NAME}_lib
        )

add_subdirectory(lib) # testing directory