cmake_minimum_required (VERSION 2.6)
project (freestream)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

# OpenMP support, based on MUSIC/CMakeLists.txt
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    message("Using Clang compiler without OpenMP parallelization... ")
else ()
    find_package(OpenMP)
    if (OpenMP_FOUND)
        string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
    endif ()
endif()

find_package (GSL REQUIRED)

string(APPEND CMAKE_CXX_FLAGS " ${GSL_CXX_FLAGS} -g -O3 -Wall")

SET(INCLUDE_DIRS ${GSL_INCLUDE_DIR})
SET(LIBS ${GSL_LIBRARIES})

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

add_subdirectory(src)

