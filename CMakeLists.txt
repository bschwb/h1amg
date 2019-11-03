cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(h1amg CXX)

set(NETGEN_BASE_DIR $ENV{NETGENDIR}/.. CACHE PATH "Netgen install dir")

find_package(NGSolve CONFIG HINTS ${NETGEN_BASE_DIR} /Applications/Netgen.app/Contents/Resources/CMake)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

include_directories(${NETGEN_INCLUDE_DIRS} "${CMAKE_SOURCE_DIR}/include")

set(CMAKE_INSTALL_PREFIX ${NGSOLVE_INSTALL_DIR})

if(USE_CCACHE)
  find_program(CCACHE_FOUND ccache)
  if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  endif(CCACHE_FOUND)
endif(USE_CCACHE)

add_subdirectory(src)

option(BUILD_TESTING "" OFF)
include(CTest)
if(BUILD_TESTING)
  include(cmake_modules/PythonLibTest.cmake)
  add_subdirectory(tests)
endif(BUILD_TESTING)
