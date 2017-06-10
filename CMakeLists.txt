cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(h1amg CXX)

find_package(Netgen CONFIG HINTS ${NETGEN_DIR} ${NETGEN_CMAKE_DIR})
find_package(NGSolve CONFIG HINTS ${NGSOLVE_DIR} ${NETGEN_CMAKE_DIR})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

include_directories(${NGSOLVE_INCLUDE_DIRS} "${CMAKE_SOURCE_DIR}/include")

set(CMAKE_INSTALL_PREFIX ${NGSOLVE_INSTALL_DIR})

if(USE_CCACHE)
  find_program(CCACHE_FOUND ccache)
  if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  endif(CCACHE_FOUND)
endif(USE_CCACHE)

add_subdirectory(src)

include(CTest)
if(BUILD_TESTING)
  include(cmake_modules/PythonLibTest.cmake)
  add_subdirectory(tests)
endif(BUILD_TESTING)
