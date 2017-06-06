cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(h1amg CXX)

set(NETGEN_BASE_DIR $ENV{NETGENDIR}/.. CACHE PATH "Netgen install dir")
set(NETGEN_LIB_DIR ${NETGEN_BASE_DIR}/lib)

find_package(NGSolve CONFIG HINTS ${NETGEN_BASE_DIR}/share/cmake ${NETGEN_BASE_DIR}/Resources/share/cmake)
find_package(Netgen CONFIG HINTS ${NETGEN_BASE_DIR}/share/cmake ${NETGEN_BASE_DIR}/Resources/share/cmake)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NGSOLVE_COMPILE_FLAGS}" )
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${NGSOLVE_LINK_FLAGS}" )
set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${NGSOLVE_LINK_FLAGS}" )
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${NGSOLVE_LINK_FLAGS}" )

set(CMAKE_INSTALL_PREFIX ${NETGEN_BASE_DIR} CACHE INTERNAL "Prefix prepended to install directories" FORCE)

include_directories(${NGSOLVE_INCLUDE_DIRS} "${CMAKE_SOURCE_DIR}/include")

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
