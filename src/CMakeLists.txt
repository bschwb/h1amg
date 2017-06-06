add_library(h1amg SHARED
  h1amg.hpp h1amg.cpp
  h1amg_helpers.hpp h1amg_helpers.cpp
  h1.hpp h1.cpp
  build_h1.hpp build_h1.cpp
  h1_helpers.hpp h1_helpers.cpp
  h1_smoothed_prol.hpp h1_smoothed_prol.cpp
  complex_mat.hpp complex_mat.cpp)

target_link_libraries(h1amg ngstd ngcomp ngla ngbla)

set_property(TARGET h1amg PROPERTY CXX_STANDARD 14)
set_property(TARGET h1amg PROPERTY CXX_STANDARD_REQUIRED ON)

install(TARGETS h1amg DESTINATION ${NETGEN_LIB_DIR})
