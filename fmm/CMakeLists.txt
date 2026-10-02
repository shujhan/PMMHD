# Builds libFMMField.a: the U_FMM field solver (FMMField.cpp, fmm_kernels_2d.hpp) plus
# the unmodified BarytreeK sources it needs, compiled with the Kokkos toolchain. PMMHD itself is not
# built here; it links this library when configured with -DENABLE_FMM=ON.
#
#   cmake -S fmm -B build_fmm -DBARYTREEK_DIR=/path/to/BarytreeK -DKokkos_ROOT=/path/to/kokkos/install \
#         -DCMAKE_INSTALL_PREFIX=/path/to/fmm_install
#   cmake --build build_fmm -j && cmake --install build_fmm
#
# CUDA backend: configure with -DCMAKE_CXX_COMPILER=${Kokkos_ROOT}/bin/nvcc_wrapper

cmake_minimum_required(VERSION 3.16)
project(PMMHD_FMM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(BARYTREEK_DIR "" CACHE PATH "BarytreeK source tree")
if(NOT EXISTS ${BARYTREEK_DIR}/src/structs.hpp)
    message(FATAL_ERROR "set -DBARYTREEK_DIR to the BarytreeK source tree")
endif()

find_package(Kokkos REQUIRED CONFIG)

# BarytreeK config header (types, max_degree)
set(PROJECT_SOURCE_DIR_SAVED ${PROJECT_SOURCE_DIR})
set(PROJECT_SOURCE_DIR ${BARYTREEK_DIR})
configure_file(${BARYTREEK_DIR}/barytreek-config.h.in ${CMAKE_CURRENT_BINARY_DIR}/barytreek-config.h)
set(PROJECT_SOURCE_DIR ${PROJECT_SOURCE_DIR_SAVED})

add_library(FMMField STATIC
    FMMField.cpp
    ${BARYTREEK_DIR}/src/tree_construction.cpp
    ${BARYTREEK_DIR}/src/interaction_list.cpp
    ${BARYTREEK_DIR}/src/upward_pass.cpp
    ${BARYTREEK_DIR}/src/downward_pass.cpp)

target_include_directories(FMMField PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${BARYTREEK_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/../inc)

target_link_libraries(FMMField PUBLIC Kokkos::kokkos)

install(TARGETS FMMField ARCHIVE DESTINATION lib)
