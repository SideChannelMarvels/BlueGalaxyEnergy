cmake_minimum_required(VERSION 3.15)

set(CMAKE_BUILD_TYPE Debug CACHE STRING "Used default build type" FORCE)

project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION})

set(CMAKE_CXX_FLAGS "-g -O3 -Wall -Wextra")
set(CMAKE_CXX_STANDARD 17)

find_package(pybind11 CONFIG REQUIRED)
find_library(PTHREAD_LIB pthread REQUIRED)
find_library(GMP_LIB gmp REQUIRED)
find_library(NTL_LIB ntl REQUIRED)

pybind11_add_module(_core MODULE src/main.cpp)

target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_link_libraries(_core PRIVATE ${NTL_LIB} ${GMP_LIB} ${PTHREAD_LIB})

install(TARGETS _core DESTINATION bluegalaxyenergy)
