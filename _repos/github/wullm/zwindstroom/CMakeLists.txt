cmake_minimum_required(VERSION 3.11)
project(zwindstroom VERSION 1.0.0 DESCRIPTION "Library for cosmic neutrino calculations")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(GSL REQUIRED)
mark_as_advanced(pkgcfg_lib_GSL_gsl pkgcfg_lib_GSL_gslcblas pkgcfg_lib_GSL_m)

file(GLOB zwindstroom_SOURCES src/*.c)

# Static library
add_library(zwindstroom_static STATIC ${zwindstroom_SOURCES})
add_library(zwindstroom::zwindstroom_static ALIAS zwindstroom_static)
set_target_properties(zwindstroom_static PROPERTIES OUTPUT_NAME zwindstroom)
set_target_properties(zwindstroom_static PROPERTIES C_STANDARD 11)
set_property(TARGET zwindstroom_static PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(zwindstroom_static PRIVATE m)
target_link_libraries(zwindstroom_static PRIVATE GSL::gsl)
target_include_directories(zwindstroom_static PUBLIC include)

# Shared library
add_library(zwindstroom_shared SHARED ${zwindstroom_SOURCES})
add_library(zwindstroom::zwindstroom_shared ALIAS zwindstroom_shared)
set_target_properties(zwindstroom_shared PROPERTIES OUTPUT_NAME zwindstroom)
set_target_properties(zwindstroom_shared PROPERTIES C_STANDARD 11)
set_property(TARGET zwindstroom_shared PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(zwindstroom_shared PRIVATE m)
target_link_libraries(zwindstroom_shared PRIVATE GSL::gsl)
target_include_directories(zwindstroom_shared PUBLIC include)

add_executable(test_zwindstroom tests/test.c)
target_link_libraries(test_zwindstroom PRIVATE zwindstroom::zwindstroom_static)
