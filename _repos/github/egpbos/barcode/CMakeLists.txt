#
# Barcode
# Copyright E.G.P. Bos and F.S. Kitaura
#
# Distributed under the terms of the MIT License.
# The full license is in the file LICENSE, distributed with this software.
#

cmake_minimum_required(VERSION 3.9.0)
# 3.8 for target_compile_features(target PRIVATE/PUBLIC cxx_std_11/17), 3.9.0 for FindOpenMP, all other features are in 3.0 at least (didn't check 2.x)

project(barcode)

#--------------------------------------- CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")
include(DownloadAndSetupFindFFTW)

#--------------------------------------- download and setup FindFFTW
DOWNLOAD_AND_SETUP_FINDFFTW()

#--------------------------------------- global settings
# warnings
add_compile_options(-Wundef -Wshadow -Wwrite-strings -Wredundant-decls -Woverloaded-virtual -Wcast-qual -Wcast-align -Wpointer-arith -Wconversion -Wold-style-cast -Wall -Wextra)

# RPATH settings, see https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

#--------------------------------------- dependencies
find_package(GSL REQUIRED)
find_package(OpenMP COMPONENTS CXX)
find_package(FFTW REQUIRED)
set(CURSES_NEED_NCURSES TRUE)
find_package(Curses REQUIRED)
#--------------------------------------- build options
include(Options)

#--------------------------------------- targets
add_subdirectory(barlib)
add_subdirectory(barcode)
add_subdirectory(test)
add_subdirectory(tools)
