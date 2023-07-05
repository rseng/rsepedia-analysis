# Copyright (C) 2018 ENKI Project
#
# Please see the AUTHORS file in the main source directory for a full list
# of contributors.
#
# This file is part of Thermocodegen
#
# Thermocodegen is free software: you can redistribute it and/or modify
# it under the terms of the MIT License
#
# Thermocodegen is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# License for more details.
#
# You should have received a copy of the MIT
# along with Thermocodegen

# Require CMake 3.0
cmake_minimum_required(VERSION 3.0) 

# Set verbose output while testing CMake
set(CMAKE_VERBOSE_MAKEFILE 1)

# Set CMake behavior
#cmake_policy(SET CMP0004 OLD)


# Default build type (can be overridden by user)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

# The project name for thermocodegen (sets up environment variables for binary and source directories)
project(THERMOCODEGEN)
# Report the build type
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

find_package(PythonInterp 3 REQUIRED)

# Set variables to help find Python library that is compatible with
# interpreter
if (PYTHONINTERP_FOUND)
  # Get Python include path from Python interpretter
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                         "import distutils.sysconfig, sys; sys.stdout.write(distutils.sysconfig.get_python_inc())"
                  OUTPUT_VARIABLE _PYTHON_INCLUDE_PATH
                  RESULT_VARIABLE _PYTHON_INCLUDE_RESULT)

  # Get Python library path from interpreter
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                        "import os, sys, inspect; sys.stdout.write(os.path.split(os.path.split(inspect.getfile(inspect))[0])[0])"
                  OUTPUT_VARIABLE _PYTHON_LIB_PATH
                  RESULT_VARIABLE _PYTHON_LIB_RESULT)

  # Set include path, if returned by interpreter
  if ("${_PYTHON_INCLUDE_RESULT}" STREQUAL "0")
    set(PYTHON_INCLUDE_DIR ${_PYTHON_INCLUDE_PATH})
  else()
    message(FATAL_ERROR "Failed to find python include path.")
  endif()

  # Add a search path for Python library based on output from
  # iterpreter
  if ("${_PYTHON_LIB_RESULT}" STREQUAL "0")
    set(CMAKE_LIBRARY_PATH ${_PYTHON_LIB_PATH})
  else()
    message(FATAL_ERROR "Failed to find python library path.")
  endif()
else()
  message(FATAL_ERROR "Could not find python interpreter.")
endif()

# Find Pythons libs
find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)

if(NOT ${PYTHONLIBS_FOUND})
  message(FATAL_ERROR "Could not find python libraries.")
endif(NOT ${PYTHONLIBS_FOUND})

list(APPEND THERMOCODEGEN_DEP_INCLUDE_DIRECTORIES "${PYTHON_INCLUDE_DIRS}")
list(APPEND THERMOCODEGEN_TARGET_LINK_LIBRARIES "${PYTHON_LIBRARIES}")
list(APPEND THERMOCODEGEN_CXX_DEFINITIONS "-DHAS_PYTHON")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
find_package(Spud REQUIRED)

if(NOT ${SPUD_FOUND})
  message(FATAL_ERROR "Could not find Spud.")
endif(NOT ${SPUD_FOUND})

list(APPEND THERMOCODEGEN_DEP_INCLUDE_DIRECTORIES "${SPUD_INCLUDE_DIRS}")
list(APPEND THERMOCODEGEN_TARGET_LINK_LIBRARIES "${SPUD_LIBRARIES}")
list(APPEND THERMOCODEGEN_CXX_DEFINITIONS "-DHAS_SPUD")


#install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

if (PYTHONINTERP_FOUND)
  set(PYTHON_LIB_DISCOVER_STR "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))")

  # Probe Python interpreter
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "${PYTHON_LIB_DISCOVER_STR}"
    RESULT_VARIABLE _THERMOCODEGEN_PYTHON_INSTALL_RESULT
    OUTPUT_VARIABLE _THERMOCODEGEN_PYTHON_INSTALL_DIR
    )
  
  if ("${_THERMOCODEGEN_PYTHON_INSTALL_RESULT}" STREQUAL "0")
    set(THERMOCODEGEN_PYTHON_INSTALL_DIR ${_THERMOCODEGEN_PYTHON_INSTALL_DIR})
  else()
    message(FATAL_ERROR "Failed to find python installation path.")
  endif()
endif()

if (THERMOCODEGEN_PYTHON_INSTALL_DIR)
  install(DIRECTORY ${PROJECT_SOURCE_DIR}/python/thermocodegen DESTINATION ${THERMOCODEGEN_PYTHON_INSTALL_DIR})
endif()

# Version information
#
# Write template header files to disk at configuration time
# (these hold the template for the git sha and version information that will
#  be replaced at build time by a suitable variable - using the script
#  below)
file(
    WRITE ${PROJECT_BINARY_DIR}/tcgversion.h.in
"#define TCG_GIT_SHA         \${GIT_SHA}
#define TCG_VERSION_RELEASE \${TCG_VERSION_RELEASE}
#define TCG_VERSION_MAJOR   \${TCG_VERSION_MAJOR}
#define TCG_VERSION_MINOR   \${TCG_VERSION_MINOR}
#define TCG_VERSION_PATCH   \${TCG_VERSION_PATCH}
#define TCG_VERSION         \"\${TCG_VERSION}\"
"
    )

file(
    WRITE ${PROJECT_BINARY_DIR}/git_sha.py.in
"git_sha = \${GIT_SHA}
"
    )

# Write a cmake script to disk at configuration time
# (this will run git - at compile time - and record the id info
#  to an environment variable which can be used to write a
#  header file for the build)
file(
    WRITE ${PROJECT_BINARY_DIR}/tcgversion.cmake
    "find_package(Git)
if(GIT_FOUND)
    execute_process (
       COMMAND \${GIT_EXECUTABLE} log --pretty=format:\"%H %ad\" -1 HEAD
       WORKING_DIRECTORY \${WORKDIR}
       OUTPUT_VARIABLE GIT_SHA
       ERROR_VARIABLE GIT_ERR
       OUTPUT_STRIP_TRAILING_WHITESPACE
       ERROR_STRIP_TRAILING_WHITESPACE
       )
    if(GIT_ERR)
       set(GIT_SHA \\\"git_log_not_found\\\")
    endif(GIT_ERR)
else(GIT_FOUND)
    set(GIT_SHA \\\"git_not_found\\\")
endif(GIT_FOUND)
find_package(PythonInterp 3 REQUIRED)
if (PYTHONINTERP_FOUND)
    execute_process(
       COMMAND \${PYTHON_EXECUTABLE} -c \"import sys; sys.path.append(\\\"\${WORKDIR}/python/thermocodegen\\\"); import version; sys.stdout.write(repr(version.version_release))\"
       OUTPUT_VARIABLE TCG_VERSION_RELEASE
       ERROR_VARIABLE PYTHON_ERR
       OUTPUT_STRIP_TRAILING_WHITESPACE
       ERROR_STRIP_TRAILING_WHITESPACE
       )
    if(PYTHON_ERR)
       set(TCG_VERSION_RELEASE \\\"version_not_found\\\")
    endif(PYTHON_ERR)
    execute_process(
       COMMAND \${PYTHON_EXECUTABLE} -c \"import sys; sys.path.append(\\\"\${WORKDIR}/python/thermocodegen\\\"); import version; sys.stdout.write(repr(version.version_major))\"
       OUTPUT_VARIABLE TCG_VERSION_MAJOR
       ERROR_VARIABLE PYTHON_ERR
       OUTPUT_STRIP_TRAILING_WHITESPACE
       ERROR_STRIP_TRAILING_WHITESPACE
       )
    if(PYTHON_ERR)
       set(TCG_VERSION_MAJOR \\\"version_not_found\\\")
    endif(PYTHON_ERR)
    execute_process(
       COMMAND \${PYTHON_EXECUTABLE} -c \"import sys; sys.path.append(\\\"\${WORKDIR}/python/thermocodegen\\\"); import version; sys.stdout.write(repr(version.version_minor))\"
       OUTPUT_VARIABLE TCG_VERSION_MINOR
       ERROR_VARIABLE PYTHON_ERR
       OUTPUT_STRIP_TRAILING_WHITESPACE
       ERROR_STRIP_TRAILING_WHITESPACE
       )
    if(PYTHON_ERR)
       set(TCG_VERSION_MINOR \\\"version_not_found\\\")
    endif(PYTHON_ERR)
    execute_process(
       COMMAND \${PYTHON_EXECUTABLE} -c \"import sys; sys.path.append(\\\"\${WORKDIR}/python/thermocodegen\\\"); import version; sys.stdout.write(repr(version.version_patch))\"
       OUTPUT_VARIABLE TCG_VERSION_PATCH
       ERROR_VARIABLE PYTHON_ERR
       OUTPUT_STRIP_TRAILING_WHITESPACE
       ERROR_STRIP_TRAILING_WHITESPACE
       )
    if(PYTHON_ERR)
       set(TCG_VERSION_PATCH \\\"version_not_found\\\")
    endif(PYTHON_ERR)
    execute_process(
       COMMAND \${PYTHON_EXECUTABLE} -c \"import sys; sys.path.append(\\\"\${WORKDIR}/python/thermocodegen\\\"); import version; sys.stdout.write(version.version)\"
       OUTPUT_VARIABLE TCG_VERSION
       ERROR_VARIABLE PYTHON_ERR
       OUTPUT_STRIP_TRAILING_WHITESPACE
       ERROR_STRIP_TRAILING_WHITESPACE
       )
    if(PYTHON_ERR)
       set(TCG_VERSION version_not_found)
    endif(PYTHON_ERR)
else(PYTHONINTERP_FOUND)
    set(TCG_VERSION_RELEASE \\\"python_not_found\\\")
    set(TCG_VERSION_MAJOR \\\"python_not_found\\\")
    set(TCG_VERSION_MINOR \\\"python_not_found\\\")
    set(TCG_VERSION_PATCH \\\"python_not_found\\\")
    set(TCG_VERSION python_not_found)
endif(PYTHONINTERP_FOUND)
configure_file(\${SRCDIR}/tcgversion.h.in \${DSTDIR}/tcgversion.h)
configure_file(\${SRCDIR}/git_sha.py.in \${DSTDIR}/git_sha.py)
"
    )

# add a target that runs the above cmake script and creates the header
# file containing the git log information
add_custom_target(
    tcgversion ALL
    ${CMAKE_COMMAND} -D SRCDIR="${PROJECT_BINARY_DIR}"
                     -D DSTDIR="${PROJECT_BINARY_DIR}"
                     -D WORKDIR="${PROJECT_SOURCE_DIR}"
                     -P "${PROJECT_BINARY_DIR}/tcgversion.cmake"
    )

install(FILES ${PROJECT_BINARY_DIR}/tcgversion.h DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/git_sha.py DESTINATION ${THERMOCODEGEN_PYTHON_INSTALL_DIR}/thermocodegen)

# Install schemas
install(DIRECTORY ${PROJECT_SOURCE_DIR}/schemas DESTINATION share/thermocodegen)

# Write a schemata file for emml files that can be installed
file(
    WRITE ${PROJECT_BINARY_DIR}/emml
    "Thermocodegen Endmember Markup Language\ndefault=${CMAKE_INSTALL_PREFIX}/share/thermocodegen/schemas/endmember.rng"
    )

# Write a schemata file for phml files that can be installed
file(
    WRITE ${PROJECT_BINARY_DIR}/phml
    "Thermocodegen Phase Markup Language\ndefault=${CMAKE_INSTALL_PREFIX}/share/thermocodegen/schemas/phase.rng"
    )
    
# Write a schemata file for rxml files that can be installed
file(
    WRITE ${PROJECT_BINARY_DIR}/rxml
    "Thermocodegen Reaction Markup Language\ndefault=${CMAKE_INSTALL_PREFIX}/share/thermocodegen/schemas/reaction.rng"
    )

# install schemata files
install(FILES ${PROJECT_BINARY_DIR}/emml DESTINATION share/thermocodegen/diamond/schemata)
install(FILES ${PROJECT_BINARY_DIR}/phml DESTINATION share/thermocodegen/diamond/schemata)
install(FILES ${PROJECT_BINARY_DIR}/rxml DESTINATION share/thermocodegen/diamond/schemata)


# install templates
install(DIRECTORY ${PROJECT_SOURCE_DIR}/templates DESTINATION share/thermocodegen)

# install examples 
install(DIRECTORY ${PROJECT_SOURCE_DIR}/examples DESTINATION share/thermocodegen)

# install tests
install(DIRECTORY ${PROJECT_SOURCE_DIR}/python/thermocodegen/tests DESTINATION share/thermocodegen)

# install notebooks
# install any binary utilities (might move to utils)
file(GLOB NOTEBOOKS ${PROJECT_SOURCE_DIR}/Notebooks/*.ipynb)
install(FILES ${NOTEBOOKS} DESTINATION share/thermocodegen/Notebooks )

#Cmake configuration files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/cmake DESTINATION share/thermocodegen)

file(GLOB SCRIPTS ${PROJECT_SOURCE_DIR}/scripts/*.py)
file(GLOB EXCLUDESCRIPTS ${PROJECT_SOURCE_DIR}/scripts/.*)
if(EXCLUDESCRIPTS)
  list(REMOVE_ITEM SCRIPTS ${EXCLUDESCRIPTS})
endif()
install(FILES ${SCRIPTS} DESTINATION bin 
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ
	OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)

# install any binary utilities (might move to utils)
file(GLOB PROGS ${PROJECT_SOURCE_DIR}/bin/tcg_*)
install(FILES ${PROGS} DESTINATION bin 
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ
	OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)


configure_file(${PROJECT_SOURCE_DIR}/cmake/templates/THERMOCODEGENConfig.cmake.in
  ${PROJECT_BINARY_DIR}/THERMOCODEGENConfig.cmake @ONLY ESCAPE_QUOTES)
configure_file(${PROJECT_SOURCE_DIR}/cmake/templates/THERMOCODEGENConfigVersion.cmake.in
  ${PROJECT_BINARY_DIR}/THERMOCODEGENConfigVersion.cmake @ONLY ESCAPE_QUOTES)
install(FILES ${PROJECT_BINARY_DIR}/THERMOCODEGENConfig.cmake DESTINATION share/thermocodegen/cmake)
install(FILES ${PROJECT_BINARY_DIR}/THERMOCODEGENConfigVersion.cmake DESTINATION share/thermocodegen/cmake)

# Just check for pybind11
find_package(pybind11 REQUIRED)

if(NOT ${PYBIND11_FOUND})
  message(FATAL_ERROR "Could not find working pybind11.")
endif(NOT ${PYBIND11_FOUND})


# Find GSL
find_package(GSL REQUIRED)



# build and install libswimdew
set(AQUEOUS_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/Thermoengine/thermoengine/thermoengine/aqueous)
set(FREESTEAM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/Thermoengine/src/FreeSteam2.1)

include_directories(${PROJECT_SOURCE_DIR}/ThermoEngine/thermoengine/thermoengine/aqueous ${PROJECT_SOURCE_DIR}/ThermoEngine/src/FreeSteam2.1 ${GSL_INCLUDE_DIRS})
file(GLOB SWIM_SOURCES ${PROJECT_SOURCE_DIR}/ThermoEngine/thermoengine/thermoengine/aqueous/*.c ${PROJECT_SOURCE_DIR}/ThermoEngine/src/FreeSteam2.1/*.c)
file(GLOB SWIM_HEADERS ${PROJECT_SOURCE_DIR}/ThermoEngine/thermoengine/thermoengine/aqueous/*.h ${PROJECT_SOURCE_DIR}/ThermoEngine/src/FreeSteam2.1/*.h)

# Set install RPATH
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
# make sure full install directory name is included in mac Dynamic library
set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib)

add_library(swimdew SHARED ${SWIM_SOURCES})
target_link_libraries(swimdew ${SWIMDEW} ${GSL_LIBRARIES})

install(TARGETS swimdew DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${SWIM_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/swimdew)

# write configuration files
file(
    WRITE ${PROJECT_BINARY_DIR}/thermocodegen.conf
"# Helper file for setting non-default ThermoCodegen environment variables

# Common Unix variables
export THERMOCODEGEN_HOME=${CMAKE_INSTALL_PREFIX}
export THERMOCODEGEN_CMAKE_PATH=${CMAKE_INSTALL_PREFIX}/share/thermocodegen/cmake
export CMAKE_MODULE_PATH=${CMAKE_INSTALL_PREFIX}/share/thermocodegen/cmake/modules:\$CMAKE_MODULE_PATH
export PYTHONPATH=${THERMOCODEGEN_PYTHON_INSTALL_DIR}:\$PYTHONPATH
export DIAMOND_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/share/thermocodegen/diamond:\$DIAMOND_CONFIG_PATH
export PATH=${CMAKE_INSTALL_PREFIX}/bin:\$PATH
export CPLUS_INCLUDE_PATH=${CMAKE_INSTALL_PREFIX}/include/swimdew:\$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:\$LD_LIBRARY_PATH
export LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:\$LIBRARY_PATH
"
    )

install(FILES ${PROJECT_BINARY_DIR}/thermocodegen.conf DESTINATION share/thermocodegen)

file(
    WRITE ${PROJECT_BINARY_DIR}/thermocodegen.configmodule
"#%Module1.0
proc ModulesHelp { } {
        global dotversion

        puts stderr \"\tSets environment for ThermoCodegen.\"
}
module-whatis \"Sets environment for Thermocodegen\"

setenv       THERMOCODEGEN_HOME        ${CMAKE_INSTALL_PREFIX}
setenv	     THERMOCODEGEN_CMAKE_PATH  ${CMAKE_INSTALL_PREFIX}/share/thermocodegen/cmake
prepend-path CMAKE_MODULE_PATH         ${CMAKE_INSTALL_PREFIX}/share/thermocodegen/cmake/modules
prepend-path PYTHONPATH                ${THERMOCODEGEN_PYTHON_INSTALL_DIR}
prepend-path DIAMOND_CONFIG_PATH       ${CMAKE_INSTALL_PREFIX}/share/thermocodegen/diamond
prepend-path PATH                      ${CMAKE_INSTALL_PREFIX}/bin
prepend-path CPLUS_INCLUDE_PATH        ${CMAKE_INSTALL_PREFIX}/include/swimdew
prepend-path LD_LIBRARY_PATH           ${CMAKE_INSTALL_PREFIX}/lib
prepend-path LIBRARY_PATH              ${CMAKE_INSTALL_PREFIX}/lib
"
    )

install(FILES ${PROJECT_BINARY_DIR}/thermocodegen.configmodule DESTINATION share/thermocodegen)


# print post-install message
add_subdirectory(cmake/post-install)
