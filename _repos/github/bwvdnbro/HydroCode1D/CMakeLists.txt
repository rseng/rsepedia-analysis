################################################################################
# This file is part of HydroCode1D
# Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#               2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# HydroCode1D is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HydroCode1D is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with HydroCode1D. If not, see <http://www.gnu.org/licenses/>.
################################################################################

cmake_minimum_required(VERSION 3.0)

project(HydroCode1D)

set(SOURCES
    HydroCode1D.cpp
    ${PROJECT_BINARY_DIR}/GitInfo.cpp
    ${PROJECT_BINARY_DIR}/Parameters.cpp

    Boundaries.hpp
    Cell.hpp
    DerivedParameters.hpp
    EOS.hpp
    HLLCRiemannSolver.hpp
    IC.hpp
    OptionNames.hpp
    PhysicalConstants.hpp
    Potential.hpp
    RiemannSolver.hpp
    SafeParameters.hpp
    Spherical.hpp
    Timer.hpp

    UserInput.hpp.in
    UserInput.hpp.sod
)

set(CONFIGURATION_NUMBER "0")

# macro that checks if a configuration value was set on the command line
# if not, the default value is set and a warning is printed
# we also print the actual value that is used
macro(check_configuration_option option_name option_default_value)
  if(NOT DEFINED ${option_name})
    message(STATUS "Using default ${option_name}")
    set(${option_name} ${option_default_value})
  endif(NOT DEFINED ${option_name})
  message(STATUS "${option_name}: ${${option_name}}")
  set(CONFIGURATION_NUMBER "${CONFIGURATION_NUMBER}+1")
  set(CONFIGURATION_KEYS "${CONFIGURATION_KEYS}\"${option_name}\",")
  set(CONFIGURATION_VALUES "${CONFIGURATION_VALUES}\"${${option_name}}\",")
endmacro(check_configuration_option)

# Configurable parameters

check_configuration_option(rmin_in_m 0.)
check_configuration_option(rmax_in_m 1.)
check_configuration_option(ncell 1000)
check_configuration_option(gamma 5./3.)
check_configuration_option(maxtime_in_s 1.)
check_configuration_option(number_of_snaps 1)
check_configuration_option(eos "EOS_IDEAL")
check_configuration_option(boundaries "BOUNDARIES_OPEN")
check_configuration_option(isothermal_temperature_in_k 500.)
check_configuration_option(potential "POTENTIAL_POINT_MASS")
check_configuration_option(g_internal 1.)
check_configuration_option(mass_point_mass_in_kg 1.)
check_configuration_option(courant_factor 0.4)
check_configuration_option(riemannsolver_type "RIEMANNSOLVER_TYPE_EXACT")
check_configuration_option(dimensionality "DIMENSIONALITY_1D")
check_configuration_option(hydro_order "HYDRO_ORDER_2")
check_configuration_option(timeline_type "TIMELINE_NORMAL")

configure_file(${PROJECT_SOURCE_DIR}/Parameters.hpp.in
               ${PROJECT_BINARY_DIR}/Parameters.hpp @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/Parameters.cpp.in
               ${PROJECT_BINARY_DIR}/Parameters.cpp @ONLY)

if(NOT DEFINED predefined_test)
  message(STATUS "No predefined test chosen, copying empty UserInput.hpp.")
  configure_file(${PROJECT_SOURCE_DIR}/UserInput.hpp.in
                 ${PROJECT_BINARY_DIR}/UserInput.hpp @ONLY)
else(NOT DEFINED predefined_test)
  set(TEST_FILENAME "${PROJECT_SOURCE_DIR}/UserInput.hpp.${predefined_test}")
  if(EXISTS ${TEST_FILENAME})
    message(STATUS "Copying predefined test ${predefined_test}.")
    configure_file(${TEST_FILENAME}
                   ${PROJECT_BINARY_DIR}/UserInput.hpp @ONLY)
  else(EXISTS ${TEST_FILENAME})
    message(ERROR "Predefined test ${predefined_test} not found!")
  endif(EXISTS ${TEST_FILENAME})
endif(NOT DEFINED predefined_test)

include_directories(${PROJECT_BINARY_DIR})

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else(OPENMP_FOUND)
  message(ERROR "OpenMP not found!")
endif(OPENMP_FOUND)

find_path(FFTW_INCLUDES fftw3.h)
find_library(FFTW_LIBRARIES NAMES fftw3)
include_directories(${FFTW_INCLUDES})
link_directories(${FFTW_LIBRARIES})

find_package(Git)

find_package(GSL)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  link_directories(${GSL_LIBRARIES})
else(GSL_FOUND)
  message(ERROR "GSL not found!")
endif(GSL_FOUND)

# Add the command and target that creates the CompilerInfo.cpp file
# Note that we need to add the 'dummy_file_that_is_not_created' to the outputs
# to force rerunnning this command every time the code is compiled
add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/GitInfo.cpp
                          dummy_file_that_is_not_created
                   COMMENT "Generating GitInfo.cpp"
                   COMMAND ${CMAKE_COMMAND}
                           -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
                           -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                           -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
                           -P ${PROJECT_SOURCE_DIR}/write_git_info.cmake
                   DEPENDS ${PROJECT_SOURCE_DIR}/.git/HEAD
                           ${PROJECT_SOURCE_DIR}/GitInfo.cpp.in
                   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_custom_target(GitInfo
                  DEPENDS ${PROJECT_BINARY_DIR}/GitInfo.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -ftree-vectorize")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wshadow")

add_executable(HydroCode1D ${SOURCES})
target_link_libraries(HydroCode1D ${FFTW_LIBRARIES})
target_link_libraries(HydroCode1D ${GSL_LIBRARIES})
