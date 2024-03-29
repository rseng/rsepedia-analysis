# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required(VERSION 3.15)

# set the project name
project(density-fitness VERSION 1.0.11 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFile)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(CTest)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

option(CCP4 "The location where ccp4 is installed" "")

if(DEFINED CCP4 AND EXISTS ${CCP4})
	set(CLIBD ${CCP4}/lib/data)
elseif(EXISTS $ENV{CCP4})
	set(CCP4 $ENV{CCP4})
	set(CLIBD ${CCP4}/lib/data)
endif()

if(EXISTS "${CCP4}")
	set(BUILD_SHARED_LIBS ON)

	list(PREPEND CMAKE_MODULE_PATH "${CCP4}/Lib")
	list(APPEND CMAKE_PREFIX_PATH ${CCP4})

	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
		set(CMAKE_INSTALL_PREFIX ${CCP4})
	endif()
endif()

if(MSVC)
	# make msvc standards compliant...
	add_compile_options(/permissive-)
	add_link_options(/NODEFAULTLIB:library)

	macro(get_WIN32_WINNT version)
		if(WIN32 AND CMAKE_SYSTEM_VERSION)
			set(ver ${CMAKE_SYSTEM_VERSION})
			string(REPLACE "." "" ver ${ver})
			string(REGEX REPLACE "([0-9])" "0\\1" ver ${ver})

			set(${version} "0x${ver}")
		endif()
	endmacro()

	get_WIN32_WINNT(ver)
	add_definitions(-D_WIN32_WINNT=${ver})
endif()

# Create a revision file, containing the current git version info
include(VersionString)
write_version_header("${CMAKE_CURRENT_SOURCE_DIR}/src")

# Optionally use mrc to create resources
find_package(Mrc QUIET)

if(MRC_FOUND)
	option(USE_RSRC "Use mrc to create resources" ON)
else()
	message(STATUS "Not using resources since mrc was not found")
endif()

if(USE_RSRC)
	message("Using resources compiled with ${MRC_EXECUTABLE}")
	mrc_write_header(${PROJECT_BINARY_DIR}/mrsrc.hpp)
endif()

# libraries
if(NOT PDB_REDO_META)
	find_package(Eigen3 REQUIRED)
	find_package(libmcfp REQUIRED)
	find_package(zeep REQUIRED)
	find_package(pdb-redo 3.0.1 REQUIRED)
endif()

add_executable(density-fitness
	${PROJECT_SOURCE_DIR}/src/density-fitness.cpp
	${PROJECT_SOURCE_DIR}/src/main.cpp)

if(USE_RSRC)
	mrc_target_resources(density-fitness ${CIFPP_SHARE_DIR}/mmcif_pdbx.dic ${CIFPP_SHARE_DIR}/mmcif_ddl.dic)
endif()

target_include_directories(density-fitness PRIVATE pdb-redo::pdb-redo ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR})
target_link_libraries(density-fitness PRIVATE pdb-redo::pdb-redo libmcfp::libmcfp zeep::zeep)

install(TARGETS ${PROJECT_NAME}
	RUNTIME DESTINATION ${BIN_INSTALL_DIR}
)

# testing
if(BUILD_TESTING)
	find_package(Boost REQUIRED)

	add_executable(unit-test-density-fitness
		${PROJECT_SOURCE_DIR}/test/unit-test-density-fitness.cpp
		${PROJECT_SOURCE_DIR}/src/density-fitness.cpp)

	target_include_directories(unit-test-density-fitness PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_BINARY_DIR})

	target_link_libraries(unit-test-density-fitness pdb-redo::pdb-redo libmcfp::libmcfp zeep::zeep Boost::headers)

	add_test(NAME unit-test-density-fitness COMMAND $<TARGET_FILE:unit-test-density-fitness> -- ${PROJECT_SOURCE_DIR}/test)

	if(USE_RSRC)
		mrc_target_resources(unit-test-density-fitness ${CIFPP_SHARE_DIR}/mmcif_pdbx.dic ${CIFPP_SHARE_DIR}/mmcif_ddl.dic)
	endif()
endif()

# manual
if(UNIX)
	install(FILES doc/density-fitness.1
		DESTINATION ${CMAKE_INSTALL_DATADIR}/man/man1)
endif()

if(EXISTS "${CCP4}/html")
	install(FILES doc/density-fitness.html
		DESTINATION ${CCP4}/html)
endif()
