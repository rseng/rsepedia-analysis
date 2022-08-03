# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

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
project(tortoize VERSION 2.0.7 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(Dart)
include(FindFilesystem)
include(GenerateExportHeader)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Filesystem REQUIRED)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

if(NOT "$ENV{CCP4}" STREQUAL "")
	set(BUILD_SHARED_LIBS ON)
	set(CCP4 $ENV{CCP4})
	list(PREPEND CMAKE_MODULE_PATH "${CCP4}/Lib")
	list(APPEND CMAKE_PREFIX_PATH ${CCP4})
	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
		set(CMAKE_PREFIX_PATH ${CCP4})
	endif()
endif()

if(MSVC)
    # make msvc standards compliant...
    add_compile_options(/permissive-)

	macro(get_WIN32_WINNT version)
		if (WIN32 AND CMAKE_SYSTEM_VERSION)
			set(ver ${CMAKE_SYSTEM_VERSION})
			string(REPLACE "." "" ver ${ver})
			string(REGEX REPLACE "([0-9])" "0\\1" ver ${ver})

			set(${version} "0x${ver}")
		endif()
	endmacro()

	get_WIN32_WINNT(ver)
	add_definitions(-D_WIN32_WINNT=${ver})

	# On Windows, do not install in the system location
	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND NOT BUILD_FOR_CCP4)
		message(STATUS "The library and auxiliary files will be installed in $ENV{LOCALAPPDATA}/${PROJECT_NAME}")
		set(CMAKE_INSTALL_PREFIX "$ENV{LOCALAPPDATA}/${PROJECT_NAME}" CACHE PATH "..." FORCE)
	endif()
endif()

if(UNIX AND NOT APPLE)
	# On Linux, install in the $HOME/.local folder by default
	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
		message(WARNING "The library and auxiliary files will be installed in $ENV{HOME}/.local")
		set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "..." FORCE)
	endif()
endif()

# Create a revision file, containing the current git version info, if any
include(VersionString)
write_version_header()

# Optionally use mrc to create resources
find_package(Mrc)

if(MRC_FOUND)
	option(USE_RSRC "Use mrc to create resources" ON)
else()
	message(WARNING "Not using resources since mrc was not found")
endif()

if(USE_RSRC STREQUAL "ON")
	set(USE_RSRC 1)
	add_compile_definitions(USE_RSRC)
endif()

# Optionally build a webservice

option(BUILD_WEBSERVICE "Build a version with a webservice daemon" OFF)

# libraries

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads)

find_package(cifpp 2.0.0 REQUIRED)
find_package(zeep 5.1.4 REQUIRED)

find_package(Boost COMPONENTS program_options REQUIRED)

add_executable(tortoize
	${PROJECT_SOURCE_DIR}/src/tortoize.cpp
	${PROJECT_SOURCE_DIR}/src/tortoize-main.cpp
	${TORTOIZE_RESOURCE})

if(BUILD_WEBSERVICE)
	target_compile_definitions(tortoize PRIVATE WEBSERVICE)
endif()

if(USE_RSRC)
	mrc_target_resources(tortoize ${CMAKE_SOURCE_DIR}/rsrc/rama-data.bin ${CMAKE_SOURCE_DIR}/rsrc/torsion-data.bin ${CIFPP_SHARE_DIR}/mmcif_pdbx_v50.dic)
endif()

target_include_directories(tortoize PRIVATE cifpp::cifpp ${CMAKE_BINARY_DIR})
target_link_libraries(tortoize cifpp::cifpp zeep::zeep std::filesystem Boost::program_options)

install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION ${BIN_INSTALL_DIR}
)

if(NOT USE_RSRC)
	install(FILES ${CMAKE_SOURCE_DIR}/rsrc/rama-data.bin ${CMAKE_SOURCE_DIR}/rsrc/torsion-data.bin
		DESTINATION ${CIFPP_SHARE_DIR})
endif()

# manual

install(FILES doc/tortoize.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

# Test code

add_executable(unit-test ${PROJECT_SOURCE_DIR}/test/unit-test.cpp ${PROJECT_SOURCE_DIR}/src/tortoize.cpp)
target_include_directories(unit-test PRIVATE cifpp::cifpp ${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR})
target_link_libraries(unit-test cifpp::cifpp zeep::zeep std::filesystem)

add_test(NAME unit-test COMMAND $<TARGET_FILE:unit-test> -- ${CMAKE_SOURCE_DIR}/test)

if(USE_RSRC)
	mrc_target_resources(unit-test ${CMAKE_SOURCE_DIR}/rsrc/rama-data.bin ${CMAKE_SOURCE_DIR}/rsrc/torsion-data.bin ${CIFPP_SHARE_DIR}/mmcif_pdbx_v50.dic)
endif()
