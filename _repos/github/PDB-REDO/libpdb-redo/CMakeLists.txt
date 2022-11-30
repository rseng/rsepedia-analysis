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

cmake_minimum_required(VERSION 3.16)

# set the project name
project(pdb-redo VERSION 3.0.4 LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Unit tests
option(ENABLE_TESTING "Build test exectuables" OFF)

include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFile)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(Dart)
include(FindPkgConfig)
include(GenerateExportHeader)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

# Build shared libraries by default (not my cup of tea, but hey)
option(BUILD_SHARED_LIBS "Build a shared library instead of a static one" OFF)

# We do not want to write an export file for all our symbols...
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Optionally build a version to be installed inside CCP4
option(BUILD_FOR_CCP4 "Build a version to be installed in CCP4" OFF)

if(BUILD_FOR_CCP4)
	if("$ENV{CCP4}" STREQUAL "" OR NOT EXISTS $ENV{CCP4})
		message(FATAL_ERROR "A CCP4 built was requested but CCP4 was not sourced")
	else()
		list(APPEND CMAKE_MODULE_PATH "$ENV{CCP4}")
		list(APPEND CMAKE_PREFIX_PATH "$ENV{CCP4}")
		set(CMAKE_INSTALL_PREFIX "$ENV{CCP4}")

		# This is the only option:
		if(WIN32)
			set(BUILD_SHARED_LIBS ON)
		endif()
	endif("$ENV{CCP4}" STREQUAL "" OR NOT EXISTS $ENV{CCP4})
endif()

# Check if CCP4 is available
if(EXISTS "$ENV{CCP4}")
	set(CCP4 $ENV{CCP4})
	set(CLIBD ${CCP4}/lib/data)
endif()

if(CCP4 AND NOT CLIBD)
	set(CLIBD ${CCP4}/lib/data)
endif()

if(MSVC)
	# make msvc standards compliant...
	add_compile_options(/permissive-)

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

# Optionally use mrc to create resources
if(WIN32 AND BUILD_SHARED_LIBS)
	message("Not using resources when building shared libraries for Windows")
else()
	find_program(MRC mrc)

	if(MRC)
		option(USE_RSRC "Use mrc to create resources" ON)
	else()
		message("Using resources not possible since mrc was not found")
	endif()

	if(USE_RSRC STREQUAL "ON")
		message("Using resources compiled with ${MRC}")
	endif()
endif()

# Libraries
find_package(CCP4 REQUIRED ccp4c clipper-core clipper-ccp4 clipper-contrib)

if(NOT PDB_REDO_META)
	find_package(newuoa REQUIRED)
	find_package(cifpp 5.0.4 REQUIRED)
endif()

pkg_check_modules(GSL gsl REQUIRED IMPORTED_TARGET)

# Create a revision file, containing the current git version info
include(VersionString)
write_version_header("${PROJECT_SOURCE_DIR}/src" "LibPDBREDO")

# Sources
set(project_sources
	${PROJECT_SOURCE_DIR}/src/AtomShape.cpp
	${PROJECT_SOURCE_DIR}/src/BondMap.cpp
	${PROJECT_SOURCE_DIR}/src/ClipperWrapper.cpp
	${PROJECT_SOURCE_DIR}/src/Compound.cpp
	${PROJECT_SOURCE_DIR}/src/DistanceMap.cpp
	${PROJECT_SOURCE_DIR}/src/MapMaker.cpp
	${PROJECT_SOURCE_DIR}/src/Ramachandran.cpp
	${PROJECT_SOURCE_DIR}/src/ResolutionCalculator.cpp
	${PROJECT_SOURCE_DIR}/src/SkipList.cpp
	${PROJECT_SOURCE_DIR}/src/Statistics.cpp
	${PROJECT_SOURCE_DIR}/src/Symmetry-2.cpp
	${PROJECT_SOURCE_DIR}/src/Minimizer.cpp
	${PROJECT_SOURCE_DIR}/src/Restraints.cpp
)

set(project_headers
	${PROJECT_SOURCE_DIR}/include/pdb-redo/AtomShape.hpp
	${PROJECT_SOURCE_DIR}/include/pdb-redo/BondMap.hpp
	${PROJECT_SOURCE_DIR}/include/pdb-redo/ClipperWrapper.hpp
	${PROJECT_SOURCE_DIR}/include/pdb-redo/Compound.hpp
	${PROJECT_SOURCE_DIR}/include/pdb-redo/DistanceMap.hpp
	${PROJECT_SOURCE_DIR}/include/pdb-redo/MapMaker.hpp
	${PROJECT_SOURCE_DIR}/include/pdb-redo/Ramachandran.hpp
	${PROJECT_SOURCE_DIR}/include/pdb-redo/ResolutionCalculator.hpp
	${PROJECT_SOURCE_DIR}/include/pdb-redo/SkipList.hpp
	${PROJECT_SOURCE_DIR}/include/pdb-redo/Statistics.hpp
	${PROJECT_SOURCE_DIR}/include/pdb-redo/Symmetry-2.hpp
	${PROJECT_SOURCE_DIR}/include/pdb-redo/Minimizer.hpp
	${PROJECT_SOURCE_DIR}/include/pdb-redo/Restraints.hpp
)

add_library(pdb-redo ${project_sources} ${project_headers})
add_library(pdb-redo::pdb-redo ALIAS pdb-redo)
set_target_properties(pdb-redo PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(pdb-redo
	PUBLIC
	"$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
	"$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
	${CCP4_INCLUDE_DIRS})
target_link_libraries(pdb-redo PUBLIC cifpp::cifpp newuoa::newuoa ${CCP4_LIBRARIES} ${GSL_LIBRARIES})

if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
	target_link_options(pdb-redo PRIVATE -undefined dynamic_lookup)
endif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")

set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(LIBRARY_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})

# Install rules
install(TARGETS pdb-redo
	EXPORT pdb-redoTargets
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT pdb-redoTargets
	FILE "pdb-redoTargets.cmake"
	NAMESPACE pdb-redo::
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pdb-redo
)

install(
	DIRECTORY include/pdb-redo
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	COMPONENT Devel
)

configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/pdb-redoConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/pdb-redo/pdb-redoConfig.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pdb-redo
	PATH_VARS INCLUDE_INSTALL_DIR LIBRARY_INSTALL_DIR
)

install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/pdb-redo/pdb-redoConfig.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/pdb-redo/pdb-redoConfigVersion.cmake"
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pdb-redo
	COMPONENT Devel
)

set(pdb-redo_MAJOR_VERSION ${CMAKE_PROJECT_VERSION_MAJOR})
set_target_properties(pdb-redo PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${pdb-redo_MAJOR_VERSION}
	INTERFACE_pdb-redo_MAJOR_VERSION ${pdb-redo_MAJOR_VERSION})

set_property(TARGET pdb-redo APPEND PROPERTY
	COMPATIBLE_INTERFACE_STRING pdb-redo_MAJOR_VERSION
)

write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/pdb-redo/pdb-redoConfigVersion.cmake"
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY AnyNewerVersion
)

# pkgconfig support
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix ${CMAKE_INSTALL_PREFIX})
set(libdir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(includedir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})

configure_file(${PROJECT_SOURCE_DIR}/cmake/libpdb-redo.pc.in
	${CMAKE_CURRENT_BINARY_DIR}/libpdb-redo.pc.in @ONLY)
file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libpdb-redo.pc
	INPUT ${CMAKE_CURRENT_BINARY_DIR}/libpdb-redo.pc.in)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libpdb-redo.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

if(ENABLE_TESTING)
	enable_testing()

	if(USE_RSRC)
		add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pdb_redo_test_rsrc.obj
			COMMAND ${MRC} -o ${CMAKE_CURRENT_BINARY_DIR}/pdb_redo_test_rsrc.obj ${CIFPP_SHARE_DIR}/mmcif_pdbx.dic ${COFF_SPEC}
		)
		set(PDB_REDO_TEST_RESOURCE ${CMAKE_CURRENT_BINARY_DIR}/pdb_redo_test_rsrc.obj)
	endif()

	list(APPEND PDB_REDO_tests pdb-redo-unit)

	if(GSL_FOUND)
		list(APPEND PDB_REDO_tests rsr)
	endif()

	foreach(PDB_REDO_TEST IN LISTS PDB_REDO_tests)
		set(PDB_REDO_TEST "${PDB_REDO_TEST}-test")
		set(PDB_REDO_TEST_SOURCE "${PROJECT_SOURCE_DIR}/test/${PDB_REDO_TEST}.cpp")

		add_executable(${PDB_REDO_TEST} ${PDB_REDO_TEST_SOURCE} ${PDB_REDO_TEST_RESOURCE})

		target_include_directories(${PDB_REDO_TEST} PRIVATE
			${PROJECT_SOURCE_DIR}/include
			${CMAKE_CURRENT_BINARY_DIR} # for config.h
		)

		target_link_libraries(${PDB_REDO_TEST} PRIVATE pdb-redo)

		if(MSVC)
			# Specify unwind semantics so that MSVC knowns how to handle exceptions
			target_compile_options(${PDB_REDO_TEST} PRIVATE /EHsc)
		endif()

		add_custom_target("run-${PDB_REDO_TEST}" DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Run${PDB_REDO_TEST}.touch ${PDB_REDO_TEST})

		add_custom_command(
			OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Run${PDB_REDO_TEST}.touch
			COMMAND $<TARGET_FILE:${PDB_REDO_TEST}> -- ${PROJECT_SOURCE_DIR}/test)

		add_test(NAME ${PDB_REDO_TEST}
			COMMAND $<TARGET_FILE:${PDB_REDO_TEST}> -- ${PROJECT_SOURCE_DIR}/test)
	endforeach()
endif()

# Tarball generation
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_SOURCE_TGZ ON)
set(CPACK_SOURCE_TBZ2 OFF)
set(CPACK_SOURCE_TXZ OFF)
set(CPACK_SOURCE_TZ OFF)
set(CPACK_SOURCE_IGNORE_FILES "/data/components.cif;/build;/.vscode;/.git;/regex/doc;/regex/example;/regex/meta;/regex/performance;/regex/test")
set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}-${PROJECT_VERSION}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME ${CPACK_PACKAGE_FILE_NAME})
include(CPack)
