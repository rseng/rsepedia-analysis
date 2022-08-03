# CMake script for the physamp package
# Author: Julien Dutheil
# Created: 24/10/2014

# Global parameters
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
PROJECT(physamp CXX)

# Compile options
SET (CMAKE_CXX_FLAGS "-std=c++11 -Wall -Weffc++ -Wshadow -Wconversion")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF()

IF(NOT NO_DEP_CHECK)
  SET(NO_DEP_CHECK FALSE CACHE BOOL
      "Disable dependencies check for building distribution only."
      FORCE)
ENDIF(NOT NO_DEP_CHECK)

IF(NO_DEP_CHECK)
  MESSAGE("-- Dependencies checking disabled. Only distribution can be built.")
ELSE(NO_DEP_CHECK)

#static linkage?
IF(NOT BUILD_STATIC)
  SET(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
ENDIF()
IF(BUILD_STATIC)
  MESSAGE("-- Static linkage requested.")
  SET(CMAKE_CXX_FLAGS "-static -static-libgcc -static-libstdc++ ${CMAKE_CXX_FLAGS}")
ENDIF()

# Check compression program
# COMPRESS_PROGRAM controls the choice of program
# COMPRESS_EXT can be used to override the file extension
if (NOT COMPRESS_PROGRAM)
  set (COMPRESS_PROGRAM gzip CACHE STRING "Set program for compressing documentation" FORCE)
endif ()
find_program (COMPRESS_BIN NAMES ${COMPRESS_PROGRAM} DOC "${COMPRESS_PROGRAM} compression program")
if (NOT COMPRESS_BIN)
  message (STATUS "${COMPRESS_PROGRAM} program not found, text doc will not be compressed")
else ()
  # Deduce COMPRESS_EXT for known compression programs if not set
  if (NOT COMPRESS_EXT)
    if (${COMPRESS_PROGRAM} STREQUAL "gzip")
      set (COMPRESS_EXT "gz")
    elseif (${COMPRESS_PROGRAM} STREQUAL "bzip2")
      set (COMPRESS_EXT "bz2")
    else ()
      set (COMPRESS_EXT "${COMPRESS_PROGRAM}") # Default: program name (works for xz/lzma)
    endif ()
  endif ()
  # Generate command line args (always add -c to output compressed file to stdout)
  if (${COMPRESS_PROGRAM} STREQUAL "gzip")
    # -n for no timestamp in files (reproducible builds)
    # -9 for maximum compression (lintian error)
    set (COMPRESS_ARGS -c9n)
  else ()
    set (COMPRESS_ARGS -c)
  endif ()
  message (STATUS "Found ${COMPRESS_BIN} compression program, using file extension .${COMPRESS_EXT}")
endif ()


# Find dependencies (add install directory to search)
IF (CMAKE_INSTALL_PREFIX)
  SET (CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF (CMAKE_INSTALL_PREFIX)

include (GNUInstallDirs)
find_package (bpp-phyl 12.0.0 REQUIRED)


# Subdirectories
ADD_SUBDIRECTORY(physamp)
ADD_SUBDIRECTORY(docs)
ADD_SUBDIRECTORY(man)

ENDIF(NO_DEP_CHECK)

# Packager
SET(CPACK_PACKAGE_NAME "physamp")
SET(CPACK_PACKAGE_VENDOR "Julien Y. Dutheil")
SET(CPACK_PACKAGE_VERSION "1.1.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The Phylogenetic Sampler package")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS")
SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL")
SET(CPACK_SOURCE_GENERATOR "TGZ")
# /!\ This assumes that an external build is used
SET(CPACK_SOURCE_IGNORE_FILES
       "/build/" 
       "/\\\\.git/" 
       "/\\\\.gitignore" 
       ${CPACK_SOURCE_IGNORE_FILES}
       )

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_DEBSOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig")
INCLUDE(CPack)

#This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

IF(NOT NO_DEP_CHECK)
IF (UNIX)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist info man)

#ADD_DEPENDENCIES(info install)

ENDIF(UNIX)

ENDIF(NOT NO_DEP_CHECK)
