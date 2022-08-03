cmake_minimum_required(VERSION 3.11)

project(rampler VERSION 2.1.1
                LANGUAGES CXX
                DESCRIPTION "Rampler is a tool for subsampling or splitting FASTA/Q files.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(FetchContent)
include(GNUInstallDirs)

find_package(bioparser 3.0.15 QUIET)
if (NOT bioparser_FOUND)
  FetchContent_Declare(
    bioparser
    GIT_REPOSITORY https://github.com/rvaser/bioparser
    GIT_TAG 3.0.15)

  FetchContent_GetProperties(bioparser)
  if (NOT bioparser_POPULATED)
    FetchContent_Populate(bioparser)
    add_subdirectory(
      ${bioparser_SOURCE_DIR}
      ${bioparser_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

find_package(biosoup 0.10.0 QUIET)
if (NOT biosoup_FOUND)
  FetchContent_Declare(
    biosoup
    GIT_REPOSITORY https://github.com/rvaser/biosoup
    GIT_TAG 0.10.0)

  FetchContent_GetProperties(biosoup)
  if (NOT biosoup_POPULATED)
    FetchContent_Populate(biosoup)
    add_subdirectory(
      ${biosoup_SOURCE_DIR}
      ${biosoup_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

add_executable(rampler
  src/main.cpp
  src/sampler.cpp)

target_link_libraries(rampler
  bioparser::bioparser
  biosoup::biosoup)

target_compile_definitions(rampler PRIVATE VERSION="${PROJECT_VERSION}")

install(TARGETS rampler DESTINATION ${CMAKE_INSTALL_BINDIR})
