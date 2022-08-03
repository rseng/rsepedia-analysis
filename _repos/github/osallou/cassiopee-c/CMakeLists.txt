cmake_minimum_required (VERSION 2.8)

enable_language(C CXX)
project (Cassiopee)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

link_directories(/usr/local/lib)


set (Cassiopee_VERSION_MAJOR 1)
set (Cassiopee_VERSION_MINOR 0)

set (PROJECT_SOURCE_DIR src)
set (PROJECT_BINARY_DIR ${CMAKE_BINARY_DIR}/bin)
set (PROJECT_DOC_DIR doc)

enable_testing()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
#add_definitions( -DBOOST_ALL_NO_LIB )
unset(Boost_INCLUDE_DIR CACHE)
unset(Boost_LIBRARY_DIRS CACHE)
set(Boost_USE_STATIC_LIBS   OFF)

find_package(Boost REQUIRED COMPONENTS iostreams serialization)


# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc ALL
${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

FIND_PACKAGE(OpenMP)
if(OPENMP_FOUND)
set(USE_OPENMP 1)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
set(USE_OPENMP 0)
endif()

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/CassiopeeConfig.h.in"
  "${PROJECT_BINARY_DIR}/CassiopeeConfig.h"
  )
 
# add the binary tree to the search path for include files
# so that we will find CassiopeeConfig.h
include_directories("${PROJECT_BINARY_DIR}")
include_directories(${Boost_INCLUDE_DIR})

add_subdirectory (src)
add_subdirectory (test)


