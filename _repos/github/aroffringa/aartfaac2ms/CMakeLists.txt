cmake_minimum_required(VERSION 2.8.8)

project(aartfaac2ms)

include(CMakeVersionInfo.txt)

if(PORTABLE)
  add_compile_options(-O3 -ggdb -Wvla -Wall -DNDEBUG)
else()
  add_compile_options(-O3 -ggdb -Wvla -Wall -DNDEBUG -march=native)
endif(PORTABLE)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

include_directories(${CMAKE_SOURCE_DIR}/aocommon/include)

find_package(AOFlagger 3 REQUIRED)
include_directories(${AOFLAGGER_INCLUDE_DIR})

find_package(CFITSIO REQUIRED)
include_directories(${CFITSIO_INCLUDE_DIRS})

find_package(Threads REQUIRED)

find_package(Boost 1.55.0 REQUIRED COMPONENTS date_time filesystem)

set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits scimath scimath_f)

include_directories(${CASACORE_INCLUDE_DIRS})

# The following stuff will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

configure_file(version.h.in version.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(aartfaac2ms main.cpp aartfaac2ms.cpp aartfaacms.cpp averagingwriter.cpp fitsuser.cpp fitswriter.cpp mswriter.cpp progressbar.cpp stopwatch.cpp threadedwriter.cpp)
target_link_libraries(aartfaac2ms
	${AOFLAGGER_LIB} ${CASACORE_LIBRARIES}
	${Boost_SYSTEM_LIBRARY} ${Boost_DATE_TIME_LIBRARY}
	Threads::Threads)

add_executable(afedit afedit.cpp)
target_link_libraries(afedit ${CASACORE_LIBRARIES} Threads::Threads)

install(TARGETS aartfaac2ms afedit DESTINATION bin)

message(STATUS "Flags passed to C++ compiler: " ${CMAKE_CXX_FLAGS})
