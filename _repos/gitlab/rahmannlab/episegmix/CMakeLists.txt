# The minimum required version of cmake
cmake_minimum_required(VERSION 3.16)

# Set required C++ standard
set(CMAKE_CXX_STANDARD, 17)

# The name of the project
project(HMM)

# Use release mode as default build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

if(NOT GTEST_SRC_DIR)
    message(STATUS "GTEST_SRC_DIR not set, disabling unit tests")
else()
    enable_testing()
    add_subdirectory(${GTEST_SRC_DIR} "gtest")
    add_subdirectory(tests)
endif()

# Configure static linking
option(BUILD_SHARED_LIBS "Build a shared (ON) or static (OFF) executable" ON)
if(APPLE)
    option(USE_LIBCPP "Use libc++ instead of libstdc++" ON)
    option(USE_CXXABI "Use cxxabi for clang" ON)
else()
    option(USE_LIBCPP "Use libc++ instead of libstdc++" OFF)
    option(USE_CXXABI "Use cxxabi for clang" OFF)
endif()


if(NOT ${BUILD_SHARED_LIBS})
    #disable -rdynamic
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
    # prefer static libraries
    if(WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    else(WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    endif(WIN32)
endif()


if(__COMPILER_GNU AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # GCC, MINGW
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${GCC_CXX_FLAGS_DEBUG}")
    if(NOT ${BUILD_SHARED_LIBS})
      set(CMAKE_EXE_LINKER_FLAGS "-static")
    endif()
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    if(${USE_LIBCPP})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CLANG_CXX_FLAGS_DEBUG}")
    if(NOT ${BUILD_SHARED_LIBS})
      set(CMAKE_EXE_LINKER_FLAGS "-static")
    endif()
elseif(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${MSVC_CXX_DEBUG_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${MSVC_CXX_RELEASE_FLAGS}")
endif()

# Find Boost
find_package(Boost 1.55.0 COMPONENTS iostreams filesystem system program_options thread REQUIRED)

if(NOT Boost_FOUND)
	message(FATAL_ERROR "Boost version 1.55.0 or higher is required")
endif()

SET(Boost_LIBRARIES
	${Boost_IOSTREAMS_LIBRARY}
	${Boost_FILESYSTEM_LIBRARY}
	${Boost_SYSTEM_LIBRARY}
	${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_THREAD_LIBRARY}
)

include_directories(${Boost_INCLUDE_DIRS})

# Find root
find_package(ROOT COMPONENTS Minuit2)
if(NOT ROOT_FOUND)
	message(FATAL_ERROR "Root is required.")
endif()

# Find OpenMP
find_package(OpenMP)
if(NOT OpenMP_CXX_FOUND)
	message(FATAL_ERROR "OpenMP is required.")
endif()

add_library(HMM SHARED
    HMM.h
    HMM.cpp
    HMM_io.cpp
    adjustableDurationHMM.h
    adjustableDurationHMM.cpp
    log_prob.h
    distribution.h
    distribution.cpp
    sample_statistics.h
    FCN.h
    matrix.h
    matrix_reader.h
    matrix_reader.cpp
    )

target_link_libraries(HMM ${Boost_LIBRARIES} ROOT::Minuit2 OpenMP::OpenMP_CXX -pthread)

# add the executables
add_executable(HMMChromSeg main.cpp)
target_link_libraries(HMMChromSeg HMM)

add_executable(TopologyHMM main_topologyHMM.cpp)
target_link_libraries(TopologyHMM HMM)

add_executable(LogLikelihood main_logL.cpp)
target_link_libraries(LogLikelihood HMM)

# link with libraries
if(NOT WIN32)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND ${USE_CXXABI})
        set(CXX_ABI c++abi)
    endif()
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        target_link_libraries(HMMChromSeg ${CXX_ABI})
   endif()
endif()