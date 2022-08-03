##############
### Header ###
##############

cmake_minimum_required( VERSION 3.10 FATAL_ERROR )
project( HiLive VERSION 2.0 LANGUAGES CXX )

###################
### Build Setup ###
###################

# Set CXX Standard
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_SANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )

# Set runtime output directory
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

# Error on unsupported compiler
if( NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU" )
    message( FATAL_ERROR
        "Compiler id '${CMAKE_CXX_COMPILER_ID}' is not supported, please \
        check the documentation." )
endif()

# Compiling options
option( VERBOSE_CONFIG "Verbose mode for quick build setup debugging." OFF )
option( CONDA "Flag for compilation in conda env." OFF )
option( LZ4_PATH "Path to manually installed LZ4 library." "")
option( SEQAN_PATH "Path to manually installed SeqAn library." "")
option( BOOST_PATH "Path to manually installed Boost library." "")
option( STRICT_WARN "Flag for strict warnings." OFF )

# Set HiLive2 compile-time variables 
add_definitions(-DHiLive_VERSION_MAJOR=2)
add_definitions(-DHiLive_VERSION_MINOR=0)

###################
### Build types ###
###################

# Release flags
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" )

# Warning flags
if ( STRICT_WARN )
	add_compile_options( -Wall -Wextra -Wshadow -Wuninitialized -Wnon-virtual-dtor
		-Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic
		-Wconversion -Wsign-conversion -Wdouble-promotion
		-Wformat=2 -Wstrict-aliasing -Wno-long-long -Wno-variadic-macros )
endif()

# SeqAn specific flags
if ( NOT CONDA )
	add_compile_options( -static -march=native )
endif()

# Set flags for compilation
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -g -pthread -W -Wall -std=gnu++14 -O0")

# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")


############################################
### Dependencies and 3rd party libraries ###
############################################

# BOOST
if ( BOOST_PATH )
	set( Boost_NO_BOOST_CMAKE TRUE )
	set( Boost_NO_SYSTEM_PATHS TRUE )
	set( BOOST_ROOT "${BOOST_PATH}" )
	set( Boost_LIBRARY_DIRS "${BOOST_PATH}/lib" ) 
endif()
set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME ON) 
find_package( Boost COMPONENTS system filesystem program_options iostreams REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

# ZLIB
find_package( ZLIB REQUIRED )

# LZ4
if ( LZ4_PATH )
	set( LZ4_LIB "${LZ4_PATH}/lib" )
endif()

include_directories(${LZ4_LIB})
link_directories(${LZ4_LIB})
set(CompressionLibs "${ZLIB_LIBRARIES};lz4")

# SEQAN
if ( SEQAN_PATH )
	set (CMAKE_MODULE_PATH "${SEQAN_PATH}/util/cmake") 
	set (SEQAN_INCLUDE_PATH "${SEQAN_PATH}/include/") 
endif()

find_package (SeqAn REQUIRED)

# Add include directories, defines, and flags for SeqAn (and its dependencies).
add_library( SeqAn INTERFACE )
include_directories (${SEQAN_INCLUDE_DIRS})
add_compile_options( ${SEQAN_DEFINITIONS} )
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")


##############################
### setup HiLive libraries ###
##############################

include_directories("${PROJECT_SOURCE_DIR}/lib")

# make a list of HiLive libraries
set (LIB_NAMES tools_static tools alnread alnstream alnout illumina_parsers kindex parallel argument_parser)
set(LIB_LIST "")
foreach (x ${LIB_NAMES})
	list(APPEND LIB_LIST "lib/${x}.cpp")
endforeach()
add_library(HiLiveLibs ${LIB_LIST})

#############################
### Build the executables ###

add_executable (hilive tools/hilive.cpp)
target_link_libraries (hilive HiLiveLibs ${CompressionLibs} ${Boost_LIBRARIES} ${SEQAN_LIBRARIES})

add_executable(hilive-build tools/build_index.cpp )
target_link_libraries(hilive-build HiLiveLibs ${CompressionLibs} ${Boost_LIBRARIES} ${SEQAN_LIBRARIES})

add_executable(hilive-out tools/hilive_out.cpp )
target_link_libraries(hilive-out HiLiveLibs ${CompressionLibs} ${Boost_LIBRARIES} ${SEQAN_LIBRARIES})

