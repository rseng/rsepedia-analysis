# Based on what I learnt from https://cliutils.gitlab.io/modern-cmake
cmake_minimum_required(VERSION 3.14...3.19)

# Create a release output by default
set(default_build_type "Release")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)


option(ALLOW_VERBOSE "Enable compile-time verbose mode" OFF)
if (ALLOW_VERBOSE OR ("${CMAKE_BUILD_TYPE}" STREQUAL Debug))
	set(ALLOW_VERBOSE ON)
	set(CMAKE_VERBOSE_MAKEFILE ON)
	message("Verbose mode enabled")
endif()


# Silence timestamp warning
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
	cmake_policy(SET CMP0135 NEW)
endif()

# Describe the project
project(udpPacketManager 	VERSION 0.9.0
							DESCRIPTION "C Library for handling CEP packet data from an international LOFAR station"
							LANGUAGES C CXX
)

set(UPM_CLI_VERSION 0.9)

message("Building udpPacketManager Version ${PROJECT_VERSION}\n\n")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native -O2 -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O2 -ffast-math")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wstrict-prototypes -Wmissing-prototypes -Wframe-larger-than=4095 -Wshadow -Wconversion -Wno-sign-conversion -Wall -Wextra -g3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wframe-larger-than=4095 -Wshadow -Wconversion -Wno-sign-conversion -Wall -Wextra -g3")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")# -fanalyzer")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")# -fanalyzer")


if("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
	if (${CMAKE_C_COMPILER_VERSION} VERSION_GREATER_EQUAL 13.0)
		set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fno-semantic-interposition")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-semantic-interposition")
	endif()
	set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fsanitize-trap=undefined")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize-trap=undefined")
else()
	# Only add GCC 10+ support for allow-store-data-races
	if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" AND ${CMAKE_C_COMPILER_VERSION} VERSION_GREATER_EQUAL 10.0)
		set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fallow-store-data-races")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fallow-store-data-races")
	endif()
	set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fno-semantic-interposition")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-semantic-interposition")
	set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wstack-usage=12287")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wstack-usage=12287")
endif()
# Ensure we have a modern python version available and prepare to install the latest build of dreamBeam
# multiprocessing.shared_memory requires >= 3.8
find_package(Python3 3.8 COMPONENTS Interpreter REQUIRED)

add_custom_command(OUTPUT pip_installed
                   COMMAND ${Python3_EXECUTABLE} -m pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt
                   COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/pip_installed
)
add_custom_target(install_python_requirements
                  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/pip_installed
)

# Find the OpenMP package on the system

message("Configuring OpenMP...")
find_package(OpenMP REQUIRED)

execute_process(COMMAND bash -c "cat /proc/cpuinfo | uniq | grep -m 2 \"siblings\" | cut -d \":\" -f 2 | sort --numeric --unique | xargs echo"
                OUTPUT_VARIABLE OMP_THREADS
)

#if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
if (OpenMP_C_LIB_NAMES MATCHES ".*(gomp).*")
	message(WARNING "The GCC OpenMP implementation has been detected (performance may be degraded, we strongly recommend avoiding the GCC OpenMP library, see compilers.md for further details)\nSetup will continue in 5 seconds.")
	set(OMP_THREADS 8)
	if (NOT("${CMAKE_BUILD_TYPE}" MATCHES ".*Debug.*"))
		execute_process(COMMAND ${CMAKE_COMMAND} -E sleep 5)
	endif()
endif()
message("")


# Install zstandard
message("Configuring Zstandard")
include(FetchContent)
set(ZSTD_BUILD_STATIC ON CACHE INTERNAL "")
set(ZSTD_BUILD_SHARED OFF CACHE INTERNAL "")
FetchContent_Declare(zstd
                     GIT_REPOSITORY https://github.com/facebook/zstd.git
                     GIT_TAG v1.5.5
)
FetchContent_MakeAvailable(zstd)
# MakeAvailable doesn't allow for libzstd_static to be used as a target
# Related? https://github.com/facebook/zstd/issues/2377
add_subdirectory(${zstd_SOURCE_DIR}/build/cmake/ ${zstd_BINARY_DIR} EXCLUDE_FROM_ALL)


#message("Configuring CSpice")
#FetchContent_Declare(cspice
#                     GIT_REPOSITORY https://github.com/tudat-team/cspice-cmake
#                     )
#FetchContent_MakeAvailable(cspice)

message("")

message("Configuring HDF5 libraries")

# zlib is "required" by hdf5. If it is present anywhere on the user's path it will be enabled and included, I cannot disable this for the life of me.
include(ExternalProject)
ExternalProject_Add(internal_zlib
                    GIT_REPOSITORY https://github.com/madler/zlib
                    GIT_TAG v1.2.12
                    CONFIGURE_COMMAND ./configure --static --prefix=./
                    BUILD_COMMAND make test install -j8
                    BUILD_IN_SOURCE TRUE # Keep source beside install prefix
                    INSTALL_COMMAND ""
                    UPDATE_COMMAND ""
)
ExternalProject_Get_Property(internal_zlib INSTALL_DIR)
set(zlib_INSTALL_DIR ${INSTALL_DIR})
add_library(libz STATIC IMPORTED)
set_property(TARGET libz PROPERTY IMPORTED_LOCATION ${zlib_INSTALL_DIR}/src/internal_zlib/libz.a)
add_dependencies(libz internal_zlib)
message("")
ExternalProject_Add(internal_hdf5
                    GIT_REPOSITORY https://github.com/HDFGroup/hdf5
                    GIT_TAG hdf5-1_12_2
                    CONFIGURE_COMMAND ./configure --disable-dependency-tracking --enable-optimization=high --disable-shared --disable-cxx --disable-hl --disable-hltools --disable-tools --enable-threadsafe --disable-fortran --enable-build-mode=production --with-zlib=${zlib_INSTALL_DIR}/src/internal_zlib/
                    BUILD_COMMAND make all -j8
                    BUILD_IN_SOURCE TRUE # Header gets moved to incorrect folder otherwise
                    INSTALL_COMMAND ""
                    UPDATE_COMMAND ""
                    DEPENDS internal_zlib
)
ExternalProject_Get_Property(internal_hdf5 INSTALL_DIR)
set(hdf5_INSTALL_DIR ${INSTALL_DIR})
add_library(libhdf5 STATIC IMPORTED)
set_property(TARGET libhdf5 PROPERTY IMPORTED_LOCATION ${hdf5_INSTALL_DIR}/src/internal_hdf5/src/.libs/libhdf5.a)
add_dependencies(libhdf5 internal_hdf5)
message("")

message("Configuring HDF5 plugins")
ExternalProject_Add(internal_bitshuffle
                    GIT_REPOSITORY https://github.com/kiyo-masui/bitshuffle
                    GIT_TAG fdfcd404ac8dcb828857a90c559d36d8ac4c2968
                    CONFIGURE_COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/src/misc/bitshuffleMakefile Makefile
                    BUILD_COMMAND CMAKE_BASE_DIR=${CMAKE_CURRENT_BINARY_DIR} make link -j4
                    BUILD_IN_SOURCE TRUE
                    INSTALL_COMMAND ""
                    UPDATE_COMMAND ""
                    DEPENDS internal_hdf5 zstd internal_zlib
)
ExternalProject_Get_Property(internal_bitshuffle INSTALL_DIR)
set(bitshuffle_INSTALL_DIR ${INSTALL_DIR})
add_library(libh5bshuf STATIC IMPORTED)
set_property(TARGET libh5bshuf PROPERTY IMPORTED_LOCATION ${bitshuffle_INSTALL_DIR}/src/internal_bitshuffle/libh5bshuf.a)
add_dependencies(libh5bshuf internal_bitshuffle)
message("")


message ("Configuring PSRDADA")
# Find the PSRDADA libraries on the system (or set NODADA if not found)
ExternalProject_Add(internal_PSRDADA
                    GIT_REPOSITORY git://git.code.sf.net/p/psrdada/code
                    GIT_TAG ba2b88
                    BUILD_IN_SOURCE TRUE
                    CONFIGURE_COMMAND ./bootstrap && ./configure --with-cuda-dir=no && cd 3rdparty && make libtimers.la # libtimers fails to compile during normal step with older Make versions
                    BUILD_COMMAND cd src/ && make libpsrdada.la dada_db dada_dbmeminfo dada_dbmetric dada_dbmonitor -j4
                    INSTALL_COMMAND ""
                    UPDATE_COMMAND ""
)
ExternalProject_Get_Property(internal_PSRDADA INSTALL_DIR)
set(psrdada_INSTALL_DIR ${INSTALL_DIR})
add_library(libpsrdada STATIC IMPORTED)
set_property(TARGET libpsrdada PROPERTY IMPORTED_LOCATION ${psrdada_INSTALL_DIR}/src/internal_PSRDADA/src/.libs/libpsrdada.a)
add_dependencies(libpsrdada internal_PSRDADA)
message("")


# Configure FFTW3 (w/ floats + omp support, optional)
message("Configuring FFTW3 (static, floats, OpenMP support enabled)")
## CMake config is outdated as of FFTW3 3.3.10, manually building to get AVX-512 support for REALTA
# I do not understand why, but the clean build process fails without the "echo $(pwd)" command.
ExternalProject_ADD(internal_FFTW3F
                    URL "https://www.fftw.org/fftw-3.3.10.tar.gz"
                    URL_HASH "MD5=8ccbf6a5ea78a16dbc3e1306e234cc5c"
                    CONFIGURE_COMMAND ./configure --disable-dependency-tracking --enable-float --enable-openmp --enable-static --enable-sse2 --enable-avx --enable-avx2 --enable-avx512 --enable-avx-128-fma
                    BUILD_COMMAND make -j4
                    BUILD_IN_SOURCE TRUE
                    INSTALL_COMMAND ""
                    UPDATE_COMMAND ""
)
ExternalProject_Get_Property(internal_FFTW3F install_dir)
set(fftw3f_INSTALL_DIR ${install_dir})
add_library(libfftw3f_omp STATIC IMPORTED)
set_property(TARGET libfftw3f_omp PROPERTY IMPORTED_LOCATION ${fftw3f_INSTALL_DIR}/src/internal_FFTW3F/threads/.libs/libfftw3f_omp.a)
add_dependencies(libfftw3f_omp internal_FFTW3F)
add_library(libfftw3f STATIC IMPORTED)
set_property(TARGET libfftw3f PROPERTY IMPORTED_LOCATION ${fftw3f_INSTALL_DIR}/src/internal_FFTW3F/.libs/libfftw3f.a)
add_dependencies(libfftw3f internal_FFTW3F)
message("")


#message("Configuring LibYAML")
#FetchContent_Declare(yaml
#                     GIT_REPOSITORY https://github.com/yaml/libyaml.git
#                     GIT_TAG 2c891fc7a770e8ba2fec34fc6b545c672beb37e6 # Git tag for 0.2.5 is erroring
#                     )
#FetchContent_MakeAvailable(yaml)


# Include compile-time parameters into the headers
configure_file(	"${CMAKE_CURRENT_SOURCE_DIR}/src/lib/lofar_udp_general.h.in"
				"${CMAKE_CURRENT_BINARY_DIR}/src/lib/lofar_udp_general.h"
				@ONLY
)

configure_file(	"${CMAKE_CURRENT_SOURCE_DIR}/src/CLI/lofar_cli_meta.h.in"
				"${CMAKE_CURRENT_BINARY_DIR}/src/CLI/lofar_cli_meta.h"
				@ONLY
)


# Setup the base library object
add_library(lofudpman STATIC
            src/lib/lofar_udp_general.c
            src/lib/lofar_udp_backends.cpp
            src/lib/lofar_udp_calibration.c
            src/lib/lofar_udp_io.c
            src/lib/lofar_udp_reader.c
            src/lib/lofar_udp_structs.c
            src/CLI/lofar_cli_meta.c
            src/lib/lofar_udp_metadata.c
            src/lib/lofar_udp_time.c)

add_dependencies(lofudpman libzstd_static libpsrdada libhdf5 libz libh5bshuf install_python_requirements) # libfftw3fomp) ##yaml) #CSpice::cspice)

# Include all of our library headers
target_include_directories(lofudpman PUBLIC
       ${CMAKE_CURRENT_SOURCE_DIR}/src/lib
       ${CMAKE_CURRENT_BINARY_DIR}/src/lib # Configured headers
       ${CMAKE_CURRENT_SOURCE_DIR}/src/CLI
       ${CMAKE_CURRENT_BINARY_DIR}/src/CLI # Configured headers
       ${CMAKE_CURRENT_SOURCE_DIR}/src/metadata

)
# Include all of our dependency headers
target_include_directories(lofudpman PUBLIC ${fftw3f_INSTALL_DIR}/src/internal_FFTW3F/api)
target_include_directories(lofudpman PUBLIC ${hdf5_INSTALL_DIR}/src/internal_hdf5/src/ ${zlib_INSTALL_DIR}/src/internal_zlib/ ${bitshuffle_INSTALL_DIR}/src/internal_bitshuffle/src/)
target_include_directories(lofudpman PUBLIC ${psrdada_INSTALL_DIR}/src/internal_PSRDADA/src/)
target_include_directories(lofudpman PUBLIC ${zstd_SOURCE_DIR}/lib/)
#target_include_directories(lofudpman PUBLIC ${yaml_SOURCE_DIR}/include/)

file(GLOB_RECURSE UDP_INCLUDE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.hpp" "${CMAKE_CURRENT_BINARY_DIR}/src/*/*.h" "${CMAKE_CURRENT_BINARY_DIR}/src/*/*.hpp")

# FFTW requires a libm (math) link
find_library(LIB_M m REQUIRED)
target_link_libraries(lofudpman PUBLIC libfftw3f_omp libfftw3f ${LIB_M})
target_link_libraries(lofudpman PUBLIC libhdf5 libz libh5bshuf ${CMAKE_DL_LIBS}) # Some cases need libdl for HDF5
target_link_libraries(lofudpman PUBLIC libpsrdada)
target_link_libraries(lofudpman PUBLIC libzstd_static)

find_library(LIB_RT rt REQUIRED) # Required for shmem functions
target_link_libraries(lofudpman PUBLIC ${LIB_RT})
#target_link_libraries(lofudpman PUBLIC yaml)

# Set a few extra compiler options
# Overwrite pre-installed libraries that we re-compile
set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "")
set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "")
set_property(TARGET lofudpman PROPERTY LINKER_LANGUAGE CXX)
set_property(TARGET lofudpman PROPERTY LINK_WHAT_YOU_USE ON)
#set_property(TARGET lofudpman PROPERTY INTERPROCEDURAL_OPTIMIZATION ON) # Static + IPO -> build failures?

# Link OpenMP
target_link_libraries(lofudpman PUBLIC OpenMP::OpenMP_CXX OpenMP::OpenMP_C)


# Setup the CLIs
add_executable(lofar_udp_extractor ${CMAKE_CURRENT_SOURCE_DIR}/src/CLI/lofar_cli_extractor.c)
add_executable(lofar_stokes_extractor ${CMAKE_CURRENT_SOURCE_DIR}/src/CLI/lofar_cli_stokes.c)
target_link_libraries(lofar_udp_extractor PUBLIC lofudpman)
target_link_libraries(lofar_stokes_extractor PUBLIC lofudpman)


include(CMakePackageConfigHelpers)
write_basic_package_version_file(
		lofudpmanConfigVersion.cmake
		VERSION ${PACKAGE_VERSION}
		COMPATIBILITY AnyNewerVersion
)

# Install everything
install(TARGETS lofudpman lofar_udp_extractor lofar_stokes_extractor
		EXPORT lofudpman
		LIBRARY DESTINATION lib
		RUNTIME DESTINATION bin
)
install(FILES ${UDP_INCLUDE_FILES} DESTINATION include)
install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/src/misc/dreamBeamJonesGenerator.py DESTINATION bin)

# Add the tests directory
enable_testing()
set(TEST_ENV_PATH ${CMAKE_CURRENT_SOURCE_DIR}/tests/lib_tests/)
add_subdirectory(tests)