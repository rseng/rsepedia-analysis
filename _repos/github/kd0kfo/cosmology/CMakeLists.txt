cmake_minimum_required(VERSION 2.8)

project(cosmology CXX)
set(COSMOLOGY_MAJOR_VERSION 2)
set(COSMOLOGY_MINOR_VERSION 10)
set(COSMOLOGY_REVISION 2)
set(AUTHOR_EMAIL david@davecoss.com)

# Old and undeveloped. Perhaps continue some day?
#AH_TEMPLATE([WEB_USAGE],[Define to use out geared towards web usage])
option(EXCLUDE_FFT "Define to not build with FFT" OFF)

set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_C_COMPILER mpicc)

# Extra Packages
find_package(BISON)
find_package(FLEX)

#libdnstd
option(LIBDNSTD_DIR "Include libdnstd library" OFF)
if(LIBDNSTD_DIR)
   include_directories(${LIBDNSTD_DIR}/include)
   link_directories(${LIBDNSTD_DIR}/lib)
endif(LIBDNSTD_DIR)

#libmygl
option(LIBMYGL_DIR "Include libmygl library" OFF)
if(LIBMYGL_DIR)
   include_directories(${LIBMYGL_DIR}/include)
   link_directories(${LIBMYGL_DIR}/lib)
endif(LIBMYGL_DIR)

#fftw3
option(FFTW3_DIR "Include FFTW3" OFF)
if(FFTW3_DIR)
   include_directories(${FFTW3_DIR}/include)
   link_directories(${FFTW3_DIR}/lib)
   message("Build with FFTW: ${FFTW3_DIR}")
endif(FFTW3_DIR)

#web conditions
#AC_ARG_ENABLE([web],AS_HELP_STRING([--enable-web],[Enable output conditions for web use (e.g. cgi). Default: disabled]),
#AC_DEFINE([WEB_USAGE],[1])
#AC_MSG_RESULT(Targeting output for web usage),)

#mpi
set(USE_MPI 1)
option(MPI_DIR "Define to compile with a specific MPI path." OFF)
if(MPI_DIR)
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -I${MPI_DIR}/include)
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_CXX_FLAGS}-I${MPI_DIR}/lib)
endif(MPI_DIR)

# netcdf
option(NETCDF_DIR "Include netcdf" OFF)
if(NETCDF_DIR)
	execute_process(COMMAND ${NETCDF_DIR}/bin/nc-config --cflags OUTPUT_VARIABLE NETCDF_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
	execute_process(COMMAND ${NETCDF_DIR}/bin/nc-config --libs OUTPUT_VARIABLE NETCDF_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
else(NETCDF_DIR)
	execute_process(COMMAND nc-config --cflags OUTPUT_VARIABLE NETCDF_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
	execute_process(COMMAND nc-config --libs OUTPUT_VARIABLE NETCDF_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(NETCDF_DIR)

# config.h
configure_file( "${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_SOURCE_DIR}/config.h" )
add_definitions(-DHAVE_CONFIG_H -I${PROJECT_SOURCE_DIR})

# version.cpp
execute_process(COMMAND git rev-parse HEAD OUTPUT_FILE ${PROJECT_BINARY_DIR}/build_git_sha.out)
file(STRINGS ${PROJECT_BINARY_DIR}/build_git_sha.out BUILD_GIT_SHA)
execute_process(COMMAND date OUTPUT_FILE ${PROJECT_BINARY_DIR}/build_time.out)
file(STRINGS ${PROJECT_BINARY_DIR}/build_time.out BUILD_TIME)
configure_file("${PROJECT_SOURCE_DIR}/src/version.cpp.in" "${PROJECT_SOURCE_DIR}/src/version.cpp")

add_executable(ray_trace_ellipse ${PROJECT_SOURCE_DIR}/src/ray_trace_ellipse.cpp ${PROJECT_SOURCE_DIR}/src/version.cpp)
add_executable(flatten ${PROJECT_SOURCE_DIR}/src/flattenmain.cpp ${PROJECT_SOURCE_DIR}/src/flatten.cpp)

BISON_TARGET(physcalc_parser ${PROJECT_SOURCE_DIR}/src/physcalc.yacc.ypp ${PROJECT_SOURCE_DIR}/src/physcalc.yacc.cpp)
FLEX_TARGET(physcalc_scanner ${PROJECT_SOURCE_DIR}/src/physcalc.lex.lpp ${PROJECT_SOURCE_DIR}/src/physcalc.lex.cpp)
ADD_FLEX_BISON_DEPENDENCY(physcalc_scanner physcalc_parser)
add_executable(physcalc ${PROJECT_SOURCE_DIR}/src/physcalc.cpp  ${PROJECT_SOURCE_DIR}/src/functions.cpp ${PROJECT_SOURCE_DIR}/src/version.cpp ${BISON_physcalc_parser_OUTPUTS} ${FLEX_physcalc_scanner_OUTPUTS})
target_link_libraries(physcalc fftw3)

add_executable(makecluster ${PROJECT_SOURCE_DIR}/src/makecluster.cpp ${PROJECT_SOURCE_DIR}/src/create_cluster.cpp)

add_executable(plane2netcdf ${PROJECT_SOURCE_DIR}/src/plane2netcdf.cpp)

set(NETCDF_BINS ray_trace_ellipse flatten physcalc makecluster plane2netcdf)
set_target_properties(${NETCDF_BINS} PROPERTIES COMPILE_FLAGS "${NETCDF_CFLAGS}")
foreach(curr_bin ${NETCDF_BINS})
		 target_link_libraries(${curr_bin} mygl dnstd ${NETCDF_LIBS})
endforeach(curr_bin)

install(TARGETS ray_trace_ellipse flatten  physcalc makecluster plane2netcdf DESTINATION bin)

# tests
enable_testing()
add_test(ray_trace_test ${CMAKE_COMMAND} -DTEST_DIR=${PROJECT_SOURCE_DIR}/test -DEXE_DIR=${PROJECT_BINARY_DIR} -P ${CMAKE_SOURCE_DIR}/test/test_schwarzschild.cmake)

