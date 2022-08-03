cmake_minimum_required(VERSION 3.8)

project(wsclean)

include(ExternalProject)
include(CMakeVersionInfo.txt)

add_compile_options(-O3 -Wall -Wzero-as-null-pointer-constant)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Debug build selected: setting linking flag --no-undefined")
  string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
else()
  add_compile_options(-DNDEBUG)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS OFF)

option(PORTABLE "Generate portable code" OFF)
option(BUILD_PACKAGES "Build Debian packages" OFF)

include(CheckCXXCompilerFlag)
if(PORTABLE)
  if(DEFINED TARGET_CPU)
    message(WARNING "You have selected to build PORTABLE binaries. "
                    "TARGET_CPU settings will be ignored.")
    unset(TARGET_CPU CACHE)
  endif()
else()
  if(DEFINED TARGET_CPU)
    add_compile_options(-march=${TARGET_CPU})
  else()
    check_cxx_compiler_flag("-march=native" COMPILER_HAS_MARCH_NATIVE)
    if(COMPILER_HAS_MARCH_NATIVE)
      add_compile_options(-march=native)
    else()
      message(
        WARNING "The compiler doesn't support -march=native for your CPU.")
    endif()
  endif()
endif()

set(ExternalSubmoduleDirectories aocommon pybind11 radler schaapcommon)
foreach(ExternalSubmodule ${ExternalSubmoduleDirectories})
  if(NOT EXISTS ${CMAKE_SOURCE_DIR}/external/${ExternalSubmodule})
    message(
      FATAL_ERROR
        "The external submodule '${ExternalSubmodule}' is missing in the external/ subdirectory. "
        "This is likely the result of downloading a git tarball without submodules. "
        "This is not supported: git tarballs do not provide the required versioning "
        "information for the submodules. Please perform a git clone of this repository."
    )
  endif()
endforeach()

# Find and include git submodules
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  # Update submodules as needed
  option(GIT_SUBMODULE "Check submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --checkout
              --depth 1
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(
        FATAL_ERROR
          "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules"
      )
    endif()
  endif()
endif()

# Include aocommon
include_directories("${CMAKE_SOURCE_DIR}/external/aocommon/include/")

# Schaapcommon
add_subdirectory("${CMAKE_SOURCE_DIR}/external/schaapcommon")
include_directories("${CMAKE_SOURCE_DIR}/external/schaapcommon/include")

include_directories(external/wgridder)

# Casacore has a separate CMake file in this directory
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

find_package(
  HDF5
  COMPONENTS C CXX
  REQUIRED)
add_definitions(${HDF5_DEFINITIONS})

set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits)
find_package(CFITSIO REQUIRED)
# Chgcentre uses a lapack function. EveryBeam provides lapack, but because
# EveryBeam is not required, chgcentre will fail to link when EveryBeam is not
# available, if LAPACK would not be explicitly required.
find_package(LAPACK REQUIRED)

find_library(FFTW3_LIB fftw3 REQUIRED HINTS ENV FFTW3_LIB)
find_library(FFTW3_THREADS_LIB fftw3_threads REQUIRED HINTS ENV FFTW3_LIB)
find_library(FFTW3F_LIB fftw3f REQUIRED HINTS ENV FFTW3_LIB)
find_library(FFTW3F_THREADS_LIB fftw3f_threads REQUIRED HINTS ENV FFTW3_LIB)
find_path(
  FFTW3_INCLUDE_DIR
  NAMES fftw3.h
  HINTS ENV FFTW3_INCLUDE)

find_package(PythonLibs 3 REQUIRED)
find_package(PythonInterp REQUIRED)
message(STATUS "Using python version ${PYTHON_VERSION_STRING}")

# Include pybind11
add_subdirectory("${CMAKE_SOURCE_DIR}/external/pybind11")
include_directories(SYSTEM ${pybind11_INCLUDE_DIRS})

#Prevent accidentally finding old BoostConfig.cmake file from casapy
set(Boost_NO_BOOST_CMAKE ON)
find_package(
  Boost
  COMPONENTS date_time filesystem system program_options
  REQUIRED)
find_library(PTHREAD_LIB pthread REQUIRED)
find_library(DL_LIB dl REQUIRED)

# Once we bump minimum CMake version to >= 3.2, we can use `find_package(GSL REQUIRED)`.
find_library(GSL_LIB NAMES gsl)
find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h)
find_library(GSL_CBLAS_LIB NAMES gslcblas)
if(NOT GSL_LIB
   OR NOT GSL_INCLUDE_DIR
   OR NOT GSL_CBLAS_LIB)
  message(FATAL_ERROR "GSL not found, but required to build WSClean!")
endif()

find_package(MPI)

if(MPI_FOUND)
  # FindMPI in CMake >= 3.10 provides MPI_CXX_COMPILE_OPTIONS, a list that can
  # be fed into add_compile_options(). In older versions of CMake, FindMPI
  # provides MPI_CXX_COMPILE_FLAGS, a white-space separated string that can be
  # fed into add_definitions().
  if(DEFINED MPI_CXX_COMPILE_OPTIONS)
    add_compile_options(${MPI_CXX_COMPILE_OPTIONS})
  else()
    add_definitions(${MPI_CXX_COMPILE_FLAGS})
  endif()
  add_definitions(-DHAVE_MPI)
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  set(MPI_CPP_FILES scheduling/mpischeduler.cpp distributed/mpibig.cpp)
else(MPI_FOUND)
  message(
    WARNING
      "MPI not found, multi-processing executable wsclean-mp will not be build."
  )
  set(MPI_LIBRARIES "")
  set(MPI_CPP_FILES)
endif(MPI_FOUND)

find_package(EveryBeam NO_MODULE QUIET)
if(${EVERYBEAM_FOUND})
  if(${EVERYBEAM_VERSION} VERSION_LESS "0.3.0" OR ${EVERYBEAM_VERSION}
                                                  VERSION_GREATER_EQUAL "0.4.0")
    message(STATUS "EveryBeam version is incompatible")
    message(
      FATAL_ERROR
        "WSClean needs EveryBeam version 0.3.x - with x >= 0 - but found version ${EVERYBEAM_VERSION}"
    )
  else()
    include_directories(${EVERYBEAM_INCLUDE_DIRS})
    add_definitions(-DHAVE_EVERYBEAM)
    message(STATUS "EveryBeam library found.")
  endif()
else(${EVERYBEAM_FOUND})
  message(
    STATUS
      "EveryBeam library not found: Functionality for calculating beams for telescopes such as LOFAR will not be available."
  )
  set(EVERYBEAM_LIB "")
endif(${EVERYBEAM_FOUND})

find_package(IDGAPI NO_MODULE QUIET)

if(IDGAPI_FOUND)
  if(${IDGAPI_VERSION} VERSION_LESS "0.8.1")
    message(STATUS "IDG version is incompatible")
    message(
      FATAL_ERROR
        "WSClean needs IDG version 0.8.x - with x >= 1 - but found version ${IDGAPI_VERSION}"
    )
  else()
    set(IDG_FILES idg/idgmsgridder.cpp)
    include_directories(${IDGAPI_INCLUDE_DIRS})
    add_definitions(-DHAVE_IDG)
    message(STATUS "Image domain gridder libraries found.")
  endif()
else()
  set(IDGAPI_LIBRARIES)
  set(IDG_FILES)
  message(
    STATUS
      "Image domain gridder libraries NOT found. This gridder will not be available."
  )
endif()

include_directories(SYSTEM ${CASACORE_INCLUDE_DIRS})
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
include_directories(SYSTEM ${CFITSIO_INCLUDE_DIR})
include_directories(SYSTEM ${FFTW3_INCLUDE_DIR})
include_directories(SYSTEM ${GSL_INCLUDE_DIR})
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})

include(CheckCXXSourceCompiles)

# GSL is required for WSClean, so always available
add_definitions(-DHAVE_GSL)

# Radler
set(RADLER_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/external/radler)
set(RADLER_LIB ${RADLER_INSTALL_PREFIX}/lib/libradler.so)
ExternalProject_Add(
  radler-module
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/radler
  BINARY_DIR ${CMAKE_BINARY_DIR}/external/radler/build
  BUILD_ALWAYS TRUE
  BUILD_BYPRODUCTS ${RADLER_LIB}
  CMAKE_ARGS
    -DCOMPILE_AS_EXTERNAL_PROJECT=On
    -DPORTABLE=${PORTABLE}
    -DCMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH}
    -DCMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH}
    -DCFITSIO_ROOT_DIR=${CFITSIO_ROOT_DIR}
    -DCASACORE_ROOT_DIR=${CASACORE_ROOT_DIR}
    -DAOCOMMON_INCLUDE_DIR=${CMAKE_SOURCE_DIR}/external/aocommon/include/
    -DSCHAAPCOMMON_SOURCE_DIR=${CMAKE_SOURCE_DIR}/external/schaapcommon/
    -DPYBIND11_SOURCE_DIR=${CMAKE_SOURCE_DIR}/external/pybind11
    -DBUILD_TESTING=OFF
    -DCMAKE_INSTALL_PREFIX=${RADLER_INSTALL_PREFIX})

set(RADLER_INCLUDE_DIR ${RADLER_INSTALL_PREFIX}/include)
include_directories(${RADLER_INCLUDE_DIR})

# The following section will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# Include GNUInstallDirs for CMAKE_INSTALL_FULL_LIBDIR
include(GNUInstallDirs)
# Use, i.e. don't skip the full RPATH for the build tree.
set(CMAKE_SKIP_BUILD_RPATH FALSE)
# When building, don't use the install RPATH already
# (but later on when installing).
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# The RPATH to be used when installing, but only if it's not a system directory.
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
     "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
endif("${isSystemDir}" STREQUAL "-1")

configure_file("${PROJECT_SOURCE_DIR}/wscversion.h.in"
               "${PROJECT_BINARY_DIR}/wscversion.h")
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add CPack directory if user wants to generate Debian packages
if(BUILD_PACKAGES)
  add_subdirectory(CPack)
endif()

add_library(
  wsclean-object OBJECT
  main/commandline.cpp
  main/progressbar.cpp
  main/stopwatch.cpp
  main/wsclean.cpp
  main/settings.cpp
  gridding/directmsgridder.cpp
  gridding/msgridderbase.cpp
  gridding/wsmsgridder.cpp
  gridding/wstackinggridder.cpp
  idg/averagebeam.cpp
  interface/wscleaninterface.cpp
  io/componentlistwriter.cpp
  io/facetreader.cpp
  io/parsetreader.cpp
  io/wscfitswriter.cpp
  math/imageoperations.cpp
  math/renderer.cpp
  model/model.cpp
  msproviders/averagingmsrowprovider.cpp
  msproviders/mscolumns.cpp
  msproviders/bdamsrowprovider.cpp
  msproviders/contiguousms.cpp
  msproviders/msdatadescription.cpp
  msproviders/directmsrowprovider.cpp
  msproviders/msprovider.cpp
  msproviders/msrowprovider.cpp
  msproviders/msrowproviderbase.cpp
  msproviders/noisemsrowprovider.cpp
  msproviders/partitionedms.cpp
  msproviders/timestepbuffer.cpp
  msproviders/synchronizedms.cpp
  msproviders/msreaders/contiguousmsreader.cpp
  msproviders/msreaders/partitionedmsreader.cpp
  msproviders/msreaders/timestepbufferreader.cpp
  msproviders/msweightcolumn.cpp
  scheduling/griddingtaskmanager.cpp
  scheduling/griddingresult.cpp
  scheduling/griddingtask.cpp
  scheduling/metadatacache.cpp
  scheduling/threadedscheduler.cpp
  structures/facetutil.cpp
  structures/imageweights.cpp
  structures/imagingtable.cpp
  structures/imagingtableentry.cpp
  structures/msselection.cpp
  structures/observationinfo.cpp
  structures/primarybeam.cpp
  structures/resources.cpp
  system/pythonfilepath.cpp
  wgridder/wgriddingmsgridder.cpp
  wgridder/wgriddinggridder_simple.cpp
  external/wgridder/ducc0/infra/threading.cc
  ${IDG_FILES}
  ${MPI_CPP_FILES})

# A number of files perform the 'core' high-performance floating point
# operations. In these files, NaNs are avoided and thus -ffast-math is
# allowed. Note that visibilities can be NaN hence this can not be turned
# on for all files.
set_source_files_properties(
  wsclean/wstackinggridder.cpp wgridder/wgriddinggridder_simple.cpp
  PROPERTIES COMPILE_FLAGS -ffast-math)

# Add radler-module dependency to wsclean-object to ensure all
# wsclean targets can link against ${RADLER_LIB}
add_dependencies(wsclean-object radler-module)

set_property(TARGET wsclean-object PROPERTY POSITION_INDEPENDENT_CODE 1)

set(WSCLEANFILES $<TARGET_OBJECTS:wsclean-object>)

set(ALL_LIBRARIES
    ${CASACORE_LIBRARIES}
    ${FFTW3_LIB}
    ${FFTW3_THREADS_LIB}
    ${FFTW3F_LIB}
    ${FFTW3F_THREADS_LIB}
    ${Boost_DATE_TIME_LIBRARY}
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${CFITSIO_LIBRARY}
    ${GSL_LIB}
    ${GSL_CBLAS_LIB}
    ${PTHREAD_LIB}
    ${IDGAPI_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${MPI_LIBRARIES}
    ${EVERYBEAM_LIB}
    schaapcommon)

# Perform the BLAS check from aocommon
include("${CMAKE_CURRENT_LIST_DIR}/external/aocommon/CMake/CheckBLAS.cmake")
check_blas(LIBRARIES ${ALL_LIBRARIES})
# Radler lib cannot be checked on BLAS consistenty as it doesn't exist
# at cmake time.
list(APPEND ALL_LIBRARIES ${RADLER_LIB})

add_library(wsclean-lib STATIC ${WSCLEANFILES})
target_link_libraries(wsclean-lib PRIVATE pybind11::embed ${ALL_LIBRARIES})
set_target_properties(wsclean-lib PROPERTIES OUTPUT_NAME wsclean)
set_target_properties(wsclean-lib PROPERTIES SOVERSION ${WSCLEAN_VERSION_SO})

add_library(wsclean-shared SHARED ${WSCLEANFILES})
target_link_libraries(wsclean-shared PRIVATE pybind11::embed ${ALL_LIBRARIES})
set_target_properties(wsclean-shared PROPERTIES SOVERSION ${WSCLEAN_VERSION_SO})

add_executable(wsclean main/main.cpp)
target_link_libraries(wsclean wsclean-lib)

if(MPI_FOUND)
  add_executable(wsclean-mp distributed/wsclean-mp.cpp distributed/slave.cpp)
  target_link_libraries(wsclean-mp wsclean-lib)
endif(MPI_FOUND)

add_executable(chgcentre chgcentre/main.cpp chgcentre/progressbar.cpp)
target_link_libraries(chgcentre ${CASACORE_LIBRARIES} ${GSL_LIB}
                      ${GSL_CBLAS_LIB} ${MPI_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(wsuvbinning EXCLUDE_FROM_ALL gridding/examples/wsuvbinning.cpp
                                            ${WSCLEANFILES})
target_link_libraries(wsuvbinning ${ALL_LIBRARIES})

add_executable(
  wspredictionexample EXCLUDE_FROM_ALL gridding/examples/wspredictionexample.cpp
                                       gridding/wstackinggridder.cpp)
target_link_libraries(wspredictionexample ${ALL_LIBRARIES})

install(TARGETS wsclean DESTINATION bin)
install(TARGETS wsclean-lib DESTINATION lib)
install(TARGETS chgcentre DESTINATION bin)
install(FILES interface/wscleaninterface.h DESTINATION include)
# Installing the radler library to ${CMAKE_INSTALL_FULL_LIBDIR}
# ensures that users do not manually have to set the LD_LIBRARY_PATH
install(FILES ${RADLER_LIB} DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR})

if(MPI_FOUND)
  install(TARGETS wsclean-mp DESTINATION bin)
endif(MPI_FOUND)

add_custom_target(
  sphinxdoc
  make html
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
  COMMENT "Generating documentation with Sphinx"
  VERBATIM)

# add target to generate API documentation with Doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                 ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(
    doxygendoc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
  add_custom_target(doc DEPENDS doxygendoc sphinxdoc)
else(DOXYGEN_FOUND)
  message(STATUS "Doxygen not found: API documentation can not compiled.")
  add_custom_target(doc DEPENDS sphinxdoc)
endif(DOXYGEN_FOUND)

if(BUILD_TESTING)
  # Boost 1.59 introduced BOOST_TEST, which several tests use.
  find_package(Boost 1.59.0 COMPONENTS unit_test_framework)
  if(Boost_FOUND)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    add_subdirectory(tests)
  else()
    message(
      STATUS
        "Boost testing framework not found (not required for wsclean: only required for running tests)."
    )
  endif()
endif()

get_directory_property(MAIN_COMPILE_OPTIONS COMPILE_OPTIONS)
string(REPLACE ";" " " MAIN_COMPILE_OPTIONS "${MAIN_COMPILE_OPTIONS}")
message(STATUS "Flags passed to C++ compiler: ${MAIN_COMPILE_OPTIONS}")
