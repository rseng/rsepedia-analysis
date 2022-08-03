cmake_minimum_required(VERSION 3.12)

project(AOFlagger)

include(CTest)
include(ExternalProject)
include(CMakeVersionInfo.txt)
set(AOFLAGGER_VERSION
    "${AOFLAGGER_VERSION_MAJOR}.${AOFLAGGER_VERSION_MINOR}.${AOFLAGGER_VERSION_SUBMINOR}"
)

# Some CMake 'find package' files are stored in this directory
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(PkgConfig)
pkg_check_modules(SIGCXX sigc++-2.0)

option(ENABLE_GUI "Build rfigui and aoqplot tools" ON)
if(ENABLE_GUI)
  pkg_check_modules(GTKMM gtkmm-3.0>=3.0.0)
  if(GTKMM_FOUND)
    add_definitions(-DHAVE_GTKMM)
  else()
    message(
      WARNING
        "The graphical user interface library GTKMM was not found; rfigui and aoqplot will not be compiled."
    )
  endif(GTKMM_FOUND)
endif()

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

# Python needs to be found before pybind11 to avoid finding different versions
# of python
find_package(
  Python
  COMPONENTS Interpreter Development
  REQUIRED)
message(STATUS "Using python version ${Python_VERSION}")
include_directories(SYSTEM ${Python_INCLUDE_DIRS})

list(APPEND ExternalSubmoduleDirectories aocommon pybind11)
foreach(ExternalSubmodule IN LISTS ExternalSubmoduleDirectories)
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

# Include aocommon/pybind11 headers
include_directories("${CMAKE_SOURCE_DIR}/external/aocommon/include")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/pybind11")
include_directories(SYSTEM ${pybind11_INCLUDE_DIR})

find_package(
  HDF5
  COMPONENTS C CXX
  REQUIRED)
add_definitions(${HDF5_DEFINITIONS})
include_directories(SYSTEM ${HDF5_INCLUDE_DIR})

set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures)
include_directories(SYSTEM ${CASACORE_INCLUDE_DIRS})

find_package(CFITSIO REQUIRED)
include_directories(SYSTEM ${CFITSIO_INCLUDE_DIR})

find_path(
  FFTW3_INCLUDE_DIR
  NAMES fftw3.h
  HINTS ENV FFTW3_INCLUDE)

find_library(GSL_LIB NAMES gsl)
find_library(GSL_CBLAS_LIB NAMES gslcblas)
find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h)
if(GSL_INCLUDE_DIR)
  include_directories(SYSTEM ${GSL_INCLUDE_DIR})
endif(GSL_INCLUDE_DIR)

find_package(Threads REQUIRED)

find_package(PNG REQUIRED)

# boost::alignment requires Boost 1.56
find_package(Boost 1.56.0 REQUIRED COMPONENTS date_time filesystem system)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

find_library(FFTW3_LIB fftw3 REQUIRED HINTS ENV FFTW3_LIB)
include_directories(SYSTEM ${FFTW3_INCLUDE_DIR})

enable_language(Fortran OPTIONAL)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Lua 5.3 REQUIRED)
include_directories(SYSTEM ${LUA_INCLUDE_DIR})

# The following stuff will set the "rpath" correctly, so that LD_LIBRARY_PATH
# doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)
# when building, don't use the install RPATH already (but later on when
# installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH which point to directories
# outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
     "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)
set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wall -Wvla -Wzero-as-null-pointer-constant -DNDEBUG -O3"
)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # GCC 8.x requires linking with stdc++fs for the filesystem library
  # https://gcc.gnu.org/onlinedocs/gcc-9.1.0/libstdc++/manual/manual/status.html#status.iso.2017
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    link_libraries(stdc++fs)
  elseif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0)
    message(
      FATAL_ERROR "The GCC version is too old, upgrade to GCC 8.0 or newer")
  endif()
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Debug build selected: setting linking flag --no-undefined")
  string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
endif()

include(CheckCXXCompilerFlag)
include(CheckSymbolExists)
include(CheckCXXSymbolExists)
include(CheckIncludeFileCXX)

check_symbol_exists(posix_fallocate "fcntl.h" HAVE_POSIX_FALLOCATE)
check_cxx_symbol_exists(exp10 "cmath" HAVE_EXP10)
if(HAVE_POSIX_FALLOCATE)
  add_definitions(-DHAVE_POSIX_FALLOCATE)
endif(HAVE_POSIX_FALLOCATE)

option(PORTABLE "Build portable binaries (with slightly decreased performance)"
       ON)
if(PORTABLE)
  message(
    STATUS
      "Portable build requested; a generic build will be created with slightly decreased performance."
  )
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif(PORTABLE)

if(GSL_LIB
   AND GSL_CBLAS_LIB
   AND GSL_INCLUDE_DIR)
  add_definitions(-DHAVE_GSL)
  link_libraries(${GSL_LIB} ${GSL_CBLAS_LIB})
  message(STATUS "GSL found.")
endif(
  GSL_LIB
  AND GSL_CBLAS_LIB
  AND GSL_INCLUDE_DIR)

configure_file(version.h.in version.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(AOQPLOT_FILES
    aoqplot/aoqplotwindow.cpp
    aoqplot/baselineplotpage.cpp
    aoqplot/datawindow.cpp
    aoqplot/grayscaleplotpage.cpp
    aoqplot/histogrampage.cpp
    aoqplot/timefrequencyplotpage.cpp
    aoqplot/twodimensionalplotpage.cpp
    aoqplot/controllers/aoqplotcontroller.cpp
    aoqplot/controllers/aoqplotpagecontroller.cpp
    aoqplot/controllers/baselinepagecontroller.cpp
    aoqplot/controllers/heatmappagecontroller.cpp
    aoqplot/controllers/histogrampagecontroller.cpp
    aoqplot/controllers/tfpagecontroller.cpp)

set(PLOT_FILES
    rfigui/maskedheatmap.cpp
    plot/colorscale.cpp
    plot/heatmap.cpp
    plot/horizontalplotscale.cpp
    plot/legend.cpp
    plot/plotbase.cpp
    plot/plotpropertieswindow.cpp
    plot/plotwidget.cpp
    plot/title.cpp
    plot/vectorimage.cpp
    plot/verticalplotscale.cpp
    plot/xyplot.cpp)

set(GUI_FILES
    rfigui/controllers/imagecomparisoncontroller.cpp
    rfigui/controllers/rfiguicontroller.cpp
    rfigui/complexplaneplotwindow.cpp
    rfigui/gotowindow.cpp
    rfigui/imagepropertieswindow.cpp
    rfigui/opendialog.cpp
    rfigui/plotframe.cpp
    rfigui/plotwindow.cpp
    rfigui/progresswindow.cpp
    rfigui/rfiguimenu.cpp
    rfigui/rfiguiwindow.cpp
    rfigui/settings.cpp
    rfigui/simulatedialog.cpp
    rfigui/strategyeditor.cpp
    util/rfiplots.cpp
    util/multiplot.cpp
    ${PLOT_FILES}
    ${AOQPLOT_FILES})

set(IMAGING_FILES imaging/uvimager.cpp imaging/model.cpp)

set(INTERFACE_FILES
    interface/aoflagger.cpp interface/flagmask.cpp interface/imageset.cpp
    interface/qualitystatistics.cpp interface/strategy.cpp)

set(LUA_FILES
    lua/datawrapper.cpp
    lua/default-strategy.cpp
    lua/functions.cpp
    lua/functionswrapper.cpp
    lua/luastrategy.cpp
    lua/optionsfunction.cpp
    lua/scriptdata.cpp
    lua/telescopefile.cpp)

set(MSIO_FILES
    msio/baselinematrixloader.cpp
    msio/baselinereader.cpp
    msio/directbaselinereader.cpp
    msio/fitsfile.cpp
    msio/memorybaselinereader.cpp
    msio/msstatreader.cpp
    msio/pngfile.cpp
    msio/reorderingbaselinereader.cpp
    msio/rspreader.cpp
    msio/singlebaselinefile.cpp
    msio/spatialtimeloader.cpp)

set(STRUCTURES_FILES
    structures/image2d.cpp
    structures/mask2d.cpp
    structures/msiterator.cpp
    structures/msmetadata.cpp
    structures/samplerow.cpp
    structures/stokesimager.cpp
    structures/timefrequencydata.cpp)

set(QUALITY_FILES
    quality/collector.cpp quality/histogramcollection.cpp
    quality/histogramtablesformatter.cpp quality/qualitytablesformatter.cpp
    quality/rayleighfitter.cpp quality/statisticscollection.cpp)

set(ALGORITHMS_FILES
    algorithms/antennaselector.cpp
    algorithms/baselineselector.cpp
    algorithms/baselinetimeplaneimager.cpp
    algorithms/combinatorialthresholder.cpp
    algorithms/fringestoppingfitter.cpp
    algorithms/fringetestcreater.cpp
    algorithms/highpassfilter.cpp
    algorithms/localfitmethod.cpp
    algorithms/morphology.cpp
    algorithms/sinusfitter.cpp
    algorithms/siroperator.cpp
    algorithms/morphologicalflagger.cpp
    algorithms/sumthreshold.cpp
    algorithms/sumthresholdmissing.cpp
    algorithms/svdmitigater.cpp
    algorithms/testsetgenerator.cpp
    algorithms/thresholdconfig.cpp
    algorithms/thresholdtools.cpp
    algorithms/timefrequencystatistics.cpp)

set(PYTHON_FILES python/pythonstrategy.cpp python/pyfunctions.cpp)

set(IMAGESETS_FILES
    imagesets/bhfitsimageset.cpp
    imagesets/filterbankset.cpp
    imagesets/fitsimageset.cpp
    imagesets/h5imageset.cpp
    imagesets/imageset.cpp
    imagesets/indexableset.cpp
    imagesets/msimageset.cpp
    imagesets/multibandmsimageset.cpp
    imagesets/parmimageset.cpp
    imagesets/pngreader.cpp
    imagesets/rfibaselineset.cpp
    imagesets/sdhdfimageset.cpp)

set(UTIL_FILES
    plot/colormap.cpp
    util/logger.cpp
    util/ffttools.cpp
    util/integerdomain.cpp
    util/plot.cpp
    util/rng.cpp
    util/stopwatch.cpp)

set(ALL_LIBRARIES
    ${CASACORE_LIBRARIES}
    ${FFTW3_LIB}
    ${CFITSIO_LIBRARY}
    ${Python_LIBRARIES}
    ${PNG_LIBRARIES}
    ${LUA_LIBRARIES}
    ${Boost_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_C_LIBRARIES})

if(GTKMM_FOUND)
  message(STATUS "GTKMM found.")
  include_directories(SYSTEM ${GTKMM_INCLUDE_DIRS})
  link_directories(${GTKMM_LIBDIR})
  set(ALL_LIBRARIES ${ALL_LIBRARIES} ${GTKMM_LIBRARIES})
endif(GTKMM_FOUND)

if(SIGCXX_FOUND)
  message(STATUS "SIGCXX found.")
  include_directories(SYSTEM ${SIGCXX_INCLUDE_DIRS})
  link_directories(${SIGCXX_LIBDIR})
  set(ALL_LIBRARIES ${ALL_LIBRARIES} ${SIGCXX_LIBRARIES})
endif(SIGCXX_FOUND)

if(GTKMM_FOUND)
  set(AOFLAGGER_PLOT_FILES ${PLOT_FILES})
endif(GTKMM_FOUND)

set(ALL_NON_GUI_FILES
    ${ALGORITHMS_FILES}
    ${AOFLAGGER_PLOT_FILES}
    ${IMAGESETS_FILES}
    ${IMAGING_FILES}
    ${INTERFACE_FILES}
    ${LUA_FILES}
    ${MSIO_FILES}
    ${PYTHON_FILES}
    ${QUALITY_FILES}
    ${STRUCTURES_FILES}
    ${UTIL_FILES})
add_library(aoflagger-lib SHARED ${ALL_NON_GUI_FILES})
set_target_properties(aoflagger-lib PROPERTIES SOVERSION 0)
set_target_properties(aoflagger-lib PROPERTIES OUTPUT_NAME aoflagger)
target_link_libraries(aoflagger-lib ${ALL_LIBRARIES})
install(TARGETS aoflagger-lib DESTINATION lib)

add_executable(aoquality applications/aoquality.cpp)
target_link_libraries(aoquality aoflagger-lib ${ALL_LIBRARIES})
install(TARGETS aoquality DESTINATION bin)

if(GTKMM_FOUND)
  add_library(aoflaggergui OBJECT ${GUI_FILES})
  set(AOFLAGGERGUI_OBJECT $<TARGET_OBJECTS:aoflaggergui>)
  add_executable(rfigui applications/rfigui.cpp ${AOFLAGGERGUI_OBJECT})
  target_link_libraries(rfigui aoflagger-lib ${ALL_LIBRARIES})
  install(TARGETS rfigui DESTINATION bin)

  add_executable(aoqplot applications/aoqplot.cpp ${AOFLAGGERGUI_OBJECT})
  target_link_libraries(aoqplot aoflagger-lib ${ALL_LIBRARIES})
  install(TARGETS aoqplot DESTINATION bin)

  add_executable(badstations applications/badstations.cpp)
  target_link_libraries(badstations aoflagger-lib ${ALL_LIBRARIES})
endif(GTKMM_FOUND)

add_executable(
  aoflagger-bin
  applications/aoflagger.cpp aoluarunner/options.cpp aoluarunner/runner.cpp
  aoluarunner/baselineiterator.cpp aoluarunner/writethread.cpp)
set_target_properties(aoflagger-bin PROPERTIES OUTPUT_NAME aoflagger)
target_link_libraries(aoflagger-bin aoflagger-lib ${ALL_LIBRARIES})
install(TARGETS aoflagger-bin DESTINATION bin)

# A number of files perform the 'core' high-performance floating point
# operations. In these files, NaNs are avoided and thus -ffast-math is allowed.
# Note that visibilities can be NaN hence this can not be turned on for all
# files.
set_source_files_properties(
  algorithms/sumthreshold.cpp algorithms/sumthresholdmissing.cpp
  PROPERTIES COMPILE_FLAGS "-ffast-math -funroll-loops")

install(FILES interface/aoflagger.h DESTINATION include)
install(DIRECTORY data/icons data/applications DESTINATION share)
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data/strategies
  DESTINATION share/aoflagger
  FILES_MATCHING
  PATTERN *.lua)

# Set-up cmake configuration files
configure_file("${PROJECT_SOURCE_DIR}/cmake/config/aoflagger-config.cmake.in"
               "${PROJECT_BINARY_DIR}/CMakeFiles/aoflagger-config.cmake" @ONLY)
configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/config/aoflagger-config-version.cmake.in"
  "${PROJECT_BINARY_DIR}/CMakeFiles/aoflagger-config-version.cmake" @ONLY)

install(FILES "${PROJECT_BINARY_DIR}/CMakeFiles/aoflagger-config.cmake"
              "${PROJECT_BINARY_DIR}/CMakeFiles/aoflagger-config-version.cmake"
        DESTINATION share/aoflagger/cmake)

add_subdirectory(python)

add_subdirectory(cpack)

add_custom_target(
  sphinxdoc
  make html
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
  COMMENT "Generating documentation with Sphinx"
  VERBATIM)

find_package(Doxygen QUIET)

if(DOXYGEN_FOUND)
  # add target to generate API documentation with Doxygen
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                 ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  message(STATUS "Doxygen found: API documentation can be compiled.")
  add_custom_target(
    doxygendoc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
  add_custom_target(doc DEPENDS doxygendoc sphinxdoc)
else()
  message(STATUS "Doxygen not found: API documentation can not compiled.")
  add_custom_target(doc DEPENDS sphinxdoc)
endif()

# This is just for development: it is to convert the default Lua stategy into a
# c++ array.
add_custom_target(
  strategy
  xxd -i data/strategies/generic-default.lua >
  ${CMAKE_CURRENT_SOURCE_DIR}/lua/default-strategy.cpp
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Converting default strategy to C code"
  VERBATIM)

# TODO(RAP-368) Make a proper integration test target.
option(ENABLE_TESTS
       "Enable the integration tests, this requires downloading test data" OFF)
if(ENABLE_TESTS)

  find_package(Boost 1.56.0 REQUIRED COMPONENTS unit_test_framework)
  add_executable(
    runtests
    test/runtests.cpp
    test/experiments/defaultstrategyspeedtest.cpp
    test/experiments/highpassfilterexperiment.cpp
    test/experiments/tthroughput.cpp
    test/lua/defaultstrategytest.cpp
    test/lua/flagnanstest.cpp
    test/lua/tmetadata.cpp
    test/lua/tscript.cpp
    test/lua/optionsfunctiontest.cpp
    test/lua/telescopefiletest.cpp
    test/interface/interfacetest.cpp
    test/quality/qualitytablesformattertest.cpp
    test/quality/statisticscollectiontest.cpp
    test/quality/statisticsderivatortest.cpp
    test/algorithms/convolutionstest.cpp
    test/algorithms/dilationtest.cpp
    test/algorithms/highpassfiltertest.cpp
    test/algorithms/medianwindow.cpp
    test/algorithms/noisestatisticstest.cpp
    test/algorithms/siroperatortest.cpp
    test/algorithms/sumthresholdmissingtest.cpp
    test/algorithms/sumthresholdtest.cpp
    test/algorithms/testtools.cpp
    test/algorithms/thresholdtoolstest.cpp
    test/algorithms/tthresholdconfig.cpp
    test/msio/tbaselinereader.cpp
    test/structures/image2dtest.cpp
    test/structures/tantennainfo.cpp
    test/structures/tearthposition.cpp
    test/structures/tfieldinfo.cpp
    test/structures/timefrequencydatatest.cpp
    test/structures/tmask2d.cpp
    test/structures/tversionstring.cpp
    test/util/numberparsertest.cpp)
  target_link_libraries(runtests aoflagger-lib ${ALL_LIBRARIES}
                        Boost::unit_test_framework)
  target_include_directories(runtests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  add_dependencies(runtests concatenate_frequency_integration_test_ms)

  add_test(
    runtests
    runtests
    -f
    JUNIT
    -k
    runtests.xml
    --catch_system_error=yes)
  add_custom_target(
    check
    COMMAND runtests
    DEPENDS runtests)

  set(TEST_DATA_DIR "${CMAKE_BINARY_DIR}/test_data")
  configure_file(${CMAKE_SOURCE_DIR}/test/config.h.in
                 ${CMAKE_BINARY_DIR}/test/config.h)

  # Test data for both unit and integration tests.
  ExternalProject_Add(
    concatenate_frequency_integration_test_ms
    EXCLUDE_FROM_ALL URL
    https://support.astron.nl/software/ci_data/aoflagger/L228163_SB150_SB151_uv.dppp.MS.tar.gz
    URL_HASH
      SHA256=1693ae6a56e22b41439246bca30edd1d32119391ac932ade45d0b748ba476e99
    SOURCE_DIR "${TEST_DATA_DIR}/concatenate_frequency"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

  # Set-up cmake configuration for pytest.
  set(INTEGRATION_DIR test/integration)
  set(BINARY_INTEGRATION_DIR ${CMAKE_BINARY_DIR}/${INTEGRATION_DIR})
  file(MAKE_DIRECTORY ${BINARY_INTEGRATION_DIR})
  configure_file(${CMAKE_SOURCE_DIR}/scripts/test/testconfig.py.in
                 ${BINARY_INTEGRATION_DIR}/testconfig.py)
  configure_file(${CMAKE_SOURCE_DIR}/scripts/test/utils.py.in
                 ${BINARY_INTEGRATION_DIR}/utils.py)

  ExternalProject_Add(
    base_integration_test_ms
    URL https://support.astron.nl/software/ci_data/aoflagger/3C196_spw5_sub1-WSRT.MS.tar.gz
    URL_HASH
      SHA256=e516c936058ba19166f5f192d2c68e83df9ca74987422eea9c158b0905ec2f85
    SOURCE_DIR "${TEST_DATA_DIR}/base/3C196_spw5_sub1-WSRT.MS"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

  execute_process(
    COMMAND python3 -m pytest
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    OUTPUT_VARIABLE PYTEST_OUTPUT)
  if(NOT ${PYTEST_OUTPUT} MATCHES "no tests ran")
    message(
      FATAL_ERROR
        "Could not run pytest.\n"
        "Please install pytest or disable testing using -DENABLE_TESTS=Off")
  endif()

  function(add_pytest)
    foreach(TEST ${ARGN})
      set(NAME "integration/${TEST}")
      add_test(
        NAME "${NAME}"
        COMMAND python3 -m pytest -v --junitxml=${CMAKE_BINARY_DIR}/${TEST}.xml
                "${CMAKE_SOURCE_DIR}/${INTEGRATION_DIR}/${TEST}.py"
        WORKING_DIRECTORY "${BINARY_INTEGRATION_DIR}")
      set_tests_properties(
        "${NAME}" PROPERTIES DEPENDS concatenate_frequency_integration_test_ms
                             DEPENDS aoflagger)
    endforeach()
  endfunction()

  add_pytest(tBase tConcatenateFrequency)

endif()
