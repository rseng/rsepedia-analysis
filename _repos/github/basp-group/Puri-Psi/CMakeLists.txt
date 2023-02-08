cmake_minimum_required(VERSION 2.8)
project(PURIPSI C CXX)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_files)

option(tests          "Enable testing"                                  on)
option(examples       "Compile examples"                                on)
#option(data          "Download measurement set for testing"            on)
option(benchmarks     "Enable benchmarking"                             on)
option(openmp         "Enable OpenMP"                                   on)
option(mpi	      "Enable MPI"                                      on)
option(logging        "Enable logging"                                  on)
option(mkl            "Use MKL for Eigen acceleration"                  on)
option(blas	      "Use BLAS for Eigen acceleration"                 on)
option(hdf5           "Use HDF5 to write/read data"                     on)
# Set version and build id of this package
include(version)

# we are doing c++11
include(AddCPP11Flags)

# sets up rpath so libraries can be found
include(rpath)

# adds logging variables
include(logging)

# include exernal dependencies
include(dependencies)

#if(tests AND data)
#    include(data_dependencies)
#endif()

if(tests)  # Adds ctest
    enable_testing()
    include(AddCatchTest)
endif()

if(benchmarks)
   include(AddBenchmark)
endif()

add_subdirectory(cpp)

# Exports PuriPsi so other packages can access it
include(export_puripsi)
