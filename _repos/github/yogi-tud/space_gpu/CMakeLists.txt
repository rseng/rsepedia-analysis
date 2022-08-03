# convenient cuda clang support was added in 3.19
cmake_minimum_required(VERSION 3.19)


    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")


# don't enable gnu extensions
set(CMAKE_CXX_EXTENSIONS OFF)


project(gpu_compressstore2 LANGUAGES CUDA CXX)

set(INCLUDES
    src
    .
)

set(PREFIX_INCLUDES
    deps/cub
    deps/thrust
)

set(SOURCES
    src/main.cu

    src/cuda_time.cuh
    src/cuda_try.cuh
    src/data_generator.cuh
    src/fast_prng.cuh

    src/kernels/data_generator.cuh
)

#add_custom_target(
 #   cub_fix
 #   sed -i "${CMAKE_CURRENT_LIST_DIR}/deps/cub/cub/iterator/tex_ref_input_iterator.cuh" -e "s|\\(^[^/]*\\)words\\[i\\] = tex1Dfetch(ref, (tex_offset \\* TEXTURE_MULTIPLE) + i);|\\1//words[i] = tex1Dfetch(ref, (tex_offset * TEXTURE_MULTIPLE) + i);|"
 #   VERBATIM
#)
add_executable(gpu_compressstore2 ${SOURCES})


option(AVXPOWER "use avx512-ish" OFF)
if(AVXPOWER)
    target_compile_definitions(gpu_compressstore2 PUBLIC AVXPOWER)
    set(CMAKE_CUDA_FLAGS "-Xcompiler -march=native")
endif()
#To include different compute capabilities, add the architecture X as following:
#PROPERTY CUDA_ARCHITECTURES 75 61 X)
set_property(TARGET gpu_compressstore2 PROPERTY CUDA_ARCHITECTURES 75 61)

target_include_directories(gpu_compressstore2 PRIVATE ${INCLUDES})
target_include_directories(gpu_compressstore2 BEFORE PRIVATE ${PREFIX_INCLUDES})
#add_dependencies(gpu_compressstore2 cub_fix)
