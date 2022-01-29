# Install h5fortran

h5fortran is typically built and installed with CMake or Meson.

## Requirements

* Fortran 2018 compiler (this project uses `submodule` and `error stop`). For example, Gfortran &ge; 7 or Intel oneAPI.
* HDF5 Fortran library (>= 1.8.7, including 1.10.x and 1.12.x)
  * MacOS / Homebrew: `brew install gcc hdf5`
  * Linux / Windows Subsystem for Linux: `apt install gfortran libhdf5-dev`
  * Windows MSYS2: `pacman -S mingw-w64-x86_64-hdf5`
  * build from source (optional): `python scripts/build_hdf5.py`

Note that some precompiled HDF5 libraries have only C / C++ without Fortran.

The library `libh5fortran` is built, link it into your program as usual along with the HDF5 libraries and include files.

## CMake

Build and self-test via:

```sh
cmake -B build
cmake --build build

# optional self-test
cd build
ctest
```

(optional) to install to a directory like ~/h5fortran:

```sh
cmake -B build -DCMAKE_INSTALL_PREFIX=~/h5fortran
cmake --install build
```

### use h5fortran from your project

```sh
cmake -B build -Dh5fortran_ROOT=~/h5fortran/
```

and in your CMakeLists.txt

```cmake
cmake_minimum_required(VERSION 3.14...3.20)
project(myProject LANGUAGES Fortran)

find_package(h5fortran)

if(NOT h5fortran_FOUND)
  include(FetchContent)

  FetchContent_Declare(H5FORTRAN
    GIT_REPOSITORY https://github.com/geospace-code/h5fortran.git
    GIT_TAG v3.6.6)
  FetchContent_MakeAvailable(H5FORTRAN)
endif()

# --- your project targets:

add_executable(myProj main.f90)
target_link_libraries(myProj PRIVATE h5fortran::h5fortran)
```

and where the main.f90 is like:

```fortran
program main

use h5fortran, only : hdf5_file
implicit none

type(hdf5_file) :: h5f

call h5f%open('h5fortran_example2.h5', action='w')
call h5f%write('/x', 123)
call h5f%close()


end program
```

### [optional] create distributable archive

If you wish to create a package archive that is usable on systems with compatible Fortran ABI, after building:

```sh
cpack --config build/CPackConfig.cmake
```

### [optional] specify a particular HDF5 library

```sh
cmake -DHDF5_ROOT=/path/to/hdf5lib -B build
```

or set environment variable `HDF5_ROOT=/path/to/hdf5lib`

## Meson

To build h5fortran as a standalone project

```sh
meson build

meson test -C build
```

### h5fortran Meson subproject

To include h5fortran as a Meson subproject, in the main project meson.build (that uses h5fortran) have like:

```meson
hdf5_proj = subproject('h5fortran')
hdf5_interface = hdf5_proj.get_variable('hdf5_interface')

my_exe = executable('myexe', 'main.f90', dependencies: hdf5_interface)
```

and have a file in the main project `subprojects/h5fortran.wrap` containing:

```ini
[wrap-git]
directory = h5fortran
url = https://github.com/geospace-code/h5fortran.git
revision = v3.6.5
```

## Standalone compiler wrapper

h5fortran can be used from the
[HDF5 compiler wrapper "h5fc"](https://support.hdfgroup.org/HDF5/Tutor/compile.html) like:

```sh
h5fc -I~/h5fortran/include myprogram.f90 ~/h5fortran/lib/libh5fortran.a
```
# h5fortran Examples

All examples assume:

```fortran
use h5fortran, only: hdf5_file
type(hdf5_file) :: h5f
```

* gzip compression may be applied for rank &ge; 2 arrays by setting `comp_lvl` to a value between 1 and 9.
  Shuffle filter is automatically applied for better compression
* string attributes may be applied to any variable at time of writing or later.
* h5f%open(..., `comp_lvl=1`) option enables GZIP compression., where comp_lvl is from 1 to 9. bigger comp_lvl gives more compression but isslower to write.

## Create new HDF5 file, with variable "value1"

```fortran
call h5f%open('test.h5', action='w')

call h5f%write('/value1', 123.)

call h5f%close()
```

## create soft links to actual variable

HDF5 soft link variables: arbitrarily many soft-linked variable names can point to an actual variable, which need not yet exist.

```fortran
call h5f%write("/x", 42)
call h5f%softlink("/x", "/y")
call h5f%softlink("/x", "/z")
```

/z and /y are soft links to /x, which need not yet exist.

## ensure all files are flushed to disk at end of program

If your program opens lots of HDF5 files and you're worried about being sure they're all flushed to disk, make this call near the very end of the program.
This flushes and closes ALL HDF5 files, even those that may be invoked directly from the HDF5 library without h5fortran.

```fortran
call hdf5_close()
```

Normally, you should be calling `%close()` on each file to flush to disk when done using a file.
If `%close()` or hdf5_close is not called, data loss can result.

```fortran
call h5f%close()
```

At any time during the program, the `%flush()` method can be called to request the operating system to write a file to disk.
This could be useful during a long-running program (say, an HPC simulation) to help ensure data isn't lost of an HDF5 file is open for a long time.
The flush request is on a per-file basis, so if multiple files are open, flush each file to protect against data loss in this case.

```fortran
call h5f%flush()
```

## read / write attributes to variable

Note that HDF5 character attributes are scalar while int32, real32, real64 attributes are 1D vectors.

Assume variable "/x" exists and then see these examples:

### write attributes

```fortran
call h%writeattr('/x', 'note','this is just a little number')
call h%writeattr('/x', 'hello', 'hi')
call h%writeattr('/x', 'life', [42])
call h%writeattr('/x', 'life_float', [42._real32, 84._real32])
call h%writeattr('/x', 'life_double', [42._real64])
```

or for the high level interface:

```fortran
call h5write_attr('myfile.h5', '/x', 'date', [2020,4,1])

call h5write_attr('myfile.h5', '/x', 'units', 'Nm^-2')
```

### read attributes

For attributes, HDF5 character values are *space-terminated* instead of null terminated.

```fortran
character(1024) :: attr_str
integer :: attr_int(1)
real(real32) :: attr32(2)
real(real64) :: attr64(1)

call h%readattr('/x', 'note', attr_str)
if (attr_str /= 'this is just a little number') error stop 'readattr value note'

call h%readattr('/x', 'life', attr_int)
call h%readattr('/x', 'life_float', attr32)
call h%readattr('/x', 'life_double', attr64)
```

or for the high level interface:

```fortran
integer :: idate(3)
character(16) :: unit_str

call h5read_attr('myfile.h5', '/x', 'date', idate)

call h5read_attr('myfile.h5', '/x', 'units', unit_str)
```

## Add/append variable "value1" to existing HDF5 file "test.h5"

* if file `test.h5` exists, add a variable to it
* if file `test.h5` does not exist, create it and add a variable to it.

```fortran
call h5f%open('test.h5', action='rw')

call h5f%write('/value1', 123.)

call h5f%close()
```

## Add gzip compressed 3-D array "value2" to existing HDF5 file "test.h5"

```fortran
real :: val2(1000,1000,3) = 0.

call h5f%open('test.h5', comp_lvl=1)

call h5f%write('/value2', val2)

call h5f%close()
```

chunk_size may optionally be set in the `%write()` method for 2-d to 7-d arrays.
compression and chunking are disabled if any element of chunk_size is less than 1
chunk_size may be manually specified in write() otherwise it will be set automatically.

Currently, data is written contiguous or compact if not compressed and is only chunked if compression is used.

## check if a variable exists

the logical method %exist() checks if a dataset (variable) exists in the opened HDF5 file.

```fortran
exists = h5f%exist("/A")
```

A convenience method that checks existence of a dataset without creating the h5 object manually is:

```fortran
exists = h5exist("my.h5", "/A")
```

## check variable shape, rank/ndims

`h5f%ndim` we didn't use `%rank` to avoid confusion with intrinsic "rank()"

```fortran
call h5f%open('test.h5', action='r')

integer :: drank
integer(hsize_t), allocatable :: dims(:)

drank = h5f%ndim('/A')
call h5f%shape('/A',dims)

if (drank /= size(dims)) error stop
```

## Read scalar, 3-D array of unknown size

```fortran
call h5f%open('test.h5', action='r')

integer(hsize_t), allocatable :: dims(:)
real, allocatable :: A(:,:,:)

call h5f%shape('/A',dims)
allocate(A(dims(1), dims(2), dims(3)))
call h5f%read('/A', A)

call h5f%close()
```

## read slice (part of) a disk array

Reading a disk HDF5 array into a variable of matching shape is done with `istart=` and `iend=` arguments, which have 1-D arguments for the start and stop index desired from each dimension.

For example, support HDF5 disk variable "/A" is shape (10,20,30) and you wish to read just part of this array like:

* dim 1: 5-7
* dim 2: 1-5
* dim 3: 2-8

then do:

```fortran
real, dimension(3,5,7) :: A

call h5f%open('test.h5', action='r')

call h5f%read('/A', A, istart=[5, 1, 2], iend=[7, 5, 8])
```

## write slice (part of) a disk array

Writing a disk HDF5 array from a variable of matching shape is done with `istart=` and `iend=` arguments, which have 1-D arguments for the start and stop index desired from each dimension.

For example, support HDF5 disk variable "/A" is shape (10,20,30) and you wish to write a slice from a variable shaped (5,7,1) with start/stop indices:

* dim 1: 3-7
* dim 2: 4-10
* dim 3: 8

then do:

```fortran
real, dimension(5,7,1) :: A

call h5f%open('test.h5')

call h5f%create('/A', H5T_NATIVE_REAL, [5,7,1])
call h5f%write('/A', A, istart=[3, 4, 8], iend=[7, 10, 8])
```

Note the h5f%create() call to open the disk variable.
This step is also needed with h5py in Python or Matlab HDF5 h5create() before h5write().

## is dataset compact, contiguous, or chunked

Assume file handle h5f was already opened, the logical status is inspected:

```fortran
is_compact = h5f%is_compact("/A")

is_contig = h5f%is_contig('/A')

is_chunked = h5f%is_chunked('/A')
```

## get chunk size

if dataset is not chunked, chunk_size == -1

```sh
call h5f%chunks('/A', chunk_size)
```

## Create group "scope"

```fortran
real :: val2(1000,1000,3) = 0.

call h5f%open('test.h5')

call h5f%write_group('/scope/')

call h5f%close()
```

## verbose / debug

set options debug and /or verbose for diagnostics

```sh
call h5f%open(..., verbose=.true., debug=.true.)
```

## Permissive syntax

We make the hdf5%open(..., action=...) like Fortran open()

* overwrite (truncate) existing file: open with `action='w'`
* append to existing file or create file: `action='rw'`
src_dir: ./src
output_dir: ./docs
project: Object-oriented Fortran HDF5 interface
project_github: https://github.com/geospace-code/h5fortran
project_website: https://geospace-code.github.io/h5fortran
summary: Object-oriented Fortran HDF5 interface
author: Michael Hirsch, Ph.D.
github: https://github.com/geospace-code
license: by
exclude: CMakeFortranCompilerId.F
         reader_lt_template.in.f90
         writer_lt_template.in.f90
         writer_template.in.f90
         writer_template_i32.in.f90
         writer_template_r64.in.f90
         writer_template_r32.in.f90
         writer_template_i64.in.f90
display: public
         protected
         private
source: false
graph: true
search: true


Simple, robust, thin HDF5 polymorphic Fortran read/write interface.
Reading or writing {real64,real32,int32,int64} from scalar to 7d is as simple as

```fortran
use h5fortran

call h5write('my.h5', '/x', x)

call h5read('my.h5', '/y', y)
```

For NetCDF4 see [nc4fortran](https://github.com/geospace-code/nc4fortran/).
h5fortran is designed for "serial" HDF5 read/write.
We don't yet implement the interface for "parallel" HDF5.

h5fortran is designed for easy use using **static** or **shared** linking from your project via:

* `cmake --install`
* CMake ExternalProject
* CMake FetchContent
* CMake + Git submodule
* Meson subproject

Uses Fortran `submodule` for clean template structure.
This easy-to-use, thin object-oriented modern Fortran library abstracts away the messy parts of HDF5 so that you can read / write various types/ranks of data with a single command.
In distinction from other high-level HDF5 interfaces, h5fortran works to deduplicate code, using polymorphism wherever feasible and extensive test suite.

Polymorphic [API](./API.md) with read/write for types int32, real32, real64 with rank scalar (0-D) through 7-D.
64-bit integers int64 are read/write from scalar through 3-D.

as well as **character (string)**.
If you need int64, we have a working example for that: src/concepts/int64.f90 that can easily be put into the h5fortran API--just make a GitHub Issue.

* HDF5 **attributes** are also supported for read/write with type character, int32, real32, real64.
* **Array slicing on read and write** is supported, that is, reading or writing part of a disk HDF5 array into a variable matching the slice shape.
* Mismatched datatypes are coerced as per standard Fortran rules. For example, reading a float HDF5 variable into an integer Fortran variable:  42.3 => 42
* Zlib (deflate) compression / decompression -- h5fortran will work without Zlib, but will save/load uncompressed data only.
* create HDF5 soft link variables--arbitrarily many soft-linked variable names can point to an actual variable, which need not yet exist.

Tested on systems with HDF5 1.8, 1.10 and 1.12 including:

* MacOS (homebrew)
* Linux (Ubuntu, CentOS)
* Windows Subsystem for Linux
* Windows MSYS2
* IBM Power with Gfortran

Compilers known to work include:

* GCC (gfortran) &ge; 7
* Intel oneAPI HPC compiler &ge; 2021

---

We welcome [contributions](https://github.com/geospace-code/.github/blob/main/CONTRIBUTING.md).
In general we hold to the geospace-code [code of conduct](https://github.com/geospace-code/.github/blob/main/CODE_OF_CONDUCT.md).

## Build

Using CMake:

```sh
git clone https://github.com/geospace-code/h5fortran.git

cd h5fortran

cmake -B build

cmake --build build
```

for more details see [Install.md](./Install.md).

For general use with non-CMake build systems, "h5fortran.pc" pkg-config file is also generated / installed.

### Autobuild HDF5

h5fortran will automatically build the HDF5 and ZLIB libraries if needed.
This is useful as many HPC have broken or ABI-incompatible HDF5 libraries installed.
Building HDF5 and ZLIB takes about a minute on a typical laptop.
To disable this autobuild behavior, use option:

```sh
cmake -B build -Dautobuild=off
```

To force building the HDF5 and ZLIB libraries, to gain better performance via optimizing for your system's CPU:

```sh
cmake -Dhdf5_external=on

cmake --build build
```

NOTE: If using Intel oneAPI on Windows, ensure that environment variable CC=icl as set manually in the command prompt:

```posh
set CC=icl
set FC=ifort
```

This is necessary to workaround techniques used by HDF5 CMake files that don't pickup the CMake `set(ENV{CC})`.
Otherwise, HDF5 build failures may result due to defaulting to icl-clang.

By default we use Zlib 2.x a.k.a. zlib-ng.
If you have a problem with Zlib-ng on your system, try the unmaintained Zlib 1.x by:

```sh
cmake -B build -Dzlib_legacy=on
```

## Usage

The simplest [example](./Examples/) h5fortran usage is like:

```fortran
use h5fortran

call h5write('golt.h5','/x', [1,2,3,4,5,6])
```

or

```fortran
use h5fortran

real :: x2

if(.not. is_hdf5('golt.h5')) error stop 'golt.h5 is not an HDF5 file'

call h5read('golt.h5', '/x', x2)
```

For detailed [examples](./Examples/) see [Examples.md](./Examples.md).

## Notes

* The first character of the filename should be a character, NOT whitespace to avoid file open/creation errors.
* Polymorphic array rank is implemented.

### Getting HDF5 library

Instead of auto-building HDF5 via H5Fortran, one may build and install the HDF5 library by:

```sh
python3 scripts/build_hdf5.py
```

### h5fortran: missing Fortran datatypes

* arrays of rank > 7: prototyped in reader_nd.f90, writer_nd.f90. Fortran 2008 arrays are up to rank 15, but only recent compilers support.
* complex64 / complex128: not natively handled in HDF5. There are performance impacts for compound datatypes. Many choose to write two datasets, one each for real and imaginary like `A_real` and `A_imag`
* non-default character kind
* logical / boolean: not supported natively by HDF5. h5py implements as [enum](https://docs.h5py.org/en/stable/faq.html#what-datatypes-are-supported).
# h5fortran API

This document provides a listing of h5fortran `public` scoped user-facing procedures and methods with a summary of their parameters.

All examples assume:

```fortran
use h5fortran, only: hdf5_file
use hdf5, only: HSIZE_T, HID_T

type(hdf5_file) :: h
```

Query HDF5 library version:

```fortran
use h5fortran, only : hdf5version
print *, hdf5version()
```

## Open / close HDF5 file reference

More than one HDF5 file can be open in a program, by declaring unique file handle (variable) like:

```fortran
type(hdf5_file) :: h1, h2, h3
```

```fortran
call h%open(filename,action,comp_lvl,verbose,debug)
!! Opens hdf5 file

character(*), intent(in) :: filename
character(*), intent(in), optional :: action  !< r, w, rw
integer, intent(in), optional      :: comp_lvl  !< 0: no compression. 1-9: ZLIB compression, higher is more compressior
logical, intent(in), optional      :: verbose, debug
```

```fortran
call h%close(close_hdf5_interface)
!! This must be called on each HDF5 file to flush buffers to disk
!! data loss can occur if program terminates before this procedure
!!
!! We don't reference count because applications might also invoke HDF5
!! directly.
!! close_hdf5_interface is when you know you have exactly one HDF5 file in your
!! application, if true it closes ALL files, even those invoked directly from HDF5.

logical, intent(in), optional :: close_hdf5_interface
```

To avoid memory leaks or corrupted files, always "finalize" all hDF5 files before STOPping the Fortran program.

```fortran
call h%flush()
!! request operating system flush data to disk.
```

## Disk variable (dataset) inquiry

To allocate variables before reading data, inquire about dataset characteristics with these procedures.

```fortran
rank = h%ndim(dataset_name)

character(*), intent(in) :: dataset_name
```

Get disk dataset shape (1D vector)

```fortran
call h%shape(dataset_name, dims)
character(*), intent(in) :: dataset_name
integer(HSIZE_T), intent(out), allocatable :: dims(:)
```

Dataset "dname" data class (i.e. integer, float, string, ...)

```fortran
integer :: class
!! H5T_INTEGER_F, H5T_FLOAT_F, H5T_STRING_F
class = h%class(dname)
character(*), intent(in) :: dname
```

Dataset "dname" datatype

```fortran
integer(HID_T) :: dtype
!! H5T_NATIVE_REAL, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_NATIVE_CHARACTER, H5T_STD_I64LE
dtype = h%dtype(dname)
character(*), intent(in) :: dname
```

Does dataset "dname" exist in this HDF5 file?

```fortran
exists = h%exist(dname)
character(*), intent(in) :: dname
```

Is dataset "dname" contiguous on disk?

```fortran
tf = h%is_contig(dname)
!! is dataset contiguous
character(*), intent(in) :: dname
```

Is dataset compact (< 64K)

```fortran
tf = h%is_compact(dname)
!! is dataset compact layout
character(*), intent(in) :: dname
```

Is dataset chunked?

```fortran
tf = h%is_chunked(dname)
!! is dataset chunked
character(*), intent(in) :: dname
```

Is this an HDF5 file?

```fortran
use h5fortran, only: is_hdf5

tf = is_hdf5('myfile.txt')  !< probably false
tf = is_hdf5('myfile.h5')  !< true if a valid HDF5 file
```

These are more advanced inquiries into the memory layout of the dataset, for advanced users:

```fortran
Layout = h%layout(dname)
!! integer :: H5D_CONTIGUOUS_F, H5D_CHUNKED_F, H5D_VIRTUAL_F, H5D_COMPACT_F
character(*), intent(in) :: dname
```

```fortran
call h%chunks(dname, chunk_size)
character(*), intent(in) :: dname
integer, intent(out) :: chunk_size(:)
```

## create dataset softlink

One of the key features of HDF5 is the ability to create dataset softlinks within an HDF5 file:

```fortran
call h%softlink(target, link)
character(*), intent(in) :: target, &  !< target path to link dataset
                            link  !< soft link path to create
```

## file write operations

```fortran
call h%write(dname,value, chunk_size, istart, iend, stride, compact)
!! write 0d..7d dataset
character(*), intent(in) :: dname
class(*), intent(in) :: value(:)  !< array to write
integer, intent(in), optional :: chunk_size(rank(value))
integer, intent(in), optional, dimension(:) :: istart, iend, stride  !< array slicing
logical, intent(in), optional :: compact  !< faster I/O for sub-64 kB datasets
```

Write dataset attribute (e.g. units or instrument)

```fortran
call h%writeattr(dname, attr, attrval)
character(*), intent(in) :: dname, attr  !< dataset name, attribute name
class(*), intent(in) :: attrval(:)  !< character, real, integer
```

## file read operations

Read data from disk to memory

```fortran
call h%read(dname, value, istart, iend, stride)
character(*), intent(in)         :: dname
class(*), intent(out) :: value(:)  !< read array to this ALLOCATED variable
integer, intent(in), optional, dimension(:) :: istart, iend, stride !< array slicing
```

Read dataset attribute into memory

```fortran
call h%readattr(dname, attr, attrval)
character(*), intent(in) :: dname, attr  !< dataset name, attribute name
class(*), intent(out) :: attrval(:)  !< character, real, integer
```

## high level operations

These are single-call operations that are slower than the object-oriented methods above.
The runtime penalty may be insignificant unless you call these functions many times, say in a for loop.

The `h5write` opens `filename` with `action='rw'` (create if not present, append if existing).

```fortran
call h5write(filename, dname, value)
character(*), intent(in) :: filename, dname
class(*), intent(in) :: value(:)
```

The `h5read` opens `filename` with `action='r'` (error if file not exist).

```fortran
call h5read(filename, dname, value)
character(*), intent(in) :: filename, dname
class(*), intent(out) :: value(:)
```
# Object-oriented Fortran HDF5 interface

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02842/status.svg)](https://doi.org/10.21105/joss.02842)
[![DOI](https://zenodo.org/badge/128736984.svg)](https://zenodo.org/badge/latestdoi/128736984)

![ci_linux](https://github.com/geospace-code/h5fortran/workflows/ci/badge.svg)
![ci_macos](https://github.com/geospace-code/h5fortran/workflows/ci_macos/badge.svg)
![ci_windows](https://github.com/geospace-code/h5fortran/workflows/ci_windows/badge.svg)
![ci_meson](https://github.com/geospace-code/h5fortran/workflows/ci_meson/badge.svg)
[![intel-oneapi](https://github.com/geospace-code/h5fortran/actions/workflows/intel-oneapi.yml/badge.svg)](https://github.com/geospace-code/h5fortran/actions/workflows/intel-oneapi.yml)

Simple, robust, thin HDF5 polymorphic Fortran read/write interface.
Reading or writing {real64,real32,int32,int64} from scalar to 7d is as simple as

```fortran
use h5fortran

call h5write('my.h5', '/x', x)

call h5read('my.h5', '/y', y)
```

For NetCDF4 see [nc4fortran](https://github.com/geospace-code/nc4fortran/).
h5fortran is designed for "serial" HDF5 read/write.
We don't yet implement the interface for "parallel" HDF5.

h5fortran is designed for easy use using static or shared linking from your project via:

* `cmake --install`
* CMake ExternalProject
* CMake FetchContent
* CMake + Git submodule
* Meson subproject

Uses Fortran `submodule` for clean template structure.
This easy-to-use, thin object-oriented modern Fortran library abstracts away the messy parts of HDF5 so that you can read / write various types/ranks of data with a single command.
In distinction from other high-level HDF5 interfaces, h5fortran works to deduplicate code, using polymorphism wherever feasible and extensive test suite.

Polymorphic [API](./API.md) with read/write for types int32, real32, real64 with rank scalar (0-D) through 7-D.
64-bit integers int64 are read/write from scalar through 3-D.

as well as **character (string)**.
If you need int64, we have a working example for that: src/concepts/int64.f90 that can easily be put into the h5fortran API--just make a GitHub Issue.

* HDF5 **attributes** are also supported for read/write with type character, int32, real32, real64.
* **Array slicing on read and write** is supported, that is, reading or writing part of a disk HDF5 array into a variable matching the slice shape.
* Mismatched datatypes are coerced as per standard Fortran rules. For example, reading a float HDF5 variable into an integer Fortran variable:  42.3 => 42
* Zlib (deflate) compression / decompression -- h5fortran will work without Zlib, but will save/load uncompressed data only.
* create HDF5 soft link variables--arbitrarily many soft-linked variable names can point to an actual variable, which need not yet exist.

Tested on systems with HDF5 1.8, 1.10 and 1.12 including:

* MacOS (homebrew)
* Linux (Ubuntu, CentOS)
* Windows Subsystem for Linux
* Windows MSYS2
* IBM Power with Gfortran

Compilers known to work (tested on CI) include:

* GCC (gfortran) &ge; 8
* Intel oneAPI &ge; 2021

---

We welcome [contributions](https://github.com/geospace-code/.github/blob/main/CONTRIBUTING.md).
In general we hold to the geospace-code [code of conduct](https://github.com/geospace-code/.github/blob/main/CODE_OF_CONDUCT.md).

## Build

Using CMake:

```sh
git clone https://github.com/geospace-code/h5fortran.git

cd h5fortran

cmake -B build

cmake --build build
```

for more details see [Install.md](./Install.md).

For general use with non-CMake build systems, "h5fortran.pc" pkg-config file is also generated / installed.

To save time, if not intended to use self-tests, you can skip the build of the test suite:

```sh
cmake -B build -DBUILD_TESTING=off
```

### Autobuild HDF5

h5fortran will automatically build the HDF5 and ZLIB libraries if needed.
This is useful as many HPC have broken or ABI-incompatible HDF5 libraries installed.
Building HDF5 and ZLIB takes about a minute on a typical laptop.
To disable this autobuild behavior, use option:

```sh
cmake -B build -Dautobuild=off
```

To force building the HDF5 and ZLIB libraries, to gain better performance via optimizing for your system's CPU:

```sh
cmake -Dhdf5_external=on

cmake --build build
```

NOTE: If using Intel oneAPI on Windows, ensure that environment variable CC=icl as set manually in the command prompt:

```posh
set CC=icl
set FC=ifort
```

This is necessary to workaround techniques used by HDF5 CMake files that don't pickup the CMake `set(ENV{CC})`.
Otherwise, HDF5 build failures may result due to defaulting to icl-clang.

By default we use Zlib 2.x a.k.a. zlib-ng.
If you have a problem with Zlib-ng on your system, try the unmaintained Zlib 1.x by:

```sh
cmake -B build -Dzlib_legacy=on
```

The user can request a specific HDF5 version like:

```sh
cmake -B build -DHDF5_VERSION=1.12.1
```

Only the HDF5 versions in cmake/libraries.json would work.

## Usage

The simplest [example](./Examples/) h5fortran usage is like:

```fortran
use h5fortran

call h5write('golt.h5','/x', [1,2,3,4,5,6])
```

or

```fortran
use h5fortran

real :: x2

if(.not. is_hdf5('golt.h5')) error stop 'golt.h5 is not an HDF5 file'

call h5read('golt.h5', '/x', x2)
```

For detailed [examples](./Examples/) see [Examples.md](./Examples.md).

## Notes

* The first character of the filename should be a character, NOT whitespace to avoid file open/creation errors.
* Polymorphic array rank is implemented.

### Getting HDF5 library

Instead of auto-building HDF5 via H5Fortran, one may build and install the HDF5 library by:

```sh
python3 scripts/build_hdf5.py
```

### h5fortran: missing Fortran datatypes

* arrays of rank > 7: prototyped in reader_nd.f90, writer_nd.f90. Fortran 2008 arrays are up to rank 15, but only recent compilers support.
* complex64 / complex128: not natively handled in HDF5. There are performance impacts for compound datatypes. Many choose to write two datasets, one each for real and imaginary like `A_real` and `A_imag`
* non-default character kind
* logical / boolean: not supported natively by HDF5. h5py implements as [enum](https://docs.h5py.org/en/stable/faq.html#what-datatypes-are-supported).
---
title: 'h5fortran: object-oriented polymorphic Fortran interface for HDF5 file IO'
tags:
authors:
  - name: Michael Hirsch
    orcid: 0000-0002-1637-6526
    affiliation: "1"
affiliations:
 - name: Boston University
   index: 1
date: 6 Nov 2020
bibliography: paper.bib
---

# Summary

h5fortran [@h5fortran] is a Fortran interface to HDF5 that abstracts away most details of a frequently-used subset of HDF5 file read/write operations.
h5fortran has object-oriented and functional interfaces that makes HDF5 use from Fortran as easy as in high-level languages.
For example, using the official HDF5 Fortran interface to write an array to disk involves over twenty subroutine calls, plus the associated program logic to call those procedures with efficient parameters.
The same array write operation with h5fortran is accomplished with a single subroutine call.
A similar reduction in user-facing complexity is achieved by h5fortran for HDF5 array reads.
h5fortran adheres to Fortran 2008 standard, working on Gfortran and Intel compilers for Linux, MacOS, Windows on Intel / AMD, ARM and IBM POWER systems.
CMake is used to build h5fortran, detecting if the HDF5 library is present and working and building HDF5 from source if necessary.
CPack can be used to generate distributable binary archives.
Conan package manager may also be used to automatically install the HDF5 library and then install h5fortran.
Meson build system is also supported by h5fortran.

h5fortran has general applicability to projects needing to do any of:

* writing variables to HDF5: scalar to 7-D, of type real32, real64 or integer
* reading variables from HDF5: scalar to 7-D, of type real32, real64 or integer
* read / write character variables to / from HDF5 file
* read / write variable attributes to / from HDF5 file
* get the shape of a disk variable to allocate a memory variable for reading that data

HDF5 does not have native support for complex numbers or booleans (`logical` Fortran datatype).
h5fortran does not yet support complex numbers or booleans.
High-level HDF5 interfaces in other code languages such as h5py have implemented these types using HDF5 struct and HDF5 enum respectively.
If the HDF5 community continues to coalesce around these *de facto* data type implementations, we may consider implementing them in h5fortran in the future.
h5fortran currently supports the serial HDF5 interface, and does not yet support the parallel MPI-based HDF5 interface.
h5fortran does not yet support extensible datasets, although we would be open to investigating adding this support upon community request.

In addition to the object-oriented interface, h5fortran provides single-command read / write procedures.
Array slicing allows reading or writing a portion of a large disk variable to/from RAM.
If the user has HDF5 with SZIP or ZLIB compression enabled, h5fortran is capable of reading and writing compressed variables, which can save over 50% disk space depending on the data lossless compressibility.
Data shuffling and Fletcher32 checksums provide better compression and a check of file integrity respectively.
h5fortran was designed for use by individual users on their laptops or embedded devices, as well as for use in HPC applications where parallel tasks need read only part of a milestone or shared HDF5 variable.

h5fortran was originally developed for the GEMINI [@gemini3d; @zettergren] ionospheric model, funded in part by NASA ROSES \#80NSSC20K0176 and DARPA Cooperative Agreement HR00112120003.
This work is approved for public release; distribution is unlimited.
The information does not necessarily reflect the position or the policy of the Government.

# Statement of need

Fortran has only raw file input-output (IO) built in to the language.
To support reproducibility of work done in any programming language and long-term usefulness of the data generated or processed, it is beneficial to use self-describing data file formats like HDF5 [@2011hdf5].
Many popular languages and libraries used for simulation and data science use the HDF5 [@hdf5] file IO library.
Most programs and libraries intended for use by practitioners such as modelers and data scientists themselves use an object-oriented HDF5 interface like h5py [@h5py; @h5pybook].

# Other programs

While other HDF5 interfaces exist, h5fortran presents a broad set of commonly used features, comprehensive test coverage and robustness across compilers and computing systems.
We have written a companion library for NetCDF4 called nc4fortran [@nc4fortran], which by design has a nearly identical user-facing API.
Other Fortran HDF5 interfaces such as HDF5_Utils [@hdf5_utils] use a functional interface mimicking the HDF5 LT functions, which require the end user to keep track of extra variables versus the single object used by h5fortran.
A package for C++ with similar priorities of using modern language features and simple commands is h5pp [@h5pp].
A template-based C++ implementation is provided in h5xx [@h5xx].

# References
# h5fortran Examples

From the h5fortran/ directory, specify the h5fortran install like:

```sh
cmake -B build -DCMAKE_INSTALL_PREFIX=~/h5fortran
cmake --build build
cmake --install build

cmake -B Examples/build -S Examples -Dh5fortran_ROOT=~/h5fortran
```

## Example 1

Example 1 shows the functional one-step interface of h5fortran

## Example 2

Example 2 shows the object-oriented interface of h5fortran, which may offer faster performance if more than one variable is being read or written.

## Example 3

Example 3 is of a C main program calling a Fortran interface to h5fortran

## Example 4

Example 4 is of a C++ main program calling a Fortran interface to h5fortran

For a C++ header-only object-oriented HDF5 library, consider [HighFive](https://github.com/BlueBrain/HighFive)


## Notes

### Non CMake build

CMake makes building much easier.
If for whatever you don't wish to use CMake, the HDF5 compiler wrapper (if available on your system) may work.
Since the HDF5 compiler wrapper is not always working or available, we strongly recommended CMake as above for any HDF5-based application.

On Ubuntu it looks like:

```sh
$ h5fc -show

gfortran -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5hl_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial
```
# Build HDF5 scripts

Here we use a dummy CMake project to reuse code from the main h5fortran project to build HDF5 and ZLIB.

```sh
cmake -B build -DCMAKE_INSTALL_PREFIX=~/local
cmake --build build
```

Optionally, build the MPI layer (parallel HDF5)

```sh
cmake -B build -DCMAKE_INSTALL_PREFIX=~/local -Dhdf5_parallel=on
cmake --build build
```

## compiler scripts

The "h5fc" and "h5cc" will be installed above under ~/local/bin.
These only appear on Linux and MacOS.
The HDF5 build system doesn't yet build these compiler wrappers on Windows.
