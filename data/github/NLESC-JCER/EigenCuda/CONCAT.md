# Change log

# [0.4.0] 10/02/2020
### Changed
  - Split the memory management (`CudaMatrix`) from the [CUBLAS](https://docs.nvidia.com/cuda/cublas/index.html) invocation (`CudaPipeline`)
  - Moved all the allocation to the smart pointers inside `CudaMatrix`
 - Removed unused headers

# [0.3.0] 26/09/2019
### Added
 - Smart pointers to handle cuda resources
 - New CudaMatrix class
 - Use Eigen::MatrixXd
 - Check available memory in the GPU before computing

### Removed
 - Template class, implementation only for double available
 - Triple tensor product
 - Shapes struct


# [0.2.0] 27/08/2019
### Added
 - Tensor matrix multiplacation using [gemmbatched](https://docs.nvidia.com/cuda/CUBLAS/index.html#CUBLAS-lt-t-gt-gemmbatched).
 - [Async calls](https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY.html#group__CUDART__MEMORY_1g85073372f776b4c4d5f89f7124b7bf79) to memory copies.
 - Properly free memory after the tensor operation is done.

# [0.1.0]

### New
 - Use a template function to perform matrix matrix multiplacation using [CUBLAS](https://docs.nvidia.com/cuda/CUBLAS/index.html).
 - Use either *pinned* (**default**) or *pageable* memory, see [cuda optimizations](https://devblogs.nvidia.com/how-optimize-data-transfers-cuda-cc/).
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3660936.svg)](https://doi.org/10.5281/zenodo.3660936)

# EigenCuda

Offload the [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) matrix matrix multiplication to an Nvidia GPU
using [cublas](https://docs.nvidia.com/cuda/cublas/index.html).

## CMake Installation

To compile execute:
```
cmake -H. -Bbuild && cmake --build build
```

To Debug compile as:
```
cmake -H. -Bbuild  -DCMAKE_BUILD_TYPE=Debug && cmake --build build
```

## Dependencies

This packages assumes that you have installed the following packages:
  
  * [Cudatoolkit](https://anaconda.org/anaconda/cudatoolkit)
  * [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page)

## Usage
### Matrix Multiplication
```cpp
#include "eigencuda.hpp"
#include "cudapipeline.hpp"

using eigencuda::CudaPipeline;
using eigencuda::CudaMatrix;

  // Call the class to handle GPU resources
  CudaPipeline cuda_pip;

Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2, 2);
Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2);
Eigen::MatrixXd C = Eigen::MatrixXd::Zero(3, 2);
Eigen::MatrixXd D = Eigen::MatrixXd::Zero(3, 2);
Eigen::MatrixXd X = Eigen::MatrixXd::Zero(3, 2);
Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(3, 2);
Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(3, 2);

// Define matrices
A << 1., 2., 3., 4.;
B << 5., 6., 7., 8., 9., 10.;
C << 9., 10., 11., 12., 13., 14.;
D << 13., 14., 15., 16., 17., 18.;
X << 23., 34., 31., 46., 39., 58.;
Y << 39., 58., 47., 70., 55., 82.;
Z << 55., 82., 63., 94., 71., 106.;

std::vector<Eigen::MatrixXd> tensor{B, C, D};
std::vector<Eigen::MatrixXd> results(3, Eigen::MatrixXd::Zero(3, 2));
CudaMatrix cuma_A{A, cuda_pip.get_stream()};
CudaMatrix cuma_B{3, 2, cuda_pip.get_stream()};
CudaMatrix cuma_C{3, 2, cuda_pip.get_stream()};

for (Index i = 0; i < 3; i++) {
  cuma_B.copy_to_gpu(tensor[i]);
  cuda_pip.gemm(cuma_B, cuma_A, cuma_C);
  results[i] = cuma_C;
}
// Expected results
bool pred_1 = X.isApprox(results[0]);
bool pred_2 = Y.isApprox(results[1]);
bool pred_3 = Z.isApprox(results[2]);
}
```
 Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
  advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at f.zapata@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org
