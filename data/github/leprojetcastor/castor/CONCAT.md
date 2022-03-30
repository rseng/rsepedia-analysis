# How to contribute

The preferred way to contribute to Castor is to fork the [main repository](https://github.com/leprojetcastor/castor) on Github.

1. *Search or open an [issue](https://github.com/leprojetcastor/castor/issues) for this project*:   
    * Search for and click into a related issue that you would like to contribute to. If the issue has not been assigned to anyone, leave a comment and assign the issue to yourself by clicking on the 'Assignees' button on the right.

    * If no such issues have been raised, open a new issue by clicking on the 'New issue' on the top right corner. Write a brief description of the issue and wait for the feedbacks from core developers before moving on to the next steps.

2. *Fork the [project repository](https://github.com/leprojetcastor/castor)*: Click on the 'Fork' button on the top right corner of the page. This creates a copy of the code under your account on the GitHub server. 

3. *Clone this copy to your local disk*: Open up a terminal on your computer, navigate to your preferred directory, and copy the following.
```
$ git clone git@github.com:<YourGithubUsername>/castor.git
$ cd castor
```

4. *Instantiate the project*: Install Castor and dependencies and compile the demos following the [user guide](https://leprojetcastor.gitlab.labos.polytechnique.fr/castor/installation.html).

5. *Create a branch to hold your changes*:
```
git checkout -b my-feature
```

6. Work on this copy on your computer using your preferred code editor such as Xcode, Visual Studio, VSCode, etc. Make sure you add the [corresponding tests] in the appropriate `demo_xxx.cpp` file. Use Git to do the version control. When you're done editing, do the following to record your changes in Git:
```
$ git add modified_files
$ git commit
```

7. *Push your changes* to Github with:
```
$ git push -u origin my-feature
```

8. Finally, go to the web page of your fork of the castor repository, and *click 'Pull request'* to send your changes to the maintainers for review. This will send an email to the committers.

If any of the above seems like magic to you, then look up the [Git documentation](https://git-scm.com/doc) on the web. It is recommended to check that your contribution complies with the following rule before submitting a pull request: all public methods should have informative docstrings with sample usage presented, compatible with the documentation process.

# THE CASTOR LIBRARY

The objective of the **castor** library is to propose **high-level semantics**, inspired by the Matlab language, allowing fast software prototyping in a low-level compiled language. It is nothing more than a matrix management layer using the tools of the **standard C++ library**, in different storage formats (full, sparse and hierarchical). Indeed, the use of IDEs 1 such as Xcode, Visual studio, Eclipse, etc. allows today to execute compiled code (C, C++, fortran, etc.) with the same **flexibility as interpreted languages** (Matlab, Python, Julia, etc.).

Full documentation is available here : http://leprojetcastor.gitlab.labos.polytechnique.fr/castor/index.html

Presentation in JOSS (Journal of Open Source Software) : [![DOI](https://joss.theoj.org/papers/10.21105/joss.03965/status.svg)](https://doi.org/10.21105/joss.03965)
---
title: 'Castor: A C++ library to code "à la Matlab"'
tags:
  - C++
  - Scientific computing
  - Fast prototyping
  - FEM/BEM simulation
authors:
  - name: Matthieu Aussal^[matthieu.aussal\@polytechnique.edu] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-2812-7578
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Marc Bakry^[marc.bakry\@polytechnique.edu] # note this makes a footnote saying 'co-first author'
    affiliation: 1
  - name: Laurent Series^[laurent.series\@polytechnique.edu]
    affiliation: 2
affiliations:
 - name: Ecole Polytechnique (CMAP), INRIA, Institut Polytechnique Paris, Route de Saclay 91128, Palaiseau, France
   index: 1
 - name: Ecole Polytechnique (CMAP), Institut Polytechnique Paris, Route de Saclay 91128, Palaiseau, France
   index: 2
date: 22 October 2021
bibliography: paper.bib
---

# Summary

The objective of the *Castor* framework is to propose high-level semantics, inspired by the Matlab language, allowing fast software prototyping in a low-level compiled language. It is nothing more than a matrix management layer using the tools of the standard C++ library (C++14 and later), in different storage formats (full, sparse and hierarchical). Linear algebra operations are built over the BLAS API and graphic rendering is performed in the VTK framework. The *Castor* framework is provided as an open source software under the LGPL 3.0, compiled and validated with clang and gcc. 

# Statement of need

Matlab is a software used worldwide in numerical prototyping, due to its particularly user-friendly semantics and its certified toolboxes. However, many use cases do not allow codes in Matlab format, for example multi-platform portability issues, proprieraty licensing and more generally code interfacing. To start meeting these needs, a header-only template library for matrix management has been developed, based on the standard C++14 and later library, by encapsulating the `std::vector` class. Many tools and algorithms are provided to simplify the development of prototypes:
 
 - dense, sparse and hierarchical matrices manipulations,
 - linear algebra computations [@blaslapack:00],
 - graphical representations [@vtk:2000].

This high-level semantic/low-level language coupling makes it possible to gain efficiency in the developpement, while ensuring performance for applications. In addition, direct access to data structures allows users to optimize the most critical parts of their code. Finally, a complete documentation is available, as well as continuous integration unit tests. All of this makes it possible to meet the needs of teaching (notebooks using a C++ interpreter such as Cling), academic research and industrial applications at the same time. 

# State of the field

For a developer accustomed to the Matlab language, it is natural to turn to prototyping tools such as Numpy or Julia, to produce open-source codes. Indeed, these tools today offer similar semantics and performance, with well-established user communities. To illustrate this similarity, the following codes perform the same tasks, with one implementation in Matlab [@MATLAB:2010] (left) and another in Julia [@bezanson2012julia] (right) :

| Matlab                            |     |     |     | Julia                                            |
| --------------------------------- | --- | --- | --- | ------------------------------------------------ |
|                                   |     |     |     | `using LinearAlgebra`                            |
| `tic`                             |     |     |     | `function test()`                                |
|                                   |     |     |     |                                                  |
| `M = [ 1  2  3 ;`                 |     |     |     | `M = [ 1  2  3 ;`                                |
| `      4  5  6 ;`                 |     |     |     | `      4  5  6 ;`                                |
| `      7  8  9 ;`                 |     |     |     | `      7  8  9 ;`                                |
| `     10 11 12];`                 |     |     |     | `     10 11 12];`                                |
| `disp(M);`                        |     |     |     | `display(M);`                                    |
| `M = (M - 1) .* ...`              |     |     |     | `M = (M .- 1) .* `                               |
| `    eye(size(M));`               |     |     |     | `    Matrix(I,size(M));`                         |
| `M(1,1) = -1;`                    |     |     |     | `M[1,1] = -1;`                                   |
| `M([2,3],1)  = -1;`               |     |     |     | `M[[2 3],1] .= -1;`                              |
| `M(4,:) = -1;`                    |     |     |     | `M[4,:] .= -1;`                                  |
| `disp(M);`                        |     |     |     | `display(M);`                                    |
|                                   |     |     |     |                                                  |
| `disp(sum(M,2));`                 |     |     |     | `display(sum(M,dims=2));`                        |
| `disp(abs(M));`                   |     |     |     | `display(abs.(M))`                               |
| `disp(sort(M,1));`                |     |     |     | `display(sort(M,dims=1));`                       |
| `disp(M*M');`                     |     |     |     | `display(M*M');`                                 |
|                                   |     |     |     |                                                  |
|                                   |     |     |     | `end`                                            |
|                                   |     |     |     |                                                  |
| `toc`                             |     |     |     |  `@time test();`                                 |
| `disp("done.");`                  |     |     |     |  `display("done.");`                             |                          

Despite the many advantages that these languages have and their high popularity, many codes are still developed natively in Fortran, C, and C++, for practical or historical reasons. Even if there are tools to automatically generate C/C++ code from a high-level language (as *Matlab Coder*), this work is often done manually by specialists. To find high-level semantics in native C++, we can turn to libraries like Eigen [@eigenweb], which offers a matrix API and efficient algebra tools. However, as the comparison below shows, the transcription from a Matlab code to an Eigen-based C++ code is not immediate: 

```c++
#include <iostream>
#include <chrono>
#include <eigen-3.4.0/Eigen/Dense>
using namespace std::chrono;
using namespace Eigen;
int main()
{
    auto tic = high_resolution_clock::now();
    
    MatrixXd  M(4,3);
    M << 1,  2,  3,
        4,  5,  6,
        7,  8,  9,
       10, 11, 12;
    std::cout << M << std::endl;
    
    M.array() -= 1;
    M.array() *= MatrixXd::Identity(M.rows(),M.cols()).array();
    M(0,0)     = -1;
    M.block<2,1>(1,0) = -MatrixXd::Ones(2,1);
    M.row(3)   = -MatrixXd::Ones(1,3);
    std::cout << M << std::endl;

    std::cout << M.rowwise().sum() << std::endl;
    std::cout << M.array().abs() << std::endl;    
    MatrixXd  Ms = M;
    for(auto col : Ms.colwise())
    {
      std::sort(col.begin(), col.end());
    }
    std::cout << Ms << std::endl;
    std::cout << M * M.transpose() << std::endl;
    
    auto toc  = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(toc - tic);
    std::cout << "Elapsed time: " << duration.count()*1e-6 << std::endl;
    std::cout << "done." << std::endl;

    return 0;
}
```
To complete this example, other references are available on this [link](https://eigen.tuxfamily.org/dox/AsciiQuickReference.txt). This is why all the features of the Castor library have been designed and developed so that the semantics at user level are as close to Matlab as what C++ allows. Moreover, to gain in portability, the manipulations of full matrices and the main algorithms depend only on the standard library which is available on the most majority of operating systems (MacOS, Linux, Windows, Android, etc.). Only advanced linear algebra tools require an external BLAS / LAPACK API, as well as graphical visualization functionality (VTK). The example below illustrates this goal:

| Matlab                            |     |     |     | Castor                            |
| -------------------------------   | --- | --- | --- | --------------------------------- |
|                                   |     |     |     | `#include "castor/matrix.hpp"`    |
|                                   |     |     |     | `using namespace castor;`         | 
|                                   |     |     |     | `int main (int argc, char* argv[])`|
|                                   |     |     |     | `{`                               |
| `tic`                             |     |     |     | `tic();`                          |
|                                   |     |     |     |                                   |
| `M = [ 1  2  3 ;`                 |     |     |     | `M = {{ 1,  2,  3},`              |
| `      4  5  6 ;`                 |     |     |     | `     { 4,  5,  6},`              |
| `      7  8  9 ;`                 |     |     |     | `     { 7,  8,  9},`              |
| `     10 11 12];`                 |     |     |     | `     {10, 11, 12}};`             |
| `disp(M);`                        |     |     |     | `disp(M);`                        |
|                                   |     |     |     |                                   |
| `M = (M - 1) .* eye(size(M));`    |     |     |     | `M = (M - 1) * eye(size(M));`     |
| `M(1,1) = -1;`                    |     |     |     | `M(0,0) = -1;`                    |
| `M([2,3],1)  = -1;`               |     |     |     | `M({1,2},0) = -1;`                |
| `M(4,:) = -1;`                    |     |     |     | `M(3,col(M)) = -1;`               |
| `disp(M);`                        |     |     |     | `disp(M);`                        |
|                                   |     |     |     |                                   |
| `disp(sum(M,2));`                 |     |     |     | `disp(sum(M,2));`                 |
| `disp(abs(M));`                   |     |     |     | `disp(abs(M));`                   |
| `disp(sort(M,1));`                |     |     |     | `disp(sort(M,1));`                |
| `disp(M*M');`                     |     |     |     | `disp(mtimes(M,transpose(M)));`   |
|                                   |     |     |     |                                   |
| `toc`                             |     |     |     | `toc();`                          |
| `disp("done.");`                  |     |     |     | `disp("done.");`                  |
|                                   |     |     |     |                                   |
|                                   |     |     |     | `return 0;`                       |
|                                   |     |     |     | `}`                               |
**Note:**
It is important to specify that the Castor library is far from offering today all the functionalities offered by Matlab and its many toolboxes. 

# Dense Matrix  

The dense matrix part of the *Castor* framework implements its own templatized class `matrix<T>` in `matrix.hpp`, where `T` can be fundamental types of C++ as well as `std::complex`. This class is built over a `std::vector<T>` which holds the values [@cpp:14]. Note that the element of a matrix is stored in row-major order and that a vector is considered as a $1\times n$ or $n\times 1$ matrix.

The class `matrix<T>` provides many useful functions and operators such as:

- builders which can be used to initialize all coefficients (`zeros`, `ones`, `eye`, etc.),
- standard algorithms over data stored in matrices (`norm`, `max`, `sort`, `argsort`, etc.),
- mathematical functions which can be applied element-wise (`cos`, `sqrt`, `conj`, etc.),
- matrix manipulations like concatenate matrices in all dimensions, find the non-zero elements or transpose them, reshape size, etc.,
- standard C++ operators which have been overloaded and work element-wise (`+`,`*`,`!`,`&&`,etc.),
- values accessors and matrix views with linear and bi-linear indexing,
- elements of linear algebra, such as the matrix product (`mtimes` or `tgemm`) and linear system resolution (multi-right-hand-side `gmres`),
- many other tools to display elements (`<<`, `disp`), save and load elements from file (ASCII or binary), etc.

The API provides more than a hundred functions and is designed such that it should feel like using Matlab. For advanced users, direct access to the data stored in the `std::vector<T>` enables all or part of an algorithm to be optimized in native C++. 

This example displays the sum of two matrices with implicit cast :

```c++
#include <iostream>
#include "castor/matrix.hpp"
using namespace castor;
int main(int argc, char* argv[])
{
    matrix<float> A = {{ 1.0,  2.0,  3.0,  4.0},
                       { 5.0,  6.0,  7.0,  8.0},
                       { 9.0, 10.0, 11.0, 12.0}};
    matrix<double> B = eye(3,4);
    auto C = A + B;
    disp(C);
    return 0;
}
```
```
Matrix 3x4 of type 'd' (96 B):
    2.00000      2.00000      3.00000      4.00000  
    5.00000      7.00000      7.00000      8.00000  
    9.00000     10.00000     12.00000     12.00000  
```

# Linear Algebra

The linear algebra part of the framework, implemented in `linalg.hpp`, provides a set of useful functions to perform linear algebra computations by linking to optimized implementations of the BLAS and LAPACK standards [@blaslapack:00] (OpenBLAS, oneAPI MKL, etc.).

The BLAS part is a straightforward overlay of the C-BLAS type III API, which is compatible with row-major ordering.  This is achieved by a template specialization of the `tgemm` function, which allows optimized implementations to take control of the computation using `sgemm`, `dgemm`, `cgemm` and `zgemm`. Thanks to this interface, naive implementations proposed in `matrix.hpp` for dense matrix-products `mtimes` and `tgemm` may be improved in term of performance, especially for large matrices. 

The LAPACK part is a direct overlay over the Fortran LAPACK API, which uses a column ordering storage convention. This interface brings new high-level functionalities, such as a linear solver (`linsolve`), matrix inversion (`inv`, `pinv`), factorizations (`qr`, `lu`), the search for eigen or singular values decompositions (`eig` ,`svd`), aca compression (`aca`), etc. It uses templatized low-level functions following the naming convention close to the LAPACK one (like `tgesdd`, `tgeqrf`, etc.).

This example displays the product of $A$ and $A^{-1}$ :

```c++
#include <iostream>
#include "castor/matrix.hpp"
#include "castor/linalg.hpp"
using namespace castor;
int main (int argc, char* argv[])
{
    matrix<> A = rand(4);
    matrix<> Am1 = inv(A); 
    disp(mtimes(A,Am1));
    return 0;
}
```
```
Matrix 4x4 of type 'd' (128 B):
 1.0000e+00   1.0408e-16  -2.7756e-17  -5.5511e-17  
          0   1.0000e+00  -5.5511e-17   1.1102e-16  
          0  -2.2204e-16   1.0000e+00  -1.1102e-16  
-2.7756e-17            0            0   1.0000e+00  
```

**Note:**
The backslash operator (`\`) not being available, the `linsolve` function allows to solve linear systems with:

  - LU decomposition with partial pivoting and row interchanges for square matrices (`[sdcz]gesv`),
  - QR or LQ factorization for overdetermined or underdetermined linear systems (`[sdcz]gels`).

In addition, an iterative multi-right-hand-side solver `gmres` is available in `matrix.hpp`, without dependency on BLAS and LAPACK. 

# 2D/3D Visualization

The graphic rendering part, provided by `graphics.hpp`, features 2D/3D customizable plotting and basic mesh generation. It is based on the well-known VTK library [@vtk:2000]. Here again, the approach tries to get as close as possible to Matlab semantics.

First, the user creates a `figure`, which is a dynamic container of data to display. The `figure` class is composed of a `vtkContextView` class, providing a view with a default interactor style, renderer, etc. Then, graphic representations can be added to the figure, using functions like `plot`, `imagesc`, `plot3`, `mesh`, etc. Options are available to customize the display of the results, such as the plotting style, legend, colorbar and others basic stuff. Finally, the `drawnow` function must be called to display all defined figures. The latters are displayed and manipulated in independent windows.

In addition, graphics exports are available in different compression formats (`png`,` jpg`, `tiff`, etc.), as well as video rendering (`ogg`). 

This example shows a basic 2D plotting of a sine function (\autoref{fig:sin}):

```c++
#include "castor/matrix.hpp"
#include "castor/graphics.hpp" 
using namespace castor;
int main (int argc, char* argv[])
{
    matrix<> X = linspace(0,10,100);
    figure fig;
    plot(fig,X,sin(X),{"r-+"},{"sin(x)"});
    plot(fig,X,cos(X),{"bx"},{"cos(x)"});
    drawnow(fig);
    return 0;
}
```
![Basic 2D plotting from Castor (using VTK).\label{fig:sin}](plot2d.png)

# Sparse matrices

Some matrices have sparse structures, with many (many) zeros that do not need to be stored [@sparse:1973]. There are adapted storage formats for this type of structure (LIL, COO, CSR, etc.), the most natural being to store the indices of rows and columns for each non-zero value, as a list of triplet $\{i,j,v\}$. For the *Castor* framework, a dedicated template class to this kind of matrix has been developed (see `smatrix.hpp`). The storage format is based on a row major sorted linear indexing. Only non-zero values and their sorted linear indices are stored in a list of pairs $\{v,l\}$:  for a $m\times n$ matrix, the following bijection is used to switch with the common bilinear indexation:

$$
\begin{aligned}
\{i,j\} &\rightarrow l = i \cdot n + j, \\
l &\rightarrow \{i=\frac{l}{n}; j= i\textrm{ mod }n\}.
\end{aligned}
$$

Accessors to all the elements are provided so that sparse matrices can be manipulated in a similar way as the dense matrices. This operation is performed by dichotomy with a convergence in $\log_2 (\text{nnz})$, where $\text{nnz}$ is the number of non-zero elements. Just like dense matrices, numerical values are stored in a templatized `std::vector<T>`. For convenience, we provide classical builders (`sparse`, `speye`, `spdiags`, etc.), standard C++ operators overloading, views, display functions (`disp`, `spy`) and some linear algebra tools (`transpose`, `mtimes`, `gmres`, etc.). 

This example displays the sum of two sparse matrices, with implicit cast and sparse to dense conversion :

```c++
#include <iostream>
#include "castor/smatrix.hpp"
using namespace castor;
int main (int argc, char* argv[])
{
    smatrix<float> As = {{0.0,  0.0,  0.0},
                         {5.0,  0.0,  7.0}};
    As(0,1) = 2.0;
    smatrix<double> Bs = speye(2,3);
    disp(As);
    disp(As(0,1)); // bilinear accessor
    disp(As(4));   // linear accessor
    disp(Bs);
    disp(full(As+Bs));
    return 0;
}
```
```
Sparse matrix 2x3 of type 'f' with 3 elements (12 B):
(0,1)  2
(1,0)  5
(1,2)  7
2
0
Sparse matrix 2x3 of type 'd' with 2 elements (16 B):
(0,0)  1
(1,1)  1
Matrix 2x3 of type 'd' (48 B):
    1.00000      2.00000            0  
    5.00000      1.00000      7.00000  
```

# Hierarchical matrices

To widen the field of applications, the $\mathcal H$-matrix format, so-called hierachical matrices [@hackbush:1999], have been added in `hmatrix.hpp`. They are specially designed for matrices with localized rank defaults. It allows a fully-populated matrix to be assembled and stored in a lighter format by compressing some parts of the original dense matrix using a low-rank representation [@rjasanow:2002]. They are constructed by binary tree subdivisions in a recursive way, with a parallel assembly of the compressed and full leaves (using the OpenMP standard). This format features a complete algebra, from elementary operations to matrix inversion. An example is given in the application section that follows.

# Application with a FEM/BEM simulation
As an application example, an acoustical scattering simulation was carried out using a boundary element method (BEM) tool, implemented with the *Castor* framework (see the *fembem* package [@fembem:21]). We consider a smooth $n$-oriented surface $\Gamma$ of some object $\Omega$, illuminated by an incident plane wave $u_i$ with wave-number $k$. The scattered field $u$ satisfies the Helmholtz equation in $\Omega$, Neumann boundary conditions (*sound-hard*) and the Sommerfeld radiation condition: 

$$
\begin{aligned}
-(\Delta u + k^2 u) &= 0 \\
-\partial_n u_i     &= 0 \\
\lim\limits_{r \to + \infty} r\textrm{ }(\partial_r u - i k u) &= 0
\end{aligned}
$$
  
The scattered field $u$ satisfies the integral representation (Neumann interior extension, see [@terasse:2007]):

$$
\label{eq1}\tag{1}
u(\textbf{x}) = - \left( \frac{1}{2}\mu(\textbf{x}) + \int_\Gamma \partial_{n_y} G(\textbf{x},\textbf{y})\mu(\textbf{y}) d_y \right) \quad \forall  \textbf{x} \in \Gamma^+,
$$

for some density $\mu$, with the Green kernel $G(\textbf{x},\textbf{y}) = \displaystyle\frac{e^{i k |x - y|}}{4 \pi |x - y| }$. Using the boundary conditions we obtain :

$$ 
\label{eq2}\tag{2}
- H\mu(\textbf{x})  = - \partial_n u_i(\textbf{x}) \quad \forall \textbf{x} \in \Gamma,
$$

where the hypersingular operator $H$ is defined by:

$$ 
\label{eq3}\tag{3}
H \mu(\textbf{x}) = \int_\Gamma \partial_{n_x} \partial_{n_y} G(\textbf{x},\textbf{y})\mu(\textbf{y}) d_y.
$$

The operator $H$ is assembled using a $P_1$ finite element discretization on a triangular mesh of the surface $\Gamma$, stored using dense matrices (`matrix.hpp`) or hierarchical matrices (`hmatrix.hpp`).

![Resonance mode at 8kHz of the human pinna (BEM with H-Matrix).\label{fig:head}](head.png)

Finally, using all the tools provided by Castor to write and solve these equations, we are able to efficiently compute the acoustic diffraction of a harmonic plane wave at 8kHz, on a human head mesh [@symare:2013]. As shown in \autoref{fig:head}, the simulation result highlights the role of the auditory pavilion as a resonator, modifying the timbre of a sound source to allow a listener's brain to precisely locate its direction. 

```c++
#include <castor/matrix.hpp>
#include <castor/smatrix.hpp>
#include <castor/hmatrix.hpp>
#include <castor/linalg.hpp>
#include <castor/graphics.hpp>
#include "fem.hpp"
#include "bem.hpp"
using namespace castor;
int main (int argc, char* argv[])
{
    // Load meshes
    matrix<double> Svtx;
    matrix<size_t> Stri;
    std::tie(Stri,Svtx) = triread("./","Head03_04kHz.ply");
    
    // Graphical representation
    figure fig;
    trimesh(fig,Stri,Svtx);
    
    // Parameters
    matrix<double> U = {0,0,-1};
    double f = 2000;
    double k = 2*M_PI*f/340;
    float tol = 1e-3;
    
    // FEM and mass matrix, sparse storage
    tic();
    femdata<double> v(Stri,Svtx,lagrangeP1,3);
    femdata<double> u(Stri,Svtx,lagrangeP1,3);
    auto Id = mass<std::complex<double>>(v);
    toc();
    
    // Left hand side '-H', equation (3), H-Matrix storage
    tic();
    auto LHSfct = [&v,&u,&k](matrix<std::size_t> Ix, matrix<std::size_t> Iy)
    {
        return -hypersingular<std::complex<double>>(v,u,k,Ix,Iy);
    };
    hmatrix<std::complex<double>> LHS(v.dof(),u.dof(),tol,LHSfct);
    toc();
    disp(LHS);
    
    // Right hand side '-dnUi', equation (2), full storage 
    auto B = - rightHandSide<std::complex<double>>(v,dnPWsource,U,k);
    
    // Solve '-H = -dnUi', equation (2), H-matrix preconditionner
    // associated to iterative solver
    hmatrix<std::complex<double>> Lh,Uh;
    tic();
    std::tie(Lh,Uh) = lu(LHS,1e-1);
    toc();
    disp(Lh);
    disp(Uh);
    auto mu = gmres(LHS,B,tol,100,Lh,Uh);
    
    // Boundary radiation, equation (1)
    tic();
    auto Dbndfct = [&v,&u,&k,&Id](matrix<std::size_t> Ix, matrix<std::size_t> Iy)
    {
        return 0.5*eval(Id(Ix,Iy)) + doubleLayer<std::complex<double>>(v,u,k,Ix,Iy);
    };
    hmatrix<std::complex<double>> Dbnd(v.dof(),u.dof(),tol,Dbndfct);
    matrix<std::complex<double>> Pbnd = mtimes(Dbnd,mu);
    toc();
    Pbnd = - gmres(Id,Pbnd,tol,100) + planeWave(v.dof(),U,k);
    
    // Graphical representation
    figure fig2;
    trimesh(fig2,Stri,Svtx,abs(Pbnd));
    
    // Export in .vtk file
    triwrite("./","head.vtk",Stri,Svtx,real(Pbnd));
    
    // Plot
    drawnow(fig);
    disp("done !");
    return 0;
}
```
  
# Acknowledgements

We thank Houssem Haddar for his precious help.

# References


.. _label-sparse-matrix:

Sparse matrix
=============

This library implements a lightweight class ``smatrix`` for the manipulation of sparse matrices, using solely the ``castor::matrix`` class without other external libraries. Many builders, functions and operators are available in a similar fashion as for the :ref:`label-class-matrix`, see the ``smatrix`` :ref:`label-smatrix-API`. Most of the manipulations should be transparent to the user and the feeling should be really `Matlab-like <https://www.mathworks.com>`_. Basic examples are given in the corresponding :ref:`label-examples-smatrix` section.

The ``smatrix`` class is installed automatically with the other headers of **castor** library.

**This is work-in-progress** 

.. _label-examples-smatrix:

Examples
........

In order to use a ``smatrix``, the corresponding header ``castor/smatrix.hpp`` needs to be included. A minimum working example would thus look like this:

.. code:: c++

    #include "castor/matrix.hpp"
    #include "castor/smatrix.hpp"

    int main()
    {
        smatrix<> s;
        // your code here
        return 0;
    }

Building a smatrix
++++++++++++++++++

A ``smatrix`` will mostly behave like a normal dense ``matrix``. For example, let us create a ``4 x 4`` ``smatrix`` of ``double``.

.. code:: c++

    smatrix<> As;
    disp(As,2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 0 elements (0 B):
    -empty-

The newly created ``smatrix`` is empty. It is rather easy to set the values of the entries. For example:

.. code:: c++

    As(1,2) = -0.5;
    disp(As,2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 1 elements (8 B):
    (1,2)  -0.5

Now ``As`` contains one non-zero element with the bilinear index ``(1,2)`` (second line, third column). However, this way of filling a ``smatrix`` **is definitely not recommended** as it involves a lot of memory management and may affect dramatically the performances. One way to do so is to first create the matrix in *triplet* format ``(I,J,VALUE)`` and only afterward build the corresponding ``smatrix`` as illustrated below:

.. code:: c++

    matrix<std::size_t> I = {0,2,3};
    matrix<std::size_t> J = {1,1,2};
    matrix<> V = {0.5, -1/3., M_PI};

    As = smatrix<>(I,J,V,4,4);
    disp(As,2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 3 elements (24 B):
    (0,1)  0.5
    (2,1)  -0.333333
    (3,2)  3.14159

Yes, we reaffected ``As`` to a new ``smatrix``. The old data is automatically discarded so one should be careful when performing such an action. As for ``matrix``, it is possible to :ref:`label-clear-smatrix` the content of a ``smatrix`` (the object is reinitialized). Let us now add a ``0``.

.. code:: c++

    As(3,3) = 0.; // but, why ?
    disp(As,2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 4 elements (32 B):
    (0,1)  0.5
    (2,1)  -0.333333
    (3,2)  3.14159
    (3,3)  0

A zero value is added. In order to clean a ``smatrix``, a simple call to ``check`` is sufficient:

.. code:: c++

    check(As);
    disp(As,2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 3 elements (24 B):
    (0,1)  0.5
    (2,1)  -0.333333
    (3,2)  3.14159

Everything went back to normal! Now, let us use one of the builders in order to create an identity sparse matrix. It is also possible to convert back to the *triplet* format.

.. code:: c++

    auto Bs = speye<>(4,4);
    disp(Bs,2);
    matrix<std::size_t> IB,JB;
    matrix<> VB;
    std::tie(IB,JB,VB) = find(Bs);
    disp(transpose(vertcat(vertcat(IB,JB),VB)),2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 4 elements (32 B):
    (0,0)  1
    (1,1)  1
    (2,2)  1
    (3,3)  1
    Matrix 4x3 of type 'd' (96 B):
              0            0      1.00000  
        1.00000      1.00000      1.00000  
        2.00000      2.00000      1.00000  
        3.00000      3.00000      1.00000 

The matrices ``IB,JB,VB`` are returned as *line* vectors. To obtain a better display, we concatenated them vertically and tranposed the result.


Basic manipulations
+++++++++++++++++++

In this section, we start with start from scratch so everything written in the previous section should be discarded from your ``main`` function. Let us create two matrices

.. code:: c++

    smatrix<> As = speye(4,4);
    
    matrix<std::size_t> I({1,1,2,2}), J({1,2,1,2});
    matrix<> V({1.,1.,1.,1.});
    smatrix<> Bs = smatrix<>(I,J,V,4,4);

``As`` is a ``4 x 4`` identity matrix and ``Bs`` is a matrix with the interior filled with ones. Here is an example of basic manipulations:

.. code:: c++

    auto Cs = 1.5*As - Bs/2.;
    disp(Cs,2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 6 elements (48 B):
    (0,0)  1.5
    (1,1)  1
    (1,2)  -0.5
    (2,1)  -0.5
    (2,2)  1
    (3,3)  1.5

What is the number of non-zero elements, again ?

.. code:: c++

    std::cout << "nnz(Cs) = " << nnz(Cs) << std::endl;

.. code:: text

    nnz(Cs) = 6

It is possible to get the value of any entry:

.. code:: c++ 

    std::cout << "Cs(1,2) = " << Cs(1,2) << std::endl;
    std::cout << "Cs(1,3) = " << Cs(1,3) << std::endl;

.. code:: text

    Cs(1,2) = -0.5
    Cs(1,3) = 0

Now, let us multiply ``Cs`` by a ``4 x 4`` dense ``matrix``:

.. code:: c++

    auto D = mtimes(Cs,ones<>(4));
    disp(D,2);  // :)

.. code:: text

    Matrix 4x4 of type 'd' (128 B):
        1.50000      1.50000      1.50000      1.50000  
        0.50000      0.50000      0.50000      0.50000  
        0.50000      0.50000      0.50000      0.50000  
        1.50000      1.50000      1.50000      1.50000

One last manipulation and we are good for this example.

.. code:: c++

    Cs(0,3) = M_PI;
    auto Es = Cs - transpose(Cs);
    check(Es); // drop the zeros
    disp(Es,2);
    disp(full(Es,2));

.. code:: text

    Sparse matrix 4x4 of type 'd' with 2 elements (16 B):
    (0,3)  3.14159
    (3,0)  -3.14159
    Matrix 4x4 of type 'd' (128 B):
              0            0            0      3.14159  
              0            0            0            0  
              0            0            0            0  
       -3.14159            0            0            0
Boids
=====

*Shared by Antoine Rideau*

| On this page you will find how to simulate using **Castor** the flocking behaviour of birds.
| Boids - contraction of "bird-oid object" - refers to an artificial life program, developed by Craig Reynolds in 1986, simulating the flocking behaviour of birds.
|
| The complex flocking behaviours arise from three simple rules describing the interaction between each boid:
|       **Separation:** steer to avoid crowding local flockmates.
|       **Alignement:** steer towards the average heading of local flockmates.
|       **Cohesion:** steer to move toward the average position of local flockmates.




The flock is composed of ``N`` boids and will fly during ``nt`` steps within a perimeter described by ``dimension``

.. code-block:: c++

    // Parameters
    int N = 100;  // Number of boids
    int nt = 300; // Number of time steps
    matrix<> dimension = {50, 50}; // Permieter lengths


| Each boids in the flock is characterized by its position :math:`(x,y)` and its velocity vector :math:`\overrightarrow{v}=(dx,dy)` . Each parameter is stored in a ``matrix<>`` and then gathered inside the ``std::vector`` ``Boids``.
| Furthermore the matrix ``Out`` contains these parameters at each steps in order to output them later on for a satisfaying visualization.

.. code-block:: c++

    // Initialization
    std::vector<matrix<>> Boids = { rand(1, N,true) * dimension(0),             // Boids x position
                                    rand(1, N) * dimension(1),                  // Boids y position
                                    rand(1, N,true) * maxSpeed - maxSpeed / 2,  // Boids x velocity : dx
                                    rand(1, N) * maxSpeed - maxSpeed / 2};      // Boids y velocity : dy
    matrix<> Out = cat(1, Boids[0], cat(1, Boids[1], cat(1, Boids[2], Boids[3])));

See :ref:`label-rand` , :ref:`label-cat` .


.. figure:: img/boids.png
    :width: 200
    :align: center
    :figclass: align-center
    
    Neighbourhood of a boid (blue) : zone of separation (red), zone of alignement (yellow) and zone of cohesion (green).

So as to apply the set of rules to each boid, it is needed to determine this boid's flockmates depending on the distance between this boid and the others boids.
This is done by the function ``boidsWithinDistance`` which returns the index of the boids who are inside the circle of ``distance`` radius  around the boid of index ``boidIndex``.

.. code-block:: c++

    matrix<size_t> boidsWithinDistance(std::vector<matrix<>> &Boids, int boidIndex, double distance)
    {
        matrix<size_t> I(boidIndex);
        return setdiff(find(sqrt((Boids[0](boidIndex) - Boids[0]) * (Boids[0](boidIndex) - Boids[0]) + (Boids[1](boidIndex) - Boids[1]) * (Boids[1](boidIndex) - Boids[1])) < distance), I);
    }

See :ref:`label-setdiff` , :ref:`label-find` .


Behaviour
---------

Borders
^^^^^^^

| There are two boundaries conditions that can be applied :

1. Boids must stay within the borders.

    .. code-block:: c++

        // Constrain boids to within the borders. If it gets too close to an edge,nudge it back in and reverse its direction.
        void stayWithinBorders(std::vector<matrix<>> &Boids, matrix<> border, double margin, double turnFactor)
        {
            auto xTooLow = find(Boids[0] < margin);
            auto xTooHigh = find(Boids[0] > border(0) - margin);
            auto yTooLow = find(Boids[1] < margin);
            auto yTooHigh = find(Boids[1] > border(1) - margin);

            if (numel(xTooLow) > 0){Boids[2](xTooLow) = eval(Boids[2](xTooLow)) + turnFactor;}
            if (numel(xTooHigh) > 0){Boids[2](xTooHigh) = eval(Boids[2](xTooHigh)) - turnFactor;}
            if (numel(yTooLow) > 0){Boids[3](yTooLow) = eval(Boids[3](yTooLow)) + turnFactor;}
            if (numel(yTooHigh > 0)){Boids[3](yTooHigh) = eval(Boids[3](yTooHigh)) - turnFactor;}
        }

See :ref:`label-find` , :ref:`label-numel` , :ref:`label-view`  .

2. Like Pacman, boids who go out by one side are wrap around to the other side.

    .. code-block:: c++

        //  Checks if boids go out of the window and if so, wraps them around to the other side.
        void wrapBorders(std::vector<matrix<>> &Boids, matrix<> border)
        {
            auto xTooLow = find(Boids[0] < 0);
            auto xTooHigh = find(Boids[0] > border(0));
            auto yTooLow = find(Boids[1] < 0);
            auto yTooHigh = find(Boids[1] > border(1));

            if (numel(xTooLow) > 0){Boids[0](xTooLow) = eval(Boids[0](xTooLow)) + border(0);}
            if (numel(xTooHigh) > 0){Boids[0](xTooHigh) = eval(Boids[0](xTooHigh)) - border(0);}
            if (numel(yTooLow) > 0){Boids[1](yTooLow) = eval(Boids[1](yTooLow)) + border(1);}
            if (numel(yTooHigh > 0)){Boids[1](yTooHigh) = eval(Boids[1](yTooHigh)) - border(1);}
        }

See :ref:`label-find` , :ref:`label-numel` , :ref:`label-view`  .


Separation
^^^^^^^^^^

For each boid, once flockmates within ``separationDistance`` are known, velocity vector is ajusted toward the opposite directions of flockmates' positions. 

.. code-block:: c++

    // Boids try to keep a small distance away from other boids
    void separation(std::vector<matrix<>> &Boids, double separationDistance, double separationFactor)
    {
        for (int i = 0; i < numel(Boids[0]); i++)
        {
            auto FlockMates = boidsWithinDistance(Boids, i, separationDistance);
            if (numel(FlockMates) > 0)
            {
                Boids[2](i) += sum(Boids[0](i) - eval(Boids[0](FlockMates))) * separationFactor;
                Boids[3](i) += sum(Boids[1](i) - eval(Boids[1](FlockMates))) * separationFactor;
            }
        }
    }

See :ref:`label-sum` , :ref:`label-view`  .



Alignement
^^^^^^^^^^

For each boid, once flockmates within ``alignementDistance`` are known, velocity vector is ajusted toward the average direction of flockmates' velocity vector. 

.. code-block:: c++

    // Boids try to match velocity with near boids
    void alignement(std::vector<matrix<>> &Boids, double alignementDistance, double alignementFactor)
    {
        for (int i = 0; i < numel(Boids[0]); i++)
        {
            auto FlockMates = boidsWithinDistance(Boids, i, alignementDistance);
            if (numel(FlockMates) > 0)
            {
                Boids[2](i) += (sum(eval(Boids[2](FlockMates))) / numel(FlockMates)) * alignementFactor;
                Boids[3](i) += (sum(eval(Boids[3](FlockMates))) / numel(FlockMates)) * alignementFactor;
            }
        }
    }

See :ref:`label-sum` , :ref:`label-numel` , :ref:`label-view`  .

Cohesion
^^^^^^^^

For each boid, once flockmates within ``cohesionDistance`` are known, velocity vector is ajusted toward the average position of flockmates.

.. code-block:: c++

    // Boids try to fly towards the centre of mass of neighbouring boids
    void cohesion(std::vector<matrix<>> &Boids, double cohesionDistance, double cohesionFactor)
    {    
        for (int i = 0; i < numel(Boids[0]); i++)
        {
            auto FlockMates = boidsWithinDistance(Boids, i, cohesionDistance);
            if (numel(FlockMates) > 0)
            {
                Boids[2](i) += ((sum(eval(Boids[0](FlockMates))) / numel(FlockMates)) - Boids[0](i)) * cohesionFactor;
                Boids[3](i) += ((sum(eval(Boids[1](FlockMates))) / numel(FlockMates)) - Boids[1](i)) * cohesionFactor;
            }
        }
    }

See :ref:`label-sum` , :ref:`label-numel` , :ref:`label-view`  .

Speed limitation
^^^^^^^^^^^^^^^^

Once the rules are applied on the boids, their speed is limited to ``maxSpeed`` .

.. code-block:: c++

    void limitSpeed(std::vector<matrix<>> &Boids, double maxSpeed)
    {
        auto Speed = sqrt(Boids[2] * Boids[2] + Boids[3] * Boids[3]);
        matrix<> I = find(Speed > maxSpeed); //Index where velocity > maxSpeed

        if (numel(I) > 0)
        {
            Boids[2](I) = (eval(Boids[2](I)) / eval(Speed(I))) * maxSpeed;
            Boids[3](I) = (eval(Boids[3](I)) / eval(Speed(I))) * maxSpeed;
        }
    }

See :ref:`label-sqrt` , :ref:`label-find` , :ref:`label-view`  .

Visualisation
-------------

| In order to have a more satisfaying animation of the flying boids Python will be used.
| Beforehand, the data stored in ``Out`` are written in text file.

.. code-block:: c++

    writetxt("../", "data.txt", Out);


.. code-block:: python

    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    import matplotlib.patches as patches
    import numpy as np

    data = np.loadtxt("./data.txt", skiprows=1)
    # print(data)

    nt = int(data.shape[0]/4)
    N = int(data.shape[1])

    xs = [data[i*4, :] for i in range(nt)]
    ys = [data[i*4+1, :] for i in range(nt)]
    dxs = [data[i*4+2, :] for i in range(nt)]
    dys = [data[i*4+3, :] for i in range(nt)]

    # Visu initialization
    fig, ax = plt.subplots()
    p = patches.Rectangle((0, 0), 50, 50, fill=True)
    ax.add_artist(p)
    ax.set_xlim(0, 50)
    ax.set_ylim(0, 50)
    ax.set_aspect('equal')

    # Length of velocity vector to then normalize its representation
    D = [np.sqrt((dxs[0][i])**2+(dys[0][i])**2) for i in range(N)]
    boids = [ax.annotate("", xy=(xs[0][i]+dxs[0][i]/D[i], ys[0][i] + dys[0][i]/D[i]), xytext=(xs[0][i], ys[0][i]), arrowprops={
                         "facecolor": "red", 'arrowstyle': 'wedge'})for i in range(N)]


    def animate(frame):
        for i in range(N):
            d = np.sqrt((dxs[frame][i])**2+(dys[frame][i])**2)
            pos = np.array([xs[frame][i], ys[frame][i]])

            boids[i].set_position(pos)
            boids[i].xy = pos + (dxs[frame][i]/d, dys[frame][i]/d)

        return boids


    # Creating the Animation object
    ani = animation.FuncAnimation(
        fig, animate, nt, interval=40, blit=True, repeat_delay=1000)
    plt.show()

.. raw:: html

    <video autoplay controls width="100%">

    <source src="./_static/boidsflight.mp4"
            type="video/mp4">

    Sorry, your browser doesn't support embedded videos.
    </video>

|                           Flight of 100 boids in a 50-by-50 square.

Code
----

Here is all the code at once, without the functions written above :


.. code-block:: c++

    int main(int argc, char const *argv[])
    {
        // Parameters
        int N = 100;  // Number of boids
        int nt = 300; // Number of time steps
        matrix<> dimension = {50, 50};

        double margin = 2.;
        double turnFactor = 1.;
        double maxSpeed = 1.5;
        double separationDistance = 2.0;
        double separationFactor = 0.05;
        double alignementDistance = 7.5;
        double alignementFactor = 0.10;
        double cohesionDistance = 8.5;
        double cohesionFactor = 0.03;

        // Initialization
        std::vector<matrix<>> Boids = { rand(1, N,true) * dimension(0),              // Boids x position
                                        rand(1, N) * dimension(1),              // Boids y position
                                        rand(1, N,true) * maxSpeed - maxSpeed / 2,   // Boids x velocity : dx
                                        rand(1, N) * maxSpeed - maxSpeed / 2};  // Boids y velocity : dy
        matrix<> Out = cat(1, Boids[0], cat(1, Boids[1], cat(1, Boids[2], Boids[3])));

        disp(Out);

        tic();
        for (int t = 0; t < nt; t++)
        {
            separation(Boids, separationDistance, separationFactor);
            alignement(Boids, alignementDistance, alignementFactor);
            cohesion(Boids, cohesionDistance, cohesionFactor);
            limitSpeed(Boids, maxSpeed);
            stayWithinBorders(Boids, dimension, margin, turnFactor);

            // Update positions
            Boids[0] += Boids[2];
            Boids[1] += Boids[3];

            // wrapBorders(Boids, dimension);

            for (int i = 0; i < 4; i++)
            {
                Out = cat(1, Out, Boids[i]);
            }
        }
        toc();

        writetxt("../", "data.txt", Out);

        return 0;
    }



References  
----------

https://eater.net/boids

https://en.wikipedia.org/wiki/Boids

http://www.red3d.com/cwr/boids/

http://www.vergenet.net/~conrad/boids/pseudocode.html.. _label-graphical-io:

Graphical input/output
++++++++++++++++++++++

.. _label-triread:

triread
-------
.. doxygenfunction:: triread(std::string const &path, std::string const &name)
    :project: castor

See :ref:`label-triwrite`.


.. _label-triwrite:

triwrite
--------
.. doxygenfunction:: triwrite(std::string const &path, std::string const &name, matrix<std::size_t> const &tri, matrix<T> const &vtx, matrix<T> const &val = {})
    :project: castor

See :ref:`label-triread`.

.. _label-basic-plot:


Basic plot
++++++++++

These functions allow to display curves or values of a ``matrix``.



.. _label-imagesc:

imagesc
-------
.. doxygenfunction:: imagesc(figure &fig, matrix<T> const &M)
    :project: castor

.. _label-plot:

plot
----
.. doxygenfunction:: plot(figure &fig, matrix<T> const &X, matrix<T> const &Y, std::vector<std::string> const &style = {""}, std::vector<std::string> const &label = {""})
    :project: castor

See :ref:`label-plot3`

.. _label-plot3:

plot3
-----
.. doxygenfunction:: plot3(figure &fig, matrix<T> const &X, matrix<T> const &Y, matrix<T> const &Z, std::string const &style = "")
    :project: castor

See :ref:`label-plot`

.. _label-spy:

spy
---
.. doxygenfunction:: spy
    :project: castor.. _label-basic:

Basics
======

In this section, the user will find informations about basic manipulations on the ``matrix`` object : creation (and destruction), etc. ``matrix`` is defined within the namespace ``castor`` meaning that in order to call the functions, their name should be preceeded by ``castor::``. In order to use directly their name, the user should add 

.. code:: c++

  using namespace castor;

in the preamble of the ``.cpp`` file. A minimum ``main`` file would then look like this:

.. code:: c++

  #include "castor/matrix.hpp"

  using namespace castor;

  int main()
  {
    matrix<double> A;
    // Write your code here
    return 0;
  }

For some more advanced features, see the corresponding section :ref:`label-advanced`. 

**The user may also refer to the examples in the** ``demo/*`` **subdirectories of the main directory of the castor project.**


Matrix creation (and destruction)
---------------------------------

There are two main paths to create a ``matrix`` object. The first way is to call one of the constructors explicitly. If no value is specified, the object will be filled with zeros. The second way it to call a *builder* (see :ref:`label-basic-builder` and :ref:`label-basic-display`).

From matrix constructor
+++++++++++++++++++++++

**Remark:** In the following, we will use a lot the :ref:`label-disp` function which is meant to produce a formatted output of the content of a ``matrix`` object.

The code below creates an empty ``matrix`` of type ``int``.

.. code:: c++

  matrix<int> A;
  disp(A,2);

.. code:: text

  Matrix 0x0 of type 'i' (0 B):
  -empty-

**Remark:** If the ``using namespace castor`` clause was not added in the preamble as explained at the beginning of this section, the code above becomes

.. code:: c++

  castor::matrix<int> A;
  castor::disp(A,2);
  // etc.

By passing the value "2" as a second argument to ``disp``, we can see that it is a ``matrix`` of size ``0x0`` of integer type whose data size is 0 Byte.

By passing as argument a single value of type *T*, a singleton ``matrix`` is created.

.. code:: c++

  matrix<int> A(M_PI);
  disp(A,2);

.. code:: text

  Matrix 1x1 of type 'i' (4 B):
  3

Here, ``A`` has been declared as an ``matrix`` of integers but a ``double`` containing the value of pi was passed as argument. As a consequence, it was cast to an ``int``, thus the value 3. Note that 4 Bytes is the size of an integer in C++ when standard compilation options are used.

Next, we initialize a matrix using an initialization-list. By passing a single list of values as argument to the constructor, a line-``matrix`` is created. By passing a list of a list, a ``matrix`` is created whose number of lines is the number of elements in the outer list and the number of columns is the number of elements in the inner lists. Please be careful that the number of elements in the inner list should be the same for all. These two options are illustrated below.

.. code:: c++

  matrix<float> A({1,2,3,4,5});       // matrix of floats
  matrix<>      B({{1,2,3},{4,5,6}}); // matrix of doubles
  disp(A,2);
  disp(B,2);

.. code:: text

  Matrix 1x5 of type 'f' (20 B):
      1.00000      2.00000      3.00000      4.00000      5.00000  
  Matrix 2x3 of type 'd' (48 B):
      1.00000      2.00000      3.00000  
      4.00000      5.00000      6.00000


Finally, it is possible to create a ``matrix`` by giving its dimensions and a fill-value. By default, the matrix is filled with 0s. In the example below, we create a ``2x3`` matrix filled with the value 4, then we modify one of the entries.

.. code:: c++

  matrix<> A(2,3,4.);
  disp(A,2);
  A(1,2) = -0.5;
  disp(A,2);

.. code:: text

  Matrix 2x3 of type 'd' (48 B):
      4.00000      4.00000      4.00000  
      4.00000      4.00000      4.00000  
  Matrix 2x3 of type 'd' (48 B):
      4.00000      4.00000      4.00000  
      4.00000      4.00000     -0.50000


Please refer to the constructors list in the :ref:`class matrix description <label-class-matrix>`. 


.. _label-basic-builder:

From builder
++++++++++++

We describe now some of the useful builders for the ``matrix`` class.

The code below creates a ``2x3`` matrix of ``double`` filled with zeros.

.. code:: c++

    matrix<> A = zeros(2,3);
    disp(A);

.. code:: text

   1.0000  1.0000  1.0000  
   1.0000  1.0000  1.0000  

**Remark:** this is equivalent to calling explicitly the ``matrix`` constructor.

The code below creates a ``1x10`` matrix of ``double`` initialized with linear spaced values :

.. code:: c++

    matrix<> A = linspace(0,1,10);
    disp(A,2);

.. code:: text

    Matrix 1x10 of type 'd' (80 B):
              0      0.11111      0.22222  ...      0.77778      0.88889      1.00000

For a ``2x2`` random matrix of ``float``, use

.. code:: c++

    matrix<float> A = rand<float>(2);
    disp(A,2);

.. code:: text

    Matrix 2x2 of type 'f' (16 B):
        0.84019      0.39438  
        0.78310      0.79844

This last result may differ depending on your random number generator.

Notes : 

- Matrices and vectors are objects of the matrix template class. A vector is considered as a (1xn) size by default or (nx1). 
- The template argument of class matrix is double by default. It is possible to specify type both for matrix constructors and builders.


Clear a matrix
++++++++++++++

If for some reason the content of a ``matrix`` needs to be cleared (for example, free some RAM), there are to possibilities. The first solution (the *clean one*) is to call the :ref:`label-clear` function.

.. code:: c++

  auto A=rand(5);
  disp(A,2);
  clear(A);
  disp(A,2);

.. code:: text

  Matrix 3x3 of type 'd' (72 B):
      0.84019      0.39438      0.78310  
      0.79844      0.91165      0.19755  
      0.33522      0.76823      0.27777  
  Matrix 0x0 of type 'd' (0 B):
  -empty-

The second one is to assign an empty ``matrix`` in place of an existing one.

.. code:: c++

  auto A=rand(3);
  disp(A,2);
  A = {};
  disp(A,2);

.. code:: text

  Matrix 3x3 of type 'd' (72 B):
      0.84019      0.39438      0.78310  
      0.79844      0.91165      0.19755  
      0.33522      0.76823      0.27777  
  Matrix 0x0 of type 'd' (0 B):
  -empty-



.. _label-basic-display:

Display
-------

A very useful function is the :ref:`label-disp` function. It produces a formatted output of the content of a ``matrix`` object with additional informations. Let us create a ``2x2`` random ``matrix``.

.. code:: c++

  auto A = rand<>(2);

The simplest call to :ref:`label-disp` displays the raw content without additional informations.

.. code:: c++

  disp(A);

.. code:: text

      0.84019      0.39438  
      0.78310      0.79844

In fact, this is equivalent to calling ``disp(A,1)``. The second (optional) argument determines the level of informations to be displayed. ``disp(A,0)`` will produce the same output as ``disp(A)`` but no end-of-line character is added to the output. ``disp(A,2)`` will add informations on the size, the type and the RAM storage of the ``matrix``, as illustrated before.

The third argument to :ref:`label-disp` is the output stream (``std::ostream``) which by default is the standard output ``std::cout``. Finally, the user may specify two additional arguments which are the number of lines and columns which need to be displayed. By default, their value is 3 meaning that the first 3 and last 3 element of each direction are displayed.

.. code:: c++

  auto A = rand(10);
  disp(A,2);

.. code:: text

  Matrix 10x10 of type 'd' (800 B):
      0.84019      0.39438      0.78310  ...      0.76823      0.27777      0.55397  
      0.47740      0.62887      0.36478  ...      0.71730      0.14160      0.60697  
      0.01630      0.24289      0.13723  ...      0.10881      0.99892      0.21826  
  ...
      0.53161      0.03928      0.43764  ...      0.73853      0.63998      0.35405  
      0.68786      0.16597      0.44010  ...      0.89337      0.35036      0.68667  
      0.95647      0.58864      0.65730  ...      0.81477      0.68422      0.91097

Now we modifiy a little bit the format.

.. code:: c++
  
  disp(A,2,std::cout,4,2);

.. code:: text

  Matrix 10x10 of type 'd' (800 B):
      0.84019      0.39438  ...      0.27777      0.55397  
      0.47740      0.62887  ...      0.14160      0.60697  
      0.01630      0.24289  ...      0.99892      0.21826  
      0.51293      0.83911  ...      0.29252      0.77136  
  ...
      0.23828      0.97063  ...      0.51254      0.66772  
      0.53161      0.03928  ...      0.63998      0.35405  
      0.68786      0.16597  ...      0.35036      0.68667  
      0.95647      0.58864  ...      0.68422      0.91097

**Note :** :ref:`label-disp` can also display the content of other variables :

.. code:: c++

  disp("pi value is :");
  disp(M_PI);

.. code:: text

    pi value is
    3.14159


Size and indexing 
-----------------

We describe now a few useful functions to begin manipulating matrices.

Size, length, numel
+++++++++++++++++++

The **size** function returns the two-element vector containing the number of rows and columns in the matrix. The result is *also* a ``matrix``.

.. code:: c++

  matrix<> A = eye(3,4);
  disp(size(A),2)

.. code:: text

  Matrix 1x2 of type 'm' (16 B):
  3  4

If a dimension is specified, :ref:`label-size` returns only the length of the specified dimensions :

.. code:: c++

  matrix<> A = eye(3,4);
  disp(size(A, 1));
  disp(size(A, 2));

.. code:: text

   3
   4

The :ref:`label-length` and :ref:`label-numel` functions returns respectively the maximum length and the number of elements in the matrix :

.. code:: c++

  matrix<> A = eye(3,4);
  disp(length(A));
  disp(numel(A));

.. code:: text

   4
   12


Accessing elements
++++++++++++++++++

The elements of a ``matrix`` can be accessed using either *linear* or *bilinear* indexing. 

*Linear* indexing consists in accessing the n-th element of the ``matrix`` in the order the data is stored. Since ``matrix`` uses a row-major layout, the rows of the ``matrix`` are concatenated one after the other. 

.. code:: c++

  matrix<> A = {{1,2,3},
                {4,5,6},
                {7,8,9}};
  disp(A(5));

.. code:: text

  6

The 5-th element of ``A`` thus holds the value 6 (the 6-th holds 7, etc.). Linear indexing is particularly useful when accessing the elements of a one-dimensional ``matrix`` (a *vector*).

**Remark:** The index numbering follows the C/C++ convention meaning that the indexes start a ``0`` and ends at ``n-1`` where ``n`` would be a dimension of the ``matrix`` (see :ref:`label-size`).

*Bilinear* indexing is the natural way to access the elements of a ``matrix``.

.. code:: c++

  disp(A(1,2));

.. code:: text

  6


Help
----

The function :ref:`label-help` allows you to display the documentation of a function at runtime. You have to give the complete path to the header file ``matrix.hpp`` in the ``documentationFiles`` variable. 

.. code:: c++

  documentationFiles =
  {
      "/complete/path/to/matrix.hpp"
  };

  help("size");

.. code:: text

  ============================ DOCUMENTATION ============================
  Help on "size":
  Size of array.

  S = size(A) for m-by-n matrix A returns the two-element vector [m,n]
  containing the number of rows and columns in the matrix.

  S = size(A,dim) returns the lengths of the specified dimensions dim.

  Example(s):
     matrix<> A = {{1,2,3},{4,5,6}};
     disp(size(A));
     disp(size(A,1));
     disp(size(A,2));

  See also:
    length, numel.
  =======================================================================


Basic operations
----------------

The ``matrix`` class is designed to be as easy of use as Matlab or Numpy arrays. As a consequence, many operators have been overloaded. We will describe here some basic manipulations with few of all of the available operators. They can be discovered at :ref:`label-operators`.

First, we create two matrices ``A`` and ``B``, we multiply the first one by ``M_PI`` and we add them.

.. code:: c++

  // create two 'double' matrices
  auto A = eye(2);
  auto B = eye(2);
  disp(A,2);
  disp(B,2);

  A *= M_PI;
  disp(A,2);

  auto C = A + B;
  disp(C,2);

.. code:: text

  Matrix 2x2 of type 'd' (32 B):
      1.00000            0  
            0      1.00000  
  Matrix 2x2 of type 'd' (32 B):
      1.00000            0  
            0      1.00000  
  Matrix 2x2 of type 'd' (32 B):
      3.14159            0  
            0      3.14159  
  Matrix 2x2 of type 'd' (32 B):
      4.14159            0  
            0      4.14159
      
Then, we create an orthogonal matrix ``L`` and we compute ``D = C - L*L'*C``. Orthogonal matrices are such that the matrix product with their transpose yields the identity matrix. Therefore, the result should be a null-matrix.

.. code:: c++

  double theta = 0.2;
  matrix<> L = {
    {std::cos(theta),-std::sin(theta)},
    {std::sin(theta),std::cos(theta)}
  };

  auto D = C - mtimes(L,mtimes(transpose(L),C));
  disp(D,2);

.. code:: text

  Matrix 2x2 of type 'd' (32 B):
          0            0  
          0            0  


We obtain the expected result. **Note that the matrix-matrix product is computed using** :ref:`label-mtimes`. Indeed, the ``*`` operator *does not* compute the matrix-matrix product but a term-by-term product. 

.. code:: c++

  auto A = ones(2);
  auto B = matrix<>(2,2,2.);
  disp(A*B,2);

.. code:: text

  Matrix 2x2 of type 'd' (32 B):
      2.00000      2.00000  
      2.00000      2.00000
.. _label-matrix-dimensions:

Matrix dimensions
+++++++++++++++++

These functions allow to recover the dimensions of a matrix : total number of elements (see :ref:`label-numel`), dimensions (see :ref:`label-size`), etc.

.. _label-length:

length
------
.. doxygenfunction:: length(matrix<T> const &A)
   :project: castor

See :ref:`label-numel`, :ref:`label-size`.

.. _label-nnz:

nnz
---
.. doxygenfunction:: nnz(matrix<T> const &A)
   :project: castor

See :ref:`label-find`, :ref:`label-size`.

.. _label-numel:

numel
-----
.. doxygenfunction:: numel(matrix<T> const &A)
   :project: castor

See :ref:`label-length`, :ref:`label-size`.

.. _label-size:

size
----
.. doxygenfunction:: size(matrix<T> const &A, int dim)
   :project: castor

.. doxygenfunction:: size(matrix<T> const &A)
   :project: castor

See :ref:`label-length`, :ref:`label-numel`.

.. _label-class-bintree:

Class bintree
+++++++++++++

This class is used by :ref:`label-class-hmatrix` for space partitioning.

.. doxygenclass:: castor::bintree
   :project: castor
   :members:
   :undoc-members:

.. _label-class-figure:

Class figure
============

The ``figure`` class implements a container which is able to display data. See :ref:`label-graphics-advanced` for some examples of use.

.. doxygenclass:: castor::figure
    :project: castor
    :members: 
Heat equation
=============

*Shared by Antoine Rideau*


On this page you will find how to code using **Castor** an example as simple as heat diffusion in 1 dimension without heat source with Dirichlet boundary described by the following equation :


.. math:: 
    
    \left\{ \begin{matrix} \displaystyle \frac{\partial u }{\partial t} = d \frac{\partial^2 u}{\partial x^2} 
    \\ u(x_{min},t) = u(x_{max}, t) = 0 
    \\ u(x,0) = \sin(x\pi )^{2}  \end{matrix} \right.


The simulation is focuses on the interval :math:`\left [ x_{min}, x_{max}\right ]=\left [ 0,2 \right ]` for the space domain which is divided between ``Nx = 1000`` points .

.. code-block:: c++
    
    //Parameters
    int Nx = 1000;
    double xmin = 0;
    double xmax = 1;
    double tf = 0.1;
    double ti = 0;

    
The space is discretized with ``dx`` steps and the time with ``dt`` steps. 

.. math:: 

    \begin{matrix} x_{i} = i \Delta x \\ t_{n} = n \Delta t \end{matrix}

The ``X`` vector stores the space grid

.. math:: 
    
    X = \begin{pmatrix} x_{0}\\ x_{1} \\ \vdots \\ x_{N_{x}-1} \end{pmatrix} ,

and the vector ``U0`` contains the initial heat distribution 

.. math:: 

    U_{0} = \begin{pmatrix} u_{0}^{0} \\ u_{1}^{0} \\ \vdots \\ u_{N_{x}-1}^{0} \end{pmatrix} \text{ with } u_{i}^{0} = \sin(x_{i})^{2} \text{ for } i=0,...,N_{x} .


.. code-block:: c++

    //Discretization
    double dx = (xmax - xmin) / (Nx - 1);
    double dt = 0.5 * (dx * dx) / d;
    matrix<> X = linspace(xmin, xmax, Nx);
    matrix<> U0 = pow(sin(X * M_PI), 2);


The second space derivative is approximated by :

.. math:: 

    \frac{\partial^2 u_{i}^{n}}{\partial x^2}\approx \frac{u_{i+1}^{n}-2u_{i}^{n}+u_{i-1}^{n}}{\Delta x^{2}} 


.. Analytical solution
.. -------------------

.. With those parameters, the analytical solution is :

.. .. math::

..     \text{ MATHS }


Explicit Euler
--------------

The specifity of the explicit Euler scheme is that the time derivative is calculated with the forward differential quotient approximation

.. math::

    \frac{\partial u_{i}^{n}}{\partial t}\approx \frac{u_{i}^{n+1}-u_{i}^{n}}{\Delta t}

So as to maintain the stability of the simulation, the time step ``dt`` has to respect the following inequality

.. math:: 

    d \frac{\Delta t}{\Delta x^{2}} \leq \frac{1}{2}

Finally, the following numerical scheme is obtained :

.. math:: 

    \frac{u_{i}^{n+1}-u_{i}^{n}}{\Delta t}- d \frac{u_{i+1}^{n}-2u_{i}^{n}+u_{i-1}^{n}}{\Delta x^{2}}=0 

When coming to the code, :math:`u_{i}^{n+1}` is expressed in function of the other dependencies of :math:`u` :

.. math:: 
    
    u_{i}^{n+1} = u_{i}^{n} + \alpha (u_{i+1}^{n}-2u_{i}^{n}+u_{i-1}^{n}) \text{ with } \alpha= \frac{d \Delta t}{\Delta x^{2}}

.. code-block:: c++

    double alpha = (dt * d / (dx * dx));
    for (double t = 0; t <= tend; t += dt)
    {
        for (int i = 1; i < Nx - 1; i++)
        {
            U(i) += alpha * (U(i - 1) - 2 * U(i) + U(i + 1));
        }
    }

Here you have all the code at once :

.. code-block:: c++

    #include "castor/matrix.hpp"
    #include "castor/graphics.hpp"

    using namespace castor;

    int main(int argc, char *argv[])
    {
        //Thermal diffusivity
        double d = 1.;

        //Parameters
        int Nx = 1000;
        double xmin = 0;
        double xmax = 2;
        double tend = 0.1;
        double ti = 0;

        //Discretization
        double dx = (xmax - xmin) / (Nx - 1);
        double dt = 0.5 * (dx * dx) / d;
        double alpha = (dt * d / (dx * dx));
        matrix<> X = linspace(xmin, xmax, Nx);
        matrix<> U0 = pow(sin(X * M_PI), 2);

        std::cout << "--- Start explicit Euler scheme ---" << endl;
        tic();
        auto U = U0;
        for (double t = 0; t <= tend; t += dt)
        {
            for (int i = 1; i < Nx - 1; i++)
            {
                U(i) += alpha * (U(i - 1) - 2 * U(i) + U(i + 1));
            }
        }
        toc();

        //Plot
        figure fig;
        plot(fig, X, U0, {"b-o"}, {"initial"});
        plot(fig, X, S, {"g-s"}, {"solution"});
        plot(fig, X, U, {"m-x"}, {"explicit"});
        drawnow(fig);

        return 0;
    }

With this code you should get these outputs :

.. code-block:: text

    --- Start explicit Euler scheme ---
    Elapsed time is 0.213486 seconds.

.. image:: img/heatexplicit.png
    :width: 400
    :align: center


Implicit Euler
--------------

The specifity of the implicit Euler scheme is that the time derivative is calculated with the backward differential quotient approximation

.. math::
    
    \frac{\partial u_{i}^{n}}{\partial t} \approx \frac{u_{i}^{n}-u_{i}^{n-1}}{\Delta t}

| This scheme is stable for any ``dt`` .
| The scheme can be written using vectors

.. math:: 

    \frac{U^{n+1}-U^{n}}{\Delta t} + \frac{d}{\Delta x}AU^{n+1} = 0 ,

where A is the :math:`N_{x}` x :math:`N_{x}` tridiagonal matrix 

.. math:: 

    A = \begin{pmatrix} -2 & 1 & 0 & \cdots  & 0 
    \\ 1 & -2 & 1 & \cdots  & 0 
    \\ \vdots & \ddots  & \ddots  & \ddots  & \vdots 
    \\ 0 & \cdots  & 1 & -2 & 1 
    \\ 0 & \cdots  & 0 & 1 & -2 \end{pmatrix} .
    

This equation leads to the following linear equation 

.. math::
    
    BU^{n+1} = U^{n} \text{ with } B = (I_{N_{x}} - \alpha A) .


.. code-block:: c++

    double alpha = (dt * d / (dx * dx));
    matrix<> e = ones(Nx, 1);
    smatrix<> A = spdiags(cat(2, cat(2, e, -2 * e), e), colon(-1, 1), Nx, Nx);
    auto B = speye(Nx) - alpha * A;
    for (double t = 0; t <= tend; t += dt)
    {
        U = linsolve(B, U);
    }

Here you have all the code at once :

.. code-block:: c++

    #include "castor/matrix.hpp"
    #include "castor/graphics.hpp"
    #include "castor/linalg.hpp"

    using namespace castor;

    int main(int argc, char *argv[])
    {
        //Thermal diffusivity
        double d = 1.;

        //Parameters
        int Nx = 1000;
        double xmin = 0;
        double xmax = 2;
        double tend = 0.1;
        double ti = 0;

        //Discretization
        double dx = (xmax - xmin) / (Nx - 1);
        double dt = 5 * (dx * dx) / d;
        double alpha = (dt * d / (dx * dx));
        matrix<> X = linspace(xmin, xmax, Nx);
        matrix<> U0 = pow(sin(X * M_PI), 2);

        std::cout << "--- Start implicit Euler scheme ---" << endl;
        auto U = transpose(U0);
        tic();
        matrix<> e = ones(Nx, 1);
        smatrix<> A = spdiags(cat(2, cat(2, e, -2 * e), e), colon(-1, 1), Nx, Nx);
        auto B = speye(Nx) - alpha * A;
        for (double t = 0; t <= tend; t += dt)
        {
            U = linsolve(B, U);
        }
        toc();

        U = transpose(U);

        //Plot
        figure fig;
        plot(fig, X, U0, {"b-o"}, {"initial"});
        plot(fig, X, S, {"g-s"}, {"solution"});
        plot(fig, X, U, {"r-+"}, {"implicit"});
        drawnow(fig);

        return 0;
    }

With this code you should get these results :

.. code-block:: text

    --- Start implicit Euler scheme ---
    Elapsed time is 4.11192 seconds.

.. image:: img/heatimplicit.png
    :width: 400
    :align: center

References
----------

https://www.ljll.math.upmc.fr/ledret/M1English/M1ApproxPDE_Chapter6-2.pdf.. _label-castorfftw:

Castor FFTW
===========

**Castor FFTW** is the **castor** wrapper for the well-known FFTW3 library. The *git* repository can be found `here <https://github.com/marcbakry/castor-fftw>`_ and the documentation at `https://marcbakry.github.io/castor-fftw/ <https://marcbakry.github.io/castor-fftw/>`_.
.. _label-class-matrix:

Class matrix
+++++++++++++

The **castor** framework implements its own templatized class ``matrix<T>`` where ``T`` can be *for example* a ``float``, ``int``, ``std::complex<double>``, etc. At its core, it is built over a ``std::vector<T>`` which holds the values. The class ``matrix`` itself provides many useful functions and operators (addition, multiplication, indexing, etc.). It is designed such that it should feel like using `Matlab <https://www.mathworks.com>`_ or `Numpy arrays <htttps://www.numpy.org>`_. The user will find here all the available constructors. Specific builders (:ref:`label-ones`, :ref:`label-eye`, etc.) may be found at :ref:`label-builders`.

.. doxygenclass:: castor::matrix
   :project: castor
   :members: 

See :ref:`label-view`, :ref:`label-cview`.


.. _label-class-smatrix:

Class smatrix
+++++++++++++

.. doxygenclass:: castor::smatrix
   :project: castor
   :members: 
.. _label-fembem:

FEM BEM
=======

The project can be found here : `<https://gitlab.labos.polytechnique.fr/leprojetcastor/fembem>`_.
.. _label-installation:

Installation
============

The simplest way to get the dense matrix part of the framework **castor** is to download the `latest version of header files <https://gitlab.labos.polytechnique.fr/leprojetcastor/castor/-/jobs/artifacts/master/download?job=deploy>`_ and include `matrix.hpp` during the compilation of your c++ program using the library, see :ref:`label-compilation` in the :ref:`label-basic` section.

For a complete installation integrating check dependencies and examples compilation, we describe below the procedure.

**MacOS** (11 and 12) and **Ubuntu** (20.04)
++++++++++++++++++++++++++++++++++++++++++++

Installing the dependencies
---------------------------

The **castor** framework depends on two external dependencies : a BLAS/LAPACK implementation in order to use optimized linear algebra, and VTK for the graphical rendering.

The BLAS/LAPACK implementation which has been tested are :

- MKL 2020.0.166 : `MKL informations <https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html>`_,   
- OpenBlas 0.3.19 : `OpenBLAS informations <https://www.openblas.net/>`_,   
- vecLib : from accelerate framework in MacOS.   

The version of VTK library which has been tested is `9.1.0 <https://vtk.org/download/>`_.

The last tool to perform installation is ``CMake``, at least version `3.18  <https://cmake.org/download/>`_.

**Note** : on macOS it is recommended to install these dependencies with `brew <https://brew.sh/>`_.

From git repository with CMake
------------------------------

You can install the **castor** library from source with ``CMake``. The source files of the library is available here `<https://gitlab.labos.polytechnique.fr/leprojetcastor/castor>`_.

.. code::

    $ git clone https://gitlab.labos.polytechnique.fr/leprojetcastor/castor.git
    $ cd castor
    $ mkdir build
    $ cd build
    $ cmake -DCMAKE_INSTALL_PREFIX=path/to/install/directory ..
    $ make install

``path/to/install/directory`` is the absolute path to the folder where **castor** header files are installed. ``CMake`` assumes this folder contains ``include`` subfolder. By default, ``CMAKE_INSTALL_PREFIX`` is set to ``/usr/local`` and header files are installed in ``/usr/local/include/castor`` directory. 

The linear algebra part and the visualization part of the library depend respectively on an optimized BLAS library (like OpenBLAS and MKL) and the VTK library. ``CMake`` tries to detect these required libraries in your system, if ``Cmake`` can not find these dependencies, you can indicate the paths to them by using ``-DCMAKE_PREFIX_PATH``:

.. code::

    $ cmake -DCMAKE_INSTALL_PREFIX=path/to/install/directory -DCMAKE_PREFIX_PATH="/path/to/optimized/blas;/path/to/vtk/" ..   

Xcode project for macOS
-----------------------

To create a Xcode project with ``CMake``, add option `-G Xcode` in the `cmake` command :
 
.. code::

    $ cmake -G Xcode -DCMAKE_INSTALL_PREFIX=path/to/install/directory 

This project is created in the `build` directory.

Windows (10)
++++++++++++

This solution has been tested on Windows 10 only (but may work on other version). Since there is no *built_in* available package manager, the different components will be installed *by-hands* using only *Windows-like* tools.

C++ compiler and CMake
----------------------

The first step is to install a suitable C++ compiler. In the following instructions we will only use the compiler provided with `Visual Studio <https://visualstudio.microsoft.com/fr/downloads/>`_, version 2019 or later (previous version may work but have not been tested). The Visual Studio framework also provides a customized command prompt named `x64 Native Tools Command Prompt`.

We will also install the `CMake <https://cmake.org>`_ tools. Download the latest binary distribution for Windows. After installing CMake, open the `x64 Native Tools Command Prompt` and execute the command `cmake-gui`. If the command fails, find the install folder of CMake and add the `CMakeInstallFolder\bin` subfolder to the Windows `%PATH%` environment variable, restart the command prompt and try again. The graphical interface of CMake should open (close it for now).

BLAS/LAPACK library
-------------------

The simplest way is probably to install `OpenBLAS <https://www.openblas.net>`_ which implements both interfaces. Compiling the library can quickly become painful as a Fortran compiler is required. Thankfully, the developers have made precompiled binaries available. Installing OpenBLAS can be done following these steps:

1. Go to `https://github.com/xianyi/OpenBLAS/releases <https://github.com/xianyi/OpenBLAS/releases>`_, look for an archive named **OpenBLAS-0.x.x-x64.zip** (or **OpenBLAS-0.x.x-x86.zip** for older architectures) in the section **Assets** corresponding to the version you wish to use and download it. The demos were tested originally with `OpenBLAS 0.3.12` so any later version should be fine.
2. Extract the downloaded archive in a folder of your choice (for example, create a folder `openblas`). This folder should now contain three subdirectories:
    - `openblas\bin` should contain a file named `libopenblas.dll`.
    - `openblas\include` should contain the header files, including `cblas.h`.
    - `openblas\lib` should contain `.lib` files (including `libopenblas.lib`) and a `openblas\lib\cmake` subfolder.
3. Add the subfolder `openblas\bin` to Windows environment variable `%PATH%`.

The BLAS and LAPACK are now ready to use.

**Remark** It is also possible to download the Intel MKL library through the framework https://github.com/oneapi-src/oneMKL or the `Intel website <https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html>`_. However, this implementation features different header names and requires a modification of the source files of **castor** (namely replace `cblas.h` by `mkl_cblas.h` wherever it appears). For this reason, we do not insist further.

VTK framework
-------------

Unfortunately, the developers of the VTK framework do not provide *ready-to-use* binaries meaning that we must compile the sources by ourselves. It is performed as follows:

1. Download the sources of VTK on the main website `https://vtk.org <https://vtk.org>`_. Choose a version of the `9.x.x` branch. Uncompress the archive in a folder of your choice.
2. Open the `x64 Native Tools Command Prompt` and move to the newly created VTK folder (use the `dir pathToFolder` command). Create a *build* folder using `mkdir build` and move to this folder.
3. Execute `cmake-gui ..` which should open the CMake graphical interface. Click on `Configure`, choose the `ninja` generator and keep the default configuration. Finally, click on `Generate`. CMake will generate the build files.
4. Go back to the command prompt and execute the command `ninja`. The compilation of VTK begins and *may* take some time (a few minutes to a few dozen of minutes depending on the computer).
5. Once the compilation is over, execute `ninja install` which will install the library in the default directory `C:\Program Files (x86)\VTK`.
6. The last step is to add the subfolder `C:\Program Files (x86)\VTK\bin` to the Windows `%PATH%` environment variable.

The installation of VTK is now completed.

Compile the demos
-----------------

In this section, we will give the instructions on how to compile the examples of castor. The steps are the following:

1. Download the sources of **castor** from the `main repo <https://gitlab.labos.polytechnique.fr/leprojetcastor/castor.git>`_.
2. Open the `x64 Native Tools Command Prompt` and got to the `castor` folder. Create a `castor\build` directory and move to it.
3. Execute `cmake-gui ..` and click on the `Configure` button. Choose the `ninja` generator on the list and let all other options by default. This last operation *should fail* as CMake cannot find BLAS/LAPACK nor VTK.
4. In the list of CMake variables, look for `VTK_DIR` and set it to the `VTK\lib\cmake\vtk-9.x` folder. Do the same for the BLAS-related variables. Look for the variable `CBLAS_INCLUDE_DIR` and set it to the `openblas\include`subfolder.
5. Click again on `Configure` then on `Generate`.
6. Finally, execute `ninja` in the command prompt to start building the demo executable. The corresponding file can then be found in the `castor\build\demo` subfolder.
N-body problem
==============

*Shared by Antoine Rideau*

| On this page you will find how to simulate using **Castor** the n-body problem with 3 celestial bodies : the Sun, Jupiter and Saturn.
|
| In physics, the n-body problem is the problem of predicting the individual motions of a group of celestial objects interacting with each other gravitationally.
|
| First of all,  Newton's law of gravity says that the gravitational force felt on a planet *P* by the Sun *S* is given by

.. math::

    \overrightarrow{F}_{S\rightarrow P} = - \overrightarrow{F}_{P\rightarrow S} = -\frac{Gm_{S}m_{P}}{d^2}\overrightarrow{u} ,

| where :
|    :math:`G` : gravitational constant,
|    :math:`m_{S}`, :math:`m_{P}` : masses of the Sun and the planet,
|    :math:`d` : distance between the Sun and the planet,
|    :math:`\overrightarrow{u}` : monic vector directed from the Sun to the planet.
|
| The system is composed of the three heavier bodies in the solar system : the Sun, Jupiter and Saturn which will be referred as *S*, *Ju* and *Sa*. This system is supposed isolated and only gravitational forces are applied on the three celestial bodies.
| According to Newton's second principle, with :math:`m` the masses, :math:`\overrightarrow{a}` the acceleration and :math:`\overrightarrow{F}` the forces :

.. math::

    \begin{matrix}
    m_{Ju}\overrightarrow{a}_{Ju} = \overrightarrow{F}_{S\rightarrow Ju}  + \overrightarrow{F}_{Sa\rightarrow Ju} ,
    \\ 
    m_{Sa}\overrightarrow{a}_{Sa} = \overrightarrow{F}_{S\rightarrow Sa}  + \overrightarrow{F}_{Ju\rightarrow Sa} ,
    \\ 
    m_{S}\overrightarrow{a}_{S} = \overrightarrow{F}_{Ju\rightarrow S}  + \overrightarrow{F}_{Sa\rightarrow S} .
    \end{matrix}


+------------+----------------------------------------------------------------+
|   Bodies   |  Masses                                                        |
|            |  (relatively to the Sun)                                       |
+======================+======================================================+
| Sun and inner planet | :math:`m_{0}` = 1.00000597682                        |
+----------------------+------------+-----------------------------------------+
|        Jupiter       | :math:`m_{1}` = 0.000954786104043                    |
+----------------------+------------+-----------------------------------------+
|        Saturn        | :math:`m_{2}` = 0.000285583733151                    |
+----------------------+------------+-----------------------------------------+
| Gravitational constant :math:`G = 2.95912208286 \times 10^{-4}`             |
+-----------------------------------------------------------------------------+


.. code-block:: c++

    // Parameters
    matrix<> m = {{1.00000597682}, {9.54786104043e-4}, {2.85583733151e-4}}; // Masses : Sun, Jupiter and Saturn
    double G = 2.95912208286e-4;                                            //Gravitation's constant


The position of a object over time is given by :math:`q(t)`.


+------------+-----------------------------+
|   Bodies   |    Initial position (AU)    |
+============+=============================+
|            |               0             |
|            +-----------------------------+
|     Sun    |               0             |
|            +-----------------------------+
|            |               0             |
+------------+-----------------------------+
|            |          −3.5023653         |
|            +-----------------------------+
|  Jupiter   |          −3.8169847         |
|            +-----------------------------+
|            |          −1.5507963         |
+------------+-----------------------------+
|            |           9.0755314         |
|            +-----------------------------+
| Saturn     |          −3.0458353         |
|            +-----------------------------+
|            |          −1.6483708         |
+------------+-----------------------------+


where AU stands for astronomical unit and :math:`1 AU = 1.495 978 707 \times 10^{11}` m.

.. code-block:: c++
    
    matrix<> qini = {0, 0, 0,                               // Sun's initial position
                    -3.5023653, -3.8169847, -1.5507963,     // Jupiter's initial position
                    9.0755314, -3.0458353, -1.6483708};     // Saturn's initial position


| Moreover, using the momentum :math:`\overrightarrow{p} = m\overrightarrow{v} = m\overrightarrow{\dot{q}}(t)` instead of the speed is more practical. 
| Indeed, the system's total momentum is preserved over time :

.. math::

    \overrightarrow{p}_{S}(t) + \overrightarrow{p}_{Ju}(t) + \overrightarrow{p}_{Sa}(t) = Constant

+------------+-----------------------------+
|   Bodies   | Initial velocity (AU/day)   |
+============+=============================+
|            |               0             |
|            +-----------------------------+
|     Sun    |               0             |
|            +-----------------------------+
|            |               0             |
+------------+-----------------------------+
|            |           0.00565429        |
|            +-----------------------------+
|  Jupiter   |           0.00565429        |
|            +-----------------------------+
|            |          −0.00190589        |
+------------+-----------------------------+
|            |           0.00168318        |
|            +-----------------------------+
| Saturn     |           0.00483525        |
|            +-----------------------------+
|            |           0.00192462        |
+------------+-----------------------------+


.. code-block:: c++

    matrix<> vini = {0, 0, 0,                                                                   // Sun's initial velocity
    0.00565429, -0.00412490, -0.00190589,                                                       // Jupiter's initial velocity
    0.00168318, 0.00483525, 0.00192462};                                                        // Saturn's initial velocity
    matrix<> pini = reshape(mtimes(transpose(m), ones(1, numel(m))), 1, numel(vini)) * vini;    // Initial momentums

See :ref:`label-reshape` , :ref:`label-mtimes`, :ref:`label-transpose`, :ref:`label-ones`, :ref:`label-numel`.

Time is discretized into ``nt`` steps 

.. math::

    t_{i} = it \times \delta t \text{ for } it = \left [ \! \left [ 0, nt-1 \right ] \! \right ]

.. code-block:: c++

    // Disretization
    int nt = 1501;
    double dt = (tend - tini) / (nt - 1);
    auto T = linspace(tini, tend, nt);

See :ref:`label-linspace`.

Scheme
------

| A symplectic Euler scheme is used in this simulation because it preserved the energy of the system unlike either forward and backward Euler scheme.
| 
| As the system is conservative the Hamiltonian can be separated in a cinetical part :math:`K(p)` and a potential part :math:`V(q)` :

.. math::

    H(q,p) = K(p) + V(q) ,

where

.. math::

    \begin{matrix}
    \displaystyle K(p) = \frac{1}{2}\frac{p^2}{m}
    & \text{ and } &
    \displaystyle V(q_{i}) = \sum_{j\neq i}- \frac{Gm_{j}m_{i}}{\left | q_{i}-q_{j} \right |} .
    \end{matrix}

With such a separation, Hamilton equation are given by 

.. math::

    \begin{matrix}
    \displaystyle \frac{\mathrm{d} q}{\mathrm{d} t} = + \frac{\mathrm{d} K}{\mathrm{d} p}
    & \text{ and } &
    \displaystyle \frac{\mathrm{d} p}{\mathrm{d} t} = - \frac{\mathrm{d} V}{\mathrm{d} q} ,
    \end{matrix}

where

.. math::

    \begin{matrix}
    \displaystyle \frac{\mathrm{d} K(p)}{\mathrm{d} p} = \frac{p}{m}
    & \text{ and } &
    \displaystyle \frac{\mathrm{d} V(q_{i})}{\mathrm{d} q} = \sum_{j\neq i} \frac{Gm_{j}m_{i}\left ( q_{i}-q_{j} \right )}{\left | q_{i}-q_{j} \right |^3}
    \end{matrix}

which result to the symplectic Euler scheme :

.. math::

    \begin{matrix}
    \displaystyle q_{n+1} = q_{n} + \frac{\mathrm{d} K}{\mathrm{d} p}(p_{n})
    & \text{ and } &
    \displaystyle p_{n+1} = p_{n} - \frac{\mathrm{d} V}{\mathrm{d} q}(q_{n+1})
    \end{matrix}

.. code-block:: c++

    // Scheme
    auto Q = zeros(nt, numel(qini));
    Q(0, col(Q)) = qini;
    auto P = zeros(nt, numel(pini));
    P(0, col(P)) = pini;
    // Symplectic Euler
    for (int it = 0; it < nt - 1; it++)
    {
        matrix<> q_n = eval(Q(it, col(Q)));
        matrix<> p_n = eval(P(it, col(P)));
        Q(it + 1, col(Q)) = q_n + dt * H_p(G, m, p_n);
        P(it + 1, col(P)) = p_n - dt * H_q(G, m, eval(Q(it + 1, col(Q))));
    }

See :ref:`label-zeros`, :ref:`label-numel`, :ref:`label-col` , :ref:`label-view`.

In the code, :math:`\displaystyle \frac{\mathrm{d} K}{\mathrm{d} p}(p)` is represented by the function ``H_p`` 

.. code-block:: c++

    matrix<> H_p(double G, matrix<> m, matrix<> p)
    {
        auto Hp = zeros(1, numel(p));
        m = reshape(mtimes(transpose(m), ones(1, numel(m))), 1, numel(p));
        Hp = p / m;
        return Hp;
    }

See :ref:`label-zeros`, :ref:`label-reshape` , :ref:`label-mtimes`, :ref:`label-transpose`, :ref:`label-ones`, :ref:`label-numel`.

and :math:`\displaystyle \frac{\mathrm{d} V}{\mathrm{d} q}(q)` by the function ``H_q``

.. code-block:: c++

    matrix<> H_q(double G, matrix<> m, matrix<> q)
    {
    
        auto q0 = eval(q(range(0, 3)));
        auto q1 = eval(q(range(3, 6)));
        auto q2 = eval(q(range(6, 9)));
        auto Hq = zeros(1, 9);
        Hq(range(0, 3)) = (G * m(0) * m(1) * ((q0 - q1) / pow(norm(q0 - q1), 3)) + G * m(0) * m(2) * ((q0 - q2) / pow(norm(q0 - q2), 3)));
        Hq(range(3, 6)) = (G * m(1) * m(0) * ((q1 - q0) / pow(norm(q1 - q0), 3)) + G * m(1) * m(2) * ((q1 - q2) / pow(norm(q1 - q2), 3)));
        Hq(range(6, 9)) = (G * m(2) * m(0) * ((q2 - q0) / pow(norm(q2 - q0), 3)) + G * m(2) * m(1) * ((q2 - q1) / pow(norm(q2 - q1), 3)));
        return Hq;
    }

See :ref:`label-range`, :ref:`label-view`, :ref:`label-zeros`, :ref:`label-norm`.


Visualisation
--------------

Simple figure with Castor
^^^^^^^^^^^^^^^^^^^^^^^^^

| The simplest method to visualize the results is to plot them using ``plot`` or ``plot3``, here ``plot3`` is used to show the motion in 3 dimensions.
| 
| For each coordinates x,y and z, the Sun's positions are subtracted in order to keep it still in the center.
| Moreover, ``transpose`` is needed because of matrix ``Q`` 's dimensions.

.. code-block:: c++

    // Visu
    figure fig;
    plot3(fig, transpose(eval(Q(row(Q), 3)) - eval(Q(row(Q), 0))), transpose(eval(Q(row(Q), 4)) - eval(Q(row(Q), 1))), transpose(eval(Q(row(Q), 5)) - eval(Q(row(Q), 2))), {"c"});
    plot3(fig, transpose(eval(Q(row(Q), 6)) - eval(Q(row(Q), 0))), transpose(eval(Q(row(Q), 7)) - eval(Q(row(Q), 1))), transpose(eval(Q(row(Q), 8)) - eval(Q(row(Q), 2))), {"b"});
    plot3(fig, zeros(1, nt), zeros(1, nt), zeros(1, nt), {"y"});

.. figure:: img/3body.png
    :width: 800
    :align: center
    :figclass: align-center
    
    Orbits of Jupiter (cyan) and Saturn (blue) around the Sun (yellow) in the center.


Video output with VTK
^^^^^^^^^^^^^^^^^^^^^

| A way to visualize the results through a video is by using C++ VTK video writer.
| However, C++ VTK video writer's behavior is very OS dependent, it works perfectly fine on MAC but can cause some issues on Linux. (Windows ?)
| 
| First of all, the source and the writer need to be initialized with name of the output file, quality, framerate and connnected together.

.. code-block:: c++

    // Initialize source and movie
    vtkNew<vtkWindowToImageFilter> source;
    vtkNew<vtkOggTheoraWriter> movie;
    movie->SetInputConnection(source->GetOutputPort());
    movie->SetFileName("nbody.avi");
    movie->SetQuality(2); // in [0,2]
    movie->SetRate(25);   // frame per seconds
    int Nplot = 150;      // < 200

Afterwards, the writer is initiated before the loop over time. 

.. code-block:: c++

    movie->Start();
    for (int it = 0; it < nt - 1; it++){...}

Then each wanted frame is plotted and added to the final movie.

.. code-block:: c++

    // Visu
    if (it % (nt / Nplot) == 0)
    {
        figure fig;
        matrix<> limits = {-10, 10};
        plot(fig, {Q(it, 3) - Q(it, 0)}, {Q(it, 4) - Q(it, 1)}, {"c"});
        plot(fig, {Q(it, 6) - Q(it, 0)}, {Q(it, 7) - Q(it, 1)}, {"b"});
        plot(fig, zeros(1), zeros(1), {"y"});
        xlim(fig, limits);
        ylim(fig, limits);
        source->SetInput(fig.GetView()->GetRenderWindow());
        source->SetInputBufferTypeToRGB();
        source->ReadFrontBufferOff();
        movie->Write();
    }

See :ref:`label-plot`, :ref:`label-xlim`, :ref:`label-ylim`.

Finally, the writer is closed after the loop over time.

.. code-block:: c++

    for (int it = 0; it < nt - 1; it++){...}
    movie->End();



.. raw:: html

    <video controls width="100%">

    <source src="./_static/3bodyvtk.mp4"
            type="video/mp4">

    Sorry, your browser doesn't support embedded videos.
    </video>

|                   Orbit of Jupiter (orange) and Saturn (green) around the Sun in the center.



Video animation with Python
^^^^^^^^^^^^^^^^^^^^^^^^^^^

| Another way to get an animation of the orbiting planets is to output our data, and then post processing those with Python. 
| To do so, first the positions in the matrix ``Q`` are stored in a .txt file using ``writetxt`` .
| So as to keep the Sun still in the center, its positions are subtracted of Jupiter's and Saturn's positions. 

.. code-block:: c++

    // Output
    writetxt("./", "dataJu.txt", cat(2, eval(Q(row(Q), 3)) - eval(Q(row(Q), 0)), eval(Q(row(Q), 4)) - eval(Q(row(Q), 1))));
    writetxt("./", "dataSa.txt", cat(2, eval(Q(row(Q), 6)) - eval(Q(row(Q), 0)), eval(Q(row(Q), 7)) - eval(Q(row(Q), 1))));

See :ref:`label-writetxt`, :ref:`label-row` .

Then the following Python code shows the beautiful animation.

.. code-block:: python

    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    import numpy as np
    from collections import deque

    # Data input
    dataJu = np.loadtxt("./build/dataJu.txt")
    dataSa = np.loadtxt("./build/dataSa.txt")

    # Parameters extraction
    nt = int(dataJu[0, 0])

    # Data processing
    dataJu = np.delete(dataJu, 0, 0)
    dataSa = np.delete(dataSa, 0, 0)

    # Visu initialization
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(autoscale_on=False, xlim=(-10, 10), ylim=(-10, 10))
    ax.set_aspect('equal')

    line, = ax.plot([], [], 'o', lw=2)
    traceJu, = ax.plot([], [], ',-', lw=1)
    traceSa, = ax.plot([], [], ',-', lw=1)
    historyJu_x, historyJu_y = deque(maxlen=nt), deque(maxlen=nt)
    historySa_x, historySa_y = deque(maxlen=nt), deque(maxlen=nt)


    def animate(i):
        # Get planets' current positions
        thisx = [0, dataJu[i, 0], dataSa[i, 0]]
        thisy = [0, dataJu[i, 1], dataSa[i, 1]]

        # Clear the trace when the animation loops
        if i == 0:
            historyJu_x.clear()
            historyJu_y.clear()
            historySa_x.clear()
            historySa_y.clear()

        # Add the current position to the trace
        historyJu_x.appendleft(thisx[1])
        historyJu_y.appendleft(thisy[1])
        historySa_x.appendleft(thisx[2])
        historySa_y.appendleft(thisy[2])

        line.set_data(thisx, thisy)  # Update planets' positions
        # Update planets' traces
        traceJu.set_data(historyJu_x, historyJu_y)
        traceSa.set_data(historySa_x, historySa_y)
        return line, traceJu, traceSa


    # Creating the Animation object
    ani = animation.FuncAnimation(
        fig, animate, nt, interval=10, blit=True)
    plt.show()



.. raw:: html

    <video controls width="100%">

    <source src="./_static/3body.mp4"
            type="video/mp4">

    Sorry, your browser doesn't support embedded videos.
    </video>

|                           Orbit of Jupiter (orange) and Saturn (green) around the Sun in the center.



Code
----

Here is all the code at once, without the functions ``H_q`` and ``H_p``  written above :

.. code-block:: c++

    #include "castor/matrix.hpp"
    #include "castor/graphics.hpp"
    #include "castor/linalg.hpp"

    using namespace castor;

        int main(int argc, char const *argv[])
    {
        // Parameters
        matrix<> m = {1.00000597682, 9.54786104043e-4, 2.85583733151e-4}; // Masses : Sun, Jupiter and Saturn
        double G = 2.95912208286e-4;                                      //Gravitation's constant

        matrix<> qini = {0, 0, 0,                                                                // Sun's initial position
                         -3.5023653, -3.8169847, -1.5507963,                                     // Jupiter's initial position
                         9.0755314, -3.0458353, -1.6483708};                                     // Saturn's initial position
        matrix<> vini = {0, 0, 0,                                                                // Sun's initial velocity
                         0.00565429, -0.00412490, -0.00190589,                                   // Jupiter's initial velocity
                         0.00168318, 0.00483525, 0.00192462};                                    // Saturn's initial velocity
        matrix<> pini = reshape(mtimes(transpose(m), ones(1, numel(m))), 1, numel(vini)) * vini; // Initial momentums

        double tini = 0.;
        double tend = 12500.;

        // Disretization
        int nt = 1501;
        double dt = (tend - tini) / (nt - 1);
        auto T = linspace(tini, tend, nt);

        // Initialize source and movie
        vtkNew<vtkWindowToImageFilter> source;
        vtkNew<vtkOggTheoraWriter> movie;
        movie->SetInputConnection(source->GetOutputPort());
        movie->SetFileName("nbody.avi");
        movie->SetQuality(2); // in [0,2]
        movie->SetRate(25);   // frame per seconds
        int Nplot = 150;      // < 200

        // Scheme
        auto Q = zeros(nt, numel(qini));
        Q(0, col(Q)) = qini;
        auto P = zeros(nt, numel(pini));
        P(0, col(P)) = pini;
        // Symplectic Euler
        tic();
        movie->Start();
        for (int it = 0; it < nt - 1; it++)
        {

            matrix<> q_n = eval(Q(it, col(Q)));
            matrix<> p_n = eval(P(it, col(P)));
            Q(it + 1, col(Q)) = q_n + dt * H_p(G, m, p_n);
            P(it + 1, col(P)) = p_n - dt * H_q(G, m, eval(Q(it + 1, col(Q))));

            // Visu
            if (it % (nt / Nplot) == 0)
            {
                figure fig;
                matrix<> L({-10, 10, -10, 10}); // Axis dimensions 
                plot(fig, Q(it, 3) - Q(it, 0) * ones(1), Q(it, 4) - Q(it, 1) * ones(1), L, {"c"});
                plot(fig, Q(it, 6) - Q(it, 0) * ones(1), Q(it, 7) - Q(it, 1) * ones(1), L, {"b"});
                plot(fig, zeros(1), zeros(1), {"y"});
                source->SetInput(fig.GetView()->GetRenderWindow());
                source->SetInputBufferTypeToRGB();
                source->ReadFrontBufferOff();
                movie->Write();
            }
        }
        movie->End();
        toc();

        // Output for Python post-processing and visualisation
        // writetxt("./", "dataJu.txt", cat(2, eval(Q(row(Q), 3)) - eval(Q(row(Q), 0)), eval(Q(row(Q), 4)) - eval(Q(row(Q), 1))));
        // writetxt("./", "dataSa.txt", cat(2, eval(Q(row(Q), 6)) - eval(Q(row(Q), 0)), eval(Q(row(Q), 7)) - eval(Q(row(Q), 1))));

        // Visu in native Castor
        // figure fig;
        // plot3(fig, transpose(eval(Q(row(Q), 3)) - eval(Q(row(Q), 0))), transpose(eval(Q(row(Q), 4)) - eval(Q(row(Q), 1))), transpose(eval(Q(row(Q), 5)) - eval(Q(row(Q), 2))), {"c"});
        // plot3(fig, transpose(eval(Q(row(Q), 6)) - eval(Q(row(Q), 0))), transpose(eval(Q(row(Q), 7)) - eval(Q(row(Q), 1))), transpose(eval(Q(row(Q), 8)) - eval(Q(row(Q), 2))), {"b"});
        // plot3(fig, zeros(1, nt), zeros(1, nt), zeros(1, nt), {"y"});

        // drawnow(fig);

        return 0;
    }





References
----------

| https://interstices.info/les-planetes-tournent-elles-rond/
|
| http://www.unige.ch/~hairer/poly/chap3.pdf.. _label-operators:

Operators
+++++++++

This section describes all the conventional ``C++`` operators which have been overloaded in ``matrix``. To the exception of ``<<``, they all work element-wise.


.. _label-operator<<:

operator<<
----------
.. doxygenfunction:: operator<<(std::ostream &flux, matrix<T> const &A)
   :project: castor

See :ref:`label-disp`.


.. _label-operator!:

operator!
---------
.. doxygenfunction:: operator!
   :project: castor

See :ref:`label-operator&&`, :ref:`label-operator||`.


.. _label-operator&&:

operator&&
----------
.. doxygenfunction:: operator&&(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator&&(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator&&(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator||`, :ref:`label-operator!`.


.. _label-operator||:

operator||
----------
.. doxygenfunction:: operator||(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator||(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator||(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator&&`, :ref:`label-operator!`.


.. _label-operator==:

operator==
----------
.. doxygenfunction:: operator==(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator==(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator==(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator!=`, :ref:`label-operator<=`, :ref:`label-operator>=`.


.. _label-operator!=:

operator!=
----------
.. doxygenfunction:: operator!=(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator!=(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator!=(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator==`, :ref:`label-operator<=`, :ref:`label-operator>=`.


.. _label-operator<=:

operator<=
----------
.. doxygenfunction:: operator<=(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator<=(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator<=(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator==`, :ref:`label-operator>=`, :ref:`label-operator<`.


.. _label-operator<:

operator<
---------
.. doxygenfunction:: operator<(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator<(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator<(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator<=`, :ref:`label-operator>`, :ref:`label-operator>=`.


.. _label-operator>=:

operator>=
----------
.. doxygenfunction:: operator>=(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator>=(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator>=(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator>`, :ref:`label-operator<=`, :ref:`label-operator<`.


.. _label-operator>:

operator>
---------
.. doxygenfunction:: operator>(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator>(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator>(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator>=`, :ref:`label-operator<`, :ref:`label-operator<=`.


.. _label-operator+:

operator+
---------
.. doxygenfunction:: operator+(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator+(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator+(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator-`, :ref:`label-operator*`, :ref:`label-operator/`.


.. _label-operator-:

operator-
---------
.. doxygenfunction:: operator-(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator-(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator-(matrix<R> const &A, matrix<S> const &B)
   :project: castor

.. doxygenfunction:: operator-(matrix<S> const &A)
   :project: castor

See :ref:`label-operator+`, :ref:`label-operator*`, :ref:`label-operator/`.


.. _label-operator*:

operator*
---------
.. doxygenfunction:: operator*(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator*(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator*(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator+`, :ref:`label-operator-`, :ref:`label-operator/`.


.. _label-operator/:

operator/
---------
.. doxygenfunction:: operator/(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator/(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator/(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator+`, :ref:`label-operator-`, :ref:`label-operator*`.
.. _label-linear-solver-func:

Linear solver 
+++++++++++++

This section regroups all of the linear algebra functions concerning linear solver which are currently implemented within **castor**. In order to access these functions, **the user must include** ``castor/linalg.hpp``. Two levels of interface are available. The high-level interface provides functions with simplified arguments while the low-level interface is much closer to the BLAS/LAPACK API.

The user will find here high-level linear algebra functions. Some examples of use may also be found at :ref:`label-linear-algebra-advanced`.


.. _label-inv:

inv
---
.. doxygenfunction:: inv(matrix<T> const &A)
   :project: castor

See :ref:`label-linsolve`, :ref:`label-pinv`, :ref:`label-gmres`.


.. _label-linsolve:

linsolve
--------
.. doxygenfunction:: linsolve(matrix<std::complex<double>> const &A, matrix<std::complex<double>> const &B)
   :project: castor
.. doxygenfunction:: linsolve(matrix<double> const &A, matrix<double> const &B)
   :project: castor
.. doxygenfunction:: linsolve(matrix<std::complex<float>> const &A, matrix<std::complex<float>> const &B)
   :project: castor
.. doxygenfunction:: linsolve(matrix<float> const &A, matrix<float> const &B)
   :project: castor

See :ref:`label-inv`, :ref:`label-pinv`, :ref:`label-gmres`.


.. _label-pinv:

pinv
----
.. doxygenfunction:: pinv(matrix<T> const &A)
   :project: castor

See :ref:`label-inv`, :ref:`label-pinv`.
Welcome to the castor library documentation
===========================================

The objective of the **castor** library is to propose **high-level semantics**, inspired by the Matlab language, allowing fast software prototyping in a low-level compiled language. It is nothing more than a **matrix management layer** using the tools of the standard **C++** library, in different storage formats (full, sparse and hierarchical). Indeed, the use of IDEs 1 such as Xcode, Visual studio, Eclipse, etc. allows today to execute compiled code (C, C++, fortran, etc.) with the **same flexibility as interpreted languages** (Matlab, Python, Julia, etc.). 

A **header-only** template library for matrix management has been developed based on the standard C++ library, notably the std::vector class. Many tools and algorithms are provided to simplify the development of scientific computing programs. Particular attention has been paid to semantics, for a simplicity of use "à la matlab", but written in C++. This high-level semantic/low-level language coupling makes it possible to gain efficiency in the prototyping phase, while ensuring performance for applications. In addition, direct access to data allows users to optimize the most critical parts of their code in native C++. Finally, **complete documentation** is available, as well as continuous integration unit tests. All of this makes it possible to meet the needs of teaching, academic issues and industrial applications at the same time. 

The **castor** library provides tools to : 

- create and manipulate dense, sparse and hierarchical matrices
- make linear algebra computations based on optimized BLAS library
- make graphical representations based on VTK library

These tools are used by applicative projects : 

- finite and boundary element method using Galerkin approximation
- analytical solutions for scattering problems 

The source files of the library is available here : `<https://gitlab.labos.polytechnique.fr/leprojetcastor/castor>`_.

As the semantics offered by **castor** library being voluntarily close to the Matlab environment, there are functions signature and their documentation inspired by it. You can refer to https://fr.mathworks.com/help/matlab/index.html.

Licensing
---------

The **castor** library is provided in open source under LGPL 3.0.

This program is free software, distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,  redistribute and/or modify it under the terms of the GNU Lesser General Public License, as published by the Free Software Foundation (version 3 or later,  http://www.gnu.org/licenses).


Gallery
-------
.. image:: img/head.png
   :width: 600
   :align: center

Boundary element computation using le projet castor. Simulation of acoustic scattering by human head, excited by plane wave at 8Khz (20.000 degrees of freedom). Hierarchical solver is used on a standard laptop.

.. toctree::
   :caption: Installation
   :maxdepth: 1
   :hidden:

   installation

.. toctree::
   :caption: Examples
   :maxdepth: 1
   :hidden:

   heat
   eigenmode
   montecarlo
   replicatormutator
   nbody
   boids

.. _label-user-guide:

.. toctree::
   :caption: User guide
   :maxdepth: 1
   :hidden:

   getting_started
   basics
   advanced
   linalg
   graphics
   sparse_matrix 
   kissfft

.. toctree::
   :caption: Dense matrix
   :maxdepth: 1
   :hidden:

   class_matrix
   class_view_cview
   algorithm
   builders
   geometry
   io
   math
   dimensions
   manipulations
   operators
   tools 
   transforms 

.. toctree::
   :caption: Linear algebra
   :maxdepth: 1
   :hidden:

   factorization
   linear_solver
   singular_eig_values
   lowlevel_linalg_func

.. toctree::
   :caption: Graphical rendering
   :maxdepth: 1
   :hidden:

   class_figure
   basic_plot
   graphical_io
   graphical_tools
   mesh_management
   mesh_plot

.. toctree::
   :caption: Sparse matrix
   :maxdepth: 1
   :hidden:

   class_smatrix
   api_smatrix 

.. toctree::
   :caption: Hierarchical matrix
   :maxdepth: 1
   :hidden:

   class_hmatrix
   class_bintree
   api_hmatrix 

.. toctree::
   :caption: Applications
   :maxdepth: 1
   :hidden:

   fembem
   analyticalscattering
   castorfftw

.. toctree::
   :caption: Contacts
   :maxdepth: 1
   :hidden:
 
   developpers

.. _label-transforms:

Transforms
++++++++++

This section describes functions concerning Fourier transform. 

.. _label-dft:

dft
---
.. doxygenfunction:: dft
   :project: castor

See :ref:`label-idft`, :ref:`label-fft`.

.. _label-fft:

fft
---
.. doxygenfunction:: fft
   :project: castor

See :ref:`label-dft`, :ref:`label-ifft`.

.. _label-idft:

idft
----
.. doxygenfunction:: idft
   :project: castor

See :ref:`label-dft`, :ref:`label-ifft`.

.. _label-ifft:

ifft
----
.. doxygenfunction:: ifft
   :project: castor

See :ref:`label-idft`, :ref:`label-fft`... _label-input-output:

Input/Output
++++++++++++

The ``matrix`` object can be saved on the drive either as a text file or a binary file.


.. _label-readbin:

readbin
-------
.. doxygenfunction:: readbin
   :project: castor

See :ref:`label-writebin`, :ref:`label-readtxt`.


.. _label-readtxt:

readtxt
-------
.. doxygenfunction:: readtxt
   :project: castor

See :ref:`label-writetxt`, :ref:`label-readbin`.


.. _label-writebin:

writebin
--------
.. doxygenfunction:: writebin
   :project: castor

See :ref:`label-readbin`, :ref:`label-writetxt`.


.. _label-writetxt:

writetxt
--------
.. doxygenfunction:: writetxt
   :project: castor

See :ref:`label-readtxt`, :ref:`label-writebin`... _label-advanced:

Advanced
========

We present here some advanced features of the ``matrix`` class.

Working with complex numbers
++++++++++++++++++++++++++++

It is possible to manipulate complex data with the ``matrix`` class. The code below creates a complex-double ``matrix`` using the predefined constant ``M_1I`` which is equal to ``std::complex<double>(0,1)`` and two matrices containing respectively the real and the imaginary part.


.. code:: c++

    matrix<> Re = {1,0,-1,0};
    matrix<> Im = {0,1,0,-1};

    disp("Imaginary number matrix ('1i') : ");
    disp(M_1I);

    disp("Complex matrix : ");
    auto Ac = Re + M_1I*Im;
    disp(Ac);


.. code:: text

  Imaginary number matrix ('1i') :
  (0,1)
  Complex matrix :
  (1.00000,0.00000)          (0.00000,1.00000)         (-1.00000,0.00000)         (0.00000,-1.00000)

We can compute the conjugate of the elements, the absolute value, or extract the imaginary part.

.. code:: c++

  disp(conj(Ac),2);
  disp(abs(Ac),2);
  disp(imag(Ac),2);

.. code:: text

  Matrix 1x4 of type 'St7complexIdE' (64 B):
        (1.00000,-0.00000)         (0.00000,-1.00000)        (-1.00000,-0.00000)          (0.00000,1.00000)  
  Matrix 1x4 of type 'd' (32 B):
      1.00000      1.00000      1.00000      1.00000  
  Matrix 1x4 of type 'd' (32 B):
            0      1.00000            0     -1.00000


Logical matrices
++++++++++++++++

Because of the special behaviour of ``std::vector<bool>``,  **matrix** does not use the native logical type ``bool`` of the C++ standard library, but uses instead ``std::uint8_t`` which corresponds to a ``unsigned short int``, or ``char``. Consequently, the ``bool`` values are converted to ``std::uint8_t`` where the value ``false`` is converted to the value ``0`` and ``true`` to the value ``1``. A logical ``matrix`` is displayed like:

.. code:: c++

    matrix<int> A = eye(3);
    matrix<int> B = ones(3);
    disp("A && B :");
    disp(A && B);

.. code:: text

    A && B :
    1  0  0
    0  1  0
    0  0  1

**WARNING (very important):** a ``matrix<logical>`` *remains intrinsically a* ``matrix<std::uint8_t>`` meaning that it behaves like one. This behavior is not natural and will be corrected in a future release.


Using some algorithms
+++++++++++++++++++++

The ``matrix`` class comes with a lot of built-in algorithms (see :ref:`label-algorithms` for a full list). First, we create some random integer data in the range ``[0,10[``.

.. code:: c++

  auto A = cast<int>(10*rand<>(1,10));
  disp(A,2,std::cout,10,10);

.. code:: text

  Matrix 1x10 of type 'i' (40 B):
  8  3  7  7  9  1  3  7  2  5

**Remark:** In the code above, we first call :ref:`label-rand` for the default ``double`` type, then we multiply by ``10`` and finally :ref:`label-cast` the result to integers. In fact, the :ref:`label-rand` function generates intrinsically numbers of type ``double`` in the range ``[0,1[`` which are converted afterward in the template type. Unfortunately, calling ``rand<int>`` will cast those numbers to zero. What we did is first generate ``double`` numbers in the range ``[0,10[``, then cast them to ``int``.

What are the minimum, average, standard deviation, and the maximum values of ``A`` ?

.. code:: c++

  disp(min(A));
  disp(mean<double>(A));
  disp(stddev<double>(A));
  disp(max(A));


.. code:: text

  1
  5.2
  2.63818
  9

**Remark:** Do not forget, in that particular case, to give the output type for the result of ``mean`` and ``stddev``. Otherwise, the result will be cast in the template type which is ``int``!

The :ref:`label-unique` elements of ``A`` are

.. code:: c++

  disp(unique(A),2,std::cout,10,10);

.. code:: text

  Matrix 1x7 of type 'i' (28 B):
  1  2  3  5  7  8  9

Now let us :ref:`label-sort` the elements.

.. code:: c++

  auto A_sorted = sort(A);
  disp(A_sorted,2,std::cout,10,10);

.. code:: text

  Matrix 1x10 of type 'i' (40 B):
  1  2  3  3  5  7  7  7  8  9  

We create a second smaller vector and we search the common elements.

.. code:: c++

  auto B = cast<int>(10*rand<>(1,4));
  disp(B,2);
  disp(intersect(A,B),2);

.. code:: text

  Matrix 1x4 of type 'i' (16 B):
  4  6  3  5  
  Matrix 1x2 of type 'i' (8 B):
  3  5  

We obtain the expected result.


View
++++

The **matrix** provides operators to extract a submatrix from a matrix instance or to assign values of a submatrix to a matrix. 
The operator ``()`` can take a list of indices to return an instance of ``class view``. 

As the ``class view`` contains a reference to the set of values corresponding to the list of indices, it is **necessary** to call the method ``eval`` to use the viewed submatrix.

It is possible to give the list of indices in linear indexing **A(L)** or bilinear indexing **A(I,J)**. The external functions all, col, range and row are useful to describe indices.

.. code:: c++

    matrix<> A = {{0, 1, 2, 4},
                  {4, 5, 6, 7},
                  {8, 9,10,11}};

    disp("Linear indexing");
    disp("  extracting :");
    matrix<> B = eval(A({1,3,5}));
    disp(B);

    disp("  assigning :"); 
    A(range(0,4)) = 0;
    disp(A);

    disp("Bilinear indexing");
    disp("  extracting :");
    B = eval(A(1,range(0,3)));;
    disp(B);

    disp("  assigning :"); 
    A({0,2}, col(A)) = 10;
    disp(A);

.. code:: text

  Linear indexing
    extracting :
      1.00000      4.00000      5.00000
    assigning :
            0            0            0            0
      4.00000      5.00000      6.00000      7.00000
      8.00000      9.00000     10.00000     11.00000
  Bilinear indexing
    extracting :
      4.00000      5.00000      6.00000
    assigning :
     10.00000     10.00000     10.00000     10.00000
      4.00000      5.00000      6.00000      7.00000
     10.00000     10.00000     10.00000     10.00000   


Internal matrix tools
+++++++++++++++++++++

In addition to constructors and operators, the matrix class has few member functions (see internal tools list in :ref:`class matrix description <label-class-matrix>`). The functions can be user-friendy for native C++ coding.

Notes:

- The method ``size``  differs from the external function ``size`` for the parameter dim=0. Indeed, in that case, the method ``size`` returns the total number of elements in the matrix while the external function ``size`` returns the two-element vector containing the number of rows and columns in the matrix.

.. code:: c++

    matrix<> A = {{1,2,3},{4,5,6}};
    disp("Method size :");
    disp(A.size(0));
    disp(A.size(1));
    disp(A.size(2));
    disp("External function size :");
    disp(size(A));
    disp(size(A,1));
    disp(size(A,2));

.. code:: text

   Method size :
   6
   2
   3
   External function size :
   2  3
   2
   3

- The method ``resize`` is more efficient than the external function ``resize`` since the methode makes the maniplulation in-place which avoids a copy.  

Low-level interface
+++++++++++++++++++

The user will find here the low-level linear algebra interface. Most of the functions here are wrappers around the CBLAS or LAPACK API and the user should be careful when using them, see :ref:`label-blaslapack-issue` to understand why.

.. _label-lpk2mat:

lpk2mat
-------
.. doxygenfunction:: lpk2mat(std::vector<S>const &V, std::size_t m, std::size_t n)
   :project: castor

See :ref:`label-mat2lpk`.


.. _label-mat2lpk:

mat2lpk
-------
.. doxygenfunction:: mat2lpk(matrix<S> const &A, std::size_t L)
   :project: castor

See :ref:`label-lpk2mat`.


.. _label-tgeev:

tgeev
-----
.. doxygenfunction:: tgeev(std::string typ, int &n, std::vector<T> &A, std::vector<T> &E, std::vector<T> &V)
   :project: castor
.. doxygenfunction:: xgeev(char &jobl, char &jobr, int &n, std::vector<zlpk> &A, std::vector<zlpk> &E, std::vector<zlpk> &V, zlpk &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgeev(char &jobl, char &jobr, int &n, std::vector<clpk> &A, std::vector<clpk> &E, std::vector<clpk> &V, clpk &wkopt, int &lwork, int &info)
   :project: castor


.. _label_tgels:

tgels
-----
.. doxygenfunction:: tgels(int &m, int &n, int &nrhs, std::vector<T> &A, std::vector<T> &B)
   :project: castor
.. doxygenfunction:: xgels(char &t, int &m, int &n, int &nrhs, std::vector<zlpk> &A, std::vector<zlpk> &B, int &l, zlpk &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgels(char &t, int &m, int &n, int &nrhs, std::vector<double> &A, std::vector<double> &B, int &l, double &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgels(char &t, int &m, int &n, int &nrhs, std::vector<clpk> &A, std::vector<clpk> &B, int &l, clpk &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgels(char &t, int &m, int &n, int &nrhs, std::vector<float> &A, std::vector<float> &B, int &l, float &wkopt, int &lwork, int &info)
   :project: castor


.. _label-tgemm-blas:

tgemm
-----
.. doxygenfunction:: tgemm(R alpha, matrix<double> const &A, matrix<double> const &B, S beta, matrix<double> &C)
   :project: castor
.. doxygenfunction:: tgemm(R alpha, matrix<std::complex<float>> const &A, matrix<std::complex<float>> const &B, S beta, matrix<std::complex<float>> &C)
   :project: castor
.. doxygenfunction:: tgemm(R alpha, matrix<std::complex<double>> const &A, matrix<std::complex<double>> const &B, S beta, matrix<std::complex<double>> &C)
   :project: castor
.. doxygenfunction:: tgemm(R alpha, matrix<float> const &A, matrix<float> const &B, S beta, matrix<float> &C)
   :project: castor

See :ref:`label-tgemm-naive`.


.. _label-tgeqrf:

tgeqrf
------
.. doxygenfunction:: tgeqrf(int &m, int &n, std::vector<T> &A, std::vector<T> &R)
   :project: castor
.. doxygenfunction:: xgeqrf(int &m, int &n, int &l, std::vector<zlpk> &A, std::vector<zlpk> &tau, zlpk &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgeqrf(int &m, int &n, int &l, std::vector<double> &A, std::vector<double> &tau, double &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgeqrf(int &m, int &n, int &l, std::vector<clpk> &A, std::vector<clpk> &tau, clpk &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgeqrf(int &m, int &n, int &l, std::vector<float> &A, std::vector<float> &tau, float &wkopt, int &lwork, int &info)
   :project: castor


.. _label-tgesdd:

tgesdd
------
.. doxygenfunction:: tgesdd(std::string typ, int &m, int &n, std::vector<T> &A, std::vector<T2> &S, std::vector<T> &U, std::vector<T> &V)
   :project: castor
.. doxygenfunction:: xgesdd(char &job, int &m, int &n, int &l, std::vector<zlpk> &A, std::vector<double> &S, std::vector<zlpk> &U, std::vector<zlpk> &V, zlpk &wkopt, int &lwork, std::vector<int> &iwork, int &info)
   :project: castor
.. doxygenfunction:: xgesdd(char &job, int &m, int &n, int &l, std::vector<double> &A, std::vector<double> &S, std::vector<double> &U, std::vector<double> &V, double &wkopt, int &lwork, std::vector<int> &iwork, int &info)
   :project: castor
.. doxygenfunction:: xgesdd(char &job, int &m, int &n, int &l, std::vector<clpk> &A, std::vector<float> &S, std::vector<clpk> &U, std::vector<clpk> &V, clpk &wkopt, int &lwork, std::vector<int> &iwork, int &info)
   :project: castor
.. doxygenfunction:: xgesdd(char &job, int &m, int &n, int &l, std::vector<float> &A, std::vector<float> &S, std::vector<float> &U, std::vector<float> &V, float &wkopt, int &lwork, std::vector<int> &iwork, int &info)
   :project: castor


.. _label-tgesv:

tgesv
-----
.. doxygenfunction:: tgesv(int &n, int &nrhs, std::vector<T> &A, std::vector<T> &B)
   :project: castor
.. doxygenfunction:: xgesv(int &n, int &nrhs, std::vector<zlpk> &A, std::vector<zlpk> &B, std::vector<int> &ipiv, int &info)
   :project: castor
.. doxygenfunction:: xgesv(int &n, int &nrhs, std::vector<double> &A, std::vector<double> &B, std::vector<int> &ipiv, int &info)
   :project: castor
.. doxygenfunction:: xgesv(int &n, int &nrhs, std::vector<clpk> &A, std::vector<clpk> &B, std::vector<int> &ipiv, int &info)
   :project: castor
.. doxygenfunction:: xgesv(int &n, int &nrhs, std::vector<float> &A, std::vector<float> &B, std::vector<int> &ipiv, int &info)
   :project: castor


.. _label-tgetrf:

tgetrf
------
.. doxygenfunction:: tgetrf(int &m, int &n, std::vector<T> &A, std::vector<T> &P, matrix<S> &U)
   :project: castor
.. doxygenfunction:: xgetrf(int &m, int &n, std::vector<zlpk> &A, std::vector<int> &ipiv, int &info)
   :project: castor
.. doxygenfunction:: xgetrf(int &m, int &n, std::vector<double> &A, std::vector<int> &ipiv, int &info)
   :project: castor
.. doxygenfunction:: xgetrf(int &m, int &n, std::vector<clpk> &A, std::vector<int> &ipiv, int &info)
   :project: castor
.. doxygenfunction:: xgetrf(int &m, int &n, std::vector<float> &A, std::vector<int> &ipiv, int &info)
   :project: castor
Replicator mutator equation
===========================

*Shared by Antoine Rideau thanks to Gael Raoul*

On this page you will find how to simulate using **Castor** the replicator mutator equation

The replicator-mutator equation is a classical model from evolutionary biology that describes how a species reacts to selection by increasing a phenotypic trait :math:`x`. The selection is represented by a reproduction rate that increases linearly with :math:`x` (the phenotype is actually the  fitness of the individual), while the population is kept constant thanks to the continuous removal of individuals, uniformly among the traits present in the population. This is completed by a mutation term: the traits constantly mutate which is described by a diffusion term. This model is used to understand mutation-selection dynamics as a whole, even though it is more directly related to experimental setups used in experimental evolutionary biology, based on chemostats. This model has been used by R. Fisher to derive the so-called *fundamental theorem of natural selection*.


.. math:: 

    \partial_{t}u = \underbrace{\sigma \Delta_{x}u}_{mutations} + \underbrace{(x - \bar{x}(t))u}_{replication} \text{ , } t > 0, 

| where :
|    :math:`x \in \mathbb{R}` : a one dimension fitness space,
|    :math:`u(t,x)` : density of a population at time :math:`t` and per unit of fitness,
|    :math:`\bar{x}(t):= \int_{\mathbb{R}}xu(x,t)dx` : mean fitness at time :math:`t` .


The population is considered constant, so

.. math::

    \int_{\mathbb{R}}u(x,t)dx = 1 .



Numeric simulation
------------------

|   ``N`` individuals are gathered within the population, each characterized by their fitness :math:`x_{i}`.
|   Ths population will evolve during ``gmax`` generations separated by ``dt`` .

.. code-block:: c++

    // Parameters
    int N = 1e3;        // Population
    int gmax = 1e4;     // Number of generations
    int Nplot = 1000;     // Number of generations plotted
    double dt = 0.01;   // Time disretization
    double sigma = 0.5; // Mutation


Initially, the population is distributed following a Gaussian distribution using ``randn`` .

.. code-block:: c++

    // Initial data
    auto parent = randn(1, N);

See :ref:`label-randn` . 

Each generation :

1. Each individual has a probability :math:`\mathbb{P} = (x_{i})_{+} \times \Delta t` ,where :math:`(x_{i})_{+}` stands for the positive part of :math:`x_{i}` , to give birth to a child

.. code-block:: c++

    // Probablity to give birth
    auto birth = dt * maximum(parent, 0);

    // Reproduction
    auto reprod = rand(size(parent));  

See :ref:`label-maximum` , :ref:`label-rand` , :ref:`label-size` .

who will inherit a fitness of :math:`x_{i} + X` with :math:`X \sim \mathcal{N}(0, \sigma^2 \Delta t)` .

.. code-block:: c++

    // Children
    auto children = parent + sigma * std::sqrt(dt) * randn(1, N);
    children = eval(children(find(reprod < birth)));

    // Update parent
    parent = cat(2, parent, children);

See :ref:`label-find` , :ref:`label-view` , :ref:`label-cat` . 

2. ``N`` individuals are uniformly choosen  to survive.

.. code-block:: c++

    // Kill parent to get N individuals
    parent = eval(parent(randperm(numel(parent), N)));

See :ref:`label-randperm` , :ref:`label-numel` .

Code
----

.. code-block:: c++

    #include <castor/matrix.hpp>
    #include <castor/graphics.hpp>

    using namespace castor;

    int main(int argc, char const *argv[])
    {
        // Parameters
        int N = 1e3;        // Population
        int gmax = 1e4;     // Number of generations
        int Nplot = 1000;     // Number of generations plotted
        double dt = 0.01;   // Time disretization
        double sigma = 0.5; // Mutation

        // Initial data
        auto parent = randn(1, N);

        // Initialize figure
        figure fig;

        // For each generation
        tic();
        for (int g = 1; g <= gmax; g++)
        {
            // Probablity to give birth
            auto birth = dt * maximum(parent, 0);

            // Reproduction
            auto reprod = rand(size(parent));

            // Children
            auto children = parent + sigma * std::sqrt(dt) * randn(1, N);
            children = eval(children(find(reprod < birth)));

            // Update parent
            parent = cat(2, parent, children);

            // Kill parent to get N individuals
            parent = eval(parent(randperm(numel(parent), N)));

            // Plot
            if (g % (gmax / Nplot) == 0)
            {
                plot(fig, parent, g * dt * ones(size(parent)), {"b"});
            }
        }
        toc();

        // Visu
        drawnow(fig);
        return 0;
    }


.. figure:: img/replicatormutator.png
    :width: 800
    :align: center
    :figclass: align-center
    
    Fitness evolution of a 1 000 individuals' population during 10 000 generations.

Reference
---------

https://openlibrary.org/books/OL7084333M/The_genetical_theory_of_natural_selection.

https://www.cirm-math.fr/RepRenc/1315/PDFfiles1315.pdf








.. _label-singular-eig-values-func:

Singular and eigen values
+++++++++++++++++++++++++

This section regroups all of the linear algebra functions concerning singular and eigen values which are currently implemented within **castor**. In order to access these functions, **the user must include** ``castor/linalg.hpp``. Two levels of interface are available. The high-level interface provides functions with simplified arguments while the low-level interface is much closer to the BLAS/LAPACK API.

The user will find here high-level linear algebra functions. Some examples of use may also be found at :ref:`label-linear-algebra-advanced`.

.. _label-eig:

eig
---
.. doxygenfunction:: eig(matrix<T> const &A, matrix<T> const &B)
   :project: castor
.. doxygenfunction:: eig(matrix<std::complex<double>> const &A, matrix<std::complex<double>> const &B, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<std::complex<float>> const &A, matrix<std::complex<float>> const &B, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<double> const &A, matrix<double> const &B, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<float> const &A, matrix<float> const &B, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: eig(matrix<std::complex<double>> const &A, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<std::complex<float>> const &A, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<double> const &A, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<float> const &A, std::string typ)
   :project: castor

See :ref:`label-qr`, :ref:`label-svd`.

.. _label-qrsvd:

qrsvd
-----
.. doxygenfunction:: qrsvd(matrix<T> const &A, matrix<T> const &B, float tol)
   :project: castor

See :ref:`label-svd`, :ref:`label-qr`.

.. _label-rank:

rank
----

.. doxygenfunction:: rank()
   :project: castor

See :ref:`label-svd`, :ref:`label-aca`.


.. _label-svd:

svd
---
.. doxygenfunction:: svd(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: svd(matrix<std::complex<double>> const &A, std::string typ)
   :project: castor
.. doxygenfunction:: svd(matrix<double> const &A, std::string typ)
   :project: castor
.. doxygenfunction:: svd(matrix<std::complex<float>> const &A, std::string typ)
   :project: castor
.. doxygenfunction:: svd(matrix<float> const &A, std::string typ)
   :project: castor

See :ref:`label-rank`, :ref:`label-eig`, :ref:`label-qr`, :ref:`label-qrsvd`, :ref:`label-aca`.
Getting started
===============

This section gives an example of use of **matrix** library and the procedure to compile it.

Example
-------

.. code::

    #include <iostream>
    #include "castor/matrix.hpp"

    using namespace castor;

    int main (int argc, char* argv[])
    {
        matrix<float> A = {{ 1.0,  2.0,  3.0,  4.0},
                           { 5.0,  6.0,  7.0,  8.0},
                           { 9.0, 10.0, 11.0, 12.0}};
    
        matrix<double> B = eye(3,4);

        auto C = A + B;

        disp(C);
    
        return 0;
    }

This example displays the sum of two matrices with implicit cast :

.. code:: text

    2.00000      2.00000      3.00000      4.00000
    5.00000      7.00000      7.00000      8.00000
    9.00000     10.00000     12.00000     12.00000


.. _label-compilation:

Compilation 
-----------

Command line
++++++++++++

The library **matrix** is a header-only library. That actually means to compile program using **matrix**, you just have to specify to the compiler the path to the directory containing the headers. For instance, with GCC, the command to compile the above example (assuming it is contained in a file named ``main.cpp``) : 

.. code:: bash

    g++ -std=c++14 -I /path/to/castor/folder main.cpp -o main

IDE
+++

Simply enter the path to the directory containing the ``castor`` folder in your favorite IDE (VSCode, Xcode, CodeBlocks). 

`Video tutorial for XCode <http://www.cmapx.polytechnique.fr/~aussal/tmp/matrix_xcode.mov>`_
.. _label-graphical-tools:

Graphical tools
+++++++++++++++

.. _label-caxis:

caxis
-----
.. doxygenfunction:: caxis(figure &fig, matrix<double> const &val)
    :project: castor

.. _label-drawnow:

drawnow
-------
.. doxygenfunction:: drawnow(figure &fig)
    :project: castor

.. _label-builders:

Builders
++++++++

In this section, the user will find all the possible builders other than the constructors for the ``matrix<>`` object. For the latter, please refer to constructors list in the :ref:`class matrix description <label-class-matrix>`. For example, :ref:`label-eye` will return the idendity matrix while :ref:`label-ones` will return a matrix filled with ones.

.. _label-colon:

colon
-----
.. doxygenfunction:: colon(U j, V k)
   :project: castor
.. doxygenfunction:: colon(U j, V i, W k)
   :project: castor

See :ref:`label-linspace`.

.. _label-diag:

diag
----
.. doxygenfunction:: diag
   :project: castor

See :ref:`label-eye`.

.. _label-eye:

eye
---
.. doxygenfunction:: eye(matrix<std::size_t> const &S)
   :project: castor
.. doxygenfunction:: eye(std::size_t m, long n = -1)
   :project: castor

See :ref:`label-ones`, :ref:`label-zeros`, :ref:`label-rand`.

.. _label-linspace:

linspace
--------
.. doxygenfunction:: linspace
   :project: castor

See :ref:`label-logspace`, :ref:`label-colon`.

.. _label-logspace:

logspace
--------
.. doxygenfunction:: logspace
   :project: castor

See :ref:`label-linspace`, :ref:`label-colon`, :ref:`label-log`.

.. _label-rand:

rand
----
.. doxygenfunction:: rand(matrix<std::size_t> const &S, bool seed = false)
   :project: castor
.. doxygenfunction:: rand(std::size_t m, long n = -1, bool seed = false)
   :project: castor

See :ref:`label-zeros`, :ref:`label-eye`, :ref:`label-ones`.

.. _label-ones:

ones
----
.. doxygenfunction:: ones(matrix<std::size_t> const &S)
   :project: castor
.. doxygenfunction:: ones(std::size_t m, long n = -1)
   :project: castor

See :ref:`label-zeros`, :ref:`label-eye`, :ref:`label-rand`.

.. _label-zeros:

zeros
-----
.. doxygenfunction:: zeros(matrix<std::size_t> const &S)
   :project: castor
.. doxygenfunction:: zeros(std::size_t m, long n = -1)
   :project: castor

See :ref:`label-ones`, :ref:`label-eye`, :ref:`label-rand`.
.. _label-geometry:

Geometry
++++++++

The functions from the **Geometry** category are useful to perform transformations between cartesian and polar or spherical coordinates. It is also possible to compute a cartesian 2-dimensional grid with :ref:`label-meshgrid` or to scatter three-dimensional nodes over a sphere with the :ref:`label-sphere` and :ref:`label-sphere2` (Fibonacci sphere) functions.


.. _label-cart2pol:

cart2pol
--------
.. doxygenfunction:: cart2pol(R const &X, S const &Y)
   :project: castor
.. doxygenfunction:: cart2pol(matrix<R> const &X, matrix<S> const &Y)
   :project: castor

See :ref:`label-pol2cart`, :ref:`label-cart2sph`, :ref:`label-sph2cart`.

.. _label-cart2sph:

cart2sph
--------
.. doxygenfunction:: cart2sph(Q const &X, R const &Y, S const &Z)
   :project: castor
.. doxygenfunction:: cart2sph(matrix<Q> const &X, matrix<R> const &Y, matrix<S> const &Z)
   :project: castor

See :ref:`label-sph2cart`, :ref:`label-cart2pol`, :ref:`label-pol2cart`.

.. _label-idx2sph:

idx2sph
--------
.. doxygenfunction:: idx2sph
   :project: castor

See :ref:`label-sph2idx`, :ref:`label-sphere`.

.. _label-meshgrid:

meshgrid
--------
.. doxygenfunction:: meshgrid
   :project: castor

See :ref:`label-kron`.

.. _label-pol2cart:

pol2cart
--------
.. doxygenfunction:: pol2cart(R const &THE, S const &RHO)
   :project: castor
.. doxygenfunction:: pol2cart(matrix<R> const &THE, matrix<S> const &RHO)
   :project: castor

See :ref:`label-cart2pol`, :ref:`label-cart2sph`, :ref:`label-sph2cart`.

.. _label-sph2cart:

sph2cart
--------
.. doxygenfunction:: sph2cart(Q const &THE, R const &PHI, S const &RHO)
   :project: castor
.. doxygenfunction:: sph2cart(matrix<Q> const &THE, matrix<R> const &PHI, matrix<S> const &RHO)
   :project: castor

See :ref:`label-cart2sph`, :ref:`label-cart2pol`, :ref:`label-pol2cart`.

.. _label-sph2idx:

sph2idx
--------
.. doxygenfunction:: sph2idx
   :project: castor

See :ref:`label-idx2sph`, :ref:`label-sphere`.

.. _label-sphere:

sphere
------
.. doxygenfunction:: sphere
   :project: castor

See :ref:`label-idx2sph`, :ref:`label-sph2idx`.

.. _label-sphere2:

sphere2
-------
.. doxygenfunction:: sphere2
   :project: castor

See :ref:`label-sphere`.
.. _label-class-hmatrix:

Class hmatrix
+++++++++++++

As part of **castor**, the ``hmatrix`` class provides a **fully abstract interface** for the manipulation of hierarchical matrices [1]_. The class itself and the related functions can be found in :ref:`label-class-hmatrix`. It uses another class :ref:`label-class-bintree` dedicated to binary space partitioning. 


**This is work-in-progress! The descriptions of the different features and examples will be added to this section progressively.**

.. [1] S. Börm, L. Grasedyck and W. Hackbusch, *Hierarchical Matrices*, 2015 (technical report)

.. doxygenclass:: castor::hmatrix
   :project: castor
   :members:
   :undoc-members:

API
===

.. _label-full-hmatrix:

full
----
.. doxygenfunction:: full(hmatrix<T> const &Ah)
    :project: castor

.. _label-gmres-hmatrix:

gmres
-----
.. doxygenfunction:: gmres(hmatrix<T> const &Ah, matrix<T> const &B, double tol, std::size_t maxit, hmatrix<T> const &Lh, hmatrix<T> const &Uh, matrix<T> const &X0 = matrix<T>())
    :project: castor
.. doxygenfunction:: gmres(hmatrix<T> const &Ah, matrix<T> const &B, double tol, std::size_t maxit, hmatrix<T> const &Ahm1, matrix<T> const &X0 = matrix<T>())
    :project: castor
.. doxygenfunction:: gmres(hmatrix<T> const &Ah, matrix<T> const &B, double tol = 1e-6, std::size_t maxit = 10, std::function<matrix<T>(matrix<T> const&)> const &Am1 = std::function<matrix<T>(matrix<T> const&)>(), matrix<T> const &X0 = matrix<T>())
    :project: castor

.. _label-inv-hmatrix:

inv
---
.. doxygenfunction:: inv(hmatrix<T> const &Ah)
    :project: castor


.. _label-linsolve-hmatrix:

linsolve
--------
.. doxygenfunction:: linsolve(hmatrix<T> const &Ah, matrix<T> const &B)
    :project: castor


.. _label-lu-hmatrix:

lu
--
.. doxygenfunction:: lu(hmatrix<T> const &Ah, double tol = 0)
    :project: castor


.. _label-mtimes-hmatrix:

mtimes
------
.. doxygenfunction:: mtimes(hmatrix<T> const &Ah, hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: mtimes(matrix<T> const &A, hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: mtimes(hmatrix<T> const &Ah, matrix<T> const &B)
    :project: castor


.. _label-size_hmatrix:

size
----
.. doxygenfunction:: size(hmatrix<T> const &Ah)
    :project: castor

See :ref:`label-length`, :ref:`label-numel`.

.. _label-spydata:

spydata
-------
.. doxygenfunction:: spydata(hmatrix<T> const &Ah)
    :project: castor

.. _label-tgeabm-hmatrix:

tgeabm
------
.. doxygenfunction:: tgeabm(T alpha, matrix<T> const &A, matrix<T> const &B, T beta, hmatrix<T> &Ch)
    :project: castor


.. _label-tgemm-hmatrix:

tgemm
-----
.. doxygenfunction:: tgemm(T alpha, hmatrix<T> const &Ah, hmatrix<T> const &Bh, T beta, hmatrix<T> &Ch)
    :project: castor


.. _label-transpose-hmatrix:

transpose
---------
.. doxygenfunction:: transpose(hmatrix<T> const &Ah)
    :project: castor


.. _label-operators-hmatrix:

operators
---------

.. _label-operator<<-hmatrix:

operator<<
++++++++++

.. doxygenfunction:: operator<<(std::ostream &flux, hmatrix<T> const &Ah)
    :project: castor


.. _label-operator+-hmatrix:

operator+
++++++++++

.. doxygenfunction:: operator+(T a, hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: operator+(hmatrix<T> const &Ah, hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: operator+(hmatrix<T> const &Ah, T b)
    :project: castor


.. _label-operator--hmatrix:

operator-
++++++++++

.. doxygenfunction:: operator-(T a, hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: operator-(hmatrix<T> const &Ah, hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: operator-(hmatrix<T> const &Ah, T b)
    :project: castor
.. doxygenfunction:: operator-(hmatrix<T> const &Ah)
    :project: castor


.. _label-operator*-hmatrix:

operator*
++++++++++

.. doxygenfunction:: operator*(T a, hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: operator*(hmatrix<T> const &Ah, T b)
    :project: castor


.. _label-operator/-hmatrix:

operator/
++++++++++

.. doxygenfunction:: operator/(hmatrix<T> const &Ah, T b)
    :project: castor


.. _label-operator+=-hmatrix:

operator+=
++++++++++
.. doxygenfunction:: operator+=(hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: operator+=(T b)
    :project: castor



.. _label-operator*=-hmatrix:

operator*=
++++++++++
.. doxygenfunction:: operator*=(T b)
    :project: castor
.. _label-mesh-management:

Mesh management
+++++++++++++++

Functions helpers to create simple volume or surface meshes from sets of nodes.

.. _label-tetboundary:

tetboundary
-----------
.. doxygenfunction:: tetboundary(matrix<std::size_t> const &tet, matrix<T> const &vtx)
    :project: castor

See :ref:`label-tetdelaunay`.


.. _label-tetdelaunay:

tetdelaunay
-----------
.. doxygenfunction:: tetdelaunay(matrix<T> const &X, matrix<T> const &Y, matrix<T> const &Z)
    :project: castor

See :ref:`label-tetboundary`.


.. _label-tridelaunay:

tridelaunay
-----------
.. doxygenfunction:: tridelaunay(matrix<T> const &X, matrix<T> const &Y, matrix<T> const &Z = {})
    :project: castor
.. _label-analyticalscattering:

Analytical Scattering
=====================

The project can be found `here <https://github.com/marcbakry/AnalyticalScattering/>`_... _label-using-kissfft:

Fast Fourier Transform
======================

The **castor** framework wraps the `kissfft library <https://github.com/mborgerding/kissfft>`_ in a standalone ``kissfft.hpp`` header file so it is not necessary to download it. Currently, only **single precision** arithmetic is supported. We describe below a basic example.

A minimum working file can be found below.

.. code:: c++

    #include "castor/matrix.hpp"
    #include "castor/kissfft.hpp"

    int main()
    {
        matrix<std::complex<float>> A = rand(3,4) + M_1I*rand(3,4);
        // your code here
        return 0;
    }

In the code above, we created a random complex single precision matrix with twelve elements. Since only the one-dimensional fft is supported, the matrix ``A`` will be treated as a one-dimensional array where the first row comes first, etc. Now we declare a second matrix ``B`` where we store the result of the forward-backward discrete Fourier transform.

.. code:: c++

    matrix<std::complex<float>> ifft(fft(A));

While ``A`` is a two-dimensional array, this is not the case for ``B``. In order to access the one-dimensional data of ``A``, we use a view on the 12 elements using ranged access (see :ref:`label-view`) and we compare it to ``B``.

.. code:: c++

    disp(norm(eval(A(range(0,12))) - B, "inf");

The output should be something more or less like

.. code:: text

    (9.31323e-08,0)

Please note that this code is part of the ``overview_kissfft.cpp`` file which can be found in the ``demo/demo_kissfft`` subfolder... _label-linear-algebra-advanced:

Linear algebra
==============

We provide a set of useful function to perform basic linear algebra manipulations available by including ``castor/linalg.hpp``. The **linear algebra** library is mostly based on the well-known BLAS and LAPACK APIs. As such, it may be interfaced with any library respecting the same naming convention like ``openblas`` or ``mkl``. 

We want to bring the attention of the user on the fact that BLAS and LAPACK are originally Fortran libraries and consequently use a column ordering storage convention, unlike ``matrix`` which uses a row-major storage. While this behavior is fully transparent for the user as long as he is using the high-level **castor** API, he may find more details at :ref:`label-blaslapack-issue`.

We give below a few examples of use. More may be found in the ``demo/demo_linalg/`` subfolder. All linear algebra functions are described at :ref:`label-singular-eig-values-func`, :ref:`label-factorization-func` and :ref:`label-linear-solver-func`



Singular Value Decomposition
----------------------------
It is very easy to compute the SVD of a given matrix. The function ``svd`` simply returns a tuple containing the singular values ``S``, the left singular vectors ``U`` and finally the transposed right singular vectors ``Vt``.

.. code:: c++

    matrix<double> A({
        {1,1,-2,1,3,-1},
        {2,-1,1,2,1,-3},
        {1,3,-3,-1,2,1},
        {5,2,-1,-1,2,1},
        {-3,-1,2,3,1,3},
        {4,3,1,-6,-3,-2}
    });
    matrix<> S, U, Vt;
    std::tie(S,U,Vt) = svd(A,"vect");
    disp(norm(A - mtimes(U,mtimes(diag(S),Vt)),"inf"));

.. code:: text

    4.44089e-15

Note that it is equally easy to compute the (eventually generalized) eigenvalues and eigenvectors (see :ref:`label-eig`) or the rank-revealing SVD (see :ref:`label-qrsvd`).


Solving a linear system
---------------------------
We provide two functions for the resolution of a linear system of equations. ``linsolve`` performs the exact inversion using the LAPACK library while ``gmres`` computes the solution iteratively using the GMRES algorithm. Note that **both** have support for multiple right-hand-sides.

The use of ``linsolve`` is straightforward. The left- and right-hand-sides should be provided in full ``matrix`` form and the function returns the ``matrix`` containing the solutions.

.. code:: c++

    // left-hand-side
    matrix<double> A({
        {1,1,-2,1,3,-1},
        {2,-1,1,2,1,-3},
        {1,3,-3,-1,2,1},
        {5,2,-1,-1,2,1},
        {-3,-1,2,3,1,3},
        {4,3,1,-6,-3,-2}
    });
    // we change the original matrix for this demo
    for(auto i=0; i<6; ++i) A(i,i) = 20.; 
    // right-hand-side, must be a column
    matrix<double> B(6,1,{4,20,-15,-3,16,-27});
    // solving
    auto X = linsolve(A,B);
    disp(X);

.. code:: text

   -0.15882  
    0.83718  
   -0.93171  
   -0.29017  
    1.15147  
   -1.31156

Compared to  ``linsolve``, one needs to specifiy the tolerance and the maximum number of iterations. The user may also provide a preconditionner, use its own definition of the matrix-vector product, etc., as specified in the documentation of :ref:`label-gmres`.

.. code:: c++

    X = gmres(A,B,1e-3,6);
    disp(X);

.. code:: text

    Start GMRES using MGCR implementation (Multiple Generalized Conjugate Residual):
    + Iteration 1 in 3.5977e-05 seconds with relative residual 0.296391.
    + Iteration 2 in 2.8044e-05 seconds with relative residual 0.0545852.
    + Iteration 3 in 2.5446e-05 seconds with relative residual 0.0119836.
    + Iteration 4 in 3.3522e-05 seconds with relative residual 0.000506301.
    GMRES converged at iteration 4 to a solution with relative residual 0.000506301.
       -0.15823  
        0.83694  
       -0.93228  
       -0.29023  
        1.15188  
       -1.31119

**Remark:** Exact convergence is always achieved when the number of iterations reaches the dimension of the matrix.


Adaptive Cross Approximation
----------------------------
The ACA algorithm performs a low-rank approximation of rank ``k`` a given matrix ``C`` of size ``m x n`` under the form ``C = mtimes(A,tranpose(B))`` where ``A`` is ``m x k`` and ``B`` is ``n x k``. The rank ``k`` is automatically obtained following a user-defined accuracy parameter ``tol``. Note that if ``tol`` is chosen too small, using the ACA algorithm may be useless ...
**Warning** : our implementation directly returns ``Bt = tranpose(B)``.

.. code:: c++

    matrix<> A, Bt;
    matrix<> C = mtimes(rand(100,20),rand(20,50));
    std::tie(A,Bt) = aca(C,1e-3);
    disp(size(A),2);
    disp(size(Bt),2);
    disp(norm(C - mtimes(A,Bt),"inf"),2);

.. code:: text

    Matrix 1x2 of type 'm' (16 B):
    100   21  
    Matrix 1x2 of type 'm' (16 B):
    21  50  
    Object of type 'd':
    3.81917e-14

**Remark:** The result may vary depending on your random number generator.



.. _label-blaslapack-issue:

Overload of some BLAS and LAPACK functions
------------------------------------------

Some functions from the BLAS and LAPACK APIs have been directly interfaced using a similar naming convention. 

The BLAS-part is in fact a straightforward overlay over the C BLAS API which enables the possibility to use row-major ordering. For example, ``tgemm`` computes the matrix-matrix product using the C BLAS API (see :ref:`label-tgemm-blas`) and provides a unique interface to the four corresponding functions (each one corresponding to one of the four types: ``float``, ``std::complex<float>``, ``double``, ``std::complex<double>``). Its main purpose is to hide the use of raw pointers to the underlying data. For information, such raw pointers can be obtained as explained below.

.. code:: c++

    matrix<double> A({1,2,3,4});
    double *pA = &A(0); // returns the adress of the first element of A

However, the LAPACK-part is a *direct* overlay over the Fortran LAPACK API. Under the hood, the functions of **castor** convert the ``matrix<>`` to a ``std::vector<>`` where the data is stored with the right ordering thanks to the ``mat2lpk`` function (see :ref:`label-mat2lpk`). The result is then converted back to a ``matrix<>`` with the ``lpk2mat`` function (see :ref:`label-lpk2mat`). The main consequence is that the functions following a naming convention close to the LAPACK one (see for example :ref:`label-tgesdd`, :ref:`label-tgeqrf`, etc.) only accept ``std::vector<>`` as input.
.. _label-tools:

Tools
+++++

In this section, we describe some useful tools provided with the **castor** framework. For example, it is possible to produce formatted output with :ref:`label-disp`. It is also possible to get :ref:`label-help` or measure the execution time with :ref:`label-tic` and :ref:`label-toc`.


.. _label-disp:

disp
----
.. doxygenfunction:: disp(matrix<T> const &A, int info = 2, std::ostream &flux = std::cout, std::size_t m = 3, std::size_t n = 3)
   :project: castor

.. doxygenfunction:: disp(T A, int info = 1, std::ostream &flux = std::cout)
   :project: castor

See :ref:`label-help`, :ref:`label-error`, :ref:`label-warning`.

.. _label-help:

help
-----
.. doxygenfunction:: help
   :project: castor

See :ref:`label-disp`, :ref:`label-error`, :ref:`label-warning`.

.. _label-error:

error
-----
.. doxygenfunction:: error
   :project: castor

See :ref:`label-warning`, :ref:`label-disp`, :ref:`label-help`.

.. _label-tic:

tic
---
.. doxygenfunction:: tic
   :project: castor

See :ref:`label-toc`.

.. _label-toc:

toc
---
.. doxygenfunction:: toc
   :project: castor

See :ref:`label-tic`.

.. _label-warning:

warning
-------
.. doxygenfunction:: warning
   :project: castor

See :ref:`label-error`, :ref:`label-disp`, :ref:`label-help`.
.. _label-smatrix-API:

API
===

.. _label-check-smatrix:

check
-----
.. doxygenfunction:: check(smatrix<T> &As)
    :project: castor

See :ref:`label-class-matrix`.


.. _label-clear-smatrix:

clear
-----
.. doxygenfunction:: clear(smatrix<T> &As)
    :project: castor

See :ref:`label-class-matrix`.


.. _label-disp-smatrix:

disp
----
.. doxygenfunction:: disp(smatrix<T> const &As, int info = 2, std::ostream &flux = std::cout, std::size_t r = 3)
    :project: castor

See :ref:`label-operator<<-smatrix`.


.. _label-find-smatrix:

find
----
.. doxygenfunction:: find(smatrix<T> const &As)
    :project: castor

See :ref:`label-index-smatrix`, :ref:`label-values-smatrix`.


.. _label-full-smatrix:

full
----
.. doxygenfunction:: full(smatrix<T> const &As)
    :project: castor

See :ref:`label-sparse-smatrix`.


.. _label-index-smatrix:

index
-----
.. doxygenfunction:: index(smatrix<T> const &As)
    :project: castor

.. _label-mtimes-smatrix:

mtimes
------
.. doxygenfunction:: mtimes(matrix<R> const &A, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: mtimes(smatrix<R> const &As, matrix<S> const &B)
    :project: castor

See :ref:`label-tgemm-naive`, :ref:`label-kron`.


.. _label-nnz-smatrix:

nnz
---
.. doxygenfunction:: nnz(smatrix<T> const &As)
    :project: castor

See :ref:`label-find`, :ref:`label-size-smatrix`.


.. _label-numel-smatrix:

numel
-----
.. doxygenfunction:: numel(smatrix<T> const &As)
    :project: castor

See :ref:`label-size-smatrix`, :ref:`label-nnz-smatrix`.


.. _label-reshape-smatrix:

reshape
-------
.. doxygenfunction:: reshape(smatrix<T> const &As, std::size_t m, std::size_t n)
    :project: castor

See :ref:`label-transpose-smatrix`.

.. _label-size-smatrix:

size
----
.. doxygenfunction:: size(smatrix<T> const &As, int dim)
    :project: castor
.. doxygenfunction:: size(smatrix<T> const &As)
    :project: castor

See :ref:`label-numel-smatrix`, :ref:`label-nnz-smatrix`.


.. _label-sparse-smatrix:

sparse
------

.. doxygenfunction:: sparse(matrix<std::size_t> const &L, matrix<T> const &V, std::size_t m, std::size_t n)
    :project: castor
.. doxygenfunction:: sparse(matrix<std::size_t> const &I, matrix<std::size_t> const &J, matrix<T> const &V)
    :project: castor
.. doxygenfunction:: sparse(matrix<std::size_t> const &I, matrix<std::size_t> const &J, matrix<T> const &V, std::size_t m, std::size_t n)
    :project: castor
.. doxygenfunction:: sparse(matrix<T> const &A)
    :project: castor

See :ref:`label-full-smatrix`.


.. _label-speye-smatrix:

speye
-----
.. doxygenfunction:: speye(matrix<std::size_t> const &S)
    :project: castor
.. doxygenfunction:: speye(std::size_t m, long n = -1)
    :project: castor

See :ref:`label-spzeros-smatrix`, :ref:`label-spones-smatrix`, :ref:`label-sprand-smatrix`, :ref:`label-eye`.


.. _label-spones-smatrix:

spones
------
.. doxygenfunction:: spones(matrix<std::size_t> const &S)
    :project: castor
.. doxygenfunction:: spones(std::size_t m, long n = -1)
    :project: castor

See :ref:`label-spzeros-smatrix`, :ref:`label-speye-smatrix`, :ref:`label-sprand-smatrix`, :ref:`label-ones`.


.. _label-sprand-smatrix:

sprand
------
.. doxygenfunction:: sprand(matrix<std::size_t> const &S, bool seed = false)
    :project: castor
.. doxygenfunction:: sprand(std::size_t m, long n = -1, bool seed = false)
    :project: castor

See :ref:`label-spzeros-smatrix`, :ref:`label-speye-smatrix`, :ref:`label-spones-smatrix`, :ref:`label-rand`.


.. _label-spzeros-smatrix:

spzeros
-------
.. doxygenfunction:: spzeros(std::size_t m, long n = -1)
    :project: castor

See :ref:`label-spones-smatrix`, :ref:`label-speye-smatrix`, :ref:`label-sprand-smatrix`, :ref:`label-zeros`.


.. _label-transpose-smatrix:

transpose
---------
.. doxygenfunction:: transpose(smatrix<T> const &As)
    :project: castor

See :ref:`label-reshape-smatrix`.


.. _label-values-smatrix:

values
------
.. doxygenfunction:: values(smatrix<T> const &As)
    :project: castor

See :ref:`label-index-smatrix`, :ref:`label-find-smatrix`.



.. _label-operator<<-smatrix:

operator<<
----------
.. doxygenfunction:: operator<<(std::ostream &flux, smatrix<T> const &As)
    :project: castor

See :ref:`label-disp-smatrix`.


.. _label-operator+-smatrix:

operator+
---------
.. doxygenfunction:: operator+(smatrix<R> const &As, S b)
    :project: castor
.. doxygenfunction:: operator+(R a, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator+(smatrix<R> const &As, matrix<S> const &B)
    :project: castor
.. doxygenfunction:: operator+(matrix<R> const &A, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator+(smatrix<R> const &As, smatrix<S> const &Bs)
    :project: castor



.. _label-operator--smatrix:

operator-
---------
.. doxygenfunction:: operator-(smatrix<R> const &As, S b)
    :project: castor
.. doxygenfunction:: operator-(R a, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator-(smatrix<R> const &As, matrix<S> const &B)
    :project: castor
.. doxygenfunction:: operator-(matrix<R> const &A, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator-(smatrix<R> const &As, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator-(smatrix<T> const &As)
    :project: castor



.. _label-operator*-smatrix:

operator*
---------
.. doxygenfunction:: operator*(smatrix<R> const &As, S Bs)
    :project: castor
.. doxygenfunction:: operator*(R As, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator*(smatrix<R> const &As, matrix<S> const &B)
    :project: castor
.. doxygenfunction:: operator*(matrix<R> const &A, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator*(smatrix<R> const &As, smatrix<S> const &Bs)
    :project: castor


.. _label-operator/-smatrix:

operator/
---------
.. doxygenfunction:: operator/(smatrix<R> const &As, S Bs)
    :project: castor
.. doxygenfunction:: operator/(matrix<R> const &As, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator/(smatrix<R> const &As, matrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator/(R As, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator/(smatrix<R> const &As, smatrix<S> const &Bs)
    :project: castor


.. _label-factorization-func:

Factorization
+++++++++++++

This section regroups all of the linear algebra functions concerning factorization which are currently implemented within **castor**. In order to access these functions, **the user must include** ``castor/linalg.hpp``. Two levels of interface are available. The high-level interface provides functions with simplified arguments while the low-level interface is much closer to the BLAS/LAPACK API.

Some examples of use may also be found at :ref:`label-linear-algebra-advanced`.


.. _label-aca:

aca
---
.. doxygenfunction:: aca(matrix<T> const &A, matrix<T> const &B, double tol = 1e-6, std::size_t rmax = 1e6)
   :project: castor
.. doxygenfunction:: aca(matrix<T> const &M, double tol = 1e-6, std::size_t rmax = 1e6)
   :project: castor
.. doxygenfunction:: aca(matrix<std::size_t> I, matrix<std::size_t> J, std::function<matrix<std::complex<double>>(matrix<std::size_t>, matrix<std::size_t>)> const &fct, double tol = 1e-6, std::size_t rmax = 1e6)
   :project: castor
.. doxygenfunction:: aca(matrix<std::size_t> I, matrix<std::size_t> J, std::function<matrix<double>(matrix<std::size_t>, matrix<std::size_t>)> const &fct, double tol = 1e-6, std::size_t rmax = 1e6)
   :project: castor

See :ref:`label-svd`, :ref:`label-rank`.


.. _label-lu:

lu
--
.. doxygenfunction:: lu(matrix<std::complex<double>> const &A)
   :project: castor
.. doxygenfunction:: lu(matrix<double> const &A)
   :project: castor
.. doxygenfunction:: lu(matrix<std::complex<float>> const &A)
   :project: castor
.. doxygenfunction:: lu(matrix<float> const &A)
   :project: castor

See :ref:`label-qr`, :ref:`label-linsolve`, :ref:`label-inv`.

.. _label-qr:

qr
--
.. doxygenfunction:: qr(matrix<std::complex<double>> const &A)
   :project: castor
.. doxygenfunction:: qr(matrix<double> const &A)
   :project: castor
.. doxygenfunction:: qr(matrix<std::complex<float>> const &A)
   :project: castor
.. doxygenfunction:: qr(matrix<float> const &A)
   :project: castor

See :ref:`label-eig`, :ref:`label-svd`, :ref:`label-lu`.
.. _label-mesh-plot:

Mesh plot
+++++++++


.. _label-edgmesh:

edgmesh
-------
.. doxygenfunction:: edgmesh(figure &fig, matrix<std::size_t> const &edg, matrix<T> const &vtx, matrix<T> const &val = {})
    :project: castor


.. _label-mesh:

mesh
----
.. doxygenfunction:: mesh(figure &fig, matrix<T> const &X, matrix<T> const &Y, matrix<T> const &Z, std::string const &options = "")
    :project: castor


.. _label-quiver:

quiver
------
.. doxygenfunction:: quiver(figure &fig, matrix<T> const &vtx, matrix<T> const &dir, matrix<T> const &val = {})
    :project: castor


.. _label-tetmesh:

tetmesh
-------
.. doxygenfunction:: tetmesh(figure &fig, matrix<std::size_t> const &tet, matrix<T> const &vtx, matrix<T> const &val = {})
    :project: castor


.. _label-trimesh:

trimesh
-------
.. doxygenfunction:: trimesh(figure &fig, matrix<std::size_t> const &tri, matrix<T> const &vtx, matrix<T> const &val = {})
    :project: castor


.. _label-vermesh:

vermesh
-------
.. doxygenfunction:: vermesh(figure &fig, matrix<std::size_t> const &ver, matrix<T> const &vtx, matrix<T> const &val = {})
    :project: castor
Monte Carlo method
==================
*Shared by Antoine Rideau*

On this page you will find how to calculate using **Castor** the value of :math:`\pi` with Monte Carlo method.

Monte Carlo methods are a broad class of computational algorithms that rely on repeated random sampling to obtain numerical results.
The underlying concept is to use randomness to solve problems that might be deterministic in principle.
They are mainly used in physical and mathematical problems related to optimization, numerical integration, and generating draws from a probability distribution.

First, draw 1 by 1 square and inscribe a quadrant within it.


Then ``n`` points are scattered randomly  with :math:`(x,y)` coordinate over the square.

.. code-block:: c++

    // Parameters
    double n = 5000.; //Number of points
    auto MX = rand(1, n); 
    auto MY = rand(1, n);


As the quadrant's surface is equal to

.. math::

    \sigma = \frac{R^2 \pi}{4} = \frac{\pi}{4} \text{ as } R = 1 ,

and the square's surface is 

.. math::

    S = R^2 = 1 ,

the probability to be in the quadrant is

.. math::

    \mathbb{P} = \frac{N_{in}}{n} = \frac{\sigma}{S} = \frac{\pi}{4} ,

| with :math:`N_{in}` number of points inside the quadrant.
| Hence,

.. math::

    \pi = 4 * \frac{N_{in}}{n} .


So, points inside the quadrant i.e points with :math:`x^2 + y^2 \leq 1` are counted.

.. code-block:: c++
    
    // Pi computation
    matrix<std::size_t> Incircle;
    Incircle = find(pow(MX, 2) + pow(MY, 2) <= 1);
    double Pi = 4. * (numel(Incircle) / n);

See :ref:`label-find-smatrix` , :ref:`label-numel` 

Another way of doing it is by using ``sum``

.. code-block:: c++

    size_t Nin = sum<size_t>(pow(MX, 2) + pow(MY, 2) <= 1);
    double Pi = 4. * (Nin / n);

See :ref:`label-sum`

Code
----

.. code-block:: c++

    #include "castor/matrix.hpp"
    #include "castor/graphics.hpp"

    using namespace castor;

    int main(int argc, char const *argv[])
    {

        // Parameters
        double n = 5000.; //Number of points
        auto MX = rand(1, n);
        auto MY = rand(1, n);

        // Pi computation
        matrix<std::size_t> Incircle;
        Incircle = find(pow(MX, 2) + pow(MY, 2) <= 1);
        double Pi = 4. * (numel(Incircle) / n);
        std::cout << "Calculated value of pi: " << Pi << endl;

        // Visu
        auto X = linspace(0, 1, 1000);
        auto Y = sqrt(1 - pow(X, 2));
        figure fig;
        plot(fig, MX, MY, {"b"});
        plot(fig, eval(MX(Incircle)), eval(MY(Incircle)), {"r"});
        plot(fig, X, Y, {"r-"});
        drawnow(fig);

        return 0;
    }

With this code you should get these outputs :

.. code-block:: text

    Calculated value of pi: 3.148


.. image:: img/montecarloresult.png
    :width: 400
    :align: center


References
----------

https://en.wikipedia.org/wiki/Monte_Carlo_method
.. _label-developpers:

Developers
==========

- Matthieu Aussal, research engineer at `CMAP <cmap.polytechnique.fr>`_, École polytechnique : matthieu.aussal at polytechnique.edu
- Marc Bakry, research engineer at `CEA <https://www.cea.fr>`_, (former post-doc at Ecole polytechnique) : marc.bakry at gmail dot com
- Laurent Series, research engineer at `CMAP <cmap.polytechnique.fr>`_, École polytechnique : laurent.series at polytechnique.edu
.. _label-mathematical-functions:

Mathematical functions
++++++++++++++++++++++

The functions described below implement common mathematical functions which can be applied element-wise to the elements of a ``matrix`` such as :ref:`label-cos`, :ref:`label-sqrt`, etc.


.. _label-abs:

abs
---
.. doxygenfunction:: abs
   :project: castor

See :ref:`label-angle`, :ref:`label-sign`.


.. _label-acos:

acos
----
.. doxygenfunction:: acos
   :project: castor

See :ref:`label-acosd`, :ref:`label-cos`.


.. _label-acosd:

acosd
-----
.. doxygenfunction:: acosd
   :project: castor

See :ref:`label-acos`, :ref:`label-cosd`.


.. _label-acosh:

acosh
-----
.. doxygenfunction:: acosh
   :project: castor

See :ref:`label-cosh`.


.. _label-angle:

angle
-----
.. doxygenfunction:: angle
   :project: castor
   
See :ref:`label-abs`.


.. _label-asin:

asin
----
.. doxygenfunction:: asin
   :project: castor

See :ref:`label-sin`, :ref:`label-asind`.


.. _label-asind:

asind
-----
.. doxygenfunction:: asind
   :project: castor

See :ref:`label-sind`, :ref:`label-asin`.


.. _label-asinh:

asinh
-----
.. doxygenfunction:: asinh
   :project: castor

See :ref:`label-sinh`.


.. _label-atan:

atan
----
.. doxygenfunction:: atan
   :project: castor

See :ref:`label-tan`, :ref:`label-atand`.


.. _label-atand:

atand
-----
.. doxygenfunction:: atand
   :project: castor

See :ref:`label-tand`, :ref:`label-atan`.


.. _label-atanh:

atanh
-----
.. doxygenfunction:: atanh
   :project: castor

See :ref:`label-tanh`.


.. _label-ceil:

ceil
----
.. doxygenfunction:: ceil
   :project: castor

See :ref:`label-floor`, :ref:`label-round`.


.. _label-conj:

conj
----
.. doxygenfunction:: conj(matrix<float> const &A)
   :project: castor
.. doxygenfunction:: conj(matrix<double> const &A)
   :project: castor
.. doxygenfunction:: conj(matrix<S> const &X)
   :project: castor

See :ref:`label-real`, :ref:`label-imag`.


.. _label-cos:

cos
---
.. doxygenfunction:: cos
   :project: castor

See :ref:`label-acos`, :ref:`label-cosd`.


.. _label-cosd:

cosd
----
.. doxygenfunction:: cosd
   :project: castor

See :ref:`label-acosd`, :ref:`label-cos`.


.. _label-cosh:

cosh
----
.. doxygenfunction:: cosh
   :project: castor

See :ref:`label-acosh`.


.. _label-deg2rad:

deg2rad
-------
.. doxygenfunction:: deg2rad(T x)
   :project: castor
.. doxygenfunction:: deg2rad(matrix<T> const &X)
   :project: castor

See :ref:`label-rad2deg`.


.. _label-exp:

exp
---
.. doxygenfunction:: exp
   :project: castor

See :ref:`label-log`, :ref:`label-log10`.


.. _label-floor:

floor
-----
.. doxygenfunction:: floor
   :project: castor

See :ref:`label-ceil`, :ref:`label-round`.


.. _label-imag:

imag
----
.. doxygenfunction:: imag(matrix<float> const &A)
   :project: castor
.. doxygenfunction:: imag(matrix<double> const &A)
   :project: castor
.. doxygenfunction:: imag(matrix<S> const &X)
   :project: castor

See :ref:`label-real`, :ref:`label-conj`, :ref:`label-angle`, :ref:`label-abs`.


.. _label-log:

log
---
.. doxygenfunction:: log
   :project: castor

See :ref:`label-log2`, :ref:`label-log10`, :ref:`label-exp`.


.. _label-log2:

log2
----
.. doxygenfunction:: log2
   :project: castor

See :ref:`label-log`, :ref:`label-log10`, :ref:`label-exp`.


.. _label-log10:

log10
-----
.. doxygenfunction:: log10
   :project: castor

See :ref:`label-log`, :ref:`label-log2`, :ref:`label-exp`.


.. _label-pow:

pow
---
.. doxygenfunction:: pow(R x, matrix<S> const &Y)
   :project: castor
.. doxygenfunction:: pow(matrix<R> const &X, S y)
   :project: castor
.. doxygenfunction:: pow(matrix<R> const &X, matrix<S> const &Y)
   :project: castor

See :ref:`label-exp`, :ref:`label-log`.


.. _label-rad2deg:

rad2deg
-------
.. doxygenfunction:: rad2deg(T x)
   :project: castor
.. doxygenfunction:: rad2deg(matrix<T> const &X)
   :project: castor

See :ref:`label-deg2rad`.


.. _label-real:

real
----
.. doxygenfunction:: real(matrix<float> const &A)
   :project: castor
.. doxygenfunction:: real(matrix<double> const &A)
   :project: castor
.. doxygenfunction:: real(matrix<S> const &X)
   :project: castor

See :ref:`label-imag`, :ref:`label-conj`, :ref:`label-angle`, :ref:`label-abs`.


.. _label-round:

round
-----
.. doxygenfunction:: round
   :project: castor

See :ref:`label-floor`, :ref:`label-ceil`.


.. _label-sign:

sign
----
.. doxygenfunction:: sign
   :project: castor

See :ref:`label-abs`.


.. _label-sin:

sin
---
.. doxygenfunction:: sin
   :project: castor

See :ref:`label-asin`, :ref:`label-sind`.


.. _label-sind:

sind
----
.. doxygenfunction:: sind
   :project: castor

See :ref:`label-asind`, :ref:`label-sin`.


.. _label-sinh:

sinh
----
.. doxygenfunction:: sinh
   :project: castor

See :ref:`label-asinh`.


.. _label-sqrt:

sqrt
----
.. doxygenfunction:: sqrt
   :project: castor

See :ref:`label-pow`.


.. _label-tan:

tan
---
.. doxygenfunction:: tan
   :project: castor

See :ref:`label-atan`, :ref:`label-tand`.


.. _label-tand:

tand
----
.. doxygenfunction:: tand
   :project: castor

See :ref:`label-atand`, :ref:`label-tan`.


.. _label-tanh:

tanh
----
.. doxygenfunction:: tanh
   :project: castor

See :ref:`label-atanh`.Eigenmodes of Helmholtz
=======================

*Shared by Antoine Rideau*

On this page you will find how to show using **Castor** the eigenmodes of the wave equation on a rectangle with fixed edges.

We start from the following D'Alembert equation on :math:`\Omega = \left [ 0, x_{0} \right ] \times \left [ 0, y_{0} \right ]` with Dirichlet boundary condition on :math:`\Gamma = \partial \Omega` 

.. math:: 

    \left\{\begin{matrix}
    - \displaystyle \frac{1}{c} \frac{\partial^2 u }{\partial t^2}(\mathbf{x}) + \Delta u (\mathbf{x}) = 0 & , & \mathbf{x} \in \Omega \setminus \Gamma
    \\ 
    u(\mathbf{x} , \cdot   ) = 0 & , & \mathbf{x} \in \Gamma 
    \end{matrix}\right.
    ,

considering

.. math::

    u (\mathbf{x},t) = V(\mathbf{x})e^{i \omega t} ,

gives the Helmholtz equation 

.. math::

    \left\{\begin{matrix}
    k^{2}V(\mathbf{x}) + \Delta V(\mathbf{x}) = 0 & , & \mathbf{x} \in \Omega \setminus \Gamma
   \\
   V(\mathbf{x}) = 0 & , & \mathbf{x} \in \Gamma
   \end{matrix}\right. 
    ,

| with :math:`k = \displaystyle \frac{\omega}{c}` .
| 
| The space domain ``L`` is discretized with ``dx`` steps which results in the meshgrid described by ``X`` and ``Y`` 

.. math:: 

    \begin{matrix} x_{i} = i \delta x & \text{ for } i = \left [ \! \left [ 0, n_{x}-1 \right ] \! \right ]\\ y_{j} = j \delta x & \text{ for } j = \left [ \! \left [ 0, n_{y}-1 \right ] \! \right ] \end{matrix}

.. code-block:: c++

    // Parameters
    matrix<> L = {1, 2}; // Dimensions
    double dx = 0.05;    // Space discretization

    // Discretization
    matrix<> X, Y;
    std::tie(X, Y) = meshgrid(colon(0, dx, L(0)), colon(0, dx, L(1)));

See :ref:`label-meshgrid`.

Laplacian
---------

The Helmholtz equation can now be written in vector form 

.. math::

    - K V(\mathbf{x}) = k^{2} V(\mathbf{x}) ,

where ``K`` stands for the matrix of the Laplacian operator.

The Laplacian operator can be approximated as 

.. math::

    \begin{matrix}
    \Delta_{\textbf{x}}u(x,y) & = & \displaystyle \frac{\partial^2 u}{\partial x^2}(x,y) + \frac{\partial^2 u}{\partial y^2}(x,y) 
    \\ 
    \Delta_{\textbf{x}}u_{i,j} & \approx & \displaystyle \frac{u_{i+1,j}+u_{i,j+1}-4u_{i,j}+u_{i-1,j}+u_{i,j-1}}{\delta_{\textbf{x}}^2} & .
    \end{matrix}
    


This expression leads to this form for ``K``

.. math::

    K = \frac{1}{\delta_{\textbf{x}}^2} \begin{pmatrix}
    -4 & 1 & 0 & \cdots & 0 & 1 & 0 & \cdots & 0\\ 
     1 & -4 & 1 & 0 & \cdots & 0 & 1 & \ddots  & \vdots \\ 
     0 & 1  & \ddots & \ddots & \ddots &  & \ddots & \ddots & 0\\ 
     \vdots& \ddots & \ddots & -4 & 1 & 0 &  & \ddots & 1\\ 
    0 &  & 0 & 1 & -4 & 1 & 0 &  & 0\\ 
     1& \ddots &  & 0 & 1 & -4 & \ddots & \ddots & \vdots \\ 
    0 & \ddots & \ddots &  & \ddots & \ddots & \ddots & 1 & 0 \\ 
     \vdots& \ddots & 1 & 0 & \cdots & 0 & 1 & -4 & 1 \\ 
    0 & \cdots & 0 & 1 & 0 & \cdots & 0 & 1 & -4
    \end{pmatrix}
    .

.. code-block:: c++

    // Laplacian
    long nx = size(X, 2), ny = size(X, 1);
    matrix<> e = ones(nx * ny, 5);
    e(row(e), 2) = -4;
    matrix<> K = full(spdiags(1. / (dx * dx) * e, {-nx, -1, 0, 1, nx}, nx * ny, nx * ny));

See :ref:`label-spdiags`

| Here, K is built as a sparse matrix using `spdiags` because it is easier to do so but I convert it into a dense matrix because **Castor** don't have a sparse solver yet.
|
| In order to take into account the homogeneous Dirichlet condition on the boundary, penalization  is used on the index where the boundaries are : index ``i`` such as ``X(i)==0``, ``X(i)==L(0)``, ``Y(i)==0`` and ``Y(i)==L(1)`` .

.. code-block:: c++

    // Penalization on boundary (Homogeneous Dirichlet condition)
    matrix<std::size_t> Ibnd;
    Ibnd = find((X == 0) || (X == L(0)) || (Y == 0) || (Y == L(1)));
    K(sub2ind(size(K), Ibnd, Ibnd)) = 1e6;

See :ref:`label-find-smatrix`, :ref:`label-sub2ind`.
    
Analytical solution
-------------------

An eigenmodes is caracterize by 2 positive integers :math:`m` and :math:`n` . Thus the eigenvalues are 

.. math:: 

    \lambda_{m,n} = c\pi \sqrt{\frac{m^2}{x_{0}}+\frac{n^2}{y_{0}}}

and the corresponding eigenmode are 

.. math::

    u_{m,n} = \sin \bigg(\frac{m\pi x}{x_{0}}\bigg) \sin \bigg(\frac{n\pi y}{y_{0}}\bigg)


.. code-block:: c++

    // Analytical
    auto Dth = zeros(nx, ny);
    for (int m = 0; m < nx; m++)
    {
        for (int n = 0; n < ny; n++)
        {
            Dth(m, n) = M_PI * sqrt(pow((m + 1) / L(0), 2) + pow((n + 1) / L(1), 2));
        }
    }

See :ref:`label-zeros` .

Eigenmodes
-----------

Once the Laplacian matrix have been built , eigenvalues are easily acquired  in the ``1`` by ``nx*ny`` vector ``D`` and eigenvectors in the ``nx*ny`` by ``nx*ny`` matrix ``V`` using the ``eig`` function

.. math:: 

    - K V(\mathbf{x}) = D V(\mathbf{x}) 

.. code-block:: c++

    // Numerical eigen values and vectors
    matrix<std::complex<double>> D, V;
    std::tie(D, V) = eig(-K, "right");

See :ref:`label-eig` .

The eigenvalues considerated are these with an imaginary part null and a real part minimal. To do so eigenvalues and eigenvectors are sorted by ascending eigenvalues.

.. code-block:: c++

    // Sort
    matrix<std::size_t> I;
    I = argsort(abs(real(D)));
    D = eval(D(I));
    V = eval(V(row(V), I));
    matrix<std::size_t> Ith;
    Ith = argsort(Dth);
    Dth = eval(Dth(Ith));

See :ref:`label-argsort` , :ref:`label-row` . 

Then for each eigenmodes ``f`` , only the real part of the corresponding eigenvector is taken.

.. code-block:: c++

    // Visu
    std::vector<figure> fig(5);
    for (int f = 0; f < fig.size(); f++)
    {
        matrix<double> Z = reshape(real(eval(V(row(V), f))), size(X, 1), size(X, 2));
        mesh(fig[f], X, Y, Z);
    }


See :ref:`label-reshape` , :ref:`label-mesh` . 

Code
----

Here you have all the code at once :

.. code-block:: c++

    #include "castor/matrix.hpp"
    #include "castor/smatrix.hpp"
    #include "castor/linalg.hpp"
    #include "castor/graphics.hpp"

    using namespace castor;

    int main(int argc, char const *argv[])
    {
        // Parameters
        matrix<> L = {1, 2}; // Dimensions
        double dx = 0.05;    // Space discretization

        // Discretization
        matrix<> X, Y;
        std::tie(X, Y) = meshgrid(colon(0, dx, L(0)), colon(0, dx, L(1)));

        // Visu mesh
        figure fig1;
        mesh(fig1, X, Y, zeros(size(X)));

        // Laplacian
        long nx = size(X, 2), ny = size(X, 1);
        matrix<> e = ones(nx * ny, 5);
        e(row(e), 2) = -4;
        matrix<> K = full(spdiags(1. / (dx * dx) * e, {-nx, -1, 0, 1, nx}, nx * ny, nx * ny));

        // Penalization on boundary (Homogeneous Dirichlet condition)
        matrix<std::size_t> Ibnd;
        Ibnd = find((X == 0) || (X == L(0)) || (Y == 0) || (Y == L(1)));
        K(sub2ind(size(K), Ibnd, Ibnd)) = 1e6;

        // Analytical
        auto Dth = zeros(nx, ny);
        for (int m = 0; m < nx; m++)
        {
            for (int n = 0; n < ny; n++)
            {
                Dth(m, n) = M_PI * sqrt(pow((m + 1) / L(0), 2) + pow((n + 1) / L(1), 2));
            }
        }

        // Numerical eigen values and vectors
        matrix<std::complex<double>> D, V;
        std::tie(D, V) = eig(-K, "right");

        // Sort
        matrix<std::size_t> I;
        I = argsort(abs(real(D)));
        D = eval(D(I));
        V = eval(V(row(V), I));
        matrix<std::size_t> Ith;
        Ith = argsort(Dth);
        Dth = eval(Dth(Ith));

        // Visu
        std::vector<figure> fig(5);
        for (int f = 0; f < fig.size(); f++)
        {
            matrix<double> Z = reshape(real(eval(V(row(V), f))), size(X, 1), size(X, 2));
            mesh(fig[f], X, Y, Z);
        }

        // Results
        std::cout << "-- Numerical eigenvalues --" << endl;
        disp(sqrt(real(eval(D(range(0, fig.size()))))), 1, fig.size());
        std::cout << "-- Analytical eigenvalues --" << endl;
        disp(eval(Dth(range(0, fig.size()))), 1, fig.size());
        std::cout << "-- Relative errors --" << endl;
        auto errRelative = abs((sqrt(real(eval(D(range(0, fig.size()))))) - eval(Dth(range(0, fig.size())))) / eval(Dth(range(0, fig.size())))) * 100;
        disp(errRelative, 1, fig.size());

        drawnow(fig1);

        return 0;
    }


With this code you should get these outputs :

.. code-block:: text

    -- Numerical eigenvalues --
    Matrix 1x5 of type 'd' (40 B):
        3.50946      4.43845      5.65288      6.45167      7.00046  
    -- Analytical eigenvalues --
    Matrix 1x5 of type 'd' (40 B):
        3.51241      4.44288      5.66359      6.47656      7.02481  
    -- Relative errors --
    Matrix 1x5 of type 'd' (40 B):
        0.08379      0.09980      0.18906      0.38427      0.34671 



.. figure:: img/results5eigenmodes.png
    :width: 1200
    :align: center
    :figclass: align-center
    
    From up left corner to bottom right corner : the meshgrid and the five first eigenmodes.




References
----------

http://ramanujan.math.trinity.edu/rdaileda/teach/s14/m3357/lectures/lecture_3_4_slides.pdf

http://www.cmap.polytechnique.fr/~jingrebeccali/frenchvietnammaster2_files/2017/LectureNotes/pde3d_mit.pdfClasses view and cview
++++++++++++++++++++++

Classes ``view`` and ``cview`` are used to extract a submatrix from matrix instance or to assign value to submatrix of matrix instance.

.. _label-view:

view
----

.. doxygenclass:: castor::view
   :project: castor
   :members:
   :undoc-members:

.. _label-cview:

cview
-----

.. doxygenclass:: castor::cview
   :project: castor
   :members:
   :undoc-members:
.. role:: red

.. _label-graphics-advanced:

Plotting and much more
======================

We describe here some aspects of the **graphics** library available within the **castor** project. It is based on the well-known `VTK library <https://www.vtk.org>`_ which must be installed by the user. Depending on the operating system, binary files *may* be directly available. However, one should pay attention to the fact that it is only compatible with version 9.x. 

After installing VTK (see :ref:`label-installation`), the user only needs to include the ``castor/graphics.hpp`` header file.

The path is very close to the one chosen by `Matlab <https://www.mathworks.com/products/matlab.html>`_ or `Matplotlib <https://matplotlib.org/>`_ (for the Python language). First, the user creates a ``figure`` which is a container which will hold the plots. Then,he adds one or multiple plots to the ``figure``. Finally, he asks the system to display the figures (one window per ``figure`` will be created) with the ``drawnow(...)`` function. Each call to ``drawnow`` displays **all** the ``figure`` objects which have been defined since the last call.


.. warning::
    Each call to ``drawnow`` is **blocking**, meaning that it will suspend the execution of the program until all the opened ``figure`` are closed. 
    This behavior cannot be modified as it is due to VTK.

    On Windows and linux, when the opened ``figure`` are closed, the program stops. A workaround consists in to call the `drawnow`` function at the end 
    of the program.

The **graphics** library features 2D/3D customizable plotting and basic triangular/tetrahedral mesh generation. We detail the use of some of these features below. The demo files may be found in the ``demo/demo_graphics/`` subfolder.

It is possible to **interact with the figures** in multiple ways:

 - "right-click and hold" : *rotate the content* in the direction of the mouse-pointer for **3D** figures. For **2D** plots, the user can use it to drag the plots in the direction of choice.

 - "middle-click and hold": *translate the content* in the direction of the mouse-pointer.

 - "left-click and hold": *zoom-in* if the pointer is in the upper window, otherwise *zoom-out*.

 - "key-pressed **e** or **q**": close window. **Warning:** in some systems, it closes *all* windows.

 - "key-pressed **r**": resets the view.

 - "key-pressed **w**": switches to *wireframe* view when applicable.

 - "key-pressed **s**": switches to *surface* view when applicable.

We advise the user to experiment with the different features mentioned above. The summary of all the available functions is available at :ref:`label-basic-plot`, :ref:`label-graphical-io`, :ref:`label-graphical-tools`, :ref:`label-mesh-management` and :ref:`label-mesh-plot`.


Basic 2D plotting
-----------------

In this example, we will plot a *sine* and a *square root* function with different styles but on the same display. First, we initialize some data and we create a ``figure`` object. The full code is available at ``demo/demo_graphics/basic.cpp``.

.. code:: c++

    matrix<> X = linspace(0,10,100);
    matrix<> Y = cos(X);
    figure fig;

Then, we add our first plot as a red solid line with diamond markers with the correct legend.

.. code:: c++

    plot(fig,X,Y,{"r-d"},{"sin(x)"});


Now we create a second set of data, reusing the already-defined variables. This is possible because the inputs are copied. This time, we plot it as a simple dotted blue line.

.. code:: c++

    X = linspace(-2,2,30);
    Y = sqrt(X);
    plot(fig,X,Y,{"b"},{"sqrt(x)"});
    // display
    drawnow(fig);

The figure should look like this:

.. image:: img/basic.png
    :width: 400
    :align: center

**Remark:** Currently, it is not possible to save the output to a graphic file.

Surface plot
------------

Surface plotting consists in plotting a surface defined by the equation ``Z = f(X,Y)``. First we create the *grid* (X,Y). The full code is available at ``demo/demo_graphics/surface_plot.cpp``.

.. code:: c++

    matrix<> X,Y;
    std::tie(X,Y) = meshgrid(linspace(-M_PI,M_PI,100));

Then, we create the surface which we want to plot, create a ``figure``, adjust the color axis and display everything.

.. code:: c++

    matrix<> Z = 2*sin(X)/X * sin(Y)/Y;
    // create the figure
    figure fig;
    caxis(fig,{-1,1}); // scaled color in the range [-1,1]
    mesh(fig,X,Y,Z);
    // display
    drawnow(fig);

The result is a 3-dimensional plot which should look like this : 

.. image:: img/surface_plot_wireframe.png
    :width: 400
    :align: center

**TIP:** It is possible to switch to a full *surface* plot by pressing the **s** key and switch back to a *wireframe* display by pressing the **w** key.

.. image:: img/surface_plot_surface.png
    :width: 400
    :align: center

Displaying nodal values
-----------------------

This feature is particularly useful if, for example, one needs to display the result of a finite element computation where the data is known at the vertices. In the following example, we create a planar mesh with triangular elements. Then we define a linearly increasing nodal data along the *x*-axis. The full code is available at ``demo/demo_graphics/nodal_values.cpp``.

.. code:: c++

    // geometric data
    matrix<> X,Y;
    std::tie(X,Y) = meshgrid(linspace(-1,1,10),linspace(-1,1,5));
    X = vertcat(X,X); 
    Y = vertcat(Y,Y);
    matrix<> Z = zeros(size(X));

    // create mesh
    matrix<> vtx;
    matrix<std::size_t> elt;
    std::tie(elt,vtx) = tridelaunay(X,Y,Z);

    // display
    figure fig;
    trimesh(fig,elt,vtx,eval(vtx(row(vtx),0)));
    drawnow(fig);

The result should look like this:

.. image:: img/nodal_values.png
    :width: 400
    :align: center


From mesh generation to file output
-----------------------------------

In this example, we create a spherical mesh and compute the normals to the triangles. We plot both on the same figure. Finally, we save the mesh to a *.ply* file. We also illustrate the use of the :ref:`label-quiver` function. The full code is available at ``demo/demo_graphics/advanced_mesh.cpp``.

First, we create the mesh using the :ref:`label-sphere2` function which creates a Fibonacci sphere.

.. code:: c++

    std::size_t nvtx=100;
    // 1. Create the mesh
    matrix<> vtx;
    matrix<std::size_t> elt;
    
    matrix<> X,Y,Z;
    std::tie(X,Y,Z) = sphere2(nvtx); // Fibonacci sphere
    // Delaunay tetrahedrisation
    std::tie(elt,vtx) = tetdelaunay(X,Y,Z);
    // Boundary extraction
    std::tie(elt,vtx) = tetboundary(elt,vtx);

Then, we compute the center ``ctr`` of the triangles and the normal vector ``nrm``. We recall that the coordinates of ``ctr`` are simply the averaged sum of the coordinates of the vertices of the triangles. The coordinates of ``nrm`` may be determined by computing the cross-product between the tangent of, two consecutive edges. In this example we normalize the length to 0.25 to get an equilibrated display.

.. code:: c++

    // 2. Compute the normal vectors and the centers of the triangles
    std::size_t nelt = size(elt,1);
    matrix<> nrm = zeros(nelt,3);
    matrix<> ctr = zeros(nelt,3);
    for(std::size_t ie=0; ie<nelt; ++ie)
    {
        // center
        for(std::size_t i=0; i<3; ++i)
        {
            ctr(ie,i) = (vtx(elt(ie,0),i)+vtx(elt(ie,1),i)+vtx(elt(ie,2),i))/3.;
        }
        // normal vector to triangle {A,B,C}
        matrix<> AB=zeros(1,3), BC=zeros(1,3), nv = zeros(1,3);
        for(std::size_t i=0; i<3; ++i)
        {
            // tangent to first edge
            AB(i) = vtx(elt(ie,1),i) - vtx(elt(ie,0),i);
            tangent to second edge
            BC(i) = vtx(elt(ie,2),i) - vtx(elt(ie,1),i);
        }
        // for single vectors, this code is faster than
        // a call to 'cross' (for the cross-product) or 'norm'.
        nv(0) = AB(1)*BC(2) - AB(2)*BC(1);
        nv(1) = AB(2)*BC(0) - AB(0)*BC(2);
        nv(2) = AB(0)*BC(1) - AB(1)*BC(0);
        double l = std::sqrt(nv(0)*nv(0)+nv(1)*nv(1)+nv(2)*nv(2));
        // normalization
        for(std::size_t i=0; i<3; ++i) 
        {
            nrm(ie,i) = nv(i)/(2*l); // arrows of length 0.5
        }
    }

Now, we plot the result.

.. code:: c++

    // 3. Plot everything
    figure fig;
    trimesh(fig,elt,vtx); // plot the mesh
    quiver(fig,ctr,nrm);  // plot the normal vectors at the centers
    drawnow(fig);

Finally, save the mesh into the current directory in the *.ply* format.

.. code:: c++

    // 4. save to .ply
    std::string path="./", name="testfile.ply";
    triwrite(path,name,elt,vtx);

The result should look like an urchin, see below. **Please note that the normal vectors may not appear when the window appears. In that case, simply clicking inside the window should do the trick.**

.. image:: img/advanced_mesh.png
    :width: 400
    :align: center
.. _label-algorithms:

Algorithms
++++++++++

Standard algorithms over data stored in matrices may be found here. Among *many* others, it is possible to :ref:`label-sort` a ``matrix``, find the :ref:`label-unique` elements, the :ref:`label-median` of the data, perform the :ref:`label-dot` product between two matrices, find the maximum value with :ref:`label-max`, etc.


.. _label-argintersect:

argintersect
------------
.. doxygenfunction:: argintersect
   :project: castor

See :ref:`label-intersect`, :ref:`label-argsetdiff`, :ref:`label-argunique`.

.. _label-argmax:

argmax
------
.. doxygenfunction:: argmax(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: argmax(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-max`, :ref:`label-maximum`, :ref:`label-argmin`, :ref:`label-argsort`.

.. _label-argmin:

argmin
------
.. doxygenfunction:: argmin(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: argmin(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-min`, :ref:`label-minimum`, :ref:`label-argmax`, :ref:`label-argsort`.

.. _label-argsetdiff:

argsetdiff
----------
.. doxygenfunction:: argsetdiff
   :project: castor

See :ref:`label-setdiff`, :ref:`label-argintersect`, :ref:`label-argunique`.

.. _label-argsort:

argsort
-------
.. doxygenfunction:: argsort
   :project: castor

See :ref:`label-sort`, :ref:`label-argmin`, :ref:`label-argmax`.

.. _label-argunique:

argunique
---------
.. doxygenfunction:: argunique
   :project: castor

See :ref:`label-unique`, :ref:`label-argsetdiff`, :ref:`label-argintersect`.

.. _label-conv:

conv
----
.. doxygenfunction:: conv
   :project: castor

See :ref:`label-fftconv`.

.. _label-cross:

cross
-----
.. doxygenfunction:: cross
   :project: castor

See :ref:`label-dot`.

.. _label-cumprod:

cumprod
-------
.. doxygenfunction:: cumprod(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: cumprod(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-cumsum`, :ref:`label-prod`.

.. _label-cumsum:

cumsum
------
.. doxygenfunction:: cumsum(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: cumsum(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-cumprod`, :ref:`label-sum`.

.. _label-diff:

diff
----
.. doxygenfunction:: diff(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: diff(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-sum`, :ref:`label-prod`.

.. _label-dot:

dot
---
.. doxygenfunction:: dot
   :project: castor

See :ref:`label-cross`.

.. _label-fftconv:

fftconv
-------
.. doxygenfunction:: fftconv
   :project: castor

See :ref:`label-conv`, :ref:`label-fft`.


.. _label-gmres:

gmres
-----
.. doxygenfunction::  gmres(matrix<T> const &A, matrix<T> const &B, double tol = 1e-6, std::size_t maxit = 10, std::function<matrix<T>(matrix<T> const&)> const &Am1 = std::function<matrix<T>(matrix<T> const&)>(), matrix<T> const &X0 = matrix<T>())
   :project: castor
.. doxygenfunction::  gmres(matrix<T> const &A, matrix<T> const &B, double tol, std::size_t maxit, matrix<T> const &Am1, matrix<T> const &X0 = matrix<T>())
   :project: castor
.. doxygenfunction:: gmres(std::function<matrix<T>(matrix<T> const&)> const &A, matrix<T> const &B, double tol = 1e-6, std::size_t maxit = 10, std::function<matrix<T>(matrix<T> const&)> const &Am1 = std::function<matrix<T>(matrix<T> const&)>(), matrix<T> const &X0 = matrix<T>())
   :project: castor



See :ref:`label-linsolve`.

.. _label-intersect:

intersect
---------
.. doxygenfunction:: intersect
   :project: castor

See :ref:`label-argintersect`, :ref:`label-setdiff`, :ref:`label-unique`, :ref:`label-union2`.

.. _label-kron:

kron
----
.. doxygenfunction:: kron(R A, matrix<S> const &B)
   :project: castor
.. doxygenfunction:: kron(matrix<R> const &A, S B)
   :project: castor
.. doxygenfunction:: kron(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-mtimes`.

.. _label-max:

max
---
.. doxygenfunction:: max(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: max(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-argmax`, :ref:`label-maximum`, :ref:`label-min`, :ref:`label-sort`.

.. _label-min:

min
---
.. doxygenfunction:: min(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: min(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-argmin`, :ref:`label-minimum`, :ref:`label-max`, :ref:`label-sort`.

.. _label-maximum:

maximum
-------
.. doxygenfunction:: maximum(R A, matrix<S> const &B)
   :project: castor
.. doxygenfunction:: maximum(matrix<R> const &A, S B)
   :project: castor
.. doxygenfunction:: maximum(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-minimum`, :ref:`label-max`, :ref:`label-min`, :ref:`label-sort`.

.. _label-mean:

mean
----
.. doxygenfunction:: mean(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: mean(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-max`, :ref:`label-min`, :ref:`label-median`, :ref:`label-variance`, :ref:`label-stddev`.

.. _label-median:

median
------
.. doxygenfunction:: median(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: median(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-max`, :ref:`label-min`, :ref:`label-variance`, :ref:`label-stddev`.

.. _label-minimum:

minimum
-------
.. doxygenfunction:: minimum(R A, matrix<S> const &B)
   :project: castor
.. doxygenfunction:: minimum(matrix<R> const &A, S B)
   :project: castor

.. doxygenfunction:: minimum(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-maximum`, :ref:`label-min`, :ref:`label-max` :ref:`label-sort`.

.. _label-mtimes:

mtimes
------
.. doxygenfunction:: mtimes(R A, matrix<S> const &B)
   :project: castor
.. doxygenfunction::  mtimes(matrix<R> const &A, S B)
   :project: castor
.. doxygenfunction:: mtimes(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-tgemm-naive`, :ref:`label-kron`.

.. _label-norm:

norm
----
.. doxygenfunction:: norm(matrix<S> const &A, std::string typ = "2")
   :project: castor
.. doxygenfunction:: norm(matrix<S> const &A, std::string typ, int dim)
   :project: castor

See :ref:`label-max`, :ref:`label-sum`.

.. _label-prod:

prod
----
.. doxygenfunction:: prod(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: prod(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-sum`, :ref:`label-diff`.

.. _label-setdiff:

setdiff
-------
.. doxygenfunction:: setdiff
   :project: castor

See :ref:`label-argsetdiff`, :ref:`label-intersect`, :ref:`label-unique`, :ref:`label-union2`.

.. _label-sort:

sort
----
.. doxygenfunction:: sort
   :project: castor

See :ref:`label-argsort`, :ref:`label-min`, :ref:`label-max`.

.. _label-stddev:

stddev
------
.. doxygenfunction:: stddev(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: stddev(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-max`, :ref:`label-min`, :ref:`label-mean`, :ref:`label-median`, :ref:`label-variance`.

.. _label-sum:

sum
---
.. doxygenfunction:: sum(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: sum(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-prod`, :ref:`label-diff`, :ref:`label-cumsum`.


.. _label-tgemm-naive:

tgemm
-----
.. doxygenfunction:: tgemm(P alpha, matrix<Q> const &A, matrix<R> const &B, S beta, matrix<T> &C)
   :project: castor


See :ref:`label-mtimes`, :ref:`label-tgemm-blas`.

.. _label-union2:

union2
------
.. doxygenfunction:: union2
   :project: castor

See :ref:`label-intersect`, :ref:`label-setdiff`, :ref:`label-unique`.

.. _label-unique:

unique
------
.. doxygenfunction:: unique
   :project: castor

See :ref:`label-argunique`, :ref:`label-intersect`, :ref:`label-setdiff`, :ref:`label-union2`.

.. _label-variance:

variance
--------
.. doxygenfunction:: variance(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: variance(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-max`, :ref:`label-min`, :ref:`label-median`, :ref:`label-mean`, :ref:`label-stddev`.
.. _label-matrix-manipulations:

Matrix manipulations
++++++++++++++++++++

Many functions are provided for the manipulation of matrices. For instance, it is possible to :ref:`label-cast` a ``matrix`` object into a different type, to concatenate matrices in all dimensions (see :ref:`label-cat`, etc.), :ref:`label-find` the non-zero elements, or to :ref:`label-transpose` it.

.. _label-all:

all
---
.. doxygenfunction:: all(matrix<T> const &A)
   :project: castor

See :ref:`label-row`, :ref:`label-col`, :ref:`label-get`, :ref:`label-set`, :ref:`label-view`, :ref:`label-cview`.


.. _label-cast:

cast
----
.. doxygenfunction:: cast
   :project: castor

See :ref:`label-class-matrix`.


.. _label-cat:

cat
---
.. doxygenfunction:: cat(int dim, R const &A, S const &B)
   :project: castor
.. doxygenfunction:: cat(int dim, R A, matrix<S> const &B)
   :project: castor
.. doxygenfunction:: cat(int dim, matrix<R> const &A, S B)
   :project: castor
.. doxygenfunction:: cat(int dim, matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-vertcat`, :ref:`label-horzcat`.


.. _label-clear:

clear
-----
.. doxygenfunction:: clear(matrix<T> &A)
   :project: castor

See :ref:`label-class-matrix`, :ref:`label-resize`.


.. _label-col:

col
---
.. doxygenfunction:: col(matrix<T> const &A)
   :project: castor

See :ref:`label-row`, :ref:`label-all`, :ref:`label-get`, :ref:`label-set`, :ref:`label-view`, :ref:`label-cview`.


.. _label-find:

find
----
.. doxygenfunction:: find(matrix<T> const &A)
   :project: castor

See :ref:`label-nnz`, :ref:`label-ind2sub`.


.. _label-get:

get
---
.. doxygenfunction:: get(matrix<T> const &A, matrix<std::size_t> const &I, matrix<std::size_t> const &J)
   :project: castor
.. doxygenfunction:: get(matrix<T> const &A, matrix<std::size_t> const &L)
   :project: castor

See :ref:`label-set`, :ref:`label-all`, :ref:`label-row`, :ref:`label-col`, :ref:`label-view`, :ref:`label-cview`.


.. _label-horzcat:

horzcat
-------
.. doxygenfunction:: horzcat
   :project: castor

See :ref:`label-vertcat`, :ref:`label-cat`.


.. _label-ind2sub:

ind2sub
-------
.. doxygenfunction:: ind2sub
   :project: castor

See :ref:`label-sub2ind`, :ref:`label-find`.


.. _label-isempty:

isempty
-------
.. doxygenfunction:: isempty
   :project: castor

See :ref:`label-isequal`, :ref:`label-isvector`.


.. _label-isequal:

isequal
-------
.. doxygenfunction:: isequal
   :project: castor

See :ref:`label-isempty`, :ref:`label-isvector`.


.. _label-isvector:

isvector
--------
.. doxygenfunction:: isvector
   :project: castor

See :ref:`label-isequal`, :ref:`label-isempty`.

.. _label-range:

range
-----
.. doxygenfunction:: range
   :project: castor

See :ref:`label-colon`, :ref:`label-get`, :ref:`label-set`, :ref:`label-view`.

.. _label-resize:

resize
------
.. doxygenfunction:: resize(std::size_t m, std::size_t n, T v = (T)NAN)
   :project: castor

See :ref:`label-reshape`.


.. _label-reshape:

reshape
-------
.. doxygenfunction:: reshape(std::size_t m, std::size_t n)
   :project: castor

See :ref:`label-resize`, :ref:`label-transpose`.


.. _label-row:

row
---
.. doxygenfunction:: row(matrix<T> const &A)
   :project: castor

See :ref:`label-all`, :ref:`label-col`, :ref:`label-get`, :ref:`label-set`, :ref:`label-view`, :ref:`label-cview`.


.. _label-set:

set
---
.. doxygenfunction:: set(matrix<T> &A, matrix<std::size_t> const &I, matrix<std::size_t> const &J, U b)
   :project: castor
.. doxygenfunction:: set(matrix<T> &A, matrix<std::size_t> const &L, matrix<U> const &B)
   :project: castor
.. doxygenfunction:: set(matrix<T> &A, matrix<std::size_t> const &I, matrix<std::size_t> const &J, matrix<U> const &B)
   :project: castor
.. doxygenfunction:: set(matrix<T> &A, matrix<std::size_t> const &L, U b)
   :project: castor

See :ref:`label-get`, :ref:`label-all`, :ref:`label-row`, :ref:`label-col`, ::ref:`label-view`, :ref:`label-cview`.


.. _label-sub2ind:

sub2ind
-------
.. doxygenfunction:: sub2ind
   :project: castor

See :ref:`label-ind2sub`, :ref:`label-find`.


.. _label-transpose:

transpose
---------
.. doxygenfunction:: transpose(matrix<T> const &A)
   :project: castor

See :ref:`label-reshape`.


.. _label-vertcat:

vertcat
-------
.. doxygenfunction:: vertcat
   :project: castor

See :ref:`label-horzcat`, :ref:`label-cat`.
