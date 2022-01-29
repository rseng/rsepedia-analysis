---
title: 'fastSF: A parallel code for computing the structure functions of turbulence'

tags:
  - C++
  - structure functions
  - turbulence
  - fluid dynamics

authors:
  - name: Shubhadeep Sadhukhan
    orcid: 0000-0002-7278-5041
    affiliation: 1
  - name: Shashwat Bhattacharya
    orcid: 0000-0001-7462-7680
    affiliation: 2
  - name: Mahendra K. Verma
    orcid: 0000-0002-3380-4561
    affiliation: 1
  

affiliations:
 - name: Department of Physics, Indian Institute of Technology Kanpur, Kanpur 208016, India
   index: 1
 - name: Department of Mechanical Engineering, Indian Institute of Technology Kanpur 208016, India
   index: 2

date: 18 February 2020

bibliography: paper.bib

---

# Summary

Turbulence is a complex phenomenon in fluid dynamics involving nonlinear interactions between multiple scales. Structure functions are popular diagnostics in the study of statistical properties properties of turbulent flows [@Kolmogorov:Dissipation; @Kolmogorov:Structure; @Frisch:book]. Some of the earlier works comprising of such analysis are those of @Gotoh:PF2002, @Kaneda:PF2003, and @Ishihara:ARFM2009 for three-dimensional (3D) hydrodynamic turbulence; @Yeung:PF2005 and @Ray:NJP2008 for passive scalar turbulence; @Biferale:NJP2004 for two-dimensional (2D) hydrodynamic turbulence; and @Kunnen:PRE2008, @Kaczorowski:JFM2013, and @Bhattacharya:PF2019 for turbulent thermal convection. Structure functions are two-point statistical quantities; thus, an accurate computation of these quantities requires averaging over many points. However, incorporation of a large number of points makes the computations very expensive and challenging. Therefore, we require an efficient parallel code for accurate computation of structure functions. In this paper, we describe the design and validation of the results of ``fastSF``, a parallel code to compute the structure functions for a given velocity or scalar field. 

``fastSF`` is a C++ application for computing the structure functions of scalar and vector fields on Cartesian grids of a 2D or 3D periodic box, stored as HDF5 files. The code employs MPI (Message Passing Interface) parallelization with equal load distribution and vectorization for efficiency on SIMD architectures. The user can select the range of the orders of the structure functions to be computed and the computed structure functions are written to HDF5 files that can be further processed by the user.

We are not aware of any other open soure or commercial packages for computing structure functions; prior studies have relied on in-house software that was never publicly released.  As an open source package, `fastSF` provides a standard high-performance implementation and thus facilitates wider use of structure functions.

``fastSF`` uses MPI [@Pacheco:book:PP] for parallelism, HDF5 [@HDF5_web] via H5SI [@H5SI_web] for reading gridded field data and writing structure functions, as well as blitz++ [@Blitz_web] for vectorized computation and yaml-cpp [@YAML_web] for reading control parameters.
 In the next section, we will briefly explain the velocity and scalar structure functions in turbulent flows.



# Velocity and scalar structure functions

We denote the velocity and scalar fields using $\boldsymbol{u}$ and $\theta$  respectively. The velocity difference between any two points $\boldsymbol{r}$ and $\boldsymbol{r}+\boldsymbol{l}$ is $\delta \boldsymbol{u} = \boldsymbol{u(r}+ \boldsymbol{l)}-\boldsymbol{u(r)}$. The difference in the parallel components of the velocity field along $\boldsymbol{l}$ is $\delta u_\parallel=\delta \boldsymbol{u}\cdot \hat{\boldsymbol{l}}$.  The corresponding difference in the perpendicular component is $\delta u_\perp= |\delta \boldsymbol{u} - \delta u_\parallel \hat{\boldsymbol{l}}|$. Assuming statistical homogeneity, we define the longitudinal velocity structure functions of order $q$ as
$$ S_q^{u_\parallel}(\boldsymbol{l}) = \langle (\delta u_\parallel)^q \rangle = \langle [\{\boldsymbol{u(r+l)}-\boldsymbol{u(r)}\}\cdot \hat{\boldsymbol{l}}]^q \rangle, \quad \quad (1)$$ 
and the transverse velocity structure functions of order 
$q$ as 
$$ S_q^{u_\perp}(\boldsymbol{l}) = \langle (\delta u_\perp)^q \rangle = \langle |\delta \boldsymbol{u} - \delta u_\parallel \hat{\boldsymbol{l}}|^q \rangle. \quad \quad (2)$$ 
Here, $\langle \cdot \rangle$ denotes spatial averaging. Similarly, we can define the scalar structure functions for the scalar field as 
$$ S_q^\theta(\boldsymbol{l}) = \langle (\delta \theta)^q\rangle = \langle [\theta (\boldsymbol{r+l}) - \theta(\boldsymbol{r})]^q \rangle. \quad \quad(3)$$

For isotropic turbulence (in addition to being homogeneous), the structure functions become functions of $l$, where $l=|\boldsymbol{l}|$. The second-order velocity structure function $S_q^{u_{\parallel}}(l)$ provides an estimate for the energy in the eddies of size $l$ or less [@Davidson:book:Turbulence]. 

![For 3D homogeneous isotropic turbulence: plots of the negative of normalized third, fifth and seventh-order longitudinal velocity structure functions vs. $l$. The negative of the normalized third-order structure function is close to $4/5$ (dashed line) in the inertial range. \label{SF_Hydro}](docs/figs/SF_hydro.png)

For 3D incompressible hydrodynamic turbulence with homegeneity and isotropy, the third-order longitudinal velocity structure function in the inertial range (scales lying between the large-scale forcing regime and the small-scale dissipation regime) is given by [@Kolmogorov:Dissipation; @Kolmogorov:Structure; @Frisch:book]
$$S_3^{u_\parallel}(l) = -\frac{4}{5} \epsilon l \sim -l, \quad \quad (4)$$
where $\epsilon$ is the viscous dissipation rate. 
For an arbitrary order $q$, @She:PRL1994 proposed that the longitudinal structure functions scale as $S_q^{u_\parallel} (l) \sim \zeta_q$, where the exponent $\zeta_q$ is given by
$$\zeta_q = \frac{q}{9} + 2\left(1 - \left( \frac{2}{3} \right)^{q/3} \right). \quad \quad (5)$$


Figure~\ref{SF_Hydro} exhibits the plots of the negative of the normalized 3rd, 5th, and 7th-order longitudinal velocity structure functions computed using the simulation data of 3D hydrodynamic turbulence [@Sadhukhan:PRF2019]. The structure functions are normalized by $(\epsilon l)^{\zeta_q}$, where $\zeta_q$ is given by Eq. (5).  In the inertial range (0.2 < l < 0.7), the normalized third-order longitudinal velocity structure function is fairly close to $4/5$ (represented by dashed line), consistent with Kolmogorov's theory. Moreover, the normalized fifth and seventh-order structure functions show a plateau for the same range of l, thus exhibiting consistency with She-Leveque's model.




In the next section, we provide a brief description of the code.

# Design of the Code
In this section, we present a sketch of the structure function computation for the velocity structure functions.  We employ vectorization and loops over $\boldsymbol{l}$, thus requiring three loops for 3D fields and two loops for 2D fields. In the following, we provide the algorithm for structure function computation for a 2D velocity field.

**Pseudo-code**

*Data*: Velocity field $\boldsymbol{u}$ in domain $(L_x, L_z)$; number of processors $P$.

*Procedure*:
 
* Divide $\boldsymbol{l}$'s among various processors. The process of data division among the processors  has been described later in this section. 
 
* For every processor:
     
    * for $\boldsymbol{l}= (l_x,l_z)$ assigned to the processor:
        
        * Compute $\delta \boldsymbol{u}(l_x,l_z)$ by taking the difference between two points with the same indices in pink and green subdomains as shown in Fig. \ref{Schematic}. This feature enables vectorized subtraction operation.
        
        * $\delta u_{\parallel}(l_x,l_z) = \delta \boldsymbol{u} \cdot \hat{\boldsymbol{l}}$ (Vectorized). 
        
        * $\delta u_{\perp}(l_x,l_z) = |\delta \boldsymbol{u} - \delta u_{\parallel} \hat{\boldsymbol{l}}$| (Vectorized). 
        
        * for order $q$:
        
            * $S_q^{u_{\parallel}}(l_x,l_z) =$ Average of $\delta u_{\parallel}^q$ (Vectorized).
            
            * $S_q^{u_{\perp}}(l_x,l_z) =$ Average of $\delta u_{\perp}^q$ (Vectorized).
            
            * Send the values of $S_q^{u_{\parallel}}(l_x,l_z)$, $S_q^{u_{\perp}}(l_x,l_z)$, $q$, $l_x$, and $l_z$ to the root process.
            
* The root process stores $S_q^{u_{\parallel}} (l_x, l_z)$ and $S_q^{u_{\perp}} (l_x, l_z)$.
            
* Stop

![The velocity difference $\delta \boldsymbol{u}(\boldsymbol{l})$ is computed by taking the difference between two points with the same indices in the pink and the green subdomains. For example, $\boldsymbol{u}(\boldsymbol{l}) - \boldsymbol{u}(0,0) = \boldsymbol{u}_B - \boldsymbol{u}_A$, where $B$ and $A$ are the origins of the green and the pink subdomains. This feature enables vecotrization of the computation. \label{Schematic}](docs/figs/Schematic.png)

Since $S_q^u(\boldsymbol{l})$ is important for intermediate scales (inertial range) only, we vary $\boldsymbol{l}$ upto half the domain size, that is, upto ($L_x/2, L_z/2$), to save computational cost. The $\boldsymbol{l}$'s are divided among MPI processors along $x$ and $z$ directions. Each MPI processor computes the structure functions for the points assigned to it and has access to the entire input data. 
After computing the structure function for a given $\boldsymbol{l}$, each processor communicates the result to the root process, which stores the $S_q^{u_\parallel}(\boldsymbol{l})$ and $S_q^{u_\perp}(\boldsymbol{l})$ arrays.

It is clear from Fig. \ref{Schematic} that the sizes of the pink or green subdomains are $(L_x-l_x)(L_z-l_z)$, which are function of $\boldsymbol{l}$'s.  This function decreases with increasing $\boldsymbol{l}$ leading to larger computational costs for small $l$ and less cost of larger $l$.   Hence, a straightforward division of the domain among the processors along $x$ and $z$ directions will lead to a load imbalance.   Therefore, we assign both large and small $\boldsymbol{l}$'s to each processor to achieve equal load distribution. We illustrate the above idea  using the following example.

Consider a one-dimensional domain of size $L=15$, for which the possible $l$'s are
$$l=\{0, 1, 2, 3 ... 15\}.$$ 
We need to compute the structure functions for $l$ ranging from 0 to 7. We divide the task among four processors, with two $l$'s assigned to each processor. The following distribution of $l$'s ensures equal load distribution:
$$\mbox{Processor 0: } \quad l=\{0,7\}, \quad \sum(L-l)=(15-0)+(15-7) = 23,$$
$$\mbox{Processor 1: } \quad l=\{1, 6\}, \quad \sum(L-l)=(15-1)+(15-6) = 23,$$
$$\mbox{Processor 2: } \quad l=\{2,5\}, \quad \sum(L-l)=(15-2)+(15-5) = 23,$$
$$\mbox{Processor 3: } \quad l=\{3, 4\}, \quad \sum(L-l)=(15-3)+(15-4) = 23.$$
Similarly, if two processors are used, then the following distribution results in load balance. 
$$\mbox{Processor 0: } \quad l=\{0, 7, 2, 5\},$$
$$\mbox{Processor 1: } \quad l=\{1, 6, 3, 4\}.$$
 This idea of load distribution has been implemented in our program and has been extended to higher dimensions. 

Note that for 2D, $l_x>0$, but $l_z$ can take both positive and negative values. However, for isotropic turbulence, the structure functions for $+l_z$ and $-l_z$ are statistically equal. Therefore, in our computations, we restrict  to $l_x>0$, $l_z>0$. For anisotropic turbulence, not discussed here, the structure functions will depend on $(l_x,l_z)$ rather than $l$; our code will be extended to such systems in future. 

For 3D turbulence, the structure functions will depend on $(l_x,l_y,l_z)$. We divide the tasks among processors over $l_x$ and $l_y$ as done  above for 2D turbulence. The aforementioned algorithm can be easily extended to the 3D case. We employ a similar method for the computation of scalar structure functions as well.
 
In the next section, we discuss the scaling of our code.

# Scaling of `fastSF`

`fastSF` is scalable over many processors due to vectorization and equal load distribution. We demonstrate the scaling of `fastSF` for the third-order longitudinal structure function for an idealized velocity field on a $128^3$ grid.  For our computation we employ a maximum of 1024 cores. We take the velocity field as
$$\boldsymbol{u} = 
\begin{bmatrix} 
x \\ y \\z
\end{bmatrix}.$$
We perform four runs on a Cray XC40 system (Shaheen II of KAUST) for this problem using a total of 16, 64, 256, and 1024 cores. We used 16 cores per node for each run. In Fig. \ref{Scaling}, we plot the inverse of time taken in seconds versus the number of cores. The best fit curve for these data points yields
$$T^{-1} \sim p^{0.986 \pm 0.002},$$
Thus, the data-points follow $T^{-1} \sim p$ curve to a good approximation. Hence, we conclude that our code exhibits good scaling. 

![Scaling of `fastSF` for the computation of third-order longitudinal velocity structure function using 16, 64, 256, and 1024 processors of Shaheen II. All the runs were conducted on a $128^3$ grid.  We observe a linear scaling. \label{Scaling}](docs/figs/SF_scaling.png)


# Conclusions

This paper provides a brief description of ``fastSF``, an efficient parallel C++ code that computes structure functions for given velocity and scalar fields. This code is shown to be scalable over many processors. An earlier version of the code was used by @Bhattacharya:PF2019 for analyzing the structure functions of turbulent convection.  We believe that ``fastSF`` will be useful to turbulence community as it facilitates wider use of structure functions.  


# Acknowledgements

We thank Roshan Samuel, Anando Chatterjee, Soumyadeep Chatterjee, and Manohar Sharma for helpful discussions during the development of ``fastSF``. We are grateful to Jed Brown, Ilja Honkonen, and Chris Green for a careful review of our work and their useful suggestions. Our computations were performed on Shaheen II at KAUST supercomputing laboratory, Saudi Arabia, under the project k1416. 

---

# References


# `fastSF`

`fastSF` is an open source hybrid parallel C++ code to compute structure functions for a given velocity or scalar field.

## Getting the Source Code

`fastSF` is hosted on GitHub. You can download the source code from the following link:

https://github.com/ShubhadeepSadhukhan1993/fastSF


## Installing `fastSF`
`fastSF` requires several libraries (listed later in this section) for compiling. These libraries can be installed in any directory as long as the said directory is added to the user $PATH. In this section, we will provide the guidelines to install these libraries in $HOME/local (which is generally recommended).

### Environment Variables
The following environment variables need to be set before compiling `fastSF`. (You may append these lines to `$HOME/.bashrc`):


`export PATH=$HOME/local/bin:$PATH`

`export PKG_CONFIG_DISABLE_UNINSTALLED=true`

`export PKG_CONFIG_PATH=$HOME/local/lib/pkgconfig:$PKG_CONFIG_PATH`

`export HDF5_ROOT=$HOME/local`

`export CPATH=$HOME/local/include/:$CPATH`

`export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH`

`export LIBRARY_PATH=$HOME/local/lib:$LIBRARY_PATH`

`export MANPATH=$HOME/local/share/man/:$MANPATH`
`

### Required Libraries

The following libraries are required for installing and running fastSF:
 
* [`Blitz++`](https://github.com/blitzpp/blitz) (Version 1.0.2)-
All array manipulations are performed using the `Blitz++` library. Download `Blitz++` from [here](https://github.com/blitzpp/blitz). After downloading, change to the `blitz-master` directory and enter the following commands

	`CC=gcc CXX=g++ ./configure --prefix=$HOME/local`
	
	`make install`

* [`YAML-cpp`](https://github.com/jbeder/yaml-cpp/releases/tag/release-0.3.0)(Version 0.3.0) - 
	The input parameters are stored in the `para.yaml` file which needs the `YAML-cpp` library to parse. Download `YAML-cpp` from [here](https://github.com/jbeder/yaml-cpp/releases/tag/release-0.3.0). Extract the zip/tar file and change the `yaml-cpp-release-0.3.0` directory. Important: Please ensure that [CMake](https://cmake.org) is installed in your system. Enter the following commands:
	
	`CC=gcc CXX=g++ cmake -DCMAKE_INSTALL_PREFIX=$HOME/local`
	
	`make install`

	
* An `MPI` (Message Passing Interface) Library - 
`fastSF` uses `MPI` for parallelism. The software was tested using [`MPICH`](https://www.mpich.org)(Version 3.3.2), however, any standard `MPI-1` implementation should be sufficient.  Here, we will provide instructions for installing `MPICH`. Download `MPICH` from [here](https://www.mpich.org/downloads/). After extraction, change to `mpich-3.3.2` folder and enter the following:

	`CC=gcc CXX=g++ ./configure --prefix=$HOME/local`
	
	`make install`

* [`HDF5`](http://turbulencehub.org/wp-content/uploads/Download_Files/hdf5-1.8.20.tar.bz2)(Version 1.8.20) -
The output files are written in HDF5 format. Download `HDF5` from [here](http://turbulencehub.org/wp-content/uploads/Download_Files/hdf5-1.8.20.tar.bz2). After extracting the tar file, change to `hdf5-1.8.20` and enter the following:

	`CC=mpicc CXX=mpicxx ./configure --prefix=$HOME/local --enable-parallel --without-zlib`

	`make install`

* [`H5SI`](https://github.com/anandogc/h5si)(Version 1.1.1) - 
This library is used for simplifying the input-output operations of `HDF5`. Download `H5SI` from [here](https://github.com/anandogc/h5si). After downloading the zip file, extract it and change to `h5si-master/trunk`. Important: Please ensure that [CMake](https://cmake.org) is installed in your system. Enter the following:

	`CXX=mpicxx cmake . -DCMAKE_INSTALL_PREFIX=$HOME/local`
	
	`make`
	
	`make install`


IMPORTANT: 

* Note that `fastSF` is not compatible with higher versions `YAML-cpp`(> 0.3.0). 


###  Compiling instruction

After downloading `fastSF`, change into `fastSF/src` directory and run the command `make` in the terminal. An executable named `fastSF.out` will be created inside the `fastSF/src` folder.

## Testing `fastSF`
`fastSF` offers an automated testing process to validate the code. The relevant test scripts can be found in the `tests/` folder of the code. To execute the tesing process, change into `fastSF` and run the command 

`bash runTest.sh`. 

The code then runs four test cases; these are as follows. 

* In the first case, the code will generate a 2D velocity field given by **u** = [*x, z*], and compute the structure functions for the given field. For this case, the longitudinal structure functions should equal *l<sup>q</sup>*. 

* In the second case, the code will generate a 2D scalar field given by *T = x + z*, and compute the structure functions for the given field. For this case, the structure functions should equal *(l<sub>x</sub> + l<sub>z</sub>)<sup>q</sup>*.

* In the third case, the code will generate a 3D velocity field given by **u** = [*x, y, z*], and compute the structure functions for the given field. For this case, the longitudinal structure functions should equal *l<sup>q</sup>*. 

* In the fourth case, the code will generate a 3D scalar field given by *T = x + y + z*, and compute the structure functions for the given field. For this case, the structure functions should equal *(l<sub>x</sub> + l<sub>y</sub> + l<sub>z</sub>)<sup>q</sup>*.

For the above cases, `fastSF` will compare the computed structure functions with the analytical results. If the percentage difference between the two values is less than 10<sup>-10</sup>, the code is deemed to have passed. 

Finally, for visualization purpose, the python script `test/test.py` is invoked. This script generates the plots of the second and third-order longitudinal structure functions versus *l*, and the density plots of the computed second-order scalar structure functions and *(l<sub>x</sub> + l<sub>z</sub>)<sup>2</sup>*. For the 3D scalar field, the density plots of the computed second-order scalar structure functions for *l<sub>y</sub> = 0.5* and *(l<sub>x</sub> + 0.5 + l<sub>z</sub>)<sup>2</sup>* are generated. These plots demonstrate that the structure functions are computed accurately. Note that the following python modules are needed to run the test script successfully:

1. `h5py`
2. `numpy`
3. `matplotlib`


## Detailed instruction for running `fastSF`

This section provides a detailed procedure to execute `fastSF` for a given velocity or scalar field.

### i) Files Required and HDF5 Schema:

All the files storing the input fields needs to be in the `hdf5` format and stored inside the `in` folder. All the input datasets are required to be precisely as per the following schema.

#### For two dimensional fields

For vector field, two datasets files are required:

1. A 2D dataset storing the x-component of the vector / velocity field. This dataset has dimensions (`Nx,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `U.V1r.h5` and the dataset as `U.V1r`.

2. A 2D dataset storing the z-component of the vector / velocity field. This dataset has dimensions (`Nx,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `U.V3r.h5` and the dataset as `U.V3r`. 

For scalar field, one dataset is required:

1. A 2D dataset storing the scalar field. This dataset has dimensions (`Nx,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `T.Fr.h5` and the dataset as `T.Fr`.

#### For three dimensional fields

For vector field, three `hdf5` files are required:

1. A 3D dataset storing the x-component of the vector / velocity field. This dataset has dimensions (`Nx,Ny,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `U.V1r.h5` and the dataset as `U.V1r`.

2. A 3D dataset storing the y-component of the vector / velocity field. This dataset has dimensions (`Nx,Ny,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `U.V2r.h5` and the dataset as `U.V2r`.

3. A 3D dataset storing the z-component of the vector / velocity field. This dataset has dimensions (`Nx,Ny,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `U.V3r.h5` and the dataset as `U.V3r`.

For scalar field, For vector field, one `hdf5` file is required:

1. A 3D dataset storing the scalar field. This dataset has dimensions (`Nx,Ny,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `T.Fr.h5` and the dataset as `T.Fr`.


**IMPORTANT:** For vector fields, it must be ensured that the dimensions of `Ux`, `Uy`, and `Uz` are the same, otherwise the code will throw an error.

### ii) `para.yaml` details

The user can specify the relevant parameters via command line or via parameters file. If no command-line options are given, the entries in the parameters file will be taken. First, we will explain how to use the parameters file, and then we will proceed to command-line options.

`fastSF` has a folder named `in`. This folder contains the input field files in `hdf5` format, and a parameters file named `para.yaml`. You need to provide the required input parameters in this file. The details of the entries are as follows:




#### `program: scalar_switch`

You can enter `true` or `false`

`true`: Calculate the structure function for a given scalar field. 

`false`: Calculate the structure function for a given velocity field. 

#### `program: 2D_switch`

You can enter `true` or `false`.

`true`: Calculate the structure function for two dimensional fields. 

`false`: Calculate the structure function for three dimensional fields.

#### `program: Only_logitudinal`

This entry is for structure functions for velocity  fields only. You can enter `true` or `false`.

`true`: Compute only the longitudinal structure function.

`false`: Compute both longitudinal and transverse structure functions.

#### `program: Processors_X`

The number of processors in x-direction. Only integer values are accepted. Note that this value should be an integer factor of the total number of processors.

#### `grid: Nx, Ny, Nz` (Only applicable if test switch is set to `true`, in which case the code generates the input fields)

The number of points along *x*, *y*, and *z* direction respectively of the  grid. Valid for both the vector and scalar fields. 
For two dimensional fields you need to provide `Nx` and `Nz`.

Only integer values will be accepted.

#### `domain_dimension: Lx, Ly, Lz`

Length of the cubical box along *x*, *y*, and *z* direction respectively. 
For two dimensional fields, you need to provide `Lx` and `Lz`. 


#### `structure_function: q1, q2`

The lower and the upper limit of the order of the structure functions to be computed.

#### `test: test_switch`

You can enter `true` or `false`

`true`: For running in test mode. Idealized velocity and scalar fields are generated internally by the code. Computed structure functions are compared with analytical results. The code is PASSED if the percentage difference between the two results is less than `1e-10`.

`false`: The "regular" mode, in which the code reads the fields from the hdf5 files in the `in` folder.

### iii) Running Instructions and Command-Line Arguments 
To run `fastSF`, change to `fastSF` directory. Ensure that the input hdf5 files follow the schema described in the previous subsection. If you want all the relevant parameters to be read from "in/para.yaml", you can simply type the following:

`mpirun -np [number of MPI processors] src/fastSF.out`
 
`fastSF` allows the user to pass the input parameters using command line arguments as well. If inputs are provided via the command-line, the corresponding inputs read from the "in/para.yaml" file get overriden. The user can also specify the input and output hdf5 file names via the command-line. The command line arguments are given as follows:

`mpirun -np [number of MPI processors] src/fastSF.out -s [scalar_switch]` 
`-d [2D_switch] -l [Only_longitudinal] -p [Processors_X] -X [Nx] -Y [Ny]`
`-Z [Nz] -x [Lx] -y [Ly] -z [Lz] -1 [q1] -2 [q2] -t [test_switch]`
`-U [Name of the hdf5 file containing the dataset storing Ux]`
`-u [Name of the dataset storing Ux]`
`-V [Name of the hdf5 file containing the dataset storing Uy]`
`-v [Name of the dataset storing Uy]`
`-W [Name of the hdf5 file containing the dataset storing Uz]`
`-w [Name of the dataset storing Uz]`
`-Q [Name of the hdf5 file containing the dataset storing T (scalar)]`
`-q [Name of the dataset storing T (scalar)]`
`-P [Name of the hdf5 file storing the transverse structure functions]`
`-L [Name of the hdf5 file storing the longitudinal structure functions]`
`-h [Help]`

The user need not give all the command line arguments; the arguments that are not provided will be read by the `in/para.yaml` file. For, if the user wants to run `fastSF` with 16 processors with 4 processors in x direction, and wants to compute only the longitudinal structure functions, the following command should be entered:

`mpirun -np 16 ./src/fastSF.out -p 4 -l true`

In this case, the number of processors in the x-direction and the longitudinal structure function switch will be taken via the command line. The rest of the parameters will be taken from the `in/para.yaml` file.

**Note:** `Nx`, `Ny`, and `Nz` should be specified only if test case is "on", in which case the code generates the input fields.

### iv) Output Information

Unless specified otherwise by the user via command-line arguments, the following output files are written by `fastSF`.

**Velocity structure functions**:

The logitudinal and transverse structure functions of order `q` are stored in the files `SF_Grid_pll.h5` and `SF_Grid_perp.h5` respectively. `SF_Grid_pll.h5` and  `SF_Grid_perp.h5` have datasets named `SF_Grid_pll`+`q`  and `SF_Grid_perp`+`q` respectively. These datasets store two/three dimensional arrays for two/three dimensional input fields.

**Scalar structure functions**:

The structure functions of order `q` are stored in the file `SF_Grid_scalar.h5` consisting of the datasets named `SF_Grid_scalar`+`q`. 

## Memory Requirements

The memory requirement per processor for running `fastSF` depends primarily on the resolution of the grid. The memory requirement also depends on the number of orders of the structure functions to be computed, number of processors *P*, and the distribution of processors in *x* and *y* (or *z*) directions. The memory requirement *M* (in bytes) can be estimated as follows:

### Two dimensional scalar field:

*M* = (20 + 4*n*)*N<sub>x</sub>N<sub>z</sub>* + 8(*N<sub>x</sub>/p<sub>x</sub> + N<sub>z</sub>p<sub>x</sub>/P*) + 32*P*.

### Three dimensional scalar field:

*M* = (16 + 2*n*)*N<sub>x</sub>N<sub>y</sub>N<sub>z</sub>* + 4*N<sub>x</sub>N<sub>y</sub>* + 8(*N<sub>x</sub>/p<sub>x</sub> + N<sub>z</sub>p<sub>x</sub>/P*) + 40*P*.

### Two dimensional vector field:


*M* = (44 + 4*n*)*N<sub>x</sub>N<sub>z</sub>* + 8(*N<sub>x</sub>/p<sub>x</sub> + N<sub>z</sub>p<sub>x</sub>/P*) + 32*P*, if only longitudinal structure functions are to be computed:

*M* = (44 + 8*n*)*N<sub>x</sub>N<sub>z</sub>* + 8(*N<sub>x</sub>/p<sub>x</sub> + N<sub>z</sub>p<sub>x</sub>/P*) + 40*P*, if both longitudinal and transverse structure functions are to be computed.

### Three dimensional vector field:

*M* = (56 + 2*n*)*N<sub>x</sub>N<sub>y</sub>N<sub>z</sub>* + 4*N<sub>x</sub>N<sub>y</sub>* + 8(*N<sub>x</sub>/p<sub>x</sub> + N<sub>z</sub>p<sub>x</sub>/P*) + 40*P*, if only longitudinal structure functions are to be computed.

*M* = (56 + 4*n*)*N<sub>x</sub>N<sub>y</sub>N<sub>z</sub>* + 4*N<sub>x</sub>N<sub>y</sub>* + 8(*N<sub>x</sub>/p<sub>x</sub> + N<sub>z</sub>p<sub>x</sub>/P*) + 48*P*, if both longitudinal and transverse structure functions are to be computed.

In the above expressions, *p<sub>x</sub>* refers to the number of processes in *x* direction and *P* refers to the total number of processors. Note that for large *N<sub>z</sub>*, the first term dominates the remaining terms; thus the memory requirement can be quickly estimated using the first term only. 

## Documentation and Validation

The documentation can be found in `fastSF/docs/index.html`. 

The validation of `fastSF` is reported [here](https://github.com/ShubhadeepSadhukhan1993/fastSF/blob/master/docs/Verification.md).

## Contributions and Bug Reports

We welcome contributions to this project. If you wish to contribute, please create a branch with a [pull request](https://github.com/ShubhadeepSadhukhan1993/fastSF/pulls) and the changes can be discussed there.

If you find a bug in the or errors in the documentation, please open a [new issue](https://github.com/ShubhadeepSadhukhan1993/fastSF/issues/new) in the Github repository and report the bug or the error. Please provide sufficient information for the bug to be reproduced.  

## License

`fastSF` is released under the terms of BSD New License.

# Testing and Validation of `fastSF`

We validate `fastSF` by comparing the numerical results with analytical results for idealized **u** and *&theta;* fields as well as with the predictions of K41 (Kolmogorov 1941a, 1941b).

### Problem 1

We consider the following 2D velocity and scalar fields: 

**u** = [*x, z*]

*&theta; = x + z*

For the above fields, it can be analytically shown that the longitudinal velocity structure functions and the scalar structure functions are

*S<sub>q</sub><sup>u<sub>ll</sub></sup> = (l<sub>x</sub><sup>2</sup> + l<sub>z</sub><sup>2</sup>)<sup>q/2</sup> = l<sup>q</sup>*,

*S<sub>q</sub><sup>&theta;</sup> = (l<sub>x</sub> + l<sub>z</sub>)<sup>q</sup>*.

We run ``fastSF`` to compute the velocity and scalar structure functions for the above fields. The resolution of the fields and the domain size are 32<sup>2</sup> and 1 x 1 respectively. We plot the second and the third-order longitudinal velocity structure functions versus *l* in Fig.1. Clearly, *S<sub>2</sub><sup>u<sub>ll</sub></sup> (l)* and *S<sub>3</sub><sup>u<sub>ll</sub></sup> (l)* equal *l<sup>2</sup>* and *l<sup>3</sup>* respectively, consistent with the analytical results. Figure 2 exhibits the density plots of the computed second-order scalar structure function *S<sub>2</sub><sup>&theta;</sup>* (**l**) along with (*l<sub>x</sub> + l<sub>z</sub>*)<sup>2</sup>. The two plots are very similar, thus showing the ``fastSF`` computes the scalar structure function correctly.


 <figure>
  <img src="figs/SF_velocity_r2D.png" alt="Trulli" style="width:70%">
  <figcaption>Fig.1: For the velocity field defined in Problem 1: plots of the second and third-order longitudinal structure functions vs. l. The second and third-order structure functions equal l<sup>2</sup> and l<sup>3</sup> respectively.</figcaption>
</figure>
 
<figure>
  <img src="figs/SF_scalar2D.png" alt="Trulli" style="width:100%">
  <figcaption>Fig.2: For the scalar field defined in Problem 1: (a) Density plot of the second-order scalar structure function as function of the displacement vector. (b) Density plot of (l<sub>x</sub> + l<sub>z</sub>)<sup>2</sup>, which is the analytical value of the second-order scalar structure function. The two density plots are very similar.</figcaption>
</figure>


### Problem 2
We consider the following 3D velocity and scalar fields: 

**u** = [*x, y, z*]

*&theta; = x + y + z*

For the above fields, it can be analytically shown that the longitudinal velocity structure functions and the scalar structure functions are

*S<sub>q</sub><sup>u<sub>ll</sub></sup> = (l<sub>x</sub><sup>2</sup> + l<sub>y</sub><sup>2</sup> + l<sub>z</sub><sup>2</sup>)<sup>q/2</sup> = l<sup>q</sup>*,

*S<sub>q</sub><sup>&theta;</sup> = (l<sub>x</sub> + l<sub>y</sub> + l<sub>z</sub>)<sup>q</sup>*.

We run ``fastSF`` to compute the velocity and scalar structure functions for the above fields. The resolution of the fields and the domain size are 32<sup>3</sup> and 1 x 1 x 1 respectively. We plot the second and the third-order longitudinal velocity structure functions versus *l* in Fig.3. Clearly, *S<sub>2</sub><sup>u<sub>ll</sub></sup> (l)* and *S<sub>3</sub><sup>u<sub>ll</sub></sup> (l)* equal *l<sup>2</sup>* and *l<sup>3</sup>* respectively, consistent with the analytical results. Figure 4 exhibits the density plots of the computed second-order scalar structure function *S<sub>2</sub><sup>&theta;</sup>* (*l<sub>x</sub>, l<sub>z</sub>*) on *l<sub>y</sub>* = 0.5 plane, along with (*l<sub>x</sub>* + 0.5 + *l<sub>z</sub>*)<sup>2</sup>. The two plots are very similar, thus showing the ``fastSF`` computes the scalar structure function correctly.


Problems 1 and 2 are used as test cases for the the code. The user is required to execute the shell script `fastSF/runTest.sh` to run the test case. On doing so, the code generates the velocity and the scalar fields as per the above relation. After computing the structure functions, the code computes the percentage difference between the theoretical and the computed values of the structure functions. If the error does not exceed 1 x 10<sup>-10</sup>, the code is deemed to have passed.

<figure>
  <img src="figs/SF_velocity_r3D.png" alt="Trulli" style="width:70%">
  <figcaption>Fig.3: For the velocity field defined in Problem 2: plots of the second and third-order longitudinal structure functions vs. l. The second and third-order structure functions equal l<sup>2</sup> and l<sup>3</sup> respectively.</figcaption>
</figure>

<figure>
  <img src="figs/SF_scalar3D.png" alt="Trulli" style="width:100%">
  <figcaption>Fig.4: For the scalar field defined in Problem 2: (a) Density plot of the second-order scalar structure function as function of the displacement vector on *l<sub>y</sub>* = 0.5 plane. (b) Density plot of (l<sub>x</sub> + 0.5 + l<sub>z</sub>)<sup>2</sup>, which is the analytical value of the second-order scalar structure function for *l<sub>y</sub>* = 0.5. The two density plots are very similar.</figcaption>
</figure>

### Problem 3

Here, we consider the classical results of Kolmogorov (Kolmogorov 1941a, 1941b) for 3D incompressible hydrodynamic turbulence with homegeneity and isotropy. In such flows, for the inertial range, which comprises of scales lying between the large-scale forcing regime and the small-scale dissipation regime, the third-order longitudinal velocity structure function is given by 

*S<sub>3</sub><sup>u<sub>ll</sub></sup> (l)* = -(4/5)*&epsilon;l*, 


where *&epsilon;* is the viscous dissipation rate (Kolmogorov1941a, 1941b; Frisch 1995; Verma 2019). For an arbitrary order *q*, She and Leveque (1994) proposed that the longitudinal structure functions scale as *S<sub>3</sub><sup>u<sub>ll</sub></sup> (l)* ~ *l<sup>&zeta;<sub>q</sub></sup>*, where the exponent *&zeta;<sub>q</sub>* is given by 


*&zeta;<sub>q</sub>* = (*q*/9) + 2 {1 - (2/3)<sup>*q*/3</sup> }.


We compute the longitudinal velocity structure functions of *q* = 3, 5, 7 using the simulation data of 3D hydrodynamic turbulence with Reynolds number (Re) of 5700. The simulation was performed using TARANG (Verma et al 2013; Chatterjee et al 2018) on a 512<sup>3</sup> grid with the domain size of (2&pi; x 2&pi; x 2&pi;). For more details on the simulation, refer to (Sadhukhan et al 2019). We run ``fastSF`` on Shaheen II to compute the structure functions, employing 4096 MPI processes. 




We normalize the third, fifth, and seventh-order longitudinal velocity structure functions with (*&epsilon;l*)<sup>*&zeta;<sub>q</sub></sup>*, where *&zeta;<sub>q</sub>* is given by She-Leveque's relation. We plot the negative of these quantities versus *l* in Fig. 5. 
The figure clearly shows that in the inertial range (0.2 < *l* < 0.8), the normalized third-order longitudinal velocity structure function is fairly close to 4/5 (represented by dashed line), consistent with Kolmogorov's theory. Moreover, the normalized fifth and seventh-order structure functions show a plateau for the same range of *l*, thus exhibiting consistency with She-Leveque's model. Note that we expect more accurate results for higher resolution simulations (Verma et al 2013).

The results obtained from Problems 1, 2, and 3 thus validate ``fastSF``. 

<figure>
  <img src="figs/SF_hydro.png" alt="Trulli" style="width:70%">
  <figcaption>Fig.5: For 3D homogeneous isotropic turbulence (Problem 3): plots of the negative of normalized third, fifth and seventh-order structure functions vs. l. The negative of the normalized third-order structure function is close to 4/5 (dashed line) in the inertial range..</figcaption>
</figure>

### Bibliography

Chatterjee, A. G., Verma, M. K., Kumar, A., Samtaney, R., Hadri, B., and Khurram, R. (2018). Scaling of a Fast Fourier Transform and a pseudo-spectral fluid solver up to 196608 cores. *J. Parallel Distrib. Comput.*, 113(3), 77-91, doi:10.1016/j.jpdc.2017.10.014

Frisch, U. (1995). *Turbulence: The Legacy of A. N. Kolmogorov*. Cambridge: Cambridge University Press. doi:10.1017/CBO9781139170666

Kolmogorov, A. N. (1941a). Dissipation of Energy in Locally Isotropic Turbulence. *Dokl Acad Nauk SSSR*, 32(1), 16–18. doi:10.1098/rspa.1991.0076

Kolmogorov, A. N. (1941b). The local structure of turbulence in incompressible viscous fluid for very large Reynolds numbers. *Dokl Acad Nauk SSSR*, 30(4), 301–305. doi:10.1098/rspa.1991.0075

Sadhukhan, S., Samuel, R., Verma, M. K., Stepanov, R., Plunian, F., and Samtaney, R. (2019). Enstrophy transfers in helical turbulence. *Phys. Rev. Fluids*, 4, 84607, doi:10.1103/PhysRevFluids.4.084607

She, Z., and Leveque, E. (1994). Universal scaling laws in fully developed turbulence. *Phys. Rev. Lett.*, 72(3), 336-339, doi:10.1103/PhysRevLett.72.336

Verma, M. K. (2019). *Energy transers in Fluid Flows: Multiscale and Spectral Perspectives*. Cambridge: Cambridge University Press. doi:10.1017/9781316810019

Verma, M. K., Chatterjee, A. G., Reddy, S., Yadav, R. K., Paul, S., Chandra, M., and Samtaney, R. (2013). Benchmarking and scaling studies of pseudospectral code Tarang for turbulence simulations. *Pramana-J. Phys*, 81(4), 617-629, doi:10.1007/s12043-013-0594-4




