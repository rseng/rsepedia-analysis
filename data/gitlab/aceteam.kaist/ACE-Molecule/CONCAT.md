# Installing ACE-Molecule

## Basic Instruction

Compiling ACE-Molecule requires CMake version 2.8 or higher. Like other project using CMake, installation is done by:
```
cd PATH/TO/ACE-Molecule
cmake [options] .
make -j4
# make install is not required.
```

ACE-Molecule requires some arguments for invoking CMake. See *Example CMake commands* and *CMake argument description* below. 
Also, we recommend parallel make (using argument -j). Optimal value of threads depend on your number of processors.
If you successfully compiled ACE-Molecule, you will get execuatable file named *ace*. ACE-Molecule is single file executable and this can be freely moved to other path.

To clean the project, invoke:
```
find -iname '*cmake*' -not -name CMakeLists.txt -not -name Find_Trilinos.cmake -not -name Version.cmake -not -name FileLists.cmake  -exec rm -rf {} \+
rm Makefile ace
```

To compile the developer API document, invoke:
```
doxygen Doxyfile
```
To clean the developer API document, invoke:
```
rm -r doxygen
```


## Dependency list
### Mandatory
 - Trilinos. Version 12.8.1 and 12.0.1 is tested. Note that the compiler for ACE-Molecule is fixed to the compiler which compiled the Trilinos library. Available at https://trilinos.org/download/public-git-repository/ .
 - Libxc Version 4.0 or higher is required. Available at http://www.tddft.org/programs/libxc/download/ .
 - GNU scientific library (GSL, will be dropped). Available at https://www.gnu.org/software/gsl/ .
### Optional
 - libxml2 for PAW support.
 - PCMSolver for PCM solvation support. Available at https://pcmsolver.readthedocs.io/en/stable/ .
 - MKL for deflation lanczos diagonalization support and performance.

## CMake argument descriptions
 - TRILINOS_PATH, LIBXC_PATH, GSL_PATH: Path to corresponding library.
 - MKL_PATH: Path to MKL. If present, deflation lanczos will be enabled by default.
 - PCMSOLVER_PATH: Path to PCMSolver. If present, PCM solvation will be enabled by default.
 - LIBXML2_PATH: Path to libxml2. If present or libxml2 is detected, PAW will be enabled by default.
 - ENABLE_CI: Compile CI routines.
 - ENABLE_EX_DIAG: Compile experimental diagonalization routines.
 - DEFLATION_LANCZOS: Compile deflation lanczos.

## Example CMake commands
### With MKL
```
cmake \
-D TRILINOS_PATH:PATH=/path/to/trilinos/ \
-D LIBXC_PATH:PATH=/path/to/libxc/ \
-D GSL_PATH:PATH=/path/to/gsl/ \
-D CMAKE_CXX_FLAGS:STRING="-g -std=c++11 -Wall -Wextra -O3" \
-D PCMSOLVER_PATH:PATH=/path/to/pcmsolver/root \
-D MKL_PATH:PATH=/path/to/mkl/ \
-D DEFLATION_LANCZOS:BOOL=TRUE \
-D ENABLE_CI:BOOL=TRUE \
-D ENABLE_EX_DIAG:BOOL=TRUE \
./   \
```

### Without MKL
```
cmake \
-D LIBXC_PATH:PATH=/path/to/libxc/ \
-D TRILINOS_PATH:PATH=/path/to/trilinos/ \
-D CMAKE_CXX_FLAGS:STRING="-g -std=c++11 -O3 -Wall -Wextra -Wno-sign-compare -I/usr/include/mpi" \
-D GSL_PATH:PATH=/path/to/gsl/ \
-D DEFLATION_LANCZOS:BOOL=FALSE \
-D ENABLE_CI:BOOL=FALSE \
./
```

### With CUDA
```
cmake \
-D LIBXC_PATH:PATH=/path/to/libxc/ \
-D CMAKE_CXX_FLAGS:STRING="-g -std=c++11 -Wall -Wextra -fPIC" \
-D CUDA_TOOLKIT_ROOT_DIR:PATH=/path/to/cuda/  \
-D CUDA_FLAGS:STRING=--cudart=shared \
-D GSL_PATH:PATH=/path/to/gsl/ \
-D DEFLATION_LANCZOS:BOOL=FALSE \
-D MKL_PATH:PATH=/path/to/mkl/ \
-D MAGMA_PATH:PATH=/path/to/magma/ \
-D TRILINOS_PATH:PATH=/path/to/trilinos/ \
./ 
```
# Contribution Guide

## How to open a bug report
Follow the template below to post a bug after you searched for similar bug.
1. General description about this bug, including what are you trying to do.
2. Expected behavior.
3. Current behavior.
4. Input files and relevant output/error log.
5. Your system environment.
6. Possible solutions.

## Coding convention
1. Use UNIX-type newline.
2. Indentation should be 4 whitespaces.
3. if-else if-else clause should be:

    ```
    if(...){
        foo();
    }
    else {
        bar();
    }
    ```
4. Avoid using bare pointer and use std library or Teuchos::RCP, unless the pointers are requested by external libraries.
5. A class should access to Teuchos::ParameterList only at constructor or dedicated initializer if it should.

## Useful manuals to read
 - ACE-Molecule developer API: http://wooyoun.kaist.ac.kr/ACE/doxygen/
 - Epetra library developer API: https://trilinos.org/docs/dev/packages/epetra/doc/html/index.html
 - Teuchos library developer API: https://trilinos.org/docs/dev/packages/teuchos/doc/html/index.html
 - MPI function guideline from KISTi (Korean): http://k-atoms.ksc.re.kr/pdf/mpi_lec.pdf

## Brief introductions to the code structure
- ACE-Molecule is a highly structured OOP code. 
- All computation routine classes, which are responsible for SCF, TDDFT, or etc. calculations are inherited from **Compute_Interface** abstract class.
- **Compute_Interface** classes are the class that directly accessed from main routine. 
**Compute_Interface** class has **compute** public method, which are called when the actual calculation is requested. 
This class is recommended starting point to inspect. For example, if you want to make changes affecting SCF routine, the **compute** method of **Scf** class is recommended starting point. In principle **Compute_Interface** classes have equal hierachy.
- **Compute** classes are the classes that computes something - like exchange correlations, Hartree potential, and Hessians.
**Compute** *directory* includes both **Compute_Interface** and **Compute** classes, since Compute_Interace classes are actually children of Compute classes.
- **Core** directory are the core modules of ACE-Molecule. Examples include the classes which performs diagonalization or provids the generalized pseudopotential. 
- **Basis** class contains all basis informations.
The **Basis** *directory* contains numerous classes, but always use **Basis** class generated by **Create_Basis** unless you exactly know what you are doing.
Usually you do not need to generate this class - Just use the object already generated.
- All calculation results are stored in the **State** class or its child class, and all basis informations are contained in the **Basis** class.
- Most classes and pointers are controled via **Teuchos::RCP** (smart reference counting pointer) for garbage collections.


### Introduction to common variables
 - **RCP&lt;Atoms&gt; atoms**: Contains the basic information about calculating molecule. The geometry, atom symbols and atomic numbers, etc. It is stored in **State** object.
   Pseudopotential-related contents, like the number of valence electrons, are not contained here.
 - **RCP&lt;Basis&gt; mesh**: Contains the information about the simulation box. 
   Mapping between flattend index (used as an index of Epetra_MultiVector) and x, y, or z-direction index is done via **decompose** and **combine** method.
   Mapping from the x, y, or z- direction index to the actual coordinate is done via **get_scaled_grid** method.
   This object also holds how orbitals or density are distributed over processors. **get_map** method returns **Epetra_Map** object, which holds such informations.
 - **RCP&lt;State&gt; state**: Contains the orbital, density, energy and its components, etc.
 - **RCP&lt;Teuchos::ParameterList&gt; parameters**: Contains input parameter informations.
   Each **Compute_Interface** classes are provided with their own parameters and the subsection **Basic_Information**.
 - Note that **orbitals, densities, potentials, and occupations** are often stored as **Array&lt;RCP&lt;Epetra_(Multi)Vector&gt; &gt;**. Array index is **spin index** if exists, and MultiVector index is orbital index, etc..

## Common errors to note
 - Wrong RCP.

```
RCP<Epetra_MultiVector> a = rcp(new Epetra_Multivector(...));
RCP<Epetra_Vector> b = rcp(a -> operator()(0));
// Segfaults when b is destroyed and a is destroyed (double free).
```
Use following code instead.
```
RCP<Epetra_Vector> b = rcp(new Epetra_Vector(*a -> operator()(0)));
```
 - Accessing Epetra_Map inside OpenMP-parallelized for loop.

Somehow following code segfaults sometimes when using OpenMP threading.
```
#pragma omp parallel for
for(int i = 0; i < size; ++i){
    int j = epetra_map.GID(i);// Somehow cause bug.
}
```
Use following code instead
```
int* gid_list = epetra_map.MyGlobalElements()
#pragma omp parallel for
for(int i = 0; i < size; ++i){
    int j = gid_list[i];
}
```
 - Computed values not properly brodcasted over processors. Make sure that all processors have same values or properly distributed over the processors by suitable object.
 - Use std::abs from <cmath> always and avoid int abs from <cstdlib> being used for float numbers.
# Advanced Computational Engine for Molecules (ACE-Molecule)

You can find more details on [wiki page](https://gitlab.com/aceteam.kaist/ACE-Molecule/wikis/home)
## **Diagonalize**

This section contains parameters related to the diagonalization for sparse matrix.

### **Tolerance**
> #### **Description**
>  Specify diagonalization tolerance. \
> Relative eigenvalue tolerance of solver. \
> If norm2(residual)/eigenvalue is lower than this variable for all eigenstates, diagonalization is considered converged.
> #### **Type**
> string

### **Solver**
> #### **Description**
> Specify Diagonalization Solver.
> #### **Type**
> string
> #### **Possible Options**
> + Pure: Use Anasazi [Default]
> + Direct: Use LAPACK default diagonalization scheme.

### **EigenSolver**
> #### **Description**
> Specify the method to diagonalize symmetric matrix
> #### **Type**
> string
> #### **Possible Options**
> + LOBPCG: Use LOBPCG (Locally Optimal Block Preconditioned Conjugated Gradient) method [Default]
> + BlockDavidson: Use BlockDavision method. It required to compile the program with USE_EX_DIAG set to true.
> + BlockKrylovSchur: Use BlockKrylovSchur method. It required to compile the program with USE_EX_DIAG set to true.

### **RandomizeInitial**
> #### **Description**
> This parameter determine how to generate initial vector for DeflationLanczos method
> #### **Type**
> int
> #### **Possible Options**
> + 0: initial vector is not randomized.
> + Otherwise: Randomly initialize initial vector [default]

### **DiagonalizeVerbosity**
> #### **Description**
> Specify verbosity of diagonalization
> It only works with Solver = Pure
> Simple, Normal, Debug
> #### **Type**
> string
> #### **Possible Options**
> + Simple
> + Normal
> + Debug

### **Locking**
> #### **Descirption**
> This parameters determine whether Locking is employed or not. It only works with *Solver* = Pure
> #### **Type**
> int
> #### **Possible Options**
> + 0: Do not use Locking.
> + Otherwise: Use Locking.

### **LockingTolerance**
> #### **Description**
> Set Locking Tolerance. \
> It only works with Solver = Pure.
> #### **Type**
> float
Description
END

### **MaxLocked**
> #### **Description**
>Set the maximum number of locking vector
> BlockSize + MaxLocked is greater than or equal to number of eigenvalues that you want to solve
> It only works with Solver = Pure and Locking != 0
> #### **Type**
> int

### **DiagonalizeSolverPrec**
> #### **Description**
> Specify types of preconditioner for diagonalization. It only works with Solver = Pure and Locking != 0
> #### **Type**
> string
> #### **Possible Options**
> + 'IC stand-alone': [experimental]
> + 'ILU stand-alone': [experimental]
> + 'ILUT stand-alone': [experimental]
> + 'Amesos stand-alone': [experimental]
> + 'ICT stand-alone': [experimental]
> +  ml: multi-level preconditioner please refer ml homepage for details
> +  None: no preconditioner is used



### **Redistribution**
> #### **Description**
> Determine whether redistribution is used or not
> #### **Type**
> int
> #### **Possible Options**
> + 0: not used [Default]
> + Otherwise:   used


### **BalancingTolerance**
> #### **Description**
> Determine how much imbalance is allowed. More information can be found in Isorropia homepage.
>#### **Type**
>float

### **PartitioningMethod**
> #### **Description**
> Criteria to determine whether eigenvectors from two different MPIGroup are the same or not. This parameter works only when Deflation Lanczos is used.
> #### **Type**
> string
> #### **Possible Options**
> + HYPERGRAPH [Default]


## CISD

This section controls parameters related to CISD

(J. Chem. Phys. 2016, 145, 224309. Jaechang Lim, Sunghwan Choi, Jaewook Kim, and Woo Youn Kim*)

#### **Subsection List**
  +  Diagonalize 

### **NumberOfEigenvalues**
> #### **Description**
> Number of eigenvalues which will be obtained from diagonalization of CISD matrix. 
> #### **Type**
> int
> #### **Possible Options**
> + key: 10 # if you want to calculate the 10 lowest states.

### **ConvergenceTolerance**
> #### **Description**
> Convergence criteria of Davidson digonalization. Relevant only if Davidson diagonalizer is used.<br>
> Defaults to 1.0E-6.
> #### **Type**
> float

### **Method**
> #### **Possible Options**
> + 0: Use Davidson diagonalization.
> + 1: Make full matrix and diagonalize it.
## **Scf**

 This section governs variables related to Scf

#### **Subsection List**
  +  Diagonalize (Optional)
  +  Exchange_Correlation (Mandatory)
  +  ISF (Optional)
  +  Mixing (Optional)
  +  ExternalField (Optional)
  +  SolvationModel (Optional)
  +  Output (Optional)


### **VarName**
> #### **Description**
> Description blah blah
> #### **Type**
> string / int / (positive) float
> #### **Possible Options**
> + key: description


### **IterateMaxCycle**
> #### **Description**
> Number of max iteration of Scf.
> #### **Type**
> (positive) int

### **ConvergenceType**
> #### **Description**
> Type to be used to calculate convergence
> #### **Possible Options**
> + Energy:  Energy difference will be used for checking convergence
> + Density: difference will be used for checking convergence [Default]
> + Potential: Potential difference will be used for checking convergence
> + HOMO: HOMO energy difference will be used for checking convergence
> + EigenvalueSum: Sum of eigenvalue difference will be used for checking convergence

### **ConvergenceTolerance**
> #### **Description**
> Scf will be end if difference of *ConvergenceType* is smaller than *ConvergenceTolerance*.
> #### **Type**
> (positive) float 

### **EnergyDecomposition**
> #### **Description**
> Calculates energy per value of *EnergyDecomposition* steps.
> If the value is 0, the program does not calculate energy at each steps; if it is 3, the program calculates energy for every 3 steps. 
> The converged total energy and its components are always calculated and printed, regardless of this option.
> Defaults to 1 if ConvergenceType is Energy, 0 otherwise.
> #### **Type**
> non-negative int


### **NumberOfEigenvalues**
>#### **Description**
> It decides number of eigenvalues to be calculated in Scf. If this option does not exists or non-positive, it defaults to the eigenvalue number of previous routine.
>#### **Type**
> int

### **DiagonalizeShouldNotBeConverged**
>#### **Description**
> This parameter decides whether diagonalization if fully converged or not. If this option is not specified, diagonalization will be fully converged.
>#### **Type**
> int
> #### **Possible Options**
> + 0: Diagonalization will be fully converged.
> + Otherwise: Diagonalization will not be fully converged. [default]

### **MakingHamiltonianMatrix**
> #### **Description**
> it decides whether diagonalization with making hamiltonian matrix or not. \
> For default, diagonalization will be done with making hamiltonian matrix.
>#### **Type**
> int 
> #### **Possible Options**
> + 0:  Diagonalization will be done without making hamiltonian matrix.
> + Otherwise:  Diagonalization will be done with making hamiltonian matrix. [default]

### **ConstructKineticMatrix**
> #### **Description**
> It decides whether calculating kinetic energy using matrix or not.
> #### **Type**
> int
> #### **Possible Options**
> + 0 Calculate kinetic energy without explicit construction of kinetic matrix.
> + Otherwise: Calculate kinetic energy using kinetic matrix. [default]



### **ComputeInitialEnergy**
> #### **Description**
> Decides whether the energy of initial state is computed or not
> #### **Type**
> int
> #### **Possible Options**
> + 0: Do not compute initial energy
> + Otherwise: Compute initial energy [default]

### **IgnoreInternalInit**
> #### **Description**
> Controls the method to obtain initial PAW atomcenter density matrix and hartree potential for PAW calculations.
> #### **Type**
> int 
> #### **Possible Options**
> 0: Construct initial PAW atomcenter density matrix from PAW dataset file. Good if the system is similar with the combination of the free atoms. Default value.
> 1: Construct initial PAW atomcenter density matrix from initial guess orbitals. Good if cube guess is used


### **Making_Hamiltonian_Matrix**
> #### **Description**
> Optional value, it decides whether diagonalization with making hamiltonian matrix or not. \
> For default, diagonalization will be done with making hamiltonian matrix.
> #### **Type**
> int
>+ 0: Diagonalization will be done without making hamiltonian matrix.
>+ Otherwise: Diagonalization will be done with making hamiltonian matrix. [default]

## **TDDFT**
 This section governes parameters related to TDDFT calculations.

#### **Subsection List**
  +  ExchangeCorrelation (Mandatory)
  +  OrbitalInfo (Mandatory if ExchangeKernel is HF_EXX) \
 Contains exchange correlation functional which used to compute the orbitals.\
 Necessary if ExchangeKernel is HF_EXX.\
 Should contain ExchangeCorrelation as subsection.

### **NumberOfStates**
> #### **Description**
> This varible controls the number of excitations computed by TDDFT.
> #### **Type** 
> (positive)int

### **TheoryLevel**
> #### **Description**
> This parameter determines the way to solve TDDFT.
> #### **Type** 
> string
> #### **Possible Options**
> + Casida: It gives oscillator_strength and negative excitations are considered. [default]
> + TammDancoff: Negative excitations are ignored.


### **SortOrbital**
> #### **Description**
> This parameter controls how to print components of
> #### **Type** 
> string
> #### **Possible Options**
> + Order: Print top few components. The number is given as MaximumOrder. [default]
> + Tolerance: Print out all components that are bigger than given criteria.

>
### **MaximumOrder**
> #### **Description**
> This determines the number of printed orbital pairs. This keyword works only with SortOrbital=Sorting case.
> #### **Type** 
> int


### **OrbitalTolerance**
> #### **Description**
> If the portion of certain orbital is larger than this critera, it will be printed.
> #### **Type** 
> float


### **GradientMatrix**
> #### **Description**
> This parameter determines how to obtain the gradients of density and orbital.
> #### **Type** 
> string
> #### **Possible Options**
> + Finite_Difference: Use the finite difference method to obtain the gradient.

### **DerivativesOrder**
> #### **Description**
> Decides derivative order for calculating GradientMatrix.\
> Only relevant when GradientMatrix is Finite_Difference.
> #### **Type** 
> int

### **ExchangeKernel**
> #### **Description**
> Determines the exchange kernel for exact exchange.
> #### **Type** 
> string
> #### **Possible Options**
> + PGG: PGG kernel. [default]
> + HF or HF_EXX: Perform TDHF or TD-HFKS calculation scheme like J. Chem. Phys. 134, 034120 (2011), TDHF(EXX). Further description can be found in DeltaCorrection.


### **DeltaCorrection**
> #### **Type** 
> int
> #### **Description**
> Specifies eigenvalue correction term. Only relevant if ExchangeKernel is HF or HF_EXX.
> #### **Possible Options**
> + 1: Calculate as J. Chem. Phys. 134, 034120 (2011). [default]
> + 2: Include only diagonal elements for occupied orbital eigenvalue correction term. Not recommended.
> + 3: Include only diagonal elements for occupied orbital eigenvalue correction term. Not implemented yet.
> + 0: Do not include eigenvalue correction term. Identical to ExchangeKernel==HF and fallbacks to it.

### **_TDClass**
> #### **Type** 
> string
> #### **Description**
> Debug input. Chooses the TDDFT calculation class. \
> If not specified (normal behavior), automatically choose appropreate routine. \
> If specified class cannot compute desired theory, calculation stops.
> #### **Possible Options**
> + C: \
> Invoke Casida equation $ CZ = -\Omega^2 Z $ computation routine. HF and KS-CI exchange kernels are not implemented yet.\
> This is chosen if the TheoryLevel is Casida and EXX exchange kernel is not HF or KS-CI.
> + ABBAL: \
> Invoke TDDFT equation $ (A* B* | B A) (X Y) = -\Omega (X Y) $ computation routine.\
> This is chosen if the TheoryLevel is Casida and HF or KS-CI exchange kernel is used.\
> This class should be removed (with this input) when HF, KS-CI casida form is implemented.
> + TDA:\
> Invoke Tamm-Dancoff equation computation routine.\ 
> This is chosen is the TheoryLevel is TDA.



### **Gradient**
> #### **Description**
> Determine how to calculate gradient of pair density.\
> #### **Type** 
> int
> #### **Possible Options**
> + 0: gradient of pair density is calculated using gradient matrix [default]
> + Otherwise: gradient of pair density is calculated using saved orbital gradient.


### **Root**
> #### **Description**
> Does not count excitations from the orbitals lower than the orbitals with this index.\
> Root-1 orbitals will be ignored.
> #### **Type** 
> int


### **MemorySaveMode**
> #### **Description**
> keeping result of kernel integration loaded on memory for acceleration
> #### **Type** 
> int
> #### **Possible Options**
> + 0: calculation of iajb-pair will use higher memory and lower communication time. [default] 
> + Otherwise: it will use lower memory but needs more communication time.

## **Guess**

 Parameters in this section govern initial guess of a target system.

#### **Subsection List**
  +  Diagonalize (Optional)
  +  Guess       (Optional, for Grid_Cutting)
  +  Scf         (Optional, for Grid_Cutting)


### **InitialGuess**
> #### **Description**
> Decides what kind of initial guess methods are used.
> #### **Type**
> int 
> #### **Possible Options**
> + -1: ramdom guess, Randomly generate guess orbitals
> + 0: Core Hamiltonian guess, although it does not really use the core Hamiltonian.\
> + 1: Atomic density guess [default], Infer initial orbitals from atomic density from UPF pseudopotential file. InitialFilenames should be given as UPF type pseudopotential. 
> + 2: Extended Hueckel guess, This guess need EHTFilename parameter
> + 3: Cube guess, This guess need EHTFilename parameter. For TDDFT calculations, Info and Info_Filetype variable is required.
> + 4: PAW guess, Infer initial orbitals from PAW data file. InitialFilename should be given as PAW-XML format.
> + 5: Restart guess, Reserved for the restart option that is currently not supported. Use cube guess option, with cube file written during the calculations.
> + 6: Grid cutting guess, Perform pre-SCF to generate guess orbitals. You need to set *Guess* *Scf* inside *Guess* section, which is related to self-consistent field calcualtions for inner spher domain.
> + 7: Cube density guess, Read and use density in cube format as initial density of SCF. Grids of calculations and cube file should be aligned.

### **InitialFilenames**
> #### **Description**
> Path to the input initial guess data file. \
> If required, InitialFilenames should be given for all atoms to be calculated, and should be given in the order of increasing atomic number. \
> This variable will be inferred from InitialFilePath and InitialFileSuffix, if those two variables exist and InitialFilenames are absent. \
> For the atomic density guess, the names of the UPF pseudopotential files should be given. \
> For the cube or cube density guess, the names of the cube files should be given. \
> For the PAW guess, the names of the PAW dataseet files should be given.
> #### **Type**
> string
>
### **EHTFilename**
> #### **Description**
> Path to the extended Hueckel parameter file.\
> See [https://github.com/greglandrum/yaehmop], which Jaewook Kim should provide. \
> Necessary for Extended Hueckel guess.
> #### **Type**
> string

### **NumberOfEigenvalues**
> #### **Description**
> Number of orbitals which is to be produced by initial guess method. \
> If not given, it is inferred from the number of total electrons.
> For the open-shell calculation, number of alpha spin orbitals will be set as this value. Same for the beta spin cases.
> 
> #### **Type**
> (positive) int
> 

### **InitialFilePath**
> #### **Description**
> Path to the directory that contains initial guess data file. \
> This is used to infer InitialFilenames. \
> For atomic density and PAW guess, it is: \
> [InitialFilenames] = [InitialFilePath]/[Atom Symbol][InitialFileSuffix] \
> For cube and cube density guess, everything in this path is used.
> #### **Type**
> string

### **InitialFileSuffix**
> #### **Description**
> This is used in conjunction to InitialFilePath to infer InitialFilenames. \
> [InitialFilenames] = [InitialFilePath]/[Atom Symbol][InitialFileSuffix]
> #### **Type**
> string

### **Info**
> #### **Description**
> Name of eigenvalue information file, which is generated from Octopus or ACE. \
> Required for TDDFT calculation with cube guess. *Info_FileT*
> #### **Type**
> string 

### **InfoType**
> #### **Description**
> Type of info file 
> #### **Type**
> string
> #### **Possible Options**
> + ACE: ACE-Molecule output which contains eigenvalue information for corresponding cube files. [default]
> + Octopus: Info file generated by Octopus program.

### **Ratio**
> #### **Description**
> Grid cutting guess will scale the full simulation box by this variable. \
> If this is not specified, smaller simulation box information should be specified in Guess section like in Basic_Information section with same grid spacing. 
>
> #### **Type**
> (positive) float
> 

### **CubeUnit**
> #### **Description**
>  Describe the unit of the input cube file. \
> This option is only relevant for Cube and Cube density guess.
> #### **Type**
> string
> #### **Possible Options**
> + bohr: Cube files are written in the bohr unit.
> + angstrom: Cube files are written in the angstrom unit.


## **ExternalField**

 Adds an additional field during the SCF run.

### **InputType**
> #### **Description**
> Specifies format. Required.
> #### **Type**
> string
> #### **Possible Options**
> + Analytic: Apply potential in analytic form (Ex, Ey, or Ez where E is constant).
> + Read: Specifies potential by numerical grid.

### **ExternalField**
> #### **Description**
> Relevant only if InputType is Analytic. Only one default option is currently available.
> Type string
> + Electric: Only and default option.

### **ExternalFieldType**
> #### **Description** 
> Relevant only if InputType is Analytic and ExternalField is Electric. Only one default option is currently available.
> #### **Type**
> string
> #### **Possible Option**
> + Static: Only and default option.

### **ExternalFieldDirection**
> #### **Description**
> Relevant only if InputType is Analytic and ExternalField is Electric.\
> Specifies the direction of electric field. \
> Potential will be Ex, Ey or Ez where E is ExternalFieldStrength.
> #### **Possible Option** 
> + x: x direction
> + y: y direction
> + z: z direction

### **ExternalFieldStrength**
> #### **Description**
> Relevant only if InputType is Analytic and ExternalField is Electric. \
> Specifies the strength of electric field. \
> Potential will be Ed where d=x,y,z is ExternalFieldDirefction.
> #### **Type**
> float

### **PotentialFilename**
> #### **Description**
> Relevant only if InputType is Read. \
> Specifies external potential on parallelpiped numerical grid. \
> Potentials are whitespace or newline-separated, and outer loop is z-index, middle loop is y-index, and inner loop is x-index.\
> Devel Note: This behavior is currently index order dependent.
> #### **Type**
> string

### **Variable Interpolation**
> #### **Description** 
> Potential interpolation scheme.
> #### **Type** 
> string
> #### **Possible Options**
> + linear: Use trilinear interpolation.
> + cubic: Use tricubic interpolation. Default value.

### **PointX**
> #### **Description**
> Number of grid points in x-direction.
> #### **Type**
> (positive) int 

### **PointY**
> #### **Description**
> Number of grid points in y-direction.
> #### **Type**
> (positive) int

### **PointZ**
> #### **Description**
> Number of grid points in z-direction.
> #### **Type**
> (positive) int


### **PotentialOffSetX**
> #### **Description**
> Position of potential file center (x) on calculation mesh. Default value: 0.0.
> #### **Type**
> float

### **PotentialOffSetY**
> #### **Description**
> Position of potential file center (y) on calculation mesh. Default value: 0.0.
> #### **Type**
> float

### **PotentialOffSetZ**
> #### **Description**
> Position of potential file center (z) on calculation mesh. Default value: 0.0.
> #### **Type**
> float


### **LatticeVecX**
> #### **Description**
> Specifies x-direction vector of the box containing the potential.
> #### **Type**
> string (ex. "1.0 0.0 0.0")

### **LatticeVecY**
> #### **Description**
> Specifies y-direction vector of the box containing the potential.
> #### **Type**
> string (ex. "0.0 1.0 0.0")

### **LatticeVecZ**
> #### **Description**
> Specifies z-direction vector of the box containing the potential.
> #### **Type**
> string (ex. "0.0 0.0 1.0")


## **BasicInformation**

 Parameters in this section govern mesh information, spin state and a geometry of target system.

#### **Subsection List**
  +  Pseudopotential (Mandatory)

### **VerboseLevel**

> #### **Description**
> Variable set level of verbosity during the calculation
>
> #### **Type**
>  int
>
> #### **Possible Options**
>+ 0: Simple - Print only simple outputs [default]
>+ 1: Normal - Print auxiliary outputs in addition to Simple. These auxiliary outputs include restating input parameters and resource used by technical routines.
>+ 2: Detail - Extremely verbose output which are not intended to be understood by non-developers.

### **Label**
> #### **Description**
> When additional results are generated during the calculation, filenames start from *Label* 
> #### **Type**
> string


### **Mode**
>#### **Description**
>Decides what kind of calculation is to be performed.
>
> #### **Type**
> string
>
> #### **Possible Options**
>+ Auto: Calculate anything specified in the input block in the specified order.
>+ Sp: Single point calculation will be performed. Calculate using input block Guess and then Scf.
>+ SpTDDFT: First, single point calculation will be performed, then TDDFT will start using orbitals from Sp calculation. For more information, please refer TDDFT section. Calculate using input block Guess, then Scf, and then TDDFT.
>+ TDDFT: TDDFT calculation will be performed. Cube guess with Info provided is highly recommended. Calculate using input block Guess and then TDDFT.
>+ Opt: Geometry optimization will be performed. This option is currently under test.

### **GeometryFilename**
> #### **Description**
> Name of the input geometry file. Now xyz and pdb fileformats are supported (Previsouly, name of this parameter is Geomtry_Filename but now old name does not work)
> #### **Type**
> string

### **GeometryFormat**
> #### **Description**
> Input geometry file format.(Previsouly, name of this parameter is Geomtry_Format but now old name does not work)
> #### **Type**
> string
> #### **Possible Options**
> + xyz: xyz format of geometry file.
> + pdb: pdb format of geometry file.

### **Grid**
> #### **Description**
> This is a keyword to select shape of simulation space.
> #### **Type**
> string
>
> #### **Possible Options**
> + Basic: The shape of simulation box will be rectangular. The size of the simulation box is governed by *Cell*
> + Sphere: The shape of simulation box will be sphere. The radius of the simulation sphere is governed by *Cell*
> + Atoms:  The simulation box is composed of points whose distance from atom is smaller than a certain value.(*Radius* see below for details) If you do not use *AbsoluteRadius* the radius will be set as *Radius* times of predefined value (van der Waals radii).
>

### **Type**
> #### **Description**
> This variable set how to construct mesh 
> #### **Type**
> string
> #### **Possible Options**
> + Points: Based on given *Cell* & *Points* information, mesh is generated.
> + Scaling: Based on given *Cell* & *Scaling* information, mesh is generated.
> + Cube: Share mesh with designated Cube file whose filename is set in *BasisCubeFilename*.

### **Cell**
> #### **Description**
> This indicates the size of cubic simulation domain in Bohr unit 
> When *Grid* is Basic, 2 times *Cell* indicate the length of simulation box.
> Otherwise, This parameter has no effect on the calculation
> #### **Type** 
> (positive) float 

### **CellDimensionX**
> #### **Description**
> This indicates the size of z axis of simulation domain in Bohr unit <br>
> When *Grid* is Basic, 2 times *Cell* indicate the length of simulation box. If *Cell* is set, this parameter does not work. <br>
> Otherwise, This parameter has no effect on the calculation <br>
> #### **Type** 
> (positive) float 

### **CellDimensionY**
> #### **Description**
> This indicates the size of z axis of simulation domain in Bohr unit  <br>
> When *Grid* is Basic, 2 times *Cell* indicate the length of simulation box. If *Cell* is set, this parameter does not work. <br>
> Otherwise, This parameter has no effect on the calculation <br>
> #### **Type** 
> (positive) float 

### **CellDimensionZ**
> #### **Description**
> This indicates the size of z axis of simulation domain in Bohr unit  <br>
> When *Grid* is Basic, 2 times *Cell* indicate the length of simulation box. If *Cell* is set, this parameter does not work. <br>
> Otherwise, This parameter has no effect on the calculation
> #### **Type** 
> (positive) float 

### **Points**
> #### **Description**
> Defines the number of grid points in all three dimension.\
> By this variable, *PointX*, *PointY*, and *PointZ* are automatically define by its value.\
> Affects basis size.
>
> #### **Type**
> (positive) int 

### **PointX**
> #### **Description**
> Defines the number of grid points in x dimension. If *Points* is set, this parameter is ignored.
>
> #### **Type**
> (positive) int

### **PointY**
> #### **Description**
> Defines the number of grid points in y dimension. If *Points* is set, this parameter is ignored.
> 
> #### **Type**
> (positive) int

### **PointZ**
> #### **Description**
> Defines the number of grid points in z dimension. If *Points* is set, this parameter is ignored.
> #### **Type**
> (positive) int 

### **Scaling** 
> #### **Description**
> This parameter defines the distance between the grid points in all three dimension in Bohr unit
> #### **Type**
> (positive) float

### **ScalingX**
> #### **Description**
> Defines the distance between the grid points in x dimension.\
> Affects basis size.\
> Should be written in Bohr unit.
> #### **Type**
> (positive) float

### **ScalingY**
> #### **Description**
> Defines the distance between the grid points in y dimension.\
> Affects basis size.\
> Should be written in Bohr unit.
> #### **Type**
> (positive) float

### **ScalingZ**
> #### **Description**
> Defines the distance between the grid points in z dimension.\
> Affects basis size.\
> Should be written in Bohr unit.
> #### **Type**
> (positive) float 

### **Radius**
> #### **Description**
> When *Grid* is set as Atoms, then the radius of simulation spheres need to be set. If *AbsoluteRadius* is set as zero, radius of each simulation sphere become covalent radius of atoms times *Radius* value. Otherwise, *Radius* value become radius in Angstrom unit
> #### **Type**
> (positive) float 

### **AbsoluteRadius**
>#### **Description**
> This parameter determines meaning of Radius; if AbsolutRadius is set as non-zero value, *Radius* become radius of sphere itslef (unit Angstrom) otherwise, the radius become *Radius* times covalent radii 
>
>#### **Type** 
> int
>
>#### **Possible Options**
>+ 0: radius of simulation sphere is set as *Radius* times covalent radii
>+ Otherwise: *Radius* values become radius of simulation sphere (unit: Angstrom) 
 
### **StoreCoreHamiltonian**
> #### **Description**
> This parameter determines whether core Hamiltonian matrix is stored in crs format or not.
>
> #### **Type**
> int
>
> #### **Possible Options**
> + 0: The core Hamiltonian matrix is not stored on a memory in crs format
> + Otherwise: The core Hamiltonian matrix is stored
>

### **BasisCubeFilename**
> #### **Description**
> Set simulation box as specified in this cube file.
> #### **Type**  
> string


### **NumElectrons**
> #### **Description**
> Number of electrons.
> #### **Type**
> (positive) float


### **Centered**
> #### **Description**
> Decides whether the center of mass of the input molecule to be translated to the simulation box origin.
> #### **Type**
> int
> #### **Possible Options**
> 0: Do not translate the molecule. [Defalut]
> 1: The center of mass will be translated to the origin.

### **ShallowCopyOrbitals**
> #### **Description**
> Shallow copy for orbitals
> #### **Type** 
> int
> #### **Possible Options**
> 0: Use deep copy [Default]
> 1: Use shallow copy. Less memory is used.

### **Polarize**
> #### **Description**
> Decides whether the spin-polarized or spin-restricted calculation to be performed.
> #### **Type**
> int
> #### **Possible Options**
> 0: Performs spin-restricted calculation.
> 1: Performs spin-polarized calculation.

### **SpinMultiplicity**
> #### **Description**
> Spin multipliciy of the moleulces. 
> 1.0 for singlet, 2.0 for doublet, etc.
> #### **Type**
> float

### **ForceCalculation**
> #### **Description**
> Decides the calculation of atomic forces after single point calculation.
> #### **Type**
> int 
> #### **Possible Options**
> No: Do not calculate atomic forces. [Default]
> Yes: Calculate atomic forces.
### **ForceDerivative**
> #### **Description**
> Decides method to evaluate force 
> #### **Possible Options**
> #### **Type**
> string
> Potential: 
> Orbital: 
### **Basis**
> #### **Description**
> Decides the basis function.
> #### **Type**
> string 
> #### **Possible Options**
> Sinc: Basis function will be used as basis function.
> FiniteDifference : Finite difference method will be used.

### **KineticMatrix**
> #### **Description**
> Decides basis for kinetic matrix.
> #### **Type**
> string
> #### **Possible Options**
> Finite_Difference:  Use finite difference method to construct the kinetic matrix. Works only with Sinc basis set. If you use this option, you should specify *DerivativesOrder*

### **DerivativesOrder**
> #### **Description**
> Decides derivative order for calculating kinetic matrix. Only relevant if *KineticMatrix* or *Basis* is set to Finite_Difference.
> #### **Type**
> int 
> #### **Possible Options**
>+ 3:  3-points central finite difference coefficients will be used. 
>+ 5:  5-points central finite difference coefficients will be used. 
>+ 7:  7-points central finite difference coefficients will be used. 
>+ 9:  9-points central finite difference coefficients will be used. 
>+ 11: 11-points central finite difference coefficients will be used. 
>+ 13: 13-points central finite difference coefficients will be used. 
>+ 21: 21-points central finite difference coefficients will be used. 

### **UseCubicBaseOnly**
> #### **Description**
> This option is relevant only if the Grid variable is set to Atoms. When reading cube file, this option is sometimes necessary.
> #### **Type**
> int
> #### **Possible Options**
>+ 0: Allow only cubic base grid.
>+ Otherwise: Allow non-cubic base grid. [Default] 

### **AllowOddPoints**
> ##### **Description**
> This parameter determines whether odd or even number of grid points will be used.
> #### **Type**
> int
> #### **Possible Options**
>+ 0: Use even number of grid points for each axis. [Default]
>+ Otherwise: Allow an odd number of points.

## DDA

This section controls parameters related to DDA (ADA)

(Sci. Rep. 2017, 7, 15775. Jaechang Lim, Sungwoo Kang, Jaewook Kim, Woo Youn Kim,* and Seol Ryu*)

#### **Subsection List**


### **InitialFilenames**
> #### **Description**
> Pseudopotential file name which will be used to calculate atomic density. This parameter is used when *ScaleUp* is minus value
> #### **Type**
> string

### **Gaussian**
> #### **Description**
> Output file of gaussian TDDFT calculation for calculation of atomic polarizability
> #### **Type**
> string

### **ACE**
> #### **Description**
> Output file of ACE-Molecule TDDFT calculation for calculation of atomic polarizability. 
> #### **Type**
> string

### **Txt**
> #### **Description**
> Text file of atomic polarizability.This file must have 3 columns for wavelength (in um), real, and imaginary value of polarizability respectively.
> #### **Type**
> string

### **StartWavelength**
> #### **Description**
> Minimum value of wavelength.
> #### **Type**
> float
> #### **Possible Options**
> + 300.0: minimum value of wavelength is 300.0nm [default]

### **EndWavelength**
> #### **Description**
> Maximum value of wavelength.
> #### **Type**
> float
> #### **Possible Options**
> + 700.0: maximum value of wavelength is 700.0nm [default]

### **ScaleUp**
> #### **Description**
> Scaling factor for atomic polarizability. The atomic polzability is multiplied by this value
> #### **Type**
> float
> #### **Possible Options**
> + 1.0: nothing change on atomic polarizability [default]
> + -1.0: atomic polarizability is scaled using local electron density, In this case, you must provide occupied molecular orbitals in Guess section.

### **Damping**
> #### **Description**
> Damping value
> #### **Type**
> float
> #### **Possible Options**
> + 1000.0: 1000nm is used for damping [default]

### **PDDamping**
> #### **Description**
> Position dependent damping value
> #### **Type**
> int
> #### **Possible Options**
> + 0: not using position dependent damping [default]

### **PDDampingMax**
> #### **Description**
> maximum value of position dependent damping value
> #### **Type**
> float
> #### **Possible Options**
> + 10000.0: 10000 nm [default]

### **PDDampingMin**
> #### **Description**
> minimum value of position dependent damping value
> #### **Type**
> float
> #### **Possible Options**
> + 1000.0: 1000 nm [default]

### **OutputBefore**
> #### **Description**
> Write atomic polarizbility of each atom without dipole-dipole interaction in file 
> #### **Type**
> int
> #### **Possible Options**
> + 0 : no write [default]
> + otherwise : write [default]

### **OutputAfter**
> #### **Description**
> Write atomic polarizbility of each atom with dipole-dipole interaction in file 
> #### **Type**
> int
> #### **Possible Options**
> + 0 : no write [default]
> + otherwise : write [default]

### **OutputDipole**
> #### **Description**
> Write dipole moment of each atom in file 
> #### **Type**
> int
> #### **Possible Options**
> + 0 : no write [default]
> + otherwise : write [default]


## **Output**

 This section controls output density, potential, orbitals and geometry.

### **Prefix**
> #### **Description**
> Decides the prefix of output density/potential/orbitals. \ 
> Output file name is: \
 [Prefix][Compute Interface].[density/potential/orbitals].s[Spin].ind[index].[txt/cube]
> #### **Type**
> string 

### **Density**
> ### **Description**
> Decides the format for the final density for corresponding Compute Interface
> #### **Type**
> string
> #### **Possible Options**
> + None: outputs for density are not generated. [default]
> + x: Print density value along x axis to file
> + y: Print density value along y axis to file
> + z: Print density value along z axis to file
> + cube: Print density as cube file

### **AllElectronDensity**
> ### **Description**
>  If density is not printed, this parameter does not affect output.
>  Decides that whether all-electron density will be printed, if available
>  This option is currently available only for the PAW pseudopotential method.
> #### **Type**
> int
> #### **Possible Options**
> + 0: Pseudo electron density will be printed. [default]
> + Otherwise: All electron density will be printed

### **Hartree**
> #### **Description**
> Decides the format for the final hartree potential output for corresponding Compute Interface
> #### **Type**
> string 
> #### **Possible Options**
> + None: outputs for Hartree potential are not generated. [default]
> + x: Print Hartree potential value along x axis to file
> + y: Print Hartree potential value along y axis to file
> + z: Print Hartree potential value along z axis to file
> + cube: Print local potential as cube file

### **Potential**
> #### **Description**
> Decides the format for the final mean field potential caused by electrons.
> #### **Type**
> string 
> #### **Possible Options**
> + None: outputs for potential are not generated. [default]
> + x: Print local potential value along x axis to file
> + y: Print local potential value along y axis to file
> + z: Print local potential value along z axis to file
> + cube: Print local potential as cube file

### **Orbitals**
>#### **Description**
> Decides the format for the final orbitals output for corresponding Compute Interface. 
>#### **Type**
> string
> #### **Possible Options**
> + None: outputs for orbital are not generated. [default]
> + x: Print orbital along x axis to file.
> + y: Print orbital along y axis to file.
> + z: Print orbital along z axis to file.
> + cube: Print orbital as cube file.

### **OrbitalStartIndex**
>#### **Description**
> Decides the lower boundary of range of orbitals to be printed. \
> Default value: 0 (Start from lowest orbital).
> #### **Type**
> int

### **OrbitalEndIndex**
> #### **Description**
> Decides the upper boundary of range of orbitals to be printed. \
> Default value: [Number of Orbitals]-1 (Print up to highest orbital calculated).
> #### **Type**
> int

### **Geometry**
> #### **Description**
> This parameter decide the name of geometry file. The format for output is fixed as "xyz" and full name of output file will be [Geometry].xyz
> #### **Type**
> string


## **SolvationModel**

Adds solvation model in SCF routine. Will be implemented to TDDFT routines, probably.

### **SolvationLibrary**
> #### **Description**
> Controls the solvation model library/routine to be used. 

> #### **Type**
> string
> #### **Possible Options**
> + PCMSolver: Currently only option. Requires compilation with PCMSolver. [default]

### **SolverType**
> #### **Description**
> Controls PCM theory.
> #### **Type**
> string
> #### **Possible Options**
> + CPCM: typical  conductor-like polarizable continuum model  
> + IEFPCM: Polarizable Continuum Model with the integral equation formalism 
> + None: Turns off the PCM routine. Default value.

### **ChargeWidth**
> #### **Description**
> Controls gaussian broadning of surface charge. See https://doi.org/10.1063/1.4932593 for further details. \
 This corresponds to alpha in equation A1 of aformentioned reference.
> #### **Type**
> float

### **Solvent**
> #### **Description**
> Defines the solvent. Currently explicit solvent definition is not supported.
 See http://pcmsolver.readthedocs.io/en/stable/users/input.html#available-solvents
> #### **Type**
> string

### **RadiiSet**
> #### **Description**
> See http://pcmsolver.readthedocs.io/en/stable/users/input.html
> ### **Type**
> string
> #### **Possible Options**
> + Bondi [Default]
> + UFF
> + Allinger

### **Area**
> #### **Description**
> See http://pcmsolver.readthedocs.io/en/stable/users/input.html<br>
 Default value is 0.3.
> #### **Type**
> float


### **Scaling**
>#### **Description**
> Determine whether radii is scaled or not. See http://pcmsolver.readthedocs.io/en/stable/users/input.html for details
>  #### **Possible Options**
> + 1: scales all radius by 1.2
> + 0: radii are not scaled [Default]

### **MinRadius**
>#### **Description**
> See http://pcmsolver.readthedocs.io/en/stable/users/input.html \
> Default value is 100, which is equivalent to turning this off.
> #### **Type**
> float


### **Correction**
> #### **Description**
> See http://pcmsolver.readthedocs.io/en/stable/users/input.html
> #### **Type**
> string

### **ProbeRadius**
> #### **Description**
> See http://pcmsolver.readthedocs.io/en/stable/users/input.html \
> Default value is 1.0.
> #### **Type**
> float

## **ExchangeCorrelation**

 This section controls the exchange correlation functional for the DFT or TDDFT calculations.<br>
 Some major exchange-correlation functional can be specified using by *FunctionalName*.<br>
 Otherwise, you can specify exchange-correlation functional using *XCFunctional* keyword, or you can specify exchange and correlation functional seperately using *XFunctional* and *CFunctional* keyword.<br>
 The hybridization of XC functional can be done with special input. See the *XFunctional* and *CFunctional* keywords.<br>

#### **Subsection List**


### **XCLibrary**
> #### **Description**
> This is a keyword to select library set of exchange correlation functional.
> #### **Type**
> string
> #### **Possible Options**
> + Libxc: Using xc library of http://www.tddft.org/programs/octopus/wiki/index.php/Libxc <br>
> The list of available functionals can be found in http://www.tddft.org/programs/octopus/wiki/index.php/Libxc_functionals or the list in http://bigdft.org/Wiki/index.php?title=XC_codes.<br>
> Meta-GGA functionals are not implemented yet.<br>
> Reference :  Miguel A. L. Marques, Micael J. T. Oliveira, and Tobias Burnus, *Comput. Phys. Commun.* **183**, 2272 (2012) [default]


### **XCFunctional**
> #### **Description**
> This is a keyword to specify exchange-correlation functional. See the list in <a href="http://bigdft.org/Wiki/index.php?title=XC_codes">BigDFT website</a> or <a href="http://www.tddft.org/programs/libxc/functionals/">libxc website</a>.<br>
> We accept the Libxc functional codes in string or integer, such as LDA_X (or 1), GGA_X_PBE (or 101), etc.<br>
> This option accepts functionals including both exchange and correlations, such as HYB_GGA_XC_B3LYP.<br>
> These functionals have **\_XC\_** in their functional code.
>
> XCFucntional with integer higher than 10000 invokes custom exchange-correlation module. Such values are listed below.<br>
> To use the range-separated hybrid, *GaussianPotential* option should be supplied. Otherwise, the program adds 100% EXX instead.<br>
> CAM-type hybrids, using both global and range-separated exact exchange, are not supported.
> #### **Type**
> int or string
> #### **Possible Options**
> + 10478: Use LC-wPBE(2Gau) exchange-correlation. This option should be used in conjunction of GaussianPotential or ErfNGaussianPreset input. See Song, J.-W., and Hirao, K., *J. Chem. Phys.* **143**, 144112 (2015).

> + Any libxc functional code or id.

### **XFunctional**
> #### **Description**
> This is a keyword to specify exchange-correlation functional. See the list in <a href="http://bigdft.org/Wiki/index.php?title=XC_codes">BigDFT website</a> or <a href="http://www.tddft.org/programs/libxc/functionals/">libxc website</a>. <br>
 We accept the Libxc functional codes in string or integer, such as LDA_X (or 1), GGA_X_PBE (or 101), etc.<br>
 > This option accepts functionals including only exchange, such as GGA_X_PBE.<br>
 > These functionals have **\_X\_** in their functional code.
 >
 > XFucntional with negative integer invokes custom exchange module. Such values are listed below.<br>
 > You may specify scaling factor, which is useful for tuning exact exchange portion. See example below for format.<br>
> #### **Type**
> int or string
> #### **Possible Options**
> + -12: Use KLI exchange.
> + "GGA_X_PBE 0.4": PBE exchange, scaled by 0.4.
> + "101 0.4": PBE exchange, scaled by 0.4.
> + Any libxc functional code or id.

### **CFunctional**
> #### **Description**
> This is a keyword to specify exchange-correlation functional. See the list in <a href="http://bigdft.org/Wiki/index.php?title=XC_codes">BigDFT website</a> or <a href="http://www.tddft.org/programs/libxc/functionals/">libxc website</a>. <br>
 We accept the Libxc functional codes in string or integer, such as LDA_X (or 1), GGA_X_PBE (or 101), etc.<br>
 > This option accepts functionals including only correlation, such as GGA_C_PBE.<br>
 > These functionals have **\_C\_** in their functional code.<br>
 > You may specify scaling factor, in quotation mark with space-separated argument. See *XFunctional*.
> #### **Type**
> int or string
> + Any libxc functional code or id.

### **GaussianPotential**
> #### **Description**
> Sets range-separated exact exchange kernel using Gaussians. This option should be supplied if range-separated hybrids.<br>
> Input should be a string, consisting two real numbers separated by space.<br>
> First number is Gaussian coefficient, and second number is Gaussian exponent. See examples below.<br>
> Multiple arguments are summed to represent short or long-range exact exchange.<br>
> #### **Type**
> string
> #### **Possible Options**
> + "0.27 0.075": Sets EXX potential to 0.27 exp(-0.075 r<sup>2</sup>).
> + "0.18 0.006": Sets EXX potential to 0.18 exp(-0.006 r<sup>2</sup>).
> + Note: Using both "0.27 0.075" and "0.18 0.006" is the setting used in Song, J.-W., and Hirao, K., *J. Chem. Phys.* **143**, 144112 (2015) to mimic erf(0.4r)/r potential.

### **ErfNGaussianPreset**
> #### **Description**
> Expand erf(wr)/r type exact exchange into N Gaussians. This option sets *GaussianPotential* from range-separation parameter for error function coefficient and the number of Gaussians.<br>
> Supply a string containing range-separation parameter and the number of Gaussians, separated by space.<br>
> Currently 2 Gaussians are supported only.
> #### **Type**
> string
> #### **Possible Options**
> + "0.4 2": Expand erf(0.4r)/r into 2 Gaussians and use it as range-separated exact exchange kernel.<br>
>           This option sets two *GaussianPotential* input, with "0.27 0.075" and "0.18 0.006".

### **CalGradientUsingDensity**
> #### **Description**
> For the GGA and meta-GGA exchange-correlation functionals, the gradient information of the electron density is needed. The gradient of the density can be obtained numerically in a two different way. One is the take derivative action to the density expression.
>
> ∇ρ(r)
>
> The other way is utilizing the derivative information of each orbitals.
>
> ∇ρ(r)=∇|φ₁|²+∇|φ₂|²+...
>       =2φ₁∇φ₁+2φ₂∇φ₂+...
>
> From the numerical reasons, the second approach is more accurate than first one.
> #### **Type**
> integer
> #### **Possible Options**
> + 0: The density gradient is calculated from the gradients of orbitals. [default]
> + Otherwise: The density gradient is calculated from the density information.
## CIS

This section controls parameters related to CIS

(Phys. Chem. Chem. Phys. 2015, 17(47), 31434-31443. Jaewook Kim+, Kwangwoo Hong+, Sunghwan Choi, Sang-Yeon Hwang, and Woo Youn Kim*)

#### **Subsection List**
  +  Diagonalize 

### **NumberOfEigenvalues**
> #### **Description**
> Number of eigenvalues which will be obtained from diagonalization of CIS matrix. 
> #### **Type**
> int
> #### **Possible Options**
> + key: 10 # if you want to calculate the 10 lowest states.

### **ConvergenceTolerance**
> #### **Description**
> Convergence criteria of Davidson digonalization. Relevant only if Davidson diagonalizer is used.<br>
> Defaults to 1.0E-6.
> #### **Type**
> float

## **Occupation**

This section governes occupation of electron. This section is defined in Scf or Guess section.

### **VarName**
> #### **Description**
> Description blah blah
> #### **Type**
> string / int / (positive) float
> #### **Possible Options**
> + key: description

### **OccupationMethod**
> #### **Description**
> Controls how to occupy electrons.
> #### **Type**
> string 
> #### **Possible Options**
> + ZeroTemp: This option will occupy electrons from orbitals with low energy.
> + Fermi: Fermi-Dirac distribution is applied.
> + Input: Get occupation from Occupation_List[_Alpha/_Beta]

### **Temperature**
> #### **Description**
> This parameter works only when *OccupationMethod* is set to Fermi. Excit some electrons near the chemical potential to make occupation number to follow Fermi-Dirac distribution.
> Note that you need to set OccupationSize to sufficiently large value.
>
> 0.1 eV (=1160.4518 K) [default]
> #### **Type**
> (non-negative) float


### **OccupationList**
> #### **Description**
> Option for OccupationMethod Input. Should supply two inputs for spin polarized calculations. \ 
 Example: OccupationList \"2.0 x4 1.5 0.5 0.0\"
> #### **Type**
> string

### **OccupationSize**
> #### **Description**
> Sets the maximum size for the electron occupation.
> For OccupationMethod is ZeroTemp, if this variable is not set or set to small values, it will fallback to the minimum occupation suitable for given NumElectrons and NumberOfEigenvalues.
> For OccupationMethod is Fermi, you need to set this value to sufficiently large value.
> For OccupationMethod is Input, this variable is ignored and can be controlled by OccupationList.
> #### **Type**
> int

### **NumberOfVirtualStates**
> #### **Description**
> Adds some unoccupied orbitals to the occupation when it is automatically determined, so it can be calculated.
> This option is ignored when OccupationSize is provided.
>
> 8 [default] when Fermi
> #### **Type**
> int
## Mixing
 This section governes parameters related to mixing methods. 
 Generates guess Hamiltonian for the next iteration from the output of previous iterations.

### **MixingType**
> #### **Description**
> Decides the field that is to be mixed.
> See the following pdf file (page 8) for the difference of thd density mixing and potential mixing.
> https://www.abinit.org/sites/default/files/oldsites/workshop_07/Presentations/Liege_Torrent2_SC_Mixing.pdf
> #### **Type**
> string
> #### **Possible Options**
>+ Density:\
> Mix density to generate guess Hamiltonian for the next step.\
> For PAW calculations, this mixes the PAW atom-centered density matrix.\
> We recommend this option for the PAW calculations.
>+ Potential:\
> Mix potential to generate guess Hamiltonian for the next step.\
> We recommend this option for the KLI calculations. [default]


### **MixingMethod**
> #### **Description**
> Decides the mixing algorithm.
> #### **Type**
> int
> #### **Possible Options**
>+ 0: Linear mixing. \
> $f^{n+1}_{in} = \alpha f^{n}_{out} + (1-\alpha) f^n_{in}$.
>+ 1 Pulay mixing (DIIS).\
> See Kresse, Phys. Rev. B 54, 11169 (1996), or http://www.hector.ac.uk/cse/distributedcse/reports/conquest/conquest/node3.html. [default]
>+ 2 Broyden mixing. \
> See Phys. Rev. C 78, 014318 (2008), section II.B, or Phys. Rev. B 38, 12807 (1998). Not recommended for PAW.


### **MixingParameter**
> #### **Description**
> Decides the mixing coefficient for the linear mixing.\
> The smaller MixingParameter, the more conservertive guess for the next iteration would be generated. \
> For the case of MixingParameter=0.0, the guess for the next iteration will be the same as the current iteration. \
> Setting this variable close to 1 will converge the calculation fast if the initial state is far from the converged state.
> However, the calculation may oscillate near the converged state, if this variable is too high. \
> Note that setting this value outside of 0-1 can work, but it slowdown convergence. \
> Defaults to 0.3.
> #### **Type** 
> float


### **PulayMixingParameter**
> #### **Description**
> Decides the mixing coefficient for the Pulay mixing.\
> This option works pretty much the same as MixingParameter. \
> Only difference is that this works for Pulay mixing.
> Defaults to MixingParameter.
> #### **Type** 
> float


### **BroydenMixingParameter**
> #### **Description**
> Decides the mixing coefficient for the Broyden mixing.\
> This option works pretty much the same as MixingParameter. \
> Only difference is that this works for Broyden mixing.
> Defaults to MixingParameter.
> #### **Type** 
> (positive) float


### **MixingStart**
> #### **Description**
> Irrelevant to the linear mixing algorithm. \
> Pulay, or Broyden mixing will start when the iterations are performed for MixingStart times. \
> Mixing algorithm is linear before this value. \
> If this variable is smaller than MixingHistory, it will be adjusted to corresponding MixingHistory. \
> Defaults to 6.
> #### **Type** 
> (positive) int

>
### **MixingHistory**
> #### **Description**
> Irrelevant to the linear mixing algorithm. \
> Pulay, or Broyden mixing will use this amount of previous iteration informations. \
> Defaults to 4.
> #### **Type** 
> (positive) int


### **UpdatePawMatrix**
> #### **Description**
> Controls whether PAW atomcenter density matrix mixing. See PAW_Density mixing also.\
> Only works with density or potential mixing.
>+ 0
> Do not update PAW atomcenter density matrix when mixing is done. Default for potential mixing. Only option for non-density non-potential mixing.
> #### **Type** 
> int
> #### **Possible Options**
>+ 1
> Update PAW atomcenter density matrix with same coefficient used to update mixing field. Default for density mixing.\
> In order to use same mixing scheme but different coefficient, use PAW_Density as mixing field.

>
## **Pseudopotential**
 This section set type of pseudopotential/PAW method.

### **Pseudopotential**
>#### **Description**
> Decides Pseudopotential/PAW method that will be used.
>#### **Type**
> int
> #### **Possible Options**
>+ 1: Pseudopotential using KB projector. [default]
>+ 3: PAW method.

### **Format**
> #### **Description**
> Decides the input format for the Pseudopotential/PAW method.
> #### **Type**
> string
> #### **Possible Options**
> +xml: PAW-XML format (https://wiki.fysik.dtu.dk/gpaw/setups/pawxml.html#pawxml).\
> This is the only option for PAW method, and exclusive to PAW method.
> + upf: UPF format (http://www.quantum-espresso.org/pseudopotentials/unified-pseudopotential-format).
> + hgh: HGH format (PRB, 58, 3641 (1998)). [default]

### **HghFormat**
> #### **Description**
> The type of HGH pseudopotential for calculations.
> #### **Type**
> string
> #### **Possible Options**
> + HGH: HGH type of pseudopotential (PRB, 58, 3641 (1998)).
> + Internal: Use hard-coded values, from same reference. [default]

### **PSFilenames**
> #### **Description**
> Path to the input pseudopotential/PAW data file. \
> PSFilenames should be given for all atoms to be calculated, and should be given in the order of increasing atomic number. \
> This variable might be inferred from PSFilePath and PSFileSuffix, if those two variables exist and PSFilenames are absent
> #### **Type**
> string

### **PSFilePath**
> #### **Description**
> Path to the directory that contains pseudopotential/PAW data file. \
> This is used to infer PSFilenames. \
> [PSFilenames] = [PSFilePath]/[Atom Symbol][PSFileSuffix]
> #### **Type**
> string

### **PSFileSuffix**
> #### **Description**
> This is used in conjunction to PSFilePath to infer PSFilenames. \
> [PSFilenames] = [PSFilePath]/[Atom Symbol][PSFileSuffix]
> #### **Type**
> string

### **UsingDoubleGrid**
> #### **Description**
> Optional argument. Exclusive to UPF and PAW pseudopotentials. \
> Use supersampling methods for accurate nuclear potential calculations.
>
> Overrides UsingFiltering.
> #### **Type**
> int
> #### **Possible Options**
> + 0: Do not use supersampling method. [default]
> + 1: Use supersampling method.

### **UsingFiltering**
> #### **Description**
> Optional argument. Exclusive to UPF pseudopotentials. \
> Use Fourier filtering method for anti-aliasing of nuclear potential representation.
> You may rafer DOI [10.1063/1.4942925] for detailed explanation.
>
> Overrided by UsingDoubleGrid.
> #### **Type**
> int
> #### **Possible Options**
> + 0: Do not use Fourier filtering method [default]
> + 1: Use Fourier filtering method.

### **NonlocalRmax**
> #### **Description**
> Optional argument. Used with both supersampling and Fourier filtering.<br>
> Decides the outer sampling region.<br>
> r_out = [NonlocalRmax] \times r_in
> Default value: 1.0
> #### **Type**
> float

### **FilterType**
> #### ** Description**
> Optional argument. Exclusive to PAW. Ignored if [FineProj] <= 1.<br>
> Decides the supersampling filter type.
> #### **Type**
> string
> #### **Possible Options**
> + Sinc: Use sinc function as supersampling filter.[default]
> + Lagrange: Use Lagrange function as supersampling filter. Currently, this option does not work with PAW method.

### **GammaLocal**
> #### ** Description**
> Optional argument. Exclusive to filtering.
> Decides cutoff radius for Bessel transformation of local potential from real space to reciprocal space
> Default value: 2.0
> cutoff_r = r_short * GammaLocal
> #### **Type**
> double

### **GammaNonlocal**
> #### ** Description**
> Optional argument. Exclusive to filtering.
> Decides cutoff radius for Bessel transformation of nonlocal potential from real space to reciprocal space
> Default value: 2.0
> cutoff_r = r_short * GammaNonlocal
> #### **Type**
> double

### **AlphaLocal**
> #### ** Description**
> Optional argument. Exclusive to filtering.
> Decides cutoff frequency for Bessel transformation of local potential from reciprocal space to real space
> Default value: 1.1
> q_cut = q_max / AlphaLocal
> #### **Type**
> double

### **AlphaNonlocal**
> #### ** Description**
> Optional argument. Exclusive to filtering.
> Decides cutoff frequency for Bessel transformation of nonlocal potential from reciprocal space to real space
> Default value: 1.1
> q_cut = q_max / AlphaNonlocal
> #### **Type**
> double

### **Eta**
> #### ** Description**
> Optional argument. Exclusive to filtering.
> Decides constant k which determine mask function exp(-k*x^2).
> Default value: q_cut*cutoff_r
> #### **Type**
> double

### **FineDimension**
> #### **Description**
> Optional argument. Exclusive to UPF and PAW pseudopotentials when UsingDoubleGrid is 1. \
> Controls upsampling rate.
> #### **Type**
> int

### **Lmax**
> #### ** Description**
> Optional argument. Exclusive to PAW.<br>
> Decides the maximum angular momentum to expand PAW partical wave.<br>
> Negative value sets the program includes all angular momentum values.<br>
> Changing this value critically affects PAW Hamiltonian correction calculation time, and may affects the total energy by order of 0.1 mHa for small molecule.
> Default value: -1
> #### **Type**
> int

### **OptimizeVectors**
> #### ** Description**
> Optional argument. Exclusive to KB projector. \
> Decides whether the local pseudopotential and nonlocal projector functions to be deleted after the Hamiltonian construction. \
> This makes less memory, but will consume additional considerable time if supersampling enabled and calculating force.  
> #### **Type**
> int
> #### **Possible Options**
>+ 0: Do not delete pseudopotential vectors. [Defalut value]
>+ 1: Delete pseudopotential vectors.

### **OverlapCutoff**
> #### ** Description**
> Controls PAW overlap matrix cutoff. Default is 1E-7, corresponds to ~1E-6 Ha energy difference. \
> Only used when PAW is selected. \
> Recommend to set this value to 0.1 $\times$ [Energy Accuracy].
> #### **Type**
> float

### **HamiltonianCutoff**
> #### ** Description**
> Controls PAW overlap matrix cutoff. Default is 1E-7, corresponds to ~1E-6 Ha energy difference. \
> Only used when PAW is selected. \
> Recommend to set this value to 0.1 $\times$ [Energy Accuracy].
> #### **Type**
> float

### **OccupancyOutput**
> #### ** Description**
> Controls PAW partial wave occupancy output, defined as the square of the dot product between orbital and PAW proejction operator.
> #### **Type**
> int
> #### **Possible Option**
> + 0: Does not print PAW partial wave occupancy at all.
> + 1: Print the sum of PAW partial wave occupancy over all orbitals.
> + 2: Print the PAW partial wave occupancy for all orbitals and atomic partial wave combinations. [default]

# How to use
1. Add ACE_reader.py to $ASE_PKG_ROOT/io and ACE_calc.py to $ASE_PKG_ROOT/calculator
   Usually $ASE_PKG_ROOT is $PYTHONPATH/site-packages/ase
2. Import and use with
    ```
    from ase.calculators.ACE_cal import ACE
    atom = ase.io.read(...)
    atom.calc = ACE(command = "mpirun -np 20 ./ace PREFIX_opt.inp > PREFIX_opt.log", ACEtemplate = "template", BasicInformation = {"ADDITIONAL": "INPUTS"}, Scf = {"SCF": "INPUTS"})
    ```
