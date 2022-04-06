C++ Implementation of SISSO with python bindings
===

## Overview
This package provides a C++ implementation of SISSO with built in Python bindings for an efficient python interface.
Future work will expand the python interface to include more postporcessing analysis tools.

For a more detailed explanation please visit our documentation [here](https://sissopp_developers.gitlab.io/sissopp/)

## Installation
The package uses a CMake build system, and compatible all versions of the C++ standard library after C++ 14.
You can access the code [here](https://gitlab.mpcdf.mpg.de/tpurcell/cpp_sisso/-/archive/master/cpp_sisso-master.tar.gz)

### Prerequisites
To install `SISSO++` the following packages are needed:

- CMake version 3.10 and up
- A C++ compiler (compatible with C++ 14 and later, e.g. gcc 5.0+ or icpc 17.0+)
- BLAS/LAPACK
- MPI

Additionally the following packages needed by SISSO++ will be installed (if they are not installed already/if they cannot be found in $PATH)

  - [Boost](https://www.boost.org) (mpi, serialization, system, filesystem, and python libraries)
  - [GTest](https://github.com/google/googletest)
  - [Coin-Clp](https://github.com/coin-or/Clp)
  - [NLopt](https://github.com/stevengj/nlopt)
  - [{fmt}](https://fmt.dev/latest/index.html) (Used for the C++ 20 [std::format](https://en.cppreference.com/w/cpp/utility/format/format) library)

To build and use the optional python bindings the following are also needed:

- [Python 3.6 or greater](https://www.python.org/)
- [numpy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [scipy](https://www.scipy.org/)
- [seaborn](https://seaborn.pydata.org/)
- [sklearn](https://scikit-learn.org/stable/index.html)
- [toml](https://pypi.org/project/toml/)

The setup of the python environment can be done using anaconda with

```bash
conda create -n sissopp_env python=3.9 numpy pandas scipy seaborn scikit-learn toml
```

### Installing `SISSO++`

`SISSO++` is installed using a cmake build system, with sample configuration files located in `cmake/toolchains/`
For example, here is `initial_config.cmake` file used to construct `SISSO++` and the python bindings using the gnu compiler.

```
###############
# Basic Flags #
###############
set(CMAKE_CXX_COMPILER g++ CACHE STRING "")
set(CMAKE_C_COMPILER gcc CACHE STRING "")
set(CMAKE_CXX_FLAGS "-O3 -march=native" CACHE STRING "")
set(CMAKE_C_FLAGS "-O3 -march=native" CACHE STRING "")

#################
# Feature Flags #
#################
set(BUILD_PYTHON ON CACHE BOOL "")
set(BUILD_PARAMS ON CACHE BOOL "")
```

Because we want to build with the python bindings in this example and assuming there is no preexisting python environment, we need to first create/activate it.
For this example we will use `conda`, but standard python installations or virtual environments are also possible.

```bash
conda create -n sissopp_env python=3.9 numpy pandas scipy seaborn scikit-learn toml
conda activate sissopp_env
```

Note if you are using a python environment with a local MKL installation then make sure the versions of all accessible MKL libraries are the same.

Now we can install `SISSO++` using `initial_config.cmake` and the following commands (this assumes gnu compiler and MKL are used, if you are using a different compiler/BLAS library change the flags to the relevant directories)

```bash
export MKLROOT=/path/to/mkl/
export BOOST_ROOT=/path/to/boost

cd ~/sissopp/
mkdir build/;
cd build/;

cmake -C initial_config.cmake ../
make
make install
```

Once all the commands are run `SISSO++` should be in the `~/SISSO++/main directory/bin/` directory.

#### Installing the Python Bindings Without Administrative Privileges

To install the python bindings on a machine where you do not have write privilege to the default python install directory (typical on most HPC systems), you must set the `PYTHON_INSTDIR` to a directory where you do have write access.
This can be done by modifying the `camke` command to:

```bash
cmake -C initial_config.cmake -DPYTHON_INSTDIR=/path/to/python/install/directory/ ../
```

A standard local python installation directory for pip and conda is `$HOME/.local/lib/python3.X/site-packages/` where X is the minor version of python.
It is important that if you do set this variable then that directory is also inside your `PYTHONPATH` envrionment variable. This can be updated with

```bash
export PYTHONPATH=$PYTHONPATH:/path/to/python/install/directory/
```

If you are using anaconda, then this can be avoided by creating a new conda environment as detailed above.

You will need to set this variable and recompile the code (remove all build files first) if you see this error

```bash

CMake Error at src/cmake_install.cmake:114 (file):
  file cannot create directory:
  ${PTYHON_BASE_DIR}/lib/python3.X/site-packages/sissopp.
  Maybe need administrative privileges.
Call Stack (most recent call first):
  cmake_install.cmake:42 (include)" 
```

#### Install the Binary Without the Python Bindings

To install only the `SISSO++` executable repeat the same commands as above but set `USE_PYTHON` in `initial_config.cmake` to `OFF`.

## Running the Code

### Input Files

To see a sample of the input files look in `~/sisso++/main directory/tests/exec_test`.
In this directory there are multiple subdirectories for different types of calculations, but the `default/` directory would be the most common application.

To use the code two files are necessary: `sisso.json` and `data.csv`.
`data.csv` stores all of the data for the calculation in a `csv` file.
The first row in the file corresponds to the feature meta data with the following format `expression (Unit)`.
For example if one of the primary features used in the set is the lattice constant of a material the header would be `lat_param (AA)`.
The first column of the file are sample labels for all of the other rows, and is used to set the sample ids in the output files.

The input parameters are stored in `sisso.json`, here is a list of all possible variables that can be set in `sisso.json`

#### `data.csv`
The data file contains all relevant data and metadata to describe the individual features and samples.
The first row of the file corresponds to the features metadata and has the following format `expression (Unit)` or `expression`.
For the cases where no `(Unit)` is included in the header then the feature is considered to be dimensionless.
For example if one of the primary features used in the set is the lattice constant of a material the header would be `lat_param (AA)`, but the number of species in the material would be `n_species` because it is a dimensionless number.


The first column provide the labels for each sample in the data file, and is used to set the sample ids in the output files.
In the simplest case, this can be just a running index.
The data describing the property vector is defined in the column with an `expression` matching the `property_key` filed in the `sisso.json` file, and will not be included in the feature space.
Additionally, an optional `Task` column whose header matches the `task_key` field in the sisso.json file can also be included in the data file.
This column maps each sample to a respective task with a label defined in the task column.
Below in a minimal example of the data file used to learn a model for a materials volume.

```
material, Structure_Type, Volume (AA^3), lat_param (AA)
C, diamond, 45.64, 3.57
Si, diamond, 163.55, 5.47
Ge, diamond, 191.39, 5.76
Sn, diamond, 293.58, 6.65
Pb, diamond, 353.84, 7.07.757
LiF, rock_salt, 67.94, 4.08
NaF, rock_salt, 103.39, 4.69
KF, rock_salt, 159.00, 5.42
RbF, rock_salt, 189.01, 5.74
CsF, rock_salt, 228.33, 6.11
```

#### `sisso.json`

All input parameters that can not be extracted from the data file are defined in the `sisso.json` file.

Here is a complete example of a `sisso.json` file where the property and task keys match those in the above data file example.

```json
{
    "data_file": "data.csv",
    "property_key": "Volume",
    "task_key": "Structure_Type",
    "opset": ["add", "sub", "mult", "div", "sq", "cb", "cbrt", "sqrt"],
    "param_opset": [],
    "calc_type": "regression",
    "desc_dim": 2,
    "n_sis_select": 5,
    "max_rung": 2,
    "max_leaves": 4,
    "n_residual": 1,
    "n_models_store": 1,
    "n_rung_store": 1,
    "n_rung_generate": 0,
    "min_abs_feat_val": 1e-5,
    "max_abs_feat_val": 1e8,
    "leave_out_inds": [0, 5],
    "leave_out_frac": 0.25,
    "fix_intercept": false,
    "max_feat_cross_correlation": 1.0,
    "nlopt_seed": 13,
    "global_param_opt": false,
    "reparam_residual": true
}
```

### Performing the Calculation
Once the input files are made the code can be run using the following command

```
mpiexec -n 2 ~/sisso++/main directory/bin/sisso++ sisso.json
```

which will give the following output for the simple problem defined above

```text
time input_parsing: 0.000721931 s
time to generate feat sapce: 0.00288105 s
Projection time: 0.00304198 s
Time to get best features on rank : 1.09673e-05 s
Complete final combination/selection from all ranks: 0.00282502 s
Time for SIS: 0.00595999 s
Time for l0-norm: 0.00260496 s
Projection time: 0.000118971 s
Time to get best features on rank : 1.38283e-05 s
Complete final combination/selection from all ranks: 0.00240111 s
Time for SIS: 0.00276804 s
Time for l0-norm: 0.000256062 s
Train RMSE: 0.293788 AA^3; Test RMSE: 0.186616 AA^3
c0 + a0 * (lat_param^3)

Train RMSE: 0.0936332 AA^3; Test RMSE: 15.8298 AA^3
c0 + a0 * ((lat_param^3)^2) + a1 * (sqrt(lat_param)^3)

```

### Analyzing the Results

Once the calculations are done, two sets of output files are generated.
Two files that summarize the results from SIS in a computer and human readable manner are stored in: `feature_space/` and every model used as a residual for SIS is stored in `models/`.
The human readable file describing the selected feature space is `feature_space/SIS_summary.txt` which contains the projection score (The Pearson correlation to the target property or model residual).
```
# FEAT_ID     Score                   Feature Expression
0             0.99997909235669924     (lat_param^3)
1             0.999036700010245471    ((lat_param^2)^2)
2             0.998534266139345261    (lat_param^2)
3             0.996929900301868899    (sqrt(lat_param)^3)
4             0.994755117666830335    lat_param
#-----------------------------------------------------------------------
5             0.0318376000648976157   ((lat_param^3)^3)
6             0.00846237838476477863  ((lat_param^3)^2)
7             0.00742498801557322716  cbrt(cbrt(lat_param))
8             0.00715447033658055554  cbrt(sqrt(lat_param))
9             0.00675695980092700429  sqrt(sqrt(lat_param))
#---------------------------------------------------------------------
```
The computer readable file file is `feature_space/selected_features.txt` and contains a the list of selected features represented by an alphanumeric code where the integers are the index of the feature in the primary feature space and strings represent the operators.
The order of each term in these expressions is the same as the order it would appear using postfix (reverse polish) notation.
```
# FEAT_ID     Feature Postfix Expression (RPN)
0             0|cb
1             0|sq|sq
2             0|sq
3             0|sqrt|cb
4             0
#-----------------------------------------------------------------------
5             0|cb|cb
6             0|cb|sq
7             0|cbrt|cbrt
8             0|sqrt|cbrt
9             0|sqrt|sqrt
#-----------------------------------------------------------------------
```

The model output files are split into train/test files sorted by the dimensionality of the model and by the train RMSE.
The model with the lowest RMSE is stored in the lowest numbered file.
For example `train_dim_2_model_0.dat` will have the best 2D model, `train_dim_2_model_1.dat` would have the second best, etc., whereas `train_dim_1_model_0.dat` will have the best 1D model.
Each model file has a large header containing information about the features selected and model generated
```
# c0 + a0 * (lat_param^3)
# Property Label: $Volume$; Unit of the Property: AA^3
# RMSE: 0.293787533962641; Max AE: 0.56084644346538
# Coefficients
# Task       a0                      c0
#  diamond,  1.000735616997855e+00, -1.551085274074442e-01,
#  rock_salt,  9.998140372873336e-01,  6.405707194855371e-02,
# Feature Rung, Units, and Expressions
# 0;  1; AA^3;                                             0|cb; (lat_param^3); $\left(lat_{param}^3\right)$; (lat_param).^3; lat_param
# Number of Samples Per Task
# Task    , n_mats_train
#  diamond, 4
#  rock_salt, 4
```
The first section of the header summarizes the model by providing a string representation of the model, defines the property's label and unit, and summarizes the error of the model.
```
# c0 + a0 * (lat_param^3)
# Property Label: $Volume$; Unit of the Property: AA^3
# RMSE: 0.293787533962641; Max AE: 0.56084644346538
```
Next the linear coefficients (as shown in the first line) for each task is listed.
```
# Coefficients
# Task       a0                      c0
#  diamond,  1.000735616997855e+00, -1.551085274074442e-01,
#  rock_salt,  9.998140372873336e-01,  6.405707194855371e-02,
```
Then a description of each feature in the model is listed, including units and various expressions.
```
# Feature Rung, Units, and Expressions
# 0;  1; AA^3;                                             0|cb; (lat_param^3); $\left(lat_{param}^3\right)$; (lat_param).^3; lat_param
```
Finally information about the number of samples in each task is given
```
# Number of Samples Per Task
# Task    , n_mats_train
#  diamond, 4
#  rock_salt, 4
```

The header for the test data files contain the same information as the training file, with an additional line at the end to list all indexes included in the test set:
```
# Test Indexes: [ 0, 5 ]
```
These indexes can be used to reproduce the results by setting `leave_out_inds` to those listed on this line.

After this header in both file the following data is stored in the file:

```
# Sample ID , Property Value        ,  Property Value (EST) ,  Feature 0 Value
```
With this data, one can plot and analyzed the model, e.g., by using the python binding.



### Using the Python Library
To see how the python interface can be used refer to the [tutorials](https://sissopp_developers.gitlab.io/sissopp/tutorial/2_python.html).
If you get an error about not being able to load MKL libraries, you may have to run `conda install numpy` to get proper linking.
Thomas A. R. Purcell# The Command Line Interface

## Running `SISSO++` from the Command Line
The `SISSO++` binary uses a `sisso.json` input file to define the parameters needed to read the data and set up the SISSO calculation.
As an example here is the `sisso.json` file we will initially use for this system (for the definition of all options refer to the [Quick Start Guide](../quick_start/code_ref):
```json
{
    "data_file": "data.csv",
    "property_key": "E_RS - E_ZB",
    "desc_dim": 4,
    "n_sis_select": 10,
    "max_rung": 2,
    "calc_type": "regression",
    "min_abs_feat_val": 1e-5,
    "max_abs_feat_val": 1e8,
    "n_residual": 10,
    "n_models_store": 1,
    "leave_out_frac": 0.0,
    "leave_out_inds": [],
    "opset": ["add", "sub", "abs_diff", "mult", "div", "inv", "abs", "exp", "log", "sin", "cos", "sq", "cb", "six_pow", "sqrt", "cbrt", "neg_exp"]
}
```
Of these parameters `n_sis_select`,  `n_residual`,  `max_rung`, and `desc_dim` are the hyperparameters that must be optimized for each calculation.
Additionally `property_key` and `task_key` both must be columns headers in the `data_file` (Here we are only using one task so `task_key` is not included).
With this input file and the provided `data.csv` file we are now able to perform SISSO with the following command
```
mpiexec -n 2 sisso++ sisso.json
```
and get the following on screen output
```
time input_parsing: 0.00142002 s
time to generate feat sapce: 1.15551 s
Projection time: 0.137257 s
Time to get best features on rank : 3.69549e-05 s
Complete final combination/selection from all ranks: 0.00072813 s
Time for SIS: 0.184525 s
Time for l0-norm: 0.025141 s
Projection time: 0.157892 s
Time to get best features on rank : 4.91142e-05 s
Complete final combination/selection from all ranks: 0.00104713 s
Time for SIS: 0.209384 s
Time for l0-norm: 0.00287414 s
Projection time: 0.165635 s
Time to get best features on rank : 9.39369e-05 s
Complete final combination/selection from all ranks: 0.000714064 s
Time for SIS: 0.218518 s
Time for l0-norm: 0.00742316 s
Projection time: 0.166316 s
Time to get best features on rank : 6.10352e-05 s
Complete final combination/selection from all ranks: 0.000627995 s
Time for SIS: 0.218724 s
Time for l0-norm: 0.0982461 s
Train RMSE: 0.125912 eV
c0 + a0 * ((E_HOMO_B / r_p_A) / (r_sigma + r_p_B))

Train RMSE: 0.0931541 eV
c0 + a0 * ((EA_B - IP_A) * (|r_sigma - r_s_B|)) + a1 * ((E_HOMO_B / r_p_A) / (r_sigma + r_p_B))

Train RMSE: 0.0725522 eV
c0 + a0 * ((r_d_B * Z_A) / (|r_p_B - r_s_B|)) + a1 * ((|EA_B - IP_A|) * (|r_sigma - r_s_B|)) + a2 * ((E_HOMO_B / r_p_A) / (r_sigma + r_p_B))

Train RMSE: 0.0588116 eV
c0 + a0 * ((E_LUMO_A / EA_A) / (r_p_B^6)) + a1 * ((|period_B - period_A|) / (r_pi * EA_B)) + a2 * ((EA_B - IP_A) * (|r_sigma - r_s_B|)) + a3 * ((E_HOMO_B / r_p_A) / (r_sigma + r_p_B))
```
The standard output provides information about what step the calculation just finished and how long it took to complete so you can see where a job failed or ran out of time.
If this:
```
Train RMSE: 0.0588116 eV
c0 + a0 * ((E_LUMO_A / EA_A) / (r_p_B^6)) + a1 * ((|period_B - period_A|) / (r_pi * EA_B)) + a2 * ((EA_B - IP_A) * (|r_sigma - r_s_B|)) + a3 * ((E_HOMO_B / r_p_A) / (r_sigma + r_p_B))
```
is not shown at the bottom of standard output then the calculation did not complete successfully.
When all calculations are complete the code prints out a summary of the best 1D, 2D, ..., {desc_dim}D models with their training RMSE/Testing RMSE (Only training if there is no test set provided as in this case).
We also see that, two additional output files are stored in `feature_space/`: `SIS_summary.txt` and `selected_features.txt`.
These files represent a human readable (`SIS_summary.txt`) and computer readable (`selected_features.txt`) summary of the selected feature space from SIS.
Below are reconstructions of both files for this calculation (To see the file click the triangle)

<details>
    <summary>feature_space/SIS_summary.txt</summary>

    # FEAT_ID     Score                   Feature Expression
    0             0.920868624862486329    ((E_HOMO_B / r_p_A) / (r_sigma + r_p_B))
    1             0.919657911026942054    ((|r_pi - r_s_A|) / (r_s_A^3))
    2             0.913575821723830561    ((r_s_B / r_p_A) / (r_pi + r_p_B))
    3             0.913335178071393861    ((IP_B / r_p_A) / (r_sigma + r_d_B))
    4             0.912838504121724292    ((IP_B / r_p_A) / (r_sigma + r_p_B))
    5             0.911915092723797449    ((IP_B / r_p_A) / (r_s_A^2))
    6             0.909420884428965848    ((r_s_B / r_p_A) / (r_p_B * r_p_A))
    7             0.908873529961838345    ((IP_B / r_p_A) / (r_sigma + r_s_B))
    8             0.908165571792444726    ((E_HOMO_B / r_p_A) / (r_s_A^2))
    9             0.90609543772524781     ((E_HOMO_B / r_p_A) / (r_sigma + r_s_B))
    #-----------------------------------------------------------------------
    10            0.477894637470377692    ((|EA_B - IP_A|) * (|r_sigma - r_s_B|))
    11            0.474277869552096332    ((EA_B - IP_A) * (|r_sigma - r_s_B|))
    12            0.45075631366629626     ((|EA_B - IP_A|) * (|r_sigma - r_p_B|))
    13            0.449774820279551457    ((EA_B - IP_A) * (|r_sigma - r_p_B|))
    14            0.429583630959742335    (|(|r_sigma|) - (|r_sigma - r_d_B|)|)
    15            0.420792823945547012    ((E_HOMO_A / r_s_A) * (|r_sigma - r_s_B|))
    16            0.413613786609261791    ((IP_A / r_s_A) * (|r_sigma - r_s_B|))
    17            0.412833543057880004    ((E_HOMO_B / r_d_A) / (|r_sigma - r_d_A|))
    18            0.40708791994290916     ((IP_B / r_d_A) / (|r_sigma - r_d_A|))
    19            0.406210357359531482    ((E_HOMO_A / r_p_A) * (|r_sigma - r_s_B|))
    #-----------------------------------------------------------------------
    20            0.39445559217196724     ((r_sigma * E_HOMO_A) / (r_pi * EA_B))
    21            0.394109435179661072    ((|period_B - period_A|) / (r_pi * EA_B))
    22            0.38222896050736499     ((|Z_B - Z_A|) / (r_pi * EA_B))
    23            0.376740768063882125    ((r_sigma / EA_B) / (|r_p_B - r_s_B|))
    24            0.372118309722369145    ((|Z_B - Z_A|) / (r_pi * Z_B))
    25            0.371428608061153243    ((r_d_B * Z_A) / (|r_p_B - r_s_B|))
    26            0.371280945046756461    ((r_sigma * IP_A) / (r_pi * EA_B))
    27            0.367930264720841615    ((r_sigma * IP_A) / (r_p_A * EA_B))
    28            0.365697893302494692    ((Z_A / r_d_A) / (Z_B^6))
    29            0.364460897568853914    ((|period_B - period_A|) / (EA_B * Z_B))
    #-----------------------------------------------------------------------
    30            0.290108702279498698    ((EA_A / r_pi) / (r_p_B^6))
    31            0.288024533781257031    ((EA_A / r_pi) / (r_s_B^6))
    32            0.287090141076753791    ((E_LUMO_A / EA_A) / (r_p_B^6))
    33            0.277661937205204767    ((E_LUMO_A / EA_A) / (r_s_B^6))
    34            0.273556648033503158    ((E_LUMO_A^6) / (r_p_B^6))
    35            0.269806696772049381    ((E_HOMO_B^6) * (EA_A / r_pi))
    36            0.267089521281145048    ((E_LUMO_A^6) / (r_s_B^6))
    37            0.265175400489099211    ((E_LUMO_A / period_B)^6)
    38            0.262777418218664849    ((E_LUMO_A^6) / (r_p_B^3))
    39            0.253659279222423484    ((E_LUMO_A / r_p_B) * (E_LUMO_B * E_LUMO_A))
    #-----------------------------------------------------------------------
</details>
This file contains the index of the selected feature space, a projection score, and a string representation of the feature.
For regression problems the score represents the Pearson correlation between the feature and target property (all feature above the first dashed line) or the highest Pearson correlation between the feature and the residual of the best `n_residual` models of the previous dimension.

<details>
    <summary>feature_space/selected_features.txt</summary>

    # FEAT_ID     Feature Postfix Expression (RPN)
    0             9|14|div|18|15|add|div
    1             19|12|abd|12|cb|div
    2             13|14|div|19|15|add|div
    3             5|14|div|18|17|add|div
    4             5|14|div|18|15|add|div
    5             5|14|div|12|sq|div
    6             13|14|div|15|14|mult|div
    7             5|14|div|18|13|add|div
    8             9|14|div|12|sq|div
    9             9|14|div|18|13|add|div
    #-----------------------------------------------------------------------
    10            7|4|abd|18|13|abd|mult
    11            7|4|sub|18|13|abd|mult
    12            7|4|abd|18|15|abd|mult
    13            7|4|sub|18|15|abd|mult
    14            18|abs|18|17|abd|abd
    15            8|12|div|18|13|abd|mult
    16            4|12|div|18|13|abd|mult
    17            9|16|div|18|16|abd|div
    18            5|16|div|18|16|abd|div
    19            8|14|div|18|13|abd|mult
    #-----------------------------------------------------------------------
    20            18|8|mult|19|7|mult|div
    21            3|2|abd|19|7|mult|div
    22            1|0|abd|19|7|mult|div
    23            18|7|div|15|13|abd|div
    24            1|0|abd|19|1|mult|div
    25            17|0|mult|15|13|abd|div
    26            18|4|mult|19|7|mult|div
    27            18|4|mult|14|7|mult|div
    28            0|16|div|1|sp|div
    29            3|2|abd|7|1|mult|div
    #-----------------------------------------------------------------------
    30            6|19|div|15|sp|div
    31            6|19|div|13|sp|div
    32            10|6|div|15|sp|div
    33            10|6|div|13|sp|div
    34            10|sp|15|sp|div
    35            9|sp|6|19|div|mult
    36            10|sp|13|sp|div
    37            10|3|div|sp
    38            10|sp|15|cb|div
    39            10|15|div|11|10|mult|mult
    #-----------------------------------------------------------------------

</details>
This files is a computer readable file used to reconstruct the selected feature space.
In these files each feature is displayed an alphanumeric string where the integers represent an index of the primary feature space, and the strings represent operations.
The order of each term matches the order of terms if the equation is written in postfix (reverse polish) notation.
In both files the change in rung is represented by the commented out dashed (--) line.

The `models/` directory is used to store the output files representing the models for each dimension:
```
ls models/
train_dim_1_model_0.dat  train_dim_2_model_0.dat  train_dim_3_model_0.dat  train_dim_4_model_0.dat
```
Each of these files represents one of the {`n_models_store`} model stored for each dimension, and can be used to reconstruct the models within python.
The file has a header that provides metadata associated with the selected features, coefficients, modeled property, and the task sizes for the calculations.
The first six lines of the header are the most important because it defines what the model is, what the error is, and the coefficients for each task.
After the header the value of the property, estimated property, and feature value for each sample is listed with the same label used in `data.csv`.
For a line by line description of the header refer to the [quick-start guide](../quick_start/code_ref.md).
An example of these files is provided here:

<details>
    <summary>models/train_dim_2_model_0.dat</summary>

    # c0 + a0 * ((EA_B - IP_A) * (|r_sigma - r_s_B|)) + a1 * ((E_HOMO_B / r_p_A) / (r_sigma + r_p_B))
    # Property Label: $E_{RS} - E_{ZB}$; Unit of the Property: eV
    # RMSE: 0.0931540779192557; Max AE: 0.356632500670745
    # Coefficients
    # Task   a0                      a1                      c0
    # all , -3.208426383962958e-02, -2.400330764129759e-01, -2.869432663392750e-01,
    # Feature Rung, Units, and Expressions
    # 0;  2; AA * eV_IP;                                       7|4|sub|18|13|abd|mult; ((EA_B - IP_A) * (|r_sigma - r_s_B|)); $\left(\left(EA_{B} - IP_{A}\right) \left(\left|r_{sigma} - r_{s, B}\right|\right)\right)$; ((EA_B - IP_A) .* abs(r_sigma - r_s_B)); EA_B,IP_A,r_sigma,r_s_B
    # 1;  2; AA^-2 * eV;                                       9|14|div|18|15|add|div; ((E_HOMO_B / r_p_A) / (r_sigma + r_p_B)); $\left(\frac{ \left(\frac{ E_{HOMO, B} }{ r_{p, A} } \right) }{ \left(r_{sigma} + r_{p, B}\right) } \right)$; ((E_HOMO_B ./ r_p_A) ./ (r_sigma + r_p_B)); E_HOMO_B,r_p_A,r_sigma,r_p_B
    # Number of Samples Per Task
    # Task, n_mats_train
    # all , 82

    # Sample ID , Property Value        ,  Property Value (EST) ,  Feature 0 Value      ,  Feature 1 Value
    AgBr        , -3.003341671137600e-02,  1.639017458701624e-02,  3.541416008482608e+00, -1.737082125260779e+00
    AgCl        , -4.279727820539900e-02,  1.221788989673778e-02,  4.414283985050924e+00, -1.836372781878862e+00
    AgF         , -1.537576731789160e-01, -1.416891410786941e-02,  7.607045749342502e+00, -2.153206644755162e+00
    AgI         ,  3.692541964119300e-02,  6.662284759163006e-02,  1.499718119621954e+00, -1.673450475112165e+00
    AlAs        ,  2.132618491086760e-01,  2.460156857673837e-01,  1.024737940396364e+00, -2.357328927365386e+00
    AlN         ,  7.294907316963900e-02,  2.456878076268038e-01,  3.482569932972883e+00, -2.684491554934336e+00
    AlP         ,  2.189583414756270e-01,  2.801436096644385e-01,  5.790749046839015e-01, -2.439939014991684e+00
    AlSb        ,  1.568687339604370e-01,  2.007443516059665e-01,  2.950349986552500e+00, -2.426113242539595e+00
    AsGa        ,  2.742777724197370e-01,  3.151675906408351e-01,  1.671179867660379e+00, -2.731829473573980e+00
    AsB         ,  8.749781837650520e-01,  7.881634596496340e-01,  3.810479887758816e+00, -4.988325717257844e+00
    BN          ,  1.712080260836270e+00,  1.606089893577266e+00,  3.161246801594280e-01, -7.928806379528248e+00
    BP          ,  1.019225161195210e+00,  1.010436354881370e+00,  4.200899478533967e+00, -5.966520988922634e+00
    BSb         ,  5.808491143689030e-01,  4.226792460264598e-01,  2.600752635485375e+00, -3.303985258847382e+00
    BaO         , -9.299855386780100e-02, -3.688421176160839e-01,  8.608614320176068e+00, -8.094809464135576e-01
    BaS         , -3.197624294261910e-01, -3.363542918017661e-01,  6.543460373759491e+00, -6.687873437791405e-01
    BaSe        , -3.434451340872330e-01, -3.321656670638599e-01,  6.165281404183006e+00, -6.356878675764703e-01
    BaTe        , -3.753868096682710e-01, -3.026651587599424e-01,  5.015472606607728e+00, -6.048993586215626e-01
    BeO         ,  6.918375772329460e-01,  5.153912242420119e-01,  6.066290623904059e+00, -4.153456575173075e+00
    BeS         ,  5.063276745431720e-01,  6.140186443576880e-01,  2.645797589979559e-01, -3.788855982192108e+00
    BeSe        ,  4.949404427752600e-01,  5.429907890755682e-01,  1.744183600985631e+00, -3.690720110296145e+00
    BeTe        ,  4.685859104938570e-01,  4.516260835242951e-01,  4.959181529624447e+00, -3.739822244098568e+00
    C2          ,  2.628603639133640e+00,  2.783574944051988e+00,  6.386751756964515e+00, -1.364575452594942e+01
    CaO         , -2.652190413191420e-01, -3.238274257878371e-01,  9.342332840592833e+00, -1.095089544385995e+00
    CaS         , -3.691331945374260e-01, -2.680040315596177e-01,  6.270424713694642e+00, -9.170452629688933e-01
    CaSe        , -3.607977344217940e-01, -2.575530743234562e-01,  5.625809745654779e+00, -8.744218061257983e-01
    CaTe        , -3.504562790767520e-01, -2.135956642585606e-01,  3.987719873696632e+00, -8.385955037322098e-01
    CdO         , -8.416135802690400e-02, -1.232460775026244e-01,  1.065231004987821e+01, -2.105829423728353e+00
    CdS         ,  7.267279591178499e-02,  1.431422536464881e-02,  4.311359919493231e+00, -1.831348860072292e+00
    CdSe        ,  8.357194908603600e-02,  4.401882316506024e-02,  2.868768109826145e+00, -1.762275469514607e+00
    CdTe        ,  1.145395321946130e-01,  1.171042634651223e-01,  3.457696716299031e-01, -1.729517037320575e+00
    BrCs        , -1.558673029940110e-01, -1.898926962890557e-01,  8.647554645764165e-01, -5.199100657183604e-01
    ClCs        , -1.503461574466200e-01, -1.571430901173148e-01,  1.238997024850578e-01, -5.573207199827011e-01
    CsF         , -1.082633186742900e-01, -8.428877325507279e-02, -1.184139604981316e+00, -6.859981467376031e-01
    CsI         , -1.623874744982460e-01, -2.139647446584292e-01,  1.354924677944221e+00, -4.851426489675994e-01
    BrCu        ,  1.524426397882050e-01,  1.751625662935415e-01,  2.324749827375701e+00, -2.235915680770981e+00
    ClCu        ,  1.562587131920740e-01,  1.703409819257667e-01,  3.357679763615902e+00, -2.353896138537911e+00
    CuF         , -1.702227234272900e-02,  1.432152572174888e-01,  6.954856486036995e+00, -2.721708123665148e+00
    CuI         ,  2.046745832631130e-01,  2.336424208281593e-01,  4.875295043074224e-02, -2.175324740635570e+00
    GaN         ,  4.334452390939990e-01,  3.544948246079955e-01,  2.884011194974658e+00, -3.057784693724841e+00
    GaP         ,  3.487517977519020e-01,  3.520992885010393e-01,  1.208441824761596e+00, -2.823837994787984e+00
    GaSb        ,  1.546252850966990e-01,  2.794927768939737e-01,  3.614065011465625e+00, -2.842902606558420e+00
    Ge2         ,  2.008525260607710e-01,  2.394501685804619e-01,  6.088560028849320e+00, -3.006837274572150e+00
    CGe         ,  8.114428802000480e-01,  4.548103795293029e-01,  1.138082099212689e+00, -3.242337197194732e+00
    GeSi        ,  2.632101701783540e-01,  2.725388874377264e-01,  6.113819992476531e+00, -3.148064336698105e+00
    AsIn        ,  1.340475751931080e-01,  1.801608579254571e-01,  4.068020252879492e-01, -2.000374593993287e+00
    InN         ,  1.537202926992900e-01,  1.448590522634359e-01,  3.816695674128193e+00, -2.309090888179985e+00
    InP         ,  1.791932872292820e-01,  2.105475238127370e-01,  7.234639577287595e-12, -2.072592650924130e+00
    InSb        ,  7.805987301981100e-02,  1.319903471163014e-01,  2.214419977436284e+00, -2.041308871200640e+00
    BrK         , -1.661759641938260e-01, -1.296547837765485e-01,  1.519640837548166e+00, -8.584026968916481e-01
    ClK         , -1.644606802110500e-01, -1.032766020647462e-01,  1.132879978116893e+00, -9.165998606480857e-01
    FK          , -1.464060984981190e-01, -3.718539436956973e-02,  5.397857467140130e-01, -1.112665405432063e+00
    IK          , -1.670391451625620e-01, -1.431459979136605e-01,  1.563489995919397e+00, -8.080581929117355e-01
    BrLi        , -3.274621288437600e-02, -2.060136493312719e-02,  2.019046121317332e+00, -1.379482839678231e+00
    ClLi        , -3.838148269915100e-02, -2.057167653941297e-03,  2.078199280986889e+00, -1.464646447821141e+00
    FLi         , -5.948831686373500e-02,  4.809501062409577e-02,  2.596776386498242e+00, -1.742901194835849e+00
    ILi         , -2.166093634150500e-02, -1.658372988181873e-02,  1.416168070359142e+00, -1.315636374733589e+00
    MgO         , -2.322747243169940e-01, -1.709637107152052e-01,  9.458655848384394e+00, -1.747482354099491e+00
    MgS         , -8.669950498824600e-02, -7.634077484966063e-02,  4.672979827143721e+00, -1.502008033928903e+00
    MgSe        , -5.530180195637500e-02, -5.633708966390212e-02,  3.594547977146389e+00, -1.441195553246026e+00
    MgTe        , -4.591286648065000e-03,  1.388515614782964e-02,  1.127931194217332e+00, -1.404045098915806e+00
    BrNa        , -1.264287278827400e-01, -1.713365131268951e-01,  2.863734255675713e+00, -8.644123624074346e-01
    ClNa        , -1.329919850813890e-01, -1.536717536375947e-01,  2.742537497136659e+00, -9.218054972107542e-01
    FNa         , -1.457881377878040e-01, -1.146505476072913e-01,  2.962752633454122e+00, -1.113806729932952e+00
    INa         , -1.148382221872450e-01, -1.700256955871654e-01,  2.461824196702322e+00, -8.161516351556487e-01
    BrRb        , -1.638205314229710e-01, -2.129364033898979e-01,  1.681775685044998e+00, -5.331156841370923e-01
    ClRb        , -1.605035540778770e-01, -1.837284810776369e-01,  1.056091935643463e+00, -5.711659393468355e-01
    FRb         , -1.355957769847010e-01, -1.206203707726129e-01,  6.544971557289886e-02, -7.016649706208032e-01
    IRb         , -1.672014421201310e-01, -2.313261599637700e-01,  1.992777843630929e+00, -4.980726751118079e-01
    Si2         ,  2.791658215483040e-01,  2.916041179034760e-01,  6.358817979658935e+00, -3.260239754057512e+00
    CSi         ,  6.690237272359810e-01,  4.822282210388984e-01,  1.101648177462507e+00, -3.351692484156472e+00
    Sn2         ,  1.696389919379700e-02,  2.567935768664463e-02,  6.363815657894407e+00, -2.153040623998902e+00
    CSn         ,  4.535379741428190e-01,  1.672793503038666e-01,  3.023496041144674e+00, -2.296472092858216e+00
    GeSn        ,  8.166336023714400e-02,  8.544807638632528e-02,  3.656280114481603e+00, -2.040137159045913e+00
    SiSn        ,  1.351087991060920e-01,  1.054171704626117e-01,  3.690378073553809e+00, -2.127887990332651e+00
    OSr         , -2.203066231741100e-01, -3.724240237555350e-01,  9.409926969679084e+00, -9.016666595490496e-01
    SSr         , -3.684341299303920e-01, -3.249109578073314e-01,  6.787670829981468e+00, -7.491039692724043e-01
    SeSr        , -3.745109517331000e-01, -3.168489728859340e-01,  6.265945830481024e+00, -7.129540446704449e-01
    SrTe        , -3.792947258625650e-01, -2.790737245833743e-01,  4.846463948497374e+00, -6.805927425088493e-01
    OZn         ,  1.019681767684230e-01,  6.602587728853203e-02,  9.268479722491282e+00, -2.709382815714799e+00
    SZn         ,  2.758133256065780e-01,  2.143486874814196e-01,  2.332991532953503e+00, -2.400270322363168e+00
    SeZn        ,  2.631368992806530e-01,  2.463580576975095e-01,  7.384497385908948e-01, -2.320488278555971e+00
    TeZn        ,  2.450012951740060e-01,  1.776248032825628e-01,  2.763715059556858e+00, -2.304848319397327e+00
</details>


## Determining the Ideal Model Complexity with Cross-Validation
While the training error always decreases with descriptor dimensionality for a given application, over-fitting can reduce the general applicability of the models outside of the training set.
In order to determine the optimal dimensionality of a model and optimize the hyperparameters associated with SISSO, we need to perform cross-validation.
The goal of cross-validation is to test how generalizable a given model is with respect to new data.
In practice, we perform cross-validation by randomly splitting the data-set into separate train/test sets and evaluate the performance of the model on the test set.

As an example we will discuss how to perform leave-out 10% using the command line.
To do this we have to modify the `sisso.json` file to automatically leave out a random sample of the training data and use that as a test set by changing `"leave_out_frac": 0.0` to `"leave_out_frac": 0.10`,
i.e. in this case SISSO will ignore 8 materials (10% of all data) during training.
In each run, this 8 materials are chosen randomly, so each SISSO run will
differ from one another.

<details>
    <summary> updated sisso.json file</summary>

    {
        "data_file": "data.csv",
        "property_key": "E_RS - E_ZB",
        "desc_dim": 4,
        "n_sis_select": 10,
        "max_rung": 2,
        "calc_type": "regression",
        "min_abs_feat_val": 1e-5,
        "max_abs_feat_val": 1e8,
        "n_residual": 10,
        "n_models_store": 1,
        "leave_out_frac": 0.10,
        "leave_out_inds": [],
        "opset": ["add", "sub", "abs_diff", "mult", "div", "inv", "abs", "exp", "log", "sin", "cos", "sq", "cb", "six_pow", "sqrt", "cbrt", "neg_exp"]
    }

</details>

Now lets make ten cross validation directories in the working directory and copy the `data.csv` and `sisso.json` into them and run separate calculations for each run.
Note the decision to begin with ten iterations is arbitrary, and not connected to the amount of data excluded from the test set.
```bash
for ii in `seq -f "%03g" 0 9`; do
    mkdir cv_$ii;
    cp sisso.json data.csv cv_$ii;
    cd cv_$ii;
    mpiexec -n 2 sisso++;
    cd ../;
done
```
Each of these directories has the same kind of output files as the non-cross-validation calculations, with the testing and training data defined in separate files in  `cv_$ii/models/`
```
ls cv_00/models/
test_dim_1_model_0.dat  test_dim_3_model_0.dat  train_dim_1_model_0.dat  train_dim_3_model_0.dat
test_dim_2_model_0.dat  test_dim_4_model_0.dat  train_dim_2_model_0.dat  train_dim_4_model_0.dat
```

A full example of the testing set output file is reproduced below:
<details>
    <summary>The test data file cv_0/models/test_dim_2_model_0.dat</summary>

    # c0 + a0 * ((E_HOMO_B / r_p_A) / (r_sigma + r_p_B))
    # Property Label: $E_{RS} - E_{ZB}$; Unit of the Property: eV
    # RMSE: 0.212994478440008; Max AE: 0.442277221520276
    # Coefficients
    # Task   a0                      c0
    # all , -2.346702839867608e-01, -3.917145321667535e-01,
    # Feature Rung, Units, and Expressions
    # 0;  2; AA^-2 * eV;                                       9|14|div|18|15|add|div; ((E_HOMO_B / r_p_A) / (r_sigma + r_p_B)); $\left(\frac{ \left(\frac{ E_{HOMO, B} }{ r_{p, A} } \right) }{ \left(r_{sigma} + r_{p, B}\right) } \right)$; ((E_HOMO_B ./ r_p_A) ./ (r_sigma + r_p_B)); E_HOMO_B,r_p_A,r_sigma,r_p_B
    # Number of Samples Per Task
    # Task,   n_mats_test
    # all,    8
    # Test Indexes: [ 2, 19, 41, 42, 50, 59, 60, 69 ]

    # Sample ID , Property Value        ,  Property Value (EST) ,  Feature 0 Value
    AgF         , -1.537576731789160e-01,  1.135790826401206e-01, -2.153206644755162e+00
    BeSe        ,  4.949404427752600e-01,  4.743878042320918e-01, -3.690720110296145e+00
    Ge2         ,  2.008525260607710e-01,  3.139008249590708e-01, -3.006837274572150e+00
    CGe         ,  8.114428802000480e-01,  3.691656586797722e-01, -3.242337197194732e+00
    FK          , -1.464060984981190e-01, -1.306050254917670e-01, -1.112665405432063e+00
    MgTe        , -4.591286648065000e-03, -6.222687007396171e-02, -1.404045098915806e+00
    BrNa        , -1.264287278827400e-01, -1.888626375989341e-01, -8.644123624074346e-01
    CSi         ,  6.690237272359810e-01,  3.948280949265375e-01, -3.351692484156472e+00

</details>

## Analyzing the Results with Python
*Note to do this part of the tutorial the python binding must also be built*

Once all of the calculations are completed the python interface provides some useful post-processing tools to easily analyze the results.
The `jackknife_cv_conv_est` tools provides a way to reasonably check the convergence of the cross-validation results with respect to the number number of calculations performed.
This tool uses [jackknife resampling](https://en.wikipedia.org/wiki/Jackknife_resampling) to calculate the mean and variance of the validation RMSEs across all cross-validation runs.
This technique essentially calculates the mean and standard error of all validation RMSEs for the system.
This data can then be used to estimate the overall validation RMSE for a given problem/set of hyper-parameters and the standard error associated with the random sampling of the test indexes.
It is important to mention that the error bars are based on the standard error of the mean of the validation RMSE, which assumes the sampling error follows a normal distribution.
Because the data set may not represent a uniform sampling of materials space, the standard error of the mean may only be a rough estimate of the true sampling error.
To visualize these results we will also use `plot_validation_rmse` at the end, and make the results easier to interpret.
```python
>>> from sissopp.postprocess.check_cv_convergence import jackknife_cv_conv_est
>>> from sissopp.postprocess.plot.cv_error_plot import plot_validation_rmse
>>> import numpy as np
>>> mean_val_rmse, var_val_rmse = jackknife_cv_conv_est("cv*")
>>> print(mean_val_rmse)
[0.16107651 0.10775548 0.10057917 0.10803806]
>>> print(np.sqrt(var_val_rmse))
[0.0159642  0.01732917 0.0118953  0.01717732]
>>> plot_validation_rmse("cv*", "cv_10._error.png").show()
```
Here is an example of the `plot_validation_rmse` output:
<details>
    <summary> Cross-Validation results </summary>

![image](./command_line/cv/cv_10_error.png)

</details>

These initial results suggest that we need to run more cross-validation samples in order to get converged results.
Using these results, we can only clearly state that the there is a significant decrease in the validation error when going from a one-dimensional model to a two-dimensional one.
However, because of the large error bars, it is impossible to determine which of the two, three, or four dimensional model is best.
To solve this lets increase the total number of samples to 100, and redo the analysis

```bash
for ii in `seq -f "%03g" 10 99`; do
    mkdir cv_$ii;
    cp sisso.json data.csv cv_$ii;
    cd cv_$ii;
    mpiexec -n 2 sisso++;
    cd ../;
done
```
```python
>>> from sissopp.postprocess.check_cv_convergence import jackknife_cv_conv_est
>>> from sissopp.postprocess.plot.cv_error_plot import plot_validation_rmse
>>> import numpy as np
>>> mean_val_rmse, var_val_rmse = jackknife_cv_conv_est("cv*")
>>> print(mean_val_rmse)
[0.15597268 0.12273297 0.10921321 0.10870643]
>>> print(np.sqrt(var_val_rmse))
[0.0051855  0.00571521 0.00398963 0.00473639]
>>> plot_validation_rmse("cv*", "cv_100._error.png").show()
```
With the additional calculations we now have relatively well converged results.
The key used in determining this is the relative size of the error bars when compared against the mean value.
For this example the estimate of the validation RMSEs for all dimensions up to the third dimension is outside the error bars of the other error bars, meaning that we can confidently say that the three-dimensional model is better than both the one and two-dimensional models.
Because the validation error for the three and four dimensional models are within each others error bars and the standard error increases when going to the fourth dimension, we can then conclude that the three-dimensional model has the ideal complexity.

<details>
    <summary> Converged cross-validation results </summary>

![image](./command_line/cv/cv_100_error.png)

</details>


## Visualizing the Cross Validation Error
The previous section illustrated how to plot the validation RMSE for each dimension of the model, but the RMSE does not give a complete picture of the model performance.
`SISSO++` also provides some utilities to plot the distribution of the average error for each sample.
To see the distributions for this system we run
```python
>>> from sissopp.postprocess.plot.cv_error_plot import plot_errors_dists
>>> plot_errors_dists("cv*", "error_cv_dist.png").show()
```
<details>
<summary> Distribution of Errors </summary>

![image](./command_line/cv/error_cv_dist.png)

</details>

These plots show the histogram of the error for each dimension with the total area normalized to one.
One thing that stands out in the plot is the large error seen in a single point for both the one and two dimensional models.
By looking at the validation errors, we find that the point with the largest error is diamond for all model dimensions, which is by far the most stable zinc-blende structure in the data set.
As a note for this setup there is a 0.22\% chance that one of the samples is never in the validation set so if `max_error_ind != 21` check if that sample is in one of the validation sets.

```python
>>> import numpy as np
>>> import pandas as pd
>>> from sissopp.postprocess.get_model_errors import get_model_errors
>>> 
>>> df = pd.read_csv("data.csv", index_col=0)
>>> 
>>> te, ve = get_model_errors("cv*", True)
>>> max_err_ind = np.nanmean(ve, axis=1).argmax(axis=0)
>>> print(df.index[max_err_ind])
Index(['C2', 'C2', 'C2', 'C2'], dtype='object', name='# Material')
```

## Optimizing the hyper-parameters of SISSO
As discussed in the previous example `desc_dim` is one of the four hyperparameters used in `SISSO++` with the others being: `n_sis_select`, `max_rung`, and `n_residual`.
Of these `n_sis_select` and `n_residual` need to be optimized together while `desc_dim` and `max_rung` can be optimized independently.
Due to the factorial increase in both computational time and required memory associated with `max_rung` only `desc_dim`, `n_sis_select`, and `n_residual` will be optimized in this exercise, but for production purposes this will also have to be studied.
Additionally the exercise will only use use relatively small SIS subspace sizes and only go up to a 3D model in order to reduce the computational time for the exercise.
The first step of this process will be setting up nine directories for each combination `n_residual` (1, 5, and 10) and `n_sis_select` (10, 50, 100) and modify the base `sisso.json` to match these new parameters (Note: the dimension of the final model will be determined in the same way as the previous example).
Here is the new base `sisso.json file`:

```json
{
    "data_file": "data.csv",
    "property_key": "E_RS - E_ZB",
    "desc_dim": 3,
    "n_sis_select": SSSS,
    "max_rung": 2,
    "calc_type": "regression",
    "min_abs_feat_val": 1e-5,
    "max_abs_feat_val": 1e8,
    "n_residual": RRRR,
    "n_models_store": 1,
    "leave_out_frac": 0.05,
    "leave_out_inds": [],
    "opset": ["add", "sub", "abs_diff", "mult", "div", "inv", "abs", "exp", "log", "sin", "cos", "sq", "cb", "six_pow", "sqrt", "cbrt", "neg_exp"]

}
```

From this base we will set up the nine directories with the following set of commands (Note: we will simply copy over the previously calculated cross-validation data from the previous exercise)
```bash
mkdir ns_10_nr_1;
cp data.csv ns_10_nr_1/;
sed "s/SSSS/10/g" sisso.json > ns_10_nr_1/sisso.json;
sed -i "s/RRRR/1/g" ns_10_nr_1/sisso.json;

mkdir ns_10_nr_5;
cp data.csv ns_10_nr_5/;
sed "s/SSSS/10/g" sisso.json > ns_10_nr_5/sisso.json;
sed -i "s/RRRR/5/g" ns_10_nr_5/sisso.json;

mkdir ns_10_nr_10;
cp data.csv ns_10_nr_10/;
sed "s/SSSS/10/g" sisso.json > ns_10_nr_10/sisso.json;
sed -i "s/RRRR/10/g" ns_10_nr_10/sisso.json;
mv cv* ns_10_nr_10

mkdir ns_50_nr_1;
cp data.csv ns_50_nr_1;
sed "s/SSSS/50/g" sisso.json > ns_50_nr_1/sisso.json;
sed -i "s/RRRR/1/g" ns_50_nr_1/sisso.json;

mkdir ns_50_nr_5;
cp data.csv ns_50_nr_5;
sed "s/SSSS/50/g" sisso.json > ns_50_nr_5/sisso.json;
sed -i "s/RRRR/5/g" ns_50_nr_5/sisso.json;

mkdir ns_50_nr_10;
cp data.csv ns_50_nr_10
sed "s/SSSS/50/g" sisso.json > ns_50_nr_10/sisso.json;
sed -i "s/RRRR/10/g" ns_50_nr_10/sisso.json;

mkdir ns_100_nr_1;
cp data.csv ns_100_nr_1;
sed "s/SSSS/100/g" sisso.json > ns_100_nr_1/sisso.json;
sed -i "s/RRRR/1/g" ns_100_nr_1/sisso.json;

mkdir ns_100_nr_5;
cp data.csv ns_100_nr_5;
sed "s/SSSS/100/g" sisso.json > ns_100_nr_5/sisso.json;
sed -i "s/RRRR/5/g" ns_100_nr_5/sisso.json;

mkdir ns_100_nr_10;
cp data.csv ns_100_nr_10;
sed "s/SSSS/100/g" sisso.json > ns_100_nr_10/sisso.json;
sed -i "s/RRRR/10/g" ns_100_nr_10/sisso.json;
```
From here we can run the same cross-validation analysis from the previous section in each of the results, and then compare the results in python
```python
>>> from sissopp.postprocess.check_cv_convergence import jackknife_cv_conv_est
>>> import numpy as np
>>> 
>>> mean_10_1, var_10_1 = jackknife_cv_conv_est("ns_10_nr_1/cv*")
>>> mean_10_5, var_10_5 = jackknife_cv_conv_est("ns_10_nr_5/cv*")
>>> mean_10_10, var_10_10 = jackknife_cv_conv_est("ns_10_nr_10/cv*")
>>> mean_50_1, var_50_1 = jackknife_cv_conv_est("ns_50_nr_1/cv*")
>>> mean_50_5, var_50_5 = jackknife_cv_conv_est("ns_50_nr_5/cv*")
>>> mean_50_10, var_50_10 = jackknife_cv_conv_est("ns_50_nr_10/cv*")
>>> mean_100_1, var_100_1 = jackknife_cv_conv_est("ns_100_nr_1/cv*")
>>> mean_100_5, var_100_5 = jackknife_cv_conv_est("ns_100_nr_5/cv*")
>>> mean_100_10, var_100_10 = jackknife_cv_conv_est("ns_100_nr_10/cv*")
>>> 
>>> print_str = f"ns:  10; nr:  1; {mean_10_1[:3]} {np.sqrt(var_10_1[:3])}\n"
>>> print_str += f"ns:  10; nr:  5; {mean_10_5[:3]} {np.sqrt(var_10_5[:3])}\n"
>>> print_str += f"ns:  10; nr: 10; {mean_10_10[:3]} {np.sqrt(var_10_10[:3])}\n"
>>> print_str += f"ns:  50; nr:  1; {mean_50_1[:3]} {np.sqrt(var_50_1[:3])}\n"
>>> print_str += f"ns:  50; nr:  5; {mean_50_5[:3]} {np.sqrt(var_50_5[:3])}\n"
>>> print_str += f"ns:  50; nr: 10; {mean_50_10[:3]} {np.sqrt(var_50_10[:3])}\n"
>>> print_str += f"ns: 100; nr:  1; {mean_100_1[:3]} {np.sqrt(var_100_1[:3])}\n"
>>> print_str += f"ns: 100; nr:  5; {mean_100_5[:3]} {np.sqrt(var_100_5[:3])}\n"
>>> print_str += f"ns: 100; nr: 10; {mean_100_10[:3]} {np.sqrt(var_100_10[:3])}"
>>> print(print_str)
ns:  10; nr:  1; [0.15680869 0.17389737 0.16029643] [0.00646652 0.04735888 0.04424095]
ns:  10; nr:  5; [0.15625644 0.12419926 0.15115378] [0.00663913 0.00631875 0.04471696]
ns:  10; nr: 10; [0.15597268 0.12273297 0.10921321] [0.0051855  0.00571521 0.00398963]
ns:  50; nr:  1; [0.15192553 0.12373729 0.13507366] [0.00513279 0.00523007 0.01695709]
ns:  50; nr:  5; [0.15262692 0.12672067 0.11011062] [0.00491465 0.00522488 0.00407753]
ns:  50; nr: 10; [0.15119692 0.13040251 0.10993919] [0.00487215 0.00506964 0.00464264]
ns: 100; nr:  1; [0.15835654 0.13728706 0.12849654] [0.00557331 0.00606579 0.01114889]
ns: 100; nr:  5; [0.15502757 0.14002783 0.12102758] [0.00489507 0.00546934 0.00612467]
ns: 100; nr: 10; [0.14996602 0.13248817 0.1070521 ] [0.00495617 0.0051647  0.00432492]
```
These results indicate that for the small SIS subspace sizes used here the validation error is stable relative to both the number of residuals and SIS subspace size, given the simliarity
However it is important to note that this will not always be the case, particularly for for larger values of `n_sis_select` .
For finding the best model over all an `n_sis_select` of 100 and `n_residual` of 10 will be used.
This choice was made because it has the lowest validation RMSE of 0.107, but all calculations that use 10 residuals will have equivalent performance (at least for the small SIS subspace size).

## Final Training Data
To get the final models we will perform the same calculation we started off the tutorial with, but with the following `sisso.json` file based on the results from the previous steps:
```
{
    "data_file": "data.csv",
    "property_key": "E_RS - E_ZB",
    "desc_dim": 3,
    "n_sis_select": 100,
    "max_rung": 2,
    "calc_type": "regression",
    "min_abs_feat_val": 1e-5,
    "max_abs_feat_val": 1e8,
    "n_residual": 10,
    "n_models_store": 1,
    "leave_out_frac": 0.0,
    "leave_out_inds": [],
    "opset": ["add", "sub", "abs_diff", "mult", "div", "inv", "abs", "exp", "log", "sin", "cos", "sq", "cb", "six_pow", "sqrt", "cbrt", "neg_exp"]
}
```
From here we can use `models/train_dim_3_model_0.dat` for all of the analysis.
In order to generate a machine learning plot for this model in matplotlib, run the following in python
```python
>>> from sissopp.postprocess.plot.parity_plot import plot_model_parity_plot
>>> plot_model_parity_plot("models/train_dim_3_model_0.dat", filename="3d_model.pdf").show()
```
The result of which is shown below:
<details>
<summary> Final 3D model </summary>

![image](./command_line/cv/3d_model.png)
</details>

Additionally you can generate a output the model as a Matlab function or a LaTeX string using the following commands.
```python
>>> from sissopp.postprocess.load_models import load_model
>>> model = load_model("models/train_dim_3_model_0.dat")
>>> print(model.latex_str)
>>> 
>>> model.write_matlab_fxn("matlab_fxn/model.m")
```

A copy of the generated matlab function is below.
<details>
<summary> Matlab function of the Final 3D model </summary>


    function P = model(X)
    % Returns the value of E_{RS} - E_{ZB} = c0 + a0 * ((r_d_B / r_d_A) * (r_p_B * E_HOMO_A)) + a1 * ((IP_A^3) * (|r_sigma - r_s_B|)) + a2 * ((IP_A / r_p_A) / (r_p_B + r_p_A))
    %
    % X = [
    %     r_d_B,
    %     r_d_A,
    %     r_p_B,
    %     E_HOMO_A,
    %     IP_A,
    %     r_sigma,
    %     r_s_B,
    %     r_p_A,
    % ]

    if(size(X, 2) ~= 8)
        error("ERROR: X must have a size of 8 in the second dimension.")
    end
    r_d_B    = reshape(X(:, 1), 1, []);
    r_d_A    = reshape(X(:, 2), 1, []);
    r_p_B    = reshape(X(:, 3), 1, []);
    E_HOMO_A = reshape(X(:, 4), 1, []);
    IP_A     = reshape(X(:, 5), 1, []);
    r_sigma  = reshape(X(:, 6), 1, []);
    r_s_B    = reshape(X(:, 7), 1, []);
    r_p_A    = reshape(X(:, 8), 1, []);

    f0 = ((r_d_B ./ r_d_A) .* (r_p_B .* E_HOMO_A));
    f1 = ((IP_A).^3 .* abs(r_sigma - r_s_B));
    f2 = ((IP_A ./ r_p_A) ./ (r_p_B + r_p_A));

    c0 = -1.3509197357e-01;
    a0 = 2.8311062079e-02;
    a1 = 3.7282871777e-04;
    a2 = -2.3703222974e-01;

    P = reshape(c0 + a0 * f0 + a1 * f1 + a2 * f2, [], 1);
    end

</details>
# Performing Classification with SISSO++

Finally, besides regression problems, `SISSO++` can be used to solve classification problems.
While we have already that this problem could solved with regression, `SISSO++` can also solve problems that can't be treated as a regression problem.
As an example of this we will adapt the previous example by replacing the property with the identifier of if the material favors the rock-salt or zinc-blende structure, and change the calculation type to be `classification`.
It is important to note that while this problem only has two classes, multi-class classification is also possible.

## The Data File
Here is the updated data file, with the property `E_RS - E_ZB (eV)` replaced with a `Class` column where any negative `E_RS - E_ZB (eV)` is replaced with 0 and any positive value replaced with 1. While this example has only one task and two classes, the method works for an arbitrary number of classes and tasks.

<details>
    <summary>Here is the full data_class.csv file for the calculation</summary>

    # Material,Class,Z_A (nuc_charge) ,Z_B (nuc_charge) ,period_A,period_B,IP_A (eV_IP) ,IP_B (eV_IP) ,EA_A (eV_IP),EA_B (eV_IP) ,E_HOMO_A (eV) ,E_HOMO_B (eV) ,E_LUMO_A (eV),E_LUMO_B (eV) ,r_s_A ,r_s_B ,r_p_A ,r_p_B ,r_d_A ,r_d_B,r_sigma ,r_pi
    AgBr,0,47,35,5,4,-8.0580997467,-12.649600029,-1.66659998894,-3.73930001259,-4.71000003815,-8.00100040436,-0.479000002146,0.708000004292,1.32000005245,0.75,1.87999999523,0.879999995232,2.97000002861,1.87000000477,1.570000052448,0.689999938012
    AgCl,0,47,17,5,3,-8.0580997467,-13.9018001556,-1.66659998894,-3.97079992294,-4.71000003815,-8.69999980927,-0.479000002146,0.574000000954,1.32000005245,0.680000007153,1.87999999523,0.759999990463,2.97000002861,1.66999995708,1.760000050064,0.63999992609
    AgF,0,47,9,5,2,-8.0580997467,-19.4043006897,-1.66659998894,-4.27349996567,-4.71000003815,-11.2939996719,-0.479000002146,1.25100004673,1.32000005245,0.409999996424,1.87999999523,0.370000004768,2.97000002861,1.42999994755,2.420000046488,0.599999934436
    AgI,1,47,53,5,5,-8.0580997467,-11.2571001053,-1.66659998894,-3.5134999752,-4.71000003815,-7.23600006104,-0.479000002146,0.212999999523,1.32000005245,0.899999976158,1.87999999523,1.07000005245,2.97000002861,1.72000002861,1.230000019072,0.730000019072
    AlAs,1,13,33,3,4,-5.78049993515,-9.26189994812,-0.3125,-1.83920001984,-2.78399991989,-5.34100008011,0.694999992847,0.0640000030398,1.09000003338,0.850000023842,1.38999998569,1.03999996185,1.94000005722,2.01999998093,0.590000033378,0.489999890318
    AlN,1,13,7,3,2,-5.78049993515,-13.5852003098,-0.3125,-1.86749994755,-2.78399991989,-7.2389998436,0.694999992847,3.0569999218,1.09000003338,0.540000021458,1.38999998569,0.509999990463,1.94000005722,1.53999996185,1.430000007149,0.329999983305
    AlP,1,13,15,3,3,-5.78049993515,-9.75059986115,-0.3125,-1.91999995708,-2.78399991989,-5.59600019455,0.694999992847,0.182999998331,1.09000003338,0.829999983311,1.38999998569,0.97000002861,1.94000005722,1.76999998093,0.680000007149,0.439999997609
    AlSb,1,13,51,3,5,-5.78049993515,-8.46829986572,-0.3125,-1.84669995308,-2.78399991989,-4.99100017548,0.694999992847,0.104999996722,1.09000003338,1,1.38999998569,1.23000001907,1.94000005722,2.05999994278,0.25,0.52999997138
    AsGa,1,31,33,4,4,-5.81820011139,-9.26189994812,-0.108099997044,-1.83920001984,-2.73200011253,-5.34100008011,0.129999995232,0.0640000030398,0.990000009537,0.850000023842,1.33000004292,1.03999996185,2.16000008583,2.01999998093,0.430000066765,0.529999971391
    AsB,1,5,33,2,4,-8.18999958038,-9.26189994812,-0.107400000095,-1.83920001984,-3.71499991417,-5.34100008011,2.24799990654,0.0640000030398,0.810000002384,0.850000023842,0.829999983311,1.03999996185,1.95000004768,2.01999998093,0.249999999997,0.209999918935
    BN,1,5,7,2,2,-8.18999958038,-13.5852003098,-0.107400000095,-1.86749994755,-3.71499991417,-7.2389998436,2.24799990654,3.0569999218,0.810000002384,0.540000021458,0.829999983311,0.509999990463,1.95000004768,1.53999996185,0.589999973774,0.050000011922
    BP,1,5,15,2,3,-8.18999958038,-9.75059986115,-0.107400000095,-1.91999995708,-3.71499991417,-5.59600019455,2.24799990654,0.182999998331,0.810000002384,0.829999983311,0.829999983311,0.97000002861,1.95000004768,1.76999998093,0.160000026226,0.160000026226
    BSb,1,5,51,2,5,-8.18999958038,-8.46829986572,-0.107400000095,-1.84669995308,-3.71499991417,-4.99100017548,2.24799990654,0.104999996722,0.810000002384,1,0.829999983311,1.23000001907,1.95000004768,2.05999994278,0.590000033375,0.249999999997
    BaO,0,56,8,6,2,-5.51569986343,-16.4332008362,0.277999997139,-3.00589990616,-3.34599995613,-9.19699954987,-2.1289999485,2.54099988937,2.15000009537,0.460000008345,2.63000011444,0.430000007153,1.35000002384,2.22000002861,3.890000194312,0.510000020262
    BaS,0,56,16,6,3,-5.51569986343,-11.7951002121,0.277999997139,-2.84489989281,-3.34599995613,-7.10599994659,-2.1289999485,0.64200001955,2.15000009537,0.740000009537,2.63000011444,0.850000023842,1.35000002384,2.36999988556,3.190000176431,0.590000033375
    BaSe,0,56,34,6,4,-5.51569986343,-10.9460000992,0.277999997139,-2.75099992752,-3.34599995613,-6.65399980545,-2.1289999485,1.31599998474,2.15000009537,0.800000011921,2.63000011444,0.949999988079,1.35000002384,2.18000006676,3.03000020981,0.629999995228
    BaTe,0,56,52,6,5,-5.51569986343,-9.86670017242,0.277999997139,-2.66599988937,-3.34599995613,-6.10900020599,-2.1289999485,0.0989999994636,2.15000009537,0.939999997616,2.63000011444,1.13999998569,1.35000002384,1.83000004292,2.700000226504,0.680000007144
    BeO,1,4,8,2,2,-9.459400177,-16.4332008362,0.630500018597,-3.00589990616,-5.59999990463,-9.19699954987,-2.09800004959,2.54099988937,1.08000004292,0.460000008345,1.21000003815,0.430000007153,2.88000011444,2.22000002861,1.400000065572,0.159999996422
    BeS,1,4,16,2,3,-9.459400177,-11.7951002121,0.630500018597,-2.84489989281,-5.59999990463,-7.10599994659,-2.09800004959,0.64200001955,1.08000004292,0.740000009537,1.21000003815,0.850000023842,2.88000011444,2.36999988556,0.700000047691,0.240000009535
    BeSe,1,4,34,2,4,-9.459400177,-10.9460000992,0.630500018597,-2.75099992752,-5.59999990463,-6.65399980545,-2.09800004959,1.31599998474,1.08000004292,0.800000011921,1.21000003815,0.949999988079,2.88000011444,2.18000006676,0.54000008107,0.279999971388
    BeTe,1,4,52,2,5,-9.459400177,-9.86670017242,0.630500018597,-2.66599988937,-5.59999990463,-6.10900020599,-2.09800004959,0.0989999994636,1.08000004292,0.939999997616,1.21000003815,1.13999998569,2.88000011444,1.83000004292,0.210000097764,0.329999983304
    C2,1,6,6,2,2,-10.8516998291,-10.8516998291,-0.87239998579,-0.87239998579,-5.41599988937,-5.41599988937,1.99199998379,1.99199998379,0.639999985695,0.639999985695,0.629999995232,0.629999995232,1.62999999523,1.62999999523,0,0.019999980926
    CaO,0,20,8,4,2,-6.4279999733,-16.4332008362,0.303900003433,-3.00589990616,-3.86400008202,-9.19699954987,-2.132999897,2.54099988937,1.75999999046,0.460000008345,2.31999993324,0.430000007153,0.680000007153,2.22000002861,3.189999908202,0.589999943972
    CaS,0,20,16,4,3,-6.4279999733,-11.7951002121,0.303900003433,-2.84489989281,-3.86400008202,-7.10599994659,-2.132999897,0.64200001955,1.75999999046,0.740000009537,2.31999993324,0.850000023842,0.680000007153,2.36999988556,2.489999890321,0.669999957085
    CaSe,0,20,34,4,4,-6.4279999733,-10.9460000992,0.303900003433,-2.75099992752,-3.86400008202,-6.65399980545,-2.132999897,1.31599998474,1.75999999046,0.800000011921,2.31999993324,0.949999988079,0.680000007153,2.18000006676,2.3299999237,0.709999918938
    CaTe,0,20,52,4,5,-6.4279999733,-9.86670017242,0.303900003433,-2.66599988937,-3.86400008202,-6.10900020599,-2.132999897,0.0989999994636,1.75999999046,0.939999997616,2.31999993324,1.13999998569,0.680000007153,1.83000004292,1.999999940394,0.759999930854
    CdO,0,48,8,5,2,-9.5813999176,-16.4332008362,0.838699996471,-3.00589990616,-5.95200014114,-9.19699954987,-1.30900001526,2.54099988937,1.23000001907,0.460000008345,1.74000000954,0.430000007153,2.59999990463,2.22000002861,2.080000013112,0.539999991662
    CdS,1,48,16,5,3,-9.5813999176,-11.7951002121,0.838699996471,-2.84489989281,-5.95200014114,-7.10599994659,-1.30900001526,0.64200001955,1.23000001907,0.740000009537,1.74000000954,0.850000023842,2.59999990463,2.36999988556,1.379999995231,0.620000004775
    CdSe,1,48,34,5,4,-9.5813999176,-10.9460000992,0.838699996471,-2.75099992752,-5.95200014114,-6.65399980545,-1.30900001526,1.31599998474,1.23000001907,0.800000011921,1.74000000954,0.949999988079,2.59999990463,2.18000006676,1.22000002861,0.659999966628
    CdTe,1,48,52,5,5,-9.5813999176,-9.86670017242,0.838699996471,-2.66599988937,-5.95200014114,-6.10900020599,-1.30900001526,0.0989999994636,1.23000001907,0.939999997616,1.74000000954,1.13999998569,2.59999990463,1.83000004292,0.890000045304,0.709999978544
    BrCs,0,55,35,6,4,-4.00619983673,-12.649600029,-0.569599986076,-3.73930001259,-2.22000002861,-8.00100040436,-0.547999978065,0.708000004292,2.46000003815,0.75,3.16000008583,0.879999995232,1.97000002861,1.87000000477,3.990000128748,0.830000042912
    ClCs,0,55,17,6,3,-4.00619983673,-13.9018001556,-0.569599986076,-3.97079992294,-2.22000002861,-8.69999980927,-0.547999978065,0.574000000954,2.46000003815,0.680000007153,3.16000008583,0.759999990463,1.97000002861,1.66999995708,4.180000126364,0.78000003099
    CsF,0,55,9,6,2,-4.00619983673,-19.4043006897,-0.569599986076,-4.27349996567,-2.22000002861,-11.2939996719,-0.547999978065,1.25100004673,2.46000003815,0.409999996424,3.16000008583,0.370000004768,1.97000002861,1.42999994755,4.840000122788,0.740000039336
    CsI,0,55,53,6,5,-4.00619983673,-11.2571001053,-0.569599986076,-3.5134999752,-2.22000002861,-7.23600006104,-0.547999978065,0.212999999523,2.46000003815,0.899999976158,3.16000008583,1.07000005245,1.97000002861,1.72000002861,3.650000095372,0.870000123972
    BrCu,1,29,35,4,4,-8.38879966736,-12.649600029,-1.6384999752,-3.73930001259,-4.85599994659,-8.00100040436,-0.64099997282,0.708000004292,1.20000004768,0.75,1.67999994755,0.879999995232,2.57999992371,1.87000000477,1.249999999998,0.609999895102
    ClCu,1,29,17,4,3,-8.38879966736,-13.9018001556,-1.6384999752,-3.97079992294,-4.85599994659,-8.69999980927,-0.64099997282,0.574000000954,1.20000004768,0.680000007153,1.67999994755,0.759999990463,2.57999992371,1.66999995708,1.439999997614,0.55999988318
    CuF,0,29,9,4,2,-8.38879966736,-19.4043006897,-1.6384999752,-4.27349996567,-4.85599994659,-11.2939996719,-0.64099997282,1.25100004673,1.20000004768,0.409999996424,1.67999994755,0.370000004768,2.57999992371,1.42999994755,2.099999994038,0.519999891526
    CuI,1,29,53,4,5,-8.38879966736,-11.2571001053,-1.6384999752,-3.5134999752,-4.85599994659,-7.23600006104,-0.64099997282,0.212999999523,1.20000004768,0.899999976158,1.67999994755,1.07000005245,2.57999992371,1.72000002861,0.909999966622,0.649999976162
    GaN,1,31,7,4,2,-5.81820011139,-13.5852003098,-0.108099997044,-1.86749994755,-2.73200011253,-7.2389998436,0.129999995232,3.0569999218,0.990000009537,0.540000021458,1.33000004292,0.509999990463,2.16000008583,1.53999996185,1.270000040536,0.370000064378
    GaP,1,31,15,4,3,-5.81820011139,-9.75059986115,-0.108099997044,-1.91999995708,-2.73200011253,-5.59600019455,0.129999995232,0.182999998331,0.990000009537,0.829999983311,1.33000004292,0.97000002861,2.16000008583,1.76999998093,0.520000040536,0.480000078682
    GaSb,1,31,51,4,5,-5.81820011139,-8.46829986572,-0.108099997044,-1.84669995308,-2.73200011253,-4.99100017548,0.129999995232,0.104999996722,0.990000009537,1,1.33000004292,1.23000001907,2.16000008583,2.05999994278,0.090000033387,0.570000052453
    Ge2,1,32,32,4,4,-7.56699991226,-7.56699991226,-0.949000000954,-0.949000000954,-4.04600000381,-4.04600000381,2.17499995232,2.17499995232,0.920000016689,0.920000016689,1.15999996662,1.15999996662,2.36999988556,2.36999988556,0,0.479999899862
    CGe,1,32,6,4,2,-7.56699991226,-10.8516998291,-0.949000000954,-0.87239998579,-4.04600000381,-5.41599988937,2.17499995232,1.99199998379,0.920000016689,0.639999985695,1.15999996662,0.629999995232,2.36999988556,1.62999999523,0.810000002382,0.249999940394
    GeSi,1,32,14,4,3,-7.56699991226,-7.75769996643,-0.949000000954,-0.992999970913,-4.04600000381,-4.16300010681,2.17499995232,0.439999997616,0.920000016689,0.939999997616,1.15999996662,1.12999999523,2.36999988556,1.88999998569,0.009999990463,0.429999947545
    AsIn,1,49,33,5,4,-5.53739976883,-9.26189994812,-0.256300002337,-1.83920001984,-2.6970000267,-5.34100008011,0.368000000715,0.0640000030398,1.12999999523,0.850000023842,1.5,1.03999996185,3.1099998951,2.01999998093,0.740000009538,0.559999942778
    InN,1,49,7,5,2,-5.53739976883,-13.5852003098,-0.256300002337,-1.86749994755,-2.6970000267,-7.2389998436,0.368000000715,3.0569999218,1.12999999523,0.540000021458,1.5,0.509999990463,3.1099998951,1.53999996185,1.579999983309,0.400000035765
    InP,1,49,15,5,3,-5.53739976883,-9.75059986115,-0.256300002337,-1.91999995708,-2.6970000267,-5.59600019455,0.368000000715,0.182999998331,1.12999999523,0.829999983311,1.5,0.97000002861,3.1099998951,1.76999998093,0.829999983309,0.510000050069
    InSb,1,49,51,5,5,-5.53739976883,-8.46829986572,-0.256300002337,-1.84669995308,-2.6970000267,-4.99100017548,0.368000000715,0.104999996722,1.12999999523,1,1.5,1.23000001907,3.1099998951,2.05999994278,0.39999997616,0.60000002384
    BrK,0,19,35,4,4,-4.43319988251,-12.649600029,-0.621299982071,-3.73930001259,-2.42600011826,-8.00100040436,-0.697000026703,0.708000004292,2.13000011444,0.75,2.44000005722,0.879999995232,1.78999996185,1.87000000477,2.940000176428,0.439999938012
    ClK,0,19,17,4,3,-4.43319988251,-13.9018001556,-0.621299982071,-3.97079992294,-2.42600011826,-8.69999980927,-0.697000026703,0.574000000954,2.13000011444,0.680000007153,2.44000005722,0.759999990463,1.78999996185,1.66999995708,3.130000174044,0.38999992609
    FK,0,19,9,4,2,-4.43319988251,-19.4043006897,-0.621299982071,-4.27349996567,-2.42600011826,-11.2939996719,-0.697000026703,1.25100004673,2.13000011444,0.409999996424,2.44000005722,0.370000004768,1.78999996185,1.42999994755,3.790000170468,0.349999934436
    IK,0,19,53,4,5,-4.43319988251,-11.2571001053,-0.621299982071,-3.5134999752,-2.42600011826,-7.23600006104,-0.697000026703,0.212999999523,2.13000011444,0.899999976158,2.44000005722,1.07000005245,1.78999996185,1.72000002861,2.600000143052,0.480000019072
    BrLi,0,3,35,2,4,-5.32910013199,-12.649600029,-0.698099970818,-3.73930001259,-2.87400007248,-8.00100040436,-0.977999985218,0.708000004292,1.64999997616,0.75,2,0.879999995232,6.92999982834,1.87000000477,2.019999980928,0.480000019072
    ClLi,0,3,17,2,3,-5.32910013199,-13.9018001556,-0.698099970818,-3.97079992294,-2.87400007248,-8.69999980927,-0.977999985218,0.574000000954,1.64999997616,0.680000007153,2,0.759999990463,6.92999982834,1.66999995708,2.209999978544,0.43000000715
    FLi,0,3,9,2,2,-5.32910013199,-19.4043006897,-0.698099970818,-4.27349996567,-2.87400007248,-11.2939996719,-0.977999985218,1.25100004673,1.64999997616,0.409999996424,2,0.370000004768,6.92999982834,1.42999994755,2.869999974968,0.390000015496
    ILi,0,3,53,2,5,-5.32910013199,-11.2571001053,-0.698099970818,-3.5134999752,-2.87400007248,-7.23600006104,-0.977999985218,0.212999999523,1.64999997616,0.899999976158,2,1.07000005245,6.92999982834,1.72000002861,1.679999947552,0.520000100132
    MgO,0,12,8,3,2,-8.03709983826,-16.4332008362,0.692499995232,-3.00589990616,-4.78200006485,-9.19699954987,-1.35800004005,2.54099988937,1.33000004292,0.460000008345,1.89999997616,0.430000007153,3.17000007629,2.22000002861,2.340000003582,0.599999934432
    MgS,0,12,16,3,3,-8.03709983826,-11.7951002121,0.692499995232,-2.84489989281,-4.78200006485,-7.10599994659,-1.35800004005,0.64200001955,1.33000004292,0.740000009537,1.89999997616,0.850000023842,3.17000007629,2.36999988556,1.639999985701,0.679999947545
    MgSe,0,12,34,3,4,-8.03709983826,-10.9460000992,0.692499995232,-2.75099992752,-4.78200006485,-6.65399980545,-1.35800004005,1.31599998474,1.33000004292,0.800000011921,1.89999997616,0.949999988079,3.17000007629,2.18000006676,1.48000001908,0.719999909398
    MgTe,0,12,52,3,5,-8.03709983826,-9.86670017242,0.692499995232,-2.66599988937,-4.78200006485,-6.10900020599,-1.35800004005,0.0989999994636,1.33000004292,0.939999997616,1.89999997616,1.13999998569,3.17000007629,1.83000004292,1.150000035774,0.769999921314
    BrNa,0,11,35,3,4,-5.22310018539,-12.649600029,-0.715699970722,-3.73930001259,-2.81900000572,-8.00100040436,-0.717999994755,0.708000004292,1.71000003815,0.75,2.59999990463,0.879999995232,6.57000017166,1.87000000477,2.679999947548,1.019999861712
    ClNa,0,11,17,3,3,-5.22310018539,-13.9018001556,-0.715699970722,-3.97079992294,-2.81900000572,-8.69999980927,-0.717999994755,0.574000000954,1.71000003815,0.680000007153,2.59999990463,0.759999990463,6.57000017166,1.66999995708,2.869999945164,0.96999984979
    FNa,0,11,9,3,2,-5.22310018539,-19.4043006897,-0.715699970722,-4.27349996567,-2.81900000572,-11.2939996719,-0.717999994755,1.25100004673,1.71000003815,0.409999996424,2.59999990463,0.370000004768,6.57000017166,1.42999994755,3.529999941588,0.929999858136
    INa,0,11,53,3,5,-5.22310018539,-11.2571001053,-0.715699970722,-3.5134999752,-2.81900000572,-7.23600006104,-0.717999994755,0.212999999523,1.71000003815,0.899999976158,2.59999990463,1.07000005245,6.57000017166,1.72000002861,2.339999914172,1.059999942772
    BrRb,0,37,35,5,4,-4.28889989853,-12.649600029,-0.590399980545,-3.73930001259,-2.3599998951,-8.00100040436,-0.704999983311,0.708000004292,2.24000000954,0.75,3.20000004768,0.879999995232,1.96000003815,1.87000000477,3.810000061988,1.090000033372
    ClRb,0,37,17,5,3,-4.28889989853,-13.9018001556,-0.590399980545,-3.97079992294,-2.3599998951,-8.69999980927,-0.704999983311,0.574000000954,2.24000000954,0.680000007153,3.20000004768,0.759999990463,1.96000003815,1.66999995708,4.000000059604,1.04000002145
    FRb,0,37,9,5,2,-4.28889989853,-19.4043006897,-0.590399980545,-4.27349996567,-2.3599998951,-11.2939996719,-0.704999983311,1.25100004673,2.24000000954,0.409999996424,3.20000004768,0.370000004768,1.96000003815,1.42999994755,4.660000056028,1.000000029796
    IRb,0,37,53,5,5,-4.28889989853,-11.2571001053,-0.590399980545,-3.5134999752,-2.3599998951,-7.23600006104,-0.704999983311,0.212999999523,2.24000000954,0.899999976158,3.20000004768,1.07000005245,1.96000003815,1.72000002861,3.470000028612,1.130000114432
    Si2,1,14,14,3,3,-7.75769996643,-7.75769996643,-0.992999970913,-0.992999970913,-4.16300010681,-4.16300010681,0.439999997616,0.439999997616,0.939999997616,0.939999997616,1.12999999523,1.12999999523,1.88999998569,1.88999998569,0,0.379999995228
    CSi,1,14,6,3,2,-7.75769996643,-10.8516998291,-0.992999970913,-0.87239998579,-4.16300010681,-5.41599988937,0.439999997616,1.99199998379,0.939999997616,0.639999985695,1.12999999523,0.629999995232,1.88999998569,1.62999999523,0.800000011919,0.199999988077
    Sn2,1,50,50,5,5,-7.04279994965,-7.04279994965,-1.03919994831,-1.03919994831,-3.86599993706,-3.86599993706,0.00800000037998,0.00800000037998,1.05999994278,1.05999994278,1.34000003338,1.34000003338,2.02999997139,2.02999997139,0,0.5600001812
    CSn,1,50,6,5,2,-7.04279994965,-10.8516998291,-1.03919994831,-0.87239998579,-3.86599993706,-5.41599988937,0.00800000037998,1.99199998379,1.05999994278,0.639999985695,1.34000003338,0.629999995232,2.02999997139,1.62999999523,1.129999995233,0.290000081063
    GeSn,1,50,32,5,4,-7.04279994965,-7.56699991226,-1.03919994831,-0.949000000954,-3.86599993706,-4.04600000381,0.00800000037998,2.17499995232,1.05999994278,0.920000016689,1.34000003338,1.15999996662,2.02999997139,2.36999988556,0.319999992851,0.520000040531
    SiSn,1,50,14,5,3,-7.04279994965,-7.75769996643,-1.03919994831,-0.992999970913,-3.86599993706,-4.16300010681,0.00800000037998,0.439999997616,1.05999994278,0.939999997616,1.34000003338,1.12999999523,2.02999997139,1.88999998569,0.329999983314,0.470000088214
    OSr,0,38,8,5,2,-6.03159999847,-16.4332008362,0.343100011349,-3.00589990616,-3.64100003242,-9.19699954987,-1.3789999485,2.54099988937,1.90999996662,0.460000008345,2.54999995232,0.430000007153,1.20000004768,2.22000002861,3.569999903442,0.669999986892
    SSr,0,38,16,5,3,-6.03159999847,-11.7951002121,0.343100011349,-2.84489989281,-3.64100003242,-7.10599994659,-1.3789999485,0.64200001955,1.90999996662,0.740000009537,2.54999995232,0.850000023842,1.20000004768,2.36999988556,2.869999885561,0.750000000005
    SeSr,0,38,34,5,4,-6.03159999847,-10.9460000992,0.343100011349,-2.75099992752,-3.64100003242,-6.65399980545,-1.3789999485,1.31599998474,1.90999996662,0.800000011921,2.54999995232,0.949999988079,1.20000004768,2.18000006676,2.70999991894,0.789999961858
    SrTe,0,38,52,5,5,-6.03159999847,-9.86670017242,0.343100011349,-2.66599988937,-3.64100003242,-6.10900020599,-1.3789999485,0.0989999994636,1.90999996662,0.939999997616,2.54999995232,1.13999998569,1.20000004768,1.83000004292,2.379999935634,0.839999973773999
    OZn,1,30,8,4,2,-10.1354999542,-16.4332008362,1.08070003986,-3.00589990616,-6.21700000763,-9.19699954987,-1.19400000572,2.54099988937,1.10000002384,0.460000008345,1.54999995232,0.430000007153,2.25,2.22000002861,1.759999960662,0.479999929672
    SZn,1,30,16,4,3,-10.1354999542,-11.7951002121,1.08070003986,-2.84489989281,-6.21700000763,-7.10599994659,-1.19400000572,0.64200001955,1.10000002384,0.740000009537,1.54999995232,0.850000023842,2.25,2.36999988556,1.059999942781,0.559999942785
    SeZn,1,30,34,4,4,-10.1354999542,-10.9460000992,1.08070003986,-2.75099992752,-6.21700000763,-6.65399980545,-1.19400000572,1.31599998474,1.10000002384,0.800000011921,1.54999995232,0.949999988079,2.25,2.18000006676,0.89999997616,0.599999904638
    TeZn,1,30,52,4,5,-10.1354999542,-9.86670017242,1.08070003986,-2.66599988937,-6.21700000763,-6.10900020599,-1.19400000572,0.0989999994636,1.10000002384,0.939999997616,1.54999995232,1.13999998569,2.25,1.83000004292,0.569999992854,0.649999916554

</details>

## Running `SISSO++` for Classification problems
For settings file the only difference between solving classification and regression is the `calc_type` key which is now `classification` instead of `regression` and we are reducing `max_rung` to be 1.
Changing `max_rung` is not a necessary step; however, if rung 2 features are included here then a one-dimensional descriptor will perfectly separate the classes.
Normally this would be a good thing, but because we also want to illustrate two-dimensional visualization tools we will restrict ourselves to a single rung.
Additionally to make it easier to visualize the model we will restrict the calculation to two dimensions, but higher dimensional models are also possible.
```json
{
    "data_file": "data_class.csv",
    "property_key": "Class",
    "desc_dim": 2,
    "n_sis_select": 20,
    "max_rung": 1,
    "calc_type": "classification",
    "min_abs_feat_val": 1e-5,
    "max_abs_feat_val": 1e8,
    "n_residual": 10,
    "n_models_store": 1,
    "leave_out_frac": 0.0,
    "leave_out_inds": [],
    "opset": ["add", "sub", "abs_diff", "mult", "div", "inv", "abs", "exp", "log", "sin", "cos", "sq", "cb", "six_pow", "sqrt", "cbrt", "neg_exp"]
}
```
With this input file and the provided data.csv file we are now able to perform SISSO with the following command
```
mpiexec -n 2 sisso++ sisso.json
```
and get the following on screen output
```
time input_parsing: 0.00104308 s
time to generate feat sapce: 0.00736403 s
Projection time: 0.0013411 s
Time to get best features on rank : 0.000585079 s
Complete final combination/selection from all ranks: 0.000355005 s
Time for SIS: 0.00245714 s
Time for l0-norm: 0.105334 s
Projection time: 0.00138497 s
Time to get best features on rank : 0.00287414 s
Complete final combination/selection from all ranks: 0.000135899 s
Time for SIS: 0.00476503 s
Time for l0-norm: 2.79099 s
Percent of training data in the convex overlap region: 2.43902%
[(r_sigma + r_p_B)]

Percent of training data in the convex overlap region: 0%
[(EA_A * Z_B), (r_sigma + r_p_B)]
```
As with the regression problems, the standard output provides information about what step the calculation just finished and how long it took to complete so you can see where a job failed or ran out of time.
However, the final summary now provides the list of features that best separate out the classes with fewest number of points inside the overlap region of the convex hulls of each class.
The two output files stored in `feature_space/` are also very similar, with the `Score` column in `feature_space/SIS_summary.txt` replaced by the 1D-convex hull overlap between all samples inside the overlap region of the previous model, modified by either how large the overlap region is (positive number) or how sperated the regions are (negative number).
<details>
    <summary>feature_space/SIS_summary.txt</summary>

    # FEAT_ID     Score                   Feature Expression
    0             2.00218777423865069     (r_sigma + r_p_B)
    1             2.0108802733799549      (|r_pi - r_p_A|)
    2             2.0108802733799549      (r_pi - r_p_A)
    3             3.00521883927864941     (r_pi * r_sigma)
    4             6.0271211617331506      (r_sigma / IP_B)
    5             6.02820376741344255     (r_sigma + r_s_B)
    6             6.03528075819598619     (r_sigma + r_p_A)
    7             6.03550443796034486     (r_sigma * r_p_A)
    8             6.04392536619566378     (r_p_A - r_s_B)
    9             6.04660180322806973     (|r_p_A - r_s_B|)
    10            7.03414291492586941     (r_sigma * r_s_A)
    11            8.05674092065169845     (r_pi + r_sigma)
    12            9.05705987306474469     (r_sigma + r_s_A)
    13            10.0641158540893727     (r_sigma / E_HOMO_B)
    14            10.0748914764239395     (r_sigma * r_s_B)
    15            10.1224773571918405     (r_s_A * EA_B)
    16            11.0919350703923243     (r_p_A * E_HOMO_B)
    17            12.0103988563879156     (r_s_A^6)
    18            12.0201816437669606     (r_p_A^6)
    19            12.0631741547198779     (r_p_A * r_s_A)
    #-----------------------------------------------------------------------
    20            -0.999999999964039432   (IP_B^6)
    21            -0.999999999608207957   (Z_A^3)
    22            -0.999999999563166098   (E_HOMO_B^6)
    23            -0.99999999731471545    (period_A^6)
    24            -0.999999995162661137   (Z_B^3)
    25            -0.999999990595021648   (IP_A^6)
    26            -0.999999978575254467   (Z_B * Z_A)
    27            -0.999999973721653945   (EA_B^6)
    28            -0.999999961553741268   (E_HOMO_A^6)
    29            -0.99999991416031242    (IP_B^3)
    30            -0.999999902601353075   (IP_B * Z_A)
    31            -0.999999878198415182   (r_d_A^6)
    32            -0.999999858299492561   (IP_A * Z_A)
    33            -0.99999985529594615    (Z_B / E_LUMO_B)
    34            -0.999999850065982798   (E_HOMO_B * Z_A)
    35            -0.999999771428597528   (period_B * Z_A)
    36            -0.999999756076769386   (E_HOMO_A * Z_A)
    37            -0.999999734902030979   (E_HOMO_B^3)
    38            -0.999999699570055189   (EA_B * Z_A)
    39            -0.99999967830096359    (EA_A * Z_B)
    #-----------------------------------------------------------------------


</details>

Additionally the model files change to better represent the classifier.
The largest changes are in the header, where the coefficients now represent the linear decision boundaries calculated using support-vector machines (SVM).
The estimated property vector in this case refers to the predicted class from SVM
<details>
    <summary>models/train_dim_2_model_0.dat</summary>

    # [(EA_A * Z_B), (r_sigma + r_p_B)]
    # Property Label: $Class$; Unit of the Property: Unitless
    # # Samples in Convex Hull Overlap Region: 0;# Samples SVM Misclassified: 0
    # Decision Boundaries
    # Task     w0                      w1                      b
    # all__0.0_1.0, -1.213391554552964e-02, -1.090594288515569e+01,  2.501183842953395e+01,
    # Feature Rung, Units, and Expressions
    # 0;  1; eV_IP * nuc_charge;                               6|1|mult; (EA_A * Z_B); $\left(EA_{A} Z_{B}\right)$; (EA_A .* Z_B); EA_A,Z_B
    # 1;  1; Unitless;                                         18|15|add; (r_sigma + r_p_B); $\left(r_{sigma} + r_{p, B}\right)$; (r_sigma + r_p_B); r_sigma,r_p_B
    # Number of Samples Per Task
    # Task, n_mats_train
    # all , 82

    # Sample ID , Property Value        ,  Property Value (EST) ,  Feature 0 Value      ,  Feature 1 Value
    AgBr        ,  0.000000000000000e+00,  0.000000000000000e+00, -5.833099961290000e+01,  2.450000047680000e+00
    AgCl        ,  0.000000000000000e+00,  0.000000000000000e+00, -2.833219981198000e+01,  2.520000040527000e+00
    AgF         ,  0.000000000000000e+00,  0.000000000000000e+00, -1.499939990046000e+01,  2.790000051256000e+00
    AgI         ,  1.000000000000000e+00,  1.000000000000000e+00, -8.832979941382000e+01,  2.300000071522000e+00
    AlAs        ,  1.000000000000000e+00,  1.000000000000000e+00, -1.031250000000000e+01,  1.629999995228000e+00
    AlN         ,  1.000000000000000e+00,  1.000000000000000e+00, -2.187500000000000e+00,  1.939999997612000e+00
    AlP         ,  1.000000000000000e+00,  1.000000000000000e+00, -4.687500000000000e+00,  1.650000035759000e+00
    AlSb        ,  1.000000000000000e+00,  1.000000000000000e+00, -1.593750000000000e+01,  1.480000019070000e+00
    AsGa        ,  1.000000000000000e+00,  1.000000000000000e+00, -3.567299902452000e+00,  1.470000028615000e+00
    AsB         ,  1.000000000000000e+00,  1.000000000000000e+00, -3.544200003135000e+00,  1.289999961847000e+00
    BN          ,  1.000000000000000e+00,  1.000000000000000e+00, -7.518000006650001e-01,  1.099999964237000e+00
    BP          ,  1.000000000000000e+00,  1.000000000000000e+00, -1.611000001425000e+00,  1.130000054836000e+00
    BSb         ,  1.000000000000000e+00,  1.000000000000000e+00, -5.477400004845000e+00,  1.820000052445000e+00
    BaO         ,  0.000000000000000e+00,  0.000000000000000e+00,  2.223999977112000e+00,  4.320000201465000e+00
    BaS         ,  0.000000000000000e+00,  0.000000000000000e+00,  4.447999954224000e+00,  4.040000200273000e+00
    BaSe        ,  0.000000000000000e+00,  0.000000000000000e+00,  9.451999902726000e+00,  3.980000197889000e+00
    BaTe        ,  0.000000000000000e+00,  0.000000000000000e+00,  1.445599985122800e+01,  3.840000212194000e+00
    BeO         ,  1.000000000000000e+00,  1.000000000000000e+00,  5.044000148776000e+00,  1.830000072725000e+00
    BeS         ,  1.000000000000000e+00,  1.000000000000000e+00,  1.008800029755200e+01,  1.550000071533000e+00
    BeSe        ,  1.000000000000000e+00,  1.000000000000000e+00,  2.143700063229800e+01,  1.490000069149000e+00
    BeTe        ,  1.000000000000000e+00,  1.000000000000000e+00,  3.278600096704400e+01,  1.350000083454000e+00
    C2          ,  1.000000000000000e+00,  1.000000000000000e+00, -5.234399914740000e+00,  6.299999952320000e-01
    CaO         ,  0.000000000000000e+00,  0.000000000000000e+00,  2.431200027464000e+00,  3.619999915355000e+00
    CaS         ,  0.000000000000000e+00,  0.000000000000000e+00,  4.862400054928000e+00,  3.339999914163000e+00
    CaSe        ,  0.000000000000000e+00,  0.000000000000000e+00,  1.033260011672200e+01,  3.279999911779000e+00
    CaTe        ,  0.000000000000000e+00,  0.000000000000000e+00,  1.580280017851600e+01,  3.139999926084000e+00
    CdO         ,  0.000000000000000e+00,  0.000000000000000e+00,  6.709599971768000e+00,  2.510000020265000e+00
    CdS         ,  1.000000000000000e+00,  1.000000000000000e+00,  1.341919994353600e+01,  2.230000019073000e+00
    CdSe        ,  1.000000000000000e+00,  1.000000000000000e+00,  2.851579988001400e+01,  2.170000016689000e+00
    CdTe        ,  1.000000000000000e+00,  1.000000000000000e+00,  4.361239981649200e+01,  2.030000030994000e+00
    BrCs        ,  0.000000000000000e+00,  0.000000000000000e+00, -1.993599951266000e+01,  4.870000123980000e+00
    ClCs        ,  0.000000000000000e+00,  0.000000000000000e+00, -9.683199763292000e+00,  4.940000116827000e+00
    CsF         ,  0.000000000000000e+00,  0.000000000000000e+00, -5.126399874684000e+00,  5.210000127556000e+00
    CsI         ,  0.000000000000000e+00,  0.000000000000000e+00, -3.018879926202800e+01,  4.720000147822000e+00
    BrCu        ,  1.000000000000000e+00,  1.000000000000000e+00, -5.734749913200000e+01,  2.129999995230000e+00
    ClCu        ,  1.000000000000000e+00,  1.000000000000000e+00, -2.785449957840000e+01,  2.199999988077000e+00
    CuF         ,  0.000000000000000e+00,  0.000000000000000e+00, -1.474649977680000e+01,  2.469999998806000e+00
    CuI         ,  1.000000000000000e+00,  1.000000000000000e+00, -8.684049868560000e+01,  1.980000019072000e+00
    GaN         ,  1.000000000000000e+00,  1.000000000000000e+00, -7.566999793080000e-01,  1.780000030999000e+00
    GaP         ,  1.000000000000000e+00,  1.000000000000000e+00, -1.621499955660000e+00,  1.490000069146000e+00
    GaSb        ,  1.000000000000000e+00,  1.000000000000000e+00, -5.513099849244000e+00,  1.320000052457000e+00
    Ge2         ,  1.000000000000000e+00,  1.000000000000000e+00, -3.036800003052800e+01,  1.159999966620000e+00
    CGe         ,  1.000000000000000e+00,  1.000000000000000e+00, -5.694000005724000e+00,  1.439999997614000e+00
    GeSi        ,  1.000000000000000e+00,  1.000000000000000e+00, -1.328600001335600e+01,  1.139999985693000e+00
    AsIn        ,  1.000000000000000e+00,  1.000000000000000e+00, -8.457900077121000e+00,  1.779999971388000e+00
    InN         ,  1.000000000000000e+00,  1.000000000000000e+00, -1.794100016359000e+00,  2.089999973772000e+00
    InP         ,  1.000000000000000e+00,  1.000000000000000e+00, -3.844500035055000e+00,  1.800000011919000e+00
    InSb        ,  1.000000000000000e+00,  1.000000000000000e+00, -1.307130011918700e+01,  1.629999995230000e+00
    BrK         ,  0.000000000000000e+00,  0.000000000000000e+00, -2.174549937248500e+01,  3.820000171660000e+00
    ClK         ,  0.000000000000000e+00,  0.000000000000000e+00, -1.056209969520700e+01,  3.890000164507000e+00
    FK          ,  0.000000000000000e+00,  0.000000000000000e+00, -5.591699838639000e+00,  4.160000175236000e+00
    IK          ,  0.000000000000000e+00,  0.000000000000000e+00, -3.292889904976300e+01,  3.670000195502000e+00
    BrLi        ,  0.000000000000000e+00,  0.000000000000000e+00, -2.443349897863000e+01,  2.899999976160000e+00
    ClLi        ,  0.000000000000000e+00,  0.000000000000000e+00, -1.186769950390600e+01,  2.969999969007000e+00
    FLi         ,  0.000000000000000e+00,  0.000000000000000e+00, -6.282899737362000e+00,  3.239999979736000e+00
    ILi         ,  0.000000000000000e+00,  0.000000000000000e+00, -3.699929845335400e+01,  2.750000000002000e+00
    MgO         ,  0.000000000000000e+00,  0.000000000000000e+00,  5.539999961856000e+00,  2.770000010735000e+00
    MgS         ,  0.000000000000000e+00,  0.000000000000000e+00,  1.107999992371200e+01,  2.490000009543000e+00
    MgSe        ,  0.000000000000000e+00,  0.000000000000000e+00,  2.354499983788800e+01,  2.430000007159000e+00
    MgTe        ,  0.000000000000000e+00,  0.000000000000000e+00,  3.600999975206400e+01,  2.290000021464000e+00
    BrNa        ,  0.000000000000000e+00,  0.000000000000000e+00, -2.504949897527000e+01,  3.559999942780000e+00
    ClNa        ,  0.000000000000000e+00,  0.000000000000000e+00, -1.216689950227400e+01,  3.629999935627000e+00
    FNa         ,  0.000000000000000e+00,  0.000000000000000e+00, -6.441299736497999e+00,  3.899999946356000e+00
    INa         ,  0.000000000000000e+00,  0.000000000000000e+00, -3.793209844826600e+01,  3.409999966622000e+00
    BrRb        ,  0.000000000000000e+00,  0.000000000000000e+00, -2.066399931907500e+01,  4.690000057220000e+00
    ClRb        ,  0.000000000000000e+00,  0.000000000000000e+00, -1.003679966926500e+01,  4.760000050067000e+00
    FRb         ,  0.000000000000000e+00,  0.000000000000000e+00, -5.313599824904999e+00,  5.030000060796000e+00
    IRb         ,  0.000000000000000e+00,  0.000000000000000e+00, -3.129119896888500e+01,  4.540000081062000e+00
    Si2         ,  1.000000000000000e+00,  1.000000000000000e+00, -1.390199959278200e+01,  1.129999995230000e+00
    CSi         ,  1.000000000000000e+00,  1.000000000000000e+00, -5.957999825478000e+00,  1.430000007151000e+00
    Sn2         ,  1.000000000000000e+00,  1.000000000000000e+00, -5.195999741550001e+01,  1.340000033380000e+00
    CSn         ,  1.000000000000000e+00,  1.000000000000000e+00, -6.235199689860000e+00,  1.759999990465000e+00
    GeSn        ,  1.000000000000000e+00,  1.000000000000000e+00, -3.325439834592000e+01,  1.479999959471000e+00
    SiSn        ,  1.000000000000000e+00,  1.000000000000000e+00, -1.454879927634000e+01,  1.459999978544000e+00
    OSr         ,  0.000000000000000e+00,  0.000000000000000e+00,  2.744800090792000e+00,  3.999999910595000e+00
    SSr         ,  0.000000000000000e+00,  0.000000000000000e+00,  5.489600181584000e+00,  3.719999909403000e+00
    SeSr        ,  0.000000000000000e+00,  0.000000000000000e+00,  1.166540038586600e+01,  3.659999907019000e+00
    SrTe        ,  0.000000000000000e+00,  0.000000000000000e+00,  1.784120059014800e+01,  3.519999921324000e+00
    OZn         ,  1.000000000000000e+00,  1.000000000000000e+00,  8.645600318880000e+00,  2.189999967815000e+00
    SZn         ,  1.000000000000000e+00,  1.000000000000000e+00,  1.729120063776000e+01,  1.909999966623000e+00
    SeZn        ,  1.000000000000000e+00,  1.000000000000000e+00,  3.674380135524000e+01,  1.849999964239000e+00
    TeZn        ,  1.000000000000000e+00,  1.000000000000000e+00,  5.619640207272000e+01,  1.709999978544000e+00


</details>

## Cross-Validation

While we won't do it here, cross-validation should also be performed for classification problems.
For those calculations the number of miscalssified points in the test set is the most important measure of the error.

## Updating the SVM Model Using `sklearn`

The final decision boundary listed in the classification model is found via [linear support vector machine (SVM) model](https://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf).
The objective function for SVM (equation 1 in the linked pdf) balances the size of the margin (distance between the points and the decision boundary) and the number of points miscalssified.
In `libsvm` the trade-off between these two components is controlled by the cost parameter, `c`, where a larger `c` prioritizes the number of miscalssified points over the margin size.
Because the basis of the classification algorithm is based on the overlap region of the convex hull, the `c` value for the SVM model is set at a fairly high value of 1000.0.
While within the scope of the SISSO algorithm this choice makes sense, it may not the best one for all models.

To account for this the python interface has the ability to refit the Linear SVM using the `svm` module of `sklearn`.
Using this functionality we can modify the `c` parameter from 1.0 to 1000.0 on a log-scale and evaluate how well each model is performing.
Importantly we could also use `sklearn` to perform cross-validation on the SVM models to help quantify that performance (we will not do that here since the data is linearly separable and a hard margin is appropriate).
To update the SVM models in python we need to run and store the updated models into separate files as shown below
```python
>>> from sissopp.postprocess.classification import update_model_svm
>>> model_1 = update_model_svm("models/train_dim_2_model_0.dat", 1.0, 100000, filename="models/train_dim_2_model_0_c_1.dat")
The updated coefficient for the decision boundaries:
[array([[ 7.34341190e-04, -8.31373991e-01,  2.06766275e+00]])]
>>> model_10 = update_model_svm("models/train_dim_2_model_0.dat", 10.0, 100000, filename="models/train_dim_2_model_0_c_10.dat")
The updated coefficient for the decision boundaries:
[array([[-2.33513390e-03, -1.83281827e+00,  4.27591385e+00]])]
>>> model_100 = update_model_svm("models/train_dim_2_model_0.dat", 100.0, 100000, filename="models/train_dim_2_model_0_c_100.dat")
The updated coefficient for the decision boundaries:
[array([[-5.91833282e-03, -4.40795646e+00,  1.01687678e+01]])]
>>> model_1000 = update_model_svm("models/train_dim_2_model_0.dat", 1000.0, 100000, filename="models/train_dim_2_model_0_c_1000.dat")
The updated coefficient for the decision boundaries:
[array([[-1.10004264e-02, -9.01866680e+00,  2.07093411e+01]])]
```
Comparing the final `c=1000.0` results to the ones found by `SISSO++` we see that the coefficients for the decision are slightly different.
These changes are a result of different SVM libraries leading to slightly different results; however, if we plot both of these models, we see that the boundaries are fairly close to each other suggesting that the changes are minor.
```python
>>> from sissopp.postprocess.plot.classification import plot_classification
>>> plot_classification("models/train_dim_2_model_0.dat", filename="sissopp.png", fig_settings={"size":{"width": 5.0, "height": 5.0}}).show()
>>> plot_classification("models/train_dim_2_model_0_c_1000.dat", filename="c_1000.png", fig_settings={"size":{"width": 5.0, "height": 5.0}}).show()
```
<details>
<summary> SISSO++ Classification </summary>

![image](./classification/sissopp.png)

</details>

<details>
<summary> sklearn SVM </summary>

![image](./classification/c_1000.png)

</details>

However as we decrease the value of `c` an increasing number of points becomes miss classified, suggesting the model is potentially over-fitting the data and would not properly classify new data points.

```python
>>> from sissopp.postprocess.plot.classification import plot_classification
>>> plot_classification("models/train_dim_2_model_0.dat", filename="sissopp.png", fig_settings={"size":{"width": 5.0, "height": 5.0}}).show()
>>> plot_classification("models/train_dim_2_model_0_c_1.dat", filename="c_1.png", fig_settings={"size":{"width": 5.0, "height": 5.0}}).show()
>>> plot_classification("models/train_dim_2_model_0_c_10.dat", filename="c_10.png", fig_settings={"size":{"width": 5.0, "height": 5.0}}).show()
>>> plot_classification("models/train_dim_2_model_0_c_100.dat", filename="c_100.png", fig_settings={"size":{"width": 5.0, "height": 5.0}}).show()
>>> plot_classification("models/train_dim_2_model_0_c_1000.dat", filename="c_1000.png", fig_settings={"size":{"width": 5.0, "height": 5.0}}).show()
```

<details>
<summary> sklearn SVM c=1.0 </summary>

![image](./classification/c_1.png)

</details>

<details>
<summary> sklearn SVM c=10.0 </summary>

![image](./classification/c_10.png)

</details>

<details>
<summary> sklearn SVM c=100.0 </summary>

![image](./classification/c_100.png)

</details>

<details>
<summary> sklearn SVM c=1000.0 </summary>

![image](./classification/c_1000.png)

</details>


# The `scikit-learn` Interface

To facilitate the integration of `SISSO++` into existing machine learning frameworks we also provide a `scikit-learn` interface via `sissopp.sklearn` submodule.
In this module we wrap all of the `SISSOSolver` objects into classes that use the same syntax as `scikit-learn` models, and can then be used in those frameworks.
Additionally we provide some utility functions to mimic the outputs of other `SISSO++` calculations.

## Using the `scikit-learn` Interface
To use the `scikit-learn` interface we must first set up the `SISSORegressor` object, which we can do either by manually specifying the parameters:

```python
>>> import numpy as np
>>> import pandas as pd
>>> from sissopp.sklearn import SISSORegressor
>>> df = pd.read_csv("data.csv", index_col=0)
>>> sisso = SISSORegressor(
...     prop_label = "E_RS - E_ZB",
...     prop_unit = "eV",
...     allowed_ops = "all",
...     n_dim = 4,
...     max_rung = 2,
...     n_sis_select = 10,
...     n_residual = 10,
... )
```

or by using a pre-existing `Inputs` object

```python
>>> import numpy as np
>>> import pandas as pd
>>> from sissopp.sklearn import SISSORegressor
>>> from sissopp.py_interface import read_csv
>>> 
>>> df = pd.read_csv("data.csv", index_col=0)
>>> inputs = read_csv(df, prop_key="E_RS - E_ZB", max_rung=2)
>>> inputs.n_sis_select = 10
>>> inputs.n_dim = 4
>>> inputs.n_residual = 10
>>> inputs.allowed_ops = ["add", "sub", "abs_diff", "mult", "div", "inv", "abs", "exp", "log", "sin", "cos", "sq", "cb", "six_pow", "sqrt", "cbrt", "neg_exp"]
>>> sisso = SISSORegressor.from_inputs(inputs)
```

From here we can fit a model to the data and evaluate its performance using the `scikit-learn` syntax

```python
>>> X = df.iloc[:, 1:]
>>> y = df.iloc[:, 0].values
>>> sisso.fit(X, y)
time to generate feat space: 10.7571 s
Projection time: 0.35614 s
Time to get best features on rank : 7.20024e-05 s
Complete final combination/selection from all ranks: 0.000281096 s
Time for SIS: 0.507253 s
Time for l0-norm: 0.000794172 s
Projection time: 0.422137 s
Time to get best features on rank : 4.1008e-05 s
Complete final combination/selection from all ranks: 0.000236034 s
Time for SIS: 0.573581 s
Time for l0-norm: 0.00115705 s
Projection time: 0.41528 s
Time to get best features on rank : 3.79086e-05 s
Complete final combination/selection from all ranks: 0.000231028 s
Time for SIS: 0.565092 s
Time for l0-norm: 0.0059731 s
Projection time: 0.412625 s
Time to get best features on rank : 3.88622e-05 s
Complete final combination/selection from all ranks: 0.000246048 s
Time for SIS: 0.565728 s
Time for l0-norm: 0.130348 s
```
and then score it using `scikit-learn` or the default metrics provided by `SISSO++` using the `get_default_model_metric` function
```python
>>> from sissopp.sklearn import regression_metric
>>> from sklearn.metrics import mean_absolute_error
>>> 
>>> y_pred = sisso.predict(X)
>>> rmse = regression_metric(y, y_pred)
>>> mae = mean_absolute_error(y, y_pred)
>>> print(rmse, mae)
0.0556078280028383, 0.044361625085390315
```

Examining the file system we see that this calculation was done in a folder based on the initial timestamp of the claculation: `SISSO_YYYY-mm-dd_hh:mm:ss`.
```bash
ls SISSO_YYYY-mm-dd_hh:mm:ss
feature_space/ models/
```
To change this to a user-defined workdir we need to set the workdir paramter of `sisso`
```python
>>> sisso.workdir = "sklearn_sisso_test"
>>> sisso.fit(X, y)
```
Now we see that SISSO was done inside of `sklearn_sisso_test`
```bash
ls sklearn_sisso_test/
feature_space/ models/
```
If we want to remove this working directory after the calculation is completed, then all we have to do is set `clean_workdir` to `True`
```python
>>> sisso.workdir = "sklearn_sisso_test_cleaned"
>>> sisso.clean_workdir = True
>>> sisso.fit(X, y)
```
and we see that `sklearn_sisso_test` was removed
```bash
ls
data.csv SISSO_YYYY-mm-dd_hh:mm:ss/ sklearn_sisso_test/
```

## Performing Cross-Validation with `scikit-learn`

Now that we can use the interface to fit models using `SISSO` we can now use the cross-validation tools of `scikit-learn`.

```python
>>> from sklearn.model_selection import cross_validate, KFold
>>> from sklearn.metrics import make_scorer
>>> from sissopp.sklearn import get_default_model_metric
>>> 
>>> kf = KFold(n_splits=10, random_state=13, shuffle=True)
>>> 
>>> metric = get_default_model_metric(sisso)
>>> scorer = make_scorer(metric, greater_is_better=False)
>>> 
>>> sisso.workdir = None
>>> sisso.clean_workdir = False
>>> 
>>> cv_results = cross_validate(sisso, X, y, cv=kf, scoring=scorer)
>>> print(cv_results["test_score"].mean())
-0.12127230691149424
```
The test score is negative because `scikit-learn` takes the negative of scores if the problem is trying to minimize a loss.
Looking inside this directory we see that we now have many timestamped directories, and that none of them contain a `test_dim{dd}_model_0.dat` file.
```bash
ls SISSO_YYYY-mm-dd_hh/models/
train_dim_1_model_0.dat  train_dim_2_model_0.dat  train_dim_3_model_0.dat  train_dim_4_model_0.dat
```
Additionally if we were to run these files in a sperate working directory
```python
>>> sisso.workdir = "cv_test"
>>> 
>>> cv_results = cross_validate(sisso, X, y, cv=kf, scoring=scorer)
>>> print(cv_results["test_score"].mean())
-0.12127230691149424
```
thme each cross-validation calculation overwrites the others
```bash
ls cv_test/
feature_space  models
```

To recover the test output files and allow for running of cross-validation insisde seperate working directories, we provide a utility function to get the scores for each split, in a seperate directory and output the results of each model into a test file
```python
>>> from sissopp.sklearn import cross_validate_from_splitter
>>> 
>>> test_scores = cross_validate_from_splitter(X, y, sisso, kf, scoring=scorer)
>>> print(np.mean(test_scores))
-0.12127230691149424
```

Now if we look into the `cv_test/` directory we see
```bash
ls cv_test/
cv_0/  cv_1/  cv_2/  cv_3/  cv_4/
cv_5/  cv_6/  cv_7/  cv_8/  cv_9/
```
with the test models stored inside each `cv_{ii}/models/` directory
```bash
ls cv_test/cv_0/models/
test_dim_1_model_0.dat  test_dim_2_model_0.dat
test_dim_3_model_0.dat  test_dim_4_model_0.dat
train_dim_1_model_0.dat  train_dim_2_model_0.dat
train_dim_3_model_0.dat  train_dim_4_model_0.dat
```
# The Python Interface

## Running `SISSO++` through the python interface
An alternative approach to using the command line is to run `SISSO++` entirely through the python interface.
This tutorial will introduce the user to how to perform SISSO calculations using the python interface, but we will not go over the postprocessing steps demonstrated in [the command line interface tutorial](1_command_line.md).
The biggest advantage to using the python bindings is that you can run and analyze the calculations within the same session; however, for larger jobs on supercomputers the command line interface is likely to be more practical.
Despite this understanding how to use the python interface is important because it allows you to test/demonstrate results in a straightforward manner.

The first step in using the python interface is creating an `Inputs` object that will then be used to construct the `FeatureSpace` and `SISSOSolver` objects.
To construct the `Inputs` object you can use the input file approach as before
```python
>>> from sissopp import Inputs
>>> inputs = Inputs("sisso.json")
```
which will then use the same functions as the command line interface to read the input files.
A second approach would be to read in the data file to create the input object and then set the non-data related properties within a script.
```python
>>> from sissopp.py_interface import read_csv
>>> inputs = read_csv("data.csv", prop_key="E_RS - E_ZB", max_rung=2, leave_out_frac=0.0)
>>> inputs.allowed_ops = [
...     "exp",
...     "neg_exp",
...     "inv",
...     "sq",
...     "cb",
...     "six_pow",
...     "sqrt",
...     "cbrt",
...     "log",
...     "abs",
...     "sin",
...     "cos",
...     "add",
...     "sub",
...     "abs_diff",
...     "mult",
...     "div"
... ]
>>> inputs.n_sis_select = 10
>>> inputs.n_dim = 4
>>> inputs.calc_type = "regression"
>>> inputs.n_residual = 10
>>> inputs.n_models_store = 1
```
Finally if there are no input files an `Inputs` object can be created entirely by setting properties of an empty `Inputs`.
```python
>>> import numpy as np
>>> import pandas as pd
>>> from sissopp import Inputs, Unit, FeatureNode
>>> df = pd.read_csv("data.csv", index_col=0)
>>> inputs = Inputs()
>>> inputs.allowed_ops = [
...     "exp",
...     "neg_exp",
...     "inv",
...     "sq",
...     "cb",
...     "six_pow",
...     "sqrt",
...     "cbrt",
...     "log",
...     "abs",
...     "sin",
...     "cos",
...     "add",
...     "sub",
...     "abs_diff",
...     "mult",
...     "div"
... ]
>>> inputs.n_sis_select = 10
>>> inputs.max_rung = 2
>>> inputs.n_dim = 4
>>> inputs.calc_type = "regression"
>>> inputs.n_residual = 10
>>> inputs.n_models_store = 1
>>> inputs.leave_out_inds = []
>>> inputs.task_names = ["all_mats"]
>>> inputs.task_sizes_train = [82]
>>> inputs.task_sizes_test = [0]
>>> inputs.prop_train = df["E_RS - E_ZB (eV)"].to_numpy()
>>> inputs.prop_test = np.array([])
>>> inputs.prop_label = "E_RS - E_ZB"
>>> inputs.prop_unit = Unit("eV")
>>> inputs.sample_ids_train = df.index.tolist()
>>> inputs.sample_ids_test = []
>>> phi_0 = []
>>> for cc, col in enumerate(df.columns[1:]):
...     expr = col.split("(")[0].strip()
...     if len(col.split("(")) == 2:
...         unit = Unit(col.split("(")[1].split(")")[0].strip())
...     else:
...         unit = Unit()
...     phi_0.append(FeatureNode(cc, expr, df[col].tolist(), [], unit))
...
>>> inputs.phi_0 = phi_0
```
Once `inputs` is created it can then be used to construct the `FeatureSpace` and `SISSOSolver` with
```python
>>> from sissopp.py_interface import get_fs_solver
>>> 
>>> feature_space, sisso = get_fs_solver(inputs, allow_overwrite=False)
```
or with
```python
>>> from sissopp import FeatureSpace, SISSORegressor
>>> 
>>> feature_space = FeatureSpace(inputs)
>>> sisso = SISSORegressor(inputs, feature_space)
```
If you use the `get_fs_solver` there is an additional keword argement `allow_overwrite` that you have to set to `True` if you want to overwrite a previous calculation.
Once created the `SISSOSolver` can then be fit via:
```
>>> sisso.fit()
Projection time: 0.332887 s
Time to get best features on rank : 3.60012e-05 s
Complete final combination/selection from all ranks: 0.000165939 s
Time for SIS: 0.478256 s
Time for l0-norm: 0.0196509 s
Projection time: 0.401822 s
Time to get best features on rank : 6.50883e-05 s
Complete final combination/selection from all ranks: 0.000201225 s
Time for SIS: 0.548018 s
Time for l0-norm: 0.00176501 s
Projection time: 0.407797 s
Time to get best features on rank : 3.79086e-05 s
Complete final combination/selection from all ranks: 0.000196934 s
Time for SIS: 0.559887 s
Time for l0-norm: 0.00716209 s
Projection time: 0.416293 s
Time to get best features on rank : 3.40939e-05 s
Complete final combination/selection from all ranks: 0.000186205 s
Time for SIS: 0.573178 s
Time for l0-norm: 0.154032 s
```
To confirm the models are the same as the one we got from the command line we can access the models directly from `sisso`
```
>>> sisso.models[-1][0]
c0 + a0 * ((E_LUMO_B^3) / sin(period_A)) + a1 * ((r_pi / r_d_A) / (E_HOMO_B^3)) + a2 * (cos(r_sigma) + exp(r_pi)) + a3 * (ln(r_pi) / exp(r_s_A))
```
Additionally the output files generated from the command line will now be present in the current working directory

## Running Cross-Validation Using the python interface
Running cross-validation with the python interface is in principle the same as doing it with [command line interface](1_command_line.md), with the main job being performing multiple calculation with different indexes left out.
To do this we adapt the previous script in the following way
```python
>>> from pathlib import Path
>>> import os
>>> 
>>> from sissopp import Inputs
>>> from sissopp.py_interface import read_csv, get_fs_solver
>>> 
>>> sisso_regs = []
>>> inputs_base = Inputs("sisso.json")
>>> 
>>> for ii in range(100):
...     work_dir = Path(f"cv_{ii:02d}")
...     work_dir.mkdir(exist_ok=True)
...     os.chdir(work_dir)
...     inputs = read_csv("../data.csv", inputs_base.prop_key, inputs=inputs_base, leave_out_frac= 0.05)
...     feature_space, sisso = get_fs_solver(inputs)
...     sisso.fit()
...     sisso_regs.append(sisso)
...     os.chdir("../")
...
```
After the calculations are finished we can then run the same analysis as we did previously using the following
```python
>>> from sissopp.postprocess.check_cv_convergence import jackknife_cv_conv_est
>>> from sissopp.postprocess.plot.cv_error_plot import plot_validation_rmse
>>> import numpy as np
>>> models = np.array([[reg.models[dim][0] for reg in sisso_regs] for dim in range(4)])
>>> mean_val_rmse, var_val_rmse = jackknife_cv_conv_est(models)
>>> print(mean_val_rmse)
[0.19977604 0.14308363 0.09543538 0.09895675]
>>> print(np.sqrt(var_val_rmse))
[0.05247252 0.02463362 0.0322683  0.01908175]
>>> plot_validation_rmse(models, "cv_100_error.pdf").show()
```
It is important to note here that, while it is not necessary to setup separate directories when using the python bindings, if you don't all output files will be overwritten reducing the reproducibility of the code.

## Using the Python Interface to Reproduce Previous Calculations
The next goal of this tutorial will be to discuss how to use the python interface to reproduce previous calculations easily.
As previous examples illustrated how the Model output files can be used to interact and extract information from the models generated by `SISSO++`, but these alone can not be used to recreate a calculation.
To increase the reproducibility of the code `SISSO++` can construct a `FeatureSpace` using a text file, which can be generated using the `phi_out_file` option in the `Inputs` object.
Depending on what we want to study there are two ways of constructing a feature space form either the `phi.txt` file or `selected_features.txt` As show below for using `selected_features.txt`

```python
>>> from sissopp import phi_selected_from_file, FeatureSpace
>>> from sissopp.py_interface import read_csv
>>> 
>>> inputs = read_csv("data.csv", prop_key="E_RS - E_ZB", max_rung=0, leave_out_frac=0.0)
>>> phi_sel = phi_selected_from_file("feature_space/selected_features.txt", inputs.phi_0)
>>> inputs.phi_0 = phi_sel
>>> feature_space = FeatureSpace(inputs)
```

and for using `phi.txt`

```python
>>> from sissopp import phi_selected_from_file, FeatureSpace
>>> from sissopp.py_interface import read_csv
>>> 
>>> inputs = read_csv("data.csv", prop_key="E_RS - E_ZB", max_rung=0, leave_out_frac=0.0)
>>> 
>>> feature_space = FeatureSpace(
>>>     "phi.txt",
>>>     inputs.phi_0,
>>>     inputs.prop_train,
>>>     [82],
>>>     project_type='regression',
>>>     cross_corr_max=1.0,
>>>     n_sis_select=100
>>> )
```

From here calculations can continue as was done in the earlier examples.

## Using the Models to Predict New Values

The previous sections all focused on how to evaluate the models using only the data initially passed to SISSO; however, the models can also be used to predict the target property for a new set of data points.
To do this we need to use the `eval` and `eval_many` functions of the models, after loading in the new dataset.
For this example we will use the same `data.csv` file, but the same procedure can be done with a new data file.

```python
>>> from sissopp.postprocess.load_models import load_model
>>> from sissopp.py_interface.import_dataframe import strip_units
>>> 
>>> model = load_model("models/train_dim_3_model_0.dat")
>>> data_predict = strip_units("data.csv")
>>> y_true = data_predict["E_RS - E_ZB"]
>>> data_predict.drop(columns=["E_RS - E_ZB"], inplace=True)
```
Here `load_model` loads the model into python, and `strip_units` will read the data file and then strip the units out of the column headers leaving only the feature's expression.
We then remove the property column to match the column numbers of the `DataFrame` with the initial `feat_ind` of the primary features.

From here we can then evaluate the model for a single point using the `eval` function

```python
>>> y_pred = model.eval(data_predict.loc["C2", :].values)
>>> print(y_pred)
2.6014256262712045
```

Additionally we can store the data inside a dictionary where the keys are the feature expressions and the values are the data for each feature.

```python
>>> data_dict = {col: data_predict.loc["C2", col] for col in data_predict.columns}
>>> y_pred = model.eval(data_dict)
>>> print(y_pred)
2.6014256262712045
```

To perform the same calculation on multiple data points, we must use the `eval_many` function instead of `eval`

```python
>>> sample_ids = ["C2", "Si2", "Ge2"]
>>> y_pred = model.eval_many(data_predict.loc[sample_ids, :].values)
>>> print(y_pred)
[2.60142563 0.25088082 0.20091549]
>>> 
>>> data_dict = {col: data_predict.loc[sample_ids, col].values for col in data_predict.columns}
>>> y_pred = model.eval_many(data_dict)
>>> print(y_pred)
[2.60142563 0.25088082 0.20091549]
```

Finally, we can store the new predictions into a separate model file, similar to the training/test data using the `prediction_to_file`

```python
>>> model.prediction_to_file(
...     "model_predict/predict_from_dict_dim_3_model_0.dat",
...     y_true[sample_ids],
...     data_dict,
...     sample_ids,
...     [], # Task ID's if applicable
... )
>>> 
>>> model.prediction_to_file(
...     "model_predict/predict_from_arr_dim_3_model_0.dat",
...     y_true[sample_ids],
...     data_predict.loc[sample_ids, :],
...     sample_ids,
...     [], # Task ID's if applicable
... )
```
Resulting in the following output files
<details>
    <summary>The resulting prediction files made from prediction_to_file.</summary>

    # c0 + a0 * ((r_sigma / EA_B) / (r_s_A^6)) + a1 * ((EA_B - IP_A) * (|r_sigma - r_s_B|)) + a2 * ((E_HOMO_B / r_p_A) / (r_sigma + r_p_B))
    # Property Label: $E_{{RS}} - E_{{ZB}}$; Unit of the Property: eV
    # RMSE: 0.0226472214238799; Max AE: 0.0282850060434131
    # Coefficients
    # Task   a0                      a1                      a2                      c0
    # all , -2.221997340952212e-01, -2.740242210520960e-02, -2.264028616363616e-01, -3.129997800655474e-01,
    # Feature Rung, Units, and Expressions
    # 0;  2; AA^-5 * eV_IP^-1;                                 18|7|div|12|sp|div; ((r_sigma / EA_B) / (r_s_A^6)); $\left(\frac{ \left(\frac{ r_{sigma} }{ EA_{B} } \right) }{ \left(r_{s, A}^6\right) } \right)$; ((r_sigma ./ EA_B) ./ (r_s_A).^6); r_sigma,EA_B,r_s_A
    # 1;  2; AA * eV_IP;                                       7|4|sub|18|13|abd|mult; ((EA_B - IP_A) * (|r_sigma - r_s_B|)); $\left(\left(EA_{B} - IP_{A}\right) \left(\left|r_{sigma} - r_{s, B}\right|\right)\right)$; ((EA_B - IP_A) .* abs(r_sigma - r_s_B)); EA_B,IP_A,r_sigma,r_s_B
    # 2;  2; AA^-2 * eV;                                       9|14|div|18|15|add|div; ((E_HOMO_B / r_p_A) / (r_sigma + r_p_B)); $\left(\frac{ \left(\frac{ E_{HOMO, B} }{ r_{p, A} } \right) }{ \left(r_{sigma} + r_{p, B}\right) } \right)$; ((E_HOMO_B ./ r_p_A) ./ (r_sigma + r_p_B)); E_HOMO_B,r_p_A,r_sigma,r_p_B
    # Number of Samples Per Task
    # Task, n_samples
    # all , 3

    # Sample ID , Property Value        ,  Property Value (EST) ,  Feature 0 Value      ,  Feature 1 Value      ,  Feature 2 Value
    C2          ,  2.628603639133640e+00,  2.601425626271205e+00, -0.000000000000000e+00,  6.386751756964515e+00, -1.364575452594942e+01
    Si2         ,  2.791658215483040e-01,  2.508808155048909e-01, -0.000000000000000e+00,  6.358817979658935e+00, -3.260239754057512e+00
    Ge2         ,  2.008525260607710e-01,  2.009154914490295e-01, -0.000000000000000e+00,  6.088560028849320e+00, -3.006837274572150e+00
</details>

# Introduction

This tutorial is based on the [Predicting energy differences between crystal structures: (Meta-)stability of octet-binary compounds](https://analytics-toolkit.nomad-coe.eu/public/user-redirect/notebooks/tutorials/descriptor_role.ipynb) tutorial created by Mohammad-Yasin Arif, Luigi Sbailò, Thomas A. R. Purcell, Luca M. Ghiringhelli, and Matthias Scheffler.
The goal of the tutorial is to teach a user how to use `SISSO++` to find and analyze quantitative models for materials properties.
In particular we will use SISSO to predict the crystal structure (rock-salt or zinc-blende) of a series of octet binaries.
The tutorial will be split into three parts: 1) explaining how to use the executable to perform the calculations and the python utilities to analyze the results and 2) How to use only python to run, analyze, and demonstrate results 3) How to perform classification problems using SISSO.

To learn more about the theory used in this tutorial please look at the following papers:

- [Ouyang, R.; et al. Phys. Rev. Mater., 2(8), 083802.](https://doi.org/10.1103/PhysRevMaterials.2.083802)
- [Ouyang, R.; et al. J. Phys. Mater., 2(2), 024002.](https://doi.org/10.1088/2515-7639/ab077b)

## Outline
The tutorials are split between solving regression problems:

- [Using the Command Line Interface](1_command_line.md)
- [Using the Python Interface](2_python.md)

And classification:

- [Classification](3_classification.md)

For several large and independent calculations using the command-line interface will be best; however, the python interface does provide a good way of performing initial tests and demonstrating the final results.

All tutorials use the octet binary dataset first described in [PRL-2015](https://doi.org/10.1103/PhysRevLett.114.105503) with the goal of predicting whether a material will crystallize in a rock-salt or zinc-blende phase.
For all applications of SISSO a data set has to be passed via a standard `csv` file where the first row represents the feature and property label and the first column are the index-label for each sample for example
```
Material, energy_diff (eV), rs_A (AA), rs_B (AA), E_HOMO_A (eV), E_HOMO_B (eV),....
AgF,-0.153758,1.32,0.41,-4.71,-11.294,...
AgCl,-0.0427973,1.32,0.68,-4.71,-8.7,...
AgBr,-0.0300334,1.32,0.75,-4.71,-8.001,...
AgI,0.0369254,1.32,0.9,-4.71,-7.236,...
...

```
The feature labels have an optional term in () that represents the units of the feature.
If no unit is passed then the feature is assumed to be unitless.

<details>
    <summary>Here is the full data.csv file for the calculation. The features describe the nuclear charge (Z); ionization potential (IP); electron affinity (EA); HOMO and LUMO energies (E_HOMO and E_LUMO); and radii of the atomic s, p and d-orbitals (r_s, r_p, and r_d) of the cation (A) and anion(B) of the materials. Additionally the radii of the \sigma and \pi orbitals of the dimer for each material is included.</summary>

    # Material,E_RS - E_ZB (eV),Z_A (nuc_charge) ,Z_B (nuc_charge) ,period_A,period_B,IP_A (eV_IP) ,IP_B (eV_IP) ,EA_A (eV_IP),EA_B (eV_IP) ,E_HOMO_A (eV) ,E_HOMO_B (eV) ,E_LUMO_A (eV),E_LUMO_B (eV) , r_s_A (AA) , r_s_B (AA) , r_p_A (AA) , r_p_B (AA) , r_d_A (AA) , r_d_B (AA), r_sigma (AA) , r_pi (AA)
    AgBr,-0.030033416711376,47,35,5,4,-8.0580997467,-12.649600029,-1.66659998894,-3.73930001259,-4.71000003815,-8.00100040436,-0.479000002146,0.708000004292,1.32000005245,0.75,1.87999999523,0.879999995232,2.97000002861,1.87000000477,1.570000052448,0.689999938012
    AgCl,-0.042797278205399,47,17,5,3,-8.0580997467,-13.9018001556,-1.66659998894,-3.97079992294,-4.71000003815,-8.69999980927,-0.479000002146,0.574000000954,1.32000005245,0.680000007153,1.87999999523,0.759999990463,2.97000002861,1.66999995708,1.760000050064,0.63999992609
    AgF,-0.153757673178916,47,9,5,2,-8.0580997467,-19.4043006897,-1.66659998894,-4.27349996567,-4.71000003815,-11.2939996719,-0.479000002146,1.25100004673,1.32000005245,0.409999996424,1.87999999523,0.370000004768,2.97000002861,1.42999994755,2.420000046488,0.599999934436
    AgI,0.036925419641193,47,53,5,5,-8.0580997467,-11.2571001053,-1.66659998894,-3.5134999752,-4.71000003815,-7.23600006104,-0.479000002146,0.212999999523,1.32000005245,0.899999976158,1.87999999523,1.07000005245,2.97000002861,1.72000002861,1.230000019072,0.730000019072
    AlAs,0.213261849108676,13,33,3,4,-5.78049993515,-9.26189994812,-0.3125,-1.83920001984,-2.78399991989,-5.34100008011,0.694999992847,0.0640000030398,1.09000003338,0.850000023842,1.38999998569,1.03999996185,1.94000005722,2.01999998093,0.590000033378,0.489999890318
    AlN,0.072949073169639,13,7,3,2,-5.78049993515,-13.5852003098,-0.3125,-1.86749994755,-2.78399991989,-7.2389998436,0.694999992847,3.0569999218,1.09000003338,0.540000021458,1.38999998569,0.509999990463,1.94000005722,1.53999996185,1.430000007149,0.329999983305
    AlP,0.218958341475627,13,15,3,3,-5.78049993515,-9.75059986115,-0.3125,-1.91999995708,-2.78399991989,-5.59600019455,0.694999992847,0.182999998331,1.09000003338,0.829999983311,1.38999998569,0.97000002861,1.94000005722,1.76999998093,0.680000007149,0.439999997609
    AlSb,0.156868733960437,13,51,3,5,-5.78049993515,-8.46829986572,-0.3125,-1.84669995308,-2.78399991989,-4.99100017548,0.694999992847,0.104999996722,1.09000003338,1,1.38999998569,1.23000001907,1.94000005722,2.05999994278,0.25,0.52999997138
    AsGa,0.274277772419737,31,33,4,4,-5.81820011139,-9.26189994812,-0.108099997044,-1.83920001984,-2.73200011253,-5.34100008011,0.129999995232,0.0640000030398,0.990000009537,0.850000023842,1.33000004292,1.03999996185,2.16000008583,2.01999998093,0.430000066765,0.529999971391
    AsB,0.874978183765052,5,33,2,4,-8.18999958038,-9.26189994812,-0.107400000095,-1.83920001984,-3.71499991417,-5.34100008011,2.24799990654,0.0640000030398,0.810000002384,0.850000023842,0.829999983311,1.03999996185,1.95000004768,2.01999998093,0.249999999997,0.209999918935
    BN,1.71208026083627,5,7,2,2,-8.18999958038,-13.5852003098,-0.107400000095,-1.86749994755,-3.71499991417,-7.2389998436,2.24799990654,3.0569999218,0.810000002384,0.540000021458,0.829999983311,0.509999990463,1.95000004768,1.53999996185,0.589999973774,0.050000011922
    BP,1.01922516119521,5,15,2,3,-8.18999958038,-9.75059986115,-0.107400000095,-1.91999995708,-3.71499991417,-5.59600019455,2.24799990654,0.182999998331,0.810000002384,0.829999983311,0.829999983311,0.97000002861,1.95000004768,1.76999998093,0.160000026226,0.160000026226
    BSb,0.580849114368903,5,51,2,5,-8.18999958038,-8.46829986572,-0.107400000095,-1.84669995308,-3.71499991417,-4.99100017548,2.24799990654,0.104999996722,0.810000002384,1,0.829999983311,1.23000001907,1.95000004768,2.05999994278,0.590000033375,0.249999999997
    BaO,-0.092998553867801,56,8,6,2,-5.51569986343,-16.4332008362,0.277999997139,-3.00589990616,-3.34599995613,-9.19699954987,-2.1289999485,2.54099988937,2.15000009537,0.460000008345,2.63000011444,0.430000007153,1.35000002384,2.22000002861,3.890000194312,0.510000020262
    BaS,-0.319762429426191,56,16,6,3,-5.51569986343,-11.7951002121,0.277999997139,-2.84489989281,-3.34599995613,-7.10599994659,-2.1289999485,0.64200001955,2.15000009537,0.740000009537,2.63000011444,0.850000023842,1.35000002384,2.36999988556,3.190000176431,0.590000033375
    BaSe,-0.343445134087233,56,34,6,4,-5.51569986343,-10.9460000992,0.277999997139,-2.75099992752,-3.34599995613,-6.65399980545,-2.1289999485,1.31599998474,2.15000009537,0.800000011921,2.63000011444,0.949999988079,1.35000002384,2.18000006676,3.03000020981,0.629999995228
    BaTe,-0.375386809668271,56,52,6,5,-5.51569986343,-9.86670017242,0.277999997139,-2.66599988937,-3.34599995613,-6.10900020599,-2.1289999485,0.0989999994636,2.15000009537,0.939999997616,2.63000011444,1.13999998569,1.35000002384,1.83000004292,2.700000226504,0.680000007144
    BeO,0.691837577232946,4,8,2,2,-9.459400177,-16.4332008362,0.630500018597,-3.00589990616,-5.59999990463,-9.19699954987,-2.09800004959,2.54099988937,1.08000004292,0.460000008345,1.21000003815,0.430000007153,2.88000011444,2.22000002861,1.400000065572,0.159999996422
    BeS,0.506327674543172,4,16,2,3,-9.459400177,-11.7951002121,0.630500018597,-2.84489989281,-5.59999990463,-7.10599994659,-2.09800004959,0.64200001955,1.08000004292,0.740000009537,1.21000003815,0.850000023842,2.88000011444,2.36999988556,0.700000047691,0.240000009535
    BeSe,0.49494044277526,4,34,2,4,-9.459400177,-10.9460000992,0.630500018597,-2.75099992752,-5.59999990463,-6.65399980545,-2.09800004959,1.31599998474,1.08000004292,0.800000011921,1.21000003815,0.949999988079,2.88000011444,2.18000006676,0.54000008107,0.279999971388
    BeTe,0.468585910493857,4,52,2,5,-9.459400177,-9.86670017242,0.630500018597,-2.66599988937,-5.59999990463,-6.10900020599,-2.09800004959,0.0989999994636,1.08000004292,0.939999997616,1.21000003815,1.13999998569,2.88000011444,1.83000004292,0.210000097764,0.329999983304
    C2,2.62860363913364,6,6,2,2,-10.8516998291,-10.8516998291,-0.87239998579,-0.87239998579,-5.41599988937,-5.41599988937,1.99199998379,1.99199998379,0.639999985695,0.639999985695,0.629999995232,0.629999995232,1.62999999523,1.62999999523,0,0.019999980926
    CaO,-0.265219041319142,20,8,4,2,-6.4279999733,-16.4332008362,0.303900003433,-3.00589990616,-3.86400008202,-9.19699954987,-2.132999897,2.54099988937,1.75999999046,0.460000008345,2.31999993324,0.430000007153,0.680000007153,2.22000002861,3.189999908202,0.589999943972
    CaS,-0.369133194537426,20,16,4,3,-6.4279999733,-11.7951002121,0.303900003433,-2.84489989281,-3.86400008202,-7.10599994659,-2.132999897,0.64200001955,1.75999999046,0.740000009537,2.31999993324,0.850000023842,0.680000007153,2.36999988556,2.489999890321,0.669999957085
    CaSe,-0.360797734421794,20,34,4,4,-6.4279999733,-10.9460000992,0.303900003433,-2.75099992752,-3.86400008202,-6.65399980545,-2.132999897,1.31599998474,1.75999999046,0.800000011921,2.31999993324,0.949999988079,0.680000007153,2.18000006676,2.3299999237,0.709999918938
    CaTe,-0.350456279076752,20,52,4,5,-6.4279999733,-9.86670017242,0.303900003433,-2.66599988937,-3.86400008202,-6.10900020599,-2.132999897,0.0989999994636,1.75999999046,0.939999997616,2.31999993324,1.13999998569,0.680000007153,1.83000004292,1.999999940394,0.759999930854
    CdO,-0.084161358026904,48,8,5,2,-9.5813999176,-16.4332008362,0.838699996471,-3.00589990616,-5.95200014114,-9.19699954987,-1.30900001526,2.54099988937,1.23000001907,0.460000008345,1.74000000954,0.430000007153,2.59999990463,2.22000002861,2.080000013112,0.539999991662
    CdS,0.072672795911785,48,16,5,3,-9.5813999176,-11.7951002121,0.838699996471,-2.84489989281,-5.95200014114,-7.10599994659,-1.30900001526,0.64200001955,1.23000001907,0.740000009537,1.74000000954,0.850000023842,2.59999990463,2.36999988556,1.379999995231,0.620000004775
    CdSe,0.083571949086036,48,34,5,4,-9.5813999176,-10.9460000992,0.838699996471,-2.75099992752,-5.95200014114,-6.65399980545,-1.30900001526,1.31599998474,1.23000001907,0.800000011921,1.74000000954,0.949999988079,2.59999990463,2.18000006676,1.22000002861,0.659999966628
    CdTe,0.114539532194613,48,52,5,5,-9.5813999176,-9.86670017242,0.838699996471,-2.66599988937,-5.95200014114,-6.10900020599,-1.30900001526,0.0989999994636,1.23000001907,0.939999997616,1.74000000954,1.13999998569,2.59999990463,1.83000004292,0.890000045304,0.709999978544
    BrCs,-0.155867302994011,55,35,6,4,-4.00619983673,-12.649600029,-0.569599986076,-3.73930001259,-2.22000002861,-8.00100040436,-0.547999978065,0.708000004292,2.46000003815,0.75,3.16000008583,0.879999995232,1.97000002861,1.87000000477,3.990000128748,0.830000042912
    ClCs,-0.15034615744662,55,17,6,3,-4.00619983673,-13.9018001556,-0.569599986076,-3.97079992294,-2.22000002861,-8.69999980927,-0.547999978065,0.574000000954,2.46000003815,0.680000007153,3.16000008583,0.759999990463,1.97000002861,1.66999995708,4.180000126364,0.78000003099
    CsF,-0.10826331867429,55,9,6,2,-4.00619983673,-19.4043006897,-0.569599986076,-4.27349996567,-2.22000002861,-11.2939996719,-0.547999978065,1.25100004673,2.46000003815,0.409999996424,3.16000008583,0.370000004768,1.97000002861,1.42999994755,4.840000122788,0.740000039336
    CsI,-0.162387474498246,55,53,6,5,-4.00619983673,-11.2571001053,-0.569599986076,-3.5134999752,-2.22000002861,-7.23600006104,-0.547999978065,0.212999999523,2.46000003815,0.899999976158,3.16000008583,1.07000005245,1.97000002861,1.72000002861,3.650000095372,0.870000123972
    BrCu,0.152442639788205,29,35,4,4,-8.38879966736,-12.649600029,-1.6384999752,-3.73930001259,-4.85599994659,-8.00100040436,-0.64099997282,0.708000004292,1.20000004768,0.75,1.67999994755,0.879999995232,2.57999992371,1.87000000477,1.249999999998,0.609999895102
    ClCu,0.156258713192074,29,17,4,3,-8.38879966736,-13.9018001556,-1.6384999752,-3.97079992294,-4.85599994659,-8.69999980927,-0.64099997282,0.574000000954,1.20000004768,0.680000007153,1.67999994755,0.759999990463,2.57999992371,1.66999995708,1.439999997614,0.55999988318
    CuF,-0.017022272342729,29,9,4,2,-8.38879966736,-19.4043006897,-1.6384999752,-4.27349996567,-4.85599994659,-11.2939996719,-0.64099997282,1.25100004673,1.20000004768,0.409999996424,1.67999994755,0.370000004768,2.57999992371,1.42999994755,2.099999994038,0.519999891526
    CuI,0.204674583263113,29,53,4,5,-8.38879966736,-11.2571001053,-1.6384999752,-3.5134999752,-4.85599994659,-7.23600006104,-0.64099997282,0.212999999523,1.20000004768,0.899999976158,1.67999994755,1.07000005245,2.57999992371,1.72000002861,0.909999966622,0.649999976162
    GaN,0.433445239093999,31,7,4,2,-5.81820011139,-13.5852003098,-0.108099997044,-1.86749994755,-2.73200011253,-7.2389998436,0.129999995232,3.0569999218,0.990000009537,0.540000021458,1.33000004292,0.509999990463,2.16000008583,1.53999996185,1.270000040536,0.370000064378
    GaP,0.348751797751902,31,15,4,3,-5.81820011139,-9.75059986115,-0.108099997044,-1.91999995708,-2.73200011253,-5.59600019455,0.129999995232,0.182999998331,0.990000009537,0.829999983311,1.33000004292,0.97000002861,2.16000008583,1.76999998093,0.520000040536,0.480000078682
    GaSb,0.154625285096699,31,51,4,5,-5.81820011139,-8.46829986572,-0.108099997044,-1.84669995308,-2.73200011253,-4.99100017548,0.129999995232,0.104999996722,0.990000009537,1,1.33000004292,1.23000001907,2.16000008583,2.05999994278,0.090000033387,0.570000052453
    Ge2,0.200852526060771,32,32,4,4,-7.56699991226,-7.56699991226,-0.949000000954,-0.949000000954,-4.04600000381,-4.04600000381,2.17499995232,2.17499995232,0.920000016689,0.920000016689,1.15999996662,1.15999996662,2.36999988556,2.36999988556,0.0,0.479999899862
    CGe,0.811442880200048,32,6,4,2,-7.56699991226,-10.8516998291,-0.949000000954,-0.87239998579,-4.04600000381,-5.41599988937,2.17499995232,1.99199998379,0.920000016689,0.639999985695,1.15999996662,0.629999995232,2.36999988556,1.62999999523,0.810000002382,0.249999940394
    GeSi,0.263210170178354,32,14,4,3,-7.56699991226,-7.75769996643,-0.949000000954,-0.992999970913,-4.04600000381,-4.16300010681,2.17499995232,0.439999997616,0.920000016689,0.939999997616,1.15999996662,1.12999999523,2.36999988556,1.88999998569,0.009999990463,0.429999947545
    AsIn,0.134047575193108,49,33,5,4,-5.53739976883,-9.26189994812,-0.256300002337,-1.83920001984,-2.6970000267,-5.34100008011,0.368000000715,0.0640000030398,1.12999999523,0.850000023842,1.5,1.03999996185,3.1099998951,2.01999998093,0.740000009538,0.559999942778
    InN,0.15372029269929,49,7,5,2,-5.53739976883,-13.5852003098,-0.256300002337,-1.86749994755,-2.6970000267,-7.2389998436,0.368000000715,3.0569999218,1.12999999523,0.540000021458,1.5,0.509999990463,3.1099998951,1.53999996185,1.579999983309,0.400000035765
    InP,0.179193287229282,49,15,5,3,-5.53739976883,-9.75059986115,-0.256300002337,-1.91999995708,-2.6970000267,-5.59600019455,0.368000000715,0.182999998331,1.12999999523,0.829999983311,1.5,0.97000002861,3.1099998951,1.76999998093,0.829999983309,0.510000050069
    InSb,0.078059873019811,49,51,5,5,-5.53739976883,-8.46829986572,-0.256300002337,-1.84669995308,-2.6970000267,-4.99100017548,0.368000000715,0.104999996722,1.12999999523,1,1.5,1.23000001907,3.1099998951,2.05999994278,0.39999997616,0.60000002384
    BrK,-0.166175964193826,19,35,4,4,-4.43319988251,-12.649600029,-0.621299982071,-3.73930001259,-2.42600011826,-8.00100040436,-0.697000026703,0.708000004292,2.13000011444,0.75,2.44000005722,0.879999995232,1.78999996185,1.87000000477,2.940000176428,0.439999938012
    ClK,-0.16446068021105,19,17,4,3,-4.43319988251,-13.9018001556,-0.621299982071,-3.97079992294,-2.42600011826,-8.69999980927,-0.697000026703,0.574000000954,2.13000011444,0.680000007153,2.44000005722,0.759999990463,1.78999996185,1.66999995708,3.130000174044,0.38999992609
    FK,-0.146406098498119,19,9,4,2,-4.43319988251,-19.4043006897,-0.621299982071,-4.27349996567,-2.42600011826,-11.2939996719,-0.697000026703,1.25100004673,2.13000011444,0.409999996424,2.44000005722,0.370000004768,1.78999996185,1.42999994755,3.790000170468,0.349999934436
    IK,-0.167039145162562,19,53,4,5,-4.43319988251,-11.2571001053,-0.621299982071,-3.5134999752,-2.42600011826,-7.23600006104,-0.697000026703,0.212999999523,2.13000011444,0.899999976158,2.44000005722,1.07000005245,1.78999996185,1.72000002861,2.600000143052,0.480000019072
    BrLi,-0.032746212884376,3,35,2,4,-5.32910013199,-12.649600029,-0.698099970818,-3.73930001259,-2.87400007248,-8.00100040436,-0.977999985218,0.708000004292,1.64999997616,0.75,2,0.879999995232,6.92999982834,1.87000000477,2.019999980928,0.480000019072
    ClLi,-0.038381482699151,3,17,2,3,-5.32910013199,-13.9018001556,-0.698099970818,-3.97079992294,-2.87400007248,-8.69999980927,-0.977999985218,0.574000000954,1.64999997616,0.680000007153,2,0.759999990463,6.92999982834,1.66999995708,2.209999978544,0.43000000715
    FLi,-0.059488316863735,3,9,2,2,-5.32910013199,-19.4043006897,-0.698099970818,-4.27349996567,-2.87400007248,-11.2939996719,-0.977999985218,1.25100004673,1.64999997616,0.409999996424,2,0.370000004768,6.92999982834,1.42999994755,2.869999974968,0.390000015496
    ILi,-0.021660936341505,3,53,2,5,-5.32910013199,-11.2571001053,-0.698099970818,-3.5134999752,-2.87400007248,-7.23600006104,-0.977999985218,0.212999999523,1.64999997616,0.899999976158,2,1.07000005245,6.92999982834,1.72000002861,1.679999947552,0.520000100132
    MgO,-0.232274724316994,12,8,3,2,-8.03709983826,-16.4332008362,0.692499995232,-3.00589990616,-4.78200006485,-9.19699954987,-1.35800004005,2.54099988937,1.33000004292,0.460000008345,1.89999997616,0.430000007153,3.17000007629,2.22000002861,2.340000003582,0.599999934432
    MgS,-0.086699504988246,12,16,3,3,-8.03709983826,-11.7951002121,0.692499995232,-2.84489989281,-4.78200006485,-7.10599994659,-1.35800004005,0.64200001955,1.33000004292,0.740000009537,1.89999997616,0.850000023842,3.17000007629,2.36999988556,1.639999985701,0.679999947545
    MgSe,-0.055301801956375,12,34,3,4,-8.03709983826,-10.9460000992,0.692499995232,-2.75099992752,-4.78200006485,-6.65399980545,-1.35800004005,1.31599998474,1.33000004292,0.800000011921,1.89999997616,0.949999988079,3.17000007629,2.18000006676,1.48000001908,0.719999909398
    MgTe,-0.004591286648065,12,52,3,5,-8.03709983826,-9.86670017242,0.692499995232,-2.66599988937,-4.78200006485,-6.10900020599,-1.35800004005,0.0989999994636,1.33000004292,0.939999997616,1.89999997616,1.13999998569,3.17000007629,1.83000004292,1.150000035774,0.769999921314
    BrNa,-0.12642872788274,11,35,3,4,-5.22310018539,-12.649600029,-0.715699970722,-3.73930001259,-2.81900000572,-8.00100040436,-0.717999994755,0.708000004292,1.71000003815,0.75,2.59999990463,0.879999995232,6.57000017166,1.87000000477,2.679999947548,1.019999861712
    ClNa,-0.132991985081389,11,17,3,3,-5.22310018539,-13.9018001556,-0.715699970722,-3.97079992294,-2.81900000572,-8.69999980927,-0.717999994755,0.574000000954,1.71000003815,0.680000007153,2.59999990463,0.759999990463,6.57000017166,1.66999995708,2.869999945164,0.96999984979
    FNa,-0.145788137787804,11,9,3,2,-5.22310018539,-19.4043006897,-0.715699970722,-4.27349996567,-2.81900000572,-11.2939996719,-0.717999994755,1.25100004673,1.71000003815,0.409999996424,2.59999990463,0.370000004768,6.57000017166,1.42999994755,3.529999941588,0.929999858136
    INa,-0.114838222187245,11,53,3,5,-5.22310018539,-11.2571001053,-0.715699970722,-3.5134999752,-2.81900000572,-7.23600006104,-0.717999994755,0.212999999523,1.71000003815,0.899999976158,2.59999990463,1.07000005245,6.57000017166,1.72000002861,2.339999914172,1.059999942772
    BrRb,-0.163820531422971,37,35,5,4,-4.28889989853,-12.649600029,-0.590399980545,-3.73930001259,-2.3599998951,-8.00100040436,-0.704999983311,0.708000004292,2.24000000954,0.75,3.20000004768,0.879999995232,1.96000003815,1.87000000477,3.810000061988,1.090000033372
    ClRb,-0.160503554077877,37,17,5,3,-4.28889989853,-13.9018001556,-0.590399980545,-3.97079992294,-2.3599998951,-8.69999980927,-0.704999983311,0.574000000954,2.24000000954,0.680000007153,3.20000004768,0.759999990463,1.96000003815,1.66999995708,4.000000059604,1.04000002145
    FRb,-0.135595776984701,37,9,5,2,-4.28889989853,-19.4043006897,-0.590399980545,-4.27349996567,-2.3599998951,-11.2939996719,-0.704999983311,1.25100004673,2.24000000954,0.409999996424,3.20000004768,0.370000004768,1.96000003815,1.42999994755,4.660000056028,1.000000029796
    IRb,-0.167201442120131,37,53,5,5,-4.28889989853,-11.2571001053,-0.590399980545,-3.5134999752,-2.3599998951,-7.23600006104,-0.704999983311,0.212999999523,2.24000000954,0.899999976158,3.20000004768,1.07000005245,1.96000003815,1.72000002861,3.470000028612,1.130000114432
    Si2,0.279165821548304,14,14,3,3,-7.75769996643,-7.75769996643,-0.992999970913,-0.992999970913,-4.16300010681,-4.16300010681,0.439999997616,0.439999997616,0.939999997616,0.939999997616,1.12999999523,1.12999999523,1.88999998569,1.88999998569,0.0,0.379999995228
    CSi,0.669023727235981,14,6,3,2,-7.75769996643,-10.8516998291,-0.992999970913,-0.87239998579,-4.16300010681,-5.41599988937,0.439999997616,1.99199998379,0.939999997616,0.639999985695,1.12999999523,0.629999995232,1.88999998569,1.62999999523,0.800000011919,0.199999988077
    Sn2,0.016963899193797,50,50,5,5,-7.04279994965,-7.04279994965,-1.03919994831,-1.03919994831,-3.86599993706,-3.86599993706,0.00800000037998,0.00800000037998,1.05999994278,1.05999994278,1.34000003338,1.34000003338,2.02999997139,2.02999997139,0.0,0.5600001812
    CSn,0.453537974142819,50,6,5,2,-7.04279994965,-10.8516998291,-1.03919994831,-0.87239998579,-3.86599993706,-5.41599988937,0.00800000037998,1.99199998379,1.05999994278,0.639999985695,1.34000003338,0.629999995232,2.02999997139,1.62999999523,1.129999995233,0.290000081063
    GeSn,0.081663360237144,50,32,5,4,-7.04279994965,-7.56699991226,-1.03919994831,-0.949000000954,-3.86599993706,-4.04600000381,0.00800000037998,2.17499995232,1.05999994278,0.920000016689,1.34000003338,1.15999996662,2.02999997139,2.36999988556,0.319999992851,0.520000040531
    SiSn,0.135108799106092,50,14,5,3,-7.04279994965,-7.75769996643,-1.03919994831,-0.992999970913,-3.86599993706,-4.16300010681,0.00800000037998,0.439999997616,1.05999994278,0.939999997616,1.34000003338,1.12999999523,2.02999997139,1.88999998569,0.329999983314,0.470000088214
    OSr,-0.22030662317411,38,8,5,2,-6.03159999847,-16.4332008362,0.343100011349,-3.00589990616,-3.64100003242,-9.19699954987,-1.3789999485,2.54099988937,1.90999996662,0.460000008345,2.54999995232,0.430000007153,1.20000004768,2.22000002861,3.569999903442,0.669999986892
    SSr,-0.368434129930392,38,16,5,3,-6.03159999847,-11.7951002121,0.343100011349,-2.84489989281,-3.64100003242,-7.10599994659,-1.3789999485,0.64200001955,1.90999996662,0.740000009537,2.54999995232,0.850000023842,1.20000004768,2.36999988556,2.869999885561,0.750000000005
    SeSr,-0.3745109517331,38,34,5,4,-6.03159999847,-10.9460000992,0.343100011349,-2.75099992752,-3.64100003242,-6.65399980545,-1.3789999485,1.31599998474,1.90999996662,0.800000011921,2.54999995232,0.949999988079,1.20000004768,2.18000006676,2.70999991894,0.789999961858
    SrTe,-0.379294725862565,38,52,5,5,-6.03159999847,-9.86670017242,0.343100011349,-2.66599988937,-3.64100003242,-6.10900020599,-1.3789999485,0.0989999994636,1.90999996662,0.939999997616,2.54999995232,1.13999998569,1.20000004768,1.83000004292,2.379999935634,0.839999973773999
    OZn,0.101968176768423,30,8,4,2,-10.1354999542,-16.4332008362,1.08070003986,-3.00589990616,-6.21700000763,-9.19699954987,-1.19400000572,2.54099988937,1.10000002384,0.460000008345,1.54999995232,0.430000007153,2.25,2.22000002861,1.759999960662,0.479999929672
    SZn,0.275813325606578,30,16,4,3,-10.1354999542,-11.7951002121,1.08070003986,-2.84489989281,-6.21700000763,-7.10599994659,-1.19400000572,0.64200001955,1.10000002384,0.740000009537,1.54999995232,0.850000023842,2.25,2.36999988556,1.059999942781,0.559999942785
    SeZn,0.263136899280653,30,34,4,4,-10.1354999542,-10.9460000992,1.08070003986,-2.75099992752,-6.21700000763,-6.65399980545,-1.19400000572,1.31599998474,1.10000002384,0.800000011921,1.54999995232,0.949999988079,2.25,2.18000006676,0.89999997616,0.599999904638
    TeZn,0.245001295174006,30,52,4,5,-10.1354999542,-9.86670017242,1.08070003986,-2.66599988937,-6.21700000763,-6.10900020599,-1.19400000572,0.0989999994636,1.10000002384,0.939999997616,1.54999995232,1.13999998569,2.25,1.83000004292,0.569999992854,0.649999916554
</details>
# Acknowledgements

`SISSO++` would not be possible without the following packages:

- [boost](https://www.boost.org/)
- [Coin-Clp](https://github.com/coin-or/Clp)
- [LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/)
- [googletest](https://github.com/google/googletest)
- [NLopt](http://github.com/stevengj/nlopt)

## How to cite these packages:

Please make sure to give credit to the right people when using `SISSO++`:
For classification problems cite:
- [How to cite Coin-Clp](https://zenodo.org/record/3748677#.YBuxDVmYVhE)
- [How to cite LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html#f203)

For use of the parameterized nodes please cite:
- [How to cite NLopt](https://nlopt.readthedocs.io/en/latest/Citing_NLopt/)
# Work that was performed using `SISSO++`

# Contributing Guidelines

When contributing to this repository, please first discuss the change you wish to make via an issue, email, or any other method with the maintainers of this repository.
This will make life easier for everyone.

## Report Issues

Please use the [issue tracker](https://gitlab.com/sissopp-developers/sissopp/-/issues) to report issues. Before posting an issue please insure that it meets the following requirements:

- The issue has not been reported previously (Have a brief look at the issues page)
- Describe the issue in terms of actual v. expected behavior
- Provide a minimal example of the issue you are seeing


## Contribute Code via Merge Request

In order to contribute code to `SISSO++`, please use a merge request (see guidelines of preparing a merge request [here](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html)).

- Please _document_ and _test_ your changes. Tests are found in `sissopp/tests` and written with [pytest](https://docs.pytest.org/en/stable/) for the python bindings and [googletest](https://github.com/google/googletest) for the C++ interface.
- If a new feature is introduced please create a minimal test of the binary file in `exec_tests`
# Input Files

To see a sample of the input files look in `/path/to/sisso++/tests/exec_test`.
In this directory there are multiple subdirectories for different types of calculations, but the `default/` directory would be the most common application.

To use the code two files are necessary: `sisso.json` and `data.csv`.
`data.csv` stores all of the data for the calculation in a `csv` file.
The first row in the file corresponds to the feature meta data with the following format `expression (Unit)`.
For example if one of the primary features used in the set is the lattice constant of a material the header would be `lat_param (AA)`.
The first column of the file are sample labels for all of the other rows, and is used to set the sample ids in the output files.

The input parameters are stored in `sisso.json`, here is a list of all possible variables that can be set in `sisso.json`

### `data.csv`
The data file contains all relevant data and metadata to describe the individual features and samples.
The first row of the file corresponds to the features metadata and has the following format `expression (Unit)` or `expression`.
For the cases where no `(Unit)` is included in the header then the feature is considered to be dimensionless.
For example if one of the primary features used in the set is the lattice constant of a material the header would be `lat_param (AA)`, but the number of species in the material would be `n_species` because it is a dimensionless number.


The first column provide the labels for each sample in the data file, and is used to set the sample ids in the output files.
In the simplest case, this can be just a running index.
The data describing the property vector is defined in the column with an `expression` matching the `property_key` filed in the `sisso.json` file, and will not be included in the feature space.
Additionally, an optional `Task` column whose header matches the `task_key` field in the sisso.json file can also be included in the data file.
This column maps each sample to a respective task with a label defined in the task column.
Below in a minimal example of the data file used to learn a model for a materials volume.

```
material, Structure_Type, Volume (AA^3), lat_param (AA)
C, diamond, 45.64, 3.57
Si, diamond, 163.55, 5.47
Ge, diamond, 191.39, 5.76
Sn, diamond, 293.58, 6.65
Pb, diamond, 353.84, 7.07.757
LiF, rock_salt, 67.94, 4.08
NaF, rock_salt, 103.39, 4.69
KF, rock_salt, 159.00, 5.42
RbF, rock_salt, 189.01, 5.74
CsF, rock_salt, 228.33, 6.11
```

### `sisso.json`

All input parameters that can not be extracted from the data file are defined in the `sisso.json` file.

Here is a complete example of a `sisso.json` file where the property and task keys match those in the above data file example.

```json
{
    "data_file": "data.csv",
    "property_key": "Volume",
    "task_key": "Structure_Type",
    "opset": ["add", "sub", "mult", "div", "sq", "cb", "cbrt", "sqrt"],
    "param_opset": [],
    "calc_type": "regression",
    "desc_dim": 2,
    "n_sis_select": 5,
    "max_rung": 2,
    "n_residual": 1,
    "n_models_store": 1,
    "n_rung_store": 1,
    "n_rung_generate": 0,
    "min_abs_feat_val": 1e-5,
    "max_abs_feat_val": 1e8,
    "leave_out_inds": [0, 5],
    "leave_out_frac": 0.25,
    "fix_intercept": false,
    "max_feat_cross_correlation": 1.0,
    "nlopt_seed": 13,
    "global_param_opt": false,
    "reparam_residual": true
}
```

A description of all fields is listed below. Anything highlighted in red should be in all `sisso.json` files, while anything highlighted in blue highlighted terms should be defined depending on the circumstances. All terms listed in black can normally remain as their default values.

#### <span style="color:red">data_file</span>

*Default: "data.csv"*

The name of the csv file where the data is stored.

#### <span style="color:red">property_key</span>

*Default: "prop"*

The expression of the column where the property to be modeled is stored.

#### <span style="color:red">task_key</span>

*Default: "task"*

The expression of the column where the task identification is stored.

#### <span style="color:red">opset</span>

A list containing the set of all operators that will be used during the feature creation step of SISSO. (If empty use all available features)

#### <span style="color:blue">param_opset</span>

A list containing the set of all operators, for which the non-linear scale and bias terms will be optimized, that will be used during the feature creation step of SISSO. (If empty none of the available features are used)

#### <span style="color:red">calc_type</span>

*Default: "regression"*

The type of calculation to run either regression, log regression, or classification

#### <span style="color:red">desc_dim</span>

The maximum dimension of the model to be created (no default value)

#### <span style="color:red">n_sis_select</span>

The number of features that SIS selects over each iteration (no default value)

#### <span style="color:red">max_rung</span>

The maximum rung of the feature (height of the tallest possible binary expression tree - 1) (no default value)

#### <span style="color:red">n_residual</span>

*Default: 1*

Number of residuals to used to select the next subset of materials in the iteration. (Affects SIS after the 1D model)

#### n_models_store

*Default: `n_residual`*

Number of models to output as file for each dimension

#### n_rung_store

*Default: `max_rung` - 1*

The number of rungs where all of the training/testing data of the materials are stored in memory.

#### <span style="color:blue">n_rung_generate</span>

*Default: 0*

The number of rungs to generate on the fly during each SIS step. Must be 1 or 0.

#### min_abs_feat_val

*Default: 1e-50*

Minimum absolute value allowed in the feature's training data

#### max_abs_feat_val

*Default: 1e50*

Maximum absolute value allowed in the feature's training data

#### leave_out_inds

The list of indexes from the data set to use as the test set. If empty and `leave_out_frac > 0` the selection will be random

#### <span style="color:blue">leave_out_frac</span>

*Default: 0.0*

Fraction (in decimal form) of the data to use as a test set. This is not used if `leave_out_inds` is set.

#### <span style="color:blue">fix_intercept</span>

*Default: false*

If true set the bias term for regression models to 0.0. For classification problems this must be

#### <span style="color:blue">max_feat_cross_correlation</span>

*Default: 1.0*

The maximum Pearson correlation allowed between selected features

#### nlopt_seed

*Default: 42*

The random seed used for seeding the pseudo-random number generator for NLopt

#### <span style="color:blue">global_param_opt</span>

*Default: false*

If true then attempt to globally optimize the non-linear scale/bias terms for the operators in `param_opset`

#### <span style="color:blue">reparam_residual</span>

*Default: false*

If true then reparameterize features based on the residuals

# Running the Code

## Perform the Calculation
Once the input files are made the code can be run using the following command

```
mpiexec -n 2 /path/to/sisso++/bin/sisso++ sisso.json
```

which will give the following output for the simple problem defined above

```text
time input_parsing: 0.000721931 s
time to generate feat sapce: 0.00288105 s
Projection time: 0.00304198 s
Time to get best features on rank : 1.09673e-05 s
Complete final combination/selection from all ranks: 0.00282502 s
Time for SIS: 0.00595999 s
Time for l0-norm: 0.00260496 s
Projection time: 0.000118971 s
Time to get best features on rank : 1.38283e-05 s
Complete final combination/selection from all ranks: 0.00240111 s
Time for SIS: 0.00276804 s
Time for l0-norm: 0.000256062 s
Train RMSE: 0.293788 AA^3; Test RMSE: 0.186616 AA^3
c0 + a0 * (lat_param^3)

Train RMSE: 0.0936332 AA^3; Test RMSE: 15.8298 AA^3
c0 + a0 * ((lat_param^3)^2) + a1 * (sqrt(lat_param)^3)

```
## Analyzing the Results

Once the calculations are done, two sets of output files are generated.
Two files that summarize the results from SIS in a computer and human readable manner are stored in: `feature_space/` and every model used as a residual for SIS is stored in `models/`.
The human readable file describing the selected feature space is `feature_space/SIS_summary.txt` which contains the projection score (The Pearson correlation to the target property or model residual).
```
# FEAT_ID     Score                   Feature Expression
0             0.99997909235669924     (lat_param^3)
1             0.999036700010245471    ((lat_param^2)^2)
2             0.998534266139345261    (lat_param^2)
3             0.996929900301868899    (sqrt(lat_param)^3)
4             0.994755117666830335    lat_param
#-----------------------------------------------------------------------
5             0.0318376000648976157   ((lat_param^3)^3)
6             0.00846237838476477863  ((lat_param^3)^2)
7             0.00742498801557322716  cbrt(cbrt(lat_param))
8             0.00715447033658055554  cbrt(sqrt(lat_param))
9             0.00675695980092700429  sqrt(sqrt(lat_param))
#---------------------------------------------------------------------
```
The computer readable file file is `feature_space/selected_features.txt` and contains a the list of selected features represented by an alphanumeric code where the integers are the index of the feature in the primary feature space and strings represent the operators.
The order of each term in these expressions is the same as the order it would appear using postfix (reverse polish) notation.
```
# FEAT_ID     Feature Postfix Expression (RPN)
0             0|cb
1             0|sq|sq
2             0|sq
3             0|sqrt|cb
4             0
#-----------------------------------------------------------------------
5             0|cb|cb
6             0|cb|sq
7             0|cbrt|cbrt
8             0|sqrt|cbrt
9             0|sqrt|sqrt
#-----------------------------------------------------------------------
```

The model output files are split into train/test files sorted by the dimensionality of the model and by the train RMSE.
The model with the lowest RMSE is stored in the lowest numbered file.
For example `train_dim_2_model_0.dat` will have the best 2D model, `train_dim_2_model_1.dat` would have the second best, etc., whereas `train_dim_1_model_0.dat` will have the best 1D model.
Each model file has a large header containing information about the features selected and model generated
```
# c0 + a0 * (lat_param^3)
# Property Label: $Volume$; Unit of the Property: AA^3
# RMSE: 0.293787533962641; Max AE: 0.56084644346538
# Coefficients
# Task       a0                      c0
#  diamond,  1.000735616997855e+00, -1.551085274074442e-01,
#  rock_salt,  9.998140372873336e-01,  6.405707194855371e-02,
# Feature Rung, Units, and Expressions
# 0;  1; AA^3;                                             0|cb; (lat_param^3); $\left(lat_{param}^3\right)$; (lat_param).^3; lat_param
# Number of Samples Per Task
# Task    , n_mats_train
#  diamond, 4
#  rock_salt, 4
```
The first section of the header summarizes the model by providing a string representation of the model, defines the property's label and unit, and summarizes the error of the model.
```
# c0 + a0 * (lat_param^3)
# Property Label: $Volume$; Unit of the Property: AA^3
# RMSE: 0.293787533962641; Max AE: 0.56084644346538
```
Next the linear coefficients (as shown in the first line) for each task is listed.
```
# Coefficients
# Task       a0                      c0
#  diamond,  1.000735616997855e+00, -1.551085274074442e-01,
#  rock_salt,  9.998140372873336e-01,  6.405707194855371e-02,
```
Then a description of each feature in the model is listed, including units and various expressions.
```
# Feature Rung, Units, and Expressions
# 0;  1; AA^3;                                             0|cb; (lat_param^3); $\left(lat_{param}^3\right)$; (lat_param).^3; lat_param
```
Finally information about the number of samples in each task is given
```
# Number of Samples Per Task
# Task    , n_mats_train
#  diamond, 4
#  rock_salt, 4
```

The header for the test data files contain the same information as the training file, with an additional line at the end to list all indexes included in the test set:
```
# Test Indexes: [ 0, 5 ]
```
These indexes can be used to reproduce the results by setting `leave_out_inds` to those listed on this line.

After this header in both file the following data is stored in the file:

```
# Sample ID , Property Value        ,  Property Value (EST) ,  Feature 0 Value
```
With this data, one can plot and analyzed the model, e.g., by using the python binding.



## Using the Python Library
To see how the python interface can be used refer to the [tutorials](../tutorial/2_python.md).
If you get an error about not being able to load MKL libraries, you may have to run `conda install numpy` to get proper linking.

Installation
---
The package uses a CMake build system, and compatible all versions of the C++ standard library after C++ 14.
You can access the code [here](http://gitlab.com/sissopp_developers/sissopp).
To download the code including the `pybind11` submodule (necessary for the python bindings) then you need to run:
```bash
git clone --recursive http://gitlab.com/sissopp_developers/sissopp
```
This will clone both the source code and `pybind11` dependency to your machine.

### Prerequisites
To install `SISSO++` the following packages are needed:

- CMake version 3.10 and up
- A C++ compiler (compatible with C++ 14 and later, e.g. gcc 5.0+ or icpc 17.0+)
- BLAS/LAPACK
- MPI

Additionally the following packages needed by SISSO++ will be installed (if they are not installed already/if they cannot be found in $PATH)

  - [Boost](https://www.boost.org) (mpi, serialization, system, and filesystem libraries)
  - [pybind11](https://pybind11.readthedocs.io/en/stable/index.html) (Provided as a git submodule)
  - [GTest](https://github.com/google/googletest)
  - [Coin-Clp](https://github.com/coin-or/Clp)
  - [NLopt](https://github.com/stevengj/nlopt)
  - [{fmt}](https://fmt.dev/latest/index.html) (Used for the C++ 20 [std::format](https://en.cppreference.com/w/cpp/utility/format/format) library)

To build and use the optional python bindings the following are also needed:

- [Python 3.6 or greater](https://www.python.org/)
- [numpy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [scipy](https://www.scipy.org/)
- [seaborn](https://seaborn.pydata.org/)
- [sklearn](https://scikit-learn.org/stable/index.html)
- [toml](https://pypi.org/project/toml/)

The setup of the python environment can be done using anaconda with

```bash
conda create -n sissopp_env python=3.9 numpy pandas scipy seaborn scikit-learn toml
```

### Installing `SISSO++`

`SISSO++` is installed using a cmake build system, with sample configuration files located in `cmake/toolchains/`
For example, here is `initial_config.cmake` file used to construct `SISSO++` and the python bindings using the gnu compiler.

```
###############
# Basic Flags #
###############
set(CMAKE_CXX_COMPILER g++ CACHE STRING "")
set(CMAKE_C_COMPILER gcc CACHE STRING "")
set(CMAKE_CXX_FLAGS "-O3 -march=native" CACHE STRING "")
set(CMAKE_C_FLAGS "-O3 -march=native" CACHE STRING "")

#################
# Feature Flags #
#################
set(BUILD_PYTHON ON CACHE BOOL "")
set(BUILD_PARAMS ON CACHE BOOL "")
```

Because we want to build with the python bindings in this example and assuming there is no preexisting python environment, we need to first create/activate it.
For this example we will use `conda`, but standard python installations or virtual environments are also possible.

```bash
conda create -n sissopp_env python=3.9 numpy pandas scipy seaborn scikit-learn toml
conda activate sissopp_env
```

Note if you are using a python environment with a local MKL installation then make sure the versions of all accessible MKL libraries are the same.

Now we can install `SISSO++` using `initial_config.cmake` and the following commands (this assumes gnu compiler and MKL are used, if you are using a different compiler/BLAS library change the flags to the relevant directories)

```bash
export MKLROOT=/path/to/mkl/
export BOOST_ROOT=/path/to/boost

cd ~/sissopp/
mkdir build/;
cd build/;

cmake -C initial_config.cmake ../
make
make install
```

Once all the commands are run `SISSO++` should be in the `~/SISSO++/main directory/bin/` directory.

#### Installing the Python Bindings Without Administrative Privileges

To install the python bindings on a machine where you do not have write privilege to the default python install directory (typical on most HPC systems), you must set the `PYTHON_INSTDIR` to a directory where you do have write access.
This can be done by modifying the `camke` command to:

```bash
cmake -C initial_config.cmake -DPYTHON_INSTDIR=/path/to/python/install/directory/ ../
```

A standard local python installation directory for pip and conda is `$HOME/.local/lib/python3.X/site-packages/` where X is the minor version of python.
It is important that if you do set this variable then that directory is also inside your `PYTHONPATH` envrionment variable. This can be updated with

```bash
export PYTHONPATH=$PYTHONPATH:/path/to/python/install/directory/
```

If you are using anaconda, then this can be avoided by creating a new conda environment as detailed above.

You will need to set this variable and recompile the code (remove all build files first) if you see this error

```bash

CMake Error at src/cmake_install.cmake:114 (file):
  file cannot create directory:
  ${PTYHON_BASE_DIR}/lib/python3.X/site-packages/sissopp.
  Maybe need administrative privileges.
Call Stack (most recent call first):
  cmake_install.cmake:42 (include)" 
```

#### Install the Binary Without the Python Bindings

To install only the `SISSO++` executable repeat the same commands as above but set `USE_PYTHON` in `initial_config.cmake` to `OFF`.

#### Testing the Compilation

To test the compilation of `SISSO++` run `make test` and you should get the following output or a subset of it depending on which parameters you are using
```
Test project /home/purcell/git/cpp_sisso/build
      Start  1: Classification
 1/17 Test  #1: Classification ...........................................   Passed    1.74 sec
      Start  2: Classification_Max_Correlation_NE_One
 2/17 Test  #2: Classification_Max_Correlation_NE_One ....................   Passed    0.55 sec
      Start  3: Classification_Generate_Project
 3/17 Test  #3: Classification_Generate_Project ..........................   Passed    0.55 sec
      Start  4: Classification_Max_Correlation_NE_One_Generate_Project
 4/17 Test  #4: Classification_Max_Correlation_NE_One_Generate_Project ...   Passed    0.54 sec
      Start  5: Regression
 5/17 Test  #5: Regression ...............................................   Passed    0.38 sec
      Start  6: Generate_Project
 6/17 Test  #6: Generate_Project .........................................   Passed    0.38 sec
      Start  7: Maximum_Correlation_NE_One
 7/17 Test  #7: Maximum_Correlation_NE_One ...............................   Passed    0.43 sec
      Start  8: Maximum_Correlation_NE_One_Generate_Project
 8/17 Test  #8: Maximum_Correlation_NE_One_Generate_Project ..............   Passed    0.39 sec
      Start  9: Log_Regression
 9/17 Test  #9: Log_Regression ...........................................   Passed    0.37 sec
      Start 10: Log_Regression_Max_Correlation_NE_One
10/17 Test #10: Log_Regression_Max_Correlation_NE_One ....................   Passed    0.37 sec
      Start 11: Log_Regression_Generate_Project
11/17 Test #11: Log_Regression_Generate_Project ..........................   Passed    0.39 sec
      Start 12: Log_Regression_Max_Correlation_NE_One_Generate_Project
12/17 Test #12: Log_Regression_Max_Correlation_NE_One_Generate_Project ...   Passed    0.40 sec
      Start 13: NL_Optimization
13/17 Test #13: NL_Optimization ..........................................   Passed    3.06 sec
      Start 14: NL_Optimization_Residuals
14/17 Test #14: NL_Optimization_Residuals ................................   Passed    9.66 sec
      Start 15: NL_Optimization_Residuals_Generate_Project
15/17 Test #15: NL_Optimization_Residuals_Generate_Project ...............   Passed    3.04 sec
      Start 16: Python_Bindings
16/17 Test #16: Python_Bindings ..........................................   Passed   12.49 sec
      Start 17: CPP_Unit_Tests
17/17 Test #17: CPP_Unit_Tests ...........................................   Passed   12.13 sec

100% tests passed, 0 tests failed out of 17

Total Test time (real) =  46.87 sec
```

Note the total test time may vary based on optimizations and computational resources used

#### Cleaning up Installation Files

To remove all previously built library and executable files run `make clean` in the build directory. If a compilation error occurs because of issues with one of the external libraries, then clean the build and recompile.

### Installing Only the Python Bindings with pip

If one only wants to use the python bindings without the executable, you can also install `SISSO++` using pip. To do this simply run

```bash
pip install .
```

inside the cloned SISSO repository. If you do not have administrative priveldges then you need to use

```bash
pip install --user .
```

This command will build all external libraries and `SISSO++` and store them inside a seperate `_sissopp` module that will be imported into `sissopp` automatically. To then test the code run
```bash
pytest test/pytest
```
And you should see the following output
```
========================================== test session starts ==========================================
platform linux -- Python 3.9.7, pytest-6.2.5, py-1.10.0, pluggy-1.0.0
rootdir: /home/purcell/git/sissopp
collected 83 items

tests/pytest/test_param.py .                                                                      [  1%]
tests/pytest/test_sisso.py .                                                                      [  2%]
tests/pytest/test_descriptor_identifier/test_class_model_from_file.py .                           [  3%]
tests/pytest/test_descriptor_identifier/test_class_model_retrain_svm.py .                         [  4%]
tests/pytest/test_descriptor_identifier/test_class_model_train_from_file.py .                     [  6%]
tests/pytest/test_descriptor_identifier/test_classifier.py .                                      [  7%]
tests/pytest/test_descriptor_identifier/test_log_reg_model_from_file.py .                         [  8%]
tests/pytest/test_descriptor_identifier/test_log_reg_train_model_from_file.py .                   [  9%]
tests/pytest/test_descriptor_identifier/test_log_regressor.py .                                   [ 10%]
tests/pytest/test_descriptor_identifier/test_reg_model_from_file.py .                             [ 12%]
tests/pytest/test_descriptor_identifier/test_reg_model_train_from_file.py .                       [ 13%]
tests/pytest/test_descriptor_identifier/test_regressor.py .                                       [ 14%]
tests/pytest/test_feature_creation/test_feat_generation/test_abs_diff_node.py .                   [ 15%]
tests/pytest/test_feature_creation/test_feat_generation/test_abs_node.py .                        [ 16%]
tests/pytest/test_feature_creation/test_feat_generation/test_add_node.py .                        [ 18%]
tests/pytest/test_feature_creation/test_feat_generation/test_cb_node.py .                         [ 19%]
tests/pytest/test_feature_creation/test_feat_generation/test_cbrt_node.py .                       [ 20%]
tests/pytest/test_feature_creation/test_feat_generation/test_cos_node.py .                        [ 21%]
tests/pytest/test_feature_creation/test_feat_generation/test_div_node.py .                        [ 22%]
tests/pytest/test_feature_creation/test_feat_generation/test_exp_node.py .                        [ 24%]
tests/pytest/test_feature_creation/test_feat_generation/test_feat_node.py .                       [ 25%]
tests/pytest/test_feature_creation/test_feat_generation/test_inv_node.py .                        [ 26%]
tests/pytest/test_feature_creation/test_feat_generation/test_log_node.py .                        [ 27%]
tests/pytest/test_feature_creation/test_feat_generation/test_model_node.py .                      [ 28%]
tests/pytest/test_feature_creation/test_feat_generation/test_mult_node.py .                       [ 30%]
tests/pytest/test_feature_creation/test_feat_generation/test_neg_exp_node.py .                    [ 31%]
tests/pytest/test_feature_creation/test_feat_generation/test_sin_node.py .                        [ 32%]
tests/pytest/test_feature_creation/test_feat_generation/test_six_pow_node.py .                    [ 33%]
tests/pytest/test_feature_creation/test_feat_generation/test_sq_node.py .                         [ 34%]
tests/pytest/test_feature_creation/test_feat_generation/test_sqrt_node.py .                       [ 36%]
tests/pytest/test_feature_creation/test_feat_generation/test_sub_node.py .                        [ 37%]
tests/pytest/test_feature_creation/test_feature_space/test_feature_space.py .                     [ 38%]
tests/pytest/test_feature_creation/test_feature_space/test_gen_feature_space_from_file.py .       [ 39%]
tests/pytest/test_feature_creation/test_feature_space/test_gen_feature_space_selected_from_file.py . [ 40%]
                                                                                                  [ 40%]
tests/pytest/test_feature_creation/test_feature_space/test_units.py .                             [ 42%]
tests/pytest/test_feature_creation/test_parameterize/test_lorentizan.py .                         [ 43%]
tests/pytest/test_feature_creation/test_parameterize/test_param_abs.py .                          [ 44%]
tests/pytest/test_feature_creation/test_parameterize/test_param_abs_diff.py .                     [ 45%]
tests/pytest/test_feature_creation/test_parameterize/test_param_add.py .                          [ 46%]
tests/pytest/test_feature_creation/test_parameterize/test_param_cb.py .                           [ 48%]
tests/pytest/test_feature_creation/test_parameterize/test_param_cbrt.py .                         [ 49%]
tests/pytest/test_feature_creation/test_parameterize/test_param_cos.py .                          [ 50%]
tests/pytest/test_feature_creation/test_parameterize/test_param_div.py .                          [ 51%]
tests/pytest/test_feature_creation/test_parameterize/test_param_exp.py .                          [ 53%]
tests/pytest/test_feature_creation/test_parameterize/test_param_inv.py .                          [ 54%]
tests/pytest/test_feature_creation/test_parameterize/test_param_log.py .                          [ 55%]
tests/pytest/test_feature_creation/test_parameterize/test_param_neg_exp.py .                      [ 56%]
tests/pytest/test_feature_creation/test_parameterize/test_param_sin.py .                          [ 57%]
tests/pytest/test_feature_creation/test_parameterize/test_param_six_pow.py .                      [ 59%]
tests/pytest/test_feature_creation/test_parameterize/test_param_sq.py .                           [ 60%]
tests/pytest/test_feature_creation/test_parameterize/test_param_sqrt.py .                         [ 61%]
tests/pytest/test_feature_creation/test_parameterize/test_param_sub.py .                          [ 62%]
tests/pytest/test_inputs/test_inputs.py .                                                         [ 63%]
tests/pytest/test_model_eval/test_model_node/test_binary_binary.py .                              [ 65%]
tests/pytest/test_model_eval/test_model_node/test_binary_unary.py .                               [ 66%]
tests/pytest/test_model_eval/test_model_node/test_op_model_nodes.py .                             [ 67%]
tests/pytest/test_model_eval/test_model_node/test_unary_binary.py .                               [ 68%]
tests/pytest/test_model_eval/test_model_node/test_unary_unary.py .                                [ 69%]
tests/pytest/test_model_eval/test_models/test_log_reg_model.py .                                  [ 71%]
tests/pytest/test_model_eval/test_models/test_reg_model.py .                                      [ 72%]
tests/pytest/test_model_eval/test_models/test_reg_model_param.py .                                [ 73%]
tests/pytest/test_model_eval/test_param_model_node/test_binary_binary_param.py .                  [ 74%]
tests/pytest/test_model_eval/test_param_model_node/test_binary_unary_param.py .                   [ 75%]
tests/pytest/test_model_eval/test_param_model_node/test_param_op_model_nodes.py .                 [ 77%]
tests/pytest/test_model_eval/test_param_model_node/test_unary_binary_param.py .                   [ 78%]
tests/pytest/test_model_eval/test_param_model_node/test_unary_unary_param.py .                    [ 79%]
tests/pytest/test_nl_opt/test_nl_opt_utils.py .                                                   [ 80%]
tests/pytest/test_sklearn/test_sisso_class_sklearn.py .....                                       [ 86%]
tests/pytest/test_sklearn/test_sisso_log_reg_sklearn.py .....                                     [ 92%]
tests/pytest/test_sklearn/test_sisso_reg_sklearn.py .....                                           [ 98%]
tests/pytest/test_sklearn/test_sklearn_interface.py .                                               [100%]

=========================================== 83 passed in 44.28s ===========================================

```
### Installation Settings for Basic Flags

Terms highlighted in red should always be set, and terms highlighted in blue should be set based on what you want to do.

#### <span style="color:red"> CMAKE_CXX_COMPILER </span>

The C++ compiler used to compile the code

#### CMAKE_C_COMPILER

The C compiler used to compile external dependencies (GTest and {fmt})

#### CMAKE_INSTALL_PREFIX

*Default the main SISSO++ directory*
The directory used to install the final binaries

#### <span style="color:red"> CMAKE_CXX_FLAGS </span>

Define the flags used by the compiler to build the C++ source files. It is recommended to use `-O2` or `-O3` level of optimizations, but it can be changed to match the compiler.

#### CMAKE_C_FLAGS

Define the flags used by the compiler to build the C source files (GTest and {fmt}). It is recommended to use `-O2` or `-O3` level of optimizations, but it can be changed to match the compiler.

#### <span style="color:blue"> CMAKE_INSTALL_PREFIX </span>

*Default: SISSO++ Source Directory*

Path to where the final library and executable files will be placed

### Installation Settings

#### <span style="color:blue"> BUILD_PARAMS </span>

*Default: OFF*

If `BUILD_PARAMS` is ON then build the operators with non-linearly optimized scale and shift parameters, and the relevant optimization files. With this flag on `NLopt`, [the parameterized operators](../cpp_api/param_node), and [the non-linear optimization](../cpp_api/nl_opt) classes will be built.

#### <span style="color:blue"> BUILD_PYTHON </span>

*Default: ON*

If `BUILD_PYTHON` is ON then build the python bindings

#### PYTHON_INSTDIR

*Default: The default python installation directory associated with the binary*

The base directory to install the python bindings to. If this is not the default value, then you must also add the directory to your `PYTHONPATH` environment variable.

#### BUILD_TESTS

*Default: OFF*

If `BULD_TESTS` is ON then build GTest based tests

#### EXTERNAL_BOOST

*Default: OFF*

If `EXTERNAL_BOOST` is ON then use the pre-built Boost Libraries currently in your path or in `$ENV{BOOST_ROOT}`
Here the `-O3` flag is for optimizations, it is recommended to stay as `-O3` or `-O2`, but it can be changed to match compiler requirements.

We recommend to use `-DEXTERNAL_BOOST=OFF` if you are building the python bindings with an anaconda environment or if you are using non-standard C++ or MPI libraries. If the boost libraries have already been built with the same C++ compiler, MPI library, and Python environment, then `EXTERNAL_BOOST=ON` can be used.

#### EXTERNAL_BUILD_N_PROCS

*Default: 1*

The number of processes passed to make when building external libraries.
---
title: 'SISSO++: A C++ Implementation of the Sure-Independence Screening and Sparsifying Operator Approach'
tags:
  - SISSO
  - Symbolic Regression
  - Physics
  - C++
  - Python
authors:
  - name: Thomas A. R. Purcell
    orcid: 0000-0003-4564-7206
    affiliation: 1
  - name: Matthias Scheffler
    affiliation: 1
  - name: Christian Carbogno
    orcid: 0000-0003-0635-8364
    affiliation: 1
  - name: Luca M. Ghiringhelli
    orcid: 0000-0001-5099-3029
    affiliation: 1
affiliations:
 - name: NOMAD Laboratory at the Fritz Haber Institute of the Max Planck Society and Humboldt University, Berlin, Germany
   index: 1
date: September 2021
bibliography: paper.bib
---

# Summary
The sure independence screening and sparsifying operator (SISSO) approach [@Ouyang2017] is an algorithm belonging to the field of artificial intelligence and more specifically a combination of symbolic regression and compressed sensing.
As a symbolic regression method, SISSO is used to identify mathematical functions, i.e. the descriptors, that best predict the target property of a data set.
Furthermore, the compressed sensing aspect of SISSO, allows it to find sparse linear models using tens to thousands of data points.
SISSO is introduced for both regression and classification tasks.
In practice, SISSO first constructs a large and exhaustive feature space of trillions of potential descriptors by taking in a set of user-provided *primary features* as a dataframe, and then iteratively applying a set of unary and binary operators, e.g. addition, multiplication, exponentiation, and squaring, according to a user-defined specification.
From this exhaustive pool of candidate descriptors, the ones most correlated to a target property are identified via sure-independence screening, from which the low-dimensional linear models with the lowest error are found via an $\ell_0$ regularization.

Because symbolic regression generates an interpretable equation, it has become an increasingly popular concept across scientific disciplines [@Wang2019a; @Neumann2020; @Udrescu2020a].
A particular advantage of these approaches are their capability to model complex phenomena using relatively simple descriptors.
SISSO has been used successfully in the past to model, explore, and predict important material properties, including the stability of different phases [@Bartel2018a; @Schleder2020]; the catalytic activity and reactivity [@Han2021; @Xu2020; @Andersen2019; @Andersen2021]; and glass transition temperatures [@Pilania2019].
Beyond regression problems, SISSO has also been used successfully to classify materials into different crystal prototypes [@Ouyang2019a], or whether a material crystallizes in its ground state as a perovskite [@Bartel2019a], or to determine whether a material is a topological insulator or not [@Cao2020].

The SISSO++ package is an open-source (Apache-2.0 licence), modular, and extensible C++ implementation of the SISSO method with Python bindings.
Specifically, SISSO++ applies this methodology for regression, log regression, and classification problems.
Additionally, the library includes multiple Python functions to facilitate the post-processing, analyzing, and visualizing of the resulting models.

# Statement of need
The main goal of the SISSO++ package is to provide a user-friendly, easily extendable version of the SISSO method for the scientific community.
While both a FORTRAN [@Ouyang] and a Matlab [@MATLAB_SISSO] implementation of SISSO exist, their lack of native Python interfaces led to the development of multiple separate Python wrappers [@SISSOKit; @pySISSO].
This package looks to rectify this situation by providing an implementation that has native Python bindings and can be used both in a massively parallel environment for discovering the descriptors and on personal computing devices for analyzing and visualizing the results.
For this reason, all computationally intensive task are written in C++ and support parallelization via MPI and OpenMP.
Additionally, the Python bindings allow one to easily incorporate the methods into computational workflows and postprocess results.
Furthermore, this enables the integration of SISSO into existing machine-learning frameworks, e.g. scikit-learn [@scikit-learn], via submodules.
The code is designed in a modular fashion, which simplifies the process of extending the code for other applications.
Finally the project's extensive documentation and tutorials provide a good access point for new users of the method.

# Features
The following features are implemented in SISSO++:

  - A C++ library for using SISSO to find analytical models for a given problem

  - Python bindings to be able to interface with the C++ objects in a Python environment

  - Postprocessing tools for visualizing models and analyzing results using Matplotlib [@Hunter:2007]

  - Access to solve an *n*-dimensional classification model using a combination of calculating the convex-hull overlap and a linear-support vector machine solver

  - The ability to include non-linear parameters within features (e.g. exp($\alpha$x) and ln(x + $\beta$)) [@Purcell2021]

  - A `scikit-learn` interface

  - Complete API documentation defining all functions of the code

  - Tutorials and Quick-Start Guides describing the basic functionality of the code

# Code Dependencies
The following libraries are used by SISSO++:

  - Boost Serialization, MPI, System, and Filesystem are used for MPI communication and file management

  - NLopt [@NLOpt] is used to optimize the non-linear bias and scale parameters within features

  - The CLP library from Coin-OR [@Coin-CLP] is used to find the number of points in the convex hull overlap region for classification problems

  - LIBSVM [@CC01a] is used to find the linear-SVM model for classification problems

  - pybind11 [@pybind11] is used to create the python bindings

  - Scikit-learn [@scikit-learn], [@sklearn_api] is used to update the SVM model for classification problems within Python

  - NumPy [@harris2020array] and pandas [@reback2020pandas; @mckinney-proc-scipy-2010] are used to represent the data structures within Python and perform array operations.

  - seaborn [@Waskom2021] and Matplotlib [@Hunter:2007] are used to generate the plots during postprocessing

# Acknowledgements
The authors would like to thank Markus Rampp and Meisam Tabriz of the MPCDF for technical support. We would also like to thank Lucas Foppa, Jingkai Quan, Aakash Naik, James Dean, Timur Bazhirov, and Luigi Sbailò for testing and providing valuable feedback. T.P. would like to thank the Alexander von Humboldt Foundation for their support through the Alexander von Humboldt Postdoctoral Fellowship Program. This project was supported by TEC1p (the European Research Council (ERC) Horizon 2020 research and innovation programme, grant agreement No. 740233), BiGmax (the Max Planck Society’s Research Network on Big-Data-Driven Materials-Science), and the FAIRmat consortium of the NFDI and FAIR-DI e.V. association.

# References
.. SISSO++ documentation master file, created by
   sphinx-quickstart on Fri Aug  6 17:42:26 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SISSO++
=======
.. toctree::
   :maxdepth: 2
   :caption: Contents:

This package provides a C++ implementation of SISSO with built in Python bindings for an efficient python interface.
Future work will expand the python interface to include more postporcessing analysis tools.

Table of Contents
^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   quick_start/QuickStart
   tutorial/tutorial
   cpp_api/cpp_api
   python_api/python_api
   development/develop
.. _tut:

Tutorial
========
.. toctree::
    :maxdepth: 2

    0_intro
    1_command_line
    2_python
    2b_sklearn
    3_classification
.. _api_mpi:

MPI
===
.. toctree::
    :maxdepth: 2

MPI_Interface
-------------
.. doxygenfile:: MPI_Interface.hpp
   :project: SISSO++

MPI_ops
-------
.. doxygenfile:: MPI_Ops.hpp
   :project: SISSO++
.. _api_models:


Model
-----
.. doxygenfile:: Model.hpp
   :project: SISSO++

ModelClassifier
---------------
.. doxygenfile:: ModelClassifier.hpp
   :project: SISSO++

ModelRegressor
--------------
.. doxygenfile:: ModelRegressor.hpp
   :project: SISSO++

ModelLogRegressor
-----------------
.. doxygenfile:: ModelLogRegressor.hpp
   :project: SISSO++

.. _api_solvers:

SISSOSolver
-----------
.. doxygenfile:: SISSOSolver.hpp
   :project: SISSO++

SISSOClassifier
---------------
.. doxygenfile:: SISSOClassifier.hpp
   :project: SISSO++

SISSORegressor
--------------
.. doxygenfile:: SISSORegressor.hpp
   :project: SISSO++

SISSOLogRegressor
-----------------
.. doxygenfile:: SISSOLogRegressor.hpp
   :project: SISSO++

.. _api_nl_opt:

Non-Linear Parameter Optimization
=================================
.. toctree::
    :maxdepth: 2

NLOptimizer
-----------
.. doxygenfile:: NLOptimizer.hpp
   :project: SISSO++

NLOptimizerRegression
---------------------
.. doxygenfile:: NLOptimizerRegression.hpp
   :project: SISSO++

NLOptimizerLogRegression
------------------------
.. doxygenfile:: NLOptimizerLogRegression.hpp
   :project: SISSO++

NLOptimizerClassification
-------------------------
.. doxygenfile:: NLOptimizerClassification.hpp
   :project: SISSO++
.. _api_sis_utils:

Sure-Independence Screening Utilities
=====================================
.. toctree::
    :maxdepth: 2

Compare Features
----------------
.. doxygenfile:: compare_features.hpp
   :project: SISSO++


Projection Functions
--------------------
.. doxygenfile:: project.hpp
   :project: SISSO++
.. _api_node:

Nodes
=====
.. toctree::
    :maxdepth: 2

Node
----
.. doxygenfile:: Node.hpp
   :project: SISSO++

FeatureNode
-----------
.. doxygenfile:: FeatureNode.hpp
   :project: SISSO++

ModelNode
---------
.. doxygenfile:: ModelNode.hpp
   :project: SISSO++


OperatorNode
------------
.. doxygenfile:: OperatorNode.hpp
   :project: SISSO++

AddNode
-------
.. doxygenfile:: add.hpp
   :project: SISSO++

SubNode
-------
.. doxygenfile:: subtract.hpp
   :project: SISSO++

AbsDiffNode
-----------
.. doxygenfile:: absolute_difference.hpp
   :project: SISSO++

MultNode
--------
.. doxygenfile:: multiply.hpp
   :project: SISSO++

DivNode
-------
.. doxygenfile:: divide.hpp
   :project: SISSO++

InvNode
-------
.. doxygenfile:: inverse.hpp
   :project: SISSO++

SqNode
------
.. doxygenfile:: square.hpp
   :project: SISSO++

CbNode
------
.. doxygenfile:: cube.hpp
   :project: SISSO++

SixPowNode
----------
.. doxygenfile:: sixth_power.hpp
   :project: SISSO++

SqrtNode
--------
.. doxygenfile:: square_root.hpp
   :project: SISSO++

CbrtNode
--------
.. doxygenfile:: cube_root.hpp
   :project: SISSO++

ExpNode
-------
.. doxygenfile:: exponential.hpp
   :project: SISSO++

NegExpNode
----------
.. doxygenfile:: negative_exponential.hpp
   :project: SISSO++

LogNode
-------
.. doxygenfile:: log.hpp
   :project: SISSO++

SinNode
-------
.. doxygenfile:: sin.hpp
   :project: SISSO++

CosNode
-------
.. doxygenfile:: cos.hpp
   :project: SISSO++

AbsNode
-------
.. doxygenfile:: absolute_value.hpp
   :project: SISSO++
.. _api_fc:

Feature Creation
================
.. toctree::
    :maxdepth: 2

Features
--------
.. toctree::
    :maxdepth: 3

    node
    node_utils

Non-Linearly Optimized Features
-------------------------------
.. toctree::
    :maxdepth: 3

    param_node

Feature Space
-------------
.. toctree::
    :maxdepth: 2

    FeatureSpace
    sis_utils
.. _api_ut:

Utilities
=========
.. toctree::
    :maxdepth: 2

    classification
    Inputs
    loss_function
    math_utils
    mpi_interface
    nl_opt
    utils
.. _api_utils:

General Utilities
=================
.. toctree::
    :maxdepth: 2

enum
----
.. doxygenfile:: enum.hpp
   :project: SISSO++

String Utilities
----------------
.. doxygenfile:: string_utils.hpp
   :project: SISSO++

Vector Utilities
----------------
.. doxygenfile:: vector_utils.hpp
   :project: SISSO++

.. _api_node_utilities:

Node Utilities
==============
.. toctree::
    :maxdepth: 2

node_utils
----------
.. doxygenfile:: node/utils.hpp
   :project: SISSO++

Node Python Utils
-----------------
.. doxygenfile:: python/py_binding_cpp_def/feature_creation/node_utils.hpp
   :project: SISSO++

Node Data Storage
-----------------
.. doxygenfile:: nodes_value_containers.hpp
   :project: SISSO++

Units
-----
.. doxygenfile:: Unit.hpp
   :project: SISSO++
.. _api_math_utils:

Math Utilities
==============
.. toctree::
    :maxdepth: 2

.. doxygenfile:: math_funcs.hpp
   :project: SISSO++

.. _api:

C++ API
=======
.. toctree::
    :maxdepth: 2

    Inputs_sec
    FeatureCreation
    DescriptorIdentification
    Utilities

ParamNodes
==========
.. toctree::
    :maxdepth: 2

AddParamNode
------------
.. doxygenfile:: parameterized_add.hpp
   :project: SISSO++

SubParamNode
------------
.. doxygenfile:: parameterized_subtract.hpp
   :project: SISSO++

AbsDiffParamNode
----------------
.. doxygenfile:: parameterized_absolute_difference.hpp
   :project: SISSO++

MultParamNode
-------------
.. doxygenfile:: parameterized_multiply.hpp
   :project: SISSO++

DivParamNode
------------
.. doxygenfile:: parameterized_divide.hpp
   :project: SISSO++

InvParamNode
------------
.. doxygenfile:: parameterized_inverse.hpp
   :project: SISSO++

SqParamNode
-----------
.. doxygenfile:: parameterized_square.hpp
   :project: SISSO++

CbParamNode
-----------
.. doxygenfile:: parameterized_cube.hpp
   :project: SISSO++

SixPowParamNode
---------------
.. doxygenfile:: parameterized_sixth_power.hpp
   :project: SISSO++

SqrtParamNode
-------------
.. doxygenfile:: parameterized_square_root.hpp
   :project: SISSO++

CbrtParamNode
-------------
.. doxygenfile:: parameterized_cube_root.hpp
   :project: SISSO++

ExpParamNode
------------
.. doxygenfile:: parameterized_exponential.hpp
   :project: SISSO++

NegExpParamNode
---------------
.. doxygenfile:: parameterized_negative_exponential.hpp
   :project: SISSO++

LogParamNode
------------
.. doxygenfile:: parameterized_log.hpp
   :project: SISSO++

SinParamNode
------------
.. doxygenfile:: parameterized_sin.hpp
   :project: SISSO++

CosParamNode
------------
.. doxygenfile:: parameterized_cos.hpp
   :project: SISSO++

AbsParamNode
------------
.. doxygenfile:: parameterized_absolute_value.hpp
   :project: SISSO++

.. _api_feat_space:

FeatureSpace
------------
.. doxygenfile:: FeatureSpace.hpp
   :project: SISSO++
.. _api_inputs:

Inputs
======
.. toctree::
    :maxdepth: 2

    Inputs
.. _api_loss:

Loss Functions
==============
.. toctree::
    :maxdepth: 2

LossFunction
------------
.. doxygenfile:: LossFunction.hpp
   :project: SISSO++

LossFunctionConvexHull
----------------------
.. doxygenfile:: LossFunctionConvexHull.hpp
   :project: SISSO++

LossFunctionPearsonRMSE
-----------------------
.. doxygenfile:: LossFunctionPearsonRMSE.hpp
   :project: SISSO++

LossFunctionLogPearsonRMSE
--------------------------
.. doxygenfile:: LossFunctionLogPearsonRMSE.hpp
   :project: SISSO++

Loss Function Utilities
-----------------------
.. doxygenfile:: loss_function/utils.hpp
   :project: SISSO++
.. _api_inputs:

InputParser
-----------
.. doxygenfile:: InputParser.hpp
   :project: SISSO++
.. _api_di:

Descriptor Identification
=========================
.. toctree::
    :maxdepth: 2

Solvers
-------
.. toctree::
    :maxdepth: 2

    solvers

Models
------
.. toctree::
    :maxdepth: 2

    model
.. _api_classification:

Classification
==============
.. toctree::
    :maxdepth: 2

ConvexHull1D
------------
.. doxygenfile:: ConvexHull1D.hpp
   :project: SISSO++

LPWrapper
---------
.. doxygenfile:: LPWrapper.hpp
   :project: SISSO++

prop_sorted_d_mat
-----------------
.. doxygenfile:: prop_sorted_d_mat.hpp
   :project: SISSO++

SVMWrapper
----------
.. doxygenfile:: SVMWrapper.hpp
   :project: SISSO++

.. _devel:

Development
===========
.. toctree::
    :maxdepth: 2

    CONTRIBUTING
    license
    Credits
    References
.. _quick_start:

Quick-Start Guide
=================
.. toctree::
    :maxdepth: 2

    Installation
    code_ref
.. _py_api_fc:

Feature Creation
================
.. toctree::
    :maxdepth: 2

Features
--------
.. toctree::
    :maxdepth: 3

    py_node
    py_node_utils

Feature Space
-------------
.. toctree::
    :maxdepth: 2

    py_FeatureSpace
.. _py_api_post:

.. currentmodule:: sissopp.postprocess

Postprocessing
==============
.. toctree::
    :maxdepth: 2

    py_post_class
    py_post_load_model
    py_post_cv_check
    py_post_parity
    py_post_map
    py_post_plot_utils
    py_post_feat_space_prev
.. _py_api_di:

Descriptor Identification
=========================
.. toctree::
    :maxdepth: 2

Solvers
-------
.. toctree::
    :maxdepth: 2

    py_solvers

Models
------
.. toctree::
    :maxdepth: 2

    py_model
.. _py_api_solvers:

SISSOSolver
-----------

.. autoclass:: sissopp.SISSOSolver
    :special-members: __init__
    :members:
    :undoc-members:


SISSOClassifier
---------------

.. autoclass:: sissopp.SISSOClassifier
    :special-members: __init__
    :members:
    :undoc-members:


SISSORegressor
--------------

.. autoclass:: sissopp.SISSORegressor
    :special-members: __init__
    :members:
    :undoc-members:


SISSOLogRegressor
-----------------

.. autoclass:: sissopp.SISSOLogRegressor
    :special-members: __init__
    :members:
    :undoc-members:


.. _py_api_cv_conv


Cross-Validation Convergence
----------------------------

.. currentmodule:: sissopp.postprocess.check_cv_convergence

.. autofunction:: jackknife_cv_conv_est

.. currentmodule:: sissopp.postprocess.plot.cv_error_plot

.. autofunction:: plot_validation_rmse

.. autofunction:: plot_errors_dists

.. autofunction:: add_violin

.. autofunction:: add_boxplot

.. autofunction:: add_pointplot
.. _py_api_load_models

Load Models
-----------

.. currentmodule:: sissopp.postprocess.load_models

.. autofunction:: sort_model_file_key
.. autofunction:: load_model
.. autofunction:: get_models
.. autofunction:: create_model_csv
.. _py_api_node:

Node
----
.. autoclass:: sissopp.Node
    :special-members: __init__
    :members:
    :undoc-members:


FeatureNode
-----------
.. autoclass:: sissopp.FeatureNode
    :special-members: __init__
    :members:
    :undoc-members:

ModelNode
---------
.. autoclass:: sissopp.ModelNode
    :special-members: __init__
    :members:
    :undoc-members:

OperatorNode1
-------------
.. autoclass:: sissopp.OperatorNode1
    :special-members: __init__
    :members:
    :undoc-members:

AddNode
-------
.. autoclass:: sissopp.AddNode
    :special-members: __init__
    :members:
    :undoc-members:

SubNode
-------
.. autoclass:: sissopp.SubNode
    :special-members: __init__
    :members:
    :undoc-members:

AbsDiffNode
-----------
.. autoclass:: sissopp.AbsDiffNode
    :special-members: __init__
    :members:
    :undoc-members:

MultNode
--------
.. autoclass:: sissopp.MultNode
    :special-members: __init__
    :members:
    :undoc-members:

DivNode
-------
.. autoclass:: sissopp.DivNode
    :special-members: __init__
    :members:
    :undoc-members:

InvNode
-------
.. autoclass:: sissopp.InvNode
    :special-members: __init__
    :members:
    :undoc-members:

SqNode
------
.. autoclass:: sissopp.SqNode
    :special-members: __init__
    :members:
    :undoc-members:

CbNode
------
.. autoclass:: sissopp.CbNode
    :special-members: __init__
    :members:
    :undoc-members:

SixPowNode
----------
.. autoclass:: sissopp.SixPowNode
    :special-members: __init__
    :members:
    :undoc-members:

SqrtNode
--------
.. autoclass:: sissopp.SqrtNode
    :special-members: __init__
    :members:
    :undoc-members:

CbrtNode
--------
.. autoclass:: sissopp.CbrtNode
    :special-members: __init__
    :members:
    :undoc-members:

ExpNode
-------
.. autoclass:: sissopp.ExpNode
    :special-members: __init__
    :members:
    :undoc-members:

NegExpNode
----------
.. autoclass:: sissopp.NegExpNode
    :special-members: __init__
    :members:
    :undoc-members:

LogNode
-------
.. autoclass:: sissopp.LogNode
    :special-members: __init__
    :members:
    :undoc-members:

SinNode
-------
.. autoclass:: sissopp.SinNode
    :special-members: __init__
    :members:
    :undoc-members:

CosNode
-------
.. autoclass:: sissopp.CosNode
    :special-members: __init__
    :members:
    :undoc-members:

AbsNode
-------
.. autoclass:: sissopp.AbsNode
    :special-members: __init__
    :members:
    :undoc-members:

AddParamNode
------------
.. autoclass:: sissopp.AddParamNode
    :special-members: __init__
    :members:
    :undoc-members:

SubParamNode
------------
.. autoclass:: sissopp.SubParamNode
    :special-members: __init__
    :members:
    :undoc-members:

AbsDiffParamNode
----------------
.. autoclass:: sissopp.AbsDiffParamNode
    :special-members: __init__
    :members:
    :undoc-members:

MultParamNode
-------------
.. autoclass:: sissopp.MultParamNode
    :special-members: __init__
    :members:
    :undoc-members:

DivParamNode
------------
.. autoclass:: sissopp.DivParamNode
    :special-members: __init__
    :members:
    :undoc-members:

InvParamNode
------------
.. autoclass:: sissopp.InvParamNode
    :special-members: __init__
    :members:
    :undoc-members:

SqParamNode
-----------
.. autoclass:: sissopp.SqParamNode
    :special-members: __init__
    :members:
    :undoc-members:

CbParamNode
-----------
.. autoclass:: sissopp.CbParamNode
    :special-members: __init__
    :members:
    :undoc-members:

SixPowParamNode
---------------
.. autoclass:: sissopp.SixPowParamNode
    :special-members: __init__
    :members:
    :undoc-members:

SqrtParamNode
-------------
.. autoclass:: sissopp.SqrtParamNode
    :special-members: __init__
    :members:
    :undoc-members:

CbrtParamNode
-------------
.. autoclass:: sissopp.CbrtParamNode
    :special-members: __init__
    :members:
    :undoc-members:

ExpParamNode
------------
.. autoclass:: sissopp.ExpParamNode
    :special-members: __init__
    :members:
    :undoc-members:

NegExpParamNode
---------------
.. autoclass:: sissopp.NegExpParamNode
    :special-members: __init__
    :members:
    :undoc-members:

LogParamNode
------------
.. autoclass:: sissopp.LogParamNode
    :special-members: __init__
    :members:
    :undoc-members:

SinParamNode
------------
.. autoclass:: sissopp.SinParamNode
    :special-members: __init__
    :members:
    :undoc-members:

CosParamNode
------------
.. autoclass:: sissopp.CosParamNode
    :special-members: __init__
    :members:
    :undoc-members:

AbsParamNode
------------
.. autoclass:: sissopp.AbsParamNode
    :special-members: __init__
    :members:
    :undoc-members:

.. _py_api:

Python API
==========
.. toctree::
    :maxdepth: 2

    py_Interface
    py_Inputs_sec
    py_FeatureCreation
    py_DescriptorIdentification
    py_Postprocessing
.. _py_api_inputs_sec:

Inputs
======
.. toctree::
    :maxdepth: 2

Inputs
------
.. toctree::
    :maxdepth: 3

    py_Inputs
.. _py_api_models:

Model
-----

.. autoclass:: sissopp.Model
    :special-members: __init__
    :members:
    :undoc-members:

ModelClassifier
---------------

.. currentmodule:: sissopp._sissopp

.. autoclass:: sissopp.ModelClassifier
    :special-members: __init__
    :members:
    :undoc-members:


ModelRegressor
--------------

.. currentmodule:: sissopp._sissopp

.. autoclass:: sissopp.ModelRegressor
    :special-members: __init__
    :members:
    :undoc-members:


ModelLogRegressor
-----------------

.. currentmodule:: sissopp._sissopp

.. autoclass:: sissopp.ModelLogRegressor
    :special-members: __init__
    :members:
    :undoc-members:


.. _py_api_read_data

Reading Datafiles
-----------------

.. currentmodule:: sissopp.py_interface.import_dataframe

.. autofunction:: read_csv
.. autofunction:: create_inputs
.. autofunction:: extract_col
.. autofunction:: get_unit
.. autofunction:: strip_units
.. _py_api_get_solvers:


Constructing Solvers
--------------------

.. currentmodule:: sissopp.py_interface.get_solver

.. autofunction:: get_fs_solver
.. _py_api_post_parity

Parity Plots
------------

.. currentmodule:: sissopp.postprocess.plot.parity_plot

.. autofunction:: plot_model_parity_plot
.. _py_api_est_feat_space:


Estimate Feature Space Size
---------------------------

.. currentmodule:: sissopp.py_interface.estimate_feat_space_size

.. autofunction:: get_max_number_feats
.. autofunction:: get_estimate_n_feat_next_rung
.. _py_api_post_an_fs


Feature Space Analysis
----------------------

.. currentmodule:: sissopp.postprocess.feature_space_analysis

.. autofunction:: get_prevelance_of_primary_features
.. _py_api_feat_space:

FeatureSpace
------------
.. autoclass:: sissopp.FeatureSpace
    :special-members: __init__
    :members:
    :undoc-members:
.. _py_api_inputs:

Inputs
------
.. autoclass:: sissopp.Inputs
    :special-members: __init__
    :members:
    :undoc-members:
.. _py_api_post_plot_utils

Plotting Utilities
------------------

.. currentmodule:: sissopp.postprocess.plot.utils

.. autofunction:: setup_plot_ax
.. autofunction:: adjust_box_widths
.. autofunction:: latexify
.. _py_api_class

Classification
--------------

.. currentmodule:: sissopp.postprocess.classification

.. autofunction:: update_model_svm

.. currentmodule:: sissopp.postprocess.plot.classification

.. autofunction:: plot_classification
.. _py_api_interface:

Python Interface
================
.. toctree::
    :maxdepth: 2

    py_reading_data_files
    py_constructing_solvers
    py_estimate_feature_space_size
.. _py_api_post_map


Plotting 2D Maps
----------------

.. currentmodule:: sissopp.postprocess.plot.maps

.. autofunction:: plot_2d_map
.. _py_api_node_utilities:

Node Utilities
==============
.. toctree::
    :maxdepth: 2

Node Python Utils
-----------------
.. autofunction:: sissopp.phi_selected_from_file

Node Data Storage
-----------------
.. autofunction:: sissopp.initialize_values_arr
.. autofunction:: sissopp.initialize_d_matrix_arr

Units
-----
.. autoclass:: sissopp.Unit
    :special-members: __init__
    :members:
    :undoc-members:
