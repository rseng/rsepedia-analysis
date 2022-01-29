<p align="center">
    <img
      src="https://github.com/pyomeca/biorbd_design/blob/main/logo_png/biorbd_full.png"
      alt="logo"
    />
</p>

BIORBD is a library to analyze biomechanical data. It provides several useful functions for the direct and inverse flow including rigid body (based on *Feathestone* equations implemented in RBDL) and muscle elements.

Biomechanical data are often analyzed using similar flow, that is inverse or direct. BIORBD implements these common analyses providing high-level and easy to use Python and MATLAB interfaces of an efficient C++ implementation. 

So, without further ado, let's begin our investigation of BIORBD!

You can get the online version of the paper for BIORBD here: [![DOI](https://joss.theoj.org/papers/10.21105/joss.02562/status.svg)](https://doi.org/10.21105/joss.02562)

# How to install
There are two main ways to install BIORBD on your computer: installing the binaries from Anaconda (easiest, but limited to C++ and Python3) or compiling the source code yourself (more versatile and up to date; for C++, Python3 and MATLAB).

## Anaconda (For Windows, Linux and Mac)
The easiest way to install BIORBD is to download the binaries from Anaconda (https://anaconda.org/) repositories (binaries are not available though for MATLAB). The project is hosted on the conda-forge channel (https://anaconda.org/conda-forge/biorbd).

After having installed properly an anaconda client [my suggestion would be Miniconda (https://conda.io/miniconda.html)] and loaded the desired environment to install BIORBD in, just type the following command:
```bash
conda install -c conda-forge biorbd
```
The binaries and includes of the core of BIORBD will be installed in `bin` and `include` folders of the environment respectively. Moreover, the Python3 binder will also be installed in the environment.

Please note that because of the way `Ipopt` is compiled on conda-forge, it was not possible to link it with `biorbd`. Therefore, the `MODULE_STATIC_OPTIM` was set to `OFF` for this particular OS.

The current building status for Anaconda release is as follow.

| License | Name | Downloads | Version | Platforms |
| --- | --- | --- | --- | --- |
|   <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/license-MIT-success" alt="License"/></a> | [![Conda Recipe](https://img.shields.io/badge/recipe-biorbd-green.svg)](https://anaconda.org/conda-forge/biorbd) | [![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/biorbd.svg)](https://anaconda.org/conda-forge/biorbd) | [![Conda Version](https://img.shields.io/conda/vn/conda-forge/biorbd.svg)](https://anaconda.org/conda-forge/biorbd) | [![Conda Platforms](https://img.shields.io/conda/pn/conda-forge/biorbd.svg)](https://anaconda.org/conda-forge/biorbd) |

## Compiling (For Windows, Linux and Mac)
The main drawback with downloading the pre-compiled version from Anaconda is that this version may be out-of-date (even if I do my best to keep the release versions up-to-date). Moreover, since it is already compiled, it doesn't allow you to modify BIORBD if you need to. Therefore, a more versatile way to enjoy BIORBD is to compile it by yourself.

The building status for the current BIORBD branches is as follow

| Name | Status |
| --- | --- |
| master | [![Build status](https://ci.appveyor.com/api/projects/status/fegqut6bypvix8ex/branch/master?svg=true)](https://ci.appveyor.com/project/pariterre/biorbd/branch/master) |
| Code coverage | [![codecov](https://codecov.io/gh/pyomeca/biorbd/branch/master/graph/badge.svg)](https://codecov.io/gh/pyomeca/biorbd) |
| DOI | [![DOI](https://zenodo.org/badge/124423173.svg)](https://zenodo.org/badge/latestdoi/124423173) |

### Dependencies
BIORBD relies on several libraries (namely eigen ([http://eigen.tuxfamily.org]) or CasADi ([https://web.casadi.org/]), rbdl-casadi (https://github.com/pyomeca/rbdl-casadi), tinyxml(http://www.grinninglizard.com/tinyxmldocs/index.html) and Ipopt (https://github.com/coin-or/Ipopt)) that one must install prior to compiling. Fortunately, all these dependencies are also hosted on the *conda-forge* channel of Anaconda. Therefore the following command will install everything you need to compile BIORBD:
```bash
conda install -c conda-forge rbdl [tinyxml] [ipopt] [pkgconfig] [cmake]
```
Please note:
- ```tinyxml``` is optional, but is required for reading VTP files;
- ```ipopt``` is optional, but is required for the *Static optimization* module;
- ```pkgconfig``` and ```cmake``` are very useful tools that can prevents lot of headaches when compiling; 


Additionnally, for the Python3 interface requires *numpy* (https://numpy.org/) and *SWIG* (http://www.swig.org/). Again, one can easily install these dependencies from Anaconda using the following command:
```bash
conda install -c conda-forge numpy swig
```

Finally, the MATLAB interface (indeed) requires MATLAB to be installed.

If ones is interested in developping BIORBD, the ```googletest``` suite is required to test your modifications. Fortunately, the CMake should download and compile the test suite for you!

### CMake
BIORBD comes with a CMake (https://cmake.org/) project. If you don't know how to use CMake, you will find many examples on Internet. The main variables to set are:

> `CMAKE_INSTALL_PREFIX` Which is the `path/to/install` BIORBD in. If you compile the Python3 binder, a valid installation of Python with Numpy should be installed relatived to this path.
>
> `BUILD_SHARED_LIBS` If you wan to build BIORBD in a shared `TRUE` or static `FALSE` library manner. Default is `TRUE`. Please note that due to the dependencies, on Windows BIORBD must be statically built.
>
> `CMAKE_BUILD_TYPE` Which type of build you want. Options are `Debug`, `RelWithDebInfo`, `MinSizeRel` or `Release`. This is relevant only for the build done using the `make` command. Please note that you will experience a slow BIORBD library if you compile it without any optimization (i.e. `Debug`), especially for all functions that requires linear algebra. 
>
> `MATH_LIBRARY_BACKEND` Choose between the two linear algebra backends, either `Eigen3` or `Casadi`. Default is `Eigen3`.
>
> `BUILD_EXAMPLE` If you want (`TRUE`) or not (`FALSE`) to build the C++ example. Default is `TRUE`.
>
> `BUILD_TESTS` If you want (`ON`) or not (`OFF`) to build the tests of the project. Please note that this will automatically download gtest (https://github.com/google/googletest). Default is `OFF`.
>
> `BUILD_DOC` If you want (`ON`) or not (`OFF`) to build the documentation of the project. Default is `OFF`.
>
> `BINDER_C` If you want (`ON`) or not (`OFF`) to build the low level C binder. Default is `OFF`. Please note that this binder is very light and will not contain most of BIORBD features.
>
> `BINDER_PYTHON3` If you want (`ON`) or not (`OFF`) to build the Python binder. Default is `OFF`.
>
> `SWIG_EXECUTABLE`  If `BINDER_PYTHON3` is set to `ON` then this variable should point to the SWIG executable. This variable should be found automatically.
>
> `BINDER_MATLAB` If you want (`ON`) or not (`OFF`) to build the MATLAB binder. Default is `OFF`. Pleaes note that `BINDER_MATLAB` can't be set to `ON` alonside to `CasADi` backend.
>
> `Matlab_ROOT_DIR` If `BINDER_MATLAB` is set to `ON` then this variable should point to the root path of MATLAB directory. Please note that the MATLAB binder is based on MATLAB R2018a API and won't compile on earlier versions. This variable should be found automatically, except on Mac where the value should manually be set to the MATLAB in the App folder. 
>
> `Matlab_biorbd_INSTALL_DIR` If `BINDER_MATLAB` is set to `ON` then this variable should point to the path where you want to install BIORBD. Typically, this is `{MY DOCUMENTS}/MATLAB`. The default value is the toolbox folder of MATLAB. Please note that if you leave the default value, you will probably need to grant administrator rights to the installer. In all cases, after the installation, you will have to add the path to the MATLAB search path by typing the following command in the MATLAB's prompt (or to add it to the `startup.m`) `addpath(genpath($Matlab_biorbd_INSTALL_DIR))`, and replacing `Matlab_biorbd_INSTALL_DIR` by your own path.
>
> `MODULE_ACTUATORS` If you want (`ON`) or not (`OFF`) to build with the actuators module. Default is `ON`. This allows to use exotic joint torques. 
>
> `MODULE_KALMAN` If you want (`ON`) or not (`OFF`) to build the Kalman filter module. Default is `ON`. The main reason to skip Kalman is that in `Debug` mode `Eigen3` will perform this very slowly and `CasADi` will always perform this slowly. 
>
> `MODULE_MUSCLES` If you want (`ON`) or not (`OFF`) to build with the muscle module. Default is `ON`. This allows to read and interact with models that include muscles.
>
> `MODULE_STATIC_OPTIM` If you want (`ON`) or not (`OFF`) to build the Static optimization module. Default is `ON` (if `ipopt` is found).
>
> `MODULE_VTP_FILES_READER` If you want (`ON`) or not (`OFF`) to build with the vtp files reader module. Default is `ON` (if `tinyxml` is found). This allows to read mesh files produced by `OpenSim`.
> 
> `SKIP_ASSERT` If you want (`ON`) or not (`OFF`) to skip the asserts in the functions (e.g. checks for sizes). Default is `OFF`. Putting this to `OFF` reduces the risks of Segmentation Faults, it will however slow down the code when using `Eigen3` backend.
>
> `SKIP_LONG_TESTS` If you want (`ON`) or not (`OFF`) to skip the tests that are long to perform. Default is `OFF`. This is useful when debugging. 


# How to use
BIORBD provides as much as possible explicit names for the filter so one can intuitively find what he wants from the library. Still, this is a C++ library and it can be sometimes hard to find what you need. Due to the varity of functions implemented in the library, minimal examples are shown here. One is encourage to have a look at the `example` and `test` folders to get a better overview of the possibility of the API. For an in-depth detail of the API, the Doxygen documentation (to come) is the way to go.

## The C++ API
The core code is written in C++, meaning that you can fully use BIORBD from C++.  Moreover, the linear algebra is using the Eigen library which makes it fairly easy to perform further computation and analyses.
The informations that follows is a basic guide that should allow you to perform everything you want to do.

### Create an empty yet valid model
To create a new valid yet empty model, just call the `biorbd::Model` class without parameter.
```C++
#include "biorbd.h"
int main()
{
    biorbd::Model myModel;
}
```
This model can thereafter be populated using the *biorbd* add methods. Even if this is not the prefered way of loading a model, one can have a look at the *src/ModelReader.cpp* in order to know what functions that must be called to populate the model manually. 

### Read and write a bioMod file
The prefered method to load a model is to read the in-house *.bioMod* format file.  To do so, one must simply call the `biorbd::Model` constructur with a valid path to the model. Afterward, one can modify manually the model and write it back to a new file. 
```C++
#include "biorbd.h"
int main()
{
    biorbd::Model myModel myModel("path/to/mymodel.bioMod");
    // Do some changes...
    biorbd::Writer::writeModel(myModel, "path/to/newFile.bioMod");
    return 0;
}
```
Please note that on Windows, the path must be `/` or `\\` separated (and not only`\`), for obvious reasons. 

### Perform some analyses
BIORBD is made to work with the RBDL functions (the doc can be found here https://rbdl.bitbucket.io/). Therefore, every functions available in RBDL is also available on BIORBD. Additionnal are of course also made available, for example the whole muscle module. 

The most obvious and probably the most used function is the forward kinematics, where one knows the configuration of the body and is interested in the resulting position of skin markers. The following code performs that task.
```C++
#include "biorbd.h"
int main()
{
    // Load the model
    biorbd::Model model("path/to/model.bioMod");
    
    // Prepare the model
    biorbd::rigidbody::GeneralizedCoordinates Q(model); 
    Q.setOnes()/10; // Set the model position
    
    // Perform forward kinematics
    std::vector<biorbd::rigidbody::NodeBone> markers(model.markers(Q));
    
    // Print the results
    for (auto marker : markers)
        std::cout << marker.name() << " is at the coordinates: " << marker.transpose() << std::endl;
    return 0;
}
```

Another common analysis to perform is to compute the effect of the muscles on the acceleration of the model. Assuming that the model that is loaded has muscles, the following code perform this task.
```C++
#include "biorbd.h"
int main()
{
    // Load the model
    biorbd::Model model("path/to/model.bioMod");
    
    // Prepare the model
    biorbd::rigidbody::GeneralizedCoordinates Q(model), Qdot(model); // position, velocity
    Q.setOnes()/10; // Set the model position
    Qdot.setOnes()/10; // Set the model velocity
    // Muscles activations
    std::vector<std::shared_ptr<biorbd::muscles::StateDynamics>> states(model.nbMuscleTotal());
    for (auto& state : states){
        state = std::make_shared<biorbd::muscles::StateDynamics>();
        state->setActivation(0.5); // Set the muscle activation
    }

    // Compute the joint torques based on muscle
    biorbd::rigidbody::GeneralizedTorque muscleTorque(
                model.muscularJointTorque(states, true, &Q, &Qdot));

    // Compute the acceleration of the model due to these torques
    biorbd::rigidbody::GeneralizedCoordinates Qddot(model);
    RigidBodyDynamics::ForwardDynamics(model, Q, Qdot, muscleTorque, Qddot);

    // Print the results
    std::cout << " The joints accelerations are: " << Qddot.transpose() << std::endl;
    return 0;
}
```

There are many other analyses and filters that are available. Please refer to the BIORBD and RBDL Docs to see what is available. 

## MATLAB
MATLAB (https://www.mathworks.com/) is a prototyping langage largely used in industry and fairly used by the biomechanical scientific community. Despite the existence of Octave as an open-source and very similar language or the growing popularity of Python as a free and open-source alternative, MATLAB remains an important player as a programming languages. Therefore BIORBD comes with a binder for MATLAB (that can theoretically used with Octave as well with some minor changes to the CMakeLists.txt file).

Most of the functions available in C++ are also available in MATLAB. Still, they were manually binded, therefore it may happen that some important one (for you) are not there. If so, do not hesitate to open an issue on GitHub to required the add of that particular function. The philosophy behind the MATLAB binder is that you open a particular model and a reference to that model is gave back to you. Thereafter, the functions can be called, assuming the pass back that model reference. That implies, however, that ones must himself deallocate the memory of the model when it is no more needed. Failing to do so results in an certain memory leak.

### Perform some analyses
Please find here the same tasks previously described for the C++ interface done in the MATLAB interface. Notice that the MATLAB interface takes advantage of the matrix nature of MATLAB and therefore can usually perform the analyses on multiple frames at once. 

Forward kinematics can be performed as follow
```MATLAB
nFrames = 10; % Set the number of frames to simulate

% Load the model
model = biorbd('new', 'path/to/model.bioMod');

% Prepare the model
Q = ones(biorbd('nQ', model), nFrames)/10; % Set the model position

% Perform the forward kinematics
markers = biorbd('markers', model, Q);

% Print the results
disp(markers);

% Deallocate the model
biorbd('delete', model);
```

The joint accelerations from muscle activations can be performed as follow
```MATLAB
nFrames = 10; % Set the number of frames to simulate

% Load the model
model = biorbd('new', 'path/to/model.bioMod');

% Prepare the model
Q = ones(biorbd('nQ', model), nFrames)/10; % Set the model position
Qdot = ones(biorbd('nQdot', model), nFrames)/10; % Set the model velocity
activations = ones(biorbd('nMuscles', model), nFrames)/2; % Set muscles activations

% Compute the joint torques based on muscle
jointTorque = biorbd('jointTorqueFromActivation', model, activations, Q, Qdot);

% Compute the acceleration of the model due to these torques
Qddot = biorbd('forwardDynamics', model, Q, Qdot, jointTorque);

% Print the results
disp(Qddot);

% Deallocate the model
biorbd('delete', model);
```

### Help
One can print all the available functions by type the `help` command
```MATLAB
biorbd('help')
```
Please note that it seems that on Windows, the command returns nothing. One must therefore look in the source code (`biorbd/binding/matlab/Matlab_help.h`) what should the command have returned.


## Python 3
Python (https://www.python.org/) is a scripting language that has taken more and more importance over the past years. So much that now it is one of the preferred language of the scientific community. Its simplicity yet its large power to perform a large variety of tasks makes it a certainty that its popularity won't decrease for the next years.

To interface the C++ code with Python, SWIG is a great tool. It creates very rapidly an interface in the target language with minimal code to write. However, the resulting code in the target language can be far from being easy to use. In effect, it gives a mixed-API not far from the original C++ language, which may not comply to best practices of the target language. When this is useful to rapidly create an interface, it sometime lacks of user-friendliness and expose the user to the possibility of the C++ such as segmentation fault (unlike the MATLAB API which won't suffer from this devil problem). 

BIORBD interfaces the C++ code using SWIG. While it has some inherent limit as discussed previously, it has the great advantage of providing almost for free the complete API. Because of that, much more of the C++ API is interfaced in Python than the MATLAB one. Again, if for some reason, part of the code which is not accessible yet is important for you, don't hesitate to open an issue asking for that particular feature!

### Perform some analyses
Please find here the same tasks previously described for the C++ interface done in the Python3 interface. Please note that the interface usually takes advantage of the numpy arrays in order to interact with the user while a vector is needed. 

Forward kinematics can be performed as follow
```Python
import numpy as np
import biorbd

# Load the model
model = biorbd.Model('path/to/model.bioMod')

# Prepare the model
Q = np.ones(model.nbQ())/10  # Set the model position

# Perform the forward kinematics
markers = model.markers(Q)

# Print the results
for marker in markers:
    print(marker.to_array())

```

The joint accelerations from muscle activations can be performed as follow
```Python
import numpy as np
import biorbd

# Load the model
model = biorbd.Model('path/to/model.bioMod')

# Prepare the model
Q = np.ones(model.nbQ())/10  # Set the model position
Qdot = np.ones(model.nbQ())/10  # Set the model velocity
states = model.stateSet()
for state in states:
    state.setActivation(0.5)  # Set muscles activations

# Compute the joint torques based on muscle
joint_torque = model.muscularJointTorque(states, Q, Qdot)

# Compute the acceleration of the model due to these torques
Qddot = model.ForwardDynamics(Q, Qdot, joint_torque)

# Print the results
print(Qddot.to_array())

```
# Model files
## *bioMod* files
The preferred method to load a model is by using a *.bioMod* file. This type of file is an in-house language that describes the segments of the model, their interactions and additionnal elements attached to them. The following section describe the structure of the file and all the tags that exists so far. 

Comments can be added to the file in a C-style way, meaning that everything a on line following a `//` will be considered as a comment and everything between `/*` and `*/` will also be ignored. 

Please note that the *bioMod* is not case dependent, so `Version`and `version` are for instance fully equivalent. The *bioMod* reader also ignore the tabulation, which is therefore only aesthetic. 

When a tag waits for multiple values, they must be separate by a space, a tabulation or a return of line. Also, anytime a tag waits for a value, it is possible to use simple equations (assuming no spaces are used) and/or variables. For example, the following snippet is a valid way to set the gravity parameter to $(0, 0, -9.81)$. 
```c
variables
    $my_useless_variable 0
endvariables
gravity 2*(1-1) -2*$my_useless_variable
        -9.81
```

### Header
#### version
The very first tag that **must** appear at the first line in file is the version of the file. The current version of the *.bioMod* files is $4$. Please note that most of the version are backward compatible, unless specified. This tag waits for $1$ value.
```c
version 4
```
From that point, the order of the tags is not important, header can even be at the end of the file. For simplicity though we suggest to put everything related to the header at the top of the file. 

#### gravity
The `gravity` tag is used to reorient and/or change the magnitude of the gravity. The default value is $(0, 0, -9.81)$. This tag waits for $3$ values.
```c
// Restate the default value
gravity 0 0 -9.81
```

#### variables / endvariables
The `variables / endvariables` tag pair allows to declare variables that can be used within the file. This allows for example to template the *bioMod* file by only changing the values in the variables. Please note that contrary to the rest of the file, the actual variables are case dependent. 

The `\$` sign is mandatory and later in the file, everything with that sign followed by the same name will be converted to the values specified in the tag. 
```c
// Restate the default value
variables
    $my_first_variable_is_an_int 10
    $my_second_variable_is_a_double 10.1
    $myThirdVariableIsCamelCase 1
    $myLastVariableIsPi pi
endvariables
```
As you may have noticed, the constant PI is defined as $3.141592653589793$.

### Definition of the model
A BIORBD model consists of a chain of segment, linked by joints with up to six DoF (3 translations, 3 rotations). It is imperative when attaching something to a segment of the model that particular segment must have been previously defined. For instance, if the `thorax` is attached to the `pelvis`, then the latter must be defined before the former in the file. 

#### Segment
The `segment xxx / endsegment`tag pair is the core of a *bioMod* file. It describes a segment of the model with the name `xxx`, that is most of the time a bone of the skeleton. For internal reasons, the name cannot be `root`. The `xxx` must be present and consists of $1$ string. The segment is composed of multiple subtags, described here. 

```c
segment default_segment
    parent ROOT
    rtinmatrix 0
    rt 0 0 0 xyz 0 0 0
    translation xyz
    rotations xyz
    com 0 0 0
    inertia 
        1 0 0
        0 1 0
        0 0 1
endsegment

segment second_segment
    parent default_segment
endsegment
```

##### parent
The `parent` tag is the name of the segment that particular segment is attached to. If no segment parent is provided, it is considered to be attached to the environment. The parent must be defined earlier in the file and is case dependent. This tag is waits for $1$ string.

##### rt
The `homogeneous matrix` of transformation (rototranslation) of the current segment relative to its parent. If `rtinmatrix` is defined to `1`, `rt` waits for a 4 by 4 matrix (the 3 x 3 matrix of rotation and the 4th column being the translation and the last row being 0 0 0 1); if it is defined to `0` it waits for the 3 rotations, the rotation sequence and the translation. The default value is the identity matrix.

##### rtinmatrix
The tag that defines if the `rt` is in matrix or not. If the `version` of the file is higher or equal than $3$, the default value is false ($0$), otherwise, it true ($1$). 

##### translations
The `translations` tag specifies the number of degrees-of-freedom in translation and their order. The possible values are `x`, `y` and/or `z` combined whatever fits the model. Please note that the vector of generalized coordinate will comply the the order wrote in this tag. If no translations are provided, then the segment has no translation relative to its parent. This tag waits for $1$ string.

##### rotations
The `rotations` tag specifies the number of degrees-of-freedom in rotation and their order. The possible values are `x`, `y` and/or `z` combined whatever fits the model. Please note that the vector of generalized coordinate will comply the the order wrote in this tag. If no rotations are provided, then the segment has no rotation relative to its parent. This tag waits for $1$ string.

##### mass
The `mass` tag specifies the mass of the segment in kilogram. This tag waits for $1$ value. The default value is $0.00001$. Please note that a pure $0$ can create a singularity. 

##### com
The $3$ values position of the `center of mass` relative to the local reference of the segment. The default values are `0 0 0`. 

##### inertia
The `inertia` tag allows to specify the matrix of inertia of the segment. It waits for $9$ values. The default values are the `identity matrix`

##### foreceplate or externalforceindex
When calculating the inverse dynamics, if force plates are used, this tag dispatch the force plates, the first force plate being $0$. If no force plate is acting on this segment, the value is $-1$. 

Warning: this tag MUST be added to a segment that has a translation and/or a rotation (i.e. that possesses at least one degree of freedom). Otherwise, it will simply be ignored

##### meshfile or ply
The path of the meshing `.bioBone`, `.ply`, `.stl` file respectively. It can be relative to the current running folder or absolute (relative being preferred) and UNIX or Windows formatted (`/` vs `\\`, UNIX being preferred).

##### mesh
If the mesh is not written in a file, it can be written directly in the segment. If so, the `mesh` tag stands for the vertex. Therefore, there are as many `mesh` tags as vertex. It waits for $3$ values being the position relative to reference of the segment. 

##### meshcolor
The color of the segment mesh given in RGV values `[0, 1]`. Default is `0.89, 0.855, 0.788`, that is bone color-ish.

##### meshscale
The scaling to apply to the provided mesh, given in `X Y Z` values. Default is `1 1 1`.

##### meshrt
The RT to apply to the provided mesh, given in `RX RY RZ seq TX TY TZ` as for RT. The default value is `0 0 0 xyz 0 0 0`. 


##### patch
The patches to define the orientation of the patches of the mesh. It waits for $3$ values being the $0-based$ of the index of the vertex defined by the `mesh`.

#### Marker
The marker with a unique name attached to a body segment. 

```c
marker my_marker
    parent segment_name
    position 0 0 0
    technical 1
    anatomical 0
endmarker
```

##### parent
The `parent` tag is the name of the segment that particular segment is attached to. The parent must be defined earlier in the file and is case dependent. This tag is waits for $1$ string.

##### position
The `position` of the marker in the local reference frame of the segment. 

##### technical
If the marker will be taged as technical (will be returned when asking technical markers). Default value is true ($1$).

##### anatomical
If the marker will be taged as anatomical (will be returned when asking anatomical markers). Default value is false ($0$).

##### axestoremove
It is possible to project the marker onto some axes, if so, write the name of the axes to project onto here. Waits for the axes in a string.

#### Imu
Same as a marker, but for inertial measurement unit. 
```c
imu my_imu
    parent segment_parent
    rtinmatrix 0
    rt 0 0 0 xyz 0 0 0
    technical 1
    anatomical 0
endimu
```
##### parent
The `parent` tag is the name of the segment that particular segment is attached to. The parent must be defined earlier in the file and is case dependent. This tag is waits for $1$ string.

##### rt
The `homogeneous matrix` of transformation (rototranslation) of the current segment relative to its parent. If `rtinmatrix` is defined to `1`, `rt` waits for a 4 by 4 matrix (the 3 x 3 matrix of rotation and the 4th column being the translation and the last row being 0 0 0 1); if it is defined to `0` it waits for the 3 rotations, the rotation sequence and the translation. The default value is the identity matrix.

##### rtinmatrix
The tag that defines if the `rt` is in matrix or not. If the `version` of the file is higher or equal than $3$, the default value is false ($0$), otherwise, it true ($1$). 

##### technical
If the marker will be taged as technical (will be returned when asking technical markers). Default value is true ($1$).

##### anatomical
If the marker will be taged as anatomical (will be returned when asking anatomical markers). Default value is false ($0$).

#### Contact
The position of a non acceleration point while computing the forward dynamics. 
```c
contact my_contact
    parent parent_segment
    position 0 0 0
    axis xyz
endcontact
```
##### parent
The `parent` tag is the name of the segment that particular segment is attached to. The parent must be defined earlier in the file and is case dependent. This tag is waits for $1$ string.

##### position
The `position` of the marker in the local reference frame of the segment. 

##### axis
The name of the `axis` that the contact acts on. If the version of the file is $1$, this tag has no effect. 

##### normal
The `normal` that the contact acts on. This tags waits for $3$ values with a norm $1$. If the version of the file is not $1$, this tag has no effect. To get the `x`, `y` and `z` axes, one must therefore define three separate contacts. 

##### acceleration
The constant `acceleration` of the contact point. The default values are `0, 0, 0`. 

#### Loopconstraint

##### 


# How to contribute

You are very welcome to contribute to the project! There are to main ways to contribute. 

The first way is to actually code new features for BIORBD. The easiest way to do so is to fork the project, make the modifications and then open a pull request to the main project. Don't forget to add your name to the contributor in the documentation of the page if you do so!

The second way is to open issues to report bugs or to ask for new features. I am trying to be as reactive as possible, so don't hesitate to do so!

# Graphical User Interface (GUI)
For now, there is no GUI for the C++ interface and the MATLAB one is so poor I decided not to release it. However, there is a Python interface that worths to have a look at. Installation procedure and documentation can be found at the GitHub repository (https://github.com/pyomeca/biorbd-viz).

# Documentation
The documentation is automatically generated using Doxygen (http://www.doxygen.org/). You can compile it yourself if you want (by setting `BUILD_DOC` to `ON`). Otherwise, you can access a copy of it that I try to keep up-to-date in the Documentation project of pyomeca (https://pyomeca.github.io/Documentation/) by selecting `biorbd`. 

# Troubleshoots
Despite my efforts to make a bug-free library, BIORBD may fails sometimes. If it does, please refer to the section below to know what to do. I will fill this section with the issue over time.

## Slow BIORBD
If you experience a slow BIORBD, you are probably using a non optimized version, that is compiled with `debug` level. Please use at least `RelWithDebInfo` level of optimization while compiling BIORBD. 

If you actually are using a released level of optimization, you may actually experiencing a bug. You are therefore welcomed to provide me with a minimal example of your slow code and I'll see how to improve the speed!

# Cite
If you use BIORBD, we would be grateful if you could cite it as follows:

```

@article{michaudBiorbd2021,
  title = {Biorbd: {{A C}}++, {{Python}} and {{MATLAB}} Library to Analyze and Simulate the Human Body Biomechanics},
  shorttitle = {Biorbd},
  author = {Michaud, Benjamin and Begon, Mickaël},
  date = {2021-01-19},
  journaltitle = {Journal of Open Source Software},
  volume = {6},
  pages = {2562},
  issn = {2475-9066},
  doi = {10.21105/joss.02562},
  url = {https://joss.theoj.org/papers/10.21105/joss.02562},
  urldate = {2021-01-19},
  abstract = {Michaud et al., (2021). biorbd: A C++, Python and MATLAB library to analyze and simulate the human body biomechanics. Journal of Open Source Software, 6(57), 2562, https://doi.org/10.21105/joss.02562},
  langid = {english},
  number = {57}
}
```

# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting one or more of the project maintainers.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.

# Contributing to `biorbd`
All contributions, bug reports, bug fixes, documentation improvements, enhancements, and ideas are welcome.
We recommend going through the list of [`issues`](https://github.com/pyomeca/biorbd/issues) to find issues that interest you, preferable those tagged with `good first issue`.
You can then get your development environment setup with the following instructions.

## Forking `biorbd

You will need your own fork to work on the code.
Go to the [biorbd project page](https://github.com/pyomeca/biorbd/) and hit the `Fork` button.
You will want to clone your fork to your machine:

```bash
git clone https://github.com/your-user-name/biorbd.git
```

## Creating and activating conda environment

Before starting any development, we recommend that you create an isolated development environment. 
The easiest and most efficient way (due to the numerous dependencies of `biorbd`) is to use an anaconda virtual environment and to create it based on the `environment_X.yml` file (where `X` is the desired backend, that is `eigen` or `casadi`). 

- Install [miniconda](https://conda.io/miniconda.html)
- `cd` to the biorbd source directory
- Install biorbd dependencies with:

```bash
conda env create -f environment_X.yml
conda activate biorbd_X
```
Where (again) `X` is the desired linear algebrea backend, that is `eigen` or `casadi`

## Implementing new features

Before starting to implement your new awesome feature, please discuss the implementation with the code owner to prevent any clashing with some other current developments. 
It is also a good idea to check the current opened pull-request to not redo something currently being developed. 
If your feature is mentioned in the issue section of GitHub, please assign it to yourself.
Otherwise, please open a new issue explaining what you are currently working on (and assign it to yourself!).

As soon as possible, you are asked to open a pull-request (see below) with a short but descriptive name. 
To tag that a pull-request is still in development, please add `[WIP]` at the beginning of the name.
Send commits that are as small as possible; 1 to 10 lines is probably a good guess, with again short but descriptive names. 
Be aware of the review done by the maintainers, they will contain useful tips and advices that should be integrated ASAP. 
Once you have responded to a specific comment, please respond `Done!` and tag it as resolved.

Make sure you create a unit test with numerical values for comparison.
As it may be tricky to create a numerical test for the CasADi's backend, you are welcome to ask for help as soon as you have created one for Eigen.  
Please note that if by accident you add a binary file in the history, your pull-request will be rejected and you will have to produce a new pull request free from the binary file. 

When you have completed the implementation of your new feature, navigate to your pull-request in GitHub and select `Pariterre` in the `Reviewers` drop menu. 
At the same time, if you think your review is ready to be merged, remove the `[WIP]` tag in the name (otherwise, your pull-request won't be merged). 
If your pull-request is accepted, there is nothing more to do. 
If changes are required, reply to all the comments and, as stated previously, respond `Done!` and tag it as resolved. 
Be aware that sometimes the maintainer can push modifications directly to your branch, so make sure to pull before continuing your work on the feature.

## Testing your code

Adding tests is required to get your development merged to the master branch. 
Therefore, it is worth getting in the habit of writing tests ahead of time so this is never an issue.
Each time you push to your pull-request, the `biorbd` test suite will run automatically.
However, we strongly encourage running the tests prior to submitting the pull-request.
If you configured your project using CMake and set `BUILD_TEST` to `ON`, the google test toolbox should download and compile itself.
You can thereafter run the test by running the binary `biorbd_tests` (the file will be at different place depending on the OS you are using).

## Convention of coding

`biorbd` tries to follow as much as possible a coherent standard, that is essentially:
  - camelCase;
  - 80 characters per line max;
  - the absence of `using namespace`
That said, if you are looking a the code, you will find unstandardize coding. 
That is the case because over the course of the past five years or so, the main coder ability has drastically changed, therefore historical non-compliant code remain, I appologize for that!
Therefore, as it is for now, unless you have good reasons to not follow a standardized approach, pull-requests are required to follow these recommendations. 

### All Submissions:

* [ ] Have you followed the guidelines in our Contributing document [doc/contribution.md]?
* [ ] Have you checked to ensure there aren't other open [Pull Requests] for the same update/change?
* [ ] Have you opened/linked the issue related to your pull request?
* [ ] Have you used the tag [WIP] for on-going changes, and removed it when the pull request was ready?

### New Feature Submissions:

1. [ ] Does your submission pass the tests for both Eigen and CasADi backends (if not please explain why this is intended)?
2. [ ] Have you linted your code locally prior to submission (using the command: `astyle -A10 -xC80 -z2 -xW --exclude="external" -n --recursive *.cpp *.h`)? 

### Changes to Core Features:

* [ ] Have you added an explanation of what your changes do and why you'd like us to include them?
* [ ] Have you written new tests for your core changes, as applicable?

---
title: 'biorbd: A C++, Python and MATLAB library to analyze and simulate the human body biomechanics'
tags:
  - Python
authors:
  - name: Benjamin Michaud
    orcid: 0000-0002-5031-1048
    affiliation: 1
  - name: Mickaël Begon
    orcid: 0000-0002-4107-9160
    affiliation: 1
affiliations:
 - name: École de Kinésiologie et de Sciences de l'Activité Physique, Université de Montréal
   index: 1
date: May 1st, 2020
bibliography: paper.bib
---

# Summary
Biomechanics is at the interface of several fields of science, such as mechanics, human physiology and robotics.
Although this transdisciplinarity encourages the emergence of new ideas, the variety of data to analyze simultaneously can be overwhelming.
Commonly biomechanical datasets are composed of skin markers trajectories (termed as markers), contact forces, electromyography (EMG) signal, inertial measurement units (IMU) kinematics, etc., which by nature are not straightforward to combine.
It is at their meeting point---the body movement---that `biorbd` steps in; `bio` standing for biomechanics and `rbd` for `rigid body dynamics`.
`biorbd` is a *feature-based development* library that targets the manipulation of biomechanical data in a comprehensive and accessible manner.
For a given musculoskeletal model, it provides functions for inverse flow---i.e., from markers to EMG---and direct flow---i.e., from EMG to markers.

Since biomechanics often requires computationally expensive or real-time computations,
the core of `biorbd` is written in C++.
Although this language provides fast computations, it lacks the flexibility of higher-level languages.
To meet the needs of the biomechanics community, Python and MATLAB binders are provided with `biorbd`.
As a result, `biorbd` can elegantly be implemented to common workflows of researchers without compromising the required speed.

Finally, biomechanical data are often multidimensional and almost always time-dependent which can be challenging to visualize.
To help with that, `bioviz` [@Michaud2018bioviz], a Python visualizer, was purposely designed.
This visualizer allows animating the model, record videos, and, for models that include muscles, plot muscular outputs against various features of the movement.

# A `biorbd` overview, the inverse and direct flow
Biomechanical analyses are usually based on one (or a mixture) of the inverse or direct flow [@kainzJointKinematicCalculation2016].
Briefly, the former uses measurements from a movement (e.g., markers) and infers its cause, while the latter assumes control (e.g., EMG) and outputs the resulting kinematics.

## Inverse flow
*Inverse kinematics*: Estimates the generalized coordinates ($q$)—i.e., the body kinematics—from body sensor measurements (e.g., markers, IMU, etc.).
The main algorithm implemented is the Extended Kalman Filter [@fohannoEstimation3DKinematics2010a] which by design facilitates the merging of multiple data sources and takes care of missing data.

*Inverse dynamics*: Estimates the generalized forces ($\tau$) producing a given generalized acceleration ($\ddot{q}$) (the second time derivative of $q$):
$$
\tau = M(q)\ddot{q} + N(q, \dot{q})
$$
where $\dot{q}$ is the generalized velocities, $M(q)$ is the mass matrix and $N(q, \dot{q})$ are the nonlinear effects.

*Static optimization*: Estimates the muscle activations ($\alpha$) producing a given $\tau$ [@andersonStaticDynamicOptimization2001].
It minimizes the muscle activation *p*-norm ($p$ usually being $2$) that matches a given $\tau$ using nonlinear optimization [Ipopt, @wachterImplementationInteriorpointFilter2006].
$$
\begin{aligned}
    & \underset{\alpha \in \mathbb{R}^m}{\text{minimize}}
    & & \left\lVert\alpha\right\rVert_p \\
    & \text{subject to}
    & & \tau_{mus_i}(\alpha ,q, \dot{q}) - \tau_{kin_i}(q, \dot{q}, \ddot{q}) = 0, &\; i=1,\ldots,n \\
    & & &  0 \leq \alpha_{t_j} \leq 1, &\; j=1,\ldots,m
\end{aligned}
$$
where $\tau_{mus_i}(\alpha ,q, \dot{q})$ and $\tau_{kin_i}(q, \dot{q}, \ddot{q})$ are $\tau$ computed from muscle forces ($F_{mus}(\alpha, q, \dot{q})$) and inverse dynamics, respectively.

## Direct flow
*Muscle activation dynamics*: Estimates the muscle activation derivative ($\dot{\alpha}$) from the muscle excitation---that is the calcium release in the muscle that triggers the muscle contraction.
Multiple activation/excitation dynamics are implemented [e.g., @manalOneparameterNeuralActivation2003;@thelenAdjustmentMuscleMechanics2003].

*Muscular joint torque*: Estimates the $\tau_{mus}$ from muscle forces ($F_{mus}(q, \dot{q}, \alpha)$) [@shermanHowComputeMuscle2010], estimated from $\alpha$ using a muscle model [e.g., @hillHeatShorteningDynamic1938; @thelenAdjustmentMuscleMechanics2003]:
$$
\tau_{mus} = J_{mus}(q)^T F_{mus}(q, \dot{q}, \alpha)
$$
where $J_{mus}(q)$ is the muscle lengths Jacobian.

*Forward dynamics*: Estimates the $\ddot{q}$ from a given $\tau$:
$$
\ddot{q} = M(q)^{-1}\tau - N(q, \dot{q})
$$
All the forward dynamics implemented in `RBDL` [@felisRBDLEfficientRigidbody2017] are available.

*Forward kinematics*: Estimates the model kinematics outputs (e.g., markers, IMU) from a given $q$, after integrating twice $\ddot{q}$.

# The dependencies
`biorbd` takes advantage of efficient back ends, especially  the `RBDL` and `CasADi` libraries.
`RBDL`, written by Martin Feliz [@felisRBDLEfficientRigidbody2017], implements Featherstone equations of spatial geometry [@featherstoneRobotDynamicsEquations2000], successfully used in the field of robotics [@macchiettoMomentumControlBalance2009; @diehlFastDirectMultiple2006; @kurfessRoboticsAutomationHandbook2018]. 
`RBDL` provides the computational core for body dynamics.
`biorbd` extends `RBDL` by giving commonly used biomechanics nomenclature, and by adding biomechanical modules, amongst others. 
`RBDL` is based on the highly efficient C++ linear algebra library `Eigen` [@guennebaud2010eigen].
Although `Eigen` is flexible and fast enough for most of the common usage, it cannot automatically provide derivatives of functions.
Therefore, `RBDL` was also augmented with the algorithmic differentiation library `CasADi` [@anderssonCasADiSoftwareFramework2019].
`CasADi` allows computing at low cost the derivatives of almost all the functions in `RBDL` and `biorbd`.
This is particularly useful when using `biorbd` in a gradient-based optimization setting.

# Statement of need
`OpenSim` [@sethOpenSimSimulatingMusculoskeletal2018] and `Anybody` [@damsgaardAnalysisMusculoskeletalSystems2006] are state-of-the-art biomechanics software that provides similar analysis flows with advanced user interface.
`Anybody` being a closed and proprietary software, the reason to create another library for the open-source community is self-explanatory.
Conversely, `OpenSim` is open-source and well established in the biomechanics community.

Nevertheless, in line with the idea that simulation software in biomechanics should be validated in multiple ways [@hicksMyModelGood2015], providing similar tools but different in their approach allows the community to cross-validate the different implementation of the algorithms.
For instance, two papers [@kimSimilaritiesDifferencesMusculoskeletal2018; @trinlerMuscleForceEstimation2019a] recently compared the outputs of `Anybody` and `OpenSim` and came to different results.
Although the authors provided plausible explanations for these differences, due to the closed-source nature of `Anybody`, they had to assume that the implementation of the algorithms are flawless in both software.
However, since a direct comparison between the actual codes is impossible, this is not verifiable.
Having multiple open source software that produces similar ends by different means is a quality assurance for the end users: "Do not put all your eggs in one basket.”
To the best of our knowledge, there is no other open-source software that provides a complete direct and inverse flow in biomechanics. 
Therefore, in our opinion, `biorbd` and `OpenSim` are complementary.

# Previous usage of `biorbd`
`biorbd` was used in most of the project of the *Laboratoire de Simulation et Modélisation du Mouvement* (S2M); particularly in analysis settings [@jacksonImprovementsMeasuringShoulder2012; @desmyttereEffect3DPrinted2020; @verdugoEffectsTrunkMotion2020] and simulation settings [@belaiseEMGmarkerTrackingOptimisation2018; @moissenetOptimizationMethodTracking2019] for a wide variety of movements (walking, piano playing, upper limb maximal exertions, etc.)
More recently, an optimal control framework for biomechanics [`bioptim`, @Michaud2018bioptim] based on Ipopt [@wachterImplementationInteriorpointFilter2006] and ACADOS [@Verschueren2019] was developed around `biorbd`.

# Acknowledgements
A huge thanks to Ariane Dang for her patience and contribution to writing the tests for `biorbd`!

# References
