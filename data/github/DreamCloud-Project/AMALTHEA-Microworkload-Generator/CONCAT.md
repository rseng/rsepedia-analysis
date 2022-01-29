# AMALTHEA-Microworkload-Generator
Command line tool generating micro workloads from an AMALTHEA application.

To get the tool you must clone this repository and its submodules. To clone the repository if you have a gite.lirmm.fr account with an SSH key registered use `git clone git@github.com:DreamCloud-Project/AMALTHEA-Microworkload-Generator.git`. Else use `git clone https://github.com/DreamCloud-Project/AMALTHEA-Microworkload-Generator.git`. Then use `cd AMALTHEA-Microworkload-Generator` followed by `git submodule init` followed by `git submodule update` to clone submodules.

## Using the tool

To ease compiling the tool a compile.py script is provided. 

The requirements are the following ones:  

- have CMake installed on your system [(https://cmake.org)](https://cmake.org/)
- have the xerces-c-dev library installed in standard includes and libs folders (using apt-get for example)
  or have xerces-c-dev library in a custom folder and define XERCES_HOME
  
## Licence

This software is made available under the  GNU Lesser General Public License v3.0

Report bugs at: mcsim-support@lirmm.fr  

(C)2016 CNRS and Université de Montpellier
