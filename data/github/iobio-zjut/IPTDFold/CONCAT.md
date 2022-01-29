# IPTDFold
### A de novo protein structure prediction by iterative partition sampling, topology adjustment, and residue-level distance deviation optimization

**Developer:**   
                Jun Liu  
                College of Information Engineering  
                University of Zhejiang University of Technology, Zhejiang  
                Email: junl@zjut.edu.cn  
		
**Contact:**  
                Guijun Zhang, Prof  
                College of Information Engineering  
                University of Zhejiang University of Technology, Zhejiang  
                Email: zgj@zjut.edu.cn  

## 1. INSTALLATION
Binaries for Linux 64 bit system has been included in the package. The Linux binary was compiled using GCC 5.4.0. Users need to have these versions of GCC compilers when using binaries.

Please Follow the below steps to install and configure IPTDFold:

- Download Rosetta3.10 source package from https://www.rosettacommons.org/software/ 
and extract it to ``"~/"`` directory.
- Compile the source code of Rosetta using the following commands:

```
 $ cd ~/Rosetta/main/source/
 $ ./scons.py -j<NumOfJobs> mode=release bin
``` 

- Copy and paste source code of ``"ClassicAbinitio.cc"``, ``"ClassicAbinitio.hh"``, ``"LJcore.cc"``, ``"LJcore.hh"``,  ``"LJAngleRotation.cc"``, ``"LJAngleRotation.hh"``, and ``"LJAngleRotation.fwd.hh"`` from ``"src/"`` folder in IPTDFold package to ``"~/Rosetta/main/source/src/protocols/abinitio/"`` folder in Rosetta. Copy and paste configuration file ``"protocols_b_6.src.settings"`` from ``"src/"`` folder in IPTDFold package to ``"~/Rosetta/main/source/src/"`` folder in Rosetta.

- Compile Rosetta source code using the following commands:

```
 $ cd ~/Rosetta/main/source/
 $ ./scons.py -j<NumOfJobs> AbinitioRelax mode=release
```

- If you want to recompile IPTDFold source code, use the following commands:

```
 $ cd ~/IPTDFold/
 $ g++ -o bin/IPTDFold src/IPTDFold.cpp
```
- If you do not recompile IPTDFold source code, you need to add executable permissions for file of ``"~/IPTDFold/bin/IPTDFold"`` using the following commands:

```
 $ chmod +x ~/IPTDFold/bin/IPTDFold/
```
## 2. INPUT
IPTDFold requires five files to generate models:

	-f	fasta			 : fasta file
	-d	dmap			 : inter-residue distance file
	-frag3	3mer_fragment_library 	 : fragment library with fragment lenth 3
	-frag9	9mer_fragment_library	 : fragment library with fragment lenth 9
	-rosetta rosetta executable file : executable file of ClassicAbinitio protocol of Rosetta

## 3. OUTPUT
The predicted model will be generated on the path where the IPTDFold program runs.

## 4. EXAMPLE
Please follow the below steps to run IPTDFold:

- Go to the ``"~/IPTDFold/example/"`` folder of IPTDFold.
  
- Run IPTDFold with the following command (the path of '-rosetta' needs to be replaced with the actual path of the Rosetta executable file.):
  
```
      $ ../bin/IPTDFold -f fasta -d distance.txt -frag3 3mer_frag_set -frag9 9mer_frag_set -rosetta ~/rosetta_src_2018.33.60351_bundle/main/source/bin/AbinitioRelax.default.linuxgccrelease
```

- Five models will generate in the ``"example/"`` folder.

## 5. DISCLAIMER
The executable software and the source code of IPTDFold is distributed free of charge 
as it is to any non-commercial users. The authors hold no liabilities to the performance 
of the program.
