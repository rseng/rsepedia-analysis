# SIMsalabim Project

SIMsalabim: A 1D drift-diffusion simulator for semiconductor devices (LEDs, solar cells, diodes, organics, perovskites). It consists of two programs that share most of their code: SimSS (simulates steady-state), and ZimT (zimulates transients).

## Table of Contents
1. [Introduction](#introduction)
2. [Quickstart guide](#quickstart-guide)
3. [Instructional videos](#instructional-videos)
4. [Running SIMsalabim](#running-simss-or-zimt)
5. [How to cite](#how-to-cite)
6. [Copyright and license](#copyright-and-license)
7. [How to contribute](#how-to-contribute)
8. [Scientific publications based on SIMsalabim](#scientific-publications-based-on-the-simsalabim-project)

## Introduction

SimSS and ZimT can be used to simulate current-voltage (JV) characteristics of semiconductor devices. This includes the effects of generation, recombination and trapping of electrons and holes, the effect of ions and dopants, and self-consistently solves the electric field that results from all charged species. It is up to the user to make sure that the physical model that underlies SIMsalabim is applicable (see Manual for more details).

The routines that make up this project are used in two different guises: SimSS and ZimT. The former is a steady-state version and is the 'main' code. The latter, ZimT, can be used to simulate transients. 
The project folder is structured as follows:
- Docs: this folder contains a manual, developer guidelines, and a log of all changes to the code
- SimSS: contains the steady-state version of this software
- Tests: contains a number of tests to assess the functioning of the codes (SimSS and ZimT)
- Units: this is where the heart of the code resides. These units are shared between SimSS and ZimT
- ZimT: the transient equivalent of SimSS


## Quickstart guide
There are two ways to get SimSS or ZimT running on your machine:
1. [Use pre compiled binaries](#Pre-compiled-binaries)
2. [Compile SIMsalabim](#Compiling-SIMsalabim)

While using the pre-compiled binaries is easiest, they will not work on every machine. If they do not work follow the instructions for compiling SIMsalabim.

### Pre-compiled binaries

SimSS and ZimT come as pre-compiled binaries for WIN or Linux (64 bits), see [releases](https://github.com/kostergroup/SIMsalabim/releases). This avoids having to download and install the FPC compiler and compiling the code which is described below. There is no guarantee that these will work, but it is worth a try. Simply copy the binary to the main SimSS/ZimT folder. The Linux one might require a change of permissions to run, something like:
<pre><code>
    chmod +x simss
 </code></pre>

### Compiling SIMsalabim

There are 3 steps to complete to run SimSS or ZimT for the first time:
  1. Install the "Free Pascal" compiler.
  2. Compile SimSS (or ZimT).
  3. Run SimSS (or ZimT). 

SimSS and ZimT should compile and run on any platform supported by the "Free Pascal" compiler. We provide detailed instructions for installing both the compiler and running the code for Windows and Linux (Ubuntu). For any other platform, please install the compiler (version 3.2.0 or newer) following the instructions from https://www.freepascal.org/. Steps 2 and 3 (compile and run) should be similar to the Ubuntu and Windows steps outlined below.

Below are the instructions for performing these steps on Ubuntu (or any other Debian based Linux distribution) and Windows.  


### Ubuntu (Linux) guide

Navigate to the SimSS or ZimT folder, where you can find simss.pas or zimt.pas. This example assumes you want to compile SimSS, but you can replace `simss` with `zimt` in the commands and run from the ZimT folder to compile and run ZimT. 

**Ubuntu 20.10 or newer:**
Steps:
  - Open a terminal in the SimSS folder
  - run the following commands to perform step 1, 2, and 3 (install, compile, run):
```
sudo apt install fpc
fpc simss
./simss
```

**Ubuntu 20.04 or older:**
We need a relatively recent version of the Pascal compiler (v 3.2.0 or higher). Therefore, manual installation of the compiler is required on Ubuntu 20.04 or older.
Steps:
  - Download the relevant version of the compiler from https://www.freepascal.org/
  - Install the compiler. This is done by running the installation script that can be found in the downloaded archive. While it can be installed without root privilege, the compiler should be manually added to path making installation a bit more tricky. To install as root run the command: `sudo ./install.sh`. 
  - Open a terminal in the SimSS folder
  - run the following commands to perform step 2 and 3 (compile and run):
```
fpc simss
./simss
```

**Installing on other Linux distributions:**
Most package managers for other Linux distributions ship the compiler package in a recent enough version (version 3.2.0 or newer of package *fpc*), so it suffices to exchange the *apt install* command with another package manager's equivalent. 

### Windows guide

Navigate to the SimSS or ZimT folder, where you can find simss.pas or zimt.pas. This example assumes you want to compile SimSS, but you can replace `simss` with `zimt` in the commands and run from the ZimT folder to compile and run ZimT. 

Steps:
  - Download and install fpc from https://www.freepascal.org/ (SimSS and ZimT require version 3.2.0 or newer)
  - Open a command prompt or powershell instance in the SimSS folder (for example by opening the folder in explorer, pressing <kbd>alt</kbd>+<kbd>D</kbd> to select the location bar, typing 'cmd' and pressing <kbd>enter</kbd>).
  - run the following commands to perform step 2 and 3 (compile, run):
  
```
fpc simss
simss.exe
```

Alternative:
  - Install Windows subsystem for Linux.
  - Follow the Ubuntu instructions above.


### Miscellaneous

**Warning message**
You may see this warning message when compiling the codes:
/usr/bin/ld.bfd: warning: link.res contains output sections; did you forget -T?
Just ignore this.

**Editing the code**
You'll probably want to use (install?) some form of integrated development environment (IDE). Most modern IDE's include syntax highlighting for the Pascal language. Popular IDEs include Lazarus, Geany, and Visual studio code (requires extension to perform syntax highlighting for pascal). Some guidelines for extending the code can be found in Docs/Developer_guidelines.md.


## Instructional videos

In order to help new users to get the most out of SimSS/ZimT, we have recorded a few instructional videos: 

1 [Installation and basic use](https://youtu.be/0mtpGJMnbFE "Installation and basic use")

Please note that at the time of recording SimSS was called SIMsalabim, so the project name and the name of one of the programs was the same.


## Running SimSS or ZimT
- all parameters are specified in device_parameters.txt (in ZimT an additional time-voltage-generation (tVG) file is present to indicate what time, generation, and voltage at every simulation step).
- on Linux run SimSS/ZimT by entering in the terminal from the folder where the compiled program is located and run using
<pre><code>
./simss
</code></pre>
or 
<pre><code>
./zimt
</code></pre>
on windows:
<pre><code>
simss.exe
</code></pre>
or 
<pre><code>
zimt.exe
</code></pre>
Alternatively, the programs can be double clicked in the file manager on both Windows and Linux to run.

- all parameters listed in device_parameters.txt can also be changed via the command line. This will override the respective parameter value in device_parameters.txt. Simply add the parameter name and value to the command line after a dash (-).
Example: change of thickness (L) and JV output file (JV_file):
<pre><code>
./simss -L 345E-9 -JV_file anotherJV.dat
</code></pre>
- multiple output files will be generated (see device_parameters) and log file.

4) Other remarks
- in general, input files and the log file have extension .txt, whereas output files that can be used for plotting (J-V curves, for example) have extension .dat.
- all changes are shown and commented on in the file Docs/Change_log.txt. Here we document not just the changes, but also explain and motivate some of the choices in naming, physical models, etc. 


## How to cite

The open-source version of the code has been published as:

M. Koopmans, V.M. Le Corre, and L.J.A. Koster, SIMsalabim: An open-source drift-diffusion simulator for semiconductor devices, J. Open Source Softw. **7**, 3727 (2022).

[The paper can be downloaded here.![DOI](https://joss.theoj.org/papers/10.21105/joss.03727/status.svg)](https://doi.org/10.21105/joss.03727)


## Copyright and license

The SIMsalabim project is licensed under the GNU Lesser General Public Licence version 3. The details of this licence can be found in the files COPYING.txt and COPYING.LESSER.txt. 
Several authors (all from the University of Groningen) have contributed to the code: 
- Dr T.S. (Tejas) Sherkar
- Dr V.M. (Vincent) Le Corre
- M. (Marten) Koopmans
- F. (Friso) Wobben
- Prof. Dr. L.J.A. (Jan Anton) Koster


## How to contribute
If you would like to would like to contribute to the SIMsalabim project or have any questions, please read the [Developer instructions](https://github.com/kostergroup/SIMsalabim/blob/master/Docs/Developer_guidelines.md). It covers how we would like to approach both questions and possible code changes.


## Scientific publications based on the SIMsalabim project

List of publications (not complete):

- V.M. Le Corre, T.S. Sherkar, M. Koopmans, and L.J.A. Koster, Identification of the Dominant Recombination Process for Perovskite Solar Cells Based on Machine Learning, Cell Rep. Phys. Sci. 2, 100346 (2021).

- Y. Firdaus, V.M. Le Corre, S. Karuthedath, W. Liu, A. Markina, W. Huang, S. Chattopadhyay, M.M. Nahid, M.I. Nugraha, A. Seitkhan, A. Basu, Y. Lin, I. McCulloch, H. Ade, J. Labram, F. Laquai, D. Andrienko, L.J.A. Koster, and T.D. Anthopoulos, Long-range exciton diffusion in molecular non-fullerene acceptors, Nature Comm. 11, 5220 (2020).

- D. Hu, Q. Yang, H. Chen, F. Wobben, V.M. Le Corre, R. Singh, L.J.A. Koster, Z. Kan, Z. Xiao, and S. Lu, 15.34% Efficiency All-Small-Molecule Organic Solar Cells with Improved Fill Factor Enabled by a Fullerene Additive, Energy Environ. Sci. 13, 2134 (2020).

- L. Hou, J. Lv, F. Wobben, V.M. Le Corre, H. Tang, R. Singh, M. Kim, F. Wang, H. Sun, W. Chen, Z. Xiao, M. Kumar, T. Xu, W. Zhang, I. McCulloch, T. Duan, H. Xie, L.J.A. Koster, S. Lu, and Z. Kan, Effects of Fluorination on Fused Ring Electron Acceptor for Active Layer Morphology, Exciton Dissociation, and Charge Recombination in Organic Solar Cells, ACS Appl. Mater. Interfaces 12, 56231 (2020).

- E.A. Duijnstee, V.M. Le Corre, M.B. Johnston, L.J.A. Koster, J. Lim, and H.J. Snaith, Understanding dark current-voltage characteristics in metal-halide perovskite single crystals, Phys. Rev. Appl. 15, 014006 (2021).

- D. Neher, J. Kniepert , A. Elimelech, and L.J.A. Koster, A new Figure of Merit for Organic Solar Cells with Transport-limited Photocurrents, Sci. Rep. 6, 24861 (2016).

- S. Shao, M. Abdu-Aguye, T.S. Sherkar, H.-H. Fang, G. ten Brink, B.J. Kooi, L.J.A. Koster, and M.A. Loi, The effect of the microstructure on trap-assisted recombination and light soaking phenomenon in hybrid perovskite solar cells, Adv. Funct. Mater. 26, 8094 (2016).

- T.S. Sherkar, C. Momblona, L. Gil-Escrig, H.J. Bolink, and L.J.A. Koster, Improving perovskite solar cells: Insights from a validated device model, Adv. Energy Mater. 1602432 (2017).
 
- T.S. Sherkar, C. Momblona, L. Gil-Escrig, J. Ávila, M. Sessolo, H. Bolink, and L.J.A. Koster, Recombination in Perovskite Solar Cells: Significance of Grain Boundaries, Interface Traps and Defect Ions, ACS Energy Lett. 2, 1214 (2017).

- V.M. Le Corre, A. Rahimi Chatri, N.Y. Doumon, and L.J.A. Koster, Charge carrier extraction in organic solar cells governed by steady-state mobilities, Adv. Energy Mater. 1701138 (2017).

- F.J.M. Colberts, M.M. Wienk, R. Heuvel, W. Li, V.M. Le Corre, L.J.A. Koster, and R.A.J. Janssen, Bilayer-Ternary Polymer Solar Cells Fabricated Using Spontaneous Spreading on Water, Adv. Energy Mater. 8, 1802197 (2018).

- N.Y. Doumon, M.V. Dryzhov, F.V. Houard, V.M. Le Corre, A. Rahimi Chatri, P. Christodoulis, and L.J.A. Koster, Photostability of Fullerene and Non-Fullerene Polymer Solar Cells: The Role of the Acceptor, ACS Appl. Mater. Interfaces 11, 8310 (2019).

- V.M. Le Corre,  M. Stolterfoht, L. Perdigón Toro, M. Feuerstein, C. Wolff, L. Gil-Escrig, H.J. Bolink, D. Neher, and L.J.A. Koster, Charge transport layers limiting the efficiency of perovskite solar cells: how to optimize conductivity, doping and thickness, ACS Appl. Energy Mater. 2, 6280 (2019). 

- E. A. Duijnstee, J.M. Ball, V.M. Le Corre, L.J.A. Koster, H.J. Snaith, and J. Lim, Towards Understanding Space-charge Limited Current Measurements on Metal Halide Perovskites, ACS Energy Lett. 5, 376 (2020).

- D. Bartesaghi and L. J. A. Koster, The effect of large compositional inhomogeneities on the performance of organic solar cells: A numerical study, Adv. Funct. Mater. 25, 2013 (2015).

- J. Kniepert, I. Lange, J. Heidbrink, J. Kurpiers, T. Brenner, L.J.A. Koster, and D. Neher, Effect of Solvent Additive on Generation, Recombination and Extraction in PTB7:PCBM Solar Cells: A conclusive Experimental and Numerical Simulation Study, J. Phys. Chem. C 119, 8310 (2015).

- D. Bartesaghi, I. del Carmen Pérez, J. Kniepert, S. Roland, M. Turbiez, D. Neher, and L.J.A. Koster, Competition between recombination and extraction of free carriers determines the fill-factor of organic solar cells, Nature Comm. 6, 7083 (2015).

- N.J. van der Kaap, I. Katsouras, K. Asadi, P.W.M. Blom, L.J.A. Koster, and D.M. de Leeuw, Charge transport in disordered semiconducting polymers driven by nuclear tunneling, Phys. Rev. B 93, 140206(R) (2016).

- D. Bartesaghi, M. Turbiez, and L. J. A. Koster, Charge Transport and Recombination in PDPP5T:[70]PCBM Organic Solar Cells: the Influence of Morphology, Org. Elec. 15, 3191 (2014).

- G. A. H. Wetzelaer, N. J. van der Kaap, L. J. A. Koster, and P. W. M. Blom, Quantifying bimolecular recombination in organic solar cells in steady-state, Adv. Energy Mater. 3, 1130 (2013).

- J. Kniepert, I. Lange, N. J. van der Kaap, L. J. A. Koster, and D. Neher, A conclusive view on charge generation, recombination and extraction in as-prepared and annealed P3HT:PCBM blends: a combined experimental-simulation work, Adv. Energy Mater. 4, 1301401 (2014).

- L. J. A. Koster, M. Kemerink, M. M. Wienk, K. Maturová, and R. A. J. Janssen, Quantifying bimolecular recombination losses in bulk heterojunction solar cells, Adv. Mater. 63, 1670 (2011).

- G. A. H. Wetzelaer, L. J. A. Koster, and P. W. M. Blom, Validity of the Einstein relation in disordered organic semiconductors, Phys. Rev. Lett. 107, 066605 (2011).

- L.J.A. Koster, S. E. Shaheen, and J. C. Hummelen, Pathways to a new efficiency regime for organic solar cells, Adv. Energy Mater. 2, 1246 (2012). Listed in Adv. Energy Mater. top-5 most accessed articles in May 2012
- D.J. Wehenkel, L.J.A. Koster, M. M. Wienk, and R. A. J. Janssen, Influence of injected charge carriers on photocurrents in polymer solar cells, Phys. Rev. B 85, 125203 (2012).
- M. Kuik, L. J. A. Koster, A. G. Dijkstra, G. A. H. Wetzelaer, and P. W. M. Blom, Non-radiative recombination losses in polymer light-emitting diodes, Org. Elec. 13, 969 (2012).

- M. Kuik, L. J. A. Koster, G. A. H. Wetzelaer, and P. W. M. Blom, Trap-assisted recombination in disordered organic semiconductors, Phys. Rev. Lett. 107, 256805 (2011).

- L. J. A. Koster, Charge carrier mobility in disordered organic blends for photovoltaics, Phys. Rev. B 81, 205318 (2010).

- M. M. Mandoc, W. Veurman, L. J. A. Koster, M. M. Koetse, J. Sweelssen, B. de Boer, and P. W. M. Blom, Charge transport in MDMO-PPV : PCNEPV all-polymer solar cells, J. Appl. Phys. 101, 104512 (2007).

- V. D. Mihailetchi, H. Xie, B. de Boer, L. J. A. Koster, and P. W. M. Blom, Charge Transport and Photocurrent Generation in Poly(3-hexylthiophene):Methanofullerene Bulk-Heterojunction Solar Cells, Adv. Funct. Mater. 16, 699 (2006).L. J. A. Koster, W. J. van Strien, W. J. E. Beek, and P. W. M. Blom, Device operation of conjugated polymer/zinc oxide bulk heterojunction solar cells, Adv. Funct. Mater. 17, 1297 (2007). 

- M. M. Mandoc, L. J. A. Koster, and P. W. M. Blom, Optimum charge carrier mobility in organic solar cells, Appl. Phys. Lett. 90, 133504 (2007).

- M. M. Mandoc, W. Veurman, L. J. A. Koster, B. de Boer, and P. W. M. Blom, Origin of the Reduced Fill Factor and Photocurrent in MDMO-PPV:PCNEPV All-Polymer Solar Cells, Adv. Funct. Mater. 17, 2167 (2007).J. D. Kotlarski, M. Lenes, L. J. A. Koster, L. H. Slooff, and P. W. M. Blom, Combined optical and electrical modeling of polymer:fullerene bulk heterojunction solar cells, J. Appl. Phys. 103, 084502 (2008).

- L. J. A. Koster, V. D. Mihailetchi, and P. W. M. Blom, Ultimate efficiency of polymer/fullerene bulk heterojunction solar cells, Appl. Phys. Lett. 88, 093511 (2006).

- M. Lenes, L. J. A. Koster, V. D. Mihailetchi, and P. W. M. Blom, Thickness dependence of the efficiency of polymer:fullerene bulk heterojunction solar cells, Appl. Phys. Lett. 88, 243502 (2006).

- L. J. A. Koster, V. D. Mihailetchi, and P. W. M. Blom, Bimolecular recombination in polymer/fullerene bulk heterojunction solar cells, Appl. Phys. Lett. 88, 052104 (2006).
 
- V. D. Mihailetchi, H. Xie, B. de Boer, L. J. A. Koster, L. M. Popescu, J. C. Hummelen, and P. W. M. Blom, Origin of the enhanced performance in poly(3-exylthiophene):methanofullerene solar cells using slow drying, Appl. Phys. Lett. 89, 012107 (2006).
  
- L. J. A. Koster, E. C. P. Smits, V. D. Mihailetchi, and P. W. M. Blom, Device model for the operation of polymer/fullerene bulk heterojunction solar cells, Phys. Rev. B 72, 085205 (2005).
  
- L. J. A. Koster, V. D. Mihailetchi, R. Ramaker, and P. W. M. Blom, Light intensity dependence of open-circuit voltage of polymer:fullerene solar cells, Appl. Phys. Lett. 86, 123509 (2005).
   
- L. J. A. Koster, V. D. Mihailetchi, H. Xie, and P. W. M. Blom, Origin of the light intensity dependence of the short-circuit current of polymer/fullerene solar cells, Appl. Phys. Lett. 87, 203502 (2005).






# SIMsalabim project: Developer Instructions

All input that can help us make SIMsalabim better and more useful for users is welcome. This document contains instructions for what to do if you would like to contribute to or tinker around with the code.

## Filing Issues
If anything is unclear or you have any feature requests, you can file an 'Issue' in the Issues section on GitHub. Please don't hesitate to do this, as it might be of help for other user or us as maintainers as well. If you have a question we can provide help in the issues section so that future users can use if for reference. 

If you have any feature requests, we would first need to discuss the how and why in detail, which we would also do in the issues section. We will decide on a plan of action on a per case basis, but expect to make many of the code changes ourselves. 

In case of bug reports we will need the exact configuration of the 'device_parameters.txt' file as many appendant bugs can be traced back to settings that result in an unphysical device. Including that, along with the terminal output (including version of the program used), will help us dial in on the problem quickly.


## Overall structure
- ZimT and SimSS share much of the same code, so try to put new code in the units, not in ZimT or SimSS
- Routines, types and constants that are specific to drift-diffusion modelling should go in DDRoutines or DDTypesAndConstants. These units should have the same version number as the ZimT/SimSS codes.
- Somewhat generic routines, types and constants should go into units TypesAndConstant, InputOutputUtils, and NumericalUtils. As yet, these units do not have a version number.
- Input parameters and parameters that are directly derived from those input parameters (like Booleans based on a 0,1 input) are stored in par (='parameters') of type TInputParameters. The order of the fields within this record (integers, reals, etc are grouped) should follow the order of the variables in the file device_parameters.txt. Note: some parameters are limited to either SimSS or ZimT. However, we do not make use of the possibility of introducing a variant part in these records as this led to unexpected behaviour. For example, when a field is only defined (through the variant part) to ZimT, it was still availalbe in SimSS and no warning or exception was given when using that part of the record (with an undefined value!).
- Other variables that remain do not change throughout the simulation are stored in stv (='static variables') of type TStaticVars.

## General remarks on changing the code
- Make sure you understand what is permitted&mdash;and what is not&mdash;under the GLPL licence.
- Describe changes in file change_log.txt: Briefly describe the changes to the code, but this is not the main point (as one can easily figure them out using a difference viewer). More importantly, motivate and/or explain why the change was made.
- Run the tests as outlined in the folder 'Tests'. These are just a few tests to assess the basic functionality so additional testing is required.
- Magic numbers used in the code are, mostly, listed as constants in DDTypesAndConstants. Here we also indicate where&mdash;in which routine&mdash;the magic number is used and what it does.
-Try to limit local variables in and parameters passed to subroutines to 32K as this is the limit of some processors. 

## Notation within the codes
- Use Linux line-endings
- Add ample comments to the code, using curly brackets {}. Describe why rather than when formulating comments.
- Naming of variables: camel case
- Naming of procedures and function: pascal case combined with snake-case: Calc_Elec_Mob, for example
- FPC/Pascal key words fully capitalized: IF ... THEN ... ELSE
- Indentation width: 4

## Parameters
- Physical units: SI only, except for work functions (eV). Work functions are given as the distance to vacuum so they are positive.
- Device parameters: try to keep the same order for ZimT and SimSS. Every parameter should have a unit (if any) and a short description. After reading the parameters, they should be checked (to some extent) whether they make sense (logical, physical) in Check_Parameters.
# SIMsalabim Tests

The literature on semiconductor physics provides many analytical expressions that hold in certain limits that can be reproduced by SimSS/ZimT. These expressions therefore serve as highly suitable validation tests for the codebase. We chose the tests in such a way that most of the functionality is covered by the test cases. In this document, we describe the tests, while we provide the test results in graph form and the used simulation parameters in the test-specific folders. In all tests, we specify in what limits SimSS/ZimT should produce the same results as the analytical expression it is tested against.

# How to run the tests
There is a script to automatically asses the SIMsalabim project that is located in the 'Tests' folder. These tests test the code against physics that apply in certain regimes of operation. For every test case there is a folder that contains the input parameters required for the test to run. In this folder, a plot showing the result of the test will be generated during the test.

A Python 3 interpretor is installed by default on most operating systems, however on Windows it should first be installed (for example from the Microsoft Store).

The test script also requires a few non-default python packages. To install these simply run the command `pip install numpy scipy matplotlib pandas` in the terminal on your computer. 

The test script also requires the free Pascal compiler if not installed already. For details on how to install the free Pascal compiler see README.txt.

To run the automated tests, run `python run_tests.py` in a terminal from the 'Tests' folder.


# SimSS tests

## Test 1: Photoconductivity of an insulator

This test is based on the work by Sokel and Hughes (R. Sokel and R.C. Hughes, J. Appl. Phys. **53**, 7414 (1982)). They provide a formula (their Eq. 49) for the current density in an insulator that is able to absorb light. It assumes that the electric field in the device is constant, i.e. there is a neglible density of electrons and holes. The file JV.dat contains the simulated and analytical result. The agreement is excellent for all voltages.

## Test 2: Space-charge-limited-current with field-dependent mobility
This test assesses the accuracy of the Poisson solver. It simulates a so-called electron-only diode (a diode where the injection of holes is blocked) with an electron mobility that depends on the electric field. The analytical result by Murgatroyd is used as a reference (see  P. N. Murgatroyd, J. Phys. D **3**, 151 (1970)). The analytical result by Murgatroyd neglects the effect of diffusion and, hence, the simulated current density should be lower than the analytical result, especially at low voltages. This is indeed the case: at high voltages, the simulation and analytical result overlap, while at low voltages the simulated current is significantly smaller.

## Test 3: Double-injection with very weak recombination

Rosenberg and Lampert (L.M. Rosenberg and M.A. Lampert, J. Appl. Phys. **41**, 508 (1970)) derive two limiting cases for double injection into an insulator, corresponding to injection of electrons and holes from Ohmic contacts. The semiconductor/insulator does not contain any trapping and recombination is of the direct type. The plasma limit is obtained when the recombination is very weak. The simulated and analytical results are quite close, especially for high voltages. Note, the analytical derivation does not consider any built-in voltage (and no diffusion) which is present in the simulations. As a result, the simulation should be below the analytical result, specifically at low voltages where diffusion and built-in voltage are relevant. 

## Test 4: Space-charge-limited diode with diffusion regime
Test 4 is another test of the space-charge-limited regime. This time, the mobility is taken as a constant (i.e. not dependent on the electric field) and we compare the high- and low-field behaviour with the analytical results as descirbed in J.A. Röhr, T. Kirchartz, and J. Nelson, J. Phys.: Condens. Matter **29**, 205901 (2017).
At high voltages, the simulation should approach the Mott-Gruney law: This is indeed the case, despite the fact that the simulation includes a built-in voltage. At low voltages, the current is linear as a consequence of background charge density that stems from the contacts. This regime is described by the moving-electrode equation. The agreement in both limits is excellent.

## Test 5: Space-charge-limited-diode with trapping
In this test, we assess the ability of SimSS to simulate a diode with a single trap level (situated mid-gap). As a reference, we also show the simulated current-voltage curve without traps. To obtain the reference, either put Bulk_tr to zero in the parameter file, or run SimSS with -Bulk_tr 0. The inclusion of traps, as expected, strongly reduces the current at low voltages. The trap-filled-limit (V<sub>TFL</sub>) is obtained as outlined in V.M. Le Corre, E.A. Duijnstee, O.  El Tambouli, J.M. Ball, H.J. Snaith, J. Lim, and L.J.A. Koster, ACS Energy Lett. **6**,  1087 (2021): two tangents are used (dashed lines in graph), one in the steepest part of the curve and one in the quadratic regime. They cross at 1.2 V, which agrees very nicely with the anticipated value for V<sub>TFL</sub> of 1 V which is based on the input parameters.


# ZimT tests

## Test 6: Open-circuit voltage in steady-state
ZimT can be used to solve for the open-circuit voltage of a solar cell, in transient cases as well as steady-state. The open-circuit in steady-state is compared with an analytical solution given in L.J.A. Koster, V.D. Mihailetchi, R. Ramaker, and P.W.M. Blom, Appl. Phys. Lett. **86**, 123509 (2005). The agreement is excellent.


## Test 7: RC-time
This tests ZimT by simulating an RC-circuit. The parameters define a capacitor with a series resistance. The voltage is changed from 0V (initial condition) to 1V. The RC-time, based on the parameters, is 5 μs. A fit to the simulated data yields an RC time of 5.36 μs, which is satisfactory.


## Test 8: Transient Photo-voltage (TPV)
This test ZimT's ability to simulate the transient decay of the open-circuit voltage (Voc) upon a small change in the light intensity. We start at a generation rate of electron-hole pairs of Gehp = 1.2E26 m<sup>-3</sup>s<sup>-1</sup> at time zero, followed by a quick reduction to Gehp = 1.0E26 m<sup>-3</sup>s<sup>-1</sup>. The Voc will thus decay a little bit. The decay rate can be estimated based on A. Rahimi Chatri, S. Torabi, V.M. Le Corre, and L.J.A. Koster, ACS Appl. Mater. Interfaces **10**, 12013 (2018). Their formula (Eq. 9) predicts a decay time of 7.09 μs for the parameters used in this test. Fitting to the simulated data yield a decay time of 6.96 μs.

## Test 9: Transient Photo-current with bulk traps (TPC)
This test assesses ZimT's ability to simulate the emptying of bulk traps: The simulation starts with a solar cell at 0 V under illumination. Quickly, the light is switched off and free carriers are swept out of the device (the mobilities are very high). Trapped electrons will be slowly released from the traps and this is what generates the current in the ms-regime. The parameters are chosen such that the detrapping time should be 0.050 s. Fitting to ZimT's results indeed yields a decay time of the current of 0.050 s as expected.

## Test 10: Transient Photo-current with interface traps at grain boundaries (TPC)
This test assesses ZimT's ability to simulate the emptying of interface traps that sit at grain-boundaries (GBs). It is similar to Test 9, but not there are 10 GBs with interface traps that trap and de-trap holes. The input parameters are chosen such that de-trapping should have a lifetime of 25 ms. Fitting to ZimT's result indeed yields a lifetime of 0.025 s as expected.

## Test 11: Transient Photo-current for estimating Urbach energy
This tests ZimT's ability to simulate detrapping of charges from bulk traps that are distributed exponentially from the conduction band defined by an Urbach energy of 100 meV. The Urbach energy can be fitted by transforming time and current output as described in MacKenzie, et al., J. Phys. Chem. C **117 (24)**, 12407-12414 (2013). First we run the steady state solver under illumination, effectively filling the traps present. Then we switch off the light and track detrapping charges through the current. We then manipulate the current and time values to estimate the DOS(E) and fit the Urbach energy, which yields 95 meV, sufficiently close to the input value.

## Test 12: Transient Photo-current for estimating Urbach energy
This tests ZimT's ability to simulate detrapping of charges from interface traps that are distributed exponentially from the valence band, defined by an Urbach of 70 meV. We use 10 grain boundaries for this purpose. The Urbach energy can be fitted by transforming time and current output as described in MacKenzie, et al., J. Phys. Chem. C **117 (24)**, 12407-12414 (2013). First we run the steady state solver under illumination, effectively filling the traps present. Then we switch off the light and track detrapping charges through the current. We then manipulate the current and time values to estimate the DOS(E) and fit the Urbach energy, which yields 63 meV, i.e. sufficiently close to the input value.












