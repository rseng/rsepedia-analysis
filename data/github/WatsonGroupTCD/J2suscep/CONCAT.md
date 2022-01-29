[![DOI](https://zenodo.org/badge/292014721.svg)](https://zenodo.org/badge/latestdoi/292014721)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02838/status.svg)](https://doi.org/10.21105/joss.02838)

# J2suscep
This package contains two standalone codes written in FORTRAN 2008 – ej_calc and suscep. 
The ej_calc code uses the spin densities obtained from DFT calculations and determines the isotropic exchange coupling between paramagnetic centres. 
The suscep code calculates the temperature dependence of magnetic susceptibility using the coupling constants.

For usage instructions and detailed documentation of the codes, please visit: https://github.com/WatsonGroupTCD/J2suscep/wiki
<html><head></head><body>


<h1>Contributing</h1>
<h2>Bugs reports and feature requests</h2>

<p>If you think you've found a bug, please report it on the 
<a href="https://github.com/WatsonGroupTCD/J2suscep/issues">Issue Tracker</a>.
You can also propose ideas for new features or ask questions about the design of J2suscep here.
Poor documentation can also be considered a bug, but please be as specific as possible when asking for improvements.</p>


<h2>Code contributions</h2>
<p>We welcome your contributions towards improving and extending the package and this is managed through Github pull requests.
For external contributions we prefer the <a href="https://guides.github.com/activities/forking/">"fork and pull"</a>
workflow while core developers use branches in the main repository:</p>

<blockquote>
<ol>
<li>First open an Issue to discuss the proposed contribution. This discussion might include how the changes fit J2suscep's scope and a
general technical approach.</li>
<li>Make your own project fork and implement the changes there.</li>
<li>Open a pull request to merge the changes into the main project. A more detailed discussion can take place there before
the changes are accepted.</li>
</ol>
</blockquote>
---
title: 'J2suscep: Calculation of magnetic exchange coupling and temperature dependence of magnetic susceptibility'
tags:
  - magnetism
  - ab initio
  - density functional theory
authors:
 - name: Swetanshu Tandon
   orcid: 0000-0002-0626-0626
   affiliation: "1, 2"
 - name: Wolfgang Schmitt
   orcid: 0000-0002-0058-9404
   affiliation: "1, 2"
 - name: Graeme W. Watson
   orcid: 0000-0001-6732-9474
   affiliation: "1"
   email: watsong@tcd.ie
affiliations:
 - name: School of Chemistry & CRANN Institute, University of Dublin, Trinity College, Dublin 2, Ireland
   index: 1
 - name: AMBER Centre, University of Dublin, Trinity College, Dublin 2, Ireland
   index: 2
date: September 2020
bibliography: paper.bib
---

# `Statement of need`

The field of molecular magnetism has attracted significant interest owing 
to the rising need to miniaturize magnets. Advances in this field 
require a deep understanding of the electronic structure of magnets, and computational 
approaches have played a key role in pushing the limits 
further [@Neese]. A deep understanding of the magnetic properties not 
only helps in making magnets smaller, but has also been instrumental for
reasoning about the structure of systems that are difficult to crystallise, and
about intermediate states during catalytic processes [@Krewald].

In molecular complexes, the magnetic interaction between paramagnetic 
centres, called the isotropic exchange interaction, has a strong influence 
on the magnetic properties. Within the framework of density functional 
theory (DFT), using what is known as the broken symmetry approach, this 
interaction can be quantified; the coupling strength is called the 
*J-value* or coupling constant [@Bsym; @Bsym1]. The broken symmetry approach requires the 
modelling of multiple spin states using DFT but the final solution obtained 
 can be dependent on the states modelled [@Cremades1; @Cremades2; 
@Rajeshkumar; @Vignesh; @Mn6]. To remove this dependency, more states 
can be modelled but this procedure results in a large number of solutions, 
which have to be averaged in some sensible way. Additionally, the problem of 
singular solutions also arises, which must be identified and removed. 
Until now, however, there has been no convenient way to do so. 

The purpose of this 
package is to provide a means to accomplish this arduous task. Additionally, 
this package can also calculate the temperature dependence of magnetic 
susceptibility, thereby enabling the comparison of computational data with
experiment. Although other codes like MAGPACK [@magpack] and PHI [@phi]
are available for the calculation of the magnetic susceptibility, 
this package provides a one-stop solution for the calculation of coupling 
constants and the determination of the magnetic susceptibility using these 
coupling constants.


# `J2suscep`

This package essentially contains two standalone codes written in FORTRAN 2008 – 
`ej_calc` and `suscep`. Both codes rely only on the LAPACK library [@Lapack]. The 
code `ej_calc` uses the data obtained from the DFT calculations and determines the 
isotropic exchange coupling between paramagnetic centres. With a completely 
programmable Hamiltonian, this code allows the calculation of any number of 
J-values and is only limited by the number of states modelled. It employs the spin 
density approach [@spin_dens], and spin densities obtained from any approach can be used for the 
calculation of J-values. The code calculates all possible solutions based on the 
Hamiltonian and the states modelled, removes any singular solutions and calculates 
an average set of coupling constants and the standard deviations. Additionally, it 
also calculates the energy of the different spin states relative to each other using 
the coupling constants, for comparison to the original DFT data.

The `suscep` code calculates the temperature dependence of magnetic susceptibility 
using the coupling constants. Similar to `ej_calc`, the Hamiltonian is flexible and 
any number of coupling constants can be provided for the calculation of the magnetic 
susceptibility. 

To illustrate the use of this package, we present the example of the Mn$_6$ complex 
shown in Figure 1 (a) [@Mn6]. The Mn atoms are arranged in an octahedron, resulting in 
each Mn being cis to 4 other Mn atoms and trans to 1 Mn atom with the Mn atoms 
interacting via Cl^-^ and phosphonate bridges. One requires 2 J-values 
to account for these cis- and trans- interactions between Mn centres. Modelling of 6 states 
for the calculation of the 2 J-values results in 60 possible solutions; using `ej_calc`, 
the final values were determined to be -1.28 (cis-coupling) and -3.48 cm^-1^ (trans-coupling) 
respectively [@Mn6]. The use of the `suscep` code calculates the temperature dependence of 
susceptibility which is shown in Figure 1 (b). This package has been used in a similar manner 
for probing the magnetic properties of a Mn$_8$ complex [@Mn8].

![(a) Structure of the Mn$_6$ complex and (b) the temperature dependence of magnetic 
susceptibility obtained using the suscep code. Colour scheme: Mn (dark blue), P (pink), Cl (green), 
C (black), N (blue) and O (red). Hydrogen atoms have been removed for clarity.](plot.tif)


# Acknowledgements

The authors are grateful to the Irish Research Council (GOIPG/2015/2952), Science Foundation 
Ireland (12/IA/1414 and 13/IA/1896) and the European Research Council (CoG 2014–647719) for funding.
The authors also  acknowledge ICHEC for computational resources.  


# References
<html><head></head><body>


<h1>Contributing</h1>
<h2>Bugs reports and feature requests</h2>

<p>If you think you've found a bug, please report it on the 
<a href="https://github.com/WatsonGroupTCD/J2suscep/issues">Issue Tracker</a>.
You can also propose ideas for new features or ask questions about the design of J2suscep here.
Poor documentation can also be considered a bug, but please be as specific as possible when asking for improvements.</p>


<h2>Code contributions</h2>
<p>We welcome your contributions towards improving and extending the package and this is managed through Github pull requests.
For external contributions we prefer the <a href="https://guides.github.com/activities/forking/">"fork and pull"</a>
workflow while core developers use branches in the main repository:</p>

<blockquote>
<ol>
<li>First open an Issue to discuss the proposed contribution. This discussion might include how the changes fit J2suscep's scope and a
general technical approach.</li>
<li>Make your own project fork and implement the changes there.</li>
<li>Open a pull request to merge the changes into the main project. A more detailed discussion can take place there before
the changes are accepted.</li>
</ol>
</blockquote>

 
</body></html>
<html><head></head><body>
<h1>Basic Details</h1>

<p>This code has been written in FORTRAN 2008. It requires one input file that specifies the following:</p>
<p>1.	Number of paramagnetic centres in the system.</p>
<p>2.	The number of unpaired electrons on each centre. </p>
<p>3.	The coupling behaviour between these centres and the coupling strength.</p>
<p>4.	The g-value.</p>
<p>5.	The strength of the applied magnetic field. </p>
<p>6.	The temperature interval and range in which the value of χ and the χT need to determined. </p>
<p>This code generates a .out file that contains the value of χ and the χT product at different temperatures. The output 
also contains details about the size of Hamiltonian matrix and the M<sub>s</sub> values associated with each paramagnetic centre.</p>
<p></p>

    
</body></html>
<html><head></head><body>

<h1>Worked Example</h1>

<p>This section describes how one can obtain the necessary data required to prepare <a href="https://github.com/WatsonGroupTCD/J2suscep/wiki/Execution#format-of-file2"> File2</a> and execute J2suscep using the available computational chemistry code <a href="https://gaussian.com/">Gaussian</a>. The example here corresponds to the system described as example1 under 'Example Input Files' subsection of <A HREF="1e_example">EJ_Calc</A> and <A HREF="2e_example">Suscep</A> in this manual. The procedure can be summarised as follows:</p>

<blockquote>
<ol>
<li>Model the ferromagnetic and broken symmetry states of the molecule under investigation.</li>
<li>Extract the spin density data and the energy of the system from the generated output files.</li>
<li>Use the spin density and energy data to prepare <a href="https://github.com/WatsonGroupTCD/J2suscep/wiki/Execution#format-of-file2"> File2</a> for J2suscep (EJ_calc).</li>
<li>Run J2suscep.</li>
</ol>
</blockquote>

<p><b>Note:</b> It is being assumed that the paths to the utilities and codes described in this example are already added to the PATH variable.</p>

<h2>Modelling Ferromagnetic and Broken Symmetry States</h2>


<p>The determination of coupling constants generally relies on the use of broken symmetry states which essentially provide a means to represent 
antiferromagnetic states using a single determinant. Most computational chemistry codes are capable of modelling such states. </p>

<p>The procedure to model broken symmetry states using the <a href="https://gaussian.com/">Gaussian</a> software is the following:</p>

<blockquote>
<ol>
<li>Obtain the crystallographic coordinates or generate a model for the system to be investigated.</li>
<li>Model the ferromagnetic state and optimise the geometry.</li>
<li>Use the ferromagnetic state to generate a guess for the antiferromagnetic state.</li>
<li>Analyse the wavefunction to ensure that the wavefunction is stable. In case of instability, optimise the wavefunction until it is stable.</li>
<li>Optimise the geometry.</li>
</ol>
</blockquote>

<p>Once the ferromagnetic and antiferromagnetic states are optimised, the wavefunction can be used to extract the spin density data. 
Mulliken charge analysis is carried out in <a href="https://gaussian.com/">Gaussian</a> by default. One can also perform the 
<a href="https://gaussian.com/population/">Hirshfeld</a> charge analysis and 
<a href="http://theory.cm.utexas.edu/henkelman/code/bader/"> Bader</a> charge analysis to determine the spin densities.</p> We personally recommend the use of Bader spin densities.
<p>(<b>Note:</b> Sometimes, the coupling constants are also determined using only the crystallographic coordinates without further optimisation. 
In such instances, only the wavefunction is optimised to model the ferromagnetic and antiferromagnetic states. The geometry is not optimised.)</p>

<p>To explain the procedure, we take the example of a dimeric Mn complex shown in Figure 1 below
<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs.htm#Libby" class="showTip Libby">Libby</a>]</font></p></p>
<p align="center"><img src="img/mn2.jpg" alt="Mn2 structure"><br /><b><i>Structure of the {Mn<sup>IV</sup><sub>2</sub>} complex. Colour scheme: Mn (dark blue), C (black), 
N (blue) and O (red). Hydrogen atoms have been removed for clarity.</i></b></p> 
<br>

<h3>Ferromagnetic state</h3>
<p>The coordinates for this complex were obtained from the CCDC database (CCDC reference: JAVYEK). 
<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs.htm#CCDC" class="showTip CCDC">CCDC</a>] </font>
The system was modelled using the PBE0 functional in conjunction with the SDDALL basis set having an effective core potential 
for the Mn atoms and 6-31g (d,p) basis set for C, O, N and H atoms. The system was divided into 3 fragments such that the Mn centres 
belong to two separate fragments and the ligand environment was included in a separate fragment. The input file for modelling the 
ferromagnetic state looks like the following:</p>
<pre>
%nprocshared=16
%mem=50000MB
%chk=ferromagnetic.chk
# opt upbe1pbe/gen geom=connectivity pseudo=read

JAVYEK

0 7 -8 1 4 4 4 4
 Mn(Fragment=2)     9.68850000    3.22950000    5.30440000
 Mn(Fragment=3)     9.68850000    3.22950000    2.55860000
 O(Fragment=1)      8.84500000    4.07300000    3.93150000
....
....
....

C O N H 0
6-31g(d,p)
_****_
Mn 0
SDDALL
_****_

Mn 0
SDDALL
</pre>
<br>


<p>The first two lines here specify the number of processors (16) and the memory (50 GB) dedicated for running the calculation. 
The third line specifies the name of the checkpoint file which stores the information about the wavefunction. 
The fourth line specifies the calculation specific terms which have the following meaning:</p>


<table class="tg">
<thead>
  <tr>
    <th class="tg-0lax">opt</th>
    <th class="tg-0lax">Optimise the geometric coordinates.</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0lax">upbe1pbe </td>
    <td class="tg-0lax">Use the PBE0 functional.</td>
  </tr>
  <tr>
    <td class="tg-0lax">gen</td>
    <td class="tg-0lax">Use the basis set described at the bottom of the input file.</td>
  </tr>
  <tr>
    <td class="tg-0lax">geom=connectivity </td>
    <td class="tg-0lax">Specify explicit atom bonding data via an additional input section.</td>
  </tr>
  <tr>
    <td class="tg-0lax">pseudo=read </td>
    <td class="tg-0lax">Use the effective core potential described at the bottom of the input file.</td>
  </tr>
</tbody>
</table>

<p>The next line specifies the title. The line below that specifies the charge and multiplicity of the overall system 
and the various fragments. This is followed by the description of the atomic coordinates, the connectivity information 
and the basis set information.</p>
<p>The complete <a href="https://github.com/WatsonGroupTCD/J2suscep/blob/master/examples/gaussian_files/ferromagnetic.com">input</a> file is provided as part of the repository.</p>
<br>
<p>To run the calculation using gaussian09, the following command can be used (assuming that the input file name is ferromagnetic.com):</p>
<pre>g09 < ferromagnetic.com > ferromagnetic.out</pre>
<p>This will generate a file called 'ferromagnetic.out' which will contain information in human readable format. Another file called 'ferromagnetic.chk' containing structure and wavefunction information in binary form will also be generated. </p>



<h3>Antiferromagnetic state</h3>

<p>The optimised ferromagnetic state is used generate a guess for the antiferromagnetic state. The <a href="https://github.com/WatsonGroupTCD/J2suscep/blob/master/examples/gaussian_files/antifer_guess.com">input</a> file for 
calculating a guess looks like the following:</p>

<pre>
%nprocshared=16
%mem=50000MB 
%chk=JAVYEK_antiferro.chk
# upbe1pbe/gen guess=(only,fragment=3) geom=connectivity pop=minimal pseudo=read

JAVYEK

0 1 -8 1 4 4 4 -4
 Mn(Fragment=2)     0.00000900    1.35724000   -0.00000800
 Mn(Fragment=3)    -0.00000900   -1.35724000   -0.00000800
 O(Fragment=1)     -1.19325800    0.00000800   -0.00001300
……
……
……
C O N H 0
6-31g(d,p)
****
Mn 0
SDDALL
****

Mn 0
SDDALL
</pre>
<br>

<p>The term ‘guess=(only,fragment=3)’ on the second line sets up the guess calculation with three fragments 
(which is the total number of fragments specified in this system). The term ‘pop=minimal’ controls the amount of data 
written in the log file and is not mandatory for the proper execution of this calculation. The execution (assuming the input file name to be 'antifer_guess.com') is performed in the same manner as for the ferromagnetic state:</p>
<pre>g09 < antifer_guess.com > antifer_guess.out</pre>


<p>Once the guess wavefunction is generated, it is analysed for any instability using the following <a href="https://github.com/WatsonGroupTCD/J2suscep/blob/master/examples/gaussian_files/antifer_stab.com">input</a> file:</p>
<pre>
%nprocshared=16
%mem=50000MB
%chk=antiferro.chk
# stable=opt upbe1pbe/gen geom=allcheck pseudo=read guess=read

N C O H 0
6-31g(d,p)
_****_
Mn 0
SDDALL
_****_

Mn 0
SDDALL
</pre>
<br>

<p>This calculation reads information from the checkpoint file antiferro.chk generated by the guess calculation. 
The term ‘stable=opt’ on the fourth line specifies that this calculation checks for the stability of the wavefunction and 
optimise it if any instability is found. The calculation updates the information in the checkpoint file antiferro.chk.</p>
<p>The execution (assuming the input file name to be 'antifer_stab.com') is performed as follows:</p>
<pre>g09 < antifer_stab.com > antifer_stab.out</pre>
<p>Once the wavefunction is optimised, it can be used to optimise the geometry using the following <a href="https://github.com/WatsonGroupTCD/J2suscep/blob/master/examples/gaussian_files/antifer_opt.com">input</a> file:</p>

<pre>
%nprocshared=16
%mem=50000MB
%chk=antiferro.chk
# opt upbe1pbe/gen geom=allcheck pseudo=read guess=read

N C O H 0
6-31g(d,p)
****
Mn 0
SDDALL
****

Mn 0
SDDALL
</pre>
<br>

<p>Again, the checkpoint file from the previous calculation is being read here for geometry optimisation.</p>
<p>The execution (assuming the input file name to be 'antifer_opt.com') is performed as follows:</p>
<pre>g09 < antifer_opt.com > antifer_opt.out</pre>
<p><b>Note:</b> Whenever running these calculations, make sure that the spin densities on the metal centres is close to what is expected. This can be done by looking up the Mullken spin densities which are printed by default in gaussian.  Unexpected spin values may indicate problems in the wavefunction optimisation.</p>



<h2>Extracting spin density and energy data</h2>


<p>Once the calculations are completed, one can start extracting the data from the output files and arrange the data in the format required for <a href="https://github.com/WatsonGroupTCD/J2suscep/wiki/Execution#format-of-file2"> File2</a>. <b> Make sure that the spin density and energy data being extracted correspond to that after the system has been optimised.</b></p>

<p>The energy information can be easily extracted from the .out/.log files. For spin denisties, if one is interested in only the Mulliken densities, they may be copied from the .out/.log files and copied to File2 in the proper format. For <a href="https://gaussian.com/population/">Hirshfeld</a> charge analysis, one needs to include 'pop=Hirshfeld' keyword in the optimisation runs for the ferromagnetic and the broken symmetry states.</p>
<p>Gaussian09 does not perform Bader analysis but it can be easily done using a code provided by the <a href="http://theory.cm.utexas.edu/henkelman/code/bader/"> Henkelman group</a>. To perform the Bader analysis, one needs to make use of the .chk files to create the formatted checkpoint (.fchk) files. For the above example, once the calculations for ferromagnetic and broken symmetry states are complete, one will have two .chk files - ferromagnetic.chk and antiferro.chk. To create the .fchk files, the <a href="https://gaussian.com/formchk/">formchk</a> utility can be used as:</p>
<pre>
formchk ferromagnetic.chk
formchk antiferro.chk
</pre>


<p>Once the formatted checkpoint files are obtained, they are used to generate the 'cube' files using the <a href="https://gaussian.com/cubegen/">cubegen</a> utility in gaussian. For both the .fchk files, two cube files - containing spin density and total electron density - will be generated as:</p>
<pre>
cubegen 2 spin=scf ferromagnetic.fchk ferromagnetic_spin.cube -4 
cubegen 2 density=scf ferromagnetic.fchk ferromagnetic_tot_dens_fine.cube -4
cubegen 2 spin=scf antiferro.fchk antiferro_spin.cube -4
cubegen 2 density=scf antiferro.fchk antiferro_tot_dens_fine.cube -4
</pre>

<p>The 'ferromagnetic_spin.cube' and 'antiferro_spin.cube' files contain information about the spin densities while the files 'ferromagnetic_tot_dens_fine.cube' and 'antiferro_tot_dens_fine.cube' contain information about the total electron densities. </p>
<p>These '.cube' files are then used to perform Bader analysis as follow:</p>
<pre>
bader ferromagnetic_spin_fine.cube -ref ferromagnetic_tot_dens_fine.cube
bader antiferro_spin_fine.cube -ref antiferro_tot_dens_fine.cube
</pre>

<p>This will generate a file called 'ACF.dat' containing the Bader spin densities for the different states. It must be kept in mind that the output file is always named ACF.dat with the bader code. Hence, it is necessary that the execution of the bader code for the various spin states in performed in a way to avoid overwriting of the output files.</p>


<p>Once the energy and spin density information is obtained, one can create <a href="https://github.com/WatsonGroupTCD/J2suscep/wiki/Execution#format-of-file2"> File2</a> which should resemble the input file <a href="https://github.com/WatsonGroupTCD/J2suscep/blob/master/examples/ej_calc_form/example1/spin"> spin</a> provided as part of the repository. This, along with <a href="https://github.com/WatsonGroupTCD/J2suscep/wiki/Execution#format-of-file1"> File1</a>, <a href="https://github.com/WatsonGroupTCD/J2suscep/wiki/Execution#format-of-file3"> File3</a> and <a href="https://github.com/WatsonGroupTCD/J2suscep/wiki/Execution#format-of-file4"> File4</a>, can be then used to execute J2suscep. The information required for <a href="https://github.com/WatsonGroupTCD/J2suscep/wiki/Execution#format-of-file1"> File1</a> and <a href="https://github.com/WatsonGroupTCD/J2suscep/wiki/Execution#format-of-file3"> File3</a> primarily depends upon the system under consideration and require no calculations as such. For this example, the details of these files are provided under the 'Example Input Files' subsection of <A HREF="1e_example">EJ_Calc</A> and <A HREF="2e_example">Suscep</A>. The information for <a href="https://github.com/WatsonGroupTCD/J2suscep/wiki/Execution#format-of-file4"> File4</a> is obtained after executing ej_calc and the details are provided in the <A HREF="1e_example">Example Input Files</A> section of Suscep.</p>


</body></html>
<html><head></head><body>
<h1>Example Input Files</h1>

<h2>Example 1</h2>
<p>To understand the construction of the input files, let us take the example of a dimeric Mn complex shown in Figure 1 below.
Here the Mn centres are in +IV oxidation state (d<sup>3</sup> spin configuration).
<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs.htm#Libby" class="showTip Libby">Libby</a>]</font></p>

<p align="center"><img src="img/mn2.jpg" alt="Mn2 structure"><br /><b><i>Structure of the {Mn<sup>IV</sup><sub>2</sub>} complex. Colour scheme: Mn (dark blue), C (black), 
N (blue) and O (red). Hydrogen atoms have been removed for clarity.</i></b></p> 
<p>The Mn centres in this complex interact with each other via the oxo groups and the Hamiltonian for this complex can be written as follows:</p>
<center> <font face="Courier New, Courier, monospace" size="-2">&nbsp;&nbsp;&nbsp;&nbsp;</font>H&#770 = -2J&lt;s<sub>1</sub>.s<sub>2</sub>&gt; </center>

<p>The input files using the Hamiltonian mentioned above can be prepared as described below. The following additional parameters are used to complete the file
 – J1: -84.0722 cm<sup>-1</sup>, g-value: 2.02, field strength: 1000 oersted, temperature range: 2-300 K, temperature step: 1 K. </p>
 <h3> File1</h3>
 <pre>
magnetic centres
2
No. of J values
1
Hamiltonian
1 2
Hamiltonian Ends
 </pre>
<br>

<p>Here the number 2 in the second line specifies that there are 2 magnetic centres. This is followed by the number of J-values and  
the details about the Hamiltonian. 

<h3> File3 (additional parameters)</h3>
<pre>
spin
1.5
1.5
g value
2.02
Field Strength
1000.0
Temperature range
2
300
Step size
1
</pre>
<br>

<p>This file starts with the specification of 
the spin on each centre. In this example, each centre has 3 unpaired electrons (Mn(IV) centres) and hence a spin value of 1.5. 
This is followed by the specification of the g-value 
which is then followed by the magnitude of the field strength in Oersted units (1000.0 in this case). The susceptibility will
be determined in the temperature range of 2-300 K with 1 K increments.</p>
<p>This example input file and the resultant output file have also been provided as separate files in the package. </p>

<h3> File4 (J-values)</h3>
<pre>
J values
-84.0722
</pre>
<br>
<p>In this file the details about the J-values are provided. In this example, 1 J-value is required and the magnitude of the J-value (-84.0722 cm<sup>-1</sup>) follows the term "J values".</p>



<h2>Example 2</h2>
This examples corresponds to a {Mn<sub>3</sub>} complex shown in Figure 2.
<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs.htm#Li" class="showTip Li">Li</a>]</font>  </p>

<p align="center"><img src="img/mn3.jpg" alt="Mn3 structure"><br /><b><i>Structure of the {Mn<sup>III</sup><sub>3</sub>} complex. Colour scheme: 
Mn (dark blue), C (black), N (blue) and O (red). Hydrogen atoms have been removed for clarity.</i></b></p> 

<p>The structure of this molecule suggests that the interaction between each pair of Mn centres would be similar and hence, the Hamiltonian 
for this system can be written as follows:</p>
<center> <font face="Courier New, Courier, monospace" size="-2">&nbsp;&nbsp;&nbsp;&nbsp;</font>H&#770 = -2J[&lt;s<sub>1</sub>.s<sub>2</sub>&gt; + 
&lt;s<sub>1</sub>.s<sub>3</sub>&gt; + &lt;s<sub>2</sub>.s<sub>3</sub>&gt;] </center>
<p>The input files using the Hamiltonian mentioned above can be prepared as described below. The following additional parameters are used to complete the file
 – J1: -9.0221 cm<sup>-1</sup>, g-value: 1.996, field strength: 1000 oersted, temperature range: 2-300 K, temperature step: 1 K. </p>

<h3> File1</h3>
<pre>
magnetic centres
3
No. of J values
1
Hamiltonian
1 2
1 3
2 3
Hamiltonian Ends
</pre>
<br>


<p>Here the number 3 in the second line specifies that there are 3 magnetic centres. This is again followed by the number of J-values (1) and  
the details about the Hamiltonian. 


<h3> File3 (additional parameters)</h3>
<pre>
spin
2.0
2.0
2.0
g value
1.996
Field Strength
1000.0
Temperature range
2
300
Step size
1
</pre>
<br>

<p>This file starts with the specification of 
the spin on each centre. In this example, each centre has 4 unpaired electrons (Mn(III) centres) and hence a spin value of 2. 
This is followed by the specification of the g-value 
which is then followed by the magnitude of the field strength in Oersted units (1000.0 in this case). The susceptibility will
be determined in the temperature range of 2-300 K with 1 K increments.</p>
<p>This example input file and the resultant output file have also been provided as separate files in the package. </p>

<h3> File4 (J-values)</h3>
<pre>
J values
-9.0221
</pre>
<br>

<p>In this file the details about the J-values are provided. In this example, 1 J-value is required and the magnitude of the J-value (-9.0221 cm<sup>-1</sup>) follows the term "J values".</p>


<h2>Example 3, 4 and 5</h2>
<p>Similar to the examples described for ej_calc code, these examples correspond to a {Mn<sub>6</sub>} coordination complex shown in the figure below.
<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs.htm#Tandon" class="showTip Tandon">Tandon</a>]</font>  </p>

<p align="center"><img src="img/mn6.jpg" alt="Mn6 structure"><br /><b><i>Structure of the {Mn<sub>6</sub>} complex. Colour scheme: Mn (dark blue), P (pink), Cl (green), C (black), 
N (blue) and O (red). Hydrogen atoms have been removed for clarity.</i></b></p>
 
<p>The Mn centres in this complex interact with each other via the phosphonate ligands and the central Cl<sup>-</sup> ion. 
Each Mn centre in this complex has 1 trans- and 4 cis- neighbours and one requires 2 J-values to account for these cis- and trans- interactions between Mn centres.
<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs.htm#Tandon" class="showTip Tandon">Tandon</a>]</font>   </p>

<p>In example 3, four of the Mn centres, Mn1, Mn2, Mn4 and Mn5, are in +IV oxidation state (d<sup>3</sup> spin configuration) while the others are in +III oxidation state. 
This decreases the overall symmetry of the complex. Therefore, the Hamiltonian for this complex can be written as follows:</p>
<center> <font face="Courier New, Courier, monospace" size="-2">&nbsp;&nbsp;&nbsp;&nbsp;</font>H&#770 = -2J<sub>1</sub>[&lt;s<sub>1</sub>.s<sub>2</sub>&gt;]
-2J<sub>2</sub>[&lt;s<sub>1</sub>.s<sub>6</sub>&gt;] -2J<sub>3</sub> [&lt;s<sub>2</sub>.s<sub>6</sub>&gt;]	
-2J<sub>4</sub>[&lt;s<sub>3</sub>.s<sub>6</sub>&gt;] -2J<sub>5</sub>[&lt;s<sub>4</sub>.s<sub>5</sub>&gt;] -2J<sub>6</sub>[&lt;s<sub>1</sub>.s<sub>3</sub>&gt; + 
&lt;s<sub>2</sub>.s<sub>3</sub>&gt;] -2J<sub>7</sub>[&lt;s<sub>1</sub>.s<sub>4</sub>&gt; + 	&lt;s<sub>2</sub>.s<sub>5</sub>&gt;]
-2J<sub>8</sub>[&lt;s<sub>1</sub>.s<sub>5</sub>&gt; + &lt;s<sub>2</sub>.s<sub>4</sub>&gt;] -2J<sub>9</sub>[&lt;s<sub>3</sub>.s<sub>4</sub>&gt;
 + &lt;s<sub>3</sub>.s<sub>5</sub>&gt;] -2J<sub>10</sub>[&lt;s<sub>4</sub>.s<sub>6</sub>&gt; + &lt;s<sub>5</sub>.s<sub>6</sub>&gt;] </center>
<br>

<p>In example 4, the Mn centres Mn1, Mn2 and Mn4 are in +IV oxidation state while the others are in +III oxidation state. The overall symmetry 
of the complex is significantly reduced and capturing the full electronic picture requires the use of the following Hamiltonian:</p>
<center> <font face="Courier New, Courier, monospace" size="-2">&nbsp;&nbsp;&nbsp;&nbsp;</font>H&#770 = -2J<sub>1</sub>[&lt;s<sub>1</sub>.s<sub>2</sub>&gt;]
-2J<sub>2</sub>[&lt;s<sub>1</sub>.s<sub>3</sub>&gt;] -2J<sub>3</sub>[&lt;s<sub>1</sub>.s<sub>4</sub>&gt;] -2J<sub>4</sub> [&lt;s<sub>1</sub>.s<sub>5</sub>&gt;]
-2J<sub>5</sub> [&lt;s<sub>1</sub>.s<sub>6</sub>&gt;] -2J<sub>6</sub> [&lt;s<sub>2</sub>.s<sub>3</sub>&gt;] -2J<sub>7</sub> [&lt;s<sub>2</sub>.s<sub>4</sub>&gt;]
-2J<sub>8</sub> [&lt;s<sub>2</sub>.s<sub>5</sub>&gt;] -2J<sub>9</sub> [&lt;s<sub>2</sub>.s<sub>6</sub>&gt;] -2J<sub>10</sub> [&lt;s<sub>3</sub>.s<sub>4</sub>&gt;]
-2J<sub>11</sub>[&lt;s<sub>3</sub>.s<sub>5</sub>&gt;] -2J<sub>12</sub> [&lt;s<sub>3</sub>.s<sub>6</sub>&gt;] -2J<sub>13</sub> [&lt;s<sub>4</sub>.s<sub>5</sub>&gt;] 
-2J<sub>14</sub> [&lt;s<sub>4</sub>.s<sub>6</sub>&gt;] -2J<sub>15</sub> [&lt;s<sub>5</sub>.s<sub>6</sub>&gt;]</center>
<br>


	
<p> In example 5, all the Mn centres are in +III oxidation state. Only 2 J-values are required to account the exchange interactions between Mn centres
and the Hamiltonian for this complex can thus be written as follows:</p>
<center> <font face="Courier New, Courier, monospace" size="-2">&nbsp;&nbsp;&nbsp;&nbsp;</font>H&#770 = -2J<sub>1</sub> [&lt;s<sub>1</sub>.s<sub>2</sub>&gt; + 	
&lt;s<sub>1</sub>.s<sub>3</sub>&gt; + &lt;s<sub>1</sub>.s<sub>5</sub>&gt; + &lt;s<sub>1</sub>.s<sub>6</sub>&gt; + &lt;s<sub>2</sub>.s<sub>3</sub>&gt; + 
&lt;s<sub>2</sub>.s<sub>4</sub>&gt; + &lt;s<sub>2</sub>.s<sub>6</sub>&gt; + &lt;s<sub>3</sub>.s<sub>4</sub>&gt; + &lt;s<sub>3</sub>.s<sub>5</sub>&gt;
 + &lt;s<sub>4</sub>.s<sub>5</sub>&gt; + &lt;s<sub>4</sub>.s<sub>6</sub>&gt; + &lt;s<sub>5</sub>.s<sub>6</sub>&gt;] -2J<sub>2</sub> [&lt;s<sub>1</sub>.s<sub>4</sub>&gt; + 	
&lt;s<sub>2</sub>.s<sub>5</sub>&gt; + &lt;s<sub>3</sub>.s<sub>6</sub>&gt;]</center>
<br>

The input files containing the Hamiltonian and the spin density data and the resultant output files for examples 3, 4 and 5 are provided in the package.



<p></p>




    
</body></html>
<html><head></head><body>
<h1>Introduction</h1>
<p>Briefly, this code calculates the temperature dependence of magnetic susceptibility for a given system for which the exchange coupling constants have been determined. 
The following Hamiltonian is used to fulfil this purpose:<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs#Boca" class="showTip Boca">Boca</a>]</font>  </p> 
<br>

<p align="center"><img src="img/one.jpg" alt="eq. 1"><br /></p>

<p>where <i>J<sub>AB</sub></i> is the coupling constant between the magnetic centres A and B, <i>I</i> and <i>J</i> represent the particular M<sub>s</sub> terms in the basis elements that the 
<i>S&#770<sub>A</sub>.S&#770<sub>B</sub></i> term operates on <i>B&#8407</i> represents the applied magnetic field, g<sub>A</sub> is the gyromagnetic tensor and μ<sub>B</sub> is the unit 
Bohr Magneton. The first term accounts for the exchange coupling between Mn centres and the second term is the Zeeman term.</p>

<p> The basic working of this code involves setting up a matrix of all possible M<sub>S</sub> (the magnetic spin quantum number) states for the given system. 
This is followed by operating each term of this matrix with eq. 1. Once this operation is complete, the matrix is diagonalised to obtain the eigenvalues which 
are then used to determine the magnetic susceptibility at different temperatures using the van Vleck equation. The details of the implementation of this 
procedure are later. </p>

<table border="0" cellpadding="2">
<tr valign="top"><td><p>J<sub>ij</sub></p></td><td>&nbsp;&nbsp;</td>
<td><p>The exchange coupling constant</p></td></tr>
<tr valign="top"><td><p>s<sub>x</sub></p></td><td>&nbsp;&nbsp;</td>
<td><p>The spin operator.</p></td></tr>
</table>
<p> Any definition of spin operators (formal spins or spin projections or spin densities or any other) can be used. Sometimes researchers exclude the pre-factor 
2 in the above formula so this needs to be kept in mind while using this code.</p>

    
</body></html>
<html><head></head><body>
<h1>For Developers</h1>

<p>The key global variables in this code are defined below:</p>

<table class="tg">
<thead>
  <tr>
    <th class="tg-amwm">  Variable </th>
    <th class="tg-amwm">  Function </th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0lax"> dimension </td>
    <td class="tg-0lax"> stores the number of magnetic centres. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> no_of_j_val </td>
    <td class="tg-0lax"> stores the total number of unique J-values defined for the Hamiltonian. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> totalspin </td>
    <td class="tg-0lax"> stores the row (or column) size of the spin hamiltonian matrix. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> spin_mat_col </td>
    <td class="tg-0lax"> stores the column size of the spin matrix. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> j_mat_col </td>
    <td class="tg-0lax"> stores the column size of the J value matrix. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> spin </td>
    <td class="tg-0lax"> stores the S value (i.e. total number of unpaired e<sup>-</sup>/2) for each metal centre. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> jval </td>
    <td class="tg-0lax"> stores the different J-values in use. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> jmatx, jmaty </td>
    <td class="tg-0lax"> store information about which J-values are associated with different spin pairs. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> hamil </td>
    <td class="tg-0lax"> stores the spin hamiltonian matrix. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> g </td>
    <td class="tg-0lax"> stores the (isotropic) g-value </td>
  </tr>
  <tr>
    <td class="tg-0lax"> B </td>
    <td class="tg-0lax"> stores the magnetic induction in Wb m<sup>-2</sup>. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> vac_perm </td>
    <td class="tg-0lax"> stores the strength of the applied field. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> field_str </td>
    <td class="tg-0lax"> permeability of free space in vacuum. </td>
  </tr>
  <tr>
  <tr>
    <td class="tg-0lax"> init_temp </td>
    <td class="tg-0lax"> stores the minimum temperature for the susceptibility measurement. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> fin_temp </td>
    <td class="tg-0lax"> stores the maximum temperature for the susceptibility measurement. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> step_size </td>
    <td class="tg-0lax"> stores the step size for temperature increments. </td>
  </tr>
    <td class="tg-0lax"> spinmat </td>
    <td class="tg-0lax"> stores all possible M<sub>s</sub> values for each metal centre. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> basis </td>
    <td class="tg-0lax"> stores all possible combinations of M<sub>s</sub> values of all metal centres. </td>
  </tr>
</tbody>
</table>
<br>

<p><b>1.	Init:</b> This subroutine reads in the input file and initialises the global variables.</p>
<p><b>2.	Printer:</b> This subroutine prints some basic details to the output file.</p>
<p><b>3.	SpinMatForm:</b> This subroutine determines the M<sub>s</sub> values for each paramagnetic centre.</p>
<p><b>4.	FormBasis:</b> The subroutine calculates the basis that will be used to span the Hamiltonian matrix. Each basis element comprises of a 
M<sub>s</sub> value for each paramagnetic centre.</p>
<p><b>5.	CheckPos and UpdatePos:</b> These subroutines are helper subroutines for the subroutine ‘FormBasis’ to cycle through all the possible 
combinations of M<sub>s</sub> values for the paramagnetic centres.</p>
<p><b>6.	HamilForm:</b> This subroutine defines the Hamiltonian matrix using the various basis elements and operates each element with the exchange operator.</p>
<p><b>7.	Kron_del:</b> This function is used to determine the value of the kronecker delta function.</p>
<p><b>8.	Diag:</b> This helper subroutine is used for matrix diagonalisation using the LAPACK library.</p>
<p><b>9.	Zeeman:</b> This subroutine adds the Zeeman term to the elements of the Hamiltonian matrix.</p>
<p><b>10.	Suscep_calc:</b> This subroutine defines the B matrix (eq. <a href="2b_theory#twenty_one" class="showTip twenty_one">21</a>), 
calculates the matrix of coefficients (C matrix, eq. <a href="2b_theory#twenty_two" class="showTip twenty_two">22</a>) 
and then calculates the temperature dependence of susceptibility.</p>

<p></p>


    
</body></html>
<html><head></head><body>
<h1>Description of Output File</h1>


<p>The format of the output file is as follows:</p>



<table class="tg">
<thead>
  <tr>
    <th class="tg-amwm">Output file sections</th>
    <th class="tg-amwm">Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-73oq">Number of Magnetic centres:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;6<br>Total number of possible interactions:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;15<br>Input obtained from spin density file<br>..........<br>..........<br>Number of spin density sets:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;6<br>Spin Matrix<br>..........<br>...........<br><br></td>
    <td class="tg-73oq">Details of the basic input parameters and <br>the information obtained from the file <br>containing the spin density information <br>are given here.</td>
  </tr>
  <tr>
    <td class="tg-73oq">No. of J value asked for = 2 <br>Interactions considered under each J value <br>J 1 &nbsp;&nbsp;1 2 &nbsp;&nbsp;1 3 &nbsp;&nbsp;1 5 &nbsp;&nbsp;1 6 &nbsp;&nbsp;2 3 &nbsp;&nbsp;2 4&nbsp;&nbsp;&nbsp;2 6&nbsp;&nbsp;&nbsp;3 4 <br>3 5&nbsp;&nbsp;&nbsp;4 5&nbsp;&nbsp;&nbsp;4 6&nbsp;&nbsp;&nbsp;5 6  <br>J 2&nbsp;&nbsp;&nbsp;1 4&nbsp;&nbsp;&nbsp;3 6&nbsp;&nbsp;&nbsp;2 5 <br>Coefficients of J values<br>.........<br>.........<br>.........<br><br></td>
    <td class="tg-73oq">The interactions accounted for by each <br>J-value are provided here. The <br>coefficients associated with each J-value <br>based on the Hamiltonian are also given.</td>
  </tr>
  <tr>
    <td class="tg-73oq">Possible combination of equations per set of equations <br>considered (i.e. total number of equations -1): <br>5 C&nbsp;&nbsp;&nbsp;2 =     10<br>Solving the equations <br>Average J values<br>Ref. eq.&nbsp;&nbsp;&nbsp; non-singular solutions&nbsp;&nbsp;&nbsp; J 1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; J 2 
	<br>1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-3.4459&nbsp;&nbsp;&nbsp;-1.2797 
	<br>2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-3.4757&nbsp;&nbsp;&nbsp;-1.2704<br>.........<br>.........<br>.........<br><br></td>
    <td class="tg-73oq">The total possible combinations after using <br>one of the equations as a reference is <br>mentioned here. This is followed by the details<br>of the average J-values obtained for each <br>set with each set differing in the equation <br>used as a reference.</td>
  </tr>
  <tr>
    <td class="tg-73oq">Removing solutions in which J-values deviate by more <br>than 3 standard deviations <br>Standard deviation (percentage standard deviation) <br>on different J-values:<br>non-singular equations&nbsp;&nbsp;&nbsp;J 1&nbsp;&nbsp;&nbsp;J 2<br>.........<br>........<br><br></td>
    <td class="tg-73oq">Information on standard deviations for each <br>J-values when cycling through all valid <br>solutions to ensure that each of them is <br>within 3 standard deviations is given here.</td>
  </tr>
  <tr>
    <td class="tg-73oq">Final results <br>Total non-singular equations&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;60 <br>Average J-values (cm-1) and standard deviations <br>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;J 1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;J 2 
	<br>J-val&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.4803&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.2787 <br>Std. Dev.     0.05( 1%)          0.04( 3%)<br><br></td>
    <td class="tg-73oq">Summary of results.</td>
  </tr>
  <tr>
    <td class="tg-73oq">Comparison of energy (cm-1) of electronic states <br>calculated by DFT <br>and energy obtained using the average J-values <br>calculated:<br>(Note: The following are the energies assuming the first <br>energy value in the spin density file as the reference value) <br>DFT Energy(cm-1)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Calculated energy(cm-1)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Absolute Difference(cm-1) <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.00000&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.00000&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.00000( NaN%)<br>........<br>........<br>........<br></td>
    <td class="tg-73oq">A comparison between the energy of the different <br>states calculated using DFT and that calculated <br>using the J-values is provided in this section. <br>This can be used to check how well the coupling <br>constants describe the different states. <br>Large differences between the two energy <br>values indicate problems in either the modelled <br>states or the Hamiltonian used.</td>
  </tr>
</tbody>
</table>




 

<p></p>




    
</body></html>
<html><head></head><body>
<h1>Theory</h1>

<p></p>

<p>For obtaining n coupling constants, at least (n+1) spin configurations need to be modelled. 
In principle, for x centres with unpaired electrons, there can be a maximum of x(x-1)/2 coupling constants and a total of 2<sup>x</sup> different spin configurations. 
Out of these 2<sup>x</sup> configurations, half would be the mirror image of the others. 
For example, for any Mn<sub>6</sub><sup>III</sup> system one can have a maximum of (6 X (6-1) / 2 =) 15 coupling constants while the total number of spin 
configurations for this complex are (2<sup>6</sup> =) 64 out of which a maximum of only 32 can be unique. </p>

<p>For any system with more than 3 paramagnetic centres, the number of configurations that are available are more than 
that required to calculate all coupling constants. This can be used to our advantage since we can model more states than required. 
This allows us to remove any dependence of the calculated coupling constants on the states modelled.<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs#Cremades" class="showTip Cremades">Cremades</a>], [<a href="refs#Cremades1" class="showTip Cremades">Cremades1</a>], [<a href="refs#Rajeshkumar" class="showTip Rajeshkumar">Rajeshkumar</a>], [<a href="refs#Vignesh" class="showTip Vignesh">Vignesh</a>], [<a href="refs#Tandon" class="showTip Tandon">Tandon</a>]</font>  
However, depending upon the size of the system and the number of states modelled, one can end up with a large number of sets of 
solutions (coupling constants) which have to be averaged in some sensible way. 
This problem has been addressed in this code.</p>


<p>When calculating r coupling constants using n (where n > r) spin configurations, then first of all, one of these configurations has 
to be used as a reference for all the other states so that we can determine the energy only related to the flipping of spin. 
With the remaining (n-1) configurations, a total of [(n-1)!⁄{r!(n-1-r)!}] (i.e. <sup>(n-1)</sup>C<sub>r</sub>) sets of solutions (i.e. r coupling constants) are obtained. 
Since any of the n modelled spin states can be used as a reference, n[(n-1)!⁄{r!(n-1-r)!}] (i.e. n(<sup>(n-1)</sup>C<sub>r</sub>) sets of solutions are obtained.</p>

<p>This code takes the Hamiltonian, the values of spin operators for each paramagnetic centre in each modelled spin state and the energy 
of each modelled spin state as input. This information is then used to calculate the coupling constants. 
As mentioned earlier, when calculating r coupling constants using n (where n > r) spin configurations, a total of n(<sup>(n-1)</sup>C<sub>r</sub>) sets of 
coupling constants are obtained. The sets which turn out to be singular are discarded. The valid solutions are averaged and the standard deviation is calculated. 
The solutions that are deviating by more than three standard deviations are discarded and the standard deviation is calculated again and this cycle is 
repeated until self-consistency is obtained. </p>


<p>The coupling constants obtained this way are then used to further determine the energy of all the possible electronic states (2<sup>0.5x</sup> for 
x paramagnetic centres). To achieve this, the Heisenberg Dirac van Vleck Hamiltonian (H&#770 = -2 &sum;<sub>i>j</sub> J<sub>ij</sub>s<sub>i</sub>s<sub>j</sub>) is used.
<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs#Heisenberg" class="showTip Heisenberg">Heisenberg</a>], [<a href="refs#Dirac" 
class="showTip Dirac">Dirac</a>], [<a href="refs#Vanvleck" class="showTip Vanvleck">Vanvleck</a>]</font>  
The spin operator value for each paramagnetic centre is determined by taking an average of the spin operator value in the different electronic states that 
is provided as input. Two electronic states are considered to be degenerate if the difference in the energy between the two is less than 1 cm<sup>-1</sup>.</p>


</body></html>
Welcome to the J2suscep wiki!

<h1>J2suscep - Documentation</h1>
<h2>Table of Contents</h2>
<li><A HREF="j2suscep"><strong>J2suscep</strong></A>   </li>
<li><A HREF="Installation"><strong>Installation Instructions</strong></A>   </li>
<li><A HREF="Execution"><strong>Code Execution</strong></A>   </li>
<li><A HREF="contribution"><strong>Contributing</strong></A>   </li>
<h3> EJ_Calc </h3>
<ul>
  <li><A HREF="1a_introduction"><strong>Introduction</strong></A>   </li>
  <li><strong><a href="1b_theory">Theory</a></strong>
  <li><strong><A HREF="1c_basic_details">Basic Details</A></strong></li>
  <li><strong><A HREF="1e_example">Example Input Files</A></strong></li>
  <li><strong><A HREF="1f_output_file_description">Description of Output File</A></strong></li>
  <li><strong><A HREF="1g_file_summary">Summary of Input and Output Files</A></strong></li>
  <li><strong><A HREF="1h_developers">For Developers</A></strong>  </li>
  </ul>
</li>

<h3>Suscep </h3>
<ul>
  <li><A HREF="2a_introduction"><strong>Introduction</strong></A>   </li>
  <li><strong><a href="2b_theory">Theory</a></strong>
  <li><strong><A HREF="2c_basic_details">Basic Details</A></strong></li>
  <li><strong><A HREF="2e_example">Example Input Files</A></strong></li>
  <li><strong><A HREF="2f_output_file_description">Description of Output File</A></strong></li>
  <li><strong><A HREF="2g_limitations">Limitations</A></strong></li>
  <li><strong><A HREF="2h_developers">For Developers</A></strong>  </li>
  </ul>
</li>
<li><A HREF="worked_example"><strong>Worked Example</strong></A>   </li>
<li><a href="refs"><strong>References</a></strong></li>
</body></html>
<html><head></head><body>
<h1>Code execution for ej_calc</h1>

<p>The execution needs to be done in two steps.</p>
<br>
<p>(i)	The first step is to execute ej_calc_form. If the input files are <a href="#format-of-file1">File1</a> and <a href="#format-of-file2">File2</a>, then the following command can be used:: </p>
<pre>ej_calc_form File1 File2</pre>
<p>It is necessary to run ej_calc_form to produce the singul file which is required to calculate the coupling constant using ej_calc_spin.</p>

<br>
<p>(ii)	This is followed by the execution of ej_calc_spin as:</p>
<pre>ej_calc_spin File1 File2</pre>
<br>

<p><a href="#format-of-file1">File1</a> contains the definition of the Hamiltonian.</p>
<p><a href="#format-of-file2">File2</a> contains the value of spin operators for each paramagnetic centre and the energy of different spin states.</p>
<p>The details about the coupling constants are printed in the files File1_form.out and File1_spin.out for formal spins and spin densities respectively.</p>

<p>The details of the format for <a href="#format-of-file1">File1</a> and <a href="#format-of-file2">File2</a> are provided later in this section.</p>
<p></p>

<h1>Code execution for suscep</h1>

<p>The execution of the suscep code requires three input files, and for input files named <a href="#format-of-file1">File1</a>, <a href="#format-of-file3">File3</a> and <a href="#format-of-file4">File4</a>, the suscep code can be executed in the following way:</p>
<pre>suscep File1 File3 File4</pre>
<p><a href="#format-of-file1">File1</a> contains the definition of the Hamiltonian.</p>
<p><a href="#format-of-file3">File3</a> contains the additional information other than the J-values required to calculate the temperature dependence of magnetic susceptibility.</p>
<p><a href="#format-of-file4">File4</a> contains the magnitude of J-values to be used for the determination of magnetic susceptibility.</p>
<p>The temperature dependence of magnetic susceptibility is printed in the file File1.out.</p>
<p>The details of the format for the input files are provided later in this section.</p>



<h1>Code execution using script</h1>
<p>For ease, a bash script (j2suscep) has also provided in the bin directory that can be used to execute the codes separately or together as described below:</p>
<h2> EJ_Calc <h2>

<p>The code ej_calc_spin can be executed using the '-s' tag for input files <a href="#format-of-file1">File1</a> and <a href="#format-of-file2">File2</a> as follows:</p>
<pre>j2suscep -s File1 File2</pre>
<p>The use of the '-s' tag leads to the execution of either ej_calc_spin (if a 'singul' file exists in the working directory), or both ej_calc_form and ej_calc_spin (if a 'singul' file does not exist in the working directory). </p>


<p>To execute the code ej_calc_form as a separate step, one can use the '-f' tag. For input files <a href="#format-of-file1">File1</a> and <a href="#format-of-file2">File2</a>, the execution can be performed as follows:
</p>
<pre>j2suscep -f File1 File2</pre>




<h2> Suscep <h2>
<p>For input files <a href="#format-of-file1">File1</a>, <a href="#format-of-file3">File3</a> and <a href="#format-of-file4">File4</a>, the code suscep can be executed as follows:
<pre>j2suscep -S File1 File3 File4</pre>

<h2> Combined execution <h2>
<p>Both the ej_calc and suscep codes can also be executed in one go. For this, again three files are required:</p>

<blockquote>
<ol>
<li><a href="#format-of-file1">File1</a>, which contains the details about the Hamiltonian.</li>
<li><a href="#format-of-file2">File2</a>, which contains information about the the value of spin operators for each paramagnetic centre and the energy of different spin states.</li>
<li><a href="#format-of-file3">File3</a>, which contains the additional information other than the J-values required to calculate the temperature dependence of magnetic susceptibility.</li>
</ol>
</blockquote>

<p> Although another file <a href="#format-of-file4">(File4)</a> is generally required for the execution of the suscep code, a combined execution of both codes generates this file automatically. The execution is performed as follows: </p>

<pre>j2suscep -a File1 File2 File3</pre>

<p>The j-values determined using the ej_calc_spin codes are written to a file named 'jvals' and are used for the calculation of the temperature dependence of magnetic susceptibility. </p>

<p>The usage information for this script can also be accessed as follows:</p>
<pre>j2suscep -h</pre>
<p>The following information is displayed with the use of the -h tag:</p>
<pre>
***J2suscep help***


For executing ej_calc_spin:
 j2suscep -s File1 File2

For executing ej_calc_form:
 j2suscep -f File1 File2

For executing suscep:
 j2suscep -S File1 File3 File4

For executing all codes together:
 j2suscep -a File1 File2 File3

For help:
 j2suscep -h

Information expected within different input files
File1 should contain the details about the Hamiltonian.
File2 should contain information about the the value of spin operators for each paramagnetic centre and the energy of different spin states.
File3 should contain the additional information other than the J-values required to calculate the temperature dependence of magnetic susceptibility.
File4 should contain the magnitude of J-values to be used for the determination of magnetic susceptibility.
</pre>


<h1>Format of input files</h1>

<h2>Format of File1</h2>
<p>File1 contains the definition of the Hamiltonian. The format of this file is as follows:</p>

<pre>
magnetic centres
&lt;number of paramagnetic centres&gt;
No. of J values
&lt;number of j-value>
Hamiltonian
&lt;First interaction associated with the first j-value&gt;
&lt;Second interaction associated with the first j-value&gt;
&lt;Third interaction associated with the first j-value&gt;
And so on
****
&lt;First interaction associated with the second j-value&gt;
&lt;Second interaction associated with the second j-value&gt;
&lt;Third interaction associated with the first j-value&gt;
And so on
****
&lt;First interaction associated with the third j-value&gt;
&lt;Second interaction associated with the third j-value&gt;
&lt;Third interaction associated with the first j-value&gt;
And so on
Hamiltonian Ends

</pre>


<p>The first line has to state 'magnetic centres' and the second line will specify how many magnetic centres are there. 
Then the number of J-values are specified. After that the Hamiltonian is specified.</p>
<p>The spin interactions under a single J-value are specified in a single block. 
The spin interactions for different J-values are separated by '****'. The line 'Hamiltonian Ends' is needed in the end and is not preceded by '****'.</p>
<p><b>Note:</b> The terms ‘magnetic centres’, ‘No. of J values’, ‘Hamiltonian’ and ‘Hamiltonian Ends’ are case-sensitive.</p>
<p>Examples of input files are provided in the <A HREF="1e_example">Example Input Files</A> subsection under the 'EJ_Calc' section and in the <A HREF="2e_example">Example Input Files</A> subsection under the 'Suscep' section.</p>

<br>
<h2>Format of File2</h2>
<p>File2 contains information about the the value of spin operators for each paramagnetic centre and the energy of different spin states. The format of this file is as follows:</p>

<pre> 
&lt;S<sub>11</sub>&gt;  &lt;S<sub>12</sub>&gt;  &lt;S<sub>13</sub>&gt; ………. &lt;E<sub>1</sub>&gt;
&lt;S<sub>21</sub>&gt;  &lt;S<sub>22</sub>&gt;  &lt;S<sub>23</sub>&gt; ………. &lt;E<sub>2</sub>&gt;
&lt;S<sub>31</sub>&gt; &lt;S<sub>32</sub>&gt; &lt;S<sub>33</sub>&gt; ………. &lt;E<sub>3</sub>&gt;
…..
….

</pre>

<p>Here &lt;S<sub>ij</sub>&gt; denotes the spin operator value for the j<sup>th</sup> paramagnetic centre in the i<sup>th</sup> electronic state and &lt;E<sub>i</sub>&gt; represents the energy of the i<sup>th</sup> state in <b>Hartrees</b>.</p>


<p>Examples of this input file are provided in the <A HREF="1e_example">Example Input Files</A> subsection under the 'EJ_Calc' section.</p>

<p></p>
<br>


<h2>Format of File3</h2>
<p>File3 contains the additional information other than the J-values required to calculate the temperature dependence of magnetic susceptibility and the format is as follows:</p>

<pre>
spin
&lt;spin on the first paramagnetic centre&gt;
&lt;spin on the second paramagnetic centre&gt;
&lt;spin on the third paramagnetic centre&gt;
And so on
g value
&lt;g value&gt;
Field Strength
&lt;Field Strength (in Oersted)&gt;
Temperature range
&lt;Minimum temperature (in Kelvin)&gt;
&lt;Maximum temperature (in Kelvin)&gt;
Step size
&lt;Interval for temperature increments&gt;
</pre>


<p>This is the format of File3. The first line has to state 'spin' and then the spin (number of unpaired electrons/2) on each paramagnetic centre is specified. This is followed by the definition of the g value which is then followed by the field strength (in Oersted) and temperature (in Kelvin). The final term is the step size for temperature (in Kelvin).</p>

<p><b>Note:</b> The terms ‘magnetic centres’, ‘spin’, ‘No. of J values’, ‘g value’, ‘Hamiltonian’, ‘Field Strength’, 'Temperature range' and 'Step size' are case-sensitive.</p>
<p>Examples of this input file are provided in the <A HREF="2e_example">Example Input Files</A> subsection under the 'Suscep' section.</p>
<br>

<h2>Format of File4</h2>
<p>This file contains the magnitude of J-values to be used for the determination of magnetic susceptibility. The format of this file is shown below:</p>
<pre>
J values
&lt;1st J-value (in cm-1)&gt;
&lt;2nd J-value (in cm-1)&gt;
&lt;3rd J-value (in cm-1)&gt;
And so on
</pre>

The first line contains the term ‘J values’ is followed by the strength of each J-value (in cm<sup>-1</sup>) specified on separate lines.
<p><b>Note:</b> The term ‘J values’ is case-sensitive.</p>

<p>Examples of this input file are provided in the <A HREF="2e_example">Example Input Files</A> subsection under the 'Suscep' section.</p>
<p></p>





</body></html>



    
</body></html>
<html><head></head><body>
<h1>Basic Details</h1>

<p>This code has been written in FORTRAN 2008 and has two executables – ej_calc_form and ej_calc_spin. 
The calculation of coupling constants using spin densities requires the use of both executables and the execution of ej_calc_form should 
be followed by execution of ej_calc_spin. Each executable requires two input files – one specifying the Hamiltonian for the system under study 
(.inp file) and another that contains the definition of spin operator (i.e. the spin density) for each paramagnetic centre and the energy 
of the system in the given spin state.</p>


<p>The executable ej_calc_form makes use of formal spin values to determine the j-values while ej_calc_spin uses the actual value of spin operators 
for the same purpose. This has been done because the value of spin operators are generally not whole numbers and it is difficult to identify singular 
sets of equations without the use of whole numbers. Since the actual values of spin operators are generally not far off from the nearest whole number, 
it is highly unlikely that the use of spin operators will remove the singularity. ej_calc_form produces 4 output files – a .out file containing 
the coupling constants and 3 other files – bigg, nrg_frm_jval_formal and singul. ej_calc_spin produces 2 output files - a .out file containing the coupling constants and nrg_frm_jval_spin. The details about the coupling constants is printed in the .out file. 
Information contained in the other files is described in the following discussion.</p>



<p>One starts by executing ej_calc_form which rounds the spin operator values to the nearest whole number and uses the corresponding whole number as the 
value of the spin operator. It identifies which sets of equations are singular, discards them and stores information about these sets in the file ‘singul’. 
In rare occasions, equations are encountered which are close to singularity thereby leading to absurd (non-infinite) solutions. To keep note of these equations, 
details are printed in the file ‘bigg’. The criteria to identify such equations has been set such that any coupling stronger than 50 cm<sup>-1</sup> is considered too big. 
These solutions are not discarded as in some cases coupling constants can be quite large. If, however, they are outliers, they will be removed eventually since at 
the end, only those solutions will be considered valid which lie within 3 standard deviations of the average value as mentioned previously. 
ej_calc_spin reads the ‘singul’ file generated by ej_calc_form and calculates the solutions for only the sets that were found to be non-singular by ej_calc_form.</p>


<p>The files nrg_frm_jval_formal and nrg_frm_jval_spin contain the approximate energy of the all possible 2<sup>0.5x</sup> states (x = number of paramagnetic centres). 
The energy is calculated using the Heisenberg Dirac van Vleck Hamiltonian (H&#770 = -2 &sum;<sub>i>j</sub> J<sub>ij</sub>s<sub>i</sub>s<sub>j</sub>). 
<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs#Heisenberg" class="showTip Heisenberg">Heisenberg</a>], [<a href="refs#Dirac" 
class="showTip Dirac">Dirac</a>], [<a href="refs#Vanvleck" class="showTip Vanvleck">Vanvleck</a>]</font> 
The spin operators for all the states are defined by taking the average of the value of spin operator for each paramagnetic centre that was provided as input. 
Since the variation of the value of spin operators in different spin states is small, the energy obtained this way serves as a reasonable approximation to the 
actual energy of these states. The true energies of all spin states can, of course, be obtained by diagonalising the full Hamiltonian. Since this diagonalisation 
becomes computationally quite expensive for large systems, the method employed in this code provides an inexpensive way to obtain a decent approximation of the 
energy of these states.</p>



<p></p>


    
</body></html>
<html><head></head><body>
<h1>Limitations</h1>


<p>The determination of the susceptibility for very large systems is limited by the available memory (RAM). This is because the size of the matrix 
containing M<sub>s</sub> elements increases rapidly with the increase in the number of unpaired electrons and paramagnetic centres. The total 
number of elements in the matrix containing the M<sub>s</sub> elements can be given as [(no. of unpaired electrons on centre 1) * (no. of unpaired 
electrons on centre 2) * (no. of unpaired electrons on centre 3) * ……...]<sup>2</sup>. For example, for a {Mn<sup>III</sup><sub>4</sub>} system, 
the matrix contains a total of (5*5*5*5)<sup>2</sup> =) 390,625 elements, for a {Mn<sup>III</sup><sub>5</sub>} system, the matrix contains a total 
of (5<sup>5</sup>*5<sup>5</sup> =) 9,765,625 elements while for a {Mn<sup>III</sup><sub>6</sub>} system, the matrix contains (5<sup>6</sup>*5<sup>6</sup> =) 
244,140,625 elements. Since the calculation of susceptibility requires the exact diagonalisation of the whole matrix, the matrix needs to be stored 
in the memory and hence the limitation.</p>

<p></p>

    
</body></html>
<html><head></head><body>
<h1>For Developers</h1>
<h2>Workflow</h2>
<p>The ej_calc code, is divided into 4 modules:</p>
<p>1.	main.F90: controls the flow of the code.</p>
<p>2.	init.F90: declares all the global variables, reads the input files and initializes the global variable using the information obtained from the input files.</p>
<p>3.	modul_lib.F90: builds up the Hamiltonian based on the definition obtained from the input file. It further forms different combinations of the equations and solves them using LAPACK routines. The average of all the valid sets of j-values and the standard deviations for each J-value is also determined using this module.</p>
<p>4.	enrg_frm_jval.F90: calculates the energy of the all spin states using the Heisenberg Dirac van Vleck Hamiltonian (H&#770 = -2 &sum;<sub>i>j</sub> J<sub>ij</sub>s<sub>i</sub>s<sub>j</sub>).
<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs#Heisenberg" class="showTip Heisenberg">Heisenberg</a>], [<a href="refs#Dirac" 
class="showTip Dirac">Dirac</a>], [<a href="refs#Vanvleck" class="showTip Vanvleck">Vanvleck</a>]</font>  </p> 
<br>

<p>A detailed breakdown of each module is provided below.</p>

<h3>init.F90</h3>
<p>This module defines all the global variables.</p>

<table class="tg">
<thead>
  <tr>
    <th class="tg-amwm"> Variable </th>
    <th class="tg-amwm"> Function </th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0lax"> num_mag_cent </td>
    <td class="tg-0lax"> stores the number of metal centres in each configuration. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> no_of_j_val </td>
    <td class="tg-0lax"> stores the total number of j-values defined for the system. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> tot_interac </td>
    <td class="tg-0lax"> stores the number of possible interactions between spins i.e. the total number of S1S2 terms. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> num_spin_dens_set </td>
    <td class="tg-0lax"> stores the number of sets of spin density available which is the total number of lines in the spin density file. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> try &nbsp;</td>
    <td class="tg-0lax"> this helps in keeping the program running just in case there is an equation set is encountered that cannot be solved. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> poss_comb </td>
    <td class="tg-0lax"> stores the possible combination of equations i.e. the <sup>(n-1)</sup>C<sub>j</sub> term. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> main_cntr </td>
    <td class="tg-0lax"> keep track of total non-singular equations.</td>
  </tr>
  <tr>
    <td class="tg-0lax"> start, finish </td>
    <td class="tg-0lax"> These in combination determine the time taken during the program execution.</td>
  </tr>
</tbody>
</table>
<br>


<table class="tg">
<thead>
  <tr>
    <th class="tg-amwm"> Array </th>
    <th class="tg-amwm"> Function </th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0lax"> hamil1, hamil2 </td>
    <td class="tg-0lax"> These in combination keep track of which centres are interacting . </td>
  </tr>
  <tr>
    <td class="tg-0lax"> jpos </td>
    <td class="tg-0lax"> keeps track of how many interactions are being included under a single j-value. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> energy </td>
    <td class="tg-0lax"> stores the energy for each electronic configuration. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> tot_avg </td>
    <td class="tg-0lax"> stores the average of all the J-value sets which are valid and which meet the standard deviation criteria. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> stddev </td>
    <td class="tg-0lax"> stores the standard deviation for the different J-values. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> std_per </td>
    <td class="tg-0lax"> stores the percentage standard deviation for the different J-values. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> spin </td>
    <td class="tg-0lax"> stores the spin operator value associated with each metal centre as given in the spin operator file. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> hamil </td>
    <td class="tg-0lax"> stores the 2S1S2 type weight terms for each interaction. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> jval </td>
    <td class="tg-0lax"> stores the total (2S1S2 type) coefficient associated with each j-value. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> jval_trkr </td>
    <td class="tg-0lax"> stores all the computed sets of j-values. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> singul </td>
    <td class="tg-0lax">   (for ej_calc_spin only) stores the location of singular sets of equations identified by ej_calc_form. It essentially stores the data in the file singul.      </td>
  </tr>
</tbody>
</table>
<br>
<br>

<h3>modul_lib.F90</h3>
<p>This module contains all the subroutines required to calculate all the possible sets of j-values and determine the valid ones based on standard deviation. 
The following subroutines form part of this module:</p>
<p><b>1.	jvals:</b> This subroutine will determine the S1S2 terms and which of these terms comes under a given j-value.</p>
<p><b>2.	solv:</b> This subroutine forms all possible sets of equations and solves them. The local variables used in this subroutine are as follows:</p>

<table class="tg">
<thead>
  <tr>
    <th class="tg-amwm"> Array </th>
    <th class="tg-amwm"> Function </th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0lax"> ctr </td>
    <td class="tg-0lax"> allows the separation and hence the evaluation of all possible combinations. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> loc_energy </td>
    <td class="tg-0lax"> stores the energy locally since that will be changed in each set of equations considered. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> enrg_tb_slvd1 </td>
    <td class="tg-0lax"> local copy of the global array ‘energy’. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> enrg_tb_slvd </td>
    <td class="tg-0lax"> stores the energy of the set of equations that will be solved in one instance. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> enrg_tb_slvd2 </td>
    <td class="tg-0lax"> (for ej_calc_form only) copy of the array ‘enrg_tb_slvd’. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> loc_jval </td>
    <td class="tg-0lax"> store the coefficients of j-values locally since that will be changed in each set of equations considered. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> jval_tb_solvd1 </td>
    <td class="tg-0lax"> local copy of the global array ‘jval’. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> jval_tb_solvd </td>
    <td class="tg-0lax"> stores the coefficients of j-values of the set of equations that will be solved in one instance. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> jval_tb_solvd2 </td>
    <td class="tg-0lax"> (for ej_calc_form only) copy of the array ‘jval_tb_solvd’. </td>
  </tr>
</tbody>
</table>
<br>

<p>In this subroutine, first of all, an electronic state is chosen as the reference state (every state is chosen as a reference state once) and a new set 
of equations is formed by subtracting the reference state from all the other states. These states are then grouped into all possible sets of r (the number 
of j-values one is looking for) and each such set is solved as simultaneous linear equations using the LAPACK library routine dgesv. The LAPACK routine 
dgesv is implemented in a different subroutine called simul_solv which also performs a few basic checks for singularity. If the solution is found to be 
singular, the value of all j-values for the given set is set to zero. If, however, the set is not found to be singular, a more thorough screening of the 
set is performed using the LAPACK library routine dgesvx which is implemented in the subroutine simul_solv1. Additionally, if any of the j-value is found to be larger than 50 cm<sup>-1</sup> for a given set of equations, details about the set and the corresponding J-values obtained are printed in the file bigg. In case of ej_calc_spin however, only those sets that are determined as non-singular by the array singul (generated from the 'singul' file created by executing ej_calc_formal) are solved.</p>
<br>


<p><b>3.	checkpos:</b> This subroutine is a helper subroutine for the subroutine ‘solv’ to cycle through all the possible combinations of r (the number of 
j-values one is looking for) spin states from n (or more accurately n-1) spin states. </p>
<p><b>4.	simul_solv:</b> This subroutine is also a helper subroutine for the subroutine ‘solv’. It checks if there are any identical pair of equations which will 
automatically make the whole set singular. If such pairs are not found, then a solution is calculated using the LAPACK library routine dgesv.</p>
<p><b>5.	simul_solv1:</b> (for ej_calc_form only) This subroutine is also a helper subroutine for the subroutine ‘solv’. It calculates the solution for the same set of equations that were considered solvable by simul_solv using the LAPACK library routine dgesvx which looks further for sets of equations close to singularity. The reason why this is used only for the sets which are found to be solvable by dgesv because it is more computationally expensive (The speed gain in keeping simul_solv and simul_solv1 has not 
been determined and considering that this code can still solve ~1,000,000 sets of equations within minutes, it was not considered worth doing). </p>
<p><b>6.	soln_chk:</b> (for ej_calc_form only) This subroutine makes sure that the code does not get stuck on any set of equations. </p>
<p><b>7.	avrg_calcs:</b> This subroutine calculates the average of all the valid set of solutions. </p>
<p><b>8.	std_dev:</b> This subroutine calculates the standard deviation for all the valid set of solutions. </p>
<p><b>9.	std_dev1:</b> This subroutine checks if all the valid j-value sets are within 3 standard deviations or not and if they are not, they are removed. </p>
<p><b>10.	avrg_calcs1:</b> This subroutine will calculate the new average for j-values since some have been eliminated by the subroutine ‘std_dev1’.  </p>
<p><b>11.	final_print:</b> This subroutine prints the final set of J-values. </p>
<p><b>12.	backtrack1:</b> This subroutine will calculate the energies using the computed J-values and compares them with the energies that are provided as input.</p>
<br>
<br>


<h3>enrg_frm_jval.F90</h3>
<p>This module calculates the energy of the spin states using the coupling constants calculated by modul_lib.F90. 
The main local variables used in this module are as follows:</p>

<table class="tg">
<thead>
  <tr>
    <th class="tg-amwm"> Array </th>
    <th class="tg-amwm"> Function </th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0lax"> unp_elec </td>
    <td class="tg-0lax"> stores the spin operator value for each metal centre. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> all_config </td>
    <td class="tg-0lax"> stores all the possible spin configurations of spin states. </td>
  </tr>
  <tr>
    <td class="tg-0lax"> new_nrg </td>
    <td class="tg-0lax"> stores the energy of all possible spin states. </td>
  </tr>
  <tr>
    <td class="tg-amwm"> Variable </td>
    <td class="tg-amwm"> Function </td>
  </tr>
  <tr>
    <td class="tg-0lax"> poss_comb1 </td>
    <td class="tg-0lax"> stores the total possible spin states. </td>
  </tr>
</tbody>
</table>
<br>

<p>The following subroutines form part of this module:</p>
<p><b>1.	avg_spin:</b> This subroutine calculates the average spin operator (in this case, the formal spin) value for each metal centre. </p>
<p><b>2.	all_poss_spins:</b> This subroutine along with the subroutine ‘update_spin’ forms all possible spin configurations, stores them in the array 
‘all_config’, determines the energy of each configuration in accordance with the specified hamiltonian and stores the energy in the array ‘new_nrg’. </p>
<p><b>3.	sort:</b> This subroutine sorts all the spin states on the basis of increasing energy using bubble sort. </p>
<p><b>4.	swap:</b> This is a helper subroutine for the subroutine ‘sort’. </p>
<p><b>5.	prnt:</b> This subroutine prints the sorted energies in the nrg_frm_jval file. </p>
<p><b>6.	update_spin:</b> This is a helper subroutine for the subroutine ‘all_poss_spins’ for cycling through all possible spin configurations.</p>
<p><b>7.	std_dev2d:</b> (for ej_calc_spin only) calculates the standard deviation for the spin operator values that are used to calculate the energy of all spin states.</p>
<br>
<br>

    
</body></html>
<html><head></head><body>
<h1>Introduction</h1>
<p>Briefly, this code calculates coupling constants using the following Hamiltonian:
<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs#Heisenberg" class="showTip Heisenberg">Heisenberg</a>], [<a href="refs#Dirac" 
class="showTip Dirac">Dirac</a>], [<a href="refs#Vanvleck" class="showTip Vanvleck">Vanvleck</a>]</font>  </p> 
<br>

<font face="Courier New, Courier, monospace" size="-2">&nbsp;&nbsp;&nbsp;&nbsp;</font>H&#770; = -2 &sum;<sub>i>j</sub> J<sub>ij</sub>s<sub>i</sub>s<sub>j</sub>	
<p>where the terms have the following meanings:</p>
<table border="0" cellpadding="2">
<tr valign="top"><td><p>J<sub>ij</sub></p></td><td>&nbsp;&nbsp;</td>
<td><p>The exchange coupling constant</p></td></tr>
<tr valign="top"><td><p>s<sub>x</sub></p></td><td>&nbsp;&nbsp;</td>
<td><p>The spin operator.</p></td></tr>
</table>
<p> Any definition of spin operators (formal spins or spin projections or spin densities or any other) can be used. Sometimes researchers exclude the pre-factor 2 in the above formula so this needs to be kept in mind while using this code.</p>

    
</body></html>
<html><head></head><body>


<h1>J2suscep</h1>

<p> This package can be used to calculate the exchange coupling constants between 
paramagnetic centres and their influence on the temperature dependence of magnetic 
susceptibility for a given system.</p>
<p>This package contains two standalone codes written in FORTRAN 2008 – ej_calc and suscep. 
The ej_calc code uses the spin densities obtained from DFT calculations and determines the 
isotropic exchange coupling between paramagnetic centres. The suscep code calculates the 
temperature dependence of magnetic susceptibility using the coupling constants. 


    
</body></html>
<html><head></head><body>
<h1>Theory</h1>

<p></p>

<p><b>Matrix Elements:</b> The matrix elements are constructed using an uncoupled basis comprising of spin elements (the magnetic spin quantum number, M<sub>s</sub>, value) 
for each magnetic centre. If a magnetic centre has N unpaired electron, then the magnetic spin quantum number values for that particular magnetic centre range from N/2 
to –N/2; the consecutive values differing from each other by unity (total (N+1) values).</p>

<p><b>Hamiltonian:</b> Only the two major contributors – the isotropic exchange term and the Zeeman term – to the Hamiltonian have been taken into consideration.</p>
<p>Isotropic Exchange term:<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs#Boca1" class="showTip Boca1">Boca1</a>]</font>  </p>

<p align="center"><img src="img/two.jpg" alt="eq. 2"><br /></p>
<p>where <i>J<sub>AB</sub></i> is the coupling constant between the magnetic centres A and B, <i>I</i> and <i>J</i> represent the particular M<sub>s</sub> terms in the basis elements 
that <i>S&#770<sub>A</sub>.S&#770<sub>B</sub></i> operate on and</p>

<p align="center"><img src="img/three.jpg" alt="eq. 3"><br /></p>
<p>where <i>S&#770<sub>A,+</sub></i> and <i>S&#770<sub>B,-</sub></i> are increment and decrement operators that are defined as follows: </p>

<p align="center"><img src="img/four_six.jpg" alt="eq. 4-6"><br /></p>
<br>

<p>Zeeman term:<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs#Boca1" class="showTip Boca1">Boca1</a>]</font> </p>
<p align="center"><img src="img/seven.JPG" alt="eq. 7"><br /></p>
<p>Here <i>B&#8407</i> represents the applied magnetic field, g<sub>A</sub> is the gyromagnetic tensor and μ<sub>B</sub> is the unit Bohr Magneton.</p>

<p align="center"><img src="img/eight.jpg" alt="eq. 8"><br /></p>
<p>where δ<sub>AB</sub> is the kroenecker delta operator which is equal to 1 if A=B. Otherwise its value is 0. <i>B<sub>i</sub></i> and g<sub>Ai</sub> represent the value of 
the magnetic field and the gyromagnetic tensor along a given direction. The effective Hamiltonian is thus given as </p>

<p align="center"><img src="img/nine.JPG" alt="eq. 9"><br /></p>
<p>The Hamiltonian matrix can thus be represented as </p>

<p align="center"><img src="img/ten.JPG" alt="eq. 10"><br /></p>
<p>where <i>A<sub>I</sub></i> and <i>A<sub>J</sub></i> represent the elements of the basis and N is the total number of elements in the basis. </p>

<p>The eigenvalues ε<sub>i</sub>  and hence the energy levels of the various magnetic states are obtained by the diagonalisation of the Hamiltonian matrix. 
Once the eigenvalues are obtained, it can be used to plot the magnetic susceptibility with respect to temperature. To calculate the magnetic susceptibility, 
the Van Vleck equation is used which is:<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs#Kahn" class="showTip Kahn">Kahn</a>]</font> </p>
<p align="center"><img src="img/eleven.JPG" alt="eq. 11"><br /></p>
<p>where <i>ϵ<sub>i</sub><sup>(1)</sup></i> and <i>ϵ<sub>i</sub><sup>(2)</sup></i> are the first and second derivative of energy with respect to the magnetic induction, μ<sub>0</sub> 
is the permeability of free space, N<sub>A</sub> is the Avogadro’s number, <i>k</i> is the Boltzmann constant and <i>T</i> is the temperature. </p>
<p>By Taylor expansion, the eigenvalues can be represented as:<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs#Boca" class="showTip Boca">Boca</a>]</font> </p>

<p align="center"><img src="img/twelve.JPG" alt="eq. 12"><br /></p>
<p>with</p>

<p align="center"><img src="img/thirteen.JPG" alt="eq. 13"><br /></p>
<p>being the coefficients to be determined. To obtain these coefficients, 5 sets of eigenvalues, <i>ε<sub>i,m</sub></i> are generated using 5 different magnetic fields <i>B<sub>m</sub></i>
, where</p>

<p align="center"><img src="img/fourteen.JPG" alt="eq. 14"><br /></p>
<p>where N = 0 and 1, and δ has been chosen to be equal B<sub>0</sub>/10. <br>Therefore,</p>

<p align="center"><img src="img/fifteen.JPG" alt="eq. 15"><br /></p>
<p>for</p>

<p align="center"><img src="img/sixteen.JPG" alt="eq. 16"><br /></p>
<p>In matrix form,</p>

<p align="center"><img src="img/seventeen.JPG" alt="eq. 17"><br /></p>
<p>or</p>


<p align="center"><img src="img/eighteen.JPG" alt="eq. 18"><br /></p>
<p>The matrix of coefficients can be calculated, provided the inverse of B exists, as follows:</p>

<p align="center"><img src="img/nineteen.JPG" alt="eq. 19"><br /></p>
<p>By assembling the rows together, we obtain,</p>

<p align="center"><img src="img/twenty.JPG" alt="eq. 20"><br /></p>
<p>i.e., </p>


<p align="center"><img src="img/twenty_one.JPG" alt="eq. 21"><br /></p>
<p>Thus all the coefficients can be determined by the following equation:</p>

<p align="center"><img src="img/twenty_two.JPG" alt="eq. 22"><br /></p>
<p>These coefficients can then be used to calculate the magnetic susceptibility using the Van Vleck equation which can now be given as:</p>

<p align="center"><img src="img/twenty_three.JPG" alt="eq. 23"><br /></p>
<p></p>



</body></html>
<html><head></head><body>
<h1>References</h1>
<table border="1">
<tr><td valign="top">Boca</td><td>Boca, R., in <i>Temperature dependence of magnetic susceptibility</i>, <i>Theoretical foundations of molecular magnetism</i>, Amsterdam; Elsevier, Oxford, 1999, 315-344.</td></tr>
<tr><td valign="top">Boca1</td><td>Boca, R., in <i>Clusters</i>, <i>Theoretical foundations of molecular magnetism</i>, Amsterdam; Elsevier, Oxford, 1999, 701-836.</td></tr>
<tr><td valign="top">CCDC</td><td>Groom, C. R., Bruno, I. J., Lightfoot, M. P. & Ward, S. C., &ldquo;The Cambridge Structural Database,&rdquo; <i>Acta Cryst. B</i>, <b>72</b> (2016) 171-179.</td></tr>
<tr><td valign="top">Cremades</td><td>Cremades, E., J. Cano, E. Ruiz, G. Rajaraman, C. J. Milios and E. K. Brechin, &ldquo;Theoretical Methods Enlighten Magnetic Properties of a Family of Mn<sub>6</sub> Single-Molecule Magnets,&rdquo; <i>Inorg. Chem.</i>, <b>48</b>(16) (2009) 8012-8019.</td></tr>
<tr><td valign="top">Cremades1</td><td>Cremades, E., T. Cauchy, J. Cano and E. Ruiz, &ldquo;Can theoretical methods go beyond the experimental data? The case of molecular magnetism,&rdquo; <i>Dalton Trans.</i>, <b>30</b> (2009) 5873-5878.</td></tr>
<tr><td valign="top">Dirac</td><td>Dirac, P. A. M., &ldquo;Quantum Mechanics of Many-Electron Systems,&rdquo; <i>Proc. R. Soc. Lond. A</i>, <b>123</b>(792) (1929) 714-733.</td></tr>
<tr><td valign="top">Heisenberg</td><td>Heisenberg, W., &ldquo;Zur Theorie des Ferromagnetismus,&rdquo; <i>Z. Physik</i>, <b>49</b>(9-10) (1928) 619-636.</td></tr>
<tr><td valign="top">Kahn</td><td>Kahn, O., in <i>Magnetization and magnetic susceptibility</i>, <i>Molecular magnetism</i>, USA; VCH Publishers, Inc, 1993, 1-8.</td></tr>
<tr><td valign="top">Libby</td><td>E. Libby, R. J. Webb, W. E. Streib, K. Folting, J. C. Huffman, D. N. Hendrickson and G. Christou, &ldquo;Crystal structure and magnetic susceptibility of the dinuclear manganese(IV) complex Mn<sub>2</sub>O<sub>2</sub>(pic)<sub>4</sub>.MeCN (picH = picolinic acid),&rdquo; <i>Inorg. Chem.</i>, <b>28</b> (1989) 4037-4040.</td></tr>
<tr><td valign="top">Rajeshkumar</td><td>Rajeshkumar, T., H. V. Annadata, M. Evangelisti, S. K. Langley, N. F. Chilton, K. S. Murray and G. Rajaraman, &ldquo;Theoretical Studies on Polynuclear {Cu<sup>II</sup><sub>5</sub>Gd<sup>III</sup><sub>n</sub>} Clusters (n=4, 2): Towards Understanding Their Large Magnetocaloric Effect,&rdquo; <i>Inorg. Chem.</i>, <b>54</b>(4) (2015) 1661-1670.</td></tr>
<tr><td valign="top">Tandon</td><td>Tandon, S., M. Venkatesan, W. Schmitt and G. W. Watson, &ldquo;Altering the nature of coupling by changing the oxidation state in a {Mn<sub>6</sub>} cage,&rdquo; <i>Dalton Trans.</i>, <b>49</b>(24) (2020) 8086-8095.</td></tr>
<tr><td valign="top">Vanvleck</td><td>Vanvleck, J. H., &ldquo;A Survey of the Theory of Ferromagnetism,&rdquo; <i>Rev. Mod. Phys.</i>, <b>17</b>(1) (1945) 27-47.</td></tr>
<tr><td valign="top">Vignesh</td><td>Vignesh, K. R., S. K. Langley, K. S. Murray and G. Rajaraman, &ldquo;What Controls the Magnetic Exchange Interaction in Mixed- and Homo-Valent Mn<sub>7</sub> Disc-Like Clusters? A Theoretical Perspective,&rdquo; <i>Chem. Eur. J.</i>, <b>21</b>(7) (2015) 2881-2892.</td></tr>
<tr> 
</table>
</body>
</html>
<h2>Table of Contents</h2>
<li><A HREF="j2suscep"><strong>J2suscep</strong></A>   </li>
<li><A HREF="Installation"><strong>Installation Instructions</strong></A>   </li>
<li><A HREF="Execution"><strong>Code Execution</strong></A>   </li>
<li><A HREF="contribution"><strong>Contributing</strong></A>   </li>
<h3> EJ_Calc </h3>
<ul>
  <li><A HREF="1a_introduction"><strong>Introduction</strong></A>   </li>
  <li><strong><a href="1b_theory">Theory</a></strong>
  <li><strong><A HREF="1c_basic_details">Basic Details</A></strong></li>
  <li><strong><A HREF="1e_example">Example Input Files</A></strong></li>
  <li><strong><A HREF="1f_output_file_description">Description of Output File</A></strong></li>
  <li><strong><A HREF="1g_file_summary">Summary of Input and Output Files</A></strong></li>
  <li><strong><A HREF="1h_developers">For Developers</A></strong>  </li>
  </ul>
</li>

<h3>Suscep </h3>
<ul>
  <li><A HREF="2a_introduction"><strong>Introduction</strong></A>   </li>
  <li><strong><a href="2b_theory">Theory</a></strong>
  <li><strong><A HREF="2c_basic_details">Basic Details</A></strong></li>
  <li><strong><A HREF="2e_example">Example Input Files</A></strong></li>
  <li><strong><A HREF="2f_output_file_description">Description of Output File</A></strong></li>
  <li><strong><A HREF="2g_limitations">Limitations</A></strong></li>
  <li><strong><A HREF="2h_developers">For Developers</A></strong>  </li>
  </ul>
</li>
<li><A HREF="worked_example"><strong>Worked Example</strong></A>   </li>
<li><a href="refs"><strong>References</a></strong></li><html><head></head><body>
<h1>Summary of Input and Output Files</h1>


<p>A summary of all input and output files is given below:</p>


<table class="tg">
<thead>
  <tr>
    <th class="tg-7btt"> Filename </th>
    <th class="tg-7btt"> File type </th>
    <th class="tg-7btt"> Function </th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0pky"> file1.inp </td>
    <td class="tg-0pky"> Input </td>
    <td class="tg-0pky">   Contains information about total number of paramagnetic centres and the Hamiltonian specifying the relationship between various metal   centres.   </td>
  </tr>
  <tr>
    <td class="tg-0pky"> file2 </td>
    <td class="tg-0pky"> Input </td>
    <td class="tg-0pky">   Contains spin densities on the paramagnetic centres in different spin states and the energy of the corresponding states.   </td>
  </tr>
  <tr>
    <td class="tg-0lax"> file1_form.out/file_spin.out </td>
    <td class="tg-0lax"> Output </td>
    <td class="tg-0lax">   Contains information about the calculated coupling constants.   </td>
  </tr>
  <tr>
    <td class="tg-0lax"> bigg </td>
    <td class="tg-0lax"> Output </td>
    <td class="tg-0lax">   Contains information about sets of equations that yield J-values &gt; 50 cm<sup>-1</sup>.   </td>
  </tr>
  <tr>
    <td class="tg-0lax"> singul </td>
    <td class="tg-0lax"> Output </td>
    <td class="tg-0lax">   Contains information about sets of equations that are singular and is produced by ej_calc_form. This file is required for the proper   execution of the executable ej_calc.    </td>
  </tr>
  <tr>
    <td class="tg-0lax"> nrg_frm_jval_form/nrg_frm_jval </td>
    <td class="tg-0lax"> Output </td>
    <td class="tg-0lax">   Contains the approximate energy of the all possible 2<sup>0.5x</sup> states (x = number of paramagnetic centres).   </td>
  </tr>
</tbody>
</table>

<p></p>
    
</body></html>
<html><head></head><body>
<h1>Description of Output File</h1>

<p>The format of the output file is as follows:</p>


<table class="tg">
<thead>
  <tr>
    <th class="tg-7btt">Output file sections</th>
    <th class="tg-7btt">Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0pky">Number of Magnetic centers:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;6 <br> g value:  1.9980 <br> Temperature range to be used (K):  2.0000 - 300.0000<br>.........<br>.........<br></td>
    <td class="tg-0pky">Details of the input parameters and the Hamiltonian that will be used.</td>
  </tr>
  <tr>
    <td class="tg-0pky">J value (cm-1) with the pairs on which they act<br>.........<br>.........<br></td>
    <td class="tg-0pky">This section specifies the interactions that are accounted for by a given J-value and should reflect the input Hamiltonian.</td>
  </tr>
  <tr>
    <td class="tg-0pky">Temperature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Chi&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Chi*Temp<br>.........<br>.........<br>.........<br></td>
    <td class="tg-0pky">This section gives the χ and χT values at different temperatures.</td>
  </tr>
</tbody>
</table>


<p></p>




    
</body></html>
<html><head></head><body>
<h1>Example Input Files</h1>

<h2>Example 1</h2>
<p>To understand the construction of the input files, let us take the example of a dimeric Mn complex shown in Figure 1 below.
Here the Mn centres are in +IV oxidation state (d<sup>3</sup> spin configuration).
<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs.htm#Libby" class="showTip Libby">Libby</a>]</font></p>

<p align="center"><img src="img/mn2.jpg" alt="Mn2 structure"><br /><b><i>Structure of the {Mn<sup>IV</sup><sub>2</sub>} complex. Colour scheme: Mn (dark blue), C (black), 
N (blue) and O (red). Hydrogen atoms have been removed for clarity.</i></b></p> 
<p>The Mn centres in this complex interact with each other via the oxo groups and the Hamiltonian for this complex can be written as follows:</p>
<center> <font face="Courier New, Courier, monospace" size="-2">&nbsp;&nbsp;&nbsp;&nbsp;</font>H&#770 = -2J&lt;s<sub>1</sub>.s<sub>2</sub>&gt; </center>

<p>This Hamiltonian can be represented using the following input file: </p>
<pre>
magnetic centres
2
No. of J values
1
Hamiltonian
1 2
Hamiltonian Ends

</pre>

<p> The spin operators were obtained using Bader population analysis and the energy was obtained from DFT calculations. The DFT calculations were 
performed using the PBE0 functional in conjunction with the SDDALL basis set having an effective core potential for the Mn atoms and 6-31g (d,p) 
basis set for C, O, N and H atoms. The first two numbers on each line represent the spin density on each of the two Mn centre in a given spin 
configuration and the last number on the line represents the energy of the corresponding state in Hartrees. The following is an example spin operator 
file for the same {Mn<sub>2</sub>} system.</p>
<pre>
2.8580703 -2.8580187 -2.102258638432228E+03
2.9222485 2.9222495 -2.102255438348413E+03
</pre>

<p>As the Mn centres are in +IV oxidation state, formally speaking, they should have three electrons in the d-orbitals. 
This is reflected in the spin density values which are close to 3. Although, Bader spin analysis was used here to determine the spin values,
it is possible to employ a different analysis scheme as well (e.g. Mulliken or Hirshfeld).</p>
<p>The energies given here are in scientific notation but they can be in non-scientific notation too.</p>
<p>These example input files and the resultant output file have also been provided as separate files in the package in the 'examples'directory.</p>
<br>

<h2>Example 2</h2>
This examples corresponds to a {Mn<sub>3</sub>} complex shown in Figure 2.
<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs.htm#Li" class="showTip Li">Li</a>]</font>  </p>

<p align="center"><img src="img/mn3.jpg" alt="Mn3 structure"><br /><b><i>Structure of the {Mn<sup>III</sup><sub>3</sub>} complex. Colour scheme: Mn (dark blue), C (black), 
N (blue) and O (red). Hydrogen atoms have been removed for clarity.</i></b></p> 

<p>The structure of this molecule suggests that the interaction between each pair of Mn centres would be similar and hence, the Hamiltonian for this system can
be written as follows:</p>
<center> <font face="Courier New, Courier, monospace" size="-2">&nbsp;&nbsp;&nbsp;&nbsp;</font>H&#770 = -2J[&lt;s<sub>1</sub>.s<sub>2</sub>&gt; + 
&lt;s<sub>1</sub>.s<sub>3</sub>&gt; + &lt;s<sub>2</sub>.s<sub>3</sub>&gt;] </center>
<p>This Hamiltonian can be represented using the following input file: </p>
<pre>
magnetic centres
3
No. of J values
1
Hamiltonian
1 2
1 3
2 3
Hamiltonian Ends
</pre>
<br>

<p>The spin operators for this complex were obtained using Bader population analysis and the energy was obtained from DFT calculations 
in a similar manner as for the {Mn<sub>2</sub>) complex described in example 1. The input file for the spins is given below: </p>
<pre>
-3.8272159 3.8433176 3.8435824 -2.973075127863509E+03
3.8610277 3.8601792 3.8611721 -2.973073907698517E+03
</pre>

<p> The first three numbers on each line represent the spin density on each of the three Mn centre in a given spin configuration and the last 
number on the line represents the energy of the corresponding state in Hartrees. As the Mn centres are in +III oxidation state, they should have three unpaired electrons  
which is reflected in the spin density values which are close to 3. </p>
<p>These example input files and the resultant output file have also been provided as separate files in the package in the 'examples'directory.</p>





<br>
<h2>Example 3, 4 and 5</h2>

<p>These examples illustrate the use of ej_calc for the calculation of coupling constants in high-nuclearity complexes 
using different oxidation states of a {Mn<sup>III</sup><sub>6</sub>} coordination complex (Figure 3) as a model.
<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs.htm#Tandon" class="showTip Tandon">Tandon</a>]</font>  </p>

<p align="center"><img src="img/mn6.jpg" alt="Mn6 structure"><br /><b><i>Structure of the {Mn<sub>6</sub>} complex. Colour scheme: Mn (dark blue), P (pink), Cl (green), C (black), 
N (blue) and O (red). Hydrogen atoms have been removed for clarity.</i></b></p>
 
<p>The Mn centres in this complex interact with each other via the phosphonate ligands and the central Cl<sup>-</sup> ion. 
Each Mn centre in this complex has 1 trans- and 4 cis- neighbours and one requires 2 J-values to account for these cis- and trans- interactions 
between Mn centres.<font face="Arial, Helvetica, sans-serif" size="-2">[<a href="refs.htm#Tandon" class="showTip Tandon">Tandon</a>]</font>  </p> 


<p>In example 3, four of the Mn centres, Mn1, Mn2, Mn4 and Mn5, are in +IV oxidation state (d<sup>3</sup> spin configuration) while the others are in +III oxidation state. 
This decreases the overall symmetry of the complex. Therefore, the Hamiltonian for this complex can be written as follows:</p>
<center> <font face="Courier New, Courier, monospace" size="-2">&nbsp;&nbsp;&nbsp;&nbsp;</font>H&#770 = -2J<sub>1</sub>[&lt;s<sub>1</sub>.s<sub>2</sub>&gt;]
-2J<sub>2</sub>[&lt;s<sub>1</sub>.s<sub>6</sub>&gt;] -2J<sub>3</sub> [&lt;s<sub>2</sub>.s<sub>6</sub>&gt;]	
-2J<sub>4</sub>[&lt;s<sub>3</sub>.s<sub>6</sub>&gt;] -2J<sub>5</sub>[&lt;s<sub>4</sub>.s<sub>5</sub>&gt;] -2J<sub>6</sub>[&lt;s<sub>1</sub>.s<sub>3</sub>&gt; + 
&lt;s<sub>2</sub>.s<sub>3</sub>&gt;] -2J<sub>7</sub>[&lt;s<sub>1</sub>.s<sub>4</sub>&gt; + 	&lt;s<sub>2</sub>.s<sub>5</sub>&gt;]
-2J<sub>8</sub>[&lt;s<sub>1</sub>.s<sub>5</sub>&gt; + &lt;s<sub>2</sub>.s<sub>4</sub>&gt;] -2J<sub>9</sub>[&lt;s<sub>3</sub>.s<sub>4</sub>&gt;
 + &lt;s<sub>3</sub>.s<sub>5</sub>&gt;] -2J<sub>10</sub>[&lt;s<sub>4</sub>.s<sub>6</sub>&gt; + &lt;s<sub>5</sub>.s<sub>6</sub>&gt;] </center>
<br>

<p>This Hamiltonian can be represented using the following input file: </p>
<pre>
magnetic centres
6
No. of J values
10
Hamiltonian
1 2
****
1 6
****
2 6
****
3 6
****
4 5
****
1 3
2 3
****
1 4
2 5
****
1 5
2 4
****
3 4
3 5
****
4 6
5 6
Hamiltonian Ends

</pre>
<br>

<p>The spin for this complex contains the spin denisty values and energies for this complex in 16 different states:</p>
<pre>
3.0571596 -3.0639469 -3.7464936 2.9855949 -2.9932085 3.7482157 -8.604041138405300E+03
3.0630259 3.0695029 3.7482020 2.9895879 2.9967993 3.7496550 -8.604040107737754E+03
3.0598187 -3.0625074 3.7479051 2.9894173 2.9952609 3.7490933 -8.604040680248647E+03
3.0619779 3.0682985 -3.7456508 2.9868525 2.9946327 3.7487343 -8.604040566607677E+03
-3.0584743 -3.0654666 3.7469864 2.9876366 2.9949252 3.7488185 -8.604040902043736E+03
3.0589748 -3.0636630 -3.7461917 2.9864062 2.9934310 3.7484142 -8.604040874943079E+03
3.0586563 -3.0627242 3.7475793 2.9888559 -2.9897365 3.7489216 -8.604041108898091E+03
3.0608172 3.0671848 -3.7463169 2.9845575 2.9917289 -3.7479541 -8.604040943152329E+03
-3.0597455 -3.0663612 -3.7469399 2.9851125 2.9926246 3.7481106 -8.604040854591891E+03
-3.0600393 -3.0659642 3.7469386 2.9873684 -2.9907475 3.7484449 -8.604041044599930E+03
3.0572042 -3.0640366 -3.7464007 2.9855883 -2.9930794 3.7481997 -8.604041137238783E+03
3.0622106 3.0684403 3.7478383 -2.9811889 2.9966250 3.7491837 -8.604040552284641E+03
3.0595895 -3.0636304 3.7475299 -2.9825245 2.9951876 3.7486818 -8.604040820396753E+03
3.0588778 -3.0640557 3.7470154 2.9869496 2.9926678 -3.7477576 -8.604040921887142E+03
-3.0596644 -3.0668741 3.7464812 2.9853693 2.9919554 -3.7480008 -8.604040857690830E+03
3.0572908 -3.0642910 3.7466734 2.9859264 -2.9933809 -3.7481567 -8.604041205874382E+03
</pre>
<br>

<p>These files along with the resultant output files are provided in the package. </p>

<p>In example 4, the Mn centres Mn1, Mn2 and Mn4 are in +IV oxidation state while the others are in +III oxidation state. The overall symmetry 
of the complex is significantly reduced and capturing the full electronic picture requires the use of the following Hamiltonian:</p>
<center> <font face="Courier New, Courier, monospace" size="-2">&nbsp;&nbsp;&nbsp;&nbsp;</font>H&#770 = -2J<sub>1</sub>[&lt;s<sub>1</sub>.s<sub>2</sub>&gt;]
-2J<sub>2</sub>[&lt;s<sub>1</sub>.s<sub>3</sub>&gt;] -2J<sub>3</sub>[&lt;s<sub>1</sub>.s<sub>4</sub>&gt;] -2J<sub>4</sub> [&lt;s<sub>1</sub>.s<sub>5</sub>&gt;]
-2J<sub>5</sub> [&lt;s<sub>1</sub>.s<sub>6</sub>&gt;] -2J<sub>6</sub> [&lt;s<sub>2</sub>.s<sub>3</sub>&gt;] -2J<sub>7</sub> [&lt;s<sub>2</sub>.s<sub>4</sub>&gt;]
-2J<sub>8</sub> [&lt;s<sub>2</sub>.s<sub>5</sub>&gt;] -2J<sub>9</sub> [&lt;s<sub>2</sub>.s<sub>6</sub>&gt;] -2J<sub>10</sub> [&lt;s<sub>3</sub>.s<sub>4</sub>&gt;]
-2J<sub>11</sub>[&lt;s<sub>3</sub>.s<sub>5</sub>&gt;] -2J<sub>12</sub> [&lt;s<sub>3</sub>.s<sub>6</sub>&gt;] -2J<sub>13</sub> [&lt;s<sub>4</sub>.s<sub>5</sub>&gt;] 
-2J<sub>14</sub> [&lt;s<sub>4</sub>.s<sub>6</sub>&gt;] -2J<sub>15</sub> [&lt;s<sub>5</sub>.s<sub>6</sub>&gt;]</center>
<br>


	
<p> In example 5, all the Mn centres are in +III oxidation state. Only 2 J-values are required to account the exchange interactions between Mn centres
and the Hamiltonian for this complex can thus be written as follows:</p>
<center> <font face="Courier New, Courier, monospace" size="-2">&nbsp;&nbsp;&nbsp;&nbsp;</font>H&#770 = -2J<sub>1</sub> [&lt;s<sub>1</sub>.s<sub>2</sub>&gt; + 	
&lt;s<sub>1</sub>.s<sub>3</sub>&gt; + &lt;s<sub>1</sub>.s<sub>5</sub>&gt; + &lt;s<sub>1</sub>.s<sub>6</sub>&gt; + &lt;s<sub>2</sub>.s<sub>3</sub>&gt; + 
&lt;s<sub>2</sub>.s<sub>4</sub>&gt; + &lt;s<sub>2</sub>.s<sub>6</sub>&gt; + &lt;s<sub>3</sub>.s<sub>4</sub>&gt; + &lt;s<sub>3</sub>.s<sub>5</sub>&gt;
 + &lt;s<sub>4</sub>.s<sub>5</sub>&gt; + &lt;s<sub>4</sub>.s<sub>6</sub>&gt; + &lt;s<sub>5</sub>.s<sub>6</sub>&gt;] -2J<sub>2</sub> [&lt;s<sub>1</sub>.s<sub>4</sub>&gt; + 	
&lt;s<sub>2</sub>.s<sub>5</sub>&gt; + &lt;s<sub>3</sub>.s<sub>6</sub>&gt;]</center>
<br>

The input files containing the Hamiltonian and the spin density data and the resultant output files for examples 4 and 5 are provided in the package.



<p></p>




    
</body></html>
<html><head></head><body>


<h1>Installation Instructions</h1>

<p> For compilation, a FORTRAN compiler, and BLAS and LAPACK libraries are required.
To install the code, specify the compiler, the relevant compiler specific options and 
the path to BLAS and LAPACK libraries in the file 'makefile.include' present in the src 
directory. The relevant options for the GNU and Intel FORTRAN compilers are provided in 
the 'arch' subdirectory. </p>
<p> <b>Note:</b> If you are using the GNU fortran compiler, then ensure that you have version 7.3.1 or newer.</p>
<p>For a fresh installation, use the following command:</p>
<pre>
make
</pre>
<p>This should compile the ej_calc 
and suscep codes and place the generated executables (ej_calc_form, ej_calc_spin, suscep) 
in the 'bin' directory outside the 'src' directory.</p>
<p>For a clean installation, use the following commands:</p>
<pre>
make clean
make
</pre>

<p> Once the installation is complete, either copy/move the contents of the bin folder to a directory in your pathway (the 'PATH' environment variable) or add the bin to the pathway. </p>
<p>In bash, the bin can be added to the pathway by adding the path to the ~/.bashrc or ~/.bash_profile file: </p>
<pre>
export PATH="$PATH:/path_to_J2suscep_directory/bin"
</pre>

<p>To load this change, either restart the terminal or use the 'source' command: </p>
<pre>
source ~/.bashrc
</pre>
<p>or</p>
<pre>
source ~/.bash_profile
</pre>

<p>To make sure that the path has been added to the PATH variable, use the following command and check if the correct path has been added to the PATH variable.</p>
<pre>
echo $PATH
</pre>


<h2>Testing the code</h2>
<p>The 'test' directory contains tests to verify the proper compilation of the code.
To run the basic tests, use the following command:</p>
<pre>
make
</pre>
<p>This will run the basic tests. Depending upon the machine and whether the compiler 
specific optimisations have been included in the compilation of the code, the execution 
of these tests may take from seconds to minutes. For most purposes, these tests are sufficient to 
verify the validity of the code.</p>
<p> To run the comprehensive tests, use the following command:</p>
<pre>
make large
</pre>
<p>This will run the longer tests. Depending upon the machine and whether the compiler 
specific optimisations have been included in the compilation of the code, the execution 
of these tests may take anytime between a few hours to a day. 
<p>To run the short and comprehensive tests together, use the following command:</p>
<pre>
make all
</pre>
<p>One can also test each executable individually by using one of the following commands depending upon which code needs to be tested:</p>
<pre>
make test_ej_calc_form
make test_ej_calc_spin
make test_suscep
</pre>

<p></p>




    
</body></html>