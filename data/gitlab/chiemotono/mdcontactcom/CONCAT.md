# Overview

## MD Contact Comparison (MDContactCom)
MDContactCom is a tool to identify differences of protein molecular dynamics from two MD simulation trajectories in terms of interresidue contact. The software

- compares residue-residue contact fluctuations of two MD trajectories, 
- quantifies the differences using similarity scores,
- identifies sites that exhibit large differences,
- then visualizes those sites on the protein structure.

## FUNCTIONS
MDContactCom consists of 4 subsystems. Each subsystem can be run individually, or all flows can be run at once.

- Contact_profile_calculator<br>
Input: a trajectory of MD simulation (multi-frame PDB format)<br>
Function: calculating a residue-residue contact map or distance map of each frame and averaging them to calculate a residue-residue contact frequency map (contact profile map) of the trajectory<br>
Output: a contact profile map in CSV format
- Contact_similarity_calculator<br>
Input: two contact profile maps (outputs of the Contact_profile_calculator)<br>
Function: calculating similarity scores, Tanimoto coefficient and Euclidean distance, of contact profile of each residue<br>
Output: a table of similarity scores in CSV format<br>
- Similarity_drawer<br>
Input: two contact profile maps (outputs of the Contact_profile_calculator), a table of similarity scores(outputs of the Contact_similarity_calculator)<br>
Function: creating PDB files to highlight residues with significant difference between two trajectories and their contacts on the 3D structure.<br>
Output: protein structure coordinate files in PDB format
- Similarity_plotter<br>
Input: a table of similarity scores (output of the Contact_similarity_calculator)<br>
Function: showing the contact similarities versus residue number as a line graph<br>
Output: a line graph in png format

# INSTALLATION
## REQUIREMENTS
This software requires the following modules:

- python
  - python 3 (tested on python 3.6.10 :: Anaconda, Inc. ) 
- The following Python modules (version numbers in brackets were used during development/testing):
  - numpy (1.14.3)
  - matplotlib (3.1.1)
  - biopandas (0.2.5)
  - pymol (2.3.5)
  - mdtraj (1.9.3)

## INSTALLATION
Download and place all scripts in “mdcontactcom” on your system. 
The following directory will be created. Set PATH of the directory.

```
mdcontactcom
├ __main__.py
├ mdcontactcom.py
├ contact_profile_calculator.py
├ contact_similarity_calculator.py
├ similarity_drawer.py
├ similarity_plotter.py
├ log.py
└ traj2pdb.py
```

# INPUT AND OUTPUT
## Input files
MDContactCom supports multi-frame PDB format as an input file format, as well as trajectory files from MD software GROMACS, Amber, CHARMM, NAMD, and Desmond.

- MD trajectory file
  - multi-frame PDB file (.pdb)
  - GROMACS (.xtc)
  - Amber (.mdcrd)
  - CHARMM (.dcd)
  - NAMD (.dcd)
  - Desmond (a directory containing all trajectory files “frame~” and “clickme.dtr”)

- Topology files<br>
Topology files are also required when you input trajectory files from MD software.
  - GROMACS (.gro)
  - Amber (.prmtop)
  - CHARMM (.psf)
  - NAMD (.psf)
  - Desmond (-out.cms)

## Output files
With two trajectory files WT.pdb and MT.pdb as input, MdContactCom outputs following files. 
See the files in bold text for main results.

|Directory|file names|description|
|:--|:--|:--|
|WT/||result of Contact_profile_calculator|
||contact_profile_calculator.log|log|
||contact_profile.csv|contact profile map|
||distance_profile.csv|averaged distance map|
||distance_map_[#frame].csv|distance map of each frame|
|MT/|result of Contact_profile_calculator|
||contact_profile_calculator.log|log|
||contact_profile.csv|contact profile map
||distance_profile.csv|averaged distance map
||distance_map_[#frame].csv|distance map of each frame|
|similarity/||result of Contact_similrarity_calculator|
||contact_similarity_calculator.log|log|
||**similarity_table.csv**|table of contact similarity|
||contact_profile_diff.csv|difference between two contact profile maps|
|drawer/||result of Similrarity_drawer|
||similarity_drawer.log|log|
||**[WT]_top[x]pct_[pos, neg, st].pdb** (x=5, 10, 20, 50, 100)|PDB file including similarity information|
||[WT]_ranktop[x]_[pos, neg, st].pdb (x=10, 20, 30, 40)|PDB file including similarity information|
||[WT]_cutoff[x]_[pos, neg, st].pdb (x=0.5, 0.7, 0.9)|PDB file including similarity information|
|plotter/||result of Similrarity_plotter
||similarity_plotter.log|log|
||**similarity_plot.png**|line graph of similarity vs residue|

With two NAMD trajectories namd_wt.dcd + namd_wt.psf and namd_mt.dcd + namd_mt.psf as input, MdContactCom outputs following files. See the files in bold text for main results.

|Directory|file names|description|
|:--|:--|:--|
|convert1/||result of executing traj2pdb.py|
||convert.log|log|
||namd_wt.pdb|trajectory file in multi PDB format|
|convert2/||result of executing traj2pdb.py|
||convert.log|log|
||namd_mt.pdb|trajectory file in multi PDB format|
|namd_wt/||result of Contact_profile_calculator|
||contact_profile_calculator.log|log|
||contact_profile.csv|contact profile map|
||distance_profile.csv|averaged distance map|
||distance_map_[#frame].csv|distance map of each frame|
|namd_mt/||result of Contact_profile_calculator|
||contact_profile_calculator.log|log|
||contact_profile.csv|contact profile map|
||distance_profile.csv|averaged distance map|
||distance_map_[#frame].csv|distance map of each frame|
|similarity/||result of Contact_similrarity_calculator|
||contact_similarity_calculator.log|log|
||**similarity_table.csv**|table of contact similarity|
||contact_profile_diff.csv|difference between two contact profile maps|
|drawer/||result of Similrarity_drawer|
||similarity_drawer.log|log|
||**[WT]_top[x]pct_[pos, neg, st].pdb** (x=5, 10, 20, 50, 100)|PDB file including similarity information|
||[WT]_ranktop[x]_[pos, neg, st].pdb (x=10, 20, 30, 40)|PDB file including similarity information|
||[WT]_cutoff[x]_[pos, neg, st].pdb (x=0.5, 0.7, 0.9)|PDB file including similarity information|
|plotter/||result of Similrarity_plotter|
||similarity_plotter.log|log|
||**similarity_plot.png**|line graph of similarity vs residue|

# RUNNING
## All flows in MdContactCom
To run all flows at once:
```
python mdcontactcom run [trajectory1] \
                        -top1 [topology1] \
                        [trajectory2] \
                        -top2 [topology2] \
                        -r [residues to be calculated] \
                        -c [contact threshold] \
                        -p [#parallelization] \
                        -s [similarity score used in “Similarity_drawer” function] \
                        -m \
                        -t [atoms to be considered in contact calculation]
```

Example: To run all flows at once with two trajectory files WT.pdb and MT.pdb as input using default parameters:
```
python mdcontactcom run WT.pdb MT.pdb
```

Example: To run all flows at once with two trajectories namd_wt.dcd + namd_wt.psf and namd_mt.dcd + namd_mt.psf as input using default parameters:
```
python mdcontactcom run namd_wt.dcd -top1 namd_wt.psf \
                        namd_mt.dcd -top2 namd_mt.psf
```

## Contact_profile_calculator
```
python mdcontactcom profile [trajectory] 
                            -top [topology] \
                            -r [residues to be calculated] \
                            -c [contact threshold] \
                            -p [#parallelization] \
                            -m \
                            -t [atoms to be considered in contact calculation]
```

Example: To run this program with a trajectory files WT.pdb as input using default parameters:
```
python mdcontactcom profile WT.pdb
```

Example:  To run this program with a trajectory namd_wt.dcd + namd_wt.psf as input using default parameters:
```
python mdcontactcom profile namd_wt.dcd -top namd_wt.psf 
```

## Contact_similarity_calculator
```
python mdcontactcom similarity [contact_profile_1] [contact_profile_2]
```
Differences between two contact profile maps are output to a file “contact_profile_diff.csv”.

Example:  To run this program with two trajectory files WT_contact_profile.csv and MT_contact_profile.csv as input using default parameters:
```
python mdcontactcom similarity WT_contact_profile.csv MT_contact_profile.csv
```

## Similarity_drawer
```
python mdcontactcom drawer [contact_profile_diff] \
                           [similarity_table] \
                           [trajecotry] \
                           -top [topology] \
                           -n [residues to be excluded] \
                           -s [similarity score used in “Similarity_drawer” function]
```

Example: To run this program with “similarity_table.csv”, ”contact_profile_diff.csv”, and a trajectory file “WT.pdb” as input using default parameters:
```
python mdcontactcom drawer contact_profile_diff.csv similarity_table.csv WT.pdb
```

Example:  To run this program with “contact_profile_diff.csv”, ”contact_profile_diff.csv”, and a trajectory file “namd.dcd” + “namd.psf” as input using default parameters:
```
python mdcontactcom drawer contact_profile_diff.csv \
                           contact_profile_diff.csv \
                           namd.dcd \
                           -top namd.psf
```

## Similarity_plotter
```
python mdcontactcom plotter [similarity_table] \
                            -b [BIN] \
                            -n [residues to be excluded]
```

Example:  To run this program with an input file “similarity_table.csv” using default parameters:
```
python mdcontactcom plotter similarity_table.csv
```
