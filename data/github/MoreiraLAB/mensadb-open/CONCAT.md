# MENSAdb

## The features extracted for the database of membrane protein dimer analysis can be replicated through this repository.

- INSTALLATION REQUIREMENTS: First you will need to run `python setup.py` in your terminal to install all the dependencies necessary for feature extraction. Excluded from these dependencies are *PSI-Blast* and *AutoDockTools* that you need to install independently. Additionally the *non-redundant (nr) database* must be downloaded from NCBI (https://ftp.ncbi.nlm.nih.gov/blast/db/).

- Before feature extraction, you should perform a PRE-PROCESSING of the PDB files. For that you need to:
	- Trim non-transmembrane residues
	- Remove heteroatoms
	- Mutate exotic amino acids
	- Model incomplete structures
	- Dimer extraction from the structure Files
	- Add hydrogens

To see additional details in how to perform data pre-processing, please see our review - "Structural Characterization of Membrane Protein Dimers" published in Methods in Molecular Biology - Protein Supersecondary Structures (https://www.springer.com/us/book/9781493991600).

### A - Obtain all the features using a single PDB file as input.

- **run.py** deploys all the below features as well as the needed libraries to attain the output files. It will look for information in the intermediate file **mensadb_fetcher.py**. To attain all the features run:
`python run.py [pdbid] [chains]`
> Example:
`python run.py 1a0t PQ`

### B - Obtain each feature individually using a single PDB file as input.

- **dssp_features.py** extracts the features from a dssp output file. Also requires the corresponding pdb file. To attain the dssp output file use the DSSP executable and run: `dssp -i [pdb_name.pdb] >[output_name.txt]`, in windows, or `mkdssp -i [pdb_name.pdb] > [output_name.txt]`, in UNIX based operating systems. To attain DSSP features, you can run `python dssp_features.py`, obtaining the following:
    - DSSP index
    - Amino acid number
    - Amino acid code
    - Chain
    - Secondary Structure
    - BP
    - ASA
    - NH-->O_1_relidx
    - O-->NH_1_relidx
    - NH-->O_1_energy
    - O-->NH_1_energy
    - TCO
    - KAPPA
    - Alpha
    - Phi
    - Psi
    - X-CA
    - Y-CA
    - Z-CA

- **features_pssm.py** extracts the pssm "jsd" features from psi-blast output file. To retrieve the pssm files needed you will require the psiblast local installation, the non-redundant (nr) database and your input file, with this, run: `psiblast -query [fasta_file.fasta] -evalue 0.001 -num_iterations 3 -db [nr] -outfmt 5 -out pssm_output_name.txt -out_ascii_pssm [output_name.pssm] -num_threads 6"`. Running this step can be very time-consuming, depending on the computer and the protein. To attain PSSM "jsd" features output, you can run: `python features_pssm.py`.

- **process_binana.py** extracts the features from the BINding ANAlyser output file (BINANA - to download go to http://rocce-vm0.ucsd.edu/data/sw/hosted/binana/#download). To attain the BINANA output, you can run: `python binana_1_2_0.py -receptor /path/to/receptor.pdbqt -ligand /path/to/ligand.pdbqt -output_file /path/to/output.pdb`, as stated in the website of this software. To use this command, you will need their binana_1_2_0.py script, as well as the ".pdbqt" input files. To attain the selected features from the BINANA output, you can run: `python process_binana.py`. A single csv will be written for each of the possible features. These features are related to a dimer, specifically.
	- Below 2.5 Angstrom residues
	- Below 4 Angstrom residues
	- Hydrogen Bonds
	- Hydrophobic contacts
	- Pi-Pi bond stack
	- T - stack
	- Cation - Pi interaction
	- Salt-bridges

- **generate_class.py** uses vmd to extract the interfacial and surface classification for each residue. Makes use of 5 other scripts that are located on the "*mensa_class*" folder. To use these scripts is required the installation of python based vmd. This can be done with: `conda install -c conda-forge vmd-python`. The whole code can be run with `generate_outputs(input_pdb).joint_call(autodock, autodock_2)`. Check the path list and replace with your locations. The possible classes are:

	- non-interface and non-surface: 0
	- non-interface and surface: 2
	- interface and surface: 3

**References**

- Preto A.J., Matos-Filipe P., Koukos P.I., Renault P., Sousa S.F., Moreira I.S. (2019) Structural Characterization of Membrane Protein Dimers. In: Kister A. (eds) Protein Supersecondary Structures. Methods in Molecular Biology, vol 1958. Humana Press, New York, NY

**Please cite**

- Matos-Filipe P., Preto A.J., Koukos P.I., Mourão J., Bonvin A.M.J.J., Moreira I.S. MENSADB: A Thorough Structural Analysis of Membrane Protein Dimers. Available in arXiv:1902.02321 (https://arxiv.org/pdf/1902.02321.pdf) 
