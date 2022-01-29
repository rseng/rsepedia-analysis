# [ICML2021] Fold2Seq: A Joint Sequence(1D)-Fold(3D) Embedding-based Generative Model for Protein Design

![Fold2Seq Architecture](/fold2seq1.png)

## Environment file: 
* [environment.yml](environment.yml)

## Data and Feature Generation:
* Go to [data/](data/) and check the README there. 

## How to train the model:
* go to [src/](src/) and run:

`python train.py --data_path $path_to_the_data_dictionary --lr $learning_rate --model_save $path_to_the_saved_model`

## How to generate sequences:
* go to [src/](src/) and run:

`python inference.py --trained_model $path_to_the_trained_model --output $path_to_the_output_file --data_path $path_to_the_data_dictionary`


## Fold2Seq generated structures against natural structures:
![Fold2Seq structures](/fold2seq3.png)
# Fold2Seq Data and Feature Generation

![Fold2Seq Architecture](/data/fold2seq2.png)

## Data
The CATH IDs of protein domains in training, validation and two test sets are in [pdb_lists/](pdb_lists/). 


## Feature Generation:
### Input File:
* In order to generate SSE density features, you need to first provide a file with all input proteins' information.  Each row describes a protein domain. The meaning of each column is:
  * Column1: The path to the PDB
  * Column2: The PDB ID
  * Column3: The chain ID
  * Column4: The starting residue ID
  * Column5: The ending residue ID
* An example of this input file is [example/domain_list.txt](example/domain_list.txt).

### Secondary Structure Assignment:
* Moreover, you need to pre-assign a secondary structure element to each residue. We provide an assignment file ([ss.txt](https://drive.google.com/file/d/1B_9JdT43-l0sVOgBJCdCRAN31tOGX8VA/view?usp=sharing)) obtained from RCSB PDB which contains most of exsiting PDBs. You can first check if your protein is in this file. If not, you can append it following the format in the file.  

### Generating features:
* To generate SSE density features, you can run:

`python fold_feat_gen.py --domain_list example/domain_list.txt  --ss ss.txt --out $path_to_the_output_dictionary`.

* It will generate a  python dictionary containing input information and fold features in `fold_features/`.

