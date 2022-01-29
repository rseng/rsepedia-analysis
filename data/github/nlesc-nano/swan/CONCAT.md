# Change Log

# 0.7.0 [Unreleased]


# 0.6.0 [21/06/2021]
## New
* Add interface to scikit regressors (#85)
* Add interface to HDF5 to store the training results (#88)
* Add interface to GPyTorch (#90)

## Changed
* Fix prediction functionality (#81)
* Return predicted and expected values scaled back (#90)

# 0.5.0 [04/05/2021]
## Changed
* Fix graph neural network implementation (#59)
* Rename the graph neural network to MPNN (message passing NN)

## New
* Interface to [se3-transformer](https://www.dgl.ai/pages/start.html) (#57)
* Calculate a guess to the molecular coordinates using a force field optimization
* Add bond distance as feature for the GNN using the optimized geometries (#59)
* Add early stopping functionality

# 0.4.0 [02/10/2020]

## Changed
* Removed duplicate CAT functionality (#36)
* Moved the properties computation and filtering to its own repo(#44, #45)

# 0.3.0 [21/08/2020]

## New
* Introduce Pipeline to filter ligands (#13, #26)
* Use [SCScore](https://pubs.acs.org/doi/10.1021/acs.jcim.7b00622)
* Use [Horovod](https://github.com/horovod/horovod) to distribute the training 
* Add [mypy](http://mypy-lang.org/) test

# 0.2.0 [25/02/2020]

## New
* Allow to  train in **GPU**.
* Use [Pytorch-geometric](https://github.com/rusty1s/pytorch_geometric) to create molecular graph convolutional networks.

## Changed
* Replace [deepchem](https://deepchem.io/) with [pytorch](https://pytorch.org)
* Use [rdkit](https://www.rdkit.org) to compute the fingerprints.

# 0.1.0

## New

* Interface to [deepchem](https://deepchem.io/) to generate regression models.
* Interface to [CAT](https://github.com/nlesc-nano/CAT) to compute solvation energies and activity coefficients.
* Train or predict values for a set of smiles using a statistical model.
---
name: Bug report
about: Something doesn't work like I expected.
title: ''
labels: bug
assignees: ''

---

Use your best judgment to provide a useful level of information. Depending on the nature of the issue, consider including, e.g.

- Which version of the software you're running
- Any logging or error output you see
---
name: Help wanted
about: I need some help with the code
title: ''
labels: help wanted
assignees: ''

---
---
name: Feature request
about: I would like a new feature to be included in the library
title: ''
labels: enhancement
assignees: ''

---

Tell us what would be a nice feature to make the **swan** library even better!
### All Submissions:

* [ ] Have you followed the guidelines in our Contributing document?
* [ ] Have you checked to ensure there aren't other open [Pull Requests](https://github.com/SCM-NV/qmflows-namd/pulls) for the same update/change?


### New Feature Submissions:

1. [ ] Does your submission pass tests?
2. [ ] Have you check your code quality (e.g. using [flake8](http://flake8.pycqa.org/en/latest/)) prior to submission?

### Changes to Core Features:

* [ ] Have you added an explanation of what your changes do and why you'd like us to include them?
* [ ] Have you written new tests for your core changes, as applicable?
* [ ] Have you successfully ran tests with your changes locally?
# Machine learning model to predict CDFT properties

The following table shows the **mean rvalue** resulting from the linear regression
between the predicted properties as function of the ground true computed properties. Therefore, the closer
these values are to 1 the better the model is to predict the properties *with
the available amount of training data*.

The mean values where computing with the results of **5** training/validation datasets picked randomly
in a 0.8/0.2 proportion, respectively.

The hyperparameters use to perform the training can be found at [model's hyperparameters](Training_hyperparameters.md).


## Amines
The [amines CDFT dataset](https://github.com/nlesc-nano/swan/blob/main/data/Amines/CDFT/all_amines.csv) has the following properties:
* Number of molecules in dataset: **3109**
* Number of labels: **15**

The [geometries file](https://github.com/nlesc-nano/swan/blob/main/data/Amines/CDFT/all_amines.csv) contains all the optimized geometries for the previous dataset.


|          Property name             | FingerprintFullyConnected | MPNN | SE3Transformer|
|:----------------------------------:|:-------------------------:|:----:|:-------------:|
| Dissocation energy (nucleofuge)    | 0.73   | 0.70   | 0.78 |
| Dissociation energy (electrofuge)  | 0.39	  | 0.32   | 0.39 |
| Electroaccepting power(w+)         | 0.42	  | 0.36   | 0.55 |
| Electrodonating power (w-)         | 0.52	  | 0.54   | 0.66 |
| Electronegativity (chi=-mu)        | 0.54	  | 0.50   | 0.61 |
| Electronic chemical potential (mu) | 0.54	  | 0.53   | 0.59 |
| Electronic chemical potential (mu+)| 0.38	  | 0.25   | 0.34 |
| Electronic chemical potential (mu-)| 0.78	  | 0.81   | 0.85 |
| Electrophilicity index (w=omega)   | 0.51	  | 0.53   | 0.62 |
| Global Dual Descriptor Deltaf+     | 0.33	  | 0.26   | 0.35 |
| Global Dual Descriptor Deltaf-     | 0.35	  | 0.26   | 0.36 |
| Hardness (eta)                     | 0.53	  | 0.53   | 0.56 |
| Hyperhardness (gamma)              | 0.39	  | 0.29   | 0.41 |
| Net Electrophilicity               | 0.47	  | 0.55   | 0.62 |
| Softness (S)                       | 0.42	  | 0.07   | 0.20 |



## Carboxylic acids
The [carboxylic acids CDFT dataset](https://github.com/nlesc-nano/swan/blob/main/data/Carboxylic_acids/CDFT/all_carboxylics.csv) has the following properties:
* Number of molecules in dataset: **11044**
* Number of labels: **15**

The [geometries file](https://github.com/nlesc-nano/swan/blob/main/data/Carboxylic_acids/CDFT/all_geometries_carboxylics.json) contains all the optimized geometries for the previous dataset.


|          Property name             | FingerprintFullyConnected | MPNN | SE3Transformer|
|:----------------------------------:|:-------------------------:|:----:|:-------------:|
| Dissocation energy (nucleofuge)    | 0.85   | 0.85 |  0.87 |
| Dissociation energy (electrofuge)  | 0.75   | 0.83 | 	0.83 |
| Electroaccepting power(w+)         | 0.63	  | 0.64 | 	0.66 |
| Electrodonating power (w-)         | 0.79	  | 0.79 | 	0.84 |
| Electronegativity (chi=-mu)        | 0.81	  | 0.85 | 	0.88 |
| Electronic chemical potential (mu) | 0.81	  | 0.80 | 	0.88 |
| Electronic chemical potential (mu+)| 0.72	  | 0.78 | 	0.77 |
| Electronic chemical potential (mu-)| 0.86	  | 0.86 | 	0.91 |
| Electrophilicity index (w=omega)   | 0.72	  | 0.75 | 	0.79 |
| Global Dual Descriptor Deltaf+     | 0.59	  | 0.59 | 	0.61 |
| Global Dual Descriptor Deltaf-     | 0.59	  | 0.59 | 	0.60 |
| Hardness (eta)                     | 0.80	  | 0.74 | 	0.83 |
| Hyperhardness (gamma)              | 0.67	  | 0.68 | 	0.70 |
| Net Electrophilicity               | 0.72	  | 0.78 | 	0.80 |
| Softness (S)                       | 0.41	  | 0.57 | 	0.44 |

# Hyperparameters used to train the NN

The hyperparameters are taken from the following references:

* [Simple Feedforward NN](https://arxiv.org/abs/1509.09292) (No convolution is performed)
* [Neural Message Passing for Quantum Chemistry](https://arxiv.org/abs/1704.01212)
* [SE(3)-Transformer](https://arxiv.org/abs/2006.10503)

## FullyConnectedFingerPrint
| Parameter       | value |
| :-------------: | :---: |
| hidden_cells    | 100   |
| batch_size      | 100   |
| learning rate   | 5e-4  |


## MPNN

|     Parameter   | value |
| :-------------: | :---: |
| num_labels      | 1     |
| num_iterations  | 3     |
| output_channels | 10    |
| batch_size      | 20    |
| learning rate   | 5e-4  |


## SE3Transformer

| Parameter       | value  |
| :-------------: | :----: |
| num_layers      | 4      |
| num_channels    | 16     |
| num_nlayers     | 0      |
| num_degrees     | 4      |
| div             | 4      |
| pooling         | 'avg'  |
| n_heads         | 1      |
| batch_size      | 32     |
| learning rate   | 1e-3   |
