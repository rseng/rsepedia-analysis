# Code for reproducible research - Fast learning of fast transforms, with guarantees

Code to reproduce experiments in "Fast learning of fast transforms, with guarantees" (Quoc-Tung Le, Léon Zheng, Elisa Riccietti, Rémi Gribonval).
Article is available at: https://hal.archives-ouvertes.fr/hal-03438881/. In the interest of reproducible research, this is exactly the version of the code used to generate the figures in the article, which has also been archived at: https://hal.inria.fr/hal-03552956.

## Getting started
Run the following command line to install the source code:
```
pip install -e .
```
To get an approximation of the Hadamard matrix as a product of butterfly factors using our hierarchical factorization 
method, run the following command:
```
python scripts/hierarchical_factorization.py \
    --num_factors 4 \
    --matrix hadamard \
    --tree balanced \
    --results_path ./results.pkl
```
The script outputs the factorization time, and the approximation error measured as the distance (in Frobenius norm)
between the target matrix and the computed approximation which involves a product of butterfly factors.

To reproduce experiments of the paper, run the scripts in the ``scripts/`` folder. A description of these scripts is 
given below. Examples on how to reproduce the figures of the paper are provided in the Jupyter notebook ``plot_figures.ipynb``.

## Structure of the repository
Source code in ``src/``.
Python scripts are in ``scripts/``.

In ``src/`` folder:
* ``hierarchical_fact.py``: implementation of the hierarchical matrix factorization method.
* ``iterative_grad_fact.py``: implementation of iterative gradient-based factorization method.
* ``tree.py``: implementation of a tree for hierarchical matrix factorization described by a tree
* ``utils.py``: helper methods for experiments

In ``scripts/`` folder:
* ``hierarchical_factorization.py``: experiment to run the hierarchical matrix factorization method on a given (noisy) target matrix
* ``hierarchical_comb_vs_balanced.py``: experiment to compare balanced hierarchical matrix factorization and 
unbalanced hierarchical matrix factorization, in terms of approximation error and computation time. 
The comparison is repeated over several sampling of the noisy matrix.
Use this script to reproduce Figure 4 in the paper.
* ``iterative_vs_hierarchical_dft.py``: experiment to compare hierarchical factorization method and 
iterative gradient-based method, in terms of approximation error and computation time. The comparison is repeated over 
several sampling of the noisy matrix.
Use this script to reproduce Figure 3 in the paper.
* ``grad_factorization.py``: experiment to compare running time and approximation error between gradient based method 
and hierarchical factorization. Use this script to reproduce Figure 2 in the paper.
  
Auxiliary source files and scripts are archived in ``archive`` folders.

The ``data/`` folder contains the MEG matrix in ``.npz`` format (use numpy to load the matrix).

## Butterfly factors and Permutation (BP) model
We refers to the following paper for detailed explanation of  BP model:
https://arxiv.org/abs/1903.05895
