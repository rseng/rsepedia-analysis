# Code for reproducible research - "Efficient Identification of Butterfly Sparse Matrix Factorization"

In the interest of reproducible research, we provide code to reproduce experiments in 
["Efficient Identification of Butterfly Sparse Matrix Factorization"](https://arxiv.org/abs/2110.01230) 
(Léon Zheng, Elisa Riccietti, Rémi Gribonval).
The hierarchical algorithm to approximate any matrix by a product of butterfly factors is also implemented in 
the [FAµST 3.25 toolbox](https://faust.inria.fr).

## Getting started

Install the requirements and the source code with the following commands:
```
pip install -e .
pip install -r requirements.txt
```

## Benchmark for hierarchical factorization
To compare the running time of the hierarchical factorization algorithm, run:
```
export USE_PROPACK=1
python benchmark_hierarchical.py --min_num_factors 2 --max_num_factors 15 --repeat 3 --results_path ./results/benchmark_hierarchical.csv
```

## Benchmark for SVD
To measure the running time of different SVD solvers from Scipy 1.8.0, run:
```
export USE_PROPACK=1
python benchmark_svd.py --matrix_type rectangle --repeat 10 --min_J 2 --max_J 8 --results_path ./results/benchmark_svd_rectangle.csv
python benchmark_svd.py --matrix_type square --repeat 10 --min_J 2 --max_J 8 --results_path ./results/benchmark_svd_square.csv
```

## Benchmark for matrix multiplication
To measure the running time of matrix-vector multiplication from Numpy, run:
```
export USE_PROPACK=1
python benchmark_matrix_multiplication.py --min_J 2 --max_J 15 --repeat_matrix 10 --repeat_vector 30 --results_path ./results/benchmark_matrix_multiplication.csv
```

## License
BSD 3-Clause License
