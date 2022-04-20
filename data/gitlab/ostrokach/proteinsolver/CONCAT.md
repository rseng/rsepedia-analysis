# ProteinSolver

[![gitlab](https://img.shields.io/badge/GitLab-main-orange?logo=gitlab)](https://gitlab.com/ostrokach/proteinsolver)
[![docs](https://img.shields.io/badge/docs-v0.1.25-blue.svg?logo=gitbook)](https://ostrokach.gitlab.io/proteinsolver/v0.1.25/)
[![poster](https://img.shields.io/static/v1?label=poster&message=html&color=yellow&logo=reveal.js)](https://ostrokach-posters.gitlab.io/2019-12-13-neurips-poster/7ad67cfdf35a4e3e8346e293dc444074/)
[![binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fmybinder%3AhTGKLsjmxRS8xNyHxRJB%40gitlab.com%2Fostrokach%2Fproteinsolver.git/v0.1.25)
[![conda](https://img.shields.io/conda/dn/ostrokach-forge/proteinsolver.svg?logo=conda-forge)](https://anaconda.org/ostrokach-forge/proteinsolver/)
[![pipeline status](https://gitlab.com/ostrokach/proteinsolver/badges/v0.1.25/pipeline.svg)](https://gitlab.com/ostrokach/proteinsolver/commits/v0.1.25/)
[![coverage report](https://gitlab.com/ostrokach/proteinsolver/badges/master/coverage.svg?job=docs)](https://ostrokach.gitlab.io/proteinsolver/v0.1.25/htmlcov/)

## Description

ProteinSolver is a deep neural network which learns to solve (ill-defined) constraint satisfaction problems (CSPs) from training data. It has shown promising results both on a toy problem of learning how to solve Sudoku puzzles and on a real-world problem of designing protein sequences that fold into a predetermined geometric shape.

## Demo notebooks

The following notebooks can be used to explore the basic functionality of `proteinsolver`.

| Notebook name               | MyBinder                                                                                                                                                                                                                        | Description                                                                                                                                                                                |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `20_sudoku_demo.ipynb`      | [![binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fmybinder%3AhTGKLsjmxRS8xNyHxRJB%40gitlab.com%2Fostrokach%2Fproteinsolver.git/v0.1.25?filepath=proteinsolver%2Fnotebooks%2F20_sudoku_demo.ipynb)      | Use a pre-trained network to solve a single Sudoku puzzle.                                                                                                                                 |
| `06_sudoku_analysis.ipynb`  | [![binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fmybinder%3AhTGKLsjmxRS8xNyHxRJB%40gitlab.com%2Fostrokach%2Fproteinsolver.git/v0.1.25?filepath=proteinsolver%2Fnotebooks%2F06_sudoku_analysis.ipynb)  | Evaluate a network trained to solve Sudoku puzzles using the validation<br>and test datasets.<br>_(This notebook is resource-intensive and is best ran on a machine with a GPU)._          |
| `20_protein_demo.ipynb`     | [![binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fmybinder%3AhTGKLsjmxRS8xNyHxRJB%40gitlab.com%2Fostrokach%2Fproteinsolver.git/v0.1.25?filepath=proteinsolver%2Fnotebooks%2F20_protein_demo.ipynb)     | Use a pre-trained network to design sequences for a single protein geometry.                                                                                                               |
| `06_protein_analysis.ipynb` | [![binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fmybinder%3AhTGKLsjmxRS8xNyHxRJB%40gitlab.com%2Fostrokach%2Fproteinsolver.git/v0.1.25?filepath=proteinsolver%2Fnotebooks%2F06_protein_analysis.ipynb) | Evaluate a network trained to reconstruct protein sequences using the<br>validation and test datasets.<br>_(This notebook is resource-intensive and is best ran on a machine with a GPU)._ |

Other notebooks in the `notebooks/` directory show how to perform more extensive validations of the networks and how to train new networks.

## Docker images

Docker images with all required dependencies are provided at: <https://gitlab.com/ostrokach/proteinsolver/container_registry>.

To evaluate a proteinsolver network from a Jupyter notebook, we can run the following:

```bash
docker run -it --rm -p 8000:8000 registry.gitlab.com/ostrokach/proteinsolver:v0.1.25 jupyter notebook --ip 0.0.0.0 --port 8000
```

## Installation

We recommend installing `proteinsolver` into a clean conda environment using the following command:

```bash
conda create -n proteinsolver -c pytorch -c conda-forge -c kimlab -c ostrokach-forge proteinsolver
conda activate proteinsolver
```

## Development

First, use `conda` to install `proteinsolver` into a new conda environment. This will also install all dependencies.

```bash
conda create -n proteinsolver -c pytorch -c conda-forge -c kimlab -c ostrokach-forge proteinsolver
conda activate proteinsolver
```

Second, run `pip install --editable .` inside the root directory of this package. This will force Python to use the development version of our code.

```bash
cd path/to/proteinsolver
pip install --editable .
```

## Pre-trained models

Pre-trained models can be downloaded using `wget` by running the following command _in the root folder of the `proteinsolver` repository_:

```bash
wget -r -nH --cut-dirs 1 --reject "index.html*" "http://models.proteinsolver.org/v0.1/"
```

For an example of how to use a pretrained ProteinSolver models in downstream applications (such as mutation ΔΔG prediction), see the [`elaspic/elaspic2`](https://gitlab.com/elaspic/elaspic2) repository, and in particular the [`src/elaspic2/plugins/proteinsolver`](https://gitlab.com/elaspic/elaspic2/-/tree/master/src/elaspic2/plugins/proteinsolver) module.

## Training and validation datasets

Data used to train and validate the "proteinsolver" network to solve Sudoku puzzles and reconstruct protein sequences can be downloaded from <http://deep-protein-gen.data.proteinsolver.org/>:

```bash
wget -r -nH --reject "index.html*" "http://deep-protein-gen.data.proteinsolver.org/"
```

The generation of the training and validation datasets was carried out in our predecessor project: [`ostrokach/protein-adjacency-net`](https://gitlab.com/ostrokach/protein-adjacency-net).

## Environment variables

- `DATAPKG_DATA_DIR` - Location of training and validation data.

## Acknowledgements

<div align="center">
<img src="docs/_static/acknowledgements.svg" width="45%" />
</div>

## References

- Strokach A, Becerra D, Corbi-Verge C, Perez-Riba A, Kim PM. _Fast and flexible protein design using deep graph neural networks_. Cell Systems (2020); 11: 1–10. doi: [10.1016/j.cels.2020.08.016](https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30327-6)
# ProteinSolver

[![gitlab](https://img.shields.io/badge/GitLab-main-orange?logo=gitlab)](https://gitlab.com/ostrokach/proteinsolver)
[![docs](https://img.shields.io/badge/docs-v0.1.25-blue.svg?logo=gitbook)](https://ostrokach.gitlab.io/proteinsolver/v0.1.25/)
[![poster](https://img.shields.io/static/v1?label=poster&message=html&color=yellow&logo=reveal.js)](https://ostrokach-posters.gitlab.io/2019-12-13-neurips-poster/7ad67cfdf35a4e3e8346e293dc444074/)
[![binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fmybinder%3AhTGKLsjmxRS8xNyHxRJB%40gitlab.com%2Fostrokach%2Fproteinsolver.git/v0.1.25)
[![conda](https://img.shields.io/conda/dn/ostrokach-forge/proteinsolver.svg?logo=conda-forge)](https://anaconda.org/ostrokach-forge/proteinsolver/)
[![pipeline status](https://gitlab.com/ostrokach/proteinsolver/badges/v0.1.25/pipeline.svg)](https://gitlab.com/ostrokach/proteinsolver/commits/v0.1.25/)
[![coverage report](https://gitlab.com/ostrokach/proteinsolver/badges/master/coverage.svg?job=docs)](https://ostrokach.gitlab.io/proteinsolver/v0.1.25/htmlcov/)

## Description

ProteinSolver is a deep neural network which learns to solve (ill-defined) constraint satisfaction problems (CSPs) from training data. It has shown promising results both on a toy problem of learning how to solve Sudoku puzzles and on a real-world problem of designing protein sequences that fold into a predetermined geometric shape.

## Demo notebooks

The following notebooks can be used to explore the basic functionality of `proteinsolver`.

| Notebook name               | MyBinder                                                                                                                                                                                                                        | Description                                                                                                                                                                                |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `20_sudoku_demo.ipynb`      | [![binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fmybinder%3AhTGKLsjmxRS8xNyHxRJB%40gitlab.com%2Fostrokach%2Fproteinsolver.git/v0.1.25?filepath=proteinsolver%2Fnotebooks%2F20_sudoku_demo.ipynb)      | Use a pre-trained network to solve a single Sudoku puzzle.                                                                                                                                 |
| `06_sudoku_analysis.ipynb`  | [![binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fmybinder%3AhTGKLsjmxRS8xNyHxRJB%40gitlab.com%2Fostrokach%2Fproteinsolver.git/v0.1.25?filepath=proteinsolver%2Fnotebooks%2F06_sudoku_analysis.ipynb)  | Evaluate a network trained to solve Sudoku puzzles using the validation<br>and test datasets.<br>_(This notebook is resource-intensive and is best ran on a machine with a GPU)._          |
| `20_protein_demo.ipynb`     | [![binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fmybinder%3AhTGKLsjmxRS8xNyHxRJB%40gitlab.com%2Fostrokach%2Fproteinsolver.git/v0.1.25?filepath=proteinsolver%2Fnotebooks%2F20_protein_demo.ipynb)     | Use a pre-trained network to design sequences for a single protein geometry.                                                                                                               |
| `06_protein_analysis.ipynb` | [![binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fmybinder%3AhTGKLsjmxRS8xNyHxRJB%40gitlab.com%2Fostrokach%2Fproteinsolver.git/v0.1.25?filepath=proteinsolver%2Fnotebooks%2F06_protein_analysis.ipynb) | Evaluate a network trained to reconstruct protein sequences using the<br>validation and test datasets.<br>_(This notebook is resource-intensive and is best ran on a machine with a GPU)._ |

Other notebooks in the `notebooks/` directory show how to perform more extensive validations of the networks and how to train new networks.

## Docker images

Docker images with all required dependencies are provided at: <https://gitlab.com/ostrokach/proteinsolver/container_registry>.

To evaluate a proteinsolver network from a Jupyter notebook, we can run the following:

```bash
docker run -it --rm -p 8000:8000 registry.gitlab.com/ostrokach/proteinsolver:v0.1.25 jupyter notebook --ip 0.0.0.0 --port 8000
```

## Installation

We recommend installing `proteinsolver` into a clean conda environment using the following command:

```bash
conda create -n proteinsolver -c pytorch -c conda-forge -c kimlab -c ostrokach-forge proteinsolver
conda activate proteinsolver
```

## Development

First, use `conda` to install `proteinsolver` into a new conda environment. This will also install all dependencies.

```bash
conda create -n proteinsolver -c pytorch -c conda-forge -c kimlab -c ostrokach-forge proteinsolver
conda activate proteinsolver
```

Second, run `pip install --editable .` inside the root directory of this package. This will force Python to use the development version of our code.

```bash
cd path/to/proteinsolver
pip install --editable .
```

## Pre-trained models

Pre-trained models can be downloaded using `wget` by running the following command _in the root folder of the `proteinsolver` repository_:

```bash
wget -r -nH --cut-dirs 1 --reject "index.html*" "http://models.proteinsolver.org/v0.1/"
```

For an example of how to use a pretrained ProteinSolver models in downstream applications (such as mutation ΔΔG prediction), see the [`elaspic/elaspic2`](https://gitlab.com/elaspic/elaspic2) repository, and in particular the [`src/elaspic2/plugins/proteinsolver`](https://gitlab.com/elaspic/elaspic2/-/tree/master/src/elaspic2/plugins/proteinsolver) module.

## Training and validation datasets

Data used to train and validate the "proteinsolver" network to solve Sudoku puzzles and reconstruct protein sequences can be downloaded from <http://deep-protein-gen.data.proteinsolver.org/>:

```bash
wget -r -nH --reject "index.html*" "http://deep-protein-gen.data.proteinsolver.org/"
```

The generation of the training and validation datasets was carried out in our predecessor project: [`ostrokach/protein-adjacency-net`](https://gitlab.com/ostrokach/protein-adjacency-net).

## Environment variables

- `DATAPKG_DATA_DIR` - Location of training and validation data.

## Acknowledgements

<div align="center">
<img src="docs/_static/acknowledgements.svg" width="45%" />
</div>

## References

- Strokach A, Becerra D, Corbi-Verge C, Perez-Riba A, Kim PM. _Fast and flexible protein design using deep graph neural networks_. Cell Systems (2020); 11: 1–10. doi: [10.1016/j.cels.2020.08.016](https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30327-6)
# Deployment

## Deploying to local machines

### Sudoku web server

```bash
docker run -d --restart unless-stopped -p 8080:8080 \
  --env PORT=8080 \
  --env NOTEBOOK_PATH=proteinsolver/notebooks/30_sudoku_dashboard.ipynb \
  registry.gitlab.com/ostrokach/proteinsolver:v0.1.25
```

### Protein design web server

```bash
docker run -d --restart unless-stopped -p 8080:8080 \
  --env PORT=8080 \
  --env NOTEBOOK_PATH=proteinsolver/notebooks/30_design_dashboard.ipynb \
  --gpus '"device=0"' \
  registry.gitlab.com/ostrokach/proteinsolver:v0.1.25
```

## Deploying to Kubernetes using Knative

### Sudoku web server

Create a Sudoku service configuration file, as given below.

```yaml
# sudoku-service.yaml
apiVersion: serving.knative.dev/v1
kind: Service
metadata:
  name: sudoku
  namespace: default
spec:
  template:
    metadata:
      name: sudoku-v1
      annotations:
        # Knative concurrency-based autoscaling (default).
        autoscaling.knative.dev/class: kpa.autoscaling.knative.dev
        autoscaling.knative.dev/metric: concurrency
        # Disable scale to zero with a minScale of 1.
        autoscaling.knative.dev/minScale: "1"
        # Limit scaling to 10 pods.
        autoscaling.knative.dev/maxScale: "10"
    spec:
      containers:
        - name: user-container
          image: registry.gitlab.com/ostrokach/proteinsolver:v0.1.25
          ports:
            - containerPort: 8080
          env:
            - name: NOTEBOOK_PATH
              value: proteinsolver/notebooks/30_sudoku_dashboard.ipynb
          resources:
            limits:
              cpu: "1.8"
              memory: 12Gi
      timeoutSeconds: 600
      containerConcurrency: 12
```

Apply the Sudoku service configuration file.

```bash
kubectl apply -f sudoku-service.yaml
```

### Protein design web server

Create a protein design service configuration file, as given below.

```yaml
# design-service.yaml
apiVersion: serving.knative.dev/v1
kind: Service
metadata:
  name: design
  namespace: default
spec:
  template:
    metadata:
      name: design-v1
      annotations:
        # Knative concurrency-based autoscaling (default).
        autoscaling.knative.dev/class: kpa.autoscaling.knative.dev
        autoscaling.knative.dev/metric: concurrency
        # Disable scale to zero with a minScale of 1.
        autoscaling.knative.dev/minScale: "1"
        # Limit scaling to 10 pods.
        autoscaling.knative.dev/maxScale: "10"
    spec:
      containers:
        - name: user-container
          image: registry.gitlab.com/ostrokach/proteinsolver:v0.1.25
          ports:
            - containerPort: 8080
          env:
            - name: NOTEBOOK_PATH
              value: proteinsolver/notebooks/30_design_dashboard.ipynb
          resources:
            limits:
              cpu: "7"
              memory: 24Gi
              nvidia.com/gpu: "1"
      timeoutSeconds: 600
      containerConcurrency: 8
```

Apply the protein design service configuration file.

```bash
kubectl apply -f design-service.yaml
```
# Usage

To use `proteinsolver` in a project:

```python
import proteinsolver
```
# Installation

## Stable release

To install `proteinsolver` into a new conda environment, run the following command:

```bash
conda create -n proteinsolver -c conda-forge -c kimlab -c ostrokach-forge proteinsolver
conda activate proteinsolver
```

This is the preferred method to install `proteinsolver`, as it will always install the most recent stable release and will also install all dependencies.

If you don't have [conda] installed, this [Python installation guide] can guide
you through the process.

[conda]: https://conda.io
[Python installation guide]: https://conda.io/docs/user-guide/install/index.html

## From sources

The sources for `proteinsolver` can be downloaded from the [GitLab repo].

You can either clone the public repository:

```bash
git clone git://gitlab.com/ostrokach/proteinsolver
```

Or download the [tarball]:

```bash
curl -OL https://gitlab.com/ostrokach/proteinsolver/repository/master/archive.tar
```

Once you have a copy of the source, you can install it with:

```bash
python setup.py install
```

***Warning: Using `pip` to install `proteinsolver` is not recommended since this method does not install the required dependencies.***

[GitLab repo]: https://gitlab.com/ostrokach/proteinsolver
[tarball]: https://gitlab.com/ostrokach/proteinsolver/repository/master/archive.tar
# Credits

## Development Lead

- Alexey Strokach

## Contributors

None yet. Why not be the first?
# Docker

```bash
docker build .
docker tag {image_id} gcr.io/awesome-dialect-184820/proteinsolver:v0.1.1
docker push
```

```bash
PORT=8080
# NOTEBOOK_PATH=proteinsolver/notebooks/30_sudoku_dashboard.ipynb
NOTEBOOK_PATH=proteinsolver/notebooks/30_design_dashboard.ipynb
docker run -d --restart unless-stopped --gpus device=0 --publish ${PORT}:${PORT} \
  registry.gitlab.com/ostrokach/proteinsolver:v0.1.25 \
  voila --no-browser \
  --MappingKernelManager.cull_interval=60 --MappingKernelManager.cull_idle_timeout=600 \
  --ExecutePreprocessor.timeout=3600 \
  --TagRemovePreprocessor.enabled=True --TagRemovePreprocessor.remove_cell_tags='{"hide"}' \
  --VoilaConfiguration.file_whitelist="['favicon.ico', '1n5uA03-distance-matrix.txt']" \
  --template=mytemplate --Voila.ip=0.0.0.0 --port=${PORT} ${NOTEBOOK_PATH}
```
# Pages

This file creates a `./public` folder containing documentation created for multiple versions (tags) of this repository.

When the repository is public, our job is easy: we simply download the `artifact.zip` file from a publicly-accessible URL (see: [downloading the latest artifacts]). However, when the repository is private, using the above-mentioned URL does not work (see: [gitlab-org/gitlab-ce#22957]). In that case, we resort to using the GitLab API instead.

If [gitlab-org/gitlab-ce#22957] is ever fixed, we would be able to specify
`--header "Private-Token: XXXXX"` or attach `&private_token=XXXXX` to the query string,
and keep using the original URL:

```bash
curl --header "Private-Token: XXXXX" \
    "https://gitlab.com/user/repo/-/jobs/artifacts/ref/download?job=job_name"
```

Good resource: <https://docs.gitlab.com/ee/api/jobs.html#download-the-artifacts-archive>.

<!-- Links -->

[downloading the latest artifacts]: https://docs.gitlab.com/ee/user/project/pipelines/job_artifacts.html#downloading-the-latest-artifacts
[gitlab-org/gitlab-ce#22957]: https://gitlab.com/gitlab-org/gitlab-ce/issues/22957
========
Contents
========

.. toctree::
   :caption: Overview
   :name: mastertoc
   :maxdepth: 2

   readme
   installation
   usage
   deployment
   authors

.. toctree::
   :caption: Modules
   :name: modules
   :maxdepth: 3

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Modules
=======

.. autosummary::
   :toctree: _modules

   proteinsolver
{{ fullname }}
{{ underline }}

.. automodule:: {{ fullname }}

   {% block functions %}
   {% if functions %}
   .. rubric:: Functions

   .. autosummary::
      :toctree: {{ objname }}
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: Classes

   .. autosummary::
      :toctree: {{ objname }}
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: Exceptions

   .. autosummary::
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:
   :no-private-members:
   :no-inherited-members:

   {% block methods %}

   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree: {{ objname }}
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
