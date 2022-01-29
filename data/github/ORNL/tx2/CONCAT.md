# TX2

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![PyPI version](https://badge.fury.io/py/tx2.svg)](https://badge.fury.io/py/tx2)
[![JOSS status](https://joss.theoj.org/papers/b7c161917e5a31af052a597bf98f0e94/status.svg)](https://joss.theoj.org/papers/b7c161917e5a31af052a597bf98f0e94)

Welcome to TX2! This library is intended to aid in the explorability and explainability of
transformer classification networks, or transformer language models with sequence classification
heads. The basic function of this library is to take a trained transformer and
test/train dataset and produce an ipywidget dashboard as seen in the screenshot below,
which can be displayed in a jupyter notebook or jupyter lab.

![screenshot]( https://raw.githubusercontent.com/ORNL/tx2/master/docs/source/screenshot.png)

NOTE: Currently this library's implementation is partially torch-dependent, and so will
not work with tensorflow/keras models - we hope to address this limitation in the future!

## Installation

You can install this package from pypi:

```bash
pip install tx2
```

NOTE: depending on the environment, it may be better to install some of the dependencies separately before
pip installing tx2, e.g. in conda:
```bash
conda install pytorch-gpu pandas scikit-learn matplotlib ipywidgets "numpy<=1.20" -c conda-forge
```

If you do not have access to a GPU on your machine, install the regular pytorch
package:
```bash
conda install pytorch pandas scikit-learn matplotlib ipywidgets "numpy<=1.20"
```

Note that `pytorch-gpu` can only be found in the `conda-forge` channel.

## Examples

Example jupyter notebooks demonstrating and testing the usage of this library can be
found in the [examples
folder](https://github.com/ORNL/tx2/tree/master/examples).

Note that these notebooks can take a while to run the first time, especially
if a GPU is not in use.

Packages you'll need to install for the notebooks to work (in addition to the
conda installs above):

```bash
pip install tqdm transformers~=4.1.1
```

Running through each full notebook will produce the ipywidget dashboard near the
end.

The tests in this repository do not depend on transformers, so raw library
functionality can be tested by running `pytest` in the project root.

## Documentation

The documentation can be viewed at [https://ornl.github.io/tx2/](https://ornl.github.io/tx2/).

The documentation can also be built from scratch with sphinx as needed.

Install all required dependencies: 
```bash
pip install -r requirements.txt
```

Build documentation:

```bash
cd docs
make html
```

The `docs/build/html` folder will now contain an `index.html`

Two notebooks demonstrating the dashboard and how to use TX2 are included
in the `examples` folder, highlighting the default and custom approaches
as discussed in the Basic Usage page of the documentation.

## Citation

To cite usage of TX2 in a publication, the DOI for this code is [https://doi.org/10.21105/joss.03652](https://doi.org/10.21105/joss.03652)

bibtex:
```
@article{Martindale2021,
  doi = {10.21105/joss.03652},
  url = {https://doi.org/10.21105/joss.03652},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {68},
  pages = {3652},
  author = {Nathan Martindale and Scott L. Stewart},
  title = {TX$^2$: Transformer eXplainability and eXploration},
  journal = {Journal of Open Source Software}
}
```
# Contributing to TX2

Help in improving TX2 is welcome! 

If you find a bug or think of an enhancement/improvement you would like to see,
feel free to fill out an appropriate
[issue](https://github.com/ORNL/tx2/issues/new/choose).

If you have a question, double check that it's not covered in our
[documentation](https://ornl.github.io/tx2).

For questions not answered by the docs or anything else that might not fit into one
of the issue templates, you can start a discussion in the [dicussions
tab](https://github.com/ORNL/tx2/discussions).

You are also welcome to contact the developers directly by emailing us at
tx2-help@ornl.gov.

## Submitting a PR

If you have added a useful feature or fixed a bug, open a new pull request with
the changes.  When submitting a pull request, please describe what the pull 
request is addressing and briefly list any significant changes made. If it's in
regards to a specific issue, please include the issue number. Please check and
follow the formatting conventions below!

## Code Formatting

This project uses the [black code formatter](https://github.com/psf/black).

Any public functions and clases should be clearly documented with 
[sphinx-style docstrings](https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html).
Local documentation can be generated with

```bash
cd docs
make html
```
---
title: 'TX$^2$: Transformer eXplainability and eXploration'
tags:
  - Python
  - explainability
  - natural language processing
  - deep networks
  - transformers
authors:
  - name: Nathan Martindale^[co-first author] ^[corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-5036-5433
    affiliation: 1
  - name: Scott L. Stewart^[co-first author] 
    orcid: 0000-0003-4320-5818
    affiliation: 1
affiliations:
 - name: Oak Ridge National Laboratory
   index: 1
date: 21 December 2021
bibliography: paper.bib
---

# Summary

The Transformer eXplainability and eXploration [@tx2], or TX$^2$ software package, is a library designed for artificial intelligence researchers to better understand the performance of transformer models [@vaswani2017attention] used for sequence classification. The tool is capable of integrating with a trained transformer model and a dataset split into training and testing populations to produce an ipywidget [@ipywidgets] dashboard with a number of visualizations to understand model performance with an emphasis on explainability and interpretability. The TX$^2$ package is primarily intended to integrate into a workflow centered around Jupyter Notebooks [@jupyternotebook], and currently assumes the use of PyTorch [@pytorch] and Hugging Face transformers library [@hf-transformers]. The dashboard includes visualization and data exploration features to aid researchers, including an interactive UMAP embedding graph [@mcinnes2018umap] to understand classification clusters, a word salience map that can be updated as researchers alter textual entries in near real time, a set of tools to understand word frequency and importance based on the clusters in the UMAP embedding graph, and a set of traditional confusion matrix analysis tools. 

# Statement of Need

Transformers, although particularly effective on a wide variety of natural language processing tasks, have the same challenge of many deep network approaches in that it is difficult to glean insight into why certain classification decisions are made [@aken2020visbert]. Various works have explored the value of analyzing the attention layers in order to provide explainability in the output of a transformer network [@vig2019multiscale]. However, analyzing attention alone can be insufficient when attempting to gain broader insight into why a transformer is performing a certain way with a specific dataset [@jain2019attention]. TX$^2$ aims to address this challenge by providing a model developer with a number of tools to explore why a certain transformer performs in a certain way for a specific dataset. This tool can help a developer determine, among other things, whether or not a specific transformer has gained a generalized understanding of the semantic meaning behind textual entries in a specific dataset. It can also help with studying the impact of language distribution shifts over time on transformer sequence classification performance.

Existing tools, such as Google PAIR's Language Interpretability Tool [@tenney2020lit], also provide a platform to use multiple visualizations to study transformer model performance. TX$^2$ differs from these tools with its emphasis on cluster analysis and easier customization of both the model interaction and dashboard itself within a Jupyter Notebook. The close integration with Jupyter Notebook is advantageous for those researchers who already rely heavily on the tools within the Jupyter ecosystem. Like the Language Interpretability Tool, TX$^2$ offers a projection map with all of the data points; however it goes further in breaking down the visual clusters and providing separate visualizations for understanding the language per cluster. Additionally, the TX$^2$ design promotes easy modification or customization depending on the researcher's needs, as researchers can completely change the presentation order of plots within the ipywidget and even add additional visualizations if desired. 

# Features

The primary visualization for the widget is a UMAP embedding graph that projects the multidimensional sequence embedding space into 2D clusters. This plots multiple controls that can be used to understand how the sequence classifier is working, including the ability to show or hide training data, highlight certain keywords, and focus on misclassifications. Below the UMAP plot, the dashboard includes a set of tools for exploring textual data including a word salience map that shows information on specific train or test data entries. The salience map serves as a proxy for word importance and is computed by recalculating the soft classifications of a particular entry in the corpus multiple times with each word individually removed. The background coloring in the map indicates the degree of impact word removal has on the classification result, with a darker background highlight corresponding to greater importance. The dashboard also includes a text entry box that is prepopulated with the text from the entry shown in the salience map. The user can use this text box to explore the impact of word addition or removal by modifying the entry. The change is reflected both in the salience map plot as well as with a change in the data point in the UMAP embedding graph.  

The dashboard also includes a set of visual clustering analysis tools. Any clustering algorithm from sklearn's [@scikit-learn] clustering module can be used to assign clusters to the data once it is projected into the UMAP embedding. The dashboard displays cluster labels, along with inter-cluster word frequency, and each words’ importance on the classification result. The salience scores for each word are calculated in aggregate for each cluster, again by iterating with the classifier while individual words are removed. There are also some sampling buttons that allow for a data example to be randomly pulled from a specific cluster so that it can be examined by the entry-specific salience map tool. Finally, it is also possible to output traditional confusion matrices as well as various evaluation scores (e.g., f1-score, accuracy, precision) as part of the dashboard.

# Integration

TX$^2$ includes two main classes: a wrapper class and a dashboard class. The wrapper class wraps around the transformer/classification model and acts as an interface between the dashboard and the transformer. The wrapper is in charge of computing and caching all the necessary data for the dashboard visualizations. The dashboard class is responsible for setting up and rendering the widget layout and handling dashboard interactivity. The flow of interactions between the TX$^2$ library and a Jupyter Notebook can be seen in \autoref{fig:interactivity}.

![Flow of interactions between a Jupyter Notebook and the tx^2 library.\label{fig:interactivity}](interaction_flow.png)

The wrapper communicates with the transformer through a set of four functions as seen in \autoref{fig:integration}. These functions include an embedding function that returns a single sequence of embeddings for each input text, a classification function that returns the predicted output class for each input text, a soft classification function that returns some output value for each class for each input text, and an encoding function that converts the text into model inputs. 

![Example of integrating a transformer with the tx^2 wrapper.\label{fig:integration}](example_interaction.png)

The default implementation for TX$^2$ assumes a huggingface pretrained model. If this use case fits the purposes of the user, they can use the default implementations for these functions. Otherwise, the user will need to redefine the functions to handle their use case while ensuring that the new functions return the necessary data and the correct format.

# Audience

The target audience for the TX$^2$ tool are machine learning or artificial intelligence researchers focused on natural language processing with transformers, and who are comfortable operating within the Jupyter ecosystem for demonstration or exploration. This open-source software is licensed under a BSD-3 clause license, is registered on [DOE Code](https://doi.org/10.11578/dc.20210129.1), and is available on [GitHub](https://github.com/ORNL/tx2). The package is also pip installable with ``pip install tx2`` with Sphinx [@sphinx] built [documentation](https://ornl.github.io/tx2/index.html). Finally, linting for this project is performed using black [@black-linter].

# Acknowledgements

The authors would like to acknowledge the US Department of Energy, National Nuclear Security Administration's Office of Defense Nuclear Nonproliferation Research and Development (NA-22) for supporting this work. 

This manuscript has been authored by UT-Battelle, LLC, under contract DE-AC05-00OR22725 with the US Department of Energy (DOE). The US government retains and the publisher, by accepting the article for publication, acknowledges that the US government retains a nonexclusive, paid-up, irrevocable, worldwide license to publish or reproduce the published form of this manuscript, or allow others to do so, for US government purposes. DOE will provide public access to these results of federally sponsored research in accordance with the DOE Public Access Plan (http://energy.gov/downloads/doe-public-access-plan).

# References
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
