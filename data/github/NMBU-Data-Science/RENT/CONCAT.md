---
title: 'RENT: A Python Package for Repeated Elastic Net Feature Selection'
tags:
  - Python
  - feature selection
authors:
  - name: Anna Jenul
    orcid: 0000-0002-6919-3483
    affiliation: 1
  - name: Stefan Schrunner
    orcid: 0000-0003-1327-4855
    affiliation: 1
  - name: Bao Ngoc Huynh
    orcid: 0000-0001-5210-132X
    affiliation: 2
  - name: Oliver Tomic
    orcid: 0000-0003-1595-9962
    affiliation: 1
affiliations:
 - name: Department of Data Science, Norwegian University of Life Sciences
   index: 1
 - name: Department of Physics, Norwegian University of Life Sciences
   index: 2
date: 10.03.2021
bibliography: paper.bib

---

# Summary
Due to modern data acquisition techniques, the number of generated features in measurement data keeps increasing. This increase can make the analysis with standard machine learning methods difficult because of underdetermined systems where the dimensionality of the feature space (number of features) exceeds the dimensionality of the object space (number of observations). A concrete example of such a situation is data acquisition in the healthcare domain, where the number of patients (observations) suffering from a specific condition may be relatively low, but a lot of measurements (number of features) are generated for each patient to acquire a good understanding of the patient's health. A very common challenge is that not all features in a high dimensional space are equally important for predictive tasks &mdash; many might even be redundant. Feature selection deals with finding the most relevant features of a dataset. With help of appropriate methodology, feature selection can reduce (a) the complexity of and (b) noise in the dataset. More importantly, data interpretation of the model becomes easier with fewer features, which is of great importance within domains such as healthcare. Even though feature selection is a well-established research topic, relatively few approaches are focusing on the stability of the selection. The important question at hand is: can we trust that the selected features are really valid or is their selection very dependent on which observations are included in the data? Providing information on the stability of feature selection is vital, especially in wide data sets where the number of features can be many times higher than the number of observations. Here, the inclusion or exclusion of a few observations can have a high impact on which features may be selected.

# Statement of Need
To get an understanding of which features are important and how stable the selection of each feature in the dataset is, a user-friendly software package is needed for this purpose.
The RENT package, implementing the feature selection method of the same name [@Jenul:2021], provides this information through an easy-to-use interface. The package includes functionalities for binary classification and regression problems. RENT is based on an ensemble of elastic net regularized models, which are trained on randomly, iid subsets of the rows of the full training data. Along with selecting informative features, the method provides information on model performance, selection stability, as well as interpretability. Compared to established feature selection packages available in `R` and Python, such as `Rdimtools` [@Rdimtools:2020] implementing Laplacian and Fisher scores or the scikit-learn feature selection module [@scikit-learn] implementing recursive feature elimination and sequential feature selection, RENT creates a deeper understanding of the data by utilizing information acquired through the ensemble. This aspect is realized through tools for post hoc data analysis, visualization, and feature selection validation provided with the package, along with an efficient and user-friendly implementation of the main methodology.

# Concept and Structure of RENT
At its core, RENT trains $K$ independent elastic net regularized models on distinct subsets of the training dataset. Each subset is generated using the scikit-learn function `train_test_split()` which delivers an iid sample from the full training dataset. The sampling processes of different subsets are mutually independent, with the condition that a single data point can appear at most once in each subset. A data point, however, can appear in multiple subsets. The framework is demonstrated in \autoref{fig:RENT}.

![Summary of RENT method [@Jenul:2021].\label{fig:RENT}](images/RENT_overview.png)

Based on three statistical cutoff criteria $\tau_1$, $\tau_2$ and $\tau_3$, relevant features are selected. While $\tau_1$ counts how often each feature was selected over $K$ models, $\tau_2$ quantifies the stability of the feature weights --- a feature where the $K$ weight signs alternate between positive and negative is less stable than a feature where all weights are of a constant sign. The third criterion $\tau_3$ deploys a Student’s $t$-test to judge whether feature weights are significantly different from zero. The presented implementation builds on an abstract class `RENT_Base` with a general skeleton for feature selection and post hoc analysis. Two inherited classes, `RENT_Classification` and `RENT_Regression`, offer target-specific methods. The constructor of `RENT_Base` initializes the different user-specific parameters such as the dataset, elastic net regularization parameters, or the number of models $K$.
After training, feature selection is conducted by use of the cutoff criteria. Deeper insights are provided by a matrix containing the cutoff criteria values of each feature, as well as a matrix comprising raw model weights of each feature throughout the $K$ elementary model. For initial analysis of the results, the package delivers multiple plotting functions, such as a barplot of $\tau_1$. Additionally, two validation studies are implemented: first, a model based on random feature selection is trained, while second, a model based on randomly permuted labels of the test dataset is obtained. Results of both validation models are compared to a model built with RENT features using Student's $t$-tests as well as empirical densities.

In addition to feature selection, RENT offers a detailed summary of prediction accuracies for the training objects. For each training object, this information can be visualized as histograms of class probabilities for classification problems or histograms of mean absolute errors for regression problems, respectively. For extended analysis,  principal component analysis reveals properties of training objects and their relation to features selected by RENT. For computation and visualization of principal components, RENT uses functionality from the `hoggorm` and `hoggormplot` packages [@Tomic:2019].

# Ongoing Research and Dissemination
The manuscript RENT - Repeated Elastic Net Technique for Feature Selection is currently under review. Further, the method and the package are used in
different master thesis projects at the Norwegian University of Life Sciences, mainly in the field of healthcare data analysis.

# Acknowledgements
We thank Runar Helin for proofreading the documentation.

# References
RENT
====

<img src="/images/RENT_logo.png" width="200"/>

RENT (Repeated Elastic Net Technique) is a feature selection method for binary classification and regression problems. At its core
RENT trains an ensemble of <img src="https://render.githubusercontent.com/render/math?math=K\in\mathbb{N}"> generalized linear models using regularized elastic net to select features. Each model <img src="https://render.githubusercontent.com/render/math?math=k=1:K"> in the ensemble is trained using a randomly, iid sampled subset of rows of the full training data. A single data point can appear at most once in each subset, but may appear in multiple subsets. From these <img src="https://render.githubusercontent.com/render/math?math=K">unique models one can acquire weight distributions for each
feature that contain rich information on the stability of feature selection and from which several adjustable classification criteria may be
defined.

More details are in the original paper published in IEEE Access: [RENT - Repeated Elastic Net Technique for Feature Selection](https://ieeexplore.ieee.org/abstract/document/9606766)

Example
-------

Below are links to Jupyter-notebooks that illustrate how to use RENT for	

* [classification](https://github.com/NMBU-Data-Science/RENT/blob/master/examples/Classification_example.ipynb) 
* [regression](https://github.com/NMBU-Data-Science/RENT/blob/master/examples/Regression_example.ipynb)
* [hyperparameter search](https://github.com/NMBU-Data-Science/RENT/blob/master/examples/Extensive_hyperparameter_search.ipynb)
* [BIC hyperparameter search](https://github.com/NMBU-Data-Science/RENT/blob/master/examples/BIC%20hyperparameter%20selection.ipynb)



Requirements
------------
Make sure that Python 3.5 or higher is installed. A convenient way to install Python and many useful packages for scientific computing is to use the [Anaconda Distribution](https://www.anaconda.com/products/individual)

* numpy >= 1.11.3
* pandas >= 1.2.3
* scikit-learn >= 0.22
* scipy >= 1.5.0
* hoggorm >= 0.13.3
* hoggormplot >= 0.13.2
* matplotlib >= 3.2.2
* seaborn >= 0.10



Installation
------------
To install the package with the pip package manager, run the following command:  
`python3 -m pip install git+https://github.com/NMBU-Data-Science/RENT.git`



Documentation
-------------

Documentation is available at [ReadTheDocs](https://rent.readthedocs.io/en/latest/). It provides detailed explanation of methods and their inputs.


Citing the RENT package
--------------

If you use RENT in a report or scientific publication, we would appreciate citations to the following paper:

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03323/status.svg)](https://doi.org/10.21105/joss.03323)

Jenul et al., (2021). RENT: A Python Package for Repeated Elastic Net Feature Selection. Journal of Open Source Software, 6(63), 3323, https://doi.org/10.21105/joss.03323 

Bibtex entry:

    @article{RENT,
    doi = {10.21105/joss.03323},
    url = {https://doi.org/10.21105/joss.03323},
    year = {2021},
    publisher = {The Open Journal},
    volume = {6},
    number = {63},
    pages = {3323},
    author = {Anna Jenul and Stefan Schrunner and Bao Ngoc Huynh and Oliver Tomic},
    title = {RENT: A Python Package for Repeated Elastic Net Feature Selection},
    journal = {Journal of Open Source Software}
    }


# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
oliver.tomic@nmbu.no.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# Contributing to RENT

Your contribution to RENT is very welcome! For the implementantion of a new feature or bug-fixing, we encourage you to send a Pull Request to https://github.com/NMBU-Data-Science/RENT. Please add a detailed and concise
description of the invented feature or the bug. In case of fixing a bug, include comments about your solution. To improve RENT even more, feel free to send us issues with bugs, you are not sure about. 
Furthermore, you can also contribute to the improvement of the Read the Docs documentation page. We are thankful for any kind of constructive criticism and suggestions.

When you are finished with your implementantions and bug-fixing, please send a Pull Request to https://github.com/NMBU-Data-Science/RENT.

## Developing RENT
RENT can easily be developed on your computer. We recommend installing RENT in a separate environment. If you use conda, the flow could be:

1. Create and activate an environment for RENT:

```
conda create -n envRENT python=3.8
conda activate envRENT
```

2. Clone a copy of RENT from source:

```
git clone https://github.com/NMBU-Data-Science/RENT.git
```

3. Install RENT requirements:

```
cd RENT
pip install -r requirements.txt
pip install -e .
```

4. Make sure that the installation works by running the entire test suite with:

```
tox .
```

## Testing RENT
Step 4 in the previous section runs the entire test suite. Navigating to the **test** folder, you can run a test for a classification and a regression problem separately with 

```
cd test
pytest test_classification.py
pytest test_regression.py
```

## Additional comments
* RENT is written in American English
* We use Semantic Versioning (https://semver.org/)

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
