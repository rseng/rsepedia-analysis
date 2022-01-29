# Change Log

## v3.1.0
- Added support for training with Datasets, generators, and other data types supported by Keras [#211](https://github.com/NLeSC/mcfly/issues/211)
- Added separate classes for each model type, also allowing easier extension by own custom models [#239](https://github.com/NLeSC/mcfly/pull/239)
- Dropped support for Python 3.5 (supported are: Python 3.6, 3.7 and 3.8)
- Extend and update documentation in readthedocs [#250](https://github.com/NLeSC/mcfly/pull/250)
- Fix broken model visualization [#256](https://github.com/NLeSC/mcfly/pull/256)

## v3.0.0
- Add ResNet architecture
- Add InceptionTime architecture
- Tensorflow dependency to 2.0
- Dropped support for Python 2.7 and added support for Python 3.7
- Early_stopping argument for train_models_on_samples() changed to early_stopping_patience
- Fix metric name issue in visualization
- Lower level functions arguments have changed by using keyword arguments dic

## v2.1.0
- Add class weight support

## v2.0.1
- Fix documentation inconsistency

## v2.0.0
- Using Tensorflow.keras instead of Keras
- Using keyword 'accuracy' instead of 'acc' in logs like latest keras versions do

## v1.0.5
- Requirements change (keras<2.3.0)

## v1.0.4
- Fixed Zenodo configuration

## v1.0.3
- Requirements change (six>=1.10.0)
- Fixed Zenodo configuration

## v1.0.2
- Small bug fixes
- Added metric option to model training
- Extended documentation
- Removed redundant dependency on Pandas

## v1.0.1
- Small bug fixes
- Compatible with Keras v2.x (no longer compatible with Keras < v2.0.0)
- Tutorial is moved to separate repository (https://github.com/NLeSC/mcfly-tutorial)

## v1.0.0
First major release.
<p align="left">
  <img src="https://raw.githubusercontent.com/NLeSC/mcfly/master/mcflylogo.png" width="200"/>
</p>

![GitHub Workflow Status](https://img.shields.io/github/workflow/status/NLeSC/mcfly/CI%20Build)
[![Coverage](https://scrutinizer-ci.com/g/NLeSC/mcfly/badges/coverage.png?b=master)](https://scrutinizer-ci.com/g/NLeSC/mcfly/statistics/)
[![PyPI](https://img.shields.io/pypi/v/mcfly.svg)](https://pypi.python.org/pypi/mcfly/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.596127.svg)](https://doi.org/10.5281/zenodo.596127)
[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org:/repo/nlesc/mcfly)
<!-- The first 12 lines are skipped while generating 'long description' (see setup.py)) -->

The goal of mcfly is to ease the use of deep learning technology for time series classification. The advantage of deep learning is that it can handle raw data directly, without the need to compute signal features. Deep learning does not require  expert domain knowledge about the data, and has been shown to be competitive with conventional machine learning techniques. As an example, you can apply mcfly on accelerometer data for activity classification, as shown in [the tutorial](https://github.com/NLeSC/mcfly-tutorial).

If you use mcfly in your research, please cite the following software paper:

D. van Kuppevelt, C. Meijer, F. Huber, A. van der Ploeg, S. Georgievska, V.T. van Hees. _Mcfly: Automated deep learning on time series._
SoftwareX,
Volume 12,
2020.
[doi: 10.1016/j.softx.2020.100548](https://doi.org/10.1016/j.softx.2020.100548)

## Installation
Prerequisites:
- Python 3.6, 3.7 or 3.8
- pip
- Tensorflow 2.0, if pip errors that it can't find it for your python/pip version

Installing all dependencies in separate conda environment:
```sh
conda env create -f environment.yml

# activate this new environment
source activate mcfly
```

To install the package, run in the project directory:

`pip install .`

### Installing on Windows
When installing on Windows, there are a few things to take into consideration. The preferred (in other words: easiest) way to install Keras and mcfly is as follows:
* Use [Anaconda](https://www.anaconda.com/download)
* Install numpy and scipy through the conda package manager (and not with pip)
* To install mcfly, run `pip install mcfly` in the cmd prompt.
* Loading and saving models can give problems on Windows, see https://github.com/NLeSC/mcfly-tutorial/issues/17

## Visualization
We build a tool to visualize the configuration and performance of the models. The tool can be found on http://nlesc.github.io/mcfly/. To run the  model visualization on your own computer, cd to the `html` directory and start up a python web server:

`python -m http.server 8888 &`

Navigate to `http://localhost:8888/` in your browser to open the visualization. For a more elaborate description of the visualization see [user manual](https://mcfly.readthedocs.io/en/latest/user_manual.html).


## User documentation
[User and code documentation](https://mcfly.readthedocs.io).

## Contributing
You are welcome to contribute to the code via pull requests. Please have a look at the [NLeSC guide](https://nlesc.gitbooks.io/guide/content/software/software_overview.html) for guidelines about software development.

We use numpy-style docstrings for code documentation.

#### Necessary steps for making a new release
* Check citation.cff using general DOI for all version (option: create file via 'cffinit')
* Create .zenodo.json file from CITATION.cff (using cffconvert)  
```cffconvert --validate```  
```cffconvert --ignore-suspect-keys --outputformat zenodo --outfile .zenodo.json```
* Set new version number in mcfly/_version.py
* Edit Changelog (based on commits in https://github.com/NLeSC/mcfly/compare/v1.0.1...master)
* Create Github release
* Upload to pypi:  
```python setup.py sdist bdist_wheel```  
```python -m twine upload --repository-url https://upload.pypi.org/legacy/ dist/*```  
(or ```python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*``` to test first)
* Check doi on zenodo
* If the visualization has changed, deploy it to github pages:
```
git subtree push --prefix html origin gh-pages
```

## Licensing
Source code and data of mcfly are licensed under the Apache License, version 2.0.
