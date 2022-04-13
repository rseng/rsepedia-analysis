# How to contribute to Rainbow

:+1::tada: First off, thanks for taking the time to contribute! :tada::+1:

## Report issues or problems with the software

* **Ensure the bug was not already reported** by searching on GitHub under [Issues](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/issues).

* If you're unable to find an open issue addressing the problem, use the relevant bug report [template](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/issues/new/choose) to create one. 

* If not possible, you can [open a blank issue](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/issues/new). Be sure to include a **title and clear description**, as much relevant information as possible, and a **code sample** or an **executable test case** demonstrating the expected behavior that is not occurring.

## Submitting Changes

* Open a new GitHub pull request with the change.

* Ensure the PR description clearly describes the change. Include the relevant issue number if applicable.

* At lease one reviewer should be added so that the code is looked over by someone else.

Always write a clear log message for your commits, for example:

    $ git commit -m "A brief summary of the commit
    > 
    > A paragraph describing what changed and its impact."
    
 ## Coding conventions
 
 Please follow the [Google Python Style Guide](https://google.github.io/styleguide/pyguide.html) when contributing to the code.
 
 ## Seek Support
 
 Please feel free to send an email to the following address for other matters:
 - 0go0vdp95@mozmail.com
# Rainbow

[![PyPI version](https://badge.fury.io/py/rainbow-optical-flow.svg)](https://pypi.org/project/rainbow-optical-flow/)
![](https://img.shields.io/badge/platform-windows%20%7C%20linux%20%7C%20macOS-lightgrey)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Python package](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/actions/workflows/python-package.yaml/badge.svg?branch=main)](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/actions/workflows/python-package.yaml)
[![GitHub license](https://img.shields.io/github/license/AlphonsG/Rainbow-Optical-Flow-For-ALI)](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/blob/main/LICENSE)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04080/status.svg)](https://doi.org/10.21105/joss.04080)

<p style="text-align:center;">
   <img src="https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/raw/main/misc/images/banner_img.png" alt="centered image" />
</p>

Software for automated Air-Liquid Interface cell culture image analysis using deep optical flow. See <a href="#software-paper">below</a> for more details.

## Table of contents
1. <a href="#installation">Installation</a>
2. <a href="#usage">Usage</a>
3. <a href="#additional-information">Additional Information</a>
4. <a href="#examples">Examples</a>
5. <a href="#community-guidelines">Community Guidelines</a>
6. <a href="#license">License</a>
7. <a href="#software-paper">Software Paper</a>
8. <a href="#our-team">Our Team</a>
9. <a href="#acknowledgements">Acknowledgements</a>

## Installation <a id="installation"></a>

Rainbow can be installed on Linux, Windows & macOS and supports Python 3.8 and above. We recommend installing and running Rainbow within a [virtual environment](https://docs.python.org/3/tutorial/venv.html). Although it is not a requirement, we also recommend installing and running Rainbow on a GPU-enabled system to minimize processing times.

1. Download and install [Python](https://www.python.org/downloads/) (Rainbow was tested using [Python version 3.8.10](https://www.python.org/downloads/release/python-3810/)).

2. Launch the terminal (*Linux* and *macOS* users) or command prompt (*Windows* users). The proceeding commands will be entered into the opened window<sup>1</sup>.

3. (Optional but recommended) Create and activate a virtual environment called 'rainbow-env' in your desired directory:

   ```python -m venv rainbow-env```

   ```. rainbow-env/bin/activate``` (*Linux* and *macOS* users) or ```rainbow-env\Scripts\activate.bat``` (*Windows* users)

   ```python -m pip install -U pip```

4. Install PyTorch by specifying your system configuration using the official [PyTorch get started tool](https://pytorch.org/get-started/locally/) and running the generated command:
   <p style="text-align:center;">
    <img src="https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/raw/main/misc/images/pytorch_get_started.png" width="750" alt="centered image" />
    </p>
   For example, Windows users without a GPU (i.e. CPU only) will run:

   ```pip install torch torchvision torchaudio```

Next, proceed wth either option A or B.
### Option A - Install from PyPI

This is the simplest and fastest way to install Rainbow, recommended for normal users.


5. Install Rainbow:

   ```pip install rainbow-optical-flow```


### Option B - Install from Source

Developers may wish to install Rainbow from source. Please ensure [Git](https://git-scm.com/downloads) and [Git LFS](https://git-lfs.github.com/) are installed before proceeding.

5. Clone this repository into your desired directory:

   ```git clone https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI.git```

6. Navigate into the cloned directory:

    ```cd Rainbow-Optical-Flow-For-ALI```

7. Install Rainbow:

   ```pip install -e .```

8. Finalize the installation by running the following commands:

   ```
   git submodule sync

   git submodule update --init --recursive
   ```

Notes:
  - <sup>1</sup>Confirm the correct python version for Rainbow has been installed using the `python -V` command in the terminal. If this command does not report the correct python version, try using the `python3 -v` command instead. If the second command produces the expected result, replace all `python` and `pip` commands in this guide with `python3` and `pip3`, respectively.

  - The virtual environment can be deactivated using:

      ```deactivate```

  - If Rainbow fails to install on Linux, it may be because `wxpython` could not be built (look for clues in the error messages printed on the terminal e.g. "Running setup.py install for wxpython ... error"). Instead, try installing `wxpython` first by following [these](https://wxpython.org/pages/downloads/) instructions (specifically "Yes, we have Linux Wheels. Sort of.") and then attempt to install Rainbow again via ```pip install rainbow-optical-flow``` (option A) or ```pip install -e .``` (option B).

## Usage <a id="usage"></a>

### Command Line Interface (CLI)

Once installed, Rainbow can be used through a CLI. Run `rainbow --help` or `rainbow -h` (within the `rainbow-env` environment if applicable) for a list of available command arguments and descriptions.

To test Rainbow using an example Air-Liquid Interface cell culture image series, follow the instructions under option B of the <a href="#installation">installation</a> procedure (except for step 7) and run the following commands in the terminal:

```
cd rainbow
rainbow ../examples/example_image_series ../misc/configs/default_config.yaml
```

After processing is finished, a folder containing similar outputs (e.g. a HTML report,  videos, images, CSV files) to those in [this](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/tree/main/examples/example_output/(tif)_191018_HNA-ALI_d14.nd2_-_191018_HNA-ALI_d14.nd2_(series_03)0000_etc) example output folder should be generated in [this](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/tree/main/examples/example_image_series) folder.
### Graphical User Interface (GUI)

Once installed, Rainbow can be be used through a GUI, which can be launched by running the command `rainbow` (within the `rainbow-env` environment if applicable).

To test Rainbow using an example Air-Liquid Interface cell culture image series, follow the instructions under option B of the <a href="#installation">installation</a> procedure (except for step 7) and run the following commands in the terminal::

 ```
 cd rainbow
 rainbow
 ```

 Then, in the GUI that opens, select [this](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/tree/main/examples/example_image_series) folder as the input image series and [this](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/blob/main/misc/configs/default_config.yaml) file as the configuration file in the GUI under 'Required Arguments' and click the 'Start' button. After processing is finished, a folder containing similar outputs (e.g. a HTML report,  videos, images, CSV files) to those in [this](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/tree/main/examples/example_output/(tif)_191018_HNA-ALI_d14.nd2_-_191018_HNA-ALI_d14.nd2_(series_03)0000_etc) example output folder should be generated in [this](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/tree/main/examples/example_image_series) folder.

## Additional Information <a id="additional-information"></a>

### Optical Flow

Rainbow uses a deep learning model called [GMA](https://arxiv.org/abs/2104.02409) to compute the optical flow in an image series. This model can be replaced with any other method for computing optical flow by writing a custom class that implements the [base_model](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/blob/main/rainbow/optical_flow/base_model.py) interface ([gma.py](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/blob/main/rainbow/optical_flow/gma.py) is an example of that).

### Analysis

Rainbow can automatically generate an analysis report after computing the optical flow in an image series. A base report file that can be modified is provided [here](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/blob/main/misc/notebooks/report.ipynb) as a Jupyter notebook. The path of a Jupyter notebook needs to specified in the config for automatic report generation (default provided).

### Scripts

The [scripts](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/raw/main/scripts) folder contains python scripts to enable additional functionality such as the ability to combine reports from multiple experiments into one file for simpler viewing and comparisons. Run `python <script-name>.py --help` in the terminal to view the usage instructions for a script.

### Automated Testing

To perform and check the status of the automated tests locally, run the command `pytest` in the terminal, with Rainbow installed, from the root directory of this repository after cloning.

## Examples <a id="examples"></a>

Examples of some of the data generated by Rainbow can be seen below.

### Raw Image Series (left) and Rainbow Optical Flow Visualisation (Right)

<img src="https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/raw/main/misc/images/raw_vs_flow.gif" width="780" />

### Magnitude Heatmaps (Left) and Quiver Plots (Right) Across Image Series

<img src="https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/raw/main/misc/images/heatmap_quiver_plot.gif" />

### Experimental Methods

Primary tracheobronchial epithelial cells were isolated through trans‐laryngeal, non‐bronchoscopic cytologic brushings via an endotracheal tube from two children (3.3 and 4.1 years), as previously described (Looi et al., 2018; Martinovich et al., 2017). The use of tracheobronchial epithelial cells for these studies were approved by the Human Research Ethics Committees of St John of God Hospital and The University of Western Australia. Cells were imaged (20x objective) on day 2 post air-lift for 2.5 hrs every 8 mins with a Nikon C2+ inverted microscope incubated at 37° C with humidified 95% air/5% CO2 using an Okolab live cell imaging chamber to generate time lapse images of maximally migrating cells as previously described (Park et al., 2015; Mitchel et al., 2020). The example image series provided in this repository contains 20 image frames at 1280 x 1024 px resolution.

#### References

Looi,K. et al. (2018) Effects of human rhinovirus on epithelial barrier integrity and function in children with asthma. Clinical & Experimental Allergy, 48, 513–524.

Martinovich,K.M. et al. (2017) Conditionally reprogrammed primary airway epithelial cells maintain morphology, lineage and disease specific functional characteristics. Scientific Reports, 7, 17971.

Mitchel,J.A. et al. (2020) In primary airway epithelial cells, the unjamming transition is distinct from the epithelial-to-mesenchymal transition. Nature Communications, 11, 5053.
Park,J.-A. et al. (2015) Unjamming and cell shape in the asthmatic airway epithelium. Nature Materials, 14, 1040–1048.

## Community guidelines <a id="community-guidelines"></a>

 Guidelines for third-parties wishing to:

- Contribute to the software
- Report issues or problems with the software
- Seek support

can be found [here](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/blob/main/CONTRIBUTING.md).

## License <a id="license"></a>

[MIT License](https://github.com/AlphonsG/Rainbow-Optical-Flow-For-ALI/blob/main/LICENSE)

## Software Paper <a id="software-paper"></a>

### Title
Rainbow: Automated Air-Liquid Interface Cell Culture Analysis Using Deep Optical Flow

### Access

https://joss.theoj.org/papers/10.21105/joss.04080

## Our Team <a id="our-team"></a>
[Learn more](https://walyanrespiratory.telethonkids.org.au/our-research/focus-areas/artifical-intelligence/)

## Acknowledgements <a id="acknowledgements"></a>

- https://github.com/philferriere/tfoptflow
- https://docs.opencv.org/3.4/d4/dee/tutorial_optical_flow.html
- https://github.com/zacjiang/GMA
---
title: 'Rainbow: Automated Air-Liquid Interface Cell Culture Analysis Using Deep Optical Flow'
tags:
  - python
  - jupyter notebook
  - bioinformatics
  - deep learning
  - biology
  - optical flow
  - image analysis
  - bioimage informatics
authors:
  - name: Alphons Gwatimba^[corresponding author]
    orcid: 0000-0002-9607-2839
    affiliation: "1, 2"
  - name: Joseph Ho
    orcid: 0000-0002-9496-506X
    affiliation: "1, 2"
  - name: Thomas Iosifidis^[co-senior author]
    orcid: 0000-0001-8462-5865
    affiliation: "1, 3, 4"
  - name: Yuliya V. Karpievitch^[co-senior author]
    orcid: 0000-0003-2586-3358
    affiliation: "1, 5"

affiliations:
  - name: Wal-yan Respiratory Research Centre, Telethon Kids Institute, University of Western Australia, Nedlands, 6009, Western Australia, Australia
    index: 1
  - name: School of Computer Science and Software Engineering, University of Western Australia, Nedlands, 6009, Western Australia, Australia
    index: 2
  - name: Centre for Cell Therapy and Regenerative Medicine, School of Medicine, University of Western Australia, Nedlands, 6009, Western Australia, Australia
    index: 3
  - name: School of Population Health, Curtin University, Bentley, 6102, Western Australia, Australia
    index: 4
  - name: School of Molecular Sciences, University of Western Australia, Nedlands, 6009, Western Australia, Australia
    index: 5
date: 19 November 2021
bibliography: paper.bib
---

# Summary

`Rainbow` is a free, open-source and cross-platform Python-based tool for Air-Liquid Interface (ALI) cell culture time-lapse image analysis. `Rainbow` accepts input time-lapse images in standard image formats, such as TIFF, or microscopy file formats, such as ND2, and then computes the optical flow; that is, the apparent motion of individual pixels in an image [@turagaAdvancesVideoBasedHuman2010; @beauchemin1995computation; @ZHAI2021107861], across multiple frames using a deep learning based pretrained optical flow model [@jiangLearningEstimateHidden2021]. `Rainbow` then applies circular data analysis to the pixel-level optical flow information to calculate the average magnitude [0, ∞] (\mu{}m) and direction [0, 360) (°) of motion between adjacent frames to quantify cell motility. Additionally, the variance of the magnitude [0, ∞] (\mu{}m) and direction [0, ∞] (°) of motion between adjacent frames is calculated to quantitatively capture the degree of heterogeneity in cell motility.

For each experiment, a CSV file with the minimum, maximum, mean, standard deviation, and variance of the magnitude and direction of cell movement between adjacent frames in an image sequence is produced. Multiple CSV files from different experiments can be combined into one file that can be analyzed for differences in cell motility across multiple experiments. `Rainbow` also includes a high-resolution and easily readable unified hue/saturation-based visualization scheme for the instantaneous vector field of motion between adjacent frames of an image sequence to qualitatively show cell motion. `Rainbow` can be used through a graphical user interface or command line interface and can generate a HTML report containing output images, videos, publication ready figures, and CSV files detailing cell dynamics (refer to Examples folder on GitHub). Importantly, our software is not limited to ALI culture image analysis, and developers can extend the software’s existing pipeline to other use cases. For example, the optical flow model can be readily substituted with different models, as we utilized the Factory Method creational software design pattern. The data analysis and report generated can be adjusted through interactive Jupyter Notebooks, allowing for a flexible and versatile system. Some of `Rainbow`’s visualizations are shown in \autoref{fig1}.

![`Rainbow` optical flow visualization. **A:** The direction of motion at any position within `Rainbow`-generated optical flow images is measured clockwise from the initial horizontal position of a unit circle (left) and is shown using hue values (right). **B:** The magnitude of motion at any position within optical flow images is shown using saturation values. High saturation (100%) corresponds to high magnitude of motion and low saturation (25%) corresponds to low magnitude of motion. **C, G:** Still frames taken from two separate ALI culture image sequences. **D, H:** Unified visualization of optical flow magnitude and direction between adjacent frames of two ALI culture image sequences using `Rainbow`. The arrow indicates the average direction of motion across the image sequence. The circles indicate three localized vortexes that the cells move around in a swirl-like motion as they change direction. **E, I:** Traditional visualization of optical flow between adjacent frames of two ALI culture image sequences using quiver plots containing vector arrows at every 70 px. **F, J:** Polar plots visualizing motion magnitude (concentric circles; µm) and direction (azimuthal angle; degrees) in the same frame of two ALI culture image sequences. Colour scale indicates the number of points migrating towards given direction. All left positioned subfigures from row 2 onwards (C-F) correspond to the same ALI culture image sequence while right positioned subfigures correspond to a different image sequence (G-J). For complete insight, refer to Examples folder on GitHub containing videos. \label{fig1}](figure_1.png){width=65%}

# Statement of need

Differentiated primary airway epithelial cells cultured at air-liquid interface (ALI) are commonly used to assess airway epithelial function in vitro in health and disease, such as asthma [@chenAirliquidInterfaceCell2019; @looiEffectsHumanRhinovirus2018; @martinovichConditionallyReprogrammedPrimary2017]. The integration of image analysis and ALI cell cultures has provided novel insights into cell dynamics, such as the recently identified unjammed-to-jammed transition of AEC characterized by changes in cell motility [@mitchelPrimaryAirwayEpithelial2020; @parkUnjammingCellShape2015]. In chronic respiratory diseases like asthma, increased cell motility and AEC unjamming have been linked to airway remodeling and disease development [@mitchelPrimaryAirwayEpithelial2020; @parkUnjammingCellShape2015]. However, the image analyses performed in these studies are limited. For example, handcrafted methods from the MATLAB Computer Vision Toolbox that compute optical flow have been used to extract cell motion information from ALI culture images [@mitchelPrimaryAirwayEpithelial2020]. This approach requires licenced software, and handcrafted optical flow estimation methods have been outperformed in terms of accuracy by deep learning methods [@savianOpticalFlowEstimation2020]. Furthermore, commonly used cell motion metrics, such as average cell speed, do not capture all unique aspects of cell motion, such as the heterogeneity of cell migration patterns across time. Cell motion is commonly visualized using vector fields, which are useful but bound by an inverse relationship between resolution and readability [@henkesDenseActiveMatter2020; @nnetuImpactJammingBoundaries2012; @osullivanIrradiationInducesEpithelial2020].

To increase understanding of lung disease mechanisms and development of new treatment options for patients, there is a need for open-source solutions for ALI culture image analyses that can be broadly implemented across cell biology laboratories. To the best of our knowledge, `Rainbow` is the first easy-to-use software tool that performs all the above analyses automatically for efficient utilization by non-programmers. `Rainbow` produces automatic cell motion quantifications, figures, and reports that are all easily transferrable into publications. We anticipate that `Rainbow` will provide cell motion characterization for each experiment and allow for easy comparisons among multiple experiments to uncover cell migration mechanisms previously undetermined in health and disease.

# Acknowledgements

We would like to thank the contribution and assistance of all the respiratory fellows, anaesthetists, nurses, hospital staff at St John of God Hospital, Subiaco (Human Research Ethics Committee study approval #901), and Western Australian Epithelial Research Program (WAERP) members. We would also like to thank the families and children participating in this project. This work was supported by the Wal-Yan Respiratory Research Centre Inspiration Award, Cystic Fibrosis Charitable Endowment Charles Bateman Charitable Trust, Western Australian Department of Health Merit Award, and BHP-Telethon Kids Blue Sky Award. Furthermore, this project relies on high quality open source Python packages: NumPy [@harris2020array], matplotlib [@Hunter:2007], pandas [@mckinney-proc-scipy-2010], [PyYaml](https://pyyaml.org/wiki/PyYAMLDocumentation), [ND2Reader](https://github.com/Open-Science-Tools/nd2reader), [Gooey](https://github.com/chriskiehl/Gooey), [Physt](https://physt.readthedocs.io/en/latest/index.html#), [imutils](https://github.com/PyImageSearch/imutils), [MoviePy](https://zulko.github.io/moviepy/), [natsort](https://github.com/SethMMorton/natsort), [PIMS](http://soft-matter.github.io/pims/v0.5/#), tqdm [@tqdmref], pytest [@pytest], FFmpeg [@tomar2006converting], OpenCV [@opencv_library], Jupyter Notebook [@jupyter], Astropy [@astropy:2013; @astropy:2018], PyTorch [@NEURIPS2019_9015], and scikit-image [@scikit-image].

# References

---
name: Bug report
about: Report issues or problems with the software
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
 - OS: [e.g. Windows]
 - Version [e.g. 10]
 - Hardware [e.g. AMD Ryzen™ 7 5800X, NVIDIA GeForce RTX™ 3090]

**Additional context**
Add any other context about the problem here.
