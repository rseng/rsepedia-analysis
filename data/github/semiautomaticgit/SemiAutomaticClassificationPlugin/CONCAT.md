# Semi Automatic Classification Plugin (SCP)

![SCP](semiautomaticclassificationplugin.png) The Semi-Automatic Classification Plugin (SCP) is a free open source plugin for QGIS that allows  for the supervised classification of remote sensing images, providing tools for the download, the preprocessing and postprocessing of images.

The overall objective of SCP is to provide a set of intertwined tools for raster processing in order to make an automatic workflow and ease the land cover classification, which could be performed also by people whose main field is not remote sensing.

Search and download is available for ASTER, GOES, Landsat, MODIS, Sentinel-1, Sentinel-2, and Sentinel-3 images. Several algorithms are available for the land cover classification. This plugin requires the installation of GDAL, OGR, Numpy, SciPy, and Matplotlib. Some tools require also the installation of SNAP (ESA Sentinel Application Platform).

## Contributing to the development

These instructions will get you a copy of the plugin up and running on your local machine for development and testing purposes.

You do not need any of these steps if you are just interested in using the plugin. 

## Before contributing

If you find a bug of if you want to add a new feature, create a new issue on GitHub to discuss it with the community. Other developers can provide valuable feedback that can improve and make your proposal or your fix even better. 

### Prerequisites

This SCP is for QGIS 3.x. It uses Python 3.x.

Check that you have installed `scipy`, `matplotlib` and `numpy`. To check:

```bash
python3
import scipy;
>>> import scipy;
>>> print(scipy.__version__);
1.0.0
```
Do the same for `matplotlib` and `numpy`.

If necessary, install the required libraries, using `easy_install` or `pip`.
Use `easy_install3` and `pip3` if you have both Python 2.x and 3.x.

```bash
sudo easy_install3 scipy
```

### Workflow overview

Create a fork of the project, using GitHub interface.

Clone your fork on your local computer. You can do it on the command line with:

```bash
mkdir -p ~/dev
git clone git@github.com:yourgithubusername/SemiAutomaticClassificationPlugin.git
```
Your fork will called `origin`. Check that with: 

```bash
cd SemiAutomaticClassificationPlugin
git remote -v
origin	git@github.com:yourgithubusername/SemiAutomaticClassificationPlugin.git (fetch)
origin	git@github.com:yourgithubusername/SemiAutomaticClassificationPlugin.git (push)
```
Add the original repository with the alias `upstream` with:

```bash
git remote add upstream https://github.com/semiautomaticgit/SemiAutomaticClassificationPlugin
```
#### To contribute

Create a new branch. 
```bash
git checkout -b mycontribution
```

Make your changes. Compile and test your changes (more details about this on the next section).

When you have done your edits, commit your local changes, with something like:

```bash
git commit -m "add some useful enhancement"
```
Push your changes back to your GitHub repository fork with:

```bash
git push origin mycontribution
```
You are now ready to issue your Pull Request. Go to your GitHub repository interface and make your Pull Request online.

### Compile and deploy on your local computer 

```bash
cd ~/dev/SemiAutomaticClassificationPlugin/
make deploy
```
`make deploy` will copy the plugin to your QGIS 3 default profile.

If you want to test your modifications in other computers or to distribute to others for testing purposes, create a zip archive with:

```bash
make package VERSION=mycontribution
```

This will create a new archive `SemiAutomaticClassificationPlugin.zip`.

In QGIS 3 you can install a plugin from the zip archive using the plugin manager interface.

![Install SCP from zip archive](docs/install%20SCP%20from%20zip%20archive.png)

## Test the SCP

Start QGIS 3 and check if the plugin is properly installed. If you are running QGIS in another computer or using another profile, install the plugin from the zip file.

and install the plugin from the zip file created.

### End user test

Test your bug fixes or new features carefully. Make sure you did not break any existing code.

Do some screen captures of the new enhancements to publish if you want to issue a pull request.

## Contributing

If the code is working as you expect, follow the steps already mentioned to issue a pull request.

1. Commit your local changes, with something like:

```bash
git commit -m "add some useful enhancement"
```
Push your changes to your GitHub repository with:

```bash
git push origin mycontribution
```
Go to your GitHub repository interface and make your Pull Request. Please be verbose on your comments.

After doing your Pull Request, make sure you are available to provide feedback to questions and comments to your contribution from other developers. In the absence of any feedback concerning your Pull Request, it will be closed. 

## Authors

* **Luca Congedo** 

See also the list of [contributors](https://github.com/semiautomaticgit/SemiAutomaticClassificationPlugin/graphs/contributors) who participated in this project.

## License

This plugin is distributed under a GNU General Public License version 3. To contribute you must accept this license.
# Semi-Automatic Classification Plugin

![SCP](semiautomaticclassificationplugin.png) The Semi-Automatic Classification Plugin (SCP) is a free open source plugin for QGIS that allows  for the supervised classification of remote sensing images, providing tools for the download, the preprocessing and postprocessing of images.

The overall objective of SCP is to provide a set of intertwined tools for raster processing in order to make an automatic workflow and ease the land cover classification, which could be performed also by people whose main field is not remote sensing.

Search and download is available for ASTER, GOES, Landsat, MODIS, Sentinel-1, Sentinel-2, and Sentinel-3 images. Several algorithms are available for the land cover classification. This plugin requires the installation of GDAL, OGR, Numpy, SciPy, and Matplotlib. Some tools require also the installation of SNAP (ESA Sentinel Application Platform).

## Plugin installation

The SCP is available for QGIS version 3.x.
The SCP is developed with Python 3 and requires the installation of GDAL (OGR), NumPy, SciPy and Matplotlib.

For the installation of QGIS and SCP on different operating systems please follow this [guide](https://semiautomaticclassificationmanual.readthedocs.io/en/latest/installation.html).


## Using the plugin

If you are new to SCP, please follow this [tutorial](https://semiautomaticclassificationmanual.readthedocs.io/en/latest/tutorial_1.html).


## Web site

All the SCP information is available from the [SCP website](https://fromgistors.blogspot.com/p/semi-automatic-classification-plugin.html).

## Documentation

Check [the user manual](https://semiautomaticclassificationmanual.readthedocs.io) or the online [tutorials](https://fromgistors.blogspot.com/search/label/Tutorial) available.

[Videos](https://www.youtube.com/user/fromgistors) are also available.

## Contributing to the development

If you find some issue that you are willing to fix, code contributions are welcome. Please read [the development notes](DEVELOPMENT.md) before contributing. 

## Authors

* **Luca Congedo** 

## License

This plugin is distributed under a GNU General Public License version 3.

## How to cite 
Congedo, Luca, (2021). Semi-Automatic Classification Plugin: A Python tool for the download and processing of remote sensing images in QGIS. Journal of Open Source Software, 6(64), 3172, https://doi.org/10.21105/joss.03172

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03172/status.svg)](https://doi.org/10.21105/joss.03172)

## Code on Zenodo
[![DOI](https://zenodo.org/badge/18038317.svg)](https://zenodo.org/badge/latestdoi/18038317)
