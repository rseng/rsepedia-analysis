# [OpenPhi: An interface to access Philips iSyntax whole slide images for computational pathology](https://gitlab.com/BioimageInformaticsGroup/openphi) [![Tweet](https://img.shields.io/twitter/url/http/shields.io.svg?style=social)](https://twitter.com/intent/tweet?text=OpenPhi:%20An%20interface%20to%20access%20Philips%20iSyntax%20whole%20slide%20images%20for%20computational%20pathology&url=https://gitlab.com/BioimageInformaticsGroup/openphi&hashtags=OpenPhi,digitalpathology,iSyntax,WSI,AI)

[![GitLab package version](https://img.shields.io/badge/Version-13.11.0-green.svg)](https://gitlab.com/BioimageInformaticsGroup/openphi/)
[![Python version](https://img.shields.io/badge/PythonVersion-v3.6.9-green.svg)](https://gitlab.com/BioimageInformaticsGroup/openphi/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://gitlab.com/BioimageInformaticsGroup/openphi/-/blob/master/LICENSE)

<p align="center">
  <img width="650" height="300" src="https://gitlab.com/BioimageInformaticsGroup/openphi/-/raw/media/WSI.png">
</p>


## Table of contents

* [Introduction](#introduction)
* [Requirements](#requirements)
* [Basic usage](#basic-usage)
* [Examples](#examples)
* [Contributing](#contributing)
* [License](#License)
* [Contact](#contact)

## Introduction

Open PatHology Interface is an Application Programming Interface (API) offering easy access to iSyntax whole slide images (WSI) produced by the Philips Ultra Fast Scanner.
The software is based on the following core functionalities:

* Reading label and macro images
* Reading regions of interest
* Reading a whole slide image
* Extracting metadata

The API is extensible and can easily be interfaced to existing applications compatible with other WSI formats.
It affords minimal effort on dealing with the intricacies of the proprietary format and adapts existing vendor-neutral code for iSyntax compatibility.

## Requirements and Installation

This software was developed and tested using SDK v2.0 on Ubuntu 18.04 and Python 3.6.9. It requires installing the Philips Software Development Kit (SDK) preliminarily, available for download from [Philips Open Pathology](https://www.openpathology.philips.com). Other required dependencies are Python libraries: [NumPy](https://numpy.org/) (>=1.18) and [Pillow](https://pillow.readthedocs.io/en/stable/) (>=8.0).

OpenPhi has no specific requirements regarding OS or Python version, however it is subject to the Philips SDK requirements. Please refer to the [Philips Open Pathology](https://www.openpathology.philips.com) portal for more information regarding SDK prerequisitions.

OpenPhi package can be directly installed through the Python Package Index (PyPI) using the following command:

    pip install openphi


## Basic usage

#### OpenPhi object

**class openphi.OpenPhi(filename)**

Parameters:

    filename(str) - the file to open.

---

#### OpenPhi object methods

**read_region(location, level, size)**

Method returns rectangular region of the WSI at a desired resolution level as PIL image.

Parameters:

    location(x, y) - top-left coordinate of the region in pixels with reference to level 0.
    level(int)     - the desired resolution level.
    size(x, y)     - width and height of the region of interest in pixels.

---

**get_thumbnail(size)**

Method returns entire WSI at desired maximum dimensions as PIL image.

Parameters:
    
    size(x, y) - maximum width and height of thumbnail WSI (default: (4000,4000)).

---

**read_wsi(level, bgvalue, channels)**

Method returns the entire WSI at a desired resolution level as PIL image.

Parameters:

    level(int)     - the desired resolution level (default: 4).
    bgvalue (int)  - background color used in the RGB channel format (default: 255).
    channels (str) - "RGB" or RGBA channel format (default: "RGBA").
	    
Note that RGBA includes the alpha channel, which can be 255 or 0 for scanned or non-scanned regions, respectively.
RGB uses bgvalue as background color for non-scanned regions.

---

**close()**

Method closes an image.

---

![OpenPhi Figure.](https://gitlab.com/BioimageInformaticsGroup/openphi/-/raw/media/OpenPhi-fig.png)

#### OpenPhi object attributes

**level_count** - number of resolution levels e.g value of 0 and level_count - 1 represent the highest and the lowest resolution level, respectively.

**dimensions** - tuple of (x, y) representing width and height of the highest resolution level (0).

**level_dimensions** - list of (x, y) tuples representing width and height of each resolution level e.g level_dimensions[n] are width and height dimensions for resolution level n.

**level_downsamples** - list of downsample factors for every resolution level e.g level_downsamples[n] is the downsample factor of level n.

**associated_images** - images associated with the slide e.g label or macro image.

---

#### OpenPhi object properties

Metadata is either extracted as Digital Imaging and Communications in Medicine (DICOM) tags or in the format of OpenSlide generic properties.

**openslide.background-color** - a slide’s background color, if any. It is represented as an RGB hex triplet.

**openslide.bounds-height** - the height of the rectangle bounding the non-empty region of the slide.

**openslide.bounds-width** - the width of the rectangle bounding the non-empty region of the slide.

**openslide.bounds-x** - the X coordinate of the rectangle bounding the non-empty region of the slide.

**openslide.bounds-y** - the Y coordinate of the rectangle bounding the non-empty region of the slide.

**openslide.comment** - a slide’s comment.

**openslide.mpp-x** - the number of microns per pixel in the X dimension of level 0.

**openslide.mpp-y** - the number of microns per pixel in the Y dimension of level 0.

**openslide.objective-power** - a slide’s objective power.

**openslide.quickhash-1** - the “quickhash-1” sum.

**openslide.vendor** - an identification of the vendor.

**DICOM_ACQUISITION_DATETIME** - DICOM acquisition time.

**DICOM_MANUFACTURERS_MODEL_NAME** - DICOM scanner model.

**DICOM_DEVICE_SERIAL_NUMBER** - DICOM scanner serial number.

## Examples

```
from openphi import OpenPhi

# Access an .isyntax image.
im = OpenPhi("myimage.isyntax")

# Get number, dimensions and downsampling factors of resolution levels.
numlevels = im.level_count
dimensions = im.level_dimensions
downsamples = im.level_downsamples

# Query metadata properties.
pixelsize = im.properties['openslide.mpp-x']
scanningtime = im.properties['DICOM_ACQUISITION_DATETIME']

# Access label and macro images.
labelim = im.associated_images['label']
macroim = im.associated_images['macro']

# Read a thumbnail of the entire WSI at a size not exceeding 2000 pixels.
thumbnail = im.get_thumbnail(size=(2000,2000))

# Read the entire WSI at resolution level 3, substituting non-scanned regions with white values.
wsi = im.read_wsi(level=3, bgvalue=255, channels="RGB")

# Read a full-res 512 x 512 tile starting from pixel x=4000, y=2000 in full-resolution coordinates.
tile = im.read_region(location=(4000,2000), level=0, size=(512,512))

# Close an image.
im.close()
```

## Tests

For running the unit tests in _/tests_, a Dockerfile is provided for building a Docker image with the necessary dependencies. Prior to building the image, the following items should be downloaded and placed in the working directory with the Dockerfile:
- _philips-sdk-ubuntu_18_04_: The Philips SDK installation files ([Philips Open Pathology](https://www.openpathology.philips.com) portal).
- _testslide.isyntax_: An .isyntax file for testing (https://zenodo.org/record/5037046/files/testslide.isyntax?download=1).

During the Docker build process, the SDK is installed and the .isyntax file is stored into /app/testslide.isyntax in the image. The image can be built by running e.g.:
 
`docker build -t openphi-test .`

A container can then be started in an interactive session by running e.g.:
 
`docker run -it --name openphi-test openphi-test /bin/bash`
 
`docker exec -it openphi-test /bin/bash`

The unit tests can then be run inside the container:
 
`cd <path_to_this_repository>`
 
`python3 -m tests.test_openphi`

## Contributing

If you would like to further contribute please open an issue or submit a pull request.

## License

This open-source software is licensed under the MIT License. Please see the [LICENSE](LICENSE) file for more details.

## Contact

* **Nita Mulliqi** (mulliqi.nita@gmail.com),
* **Kimmo Kartasalo** (kimmo.kartasalo@gmail.com).

If you find this code useful in your research, please cite:

```
@article{mulliqi2021openphi,
  title={OpenPhi: an interface to access Philips iSyntax whole slide images for computational pathology},
  author={Mulliqi, Nita and Kartasalo, Kimmo and Olsson, Henrik and Ji, Xiaoyi and Egevad, Lars and Eklund, Martin and Ruusuvuori, Pekka},
  journal={Bioinformatics},
  volume={37},
  number={21},
  pages={3995--3997},
  year={2021},
  publisher={Oxford University Press}
}
```

