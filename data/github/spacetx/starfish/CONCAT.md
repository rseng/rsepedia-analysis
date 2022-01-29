## [0.2.2] - 2021-04-29
- Updates requirements
- Updates to documentation
- Import match_histograms from skimage.exposure
- Add necessary coords for IntensityTable when using Nearest Neighbors strategy (#1928)
- Fix localmaxpeakfinder spot_props filter (#1839)

## [0.2.1] - 2020-06-08
- Bump napari to 0.3.4 (#1889)
- fix how spot_ids are handled by build_traces_sequential and Label._assign() (#1872)
- reorganized examples gallery and made clarifications to example pipelines and formatting (#1880)
- added image registration tutorial (#1874)
- Add assigning spots to cells docs (#1832)
- Add a Quick Start tutorial (#1869)
- Update starfish installation guide for JOSS submission (#1868)
- Add image segmentation docs (#1821)
- Changing return value of PixelDecoding to DecodedIntensityTable (#1823)
- Ensure that LocalMaxPeakFinder works in 3D (#1822)
- Deprecate is_volume parameter with trackpy (#1820)
- Fix on-demand calculation of BinaryMaskCollection's regionprops (#1819)
- Remove workaround for non-3D images (#1808)
- improve from_code_array validation (#1806)
- Add group_by for tilefetcher-based ImageStack construction (#1796)

## [0.2.0] - 2020-01-31
- Add level_method to the clip filters. (#1758)
- adding method to use installed ilastik instance (#1740)
- Create a TileFetcher-based constructor for ImageStack (#1737)
- adding mouse v human example to starfish.data (#1741)
- adding method to binary mask collection that imports labeled images from external sources like ilastik (#1731)
- Remove starfish.types.Clip (#1729)
- Move watershed segmentation from morphology.Binarize to morphology.Segment (#1720)
- Link to the available datasets in "loading data" section (#1722)
- Document workaround for python3.8 (#1705)
- Wrap skimage's watershed (#1700)
- Add 3D support to target assignment. (#1699)
- Pipeline component and implementation for merging BinaryMaskCollections (#1692)
- Mechanism to reduce multiple masks into one (#1684)

## [0.1.10] - 2019-12-13
- Bump slicedimage to 4.1.1 (#1697)
- Make map/reduce APIs more intuitive (#1686)
- updates roadmap to reflect 2020H1 plans
- adding aws scaling vignette (#1638)
- Use thresholded binarize and mask filtering in existing watershed code. (#1671)
- adding spot ids to pixel results (#1687)
- Implement Labeling algorithms (#1680)
- Thresholded binarize conversion algorithm (#1651)
- Area filter for binary masks (#1673)
- Fix stain generation in watershed (#1670)
- Use the new levels module. (#1669)
- Linear image leveling (#1666)
- add axis labels to display() (#1682)
- Clip method for Richardson Lucy (#1668)
- Filters for mask collections (#1659)
- Provide an apply method to binary mask collections. (#1655)
- adding convience method for slicing codebook data (#1626)
- Fix display tests and code (#1664)
- Additional builders for BinaryMaskCollection (#1637)
- Methods for uncropping binary masks. (#1647)
- Improve coordinate handling code for BinaryMaskCollection and LabelImage (#1632)

## [0.1.9] - 2019-11-18
- Create an ArrayLike type (#1649)
- Verify that binary masks can be generated from empty label images (#1634)
- Add a morphology package to hold BinaryMaskCollection, LabelImage, and their respective operators (#1631)
- fixing travis (#1648)
- Support multiple codewords for the same target (#1646)
- Update data model for BinaryMaskCollection (#1628)
- Test for Codebook.to_json / open_json (#1645)
- Simplify Dockerfile (#1642)
- Switch to version exclusion for scikit-image workaround (#1629)
- Clean up binary mask (#1622)
- adding an extras feild to SpotFindingResults (#1615)
- deleting Decode and Detect modules in lieu of spot finding refactor (#1598)
- Fix install issues (#1641)
- Upgrade to slicedimage 4.1.0 (#1639)
- Update vocabulary for LabelImage I/O operations. (#1630)
- Add a label image data type (#1619)
- Remove deprecated code (#1621)
- fixing bug with codebook.to_json (#1625)
- Don't fill a new ImageStack with NaN (#1609)
- Rename SegmenationMaskCollection to BinaryMaskCollection (#1611)
- Remove hack to force anonymous memory mapping on osx (#1618)

## [0.1.8] - 2019-10-18
- Logging improvements (#1617)
- Make regionprops available per mask (#1610)
- Don't use mypy 0.740 (#1616)
- changing test code to use new spot finding modules (#1597)
- refactoring allen smFish with new spot finding (#1593)
- clean up max projection (#1379)
- Use masked fill to produce labeled images (#1582)
- Replace most instances of starfish.image.Filter.Reduce with imagestack.reduce (#1548)
- implementing starMap spot finding refactor (#1592)
- Add __slots__ to classes that subclass xr.DataArray (#1607)
- Convert SegmentationMaskCollection to a dict-like object (#1579)
- Test case for multiprocessing + imagestack (#1589)
- Masked fill method (#1581)
- Add map/reduce methods to ImageStack (#1539)
- Unify FunctionSource in Map and Reduce (#1540)

## [0.1.7] - 2019-10-09
- ISS refactored with new spot finding path (#1518)
- Fix bugs in per-round-max-decoder (#1602)
- Fix dimension ordering on Codebook and IntensityTable (#1600)
- provanance logging refactor and support for SpotFindingResults (#1517)
- napari 0.2.0 release (#1599)
- starfish.display: unpin napari version, add tests, view masks separately (#1570)
- adding coordinate support to SpotFindingResults (#1516)
- adding new SpotFindingResults data structure and new packages (#1515)

## [0.1.6] - 2019-09-18
- Switch to python multithreading (#1544)
- Don't waste memory/compute in preserve_float_range (#1545)
- Get rid of shared state for LocalMaxPeakFinder (#1541)
- map filter (#1520)
- funcs passed to apply and transform can use positional arguments (#1519)
- import SegmentationMaskCollection in main starfish (#1527)
- Enable Windows builds on master (#1538)
- Throw a warning when the data size is unusual. (#1525)


## [0.1.5] - 2019-08-12
- Update the documentation for data formatters (#1476)
- add ability to convert segmentation masks to a label image
- If in_place=True, we should return None (#1473)
- downgrade pyparsing (#1467)
- fixes unicode in issue template (#1464)
- Adds issue templates (#1460)
- Updating requirements. (#1461)
- Bump to slicedimage 4.0.1 (#1458)
- on-demand loading of data. (#1456)
- Get rid of the check-requirements cron job. (#1448)
- Fixing travis build  (#1457)
- removing duplicate file (#1455)
- Remove Cli (#1444)


## [0.1.4] - 2019-07-16
- Update in-place experiment writing to use the new WriterContract API in slicedimage 4.0.0 (#1447)
- data set formatter with fixed filenames (#1421)

## [0.1.3] - 2019-07-09
- Instantiate the multiprocessing pool using `with` (#1436)
- Slight optimization of pixel decoding  (#1412)
- [easy] point starfish.data.osmFISH() to new dataset (#1425)
- [easy] Warn about the deprecation of the MaxProject filter (#1390)

## [0.1.2] - 2019-06-19
- Refactor reduce to take an optional module and only a function name. (#1386)
- Codify the expectation that in-place experiment construction does not rely on TileFetcher data (#1389)
- Warn and return empty SpotAttributes when PixelDecoding finds 0 spots (#1400)
- updating data.merfish link to full dataset (#1406)
- Rename tile_coordinates to tile_identifier (#1401)
- Support for irregular images in the builder (#1382)
- Fix how we structure the run notebook rules. (#1384)
- updated loading data docs and added image of napari viewer (#1387)
- Format complete ISS experiment and expose in starfish.data (#1316)
- Add concatenate method for ExpressionMatrix (#1381)
- Add TransformsList __repr__ (#1380)
- Fix 3d smFISH notebook as well. (#1385)
- Add custom clip Filter classes (#1376)
- Fix smFISH notebook. (#1383)
- Add Filter.Reduce (general dimension reduction for ImageStack) (#1342)
- Handle denormalized numbers when normalizing intensities/codebooks (#1371)
- TileFetcher formats complete 496 fov MERFISH dataset (#1341)
- Refactor fov.getImage() to fov.getImages() (#1346)
- Add the ability to write labeled experiments (#1374)
- Add inplace TileFetcher module back to public builder API (#1375)
- Always create Z coordinates, even on 4D datasets. (#1358)
- Create an all-purpose ImageStack factory (#1348)
- Remove physical_coordinate_calculator.py (#1352)
- ImageStack parsers should provide coordinates as an array (#1351)
- bump to slicedimage 3.1.1 (#1343)
- Creating a standard starfish.wdl that can be run with any recipe file  (#1364)

## [0.1.1] - 2019-05-16
- [Easy] Fixing spot detection for labeled axes (#1347)
- Schema versioning (#1278)
- Add a missing parameter to wrapper for trackpy local max peak finder (#1300)
- Fix physical coordinate calculator (#1350)
- Fix spot detection for labeled data. (#1349)
- Adding back ability to crop on fov.get_image() (#1329)
- RFC: Base calling filter for in situ sequencing (#1281)
- Preserve subpixel offsets from spot detection (#1330)

## [0.1.0] - 2019-04-23
- public/private separation (#1244)
- Recipe and recipe execution (#1192)
- 3d smFISH notebook (#1238)
- SeqFISH notebook (#1239)
- Adding windows install instructions (#1227)
- vectorize labeling spot lookups (#1215)
- vectorize imagestack -> intensity_table coordinate transfer (#1212)
- Fix the restoration of non-indexed axes. (#1189)
- Allow for intensity tables with labeled axes (#1181)
- ImageStack select on Physical Coordinates (#1147)
- fixing Clip.SCALE_BY_IMAGE (#1193)
- Update BaristaSeq text, fix LinearUnmixing (#1188)
- Update STARmap notebook for SpaceJam (#1199)
- replace label images with segmentation masks (#1135)
- BaristaSeq + Plot tools update (#1171)
- Intensity Table Concat Processing (#1118)

## [0.0.36] - 2019-04-10
- Update strict requirements (#1142)
- High level goal: detect spots should accept imagestacks and not numpy arrays. (#1143)
- Remove cropping from PixelSpotDetector, (#1120)
- Add LocalSearchBlobDetector to support BaristaSeq, SeqFISH, STARmap (#1074)
- Indirect File click types (#1124)
- Move the registration tests next to their sources. (#1134)
- Test to verify that inplace experiment construction works. (#1131)
- Additional support code for building experiments in-place. (#1127)

## [0.0.35] - 2019-04-03
- Transfer physical Coords to Expression Matrix (#965)
- Support for hierarchical directory structures for experiments. (#1126)
- Pipeline Components: LearnTransform and ApplyTransform (#1083)
- Restructure the relationship between PipelineComponent and AlgorithmBase (#1095)


## [0.0.34] - 2019-03-21
- Adding ability to pass aligned group to Imagestack.from_path_or_url (#1069)
- Add Decoded Spot Table (#1087)
- Enable appending to existing napari viewer in display() (#1093)
- Change tile shape to a dict by default (#1072)
- Add ElementWiseMult Filter Pipeline Component (#983)
- Add linear unmixing pipeline component (#1056)
- Spiritual Bag of Images Refactor: Part 1 (#986)
- Add to provenance log   (#968)

## [0.0.33] - 2019.02.14
- Last release without a changelog!

[0.2.2]: https://github.com/spacetx/starfish/releases/tag/0.2.2
[0.2.1]: https://github.com/spacetx/starfish/releases/tag/0.2.1
[0.2.0]: https://github.com/spacetx/starfish/releases/tag/0.2.0
[0.1.10]: https://github.com/spacetx/starfish/releases/tag/0.1.10
[0.1.9]: https://github.com/spacetx/starfish/releases/tag/0.1.9
[0.1.8]: https://github.com/spacetx/starfish/releases/tag/0.1.8
[0.1.7]: https://github.com/spacetx/starfish/releases/tag/0.1.7
[0.1.6]: https://github.com/spacetx/starfish/releases/tag/0.1.6
[0.1.5]: https://github.com/spacetx/starfish/releases/tag/0.1.5
[0.1.4]: https://github.com/spacetx/starfish/releases/tag/0.1.4
[0.1.3]: https://github.com/spacetx/starfish/releases/tag/0.1.3
[0.1.2]: https://github.com/spacetx/starfish/releases/tag/0.1.2
[0.1.1]: https://github.com/spacetx/starfish/releases/tag/0.1.1
[0.1.0]: https://github.com/spacetx/starfish/releases/tag/0.1.0
[0.0.36]: https://github.com/spacetx/starfish/releases/tag/0.0.36
[0.0.35]: https://github.com/spacetx/starfish/releases/tag/0.0.35
[0.0.34]: https://github.com/spacetx/starfish/releases/tag/0.0.34
[0.0.33]: https://github.com/spacetx/starfish/releases/tag/0.0.33
---
title: 'starfish: scalable pipelines for image-based transcriptomics'
tags:
  - Python
  - xarray
  - skimage
  - microscopy
  - imaging
  - biology
  - single cell biology
  - spatial transcriptomics
  - MERFISH
  - In Situ Sequencing
  - osmFISH
  - smFISH
  - BaristaSeq
  - dartfish
  - starmap
  - seqFISH
authors:
  - name: Shannon Axelrod
    affiliation: 1
  - name: Matthew Cai
    orcid: 0000-0003-4998-6328
    affiliation: 1
  - name: Ambrose J. Carr
    orcid: 0000-0002-8457-2836
    affiliation: 1
  - name: Jeremy Freeman
    orcid: 0000-0001-7077-7972
    affiliation: 1
  - name: Deep Ganguli
    affiliation: 1
  - name: Justin T. Kiggins
    orcid: 0000-0002-4638-7015
    affiliation: 1
  - name: Brian Long
    orcid: 0000-0002-7793-5969
    affiliation: 2
  - name: Tony Tung
    affiliation: 1
  - name: Kevin A. Yamauchi
    orcid: 0000-0002-7818-1388
    affiliation: 3
affiliations:
 - name: Chan Zuckerberg Initiative
   index: 1
 - name: Allen Institute for Brain Science
   index: 2
 - name: Chan Zuckerberg Biohub
   index: 3
date: 24 April 2020
bibliography: paper.bib
---

# Summary

The exploding field of single cell transcriptomics has begun to enable deep analysis of gene expression and cell types, but spatial context is lost in the preparation of tissue for these assays.
Recent developments in biochemistry, microfluidics, and microscopy have come together to bring about an “alphabet soup” of technologies that enable sampling gene expression in situ, with varying levels of spatial resolution, sensitivity, and genetic depth.
These technologies promise to permit biologists to ask new questions about the spatial relationships between cell type and interactions between gene expression and cell morphology.
However, these assays generate very large microscopy datasets which are challenging to process using general microscopy analysis tools.
Furthermore, many of these assays require specialized analysis to decode gene expression from multiplexed experimental designs.

# Statement of Need

*starfish* is a Python library for processing images generated by microscopy-based spatial transcriptomics assays.
It lets biologists build scalable pipelines that localize and quantify RNA transcripts in image data generated by any hybridization- or sequencing-based *in situ* transcriptomics method, from classic RNA single-molecule FISH to combinatorial barcoded assays.
Image processing of an experiment is divided into fields of view (FOV) that correspond to the data produced by a microscope at a single location on a microscope slide.
*starfish* lets users register and pre-process images in each FOV, localize spots representing tagged RNA molecules in 3D, decode the identity of those molecules according to the experimental design, segment cells, assign the spots to cells, then aggregate spots into a cell x gene expression matrix.
This spatially-annotated gene expression matrix can then be analyzed and visualized in downstream tools for single-cell biology, such as Seurat [@Seurat], Bioconductor [@Bioconductor], Scanpy [@ScanPy], and cellxgene [@cellxgene].

To enable large scale processing of these data, *starfish* leverages a 5-dimensional imaging data model (x, y, z, round, channel) backed by the cloud-friendly `spacetx-format` file format, and [slicedimage](https://github.com/spacetx/slicedimage/), an interface for lazy, distributed loading of `spacetx-format` datasets.  
Furthermore, *starfish* implements comprehensive logging of all data processing steps for provenance tracking, reproducibility, and transparency.
*starfish* is built on top of popular Python tools like xarray [@xarray] and scikit-image [@scikit-image].

There are a number of other tools which support localization and quantification of spots in fluorescent microscopy images, including ImageJ and CellProfiler, however these tools do not support multiplexed decoding of gene targets necessary for many assays. Other tools which are designed more specifically to handle the kinds of assays that starfish supports include [dotdotdot](https://github.com/LieberInstitute/dotdotdot) [@dotdotdot], a MATLAB toolbox designed for RNAscope assays;  [pysmFISH](https://github.com/linnarsson-lab/pysmFISH/), a Python package designed for smFISH assays; and [SMART-Q](https://github.com/shenlab-ucsf/SMART-Q), a fork from an earlier development release of starfish adding support for immunostaining and other features [@SMART-Q].

starfish requires a working knowledge of Python and fluorescent image analysis for a user to create an analysis pipeline.
To help new users get started and support the broader single cell biology community in learning how to work with these data, *starfish* maintains example datasets and reference implementations ported from published assays, including
MERFISH [@MERFISH],
In Situ Sequencing [@ISS],
osmFISH [@osmFISH],
BaristaSeq [@BaristaSeq],
smFISH [@smFISH],
DARTFISH [@DARTFISH],
STARmap [@starmap],
and seqFISH [@seqFISH].
To take advantage of starfish's support for large scale processing, users must have familiarity with cluster or cloud computing.

*starfish* was developed alongside the [SpaceTx project](https://spacetx-starfish.readthedocs.io/en/stable/about/index.html), a CZI-funded effort to compare spatial transcriptomics methods in the context of determining cell types in the brain [@starfish-sfn].
*starfish* is currently in use by multiple research groups, including the [Allen Institute for Brain Science](https://alleninstitute.org), the [Chan Zuckerberg Biohub](https://www.czbiohub.org/), and the [Zhang Lab at UC San Diego](http://jinzhanglab.ucsd.edu/).
These groups support multiple large-scale projects profiling *in situ* gene expression, including the SpaceTx consortium, the [Human Cell Atlas](https://www.humancellatlas.org/), the [BRAIN Initiative Cell Census Network](https://biccn.org/), and the [HuBMAP Consortium](https://hubmapconsortium.org/).


# Acknowledgements

We would like to acknowledge the following contributions, which have been invaluable to the development of starfish.

For direct contributions of code and documentation, we thank
Olga Botvinnik,
Gökçen Eraslan,
Kira Evans,
Marcus Kinsella,
Nicholas Mei,
Josh Moore,
Nicholas Sofroniew,
and Ola Tarkowska.

We would like to thank the members of the SpaceTx consortium for working closely with us to understand their image processing steps, port pipelines into starfish, and contribute example datasets.
We especially would like to thank
Jeff Moffitt and Xiaowei Zhuang for their help with the MERFISH pipeline and example data,
Marco Mignardi and Mats Nilsson for their help with the In Situ Sequencing pipeline and example data,
Simone Codeluppi and Sten Linnarsson for their help with the osmFISH pipeline and example data,
Xioayin Chen and Anthony Zador for their help with the BaristaSeq pipeline and example data,
Richard Que and Kun Zhang for their help with the DARTFISH pipeline and example data,
Xiao Wang, William Allen, and Karl Deisseroth for their help with the STARmap pipeline and example data,
Dan Goodwin and Ed Boyden for their help with the ExFISH pipeline and example data,
Nico Pierson, Sheel Shah, and Long Cai for their help with the seqFISH pipeline and example data,
and Denis Shapiro for their help with the Imaging Mass Cytometry example data.
Brian Long would like to thank the Allen Institute founder, Paul G. Allen, for his vision, encouragement and support.
Finally, we would like to thank the Science and Technology teams and the Single Cell Program at the Chan Zuckerberg Initiative for their support of this work.

# References
---
name: "\U0001F4D8 User Story"
about: This is a template for a user story
title: ''
labels: feature
assignees: ''

---

**As a** _______________________,
**I want** _______________________ **so that** _______________________.

### Acceptance Criteria

1. [If I do A.]
1. [B should happen.]

[
Also, here are a few points that need to be addressed:

1. Constraint 1;
1. Constraint 2;
1. Constraint 3.
]


### Resources:

* Sketches: [Here goes a URL to sketches]


### Notes

[Some complementary notes if necessary:]

* > Here goes a quote from an email
* Here goes whatever useful information can exist…
---
name: "\U0001F41B Bug report"
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

#### Description
<!-- Example: Joblib Error thrown when calling fit on LatentDirichletAllocation with evaluate_every > 0-->

#### Steps/Code to Reproduce
<!--
If the code is too long, feel free to put it in a public gist and link
it in the issue: https://gist.github.com
-->

#### Expected Results
<!-- Example: No error is thrown. Please paste or describe the expected results.-->

#### Actual Results
<!-- Please paste or specifically describe the actual output or traceback. -->

#### Versions
<!--
import platform; print(platform.platform())
import sys; print("Python", sys.version)
import numpy; print("NumPy", numpy.__version__)
import scipy; print("SciPy", scipy.__version__)
import skimage; print("scikit-image", skimage.__version__)
import pandas; print("pandas", pandas.__version__)
import sklearn; print("sklearn", sklearn.__version__)
import xarray; print("xarray", xarray.__version__)
import sympy; print("sympy", sympy.__version__)
import starfish; print("starfish", starfish.__version__)
-->


<!-- Thanks for contributing! -->
---
name: "\U00002728 Feature request"
about: Suggest an idea for this project
title: ''
labels: feature
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
# SpaceTx Image Format Specification

## Introduction

This document describes the SpaceTx file format specification for image-based biological assays.
The SpaceTx format is designed to support the construction of a data tensor, which generalizes
the approaches taken by sequential single-molecule FISH and assays that identify targets
by building codes over multiple imaging rounds. Each of theses assays produce images that can
form a data tensor. The data tensor contains a series of (x, y) planar image planes that
represent specific z-planes, imaging channels (c), and imaging rounds (r). Together these form a
5-Dimensional tensor (r, c, z, y, x) that serves as a general representation of an image-based 
transcriptomics or proteomics assay, and is the substrate of the starfish package.

This format should be self-describing and it should specify both how a set of 2D images form a 
field of view and how multiple fields of view interact to form a larger experiment. The SpaceTx
format accomplishes this by combining these images, stored in 2-Dimensional image formats, 
with a series of JSON files that describe how to organize each image file into the 5-Dimensional 
imaging tensor. Combined with imaging metadata and a pipeline recipe, both of which are defined
elsewhere, these files enable a pipeline to generate the desired outputs of a spatial assay: a gene 
expression matrix augmented with spatial locations of transcripts and cells.

## Format Specification

Here, we tabulate the minimum set of required json files that can describe the
imaging data of a spaceTx experiment with brief descriptions of their purpose:

| Type            | Description                                                                                              |
|:----------------|:---------------------------------------------------------------------------------------------------------|
| Experiment      | links the data manifests and codebook together                                                           |
| Manifest        | file locations of each field of view                                                                     |
| Field of View   | describes how individual 2-d image planes form an image tensor                                           |
| Codebook        | maps patterns of intensity in the channels and rounds of a field of view to target molecules             |

Each of these input types and their file formats are described in detail in the following sections.

## Experiment

The data manifest is a JSON file that ties together all information about the data images, auxiliary images (like nuclear stains), and the codebook needed to decode the experiment.
It is the file read by starfish to load data into the analysis environment.

Example:
```json
{
  "version": "0.0.0",
  "images": {
    "primary": "primary_images.json",
    "nuclei": "nuclei.json"
  },
  "codebook": "codebook.json",
  "extras": {
    "is_space_tx_cool": true
  }
}
```

## Manifest

Both the `primary_images.json` and `nuclei.json` files referenced by the above `experiment.json` may contain links to Field of View Manifests (for simple experiments with only one field of view, these fields may also directly reference a field of view).
The Manifest is a simple association of a field of view name with the json file that defines the field of view.
In this example, we demonstrate a primary images manifest with three fields of view.
Such an experiment would likely also have a nuclei manifest, which would _also_ contain three fields of view.

```json
{
  "version": "0.0.0",
  "contents": {
    "fov_000": "primary-images-fov_000.json",
    "fov_001": "primary-images-fov_001.json",
    "fov_00N": "primary-images-fov_002.json"
  },
  "extras": null
}
```

## Field of View

The field of view is the most complex file in the spaceTx format, and must be created for each data tensor and auxiliary image tensor in an experiment.
It provides two key types of information: information about the field of view, and information about each tile contained in it.

The field_of_view.json file specifies the shape of the image tensor, including the size of the (X, Y) image in pixels, and the number of z-planes, imaging channels, and imaging rounds in the experiment.
Thus, an image tensor has shape (r, c, z, y, x), though y and x are limited to at most 3000 pixels.
For experiments that do not leverage all of these concepts, the values can simply be set to one, and that dimension of the tensor will be ignored.
For example, smFISH experiments that do not leverage multiple imaging rounds have shape (1, c, z, y, x).

For each individual tile, the Field of View specifies the portion of the tensor the tile corresponds to by providing the indices of the tile in (r, c, z), the location of the tile, and the sha256 hash of the file data, to guard against corruption.

Finally, each tile also specifies the coordinates of the image in physical space, relative to some experiment-wide reference point specified in micrometers.

The below example describes a 4-channel, 3-round barcoded experiment that samples a tissue section
using 1 discrete z-plane. For conciseness, the tile data is truncated, and shows only the
information for two tiles, while in practice there would be 4 * 3 * 1 tiles.

```json
{
    "default_tile_format": "TIFF",
    "dimensions": [
        "z",
        "xc",
        "x",
        "yc",
        "y",
        "zc",
        "c",
        "r"
    ],
    "extras": {},
    "shape": {
        "c": 4,
        "r": 3,
        "z": 1
    },
    "tiles": [
        {
            "coordinates": {
                "xc": [
                    0.0,
                    0.112
                ],
                "yc": [
                    0.0,
                    0.0539
                ],
                "zc": [
                    0.0,
                    0.0001
                ]
            },
            "file": "primary-fov_000-c0-r0-z0.tiff",
            "indices": {
                "c": 0,
                "r": 0,
                "z": 0
            },
            "sha256": "bbd9098fa11918ba4e09672789000fa94c0fec4128b071a6b5dfb3b2f4d04df8",
            "tile_format": "TIFF",
            "tile_shape": {
                "x": 1120,
                "y": 539
            }
        },
        {
            "coordinates": {
                "xc": [
                    0.0,
                    0.112
                ],
                "yc": [
                    0.0,
                    0.0539
                ],
                "zc": [
                    0.0,
                    0.0001
                ]
            },
            "file": "primary-fov_000-c1-r0-z0.tiff",
            "indices": {
                "c": 1,
                "r": 0,
                "z": 0
            },
            "sha256": "4692207a483e6c482db37239af9ba4867d05caf80795aabc5b66a943cb0d60df",
            "tile_format": "TIFF",
            "tile_shape": {
                "x": 1120,
                "y": 539
            }
        }
    ],
    "version": "0.1.0"
}
```

.. _sptx_codebook_format:

## Codebook

The final part of the spaceTx specification, the codebook describes how intensities detected in the image tensor correspond to the targets of the assay.
The codebook is an array, where each object in the array lists a codeword and the target it corresponds to.
Each codeword is made up of one or more json objects, each of which describe the expected intensity value for tiles of specific (channel, round) combinations.

For smFISH experiments where each channel corresponds to a different target and there is only one imaging round, the codebook is very simple:

```json
{
  "version": "0.0.0",
  "mappings": [
    {
      "codeword": [
        {"c": 0, "r": 0, "v": 1}
      ],
      "target": "SCUBE2"
    },
    {
      "codeword": [
        {"c": 1, "r": 0, "v": 1}
      ],
      "target": "BRCA"
    },
    {
      "codeword": [
        {"c": 2, "r": 0, "v": 1}
      ],
      "target": "ACTB"
    }
  ]
}
```
In this example, channels 0, 1, and 2 correspond to `SCUBE2`, `BRCA`, and `ACTB`, respectively.
In contrast, a barcoded experiment may have a more complex codebook:

```json
{
  "version": "0.0.0",
  "mappings": [
    {
      "codeword": [
        {"r": 0, "c": 0, "v": 1},
        {"r": 0, "c": 1, "v": 1}
      ],
      "target": "SCUBE2"
    },
    {
      "codeword": [
        {"r": 0, "c": 0, "v": 1},
        {"r": 1, "c": 1, "v": 1}
      ],
      "target": "BRCA"
    },
    {
      "codeword": [
        {"r": 0, "c": 1, "v": 1},
        {"r": 1, "c": 0, "v": 1}
      ],
      "target": "ACTB"
    }
  ]
}
```

The above example describes the coding scheme of an experiment with 2 rounds and 2 channels, where each code expects exactly two images out of four to produce signal for a given target.
In the above example, a spot in the image tensor would decode to `SCUBE2` if the spot was detected in (round=0, channel=0) and (round=0, channel=1).

Note that the codebook only states opinions about non-zero expected fluorescence values, and is not designed to distinguish between a literal `0` and a missing value, which would imply "any fluorescence level". For example, the following two codewords are treated identically and cannot be distinguished at the level of the codebook.

```json
[
  {"r": 0, "c": 1, "v": 1},
  {"r": 1, "c": 0, "v": 0}
]

[
  {"r": 0, "c": 1, "v": 1},
]
```

Libraries leveraging SpaceTx Format can implement logic to support distinctions between these types by implementing a flag that allows the user to specify whether missing or zero `(round, channel)` combinations should be instantiated as `0` or `null`. The latter is needed to support smFISH experiments for which each `(round, channel)` pair codes for a separate gene, and for which fluorescence intensities in other channels are _irrelevant_ and should not be considered by a decoder function.

The flexibility of this codebook format to describe any coding scheme also increases the
potential that errors in making the codebook go undetected. The `starfish validate codebook` command
should be used to check for egregious mistakes, but the only way to ensure the codebook is
correct is to have a good understanding of the format and inspect the json file.## JSON Schemas

Each of the json files that comprise a SpaceTx image fileset can be validated
against one of several jsonschemas ([Draft 4](http://json-schema.org/specification-links.html#draft-4)).

| Schema location                        | Description                                                                            |
|:---------------------------------------|:---------------------------------------------------------------------------------------|
| codebook/codebook.json                 | maps patterns of intensity in the channels and rounds of a field of view to target molecules |
| codebook/codeword.json                 | describes the individiual codes contained in a codebook                                | 
| experiment.json                        | top-level object mapping manifests and codebook together                               |
| extras.json                            | extension point used in multiple schemas for storing key/value pairs                   |
| field_of_view/tiles/indices.json       | describes the categorical indices (channel, round, and z-section) of a tile            |
| field_of_view/tiles/coordinates.json   | physical coordinates of a tile                                                         |
| field_of_view/tiles/tiles.json         | specification of a 2-D image tile                                                      |
| field_of_view/field_of_view.json       | 5-D image consisting of multiple 2-D image tiles                                       |
| fov_manifest.json                      | manifest listing one or more json files that each describe a field of view             |
| version.json                           | general purpose version specification used in multiple schemas                         |

See the API and Usage documentation for how to validate your documents against the jsonschema from your
code or the command-line, respectively.
