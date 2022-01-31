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
Contributing to Starfish
========================

If you're reading this section, you're probably interested in contributing to Starfish.  Welcome and thanks for your interest in contributing!

How can I contribute?
---------------------

Starfish is designed to specify pipeline recipes for image processing. To support this use, the library is composed as a series of `pipeline_component` modules.
The objects that sub-class `pipeline_component` spawn a command-line interface that should identically track the API of the python library.

A typical starfish run consists of running one or more image processing filter stages, and identifying features through either a spot- or pixel-based approach.
The identified features are then decoded into the genes that they correspond to by mapping the fluorescence channel (and optionally hybridization round) using a codebook.
Finally, the filtered data are segmented, identifying which cell each feature belongs to.

Creating a new algorithm for an existing `pipeline_component`
-------------------------------------------------------------

For example, to add a new image filter, one would:

1. Create a new python file `new_filter.py` in the `starfish/core/image/Filter/` directory.
2. Find the corresponding `AlgorithmBase` for your component.
   For filters, this is `FilterAlgorithm`, which is found in `starfish/core/image/Filter/_base.py`.
   Import that base into `new_filter.py`, and have your new algorithm subclass it,
   e.g. create `class NewFilter(FilterAlgorithm)`
3. Implement all required methods from the base class.

That's it! If at any point something gets confusing, it should be possible to look at existing pipeline components of
the same category for guidance on implementation.

Reporting bugs
--------------

- Bugs can be contributed as issues in the starfish repository.
  Please check to make sure there is no existing issue that describes the problem you
  have identified before adding your bug.
- When reporting issues please include as much detail as possible about your operating system,
  starfish version, slicedimage version, and python version. Whenever possible, please also include a brief,
  self-contained code example that demonstrates the problem, including a full traceback.

Code contributions
------------------

- Don't break the build: pull requests are expected to pass all automated CI checks.
  You can run those checks locally by running `make all` in starfish repository root.
- All Pull Request comments must be addressed, even after merge.
- All code must be reviewed by at least 1 other team member.
- All code must have typed parameters, and when possible, typed return calls (see
  `PEP484 <https://www.python.org/dev/peps/pep-0484>`_).
  We also encourage rational use of in-line variable annotation when the type of a newly defined object is not clear
  (see `PEP526 <https://www.python.org/dev/peps/pep-0526/>`_.
- All code must be documented according to `numpydoc <https://numpydoc.readthedocs.io/en/latest/>`_ style guidelines.
- Numpy provides an excellent `development workflow <https://docs.scipy.org/doc/numpy/dev/gitwash/development_workflow.html>`_
  that we encourage you to follow when developing features for starfish!
- All commits must have informative summaries; try to think about what will still make sense looking back on them next year.
- All pull requests should describe your testing plan.
- When merging a pull request, squash commits down to the smallest logical number of commits. In cases where a single commit
  suffices, use the "Squash and Merge" strategy, since it adds the PR number to the commit name. If multiple commits remain,
  use "Rebase and Merge".

Notebook contributions
----------------------

- All `.ipynb` files should have a corresponding `.py` file.
  Use `nbencdec <https://github.com/ttung/nbencdec>`_ to generate the corresponding `.py` file.
- The `.py` files allow refactor commands in the codebase to find code in the `.py` files,
  which is an important to keep the notebooks working as starfish evolves.


Project Tracking
-----------------
- The starfish team uses `zenhub <https://app.zenhub.com/workspaces/starfish-dev-5b4e05b4c93e4717b2160fdb/board>`_ to track the team's progress
and upcoming engineering raodmap.
- If you would like to see the team's board in github and the correct rendering of epic issues it's recommended that you download the chrome extension
`here <https://www.zenhub.com/extension>`_
starfish: scalable pipelines for image-based transcriptomics
========


.. image:: docs/source/_static/design/logo.png
    :scale: 50 %

-----------------------------

.. image:: https://img.shields.io/badge/dynamic/json.svg?label=forum&url=https%3A%2F%2Fforum.image.sc%2Ftags%2Fstarfish.json&query=%24.topic_list.tags.0.topic_count&colorB=brightgreen&suffix=%20topics&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAA4AAAAOCAYAAAAfSC3RAAABPklEQVR42m3SyyqFURTA8Y2BER0TDyExZ+aSPIKUlPIITFzKeQWXwhBlQrmFgUzMMFLKZeguBu5y+//17dP3nc5vuPdee6299gohUYYaDGOyyACq4JmQVoFujOMR77hNfOAGM+hBOQqB9TjHD36xhAa04RCuuXeKOvwHVWIKL9jCK2bRiV284QgL8MwEjAneeo9VNOEaBhzALGtoRy02cIcWhE34jj5YxgW+E5Z4iTPkMYpPLCNY3hdOYEfNbKYdmNngZ1jyEzw7h7AIb3fRTQ95OAZ6yQpGYHMMtOTgouktYwxuXsHgWLLl+4x++Kx1FJrjLTagA77bTPvYgw1rRqY56e+w7GNYsqX6JfPwi7aR+Y5SA+BXtKIRfkfJAYgj14tpOF6+I46c4/cAM3UhM3JxyKsxiOIhH0IO6SH/A1Kb1WBeUjbkAAAAAElFTkSuQmCC
    :target: https://forum.image.sc/tag/starfish
    :alt: Image.sc forum

.. image:: https://img.shields.io/pypi/v/starfish
    :target: https://pypi.org/project/starfish/
    :alt: PyPI

.. image:: https://img.shields.io/pypi/dm/starfish
   :target: https://pypistats.org/packages/starfish
   :alt: PyPI - Downloads

.. image:: https://readthedocs.org/projects/spacetx-starfish/badge/?version=latest
    :target: https://spacetx-starfish.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://travis-ci.com/spacetx/starfish.svg?branch=master
    :target: https://travis-ci.com/spacetx/starfish

.. image:: https://codecov.io/gh/spacetx/starfish/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/spacetx/starfish

-------------------------------

*starfish* is a Python library for processing images of image-based spatial transcriptomics.
It lets you build scalable pipelines that localize and quantify RNA transcripts in image data
generated by any FISH method, from simple RNA single-molecule FISH to combinatorial barcoded
assays.

- `Starfish enterprise: finding RNA patterns in single cells <https://www.nature.com/articles/d41586-019-02477-9>`_, Nature Technology Feature, Aug 2019
- `Developing a Computational Pipeline and Benchmark Datasets for Image-Based Transcriptomics <https://www.ascb.org/science-news/developing-a-computational-pipeline-and-benchmark-datasets-for-image-based-transcriptomics/>`_, ASCB Science News, Dec 2018

Documentation
-------------

See `spacetx-starfish.readthedocs.io <https://spacetx-starfish.readthedocs.io/en/latest/>`_ for the `quickstart <https://spacetx-starfish.readthedocs.io/en/latest/gallery/quick_start/plot_quick_start.html>`_, `user guide <https://spacetx-starfish.readthedocs.io/en/latest/user_guide/>`_, `examples <https://spacetx-starfish.readthedocs.io/en/latest/gallery/>`_, and `API <https://spacetx-starfish.readthedocs.io/en/latest/api/>`_.

Installation
------------

starfish supports python 3.7 and above and can easily be installed from PyPI:

.. code-block:: bash

    $ pip install starfish[napari]

`For more detailed installation instructions, see here <https://spacetx-starfish.readthedocs.io/en/latest/installation/>`_.

Python Version Notice
---------------------

starfish will be dropping support for python 3.6 in the next release due to
minimum python=3.7 version requirements in upstream dependencies.

Contributing
------------

We welcome contributions from our users! See our contributing.rst_ and `developer guide`_ for more information.

.. _`developer guide`: https://spacetx-starfish.readthedocs.io/en/latest/developer_guide/
.. _contributing.rst: https://github.com/spacetx/starfish/blob/master/CONTRIBUTING.rst

Help, support, and questions
----------------

- Forum: ask on the `Image.sc forum <https://forum.image.sc/tag/starfish>`_
- Email: `starfish@chanzuckerberg.com <mailto:starfish@chanzuckerberg.com>`_
.. raw:: html

    <style>
    #banner_img {
        background: url("_static/design/ivana-cajina-asuyh-_ZX54-unsplash.png") no-repeat center;
        /*Photo by Ivana Cajina on Unsplash*/
        background-size: cover;
    }
    </style>

    <div class="jumbotron jumbotron-fluid" id="banner_img">
    <div class="container-fluid">

    <h1>starfish: scalable pipelines for image-based transcriptomics</h1>

    <div style="clear: both"></div>
    <div class="container-fluid">
      <div class="row">
        <div class="col-lg-7">
          <div class="bs-component">
              <a type="button" class="btn btn-primary" href="gallery/quick_start/plot_quick_start.html">quickstart</button></a>
              <a type="button" class="btn btn-primary" href="https://github.com/spacetx/starfish">github</button></a>
              <a type="button" class="btn btn-primary" href="https://forum.image.sc/tag/starfish">user support</button></a>
          </div>
        </div>
      </div>
    </div>
    </div>
    </div>
    <br>

    <div class="container-fluid">
    <div class="row">
    <div class="col-md-7">
    <br>

*starfish* is a Python library for processing images of image-based spatial transcriptomics.
It lets you build scalable pipelines that localize and quantify RNA transcripts in image data
generated by any FISH method, from simple RNA single-molecule FISH to combinatorial barcoded
assays. Image processing of an experiment is divided into fields of view (FOV) that correspond to
the data produced by a microscope at a single location on a microscope slide. Each FOV is
processed into a table of spots localized in 3D and then spots are aggregated into a cell x gene
expression matrix by comparing physical positions of spots and cells.

If you want to give it a try, the :ref:`quick start tutorial <quick start>` will guide you from
installation to running a pipeline in under 10 minutes. For more comprehensive instructions on
how to use starfish, see the user guide on :ref:`creating image-based transcriptomics processing
pipelines <user_guide>`, which organizes and contextualizes the
tutorials on running starfish using the API. Finally, advanced users can examine the
:ref:`Data Structures <data structures>` and :ref:`Help & Reference <help and reference>`
sections to learn more details about starfish and its object models.

In addition to the library of image processing functions, starfish introduces a standardized data
format for image-based spatial transcriptomic assays. Examples and tutorials for formatting
your data into the SpaceTx format can be found in :ref:`section_formatting_data`.

.. raw:: html

  </div>
  <div class="col-md-5">
  <h2>Features</h2>

* Formatting Data: :ref:`API<experiment_builder>` | :ref:`Tutorial<section_formatting_data>`
* Registering Images: :ref:`API<learn_transform>` | :ref:`Tutorial<tutorial_image_registration>`
* Filtering Images: :ref:`API<filtering>` | :ref:`Tutorial<section_correcting_images>`
* Finding Spots: :ref:`API<spot_finding>` | :ref:`Tutorial<section_finding_and_decoding>`
* Decoding Spots: :ref:`API<decode_spots>` | :ref:`Tutorial<section_finding_and_decoding>`
* Segmenting Cells: :ref:`API<segmentation>` | :ref:`Tutorial<section_segmenting_cells>`
* Assigning Spots: :ref:`API<target_assignment>` | :ref:`Tutorial<tutorial_assigning_spots>`
* Example Pipelines: :ref:`Gallery<pipeline_examples>`

.. raw:: html

  </div>
  <div class="col-md-5">
  <h2>Contents</h2>

.. toctree::
   :maxdepth: 1

   installation/index
   user_guide/index
   gallery/index
   help_and_reference/index
   api/index
   developer_guide/index
   about/index

.. raw:: html

   </div>
   </div>
   </div>
.. _contributing:

Developer Guide
===============

Welcome and thanks for your interest in contributing!

How can I contribute?
---------------------

Starfish is designed to specify pipeline recipes for image processing. To support this use, the library is composed as a series of `pipeline_component` modules.
The objects that sub-class `pipeline_component` spawn a command-line interface that should identically track the API of the python library.

A typical starfish run consists of running one or more image processing filter stages, and identifying features through either a spot- or pixel-based approach.
The identified features are then decoded into the genes that they correspond to by mapping the fluorescence channel (and optionally hybridization round) using a codebook.
Finally, the filtered data are segmented, identifying which cell each feature belongs to.

Installing *starfish* for developers
------------------------------------

If you are on a mac, make sure you have the `XCode CommandLine Tools`_
installed.  Check out the code for starfish and set up a `virtual environment`_.

.. _`XCode CommandLine Tools`: https://developer.apple.com/library/archive/technotes/tn2339/_index.html
.. _`virtual environment`: #using-virtual-environments

.. code-block:: bash

    $ git clone git://github.com/spacetx/starfish.git
    $ cd starfish

Then install starfish:

.. code-block:: bash

    $ make install-dev

Creating a new algorithm for an existing `pipeline_component`
-------------------------------------------------------------

For example, to add a new image filter, one would:

1. Create a new python file `new_filter.py` in the `starfish/core/image/Filter/` directory.
2. Find the corresponding `AlgorithmBase` for your component.
   For filters, this is `FilterAlgorithm`, which is found in `starfish/core/image/Filter/_base.py`.
   Import that base into `new_filter.py`, and have your new algorithm subclass it,
   e.g. create `class NewFilter(FilterAlgorithm)`
3. Implement all required methods from the base class.

That's it! If at any point something gets confusing, it should be possible to look at existing pipeline components of
the same category for guidance on implementation.

Reporting bugs
--------------

- Bugs can be contributed as issues in the starfish repository.
  Please check to make sure there is no existing issue that describes the problem you
  have identified before adding your bug.
- When reporting issues please include as much detail as possible about your operating system,
  starfish version, slicedimage version, and python version. Whenever possible, please also include a brief,
  self-contained code example that demonstrates the problem, including a full traceback.

Code contributions
------------------

- Don't break the build: pull requests are expected to pass all automated CI checks.
  You can run those checks locally by running `make all` in starfish repository root.
- All Pull Request comments must be addressed, even after merge.
- All code must be reviewed by at least 1 other team member.
- All code must have typed parameters, and when possible, typed return calls (see
  `PEP484 <https://www.python.org/dev/peps/pep-0484>`_).
  We also encourage rational use of in-line variable annotation when the type of a newly defined object is not clear
  (see `PEP526 <https://www.python.org/dev/peps/pep-0526/>`_.
- All code must be documented according to `numpydoc <https://numpydoc.readthedocs.io/en/latest/>`_ style guidelines.
- Numpy provides an excellent `development workflow <https://docs.scipy.org/doc/numpy/dev/gitwash/development_workflow.html>`_
  that we encourage you to follow when developing features for starfish!
- All commits must have informative summaries; try to think about what will still make sense looking back on them next year.
- All pull requests should describe your testing plan.
- When merging a pull request, squash commits down to the smallest logical number of commits. In cases where a single commit
  suffices, use the "Squash and Merge" strategy, since it adds the PR number to the commit name. If multiple commits remain,
  use "Rebase and Merge".

Notebook contributions
----------------------

- All `.ipynb` files should have a corresponding `.py` file.
  Use `nbencdec <https://github.com/ttung/nbencdec>`_ to generate the corresponding `.py` file.
- The `.py` files allow refactor commands in the codebase to find code in the `.py` files,
  which is an important to keep the notebooks working as starfish evolves.

Debugging Errors
----------------

First, thank you for using Starfish and SpaceTx-Format! Feedback you provide on features and the
user experience is critical to making Starfish a successful tool. Because we iterate quickly on this
feedback to add new features, things change often, which can result in your code getting out of sync
with your data. When that happens, you may observe errors.

Most of the time, you can fix this problem by pulling the most recent version of the code,
reinstalling starfish, and restarting your environment. If you're using starfish with datasets from
spaceTx located on our cloudfront distribution, we're committed to keeping that data up to date.
Updated versions of the notebook will reference the correct data version, and copying over the
new link should fix any issues.

For example, if a notebook references in-situ sequencing data from August 23rd, and a breaking
change occurs on September 26th, it would be necessary to replace the experiment link to point at
data that was updated to work post-update:

.. code-block:: diff

    - http://spacetx.starfish.data.public.s3.amazonaws.com/browse/formatted/20180823/iss_breast/experiment.json
    + http://spacetx.starfish.data.public.s3.amazonaws.com/browse/formatted/20180926/iss_breast/experiment.json

If you're using your own data with starfish, you may need to re-run your data ingestion workflow
based on :py:class:`starfish.experiment.builder.providers.TileFetcher` and
:py:class:`starfish.experiment.builder.providers.FetchedTile` to generate up-to-date versions of spaceTx-format.

Upgrading to a new version
--------------------------

If you've installed from pypi, upgrading is as simple as reinstalling starfish.

.. code-block:: bash

    pip install --upgrade starfish

If you've installed our development version to take advantage of new features in real time, you'll
need to fetch changes and reinstall. Assuming you've cloned the respository into ``./starfish``,
you can install the newest version as follows:

.. code-block:: bash

    cd ./starfish
    git checkout master
    git pull
    pip3 install .

Reporting bugs
--------------

Bugs can be contributed as issues in the starfish repository. Please check to make sure there
is no existing issue that describes the problem you have identified before adding your bug.

When reporting issues please include as much detail as possible about your operating system,
starfish version, slicedimage version, and python version. Much of this can be accomplished by
sending us the output of ``pip freeze``:

.. code-block:: bash

    pip freeze > environment.txt

Whenever possible, please also include a brief, self-contained code example that demonstrates the
problem, including a full traceback.
.. _installation:

Installation
============

Starfish supports python 3.6 and above (python 3.7 recommended). To install the starfish package,
first verify that your python version is compatible. You can check this by running :code:`python
--version`.

The output should look similar to this:

.. code-block:: bash

   $ python --version
   Python 3.7.7

.. warning::
    While starfish itself has no known issues with python 3.8, scikit-image is not fully
    compatible with python 3.8. As such, installation of scikit-image, as part of starfish
    installation, may unexpectedly fail. The workaround is to install numpy first before
    installing starfish or scikit-image.


Using virtual environments
--------------------------

Starfish lists minimum versions for its dependencies for access to new features and algorithms.
These more up-to-date packages may create conflicts in your existing scripts or other packages,
so we recommend using a virtualenv_. You can create a work folder and set up the virtual
environment like:

.. _virtualenv: https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments

.. code-block:: bash

    $ mkdir starfish
    $ cd starfish
    $ python -m venv .venv
    $ source .venv/bin/activate

Conda_ users can set one up like so:

.. _Conda: https://www.anaconda.com/distribution/

.. code-block:: bash

    $ conda create -n starfish "python=3.7"
    $ conda activate starfish

Installing *starfish*
---------------------

Starfish can easily be installed using pip:

.. code-block:: bash

    $ pip install starfish

.. note::
    If using python 3.8, first install numpy using pip before installing starfish.

To use napari for interactive image visualization via :py:func:`.display` you must also
install napari:

.. code-block:: bash

    $ pip install starfish[napari]

Interactive visualization with napari also requires using Qt (e.g. by running the magic command
`%gui qt` in a jupyter notebook or ipython shell.)

Installing *starfish* on Windows
--------------------------------

Windows users can install starfish in the same way. Again, we recommend using a conda or virtual
environment with python 3.7. Here is how you would install starfish in a virtual environment
created with python's ``venv`` module:

.. code-block:: bat

    > mkdir starfish
    > cd starfish
    > python -m venv .venv
    > .venv\Scripts\activate.bat
    > pip install starfish
    > pip install starfish[napari]

.. note::
    Python 3.8 has trouble installing scikit-image v0.15.0 and the ``pip install numpy``
    workaround does not solve this issue on Windows.

Jupyter notebook
----------------

To run starfish in a jupyter notebook (recommended for creating an image processing pipeline) add
the virtualenv kernel to jupyter by activating your virtual environment and then:

.. code-block:: bash

    $ python -m ipykernel install --user --name=<venv_name>

Now you should be able to select ``venv_name`` as the kernel in a jupyter notebook to have access
to the starfish library.
.. _about:

About
=====

Starfish is developed and maintained by the Science team at the Chan Zuckerberg Initiative.
It began as a project to support the SpaceTx consortium.
SpaceTx is a consortium effort to benchmark image based transcriptomic methods by applying 10 different methods on a common tissue source, standardizing the raw data formats and using standardized analysis pipelines.

.. _citing:

Citing Starfish
---------------

to cite starfish, please use::

    Axelrod S, Carr AJ, Freeman J, Ganguli D, Long B, Tung T, and others.
    Starfish: Open Source Image Based Transcriptomics and Proteomics Tools, 2018-,
    http://github.com/spacetx/starfish [Online; accessed <date>].

Here’s an example of a BibTeX entry::

    @misc{,
    author = {
        Shannon Axelrod, Ambrose J Carr, Jeremy Freeman, Deep Ganguli, Brian Long, Tony Tung,
        and others
    },
    title = {{Starfish}: Open Source Image Based Transcriptomics and Proteomics Tools},
    year = {2018--},
    url = "http://github.com/spacetx/starfish",
    note = {[Online; accessed <date>]}
    }

.. _license:

License
-------

.. literalinclude:: ../../../LICENSE
.. _user_guide:

User Guide
==========

Welcome to the user guide for building an image processing pipeline using starfish! This tutorial
will cover all the steps necessary for going from raw images to a single cell gene expression
matrix. If you only have a few minutes to try out starfish, check out the
:ref:`Quick Start<quick start>` to see a demonstration of how starfish works. starfish
requires a working knowledge of Python and fluorescent image analysis to create an analysis pipeline.
If you are ready to learn how to build your own image processing pipeline using starfish then read on!

This part of the tutorial goes into more detail about why each of the stages in the example are
needed, and provides some alternative approaches that can be used to build similar pipelines.

The core functions of starfish pipelines are the detection (and :term:`decoding<Decoding>`)
of spots, and the segmentation of cells. Each of the other approaches are designed to address
various characteristics of the imaging system, or the optical characteristics of the tissue
sample being measured, which might bias the resulting spot calling, decoding, or cell
segmentation decisions. Not all parts of image processing are always needed; some are dependent
on the specific characteristics of the tissues. In addition, not all components are always found
in the same order. *Starfish* is flexible enough to omit some pipeline stages or disorder them,
but the typical order might match the following. The links show how and when to use each
component of *starfish*, and the final section demonstrates putting together a "pipeline recipe"
and running it on an experiment.

.. _section_formatting_data:

Formatting Data
---------------

In order to load the experiment into a starfish pipeline the data must be in
:ref:`SpaceTx Format<sptx_format>`, which is a standardized format that utilizes json
files to organize single-plane tiffs for image-based spatial transcriptomics data. If the data
you want to process isn't already in SpaceTx Format, there are a few methods to convert
your data.

.. note::

    When converting data to SpaceTx Format is too costly, images can be loaded directly without
    formatting by :ref:`tilefetcher_loader`. This is a workaround and only recommended if
    reading and writing all the images is infeasible. The experiment JSON files like the codebook
    will still need to be created.

The primary method is to use :py:func:`.format_structured_dataset`, a conversion tool, on
data that is structured as 2D image tiles with specific filenames and a CSV
file containing the physical coordinates of each image tile. This method requires minimal Python
knowledge. You can manually organize your images, but for large datasets you will want to use a
script (e.g. Matlab) to move and rename files into the structured data format. The structured
data must be 2D image tiles in TIFF, PNG, or NPY file format.

Users who are familiar with Python and starfish also have the option of using
:py:class:`.TileFetcher` and :py:class:`.FetchedTile`, a set of user-defined interfaces the
experiment builder uses for obtaining the data corresponding to each tile location. Any data
format that can be read as a numpy array can be used.

Lastly, there is a 3rd party `spacetx writer`_ which writes SpaceTx-Format experiments using the
`Bio-Formats`_ converter. Bio-Formats can read a variety of input formats, so might be a
relatively simple approach for users familiar with those tools.

.. _spacetx writer: https://github.com/spacetx/spacetx-writer
.. _Bio-Formats: https://www.openmicroscopy.org/bio-formats/

After converting, you can use :ref:`starfish validate<cli_validate>` to ensure that the experiment
files meet the format specifications before loading.

Your first time applying these generalized tools to convert your data can be time-consuming. If
you just want to try starfish before investing the time to format your data, you can use one of the
:ref:`formatted example datasets <datasets>` included in the starfish library.

* Tutorial: :ref:`Formatting structured data<format_structured_data>`
* Tutorial: :ref:`Formatting with TileFetcher<format_tilefetcher>`

.. _section_loading_data:

Loading Data
------------

Once the data is in :ref:`SpaceTx Format<sptx_format>`, loading the whole experiment into starfish
is simple. The only options are for selecting which :term:`FOVs <Field of View (FOV)>` and
subsets to load into memory.

As mentioned in the previous section, it is also possible to
:ref:`directly load data <tilefetcher_loader>` that has not been formatted, although
there may be performance implications in doing so. This method is also more complicated.

* Tutorial: :ref:`Loading SpaceTx Formatted Data <loading_data>`
* Tutorial: :ref:`Loading Data Without Formatting <tilefetcher_loader>`

.. _section_manipulating_images:

Manipulating Images
-------------------

Sometimes it can be useful subset the images by, for example, excluding out-of-focus images or
cropping out edge effects. For sparse data, it can be useful to project the z-volume into a single
image, as this produces a much faster processing routine. Starfish supports the cropping and
projecting of :py:class:`.ImageStack`\s with the :py:meth:`.sel` and :py:meth:`.reduce` methods.

* Tutorial: :ref:`Cropping <tutorial_cropping>`
* Tutorial: :ref:`Projecting <tutorial_projection>`

.. _section_correcting_images:

Correcting Images
-----------------

These stages are typically specific to the microscope, camera, filters, chemistry, and any tissue
handling or microfluidices that are involved in capturing the images. These steps are typically
*independent* of the assay. *Starfish* enables the user to design a pipeline that matches their
imaging system and provides some basic image correction methods.

* Tutorial: :ref:`Illumination Correction <tutorial_illumination_correction>`
* Tutorial: :ref:`Image Registration <tutorial_image_registration>`

.. _section_improving_snr:

Enhancing Signal & Removing Background Noise
--------------------------------------------

These stages are usually specific to the sample being analyzed. For example, tissues often have
some level of autofluorescence which causes cellular compartments to have more background noise than
intracellular regions. This can confound spot finders, which look for local intensity differences.
These approaches ameliorate these problems.

* Tutorial: :ref:`Removing Autofluorescence <tutorial_removing_autoflourescence>`

.. _section_normalizing_intensities:

Normalizing Intensities
-----------------------

Most assays are designed such that intensities need to be compared between :term:`rounds<Imaging
Round>` and/or :term:`channels<Channel>` in order to :term:`decode<Decoding>` spots. As a basic
example, smFISH spots are labeled by the channel with the highest intensity value. But because
different channels use different fluorophores, excitation sources, etc. the images have different
ranges of intensity values. The background intensity values in one channel might be as high as
the signal intensity values of another channel. Normalizing image intensities corrects for these
differences and allows comparisons to be made.

Whether to normalize
^^^^^^^^^^^^^^^^^^^^

The decision of whether to normalize depends on your data and decoding method used in the next
step of the pipeline.
If your :py:class:`.ImageStack` has approximately the same
range of intensities across rounds and
channels then normalizing may have a trivial effect on pixel values. Starfish provides utility
functions :ref:`imshow_plane<tutorial_imshow_plane>` and
:ref:`intensity_histogram<tutorial_intensity_histogram>` to visualize images and their intensity
distributions.

Accurately normalized images is important if you plan to decode features with
:py:class:`.MetricDistance` or :py:class:`.PixelSpotDecoder`. These two algorithms use the
:term:`feature trace<Feature (Spot, Pixel) Trace>` to construct a vector whose distance from
other vectors is used decode the feature. Poorly normalized images with some systematic or random
variation in intensity will bias the results of decoding.

However if you decode with :py:class:`.PerRoundMaxChannel`, which only compares intensities
between channels of the same round, precise normalization is not necessary. As long the intensity
values of signal in all three channels are greater than background in all three channels the
features will be decoded correctly.

How to normalize
^^^^^^^^^^^^^^^^

How to normalize depends on your data and a key assumption. There are two approaches for
normalizing images in starfish:

Normalizing Intensity Distributions
"""""""""""""""""""""""""""""""""""

If you know a priori that image volumes acquired for every channel and/or every round should have
the same distribution of intensities then the intensity *distributions* of image volumes can be
normalized with :py:class:`.MatchHistograms`. Typically this means the number of spots and amount of
background autofluorescence in every image volume is approximately uniform across channels and/or
rounds.

* Tutorial: :ref:`Normalizing Intensity Distributions<tutorial_normalizing_intensity_distributions>`

Normalizing Intensity Values
""""""""""""""""""""""""""""

In most data sets the differences in gene expression leads to too much variation in number of
spots between channels and rounds. Normalizing intensity distributions would incorrectly skew the
intensities. Instead you can use :py:class:`.Clip`, :py:class:`.ClipPercentileToZero`, and
:py:class:`.ClipValueToZero` to normalize intensity *values* by clipping extreme values and
rescaling.

* Tutorial: :ref:`Normalizing Intensity Values <tutorial_normalizing_intensity_values>`

.. _section_finding_and_decoding:

Finding and Decoding Spots
--------------------------

Finding and decoding bright spots is the unique core functionality of starfish and is necessary in
every image-based transcriptomics processing pipeline. The inputs are all the images from a
:term:`FOV <Field of View (FOV)>` along with a :term:`codebook <Codebook>` that describes the
experimental
design. The output after decoding is a :term:`DecodedIntensityTable` that contains the
location, intensity values, and mapped :term:`target <Target>` of every detected
:term:`feature <Feature>`.

Every assay uses a set of rules that the :term:`codewords <Codeword>` in the codebook
must follow (e.g. each target has one hot channel in each round). These rules determine which
decoding methods in starfish should be used. See :ref:`section_which_decoding_approach` to
learn about different codebook designs and how to decode them.

There are two divergent decoding approaches, spot-based and pixel-based, used in the image-based
transcriptomics community when it comes to analyzing spots in images:

.. image:: /_static/design/decoding_flowchart.png
   :scale: 50 %
   :alt: Decoding Flowchart
   :align: center

Spot-Based Decoding
^^^^^^^^^^^^^^^^^^^

The spot-based approach finds spots in each image volume based on the brightness of regions
relative to their surroundings and then builds a :term:`spot trace<Feature (Spot, Pixel) Trace>`
using the appropriate :ref:`TraceBuildingStrategies<howto_tracebuildingstrategies>`. The spot
traces can then be mapped, or *decoded*, to codewords in the codebook using a
:py:class:`.DecodeSpotsAlgorithm`.

.. list-table::
   :widths: auto
   :header-rows: 1

   * - When to Use
     - How-To
   * - Images are amenable to spot
       detection methods
     - :ref:`howto_spotfindingresults`
   * - Data is from sequential methods
       like smFISH
     - :ref:`howto_simplelookupdecoder`
   * - Spots are sparse and may not be
       aligned across all rounds
     - :ref:`Use TraceBuildingStrategies.NEAREST_NEIGHBOR <howto_tracebuildingstrategies>`

* Tutorial: :ref:`Spot-Based Decoding with FindSpots and DecodeSpots <tutorial_spot_based_decoding>`

Pixel-Based Decoding
^^^^^^^^^^^^^^^^^^^^

The pixel-based approach first treats every pixel as a :term:`feature <Feature>` and constructs a
corresponding :term:`pixel trace<Feature (Spot, Pixel) Trace>` that is mapped to codewords.
Connected component analysis is then used to label connected pixels with the same codeword as an RNA
spot.

* Tutorial: :ref:`Pixel-Based Decoding with DetectPixels <tutorial_pixel_based_decoding>`

.. _section_which_decoding_approach:

What Decoding Pipeline Should I Use?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are unsure which spot finding and decoding methods are compatible with your data here is a
handy table that summarizes the three major :term:`codebook <Codebook>` designs and what methods
can be used to decode each of them. If your codebook doesn't fall into any of these categories,
`make a feature request on github <https://github.com/spacetx/starfish/issues/new/choose>`_, we
would love to hear about unique codebook designs!

.. _tab-codebook-designs:

.. table::
   :class: "table-bordered"

   +-----------------+---------------------------+-------------------------+--------------------------+
   | Name            | Linearly Multiplexed      | One Hot Exponentially   | Exponentially Multiplexed|
   |                 |                           | Multiplexed             |                          |
   +=================+===========================+=========================+==========================+
   | Assays          | - sequential smFISH       | - In Situ Sequencing    | - MERFISH                |
   |                 | - RNAscope                | - seqFISH               | - DARTFISH               |
   |                 | - osmFISH                 | - FISSEQ                | - seqFISH+               |
   |                 |                           | - STARmap               |                          |
   |                 |                           | - BaristaSeq            |                          |
   +-----------------+---------------------------+-------------------------+--------------------------+
   | Example 7-round | |linear1|                 | |onehot1|               | |multiplex1|             |
   | Codebook        |                           |                         |                          |
   | Diagrams        | |linear2|                 | |onehot2|               | |multiplex2|             |
   +-----------------+---------------------------+-------------------------+--------------------------+
   | Description     | Codewords have only one   | Codewords are one hot   | Each codeword is a       |
   |                 | round and channel with    | in each round           | combination of signals   |
   |                 | signal                    |                         | over multiple rounds     |
   +-----------------+---------------------------+-------------------------+--------------------------+
   | Reference Image | No                        | Yes                     | Yes                      |
   | Needed?         |                           |                         |                          |
   +-----------------+---------------------------+-------------------------+--------------------------+
   | starfish        | - SimpleLookup            | - Exact_Match or        | - Pixel-based            |
   | Pipeline        | - Sequential +            |   Nearest_Neighbor      | - Exact_Match +          |
   | Options         |   PerRoundMaxChannel      | - PerRoundMaxChannel or |   MetricDistance         |
   |                 |                           |   MetricDistance        | - Nearest_Neighbor +     |
   |                 |                           |                         |   MetricDistance         |
   +-----------------+---------------------------+-------------------------+--------------------------+

.. |linear1| image:: /_static/design/linear_codebook_1.png
   :scale: 10%
   :align: middle
.. |linear2| image:: /_static/design/linear_codebook_2.png
   :scale: 10%
   :align: middle
.. |onehot1| image:: /_static/design/onehot_codebook_1.png
   :scale: 10%
   :align: middle
.. |onehot2| image:: /_static/design/onehot_codebook_2.png
   :scale: 10%
   :align: middle
.. |multiplex1| image:: /_static/design/multiplex_codebook_1.png
   :scale: 10%
   :align: middle
.. |multiplex2| image:: /_static/design/multiplex_codebook_2.png
   :scale: 10%
   :align: middle

.. _section_segmenting_cells:

Segmenting Cells
----------------

Unlike single-cell RNA sequencing, image-based transcriptomics methods do not physically separate
cells before acquiring RNA information. Therefore, in order to characterize cells, the RNA must be
assigned into single cells by partitioning the image volume. Accurate unsupervised cell-segmentation
is an `open problem <https://www.kaggle.com/c/data-science-bowl-2018>`_ for all biomedical imaging
disciplines ranging from digital pathology to neuroscience.

The challenge of segmenting cells depends on the structural complexity of the sample and quality
of images available. For example, a sparse cell mono-layer with a strong cytosol stain would be
trivial to segment but a dense heterogeneous population of cells in 3D tissue with only a DAPI stain
can be impossible to segment perfectly. On the experimental side, selecting good cell stains and
acquiring images with low background will make segmenting a more tractable task.

There are many approaches for segmenting cells from image-based transcriptomics assays. Below are
a few methods that are implemented or integrated with starfish to output a
:py:class:`.BinaryMaskCollection`, which represents a collection of labeled objects. If you do not
know which segmentation method to use, a safe bet is to start with thresholding and watershed. On
the other hand, if you can afford to manually define :term:`ROI <Region of Interest (ROI)>` masks
there is no better way to guarantee accurate segmentation.

.. note::
    While there is no "ground truth" for cell segmentation, the closest approximation is manual
    segmentation by an expert in the tissue of interest.

Thresholding and Watershed
^^^^^^^^^^^^^^^^^^^^^^^^^^

The traditional method for segmenting cells in fluorescence microscopy images is to threshold the
image into foreground pixels and background pixels and then label connected foreground as
individual cells. Common issues that affect thresholding such as background noise can be corrected
by preprocessing images before thresholding and filtering connected components after. There are
`many automated image thresholding algorithms <https://imagej.net/Thresholding>`_ but currently
starfish requires manually selecting a global threshold value in :py:class:`.ThresholdBinarize`.

When overlapping cells are labeled as one connected component, they are typically segmented by
using a `distance transformation followed by the watershed algorithm <https://www.mathworks
.com/company/newsletters/articles/the-watershed-transform-strategies-for-image-segmentation
.html>`_. Watershed is a classic image processing algorithm for separating objects in images and
can be applied to all types of images. Pairing it with a distance transform is particularly
useful for segmenting convex shapes like cells.

A segmentation pipeline that consists of thresholding, connected component analysis, and watershed
is simple and fast to implement but its accuracy is highly dependent on image quality.
The signal-to-noise ratio of the cell stain must be high enough for minimal errors after
thresholding and binary operations. And the nuclei or cell shapes must be convex to meet the
assumptions of the distance transform or else it will over-segment. Starfish includes the basic
functions to build a watershed segmentation pipeline and a predefined :py:class:`.Watershed`
segmentation class that uses the :term:`primary images<Primary Images>` as the cell stain.

* Tutorial: :ref:`Ways to segment by thresholding and watershed<tutorial_watershed_segmentation>`

Manually Defining Cells
^^^^^^^^^^^^^^^^^^^^^^^

The most accurate but time-consuming approach is to manually segment images using a tool such as
`ROI manager <https://imagej.net/docs/guide/146-30.html#fig:The-ROI-Manager>`_ in FIJI (ImageJ). It
is a straightforward process that starfish supports by importing
:term:`ROI <Region of Interest (ROI)>` sets stored in ZIP archives to be imported as a
:py:class:`.BinaryMaskCollection`. These masks can then be integrated into the pipeline for
visualization and assigning spots to cells.

* Tutorial: :ref:`Loading ImageJ ROI set<tutorial_manual_segmentation>`

Machine-Learning Methods
^^^^^^^^^^^^^^^^^^^^^^^^

Besides the two classic cell segmentation approaches mentioned above, there are machine-learning
methods that aim to replicate the accuracy of manual cell segmentation while reducing the labor
required. Machine-learning algorithms for segmentation are continually improving but there is no
perfect solution for all image types yet. These methods require training data (e.g. stained
images with manually defined labels) to train a model to predict cell or nuclei locations in test
data. There are `exceptions that don't require training on your specific data <http://www.cellpose
.org/>`_ but generally training the model is something to consider when evaluating how much time
each segmentation approach will require.

Starfish currently has built-in functionality to support `ilastik <https://www.ilastik.org/>`_, a
segmentation toolkit that leverages machine-learning. Ilastik has a Pixel Classification
workflow that performs semantic segmentation of the image, returning probability maps for each
label such as cells and background. To transform the images of pixel probabilities to binary
masks, you can use the same thresholding and watershed methods in starfish that are used for
segmenting images of stained cells.

* Tutorial: :ref:`Using ilastik in starfish<tutorial_ilastik_segmentation>`

.. _section_assigning_spots:

Assigning Spots to Cells
------------------------

After segmenting images to find cell boundaries, RNA spots in the :py:class:`.DecodedIntensityTable`
can be assigned to cells and then the table can be reorganized to create a single cell gene
:py:class:`.ExpressionMatrix`. These matrices are the data structure most often generated and used
by single-cell RNAseq analysis packages (e.g. `scanpy <https://icb-scanpy.readthedocs-hosted
.com/en/stable/>`_) to cluster and classify cell types. Compared to single-cell RNAseq, image-based
transcriptomic methods provide additional information about the cell, such as its location, size,
and morphology. The :py:class:`.ExpressionMatrix` holds both the 2-Dimensional matrix and cell
metadata produced by these image-based methods. This data is what links the histological context of
single cells to their transcriptomes.

In a starfish pipeline, the first step to creating a gene expression matrix is assigning spots,
aka :term:`features <Feature>`, to cells defined in a :py:class:`.BinaryMaskCollection` as cell
masks. This is done by using :py:class:`.Label` to label features with ``cell_id``\s. Currently,
:py:class:`.Label` assumes every cell mask created by
:ref:`cell segmentation<section_segmenting_cells>` encompasses a whole cell. RNA spots
with spatial coordinates that are within a cell mask are assigned to that cell and spots that do
not fall within any cell mask are not assigned a ``cell_id``. Therefore, the accuracy and
percent yield of assigned spots is largely dependent on the quality and completeness of cell
segmentation.

For data without well segmented cells, such as when no cell stain images are available, there is
potential for more sophisticated methods to assign spots to cells. For example, there are a
number of segmentation-free approaches for grouping spots into cells that starfish would like to
support in the `future <https://github.com/spacetx/starfish/issues/1675>`_.

* Tutorial: :ref:`tutorial_assigning_spots`


.. _section_working_with_starfish_outputs:

Working with starfish outputs
-----------------------------

Once you've processed your data with starfish, you are ready to load the output
files into tools like Seurat and ScanPy for further analysis. Starfish lets you
save expression matrices and segmentation masks in a variety of data formats.

* Tutorial: :ref:`working_with_starfish_outputs`

.. _section_processing_at_scale:

Processing at Scale with AWS
----------------------------

When you are ready to scale up your analysis to a full experiment or multiple
experiments, starfish can be deployed on the cloud or an an institutional
high performance computing cluster for efficient analysis of large datasets.
Implementation details will vary based on the compute resources at your disposal,
but below we demonstrate how you can analysis a full dataset on AWS.

* Tutorial: :ref:`processing_at_scale`

.. _section_further_reading:

Further Reading
---------------

Additional resources are available in :ref:`help and reference`.



.. toctree::
   :hidden:

   working_with_starfish_outputs/index
   processing_at_scale/index
.. _processing_at_scale:

Processing With AWS
===================

This tutorial will walk you through how to create two AWS Batch jobs that will use your existing starfish pipeline and
apply it in parallel to your entire experiment. Before we begin make sure you've completed the following prerequisites:

Prerequisites
-------------
- Download the template files and script needed to set up and run your aws job :download:`here </_static/starfish-aws-templates.zip>`
- Create an aws account with access to the console `Create Account <https://aws.amazon.com/premiumsupport/knowledge-center/create-and-activate-aws-account/>`__
- Create a starfish pipeline for processing a singe field of view from :ref:`ImageStack` to :ref:`DecodedIntensityTable`
- Convert your dataset to SpaceTx format. :ref:`section_formatting_data`
- Make sure you have the awscli installed ``pip install awscli``

For this tutorial we will be working with a 15 field of view `ISS dataset <https://s3.amazonaws.com/spacetx.starfish.data.public/browse/formatted/iss/20190506/experiment.json>`_

Alright let's begin!

Set up your data
-----------------

Upload your Data
+++++++++++++++++

If your dataset is already uploaded to an s3 bucket you can skip this section.

Open the aws console and navigate to `s3 <https://console.aws.amazon.com/s3/home>`_.  Select `Create Bucket` and enter the name of your bucket.
We name ours `aws-processing-example`. Once your bucket is created navigate into it and create a new folder for your spacetx-dataset. We call ours
`iss-spacetx-formatted`:

.. figure:: /_static/images/aws-tutorial-figure1.png
   :align: center

We also create a folder to hold the results from our processing called `iss-spacetx-formatted-results`.

Now that our buckets are created it's time to sync our dataset to it. Open up terminal and navigate to the folder with your spacetx formatted dataset. Run the following command:

``aws s3 sync . s3://<PATH_TO_DATASET_FOLDER>``

Ours looks like this:

``aws s3 sync . s3://aws-processing-example/iss-spacetx-formatted/``

To test that our experiment has been properly uploaded we try loading it up with starfish:


>>> from starfish import Experiment
>>> e = Experiment.from_json("https://s3.amazonaws.com/aws-processing-example/iss-spacetx-formatted/experiment.json")
>>> e

::

    {fov_000: <starfish.FieldOfView>
      Primary Image: <slicedimage.TileSet (z: 1, c: 4, r: 4, x: 1390, y: 1044)>
      Auxiliary Images:
        nuclei: <slicedimage.TileSet (z: 1, c: 1, r: 1, x: 1390, y: 1044)>
        dots: <slicedimage.TileSet (z: 1, c: 1, r: 1, x: 1390, y: 1044)>
    fov_001: <starfish.FieldOfView>
      Primary Image: <slicedimage.TileSet (z: 1, c: 4, r: 4, x: 1390, y: 1044)>
      Auxiliary Images:
        nuclei: <slicedimage.TileSet (z: 1, c: 1, r: 1, x: 1390, y: 1044)>
        dots: <slicedimage.TileSet (z: 1, c: 1, r: 1, x: 1390, y: 1044)>
    fov_002: <starfish.FieldOfView>
      Primary Image: <slicedimage.TileSet (z: 1, c: 4, r: 4, x: 1390, y: 1044)>
      Auxiliary Images:
        nuclei: <slicedimage.TileSet (z: 1, c: 1, r: 1, x: 1390, y: 1044)>
        dots: <slicedimage.TileSet (z: 1, c: 1, r: 1, x: 1390, y: 1044)>
    fov_003: <starfish.FieldOfView>
      Primary Image: <slicedimage.TileSet (z: 1, c: 4, r: 4, x: 1390, y: 1044)>
      Auxiliary Images:
        nuclei: <slicedimage.TileSet (z: 1, c: 1, r: 1, x: 1390, y: 1044)>
        dots: <slicedimage.TileSet (z: 1, c: 1, r: 1, x: 1390, y: 1044)>
      ...,
    }


Looks good! Let's move on.


Create a Recipe File
+++++++++++++++++++++

From the template files downloaded open the file named `recipe.py`. The file should include the method signature ``def process_fov(fov: FieldOfView, codebook: Codebook) -> DecodedIntensityTable:``.
Within this method add your starfish pipeline code for processing a single field of view. The return value should be a :ref:`DecodedIntensityTable`. Here's what our recipe.py file looks like.

.. code-block:: python

    from starfish import Codebook, DecodedIntensityTable, FieldOfView
    from starfish.image import ApplyTransform, Filter, LearnTransform
    from starfish.spots import DecodeSpots, FindSpots
    from starfish.types import Axes, FunctionSource


    def process_fov(fov: FieldOfView, codebook: Codebook) -> DecodedIntensityTable:
        """Process a single field of view of ISS data
        Parameters
        ----------
        fov : FieldOfView
            the field of view to process
        codebook : Codebook
            the Codebook to use for decoding

        Returns
        -------
        DecodedSpots :
            tabular object containing the locations of detected spots.
        """

        # note the structure of the 5D tensor containing the raw imaging data
        imgs = fov.get_image(FieldOfView.PRIMARY_IMAGES)
        dots = fov.get_image("dots")
        nuclei = fov.get_image("nuclei")

        print("Learning Transform")
        learn_translation = LearnTransform.Translation(reference_stack=dots, axes=Axes.ROUND, upsampling=1000)
        transforms_list = learn_translation.run(imgs.reduce({Axes.CH, Axes.ZPLANE}, func="max"))

        print("Applying transform")
        warp = ApplyTransform.Warp()
        registered_imgs = warp.run(imgs, transforms_list=transforms_list, verbose=True)

        print("Filter WhiteTophat")
        filt = Filter.WhiteTophat(masking_radius=15, is_volume=False)

        filtered_imgs = filt.run(registered_imgs, verbose=True)
        filt.run(dots, verbose=True, in_place=True)
        filt.run(nuclei, verbose=True, in_place=True)

        print("Detecting")
        detector = FindSpots.BlobDetector(
            min_sigma=1,
            max_sigma=10,
            num_sigma=30,
            threshold=0.01,
            measurement_type='mean',
        )
        dots_max = dots.reduce((Axes.ROUND, Axes.ZPLANE), func="max", module=FunctionSource.np)
        spots = detector.run(image_stack=filtered_imgs, reference_image=dots_max)

        print("Decoding")
        decoder = DecodeSpots.PerRoundMaxChannel(codebook=codebook)
        decoded = decoder.run(spots=spots)
        return decoded

Upload your recipe to s3. To make things easy we upload our recipe file to the same directory our experiment dataset lives in.

``aws s3 cp recipe.py s3://aws-processing-example/iss-spacetx-formatted/``

Set up your Batch Jobs
----------------------

So now we have our data and recipe uploaded and ready to go in s3, let's move on to actually creating our processing jobs.
Our final workflow will be composed of two jobs:

- Process each Field of View in parallel using an AWS Batch Array Job
- Combine the results from each Field of View into one large DecodedIntensityTable using an AWS Batch Job


Create a custom IAM Role
+++++++++++++++++++++++++

Before we can register our jobs we need to set up an IAM role that has access to AWSBatchServices and
our newly created s3 bucket. Navigate to the `IAM console <https://console.aws.amazon.com/iam/home>`_ an select *Roles* from the left panel.
Click *Create Role* we've called ours `spacetx-batch-uploader`. From the list of available services to prevision your role with select *batch*. Then click through the rest of the wizard
using the default settings and create the role.

We also need to give this role read and write access to our newly created s3 bucket. To do this we make a new policy and attach it to the `spacetx-batch-uploader` role.

Select *policies* from the left hand panel and click *create policy*. Click on the JSON editor and paste in the following code:


::

    {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Sid": "ListObjectsInBucket",
                "Effect": "Allow",
                "Action": [
                    "s3:ListBucket"
                ],
                "Resource": [
                    "arn:aws:s3:::<YOUR BUCKET>"
                ]
            },
            {
                "Sid": "AllObjectActions",
                "Effect": "Allow",
                "Action": "s3:*Object",
                "Resource": [
                    "arn:aws:s3:::<YOUR BUCKET>/*"
                ]
            }
        ]
    }


Here's what ours looks like:

::

    {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Sid": "ListObjectsInBucket",
                "Effect": "Allow",
                "Action": [
                    "s3:ListBucket"
                ],
                "Resource": [
                    "arn:aws:s3:::aws-processing-example"
                ]
            },
            {
                "Sid": "AllObjectActions",
                "Effect": "Allow",
                "Action": "s3:*Object",
                "Resource": [
                    "arn:aws:s3:::aws-processing-example/*"
                ]
            }
        ]
    }


Name your policy and save it, we named our `spacetx-batch-uploader`. Now navigate back to your new role and attach your s3 uploader policy.
Our `spacetx-batch` role summery now looks like this:

.. figure:: /_static/images/aws-tutorial-figure2.png
   :align: center

Note the ARN of your new role (circled in the image). You'll need it in the next few steps.


Register your Jobs
+++++++++++++++++++

Follow the `Getting Started Guide <http://docs.aws.amazon.com/batch/latest/userguide/Batch_GetStarted.html>`_ and ensure you have a valid job queue and compute environment. For this tutorial
we used the default parameters (our job queue is still called first-run-job-queue).

Here's what out Batch Dashboard looks like:

.. figure:: /_static/images/aws-tutorial-figure5.png
   :align: center

Alright now it's time to register our batch jobs. From the template files open up the file named `register-process-fov-job.json`. This file describes a batch job that will create an ec2 instance using the docker container `spacetx/process-fov`
that processes a specified single field of view using your recipe.py file. Replace the string "ADD ARN" with the aws ARN of the role you just created in the last step. Our file looks like this:

::

    {
      "jobDefinitionName": "process-fov",
      "type": "container",
      "containerProperties": {
        "jobRoleArn": "arn:aws:iam::422553907334:role/spacetx-batch",
        "image": "spacetx/process-fov",
        "vcpus": 1,
        "memory": 2500
      }
    }

NOTE: if your starfish processing is memory expensive you can adjust the allocated memory for each created instance using the `memory` parameter.

Then from the directory where this file lives run the following command:

``aws batch submit-job --cli-input-json file://register-process-fov-job.json``

You can check that your jobs had been successfully registered by navigating to the `Job Definitions page <https://console.aws.amazon.com/batch/home>`_.

Here's what our's looks like:

.. figure:: /_static/images/aws-tutorial-figure3.png
   :align: center

Now open the file named `register-merge-job.json`. This file describes a batch job that will create an ec2 instance using the docker container `spacetx/merge-batch-job` that merges together all your processed results into
one `DecodedIntensityTable`. Again replace the string "ADD ARN" with the aws ARN of your batch processing role. Our file looks like this:

::

    {
      "jobDefinitionName": "merge-job",
      "type": "container",
      "containerProperties": {
        "jobRoleArn": "arn:aws:iam::422553907334:role/spacetx-batch",
        "image": "spacetx/merge-batch-job",
        "vcpus": 1,
        "memory": 2500
      }
    }

Then from the directory where this file lives run the following command:

``aws batch submit-job --cli-input-json file://register-merge-job.json``

Again, check that your job has been successfully registered from the job console, our two jobs are ready to go!

.. figure:: /_static/images/aws-tutorial-figure4.png
   :align: center


Run your Batch Jobs
-------------------

Now that we've set everything up it's time to run our jobs! The script `starfish-workflow.py` will handle submitting the process-fov array job
then the merge job with a dependency on the first job to finish. All you'll need to do is run the script with a few parameters:

::

    --experiment-url: The path to your experiment.json file. Our is "s3://aws-processing-example/iss-spacetx-formatted/experiment.json"

    --num-fovs: The number of fields of view in the experiment. We have 15

    --recipe-location: The path to your recipe file in s3. Ours is "s3://aws-processing-example/aws-processing-example/iss-spacetx-formatted/recipe.py"

    --results-location: The s3 bucket to copy the results from the job to. Ours is "s3://aws-processing-example/iss-spacetx-formatted-results/"

    --job-queue: The name of your job queue to run your jobs. Ours is "first-run-job-queue"


Now we run our script:

::

    $ python3 starfish-workflow.py \
    >     --experiment-url "s3://aws-processing-example/iss-spacetx-formatted/experiment.json" \
    >     --num-fovs 15 \
    >     --recipe-location "s3://aws-processing-example/aws-processing-example/iss-spacetx-formatted/recipe.py" \
    >     --results-bucket "s3://aws-processing-example/iss-spacetx-formatted-results/" \
    >     --job-queue "first-run-job-queue"
    Process fovs array job 39a13edd-8cca-4e7e-9379-aa3cf757c72e successfully submitted.
    Merge results job ac5d49f5-a12e-4176-96e2-f697c6cf0a12 successfully submitted.

To monitor the status of both jobs navigate to the `AWS Batch Dashboard <https://console.aws.amazon.com/batch/home>`_. You should see 2
jobs under PENDING

.. figure:: /_static/images/aws-tutorial-figure6.png
   :align: center

From here you should be able to click on the jobs and track their movement through the RUNNABLE -> RUNNING -> SUCCEEDED states.
NOTE: Batch jobs may take up to 10 minutes to move from PENDING to RUNNABLE. When both jobs have reached the SUCCEEDED state check
that everything worked by navigating to your results bucket. The bucket should include the processed results from
each field of view as well as the concatenated results called `merged_decoded_fovs.nc`. Here's what our bucket contains:

.. figure:: /_static/images/aws-tutorial-figure7.png
   :align: center

And that's it! You have successfully set up and processed your experiment using aws. As long as you keep your job definitions you can rerun the jobs
using the same command anytime.
.. _working_with_starfish_outputs:

Working with Starfish Outputs
-----------------------------

Starfish's output_formats are serialized as :code:`netcdf` or :code:`csv` files. These files are
easy to work with in both Python_ and R_.

.. _Python: https://www.python.org/
.. _R: https://www.r-project.org/about.html

To work with the :py:class:`IntensityTable` in python, it's as simple as using that object's open_netcdf
command:

.. code-block:: python

    import starfish
    example_netcdf_file: str = "docs/source/_static/example_data_files/decoded.nc"
    intensity_table: starfish.IntensityTable = starfish.IntensityTable.open_netcdf(example_netcdf_file)
    print(intensity_table)

in R_, the ncdf4 library allows the :code:`.nc` archive, which is based on hdf5, to be opened.
It will contain a number of variables, each of which can be accessed by name. Alternative
installation instructions can be accessed here. Alternative installation instructions can be
accessed here_:

.. _here: http://cirrus.ucsd.edu/~pierce/ncdf/

.. code-block:: R

    install.packages("ncdf4")
    library("ncdf4")
    example_netcdf_file <- "docs/source/_static/example_data_files/decoded.nc"
    netcdf_connection <- nc_open(example_netcdf_file)

    # access the z-coordinate vector
    zc <- ncvar_get(netcdf_connection, "zc")
    head(zc)

    # access the 3-dimensional data structure containing intensty information
    # this variable has a special name, the rest are accessible with the string
    # constants you would expect from the starfish python API.
    data <- ncvar_get(netcdf_connection, "__xarray_dataarray_variable__")
    head(data)

To work with the decoded table is even simpler, as they are stored as :code:`.csv` files, and can
be read natively by pandas in Python and natively in R.

Python:

.. code-block:: python

    import pandas as pd
    example_decoded_spots_file: str = "docs/source/_static/example_data_files/decoded.csv"
    table: pd.DataFrame = pd.read_csv(example_decoded_spots_file, index_col=0)
    table.head()

R:

.. code-block:: R

    example_decoded_spots_file <- "docs/source/_static/example_data_files/decoded.csv"
    table <- read.csv(file=example_decoded_spots_file, header=TRUE, sep=',', row.names=1)
    head(table)

Output Formats
^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1
   :caption: Contents:

.. toctree::
    IntensityTable/index.rst

.. toctree::
    ExpressionMatrix/index.rst

.. toctree::
    DecodedSpots/index.rst

.. toctree::
    SegmentationMask/index.rst.. _ExpressionMatrixSpecification:

ExpressionMatrix
================
The :py:class:`ExpressionMatrix` is a 2-dimensional :code:`cells (x)` by :code:`genes (y)` array
whose values contain the expression of a gene in a particular cell. The :py:class:`ExpressionMatrix`
is additionally annotated with the :code:`x, y, z` pixel coordinates of the centroid of the cell in
a pixel space, and :code:`xc, yc, zc` in physical coordinate space. Additional metadata may be added
at the user's convenience to either the cells or genes.

Data
----
Gene expression are stored as numeric values, typically as integers in Image-based transcriptomics
experiments, since they represent counted fluorescent spots, each corresponding to a single detected
RNA molecule.

Metadata
--------

:code:`cells; cell_id (int):` cell identifier

:code:`cells; x, y, z (int):` coordinates of cell centroid in pixel space

:code:`cells; xc, yc, zc (int):` coordinates of cell centroid in global coordinate space (um)

:code:`genes; gene_id (int):` GENCODE gene ID

:code:`genes; gene_name (int):` Human-readable gene symbol (e.g. HGNC gene symbol for human data)

Implementation
--------------
Starfish Implements the :py:class:`ExpressionMatrix` as an :code:`xarray.DataArray` object to take
advantage of `xarray's`_ high performance, flexible metadata storage capabilities, and serialization
options

.. _`xarray's`: http://xarray.pydata.org/en/stable/

Serialization
-------------
The :py:class:`ExpressionMatrix` can leverage any of the :code:`xarray` serialization features,
including csv, zarr, and netcdf. We choose netcdf as it currently has the strongest support and
interoperability between R and python. Users can load and manipulate :py:class:`ExpressionMatrix`
using R by loading them with the `ncdf4`_ package. In the future, NetCDF serialization may be
deprecated if R gains Zarr support.

.. _ncdf4: https://cran.r-project.org/web/packages/ncdf4/index.html
.. _SegmentationMask:

SegmentationMask
================
The :py:class:`SegmentationMask` is an integer array the same size as the :code:`x-y` plane of the
input :py:class:`ImageStack`. The values of the :py:class:`SegmentationMask` indicate which cell
each pixel corresponds to. At the moment, only hard assignment is supported.

Implementation
--------------
The :py:class:`SegmentationMask` is currently implemented as a :code:`numpy.ndarray` object. We are
revisiting the optimal object for this information, and are aware of the pending need to use
segmented instances of cells to label points in 3d, and also to support non-cell instances. The
current implementation is expected to change.

Serialization
-------------
numpy arrays can be saved in a variety of formats of which the most common is npy_

.. _npy: https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.save.htmlIntensityTable
==============

The :py:class:`IntensityTable` summarizes information about RNA features that are detected in
image-based transcriptomics experiments. The most common features are spots, pixels, or connected
pixels that result from the fluorescence of a probe that has been experimentally attached to an RNA
molecule.

To gather information on the RNA present in a tissue slice, image-based transcriptomics experiments
require the imaging of the same RNA molecules over multiple rounds and across multiple fluorescence
channels. In multiplex assays, the identity of these spots is only determined after measuring the
spot's intensity across all rounds and channels. In sequential assays, each round and channel
identifies all detected RNA molecules for a specific target.

The :py:class:`IntensityTable` stores a summary of the intensity of each feature across each round
and channel, forming a 3-dimensional array with dimensions (feature, round, channel). It also stores
metadata that declare how the features were measured and are defined, and metadata for each
individual feature. Finally, it supports addition of arbitrary metadata for each feature, round,
channel, or on the detection of those features.

Data
----

IntensityTables are stored as three dimensional arrays of shape :code:`(f, r, c)`, where :code:`f`
is the number of detected features in the experiment, :code:`r` is the number of rounds, and
:code:`c` is the number of channels. Entries in the array are floats in the range of :code:`[0, 1]`,
where :code:`1` represents the maximum intensity and :code:`0`, the minimum intensity.

Table Metadata
--------------
The IntensityTable stores some standard metadata that record how features are summarized.

:code:`intensity_measurement_type (string):`

In the case where features are composed of spots or connected pixels, the intensity of each pixel
in the feature must be integrated into a single measurement. The :code:`intensity_measurement_type`
field stores the calculation made to carry out this integration. This entry is filled by the spot
detector. The default value is :code:`max` but other common methods might include :code:`mean` or
:code:`median`.

:code:`area_measurment_type (string):`

This size of each feature can be summarize in different ways that often depend on how they're
detected. For example, it is natural to summarize spots as circles, while pixel-based
approaches can produce irregular features that are difficult to summarize parametrically. This
entry is filled by the spot detector and describes how areas are calculated.


Feature Metadata
----------------
:py:class:`IntensityTable` Features are annotated with several types of metadata. These include the
:code:`x, y, z` position of each feature in pixel and coordinate space and the area of the feature.
Additionally, there are entries that can be filled by a decoder which determines the target that a
spot corresponds to, or a Segmentation mask, which determines which cell an object corresponds to.
Users may add additional metadata columns with names that do not overlap this list as they see fit.
The following metadata are always present on IntensityTables. In some cases, these values may be
set to NaN or null if the IntensityTable has not been decoded or merged with a segmentation result.

:code:`x, y, z (int):` coordinates of the centroid of the feature in pixel space

:code:`xc, yc, zc (float):` coordinates of the centroid of the feature in global coordinates (um)

:code:`area (float):` area of the spot in pixels  # TODO: should this be in um^2?

:code:`cell (int):` the id for the cell that this feature belongs to

:code:`gene (str):` gene name of target RNA, (gene symbol or GENCODE gene ID)

Implementation
--------------
Starfish Implements the :py:class:`ExpressionMatrix` as an :code:`xarray.DataArray` object to take
advantage of `xarray's`_ high performance, flexible metadata storage capabilities, and Serialization
options

.. _`xarray's`: http://xarray.pydata.org/en/stable/

Serialization
-------------
The :py:class:`IntensityTable` can leverage any of the :code:`xarray` serialization features,
including csv, zarr, and netcdf. We choose netcdf as it currently has the strongest support and
interoperability between R and python. Users can load and manipulate :py:class:`IntensityTables`
using R by loading them with the `ncdf4`_ package. In the future, NetCDF serialization may be
deprecated if R gains Zarr support.

.. _ncdf4: https://cran.r-project.org/web/packages/ncdf4/index.html
.. _DecodedSpotsSpecification:

DecodedSpots
================
The :py:class:`DecodedSpots` is a 2-dimensional tabular data structure where each record represents
a spot, and each record contains, at minimum, columns that specify the `x` and `y` coordinates of
the spot in physical space and the genes that it targets.

Additional columns are not validated, however there are several optional columns with standard
names. If those data are available, methods that take :py:class:`DecodedSpots` objects may make
use of them, so it is in users interests to name those columns appropriately.

Required Columns
----------------

These columns must be present for an object to be constructed.

:code:`target (str):` the gene target for the spot

:code:`x (float):` the x coordinate of the spot in physical space

:code:`y (float):` the y coordinate of the spot in physical space

Optional Columns
----------------

These columns are not validated, but have special meaning to :py:class:`DecodedSpots`

:code:`z (float):` the z coordinate of the spot in physical space

:code:`target_probability (float):` the quality of the decoding, or probability that the spot is
associated with the listed target.

:code:`cell (int):` the identifier of the cell that contains this spot.

:code:`cell_probability (float):` the quality of a cell association, or probability that the spot
is associated with the listed cell.

Implementation
--------------
Starfish Implements the :py:class:`DecodedSpots` as a wrapper for the :code:`pd.DataFrame` object

Serialization
-------------
:py:class:`DecodedSpots` defines methods to save and load the spot file as a csv file.
.. _help and reference:

Help & Reference
================

.. toctree::
    :maxdepth: 1

    data_model/index
    spacetx-format/index
    available_datasets/index
    configuration/index
    glossary/glossary
    spacetx_consortium/index
    request_support/index
.. _spacetx consortium:

SpaceTx Consortium
==================

The SpaceTx consortium is a group led by Ed Lein and the Allen Institute for Brain Science that
is engaged in benchmarking image-based transcriptomics assays. In collaboration with the lead
developers of 10+ image-based transcriptomics assays, they are applying each technology to adjacent
sections of mouse and human primary visual cortex. The goal of this experiment is to understand
what assay(s) are best applied to (a) spatially localizing cells by type based on marker genes and
(b) obtaining robust phenotypes of cells by measuring RNA from thousands of genes.

Consortium members include:

- Ed Boyden, MIT
- Long Cai, Caltech
- Fei Chen, Broad Institute
- Karl Dieserroth, Stanford
- Ed Lein, Allen Institute for Brain Science
- Sten Linnarrson, Karolinska Institutet
- Joakim Lundeberg, SciLifeLab
- Jeffrey Moffit, Harvard
- Mats Nilsson, SciLifeLab
- Aviv Regev, Broad Institute
- Anthony Zador, Cold Spring Harbor Laboratory
- Kun Zhang, University of California San Diego
- Xiaowei Zhuang, Harvard

Pipelines
---------

Follow the links in the table below to see pipelines from each of the SpaceTx groups ported to starfish.
All of the example pipelines can be found in the `notebooks directory <https://github.com/spacetx/starfish/tree/master/notebooks/>`_.

====================  ==========  ===================  ==================
 Assay                Loads Data  Single-FoV Pipeline  Multi-FoV Pipeline
--------------------  ----------  -------------------  ------------------
 MERFISH              [x]         [x] mer_             in process
 ISS                  [x]         [x] iss_             [x]
 osmFISH              [x]         [x] osm_             [ ]
 smFISH               [x]         [x] 3ds_             [x]
 BaristaSeq           [x]         [x] bar_             [ ]
 DARTFISH             [x]         [x] dar_             [ ]
 ex-FISH              [x]         [ ]                  [ ]
 StarMAP              [x]         [x] str_             [ ]
 seq-FISH             [x]         in process: seq_     [ ]
 FISSEQ               no data     no pipeline          [ ]
====================  ==========  ===================  ==================

.. _mer: https://github.com/spacetx/starfish/blob/master/notebooks/MERFISH.ipynb
.. _iss: https://github.com/spacetx/starfish/blob/master/notebooks/ISS.ipynb
.. _osm: https://github.com/spacetx/starfish/blob/master/notebooks/osmFISH.ipynb
.. _bar: https://github.com/spacetx/starfish/blob/master/notebooks/BaristaSeq.ipynb
.. _dar: https://github.com/spacetx/starfish/blob/master/notebooks/DARTFISH.ipynb
.. _str: https://github.com/spacetx/starfish/blob/master/notebooks/STARmap.ipynb
.. _seq: https://github.com/spacetx/starfish/blob/master/notebooks/SeqFISH.ipynb
.. _3ds: https://github.com/spacetx/starfish/blob/master/notebooks/smFISH.ipynb
Available datasets
==================

Starfish makes a small number of datasets :ref:`easily accessible<datasets>`.

To load them: :code:`import starfish.data`... _cli_config:

Configuration
=============

`starfish` commands can all be configured in a unified way.

If the :ref:`STARFISH_CONFIG <env_config>` environment variable is set, then it will
be loaded, either (1) as a JSON file if the value starts with ``@`` or (2) as a JSON
string for all other values. For example:

::

    export STARFISH_CONFIG=@~/.my.starfish

::

    export STARFISH_CONFIG='{"verbose": false}'

Otherwise, ``~/.starfish/config`` will be loaded if it exists and is a valid JSON file.
If neither is true, then only the :ref:`default values <env_defaults>` will apply.


Additionally, the individual properties from the configuration JSON can be set by
environment variable. These values take precedence if also set in the configuration
file. For example:

::

 {
     "validation": {
         "strict": true
     }
 }

can also be specified as:

::

    export STARFISH_VALIDATION_STRICT=true

Other valid values for "true" are: "TRUE", "True", "yes", "y", "1", "on", and "enabled".

.. _env_defaults:

Default values
--------------

To not require any mandatory configuration, the following values will be
assumed if no configuration is available.

::

 {
     "slicedimage": {
         "caching": {
             "debug": false,
             "directory": "~/.starfish/cache",
             "size_limit": 5e9
         },
     },
     "validation": {
         "strict": false
     },
     "verbose": true
 }

.. _env_main:

Main environment variables
--------------------------

.. _env_config:

``STARFISH_CONFIG``
~~~~~~~~~~~~~~~~~~~

The primary configuration variable is ``STARFISH_CONFIG``.
By default, it is equivalent to ``~/.starfish/config`` which need not exist.

.. _env_validate_strict:

``STARFISH_VALIDATION_STRICT``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Strict validation will run the equivalent of :ref:`starfish validate <cli_validate>`
on every experiment that is loaded. Execution will halt with an exception. The initial
loading of an experiment will also take longer since files must be downloaded for
validation. If caching is enabled, the overall impact should not be significant.

.. _env_verboase:

``STARFISH_VERBOSE``
~~~~~~~~~~~~~~~~~~~~

Whether or not various commands should should print internal status messages.
By default, true.

.. _env_backend:

Backend environment variables
-----------------------------

Starfish currently uses the slicedimage library as a backend for storing large image sets.
Configuration values starting with ``SLICEDIMAGE_`` (or optionally, ``STARFISH_SLICEDIMAGE_``)
will be passed to the backend without modification.

.. _env_slicedimage_caching_size_limit:

``SLICEDIMAGE_CACHING_SIZE_LIMIT``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Maximum size of the :ref:`cache directory <env_slicedimage_caching_directory>`.
By default, size_limit is 5GB. Setting size_limit to 0 disables caching.

.. _env_slicedimage_caching_directory:

``SLICEDIMAGE_CACHING_DIRECTORY``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Directory where the cache will be stored.
By default, the directory ``~/.starfish/cache``.

.. _env_slicedimage_debug:

``SLICEDIMAGE_CACHING_DEBUG``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whether or not to print which files are being cached.
By default, false.
SpaceTx Format
--------------

.. toctree::
   :maxdepth: 1
   :caption: Contents:

.. toctree::
    SpaceTxFormat/index.rst

.. toctree::
    Validation/index.rst
.. _sptx_format:

.. mdinclude:: ../../../../../starfish/spacetx_format/README.md
.. _schema:

.. mdinclude:: ../../../../../starfish/spacetx_format/schema/README.md

.. _cli_validate:

Validation
==========

The `starfish validate` command provides a way to check that a fileset based on the
:ref:`sptx_format` is valid. One of the schema requirements is that the codebook is version 0.0.0
and the experiment is 4.0.0 or 5.0.0.

Usage
^^^^^

starfish validate --help will provide instructions on how to use the tool:

.. program-output:: env MPLBACKEND=Agg starfish validate --help


Examples
^^^^^^^^

.. code-block:: bash

    $ starfish validate experiment tmp/experiment.json > /dev/null && echo ok

Validating the experiment, validates all of the included files. These files can also be individually validated:

.. code-block:: bash

    $ starfish validate codebook tmp/codebook.json.. _data_model:

Data Model
==========

The starfish package is designed to allow users to create pipelines to count spots and
assign those spots to cells. In order to accomplish this, starfish defines a series of
abstractions, and those abstractions in turn define how we require the data to be formatted,
and how the data is processed.

Imaging experiments that are compatible with starfish typically cannot capture an entire tissue
slide within the field of view of a microscope. As a result, these experiments are often captured as
collages of sometimes-overlapping images that are stitched back into panoramas with software.

Because the complete images are often tens of thousands of pixels in :code:`x` and :code:`y`, it
can be challenging to process panoramas. In contrast, fields of view are often small enough
that they can be processed on laptops. As a result, panoramas are stored in :code:`Experiment`
objects, which are composed of all the fields of view captured by the microscope.

Note that at the current time, this means there is some up-front cost to users, who will need to
determine, using our documentation, how to reformat their microscopy data. We are working with
`Bio-Formats <bio_formats>`_ to automate this conversion and reduce this up-front cost.

Field of View
-------------

A field of view in a starfish-compatible experiment contains images of several types. Each field
of view **must** contain images containing spots. These are the "primary" images. Additional images
with the same :code:`(y, x)` size of nuclei (for segmenting cells) or fiduciary markers (for
registering images within the field of view) can also be associated with a field of view, if the
user captures them. It is important to note that all images regardless of the fluorescence channel,
time point/round, or z-plane that are taken of the specific :code:`(y, x)` coordinate area should
be included in a field of view:

.. image:: /_static/design/imagestack.png

This data structure is general across the types of data that we've observed in image-based
transcriptomics and proteomics studies, which have variables numbers of fluorescence channels
and imaging rounds. Starfish simply builds up a set of 5D tensor for each field of view, one for
each image type (primary, dots, nuclei, ...). The dimensions are :code:`round (time), channel,
z, y, x`.

.. image:: /_static/design/field_of_view.png

Processing Model
----------------

Starfish breaks up data by field of view to enable efficient processing. Unlike stitched images,
fields of view are often small enough to work with on laptops, and are efficient to pass around to
cloud compute instances or HPC systems. They're also more easily broken up into pieces that fit on
modern GPUs, which while not yet integrated in starfish, show promise for speeding early image
processing tasks.

This processing model is particularly amenable to processing on the cloud, where there is no joint
file system that all the compute instances can access. In these ecosystems, input data for pipelines
must be localized and downloaded to the machine, and uploaded back to the data store when processing
is completed. By working with small 2-dimensional planes, starfish is able to exercise granular
control over data upload, producing an efficient processing system for large imaging experiments
that works both on local HPC clusters and on the cloud.

.. image:: /_static/design/processing_model.png

Dual Coordinate Systems
-----------------------

Because starfish requires that microscopy experiments are broken up into fields of view, it is
important to keep track of the physical bounding box of each field of view. This information will be
needed to convert data into SpaceTx format, and outputs that starfish produces will track objects
in this coordinate space.

Internally, starfish treats images as tensors, and will primarily operate on pixel coordinates,
bringing the physical coordinates along for the ride. This duality will pervade the remainder of
the documentation.


Next Steps
----------

At this point, the tutorial forks. You can either dive into formatting your data in SpaceTx-Format,
or play skip forward and play with starfish using our pre-constructed example datasets.
We suggest the latter, as it will give you a sense of starfish's capabilities before you put work
into reformatting your data.
Here, we clearly define the relevant terms needed to understand the spaceTx pipeline specification

Glossary
--------

.. glossary::

    Channel
        An imaging mode that captures a continuous-valued feature from a field of view. Examples of channels include the read-out from a fluorescent dye, such as Cy3, or a the abundance of an isotope captured from a mass spectrometer.

    Imaging Round
        Several image-based transcriptomics and proteomics approaches will image the same tissue multiple times. Each time the tissue is imaged is a discrete imaging round.

    Target
        A feature that is the target of quantification by an image-based assays. Common targets include mRNA transcripts or proteins.

    IntensityTable
        An intensity Table contains the features identified in an ImageStack. It can be thought of as an array whose entries are the intensities of each feature across the imaging rounds and channels of a field of view. Starfish exposes several processing tools to decode the features of the table, estimate their qualities, and assign features to cells.

    DecodedIntensityTable
        A representation of a decoded intensity table. Contains the features identified in an ImageStack as well as their associated target values.

    Codeword
        A codeword maps expected intensities across multiple image tiles within a field of view to the target that is encoded by the codeword.

    Codebook
        A codebook contains all the codewords needed by an experiment to decode an IntensityTable. It also contains a mapping of channels to the integer indices that are used by starfish to represent them internally.

    Pipeline
        A sequence of data processing steps to process the inputs into the desired outputs.

    Pipeline Component
        A single data processing step in the pipeline, as defined by its input and output file formats, e.g., the spot-detection component takes as input an image and outputs a table of spot locations, shapes, and intensities.

    Pipeline Component Algorithm
        A specific algorithm type that adheres to the specified input and output file formats required by the component it belongs to. For example, a spot-detection component algorithm can be realized as a Gaussian blob detector or a connected components labeller. Both find spots and accept the same inputs and produce the same outputs, hence belong to the same component. However, the underlying properties of the algorithms (and parameterizations) may be quite different.

    Pipeline Specification
        A document describing the pipeline in detail, including an ordered list of all the pipeline components, and expected input/output file formats at each step of computation.

    Pipeline Implementation
        Actual code for the pipeline. This code will be packaged as a well-documented Python library and corresponding command line tool for use by consortium members to facilitate easy sharing and comparison of results across labs/methods.

    Manifest
        The data manifest is a file that includes the locations of all fields of view for either primary or auxiliary images.

    Field of View (FOV)
        A collection of Image Tiles corresponding to a specific volume or plane of the sample, under which the signal for all channels and all imaging rounds were acquired. All tiles within this FOV are the same size, but the manifest allows for different spatial coordinates for different imaging rounds or channels (to accommodate slight movement between rounds, for example).
        In microscopy, a field of view corresponds to the camera sensor mapped to the sample plane, and many such fields of view are expected to be taken per tissue slice.

    Region of Interest (ROI)
        Areas of an image identified for a particular purpose, such as to define the boundaries of a cell.

    Image Tile
        A single plane, single channel, single round 2D image. In the manifest, each tile has information about its (X,Y,Z) coordinates in space, and information about which imaging round (R) and/or fluorescence channel (C) it was acquired under.

    ImageSlice
        The image volume corresponding to a single round and single channel of the Field of View.

    Coordinates (Tile)
        Coordinates refer to the physical location of a Tile with respect to some independent reference.  If a pair of values are provided, it corresponds to the physical coordinates of the edges.  If a single value is provided, it corresponds to the center of the tile.  For x and y, two values are required.  For z, both a single value and a pair of values are valid.

    Axes (Tile)
        Each pixel is located along five axes, which are round, channel, z-plane, y, and x.

    Labels (Tile)
        Labels are the valid set of values for each axis.

    Index (Tile)
        An index indicates the label or range of labels for a given axis.  These should be a whole number (non-negative integers) or a python contiguous slice representing a range.

    Selectors (Tile)
        A mapping of axes (round, channel, and z-plane) to their respective index.  These are expressed as a mapping from Axis to index.

    Primary Images
        The primary image data for an experiment. Primary images contain information on imaging targets. primary images build fields of view that usually contain multiple channels and may contain multiple imaging rounds. Primary images can be decoded to identify the abundance of transcript or protein targets.

    Auxiliary Images
        Any user-submitted additional images for analysis beyond the primary images. These images may be of lower dimension than the primary images (e.g., single channel images), but should span the same spatial extent as the primary images acquired under the same FOV. Auxiliary images are used to aid the image processing of the primary images.

        Examples of such data may include:

        Nuclei (DAPI or similar nuclear stain): this required image shows cell nuclei and is crucial for cell segmentation further on down the pipeline.

        Dots: an image containing the locations of imaging features across a field of view.

        Other stains or labels: these optional (but recommended) image(s), including but not limited to antibody stains, may capture additional information about cell boundaries or subcellular structure that will be useful for cell segmentation and/or additional spatial analyses.

    Registration
        Refers to the process of aligning multiple images of the same spatial location, mostcommonly across multiple rounds of imaging within a FOV.

    Stitching
        The process of combining images from multiple fields of view into a larger image thatspans the extent of the sample.

    Feature
        The value of a spot (aggregated across all pixel values circumscribed by that spot) or the value of a single pixel.

    Feature (Spot, Pixel) Trace
        Feature intensity values across all imaging rounds and/or color channels. These map to codewords in a codebook.

    Decoding
        Matching putative barcodes to codewords in a codebook to read out the corresponding target believed to be associated with that barcode.

    Rolony
        A rolling-circle amplified "colony", or rolony, is an amplicon produced by image-based transcriptomics assays that use circular probes to increase signal... _request support:

Request Support
---------------

Starfish user support is hosted on the `image.sc forum <https://forum.image.sc/tag/starfish>`_ where
questions and discussion will serve as a continually expanding resource for starfish users. Please
post any help or support requests on `image.sc <https://forum.image.sc/tag/starfish>`_ with the
``starfish`` tag.

To see the code or report a bug, please visit the `github repository <https://github
.com/spacetx/starfish>`_.
Example Data Files
==================

This folder contains some example data files for use with the starfish documentation. The files are
as follows:

1. decoded.csv: A :code:`csv` file containing the :py:class:`DecodedSpots` output from the ISS
pipeline run on testing data and serialized using :py:meth:`DecodedSpots.save`.

2. decoded.nc: A :code:`netcdf` file containing the :py:class:`IntensityTable` output from the ISS
pipeline run on testing data and serialized using :py:meth:`IntensityTable.to_netcdf`.
.. _API:

API
===

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. toctree::
   data_structures/index

.. toctree::
   slicedimage/index

.. toctree::
   experiment_builder/index

.. toctree::
   image/index.rst

.. toctree::
   spots/index.rst

.. toctree::
   morphology/index.rst

.. toctree::
   types/index.rst

.. toctree::
   validation/index.rst

.. toctree::
   utils/index.rst

.. toctree::
   display/index.rst

.. toctree::
   plot_utils/index.rst

.. toctree::
   datasets/index.rst
.. _experiment_builder:

Data Formatting
===============

Classes and functions for converting data to SpaceTx Format and building experiments.

Converting Structured Data
--------------------------

:ref:`format_structured_data` tutorial.

.. code-block:: python

    from starfish.experiment.builder import format_structured_dataset

.. automodule:: starfish.core.experiment.builder.structured_formatter
   :members: format_structured_dataset

Tile Fetcher Interface
----------------------

:ref:`format_tilefetcher` tutorial and :ref:`tilefetcher_loader` tutorial.

.. code-block:: python

    from starfish.experiment.builder import FetchedTile, TileFetcher, write_experiment_json

.. automodule:: starfish.core.experiment.builder.providers
   :members: FetchedTile, TileFetcher

.. automodule:: starfish.core.experiment.builder.builder
   :members: write_experiment_json
.. _slicedimage:

slicedimage
===========

.. automodule:: slicedimage._formats
   :members:
.. _validation:

Validation
==========

Validators are provided for validating a SpaceTx fileset against the :ref:`schema`.

Validators
----------

.. autoclass:: starfish.core.spacetx_format.util.SpaceTxValidator
   :members:
   :exclude-members: fuzz_object

Helpers
-------

In addition, the starfish.spacetx_format.validate_sptx module contains helpers to simplify
iterating over the tree of json files and their respective schemas.


.. automodule:: starfish.core.spacetx_format.validate_sptx
   :members:


Error messages
--------------

Descriptive error messages are printed as warnings while validation takes place.
For example:

::

    starfish/starfish/core/spacetx_format/util.py:82: UserWarning:
     'contents' is a required property
            Schema:                 https://github.com/spacetx/starfish/starfish/file-format/schema/fov-manifest.json
            Subschema level:        0
            Path to error:          required
            Filename:               ../field_of_view/field_of_view.json
    
      warnings.warn(message)


This message tells you which schema has failed validation (``fov-manifest.json``), what type of error
has been encountered (``a required field is missing``), and the name of the file which is invalid
(``field_of_view.json``) if it has been provided. Validation of a json object will simply omit the
``Filename:`` field.
.. _config:

Configuration
=============

.. _starfishconfig:

StarfishConfig
--------------

.. autoclass:: starfish.config.StarfishConfig
       :members:

.. _environ:

environ
-------

.. autoclass:: starfish.config.environ
       :members:
.. _logging:

Provenance Logging
===================

Every method run on an ImageStack is recorded on the stacks log attribute. Each entry includes

- Method name
- Arguments supplied to the method
- os information
- starfish version
- release tag


Example
--------
This is an example of the formatted provenance log after processing an ImageStack through the ISS pipeline.

.. code-block:: python

    >>>pprint(stack.log)

    [{'arguments': {'clip_method': '"<Clip.CLIP: \'clip\'>"',
                    'is_volume': False,
                    'masking_radius': 15},
      'dependencies': {'numpy': '1.16.1',
                       'pandas': '0.24.1',
                       'scikit-image': '0.14.2',
                       'scikit-learn': '0.20.2',
                       'scipy': '1.2.1',
                       'sympy': '1.3',
                       'xarray': '0.11.3'},
      'method': 'WhiteTophat',
      'os': {'Platform': 'Darwin',
             'Python Version': '3.6.5',
             'Version:': 'Darwin Kernel Version 17.7.0: Thu Jun 21 22:53:14 PDT '
                         '2018; root:xnu-4570.71.2~1/RELEASE_X86_64'},
      'release tag': 'Running starfish from source',
      'starfish version': '0.0.36+4.g169fe89b.dirty'},
     {'arguments': {},
      'dependencies': {'numpy': '1.16.1',
                       'pandas': '0.24.1',
                       'scikit-image': '0.14.2',
                       'scikit-learn': '0.20.2',
                       'scipy': '1.2.1',
                       'sympy': '1.3',
                       'xarray': '0.11.3'},
      'method': 'Warp',
      'os': {'Platform': 'Darwin',
             'Python Version': '3.6.5',
             'Version:': 'Darwin Kernel Version 17.7.0: Thu Jun 21 22:53:14 PDT '
                         '2018; root:xnu-4570.71.2~1/RELEASE_X86_64'},
      'release tag': 'Running starfish from source',
      'starfish version': '0.0.36+4.g169fe89b.dirty'},
     {'arguments': {},
      'dependencies': {'numpy': '1.16.1',
                       'pandas': '0.24.1',
                       'scikit-image': '0.14.2',
                       'scikit-learn': '0.20.2',
                       'scipy': '1.2.1',
                       'sympy': '1.3',
                       'xarray': '0.11.3'},
      'method': 'Warp',
      'os': {'Platform': 'Darwin',
             'Python Version': '3.6.5',
             'Version:': 'Darwin Kernel Version 17.7.0: Thu Jun 21 22:53:14 PDT '
                         '2018; root:xnu-4570.71.2~1/RELEASE_X86_64'},
      'release tag': 'Running starfish from source',
      'starfish version': '0.0.36+4.g169fe89b.dirty'},
     {'arguments': {'detector_method': '"<function blob_log at 0x1233208c8>"',
                    'is_volume': True,
                    'max_sigma': 10,
                    'measurement_function': '"<function mean at 0x10c41c378>"',
                    'min_sigma': 1,
                    'num_sigma': 30,
                    'overlap': 0.5,
                    'threshold': 0.01},
      'dependencies': {'numpy': '1.16.1',
                       'pandas': '0.24.1',
                       'scikit-image': '0.14.2',
                       'scikit-learn': '0.20.2',
                       'scipy': '1.2.1',
                       'sympy': '1.3',
                       'xarray': '0.11.3'},
      'method': 'BlobDetector',
      'os': {'Platform': 'Darwin',
             'Python Version': '3.6.5',
             'Version:': 'Darwin Kernel Version 17.7.0: Thu Jun 21 22:53:14 PDT '
                         '2018; root:xnu-4570.71.2~1/RELEASE_X86_64'},
      'release tag': 'Running starfish from source',
      'starfish version': '0.0.36+4.g169fe89b.dirty'}]

.. _utils:

Utilities
=========

A number of utilities exist for simplifying work with the starfish library.

The :ref:`StarfishConfig` object can be instantiated at any point to provide configuration
regarding caching, validation, and similar low-level concerns. This is especially important
when directly calling out to IO backends like slicedimage. If a configuration setting needs
to be temporarily modified, use the :ref:`environ` context manager to set individual values.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. toctree::
   config.rst

.. toctree::
    logging.rst.. _image:

Image Manipulation
==================


starfish provides a variety of image manipulation methods that aid in the quantification of image-based transcriptomics
experiments. These include :py:class:`~starfish.image.Filter`, which remove background fluorescence and enhance spots,
:py:class:`~starfish.image.LearnTransform`, which learn transforms to align images across rounds and channels,
:py:class:`~starfish.image.ApplyTransform`, which apply learned transforms to images, and finally,
:py:class:`~starfish.image.Segmentation`, to identify the locations of cells.


.. _filtering:

Filtering
---------

Filters can be imported using ``starfish.image.Filter``, which registers all classes that subclass
``FilterAlgorithm``:

.. code-block:: python

    from starfish.image import Filter

.. automodule:: starfish.image.Filter
   :members:


.. _learn_transform:

Learn Transform
---------------

LearnTransform can be imported using ``starfish.image.LearnTransform``, the subclasses of
``LearnTransformAlgorithm`` are available for transform learning.

.. code-block:: python

    from starfish.image import LearnTransform

.. automodule:: starfish.image.LearnTransform
   :members:


.. _apply_transform:

Apply Transform
---------------

ApplyTransform can be imported using ``starfish.image.ApplyTransform``, the subclasses of
``ApplyTransformAlgorithm`` are available for transform learning.

.. code-block:: python

    from starfish.image import ApplyTransform

.. automodule:: starfish.image.ApplyTransform
   :members:


.. _segmentation:

Segmentation
------------

Segmentation can be imported using ``starfish.image.Segment``, which registers all classes that subclass
``SegmentAlgorithm``:

.. code-block:: python

    from starfish.image import Segment

.. automodule:: starfish.image.Segment
   :members:
.. _plot_utils:

Plotting Utilities
==================

This module contains a series of utilities for creating two dimensional plots that are useful for
generating documentation and vignettes. We suggest that users leverage :py:func:`starfish.display`
for their plotting needs, as the interactive viewer is better able to handle the array of features
that starfish needs.

.. code-block:: python

    from starfish.util import plot

.. automodule:: starfish.util.plot
    :members:


.. _spots:

Spots
=====

Starfish provides a number of methods for which spots (or other regions of interest) are the main substrate.
These include :py:class:`starfish.spots.DetectPixels`, which exposes methods that identify which target code best corresponds to each pixel, and merges adjacent pixels into ROIs,
:py:class:`starfish.spots.FindSpots`, which exposes methods that find bright spots against dark backgrounds,
:py:class:`starfish.spots.DecodeSpots`, which exposes methods that match patterns of spots detected across rounds and channels in the same spatial positions with target codes, and
:py:class:`starfish.spots.AssignTargets`, which exposes methods to assign spots to cells.

.. _detect_pixels:

Detecting Pixels
----------------

Pixel Detectors can be imported using ``starfish.spots.DetectPixels``, which registers all classes that subclass ``DetectPixelsAlgorithm``:

.. code-block:: python

    from starfish.spots import DetectPixels

.. automodule:: starfish.spots.DetectPixels
    :members:

.. _spot_finding:

Finding Spots
---------------

Spot Finders can be imported using ``starfish.spots.FindSpots``, which registers all classes that subclass ``FindSpotsAlgorithm``:


.. _`Spot Finding Refactor Plan`: https://github.com/spacetx/starfish/issues/1514

.. code-block:: python

    from starfish.spots import FindSpots

.. automodule:: starfish.spots.FindSpots
    :members:

.. _decode_spots:

Decoding Spots
---------------

Spot Decoders can be imported using ``starfish.spots.DecodeSpots``, which registers all classes that subclass ``DecodeSpotsAlgorithm``:

.. code-block:: python

    from starfish.spots import DecodeSpots

.. automodule:: starfish.spots.DecodeSpots
   :members:


.. _target_assignment:

Target Assignment
-----------------

Target Assignment can be imported using ``starfish.spots.AssignTargets``, which registers all classes that subclass ``AssignTargetsAlgorithm``:

.. code-block:: python

    from starfish.spots import AssignTargets

.. automodule:: starfish.spots.AssignTargets
   :members:
.. _datasets:

Datasets
========

We maintain an API to let users access a handful of curated datasets for learning how to use
starfish, and benchmarking algorithms and pipelines. See how to load these example datasets
:ref:`here<loading_data>`.


.. automodule:: starfish.data
   :members:
.. _ImageStack:

ImageStack
==========

.. autoclass:: starfish.core.imagestack.imagestack.ImageStack
   :members:
.. _ExpressionMatrix:

ExpressionMatrix
================

Expression Matrix is a wrapper for xarray.DataArray that provides serialization adapters to netCDF,
loom, and AnnData to enable it to be used with single-cell analysis software packages such as
scanpy (python) and seurat (R).

.. autoclass:: starfish.core.expression_matrix.expression_matrix.ExpressionMatrix
   :members:
.. _FieldOfView:

Field of View
=============

.. autoclass:: starfish.core.experiment.experiment.FieldOfView
    :members:
.. _binary_mask:

BinaryMaskCollection
====================

.. autoclass:: starfish.morphology.BinaryMaskCollection
   :members:
.. _Codebook:

Codebook
========

.. autoclass:: starfish.core.codebook.codebook.Codebook
   :members:
.. _IntensityTable:

IntensityTable
==============

.. autoclass:: starfish.core.intensity_table.intensity_table.IntensityTable
   :members:

.. automodule:: starfish.core.intensity_table.intensity_table_coordinates
   :members:
.. _Experiment:

Experiment
==========

.. autoclass:: starfish.core.experiment.experiment.Experiment
    :members:
.. _data structures:

Data Structures
===============

The top-level object in a starfish workflow is the :ref:`Experiment`. It is composed of one or more
:ref:`FieldOfView` objects, and a :ref:`Codebook`, which maps detected spots to the entities they
target.

Each :ref:`FieldOfView` consists of a set of Primary Images and optionally, Auxiliary images that
may contain information on nuclei (often used to seed segmentation) or fiduciary beads (often used
to enable fine registration).

Both Primary and Auxiliary Images are referenced by slicedimage_ TileSet_ objects, which map
two dimensional image tiles stored on disk into a 5-dimensional Image Tensor that labels each
``(z, y, x)`` tile with the ``round`` and ``channel`` that it corresponds to. When loaded into
memory, these Image Tensors are stored in :ref:`ImageStack` objects. The :ref:`ImageStack` is what
starfish uses to execute image pre-processing, and serves as the substrate for spot finding.

Identified spots are stored in the :ref:`IntensityTable`, which stores the intensity of the spot
across each of the rounds and channels that it is detected in. It also stores assigned genes when
decoded with a :ref:`Codebook` and assigned cells when combined with a segmentation results.

Finally, the :ref:`IntensityTable` can be converted into an :ref:`ExpressionMatrix` by summing all
of the spots detected for each gene across each cell. The ExpressionMatrix provides conversion and
serialization for use in single-cell analysis environments such as Seurat_ and Scanpy_.

.. TODO ambrosejcarr: think about removing PipelineComponent to another part of the docs

.. _slicedimage: https://github.com/spacetx/slicedimage

.. _TileSet: https://github.com/spacetx/slicedimage/blob/master/slicedimage/_tileset.py

.. _Seurat: https://satijalab.org/seurat/

.. _Scanpy: https://scanpy.readthedocs.io/en/latest/

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. toctree::
   experiment.rst

.. toctree::
   field_of_view.rst

.. toctree::
   image_stack.rst

.. toctree::
   codebook.rst

.. toctree::
   expression_matrix.rst

.. toctree::
   intensity_table.rst

.. toctree::
   decoded_intensity_table.rst

.. toctree::
   binary_mask.rst

.. toctree::
   label_image.rst
.. _DecodedIntensityTable:

DecodedIntensityTable
=====================

.. autoclass:: starfish.core.intensity_table.decoded_intensity_table.DecodedIntensityTable
   :members:
.. _label_image:

LabelImage
==========

.. autoclass:: starfish.morphology.LabelImage
   :members:
.. _types:

Types
=====

.. contents::
    :local:

Starfish uses a series of constants and container types to enable type checking, and provides a set
of useful hints about how to interact with starfish objects.

Coordinates
-----------

Coordinates holds constants that store with the physical coordinates of a field of view. They define
a field of view's relative location to some global scale parameter, and identify how to stitch or
combine multiple fields of view.

.. autoclass:: starfish.core.types.Coordinates
    :members:
    :undoc-members:

Physical Coordinates
---------------------
.. autoclass:: starfish.core.types.PhysicalCoordinateTypes
    :members:
    :undoc-members:

Axes
----
Axes holds constants that represent indexers into the dimensions of the :py:class:`ImageStack`
5-d image tensor. They are re-used by objects that inherit subsets of these Axes, such as:

1. :py:class:`IntensityTable`, which stores spot coordinates and pixel traces across
   rounds and channels
2. :py:class:`Codebook`, which stores expected image intensities across imaging rounds and
   channels

.. autoclass:: starfish.core.types.Axes
    :members:
    :undoc-members:

Features
--------

Features holds constants that represent characteristics of detected image features (most often
spots, but sometimes also individual pixels).

.. autoclass:: starfish.core.types.Features
    :members:
    :undoc-members:


SpotAttributes
--------------

SpotAttributes defines the minimum amount of information required by starfish to describe
a spot. It also contains methods to save these attributes to files that can be used to visualize
detected spots.

.. autoclass:: starfish.core.types.SpotAttributes
    :members:

Levels
------

.. autoclass:: starfish.types.Levels
    :members:
    :undoc-members:
.. _display:

Interactive Image Viewer: display()
===================================

For fast and interactive visualization of an :py:class:`.ImageStack`,
:py:class:`.IntensityTable`, and :py:class:`.BinaryMaskCollection`, starfish utilizes napari.
Detailed usage instructions can be found in the :ref:`plot_display` tutorial.

.. code-block:: python

    from starfish import display

.. automodule:: starfish.core._display
    :members:.. _morphology:

Morphology Transformations
==========================

starfish provides a variety of methods to perform transformations on morphological data.  These include:

* :py:class:`~starfish.morphology.Binarize`, which transforms image data into morphological data.
* :py:class:`~starfish.morphology.Filter`, which performs filtering operations on morphological data.
* :py:class:`~starfish.morphology.Merge`, which combines different sets of morphological data.
* :py:class:`~starfish.morphology.Segment`, which performs segmentation operations to yield morphological data.

.. _binarize:

Binarize
--------

Binarizing operations can be imported using ``starfish.morphology.Binarize``, which registers all classes that subclass :py:class:`~starfish.morphology.Binarize.BinarizeAlgorithm`:

.. code-block:: python

    from starfish.morphology import Binarize

.. automodule:: starfish.morphology.Binarize
   :members:

.. _morphological_filter:

Filter
------

Filtering operations can be imported using ``starfish.morphology.Filter``, which registers all classes that subclass :py:class:`~starfish.morphology.Filter.FilterAlgorithm`:

.. code-block:: python

    from starfish.morphology import Filter

.. automodule:: starfish.morphology.Filter
   :members:

.. _merge:

Merge
-----

Filtering operations can be imported using ``starfish.morphology.Merge``, which registers all classes that subclass :py:class:`~starfish.morphology.Merge.MergeAlgorithm`:

.. code-block:: python

    from starfish.morphology import Merge

.. automodule:: starfish.morphology.Merge
   :members:

Segment
-------

Filtering operations can be imported using ``starfish.morphology.Segment``, which registers all classes that subclass :py:class:`~starfish.morphology.Segment.SegmentAlgorithm`:

.. code-block:: python

    from starfish.morphology import Segment

.. automodule:: starfish.morphology.Segment
   :members:
