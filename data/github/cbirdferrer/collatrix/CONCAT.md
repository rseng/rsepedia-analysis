# CollatriX
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02328/status.svg)](https://doi.org/10.21105/joss.02328)  [![DOI](https://zenodo.org/badge/243385218.svg)](https://zenodo.org/badge/latestdoi/243385218) [![Anaconda-Server Badge](https://anaconda.org/cbird/collatrix/badges/version.svg)](https://anaconda.org/cbird/collatrix)

## Background
This function collates the csv outputs from the MorphoMetriX photogrammetry GUI (https://github.com/wingtorres/morphometrix) into one large single data frame containing the image, animal ID, measurements, and notes.  
CollatriX was designed with several add-ons. A figure showing the different routes available is included below:  
![alt text](https://github.com/cbirdferrer/collatrix/blob/master/images/Figure1.png)

The altitude calibration function (`collatrix.altitude_calib`) can be used to calculate corrected altitudes using images of an object of known length. If used, this function should be used before the main function. The output can be used to create the safety input file for the main `collatrix` function. Note that the altitude calibration function is not required, the user can start the workflow using the main `collatrix` function. The output of this main function can then be used to calculate metrics of whale body condition (`collatrix.whale_bc`) if desired.

## Documentation
Information on how to install and use `collatrix` and example code can be found in our [wiki](https://github.com/cbirdferrer/collatrix/wiki)

## Demo
A demonstration is available in the [demo](https://github.com/cbirdferrer/collatrix/tree/master/demo) directory. The demo includes a separate README file with instructions for what inputs to use.

### Automated Testing
A GitHub Action has been set up to test the main `collatrix` function when an update is pushed to the master branch using `pytest`. If you are working on editing `collatrix` and would like to run the automated tests locally, open a terminal or command prompt window, then change the directory (`cd`) to the folder where you have cloned `collatrix` to, then type `pytest`, the test should then run.

For example:

```bash
(base) :~ user$ cd github/collatrix
(base) :collatrix user$ pytest
```

# Attribution
If you use this software please cite our paper:  
*Bird C.N. & Bierlich K.C., (2020). CollatriX: A GUI to collate MorphoMetriX outputs. Journal of Open Source Software, 5(51), 2328, https://doi.org/10.21105/joss.02328*

Additionally:
* if you used `collatrix.whale_bc` to calculate Body Volume cite:
*Christiansen, F., Vivier, F., Charlton, C., Ward, R., Amerson, A., Burnell, S., & Bejder, L. Maternal body size and condition determine calf growth rates in southern right whales (2018). Maternal body size and condition determine calf growth rates in southern right whales. Marine Ecology Progress Series, 592, 267–281. http://doi.org/10.3354/meps12522*

* if you used `collatrix.whale_bc ` to calculate BAI or `collatrix.altitude_calib` for the altitude calibration cite: 
*Burnett, Jonathan D., Leila Lemos, Dawn Barlow, Michael G. Wing, Todd Chandler, and Leigh G. Torres. 2018. “Estimating Morphometric Attributes of Baleen Whales with Photogrammetry from Small UASs: A Case Study with Blue and Gray Whales.” Marine Mammal Science 35 (1): 108–39. https://doi.org/10.1111/mms.12527.*

# Contributing
We designed CollatriX with future collaborations in mind and we'd love for you to contribute! If you'd like to contribute please see our [contributing guidelines](https://github.com/cbirdferrer/collatrix/blob/master/CONTRIBUTING.md)

# Code of Conduct
See our [code of conduct](https://github.com/cbirdferrer/collatrix/blob/master/CODE_OF_CONDUCT.md)

# License
[![Anaconda-Server Badge](https://anaconda.org/cbird/collatrix/badges/license.svg)](https://anaconda.org/cbird/collatrix)

Copyright (c) 2020 Clara Bird, KC Bierlich

`Collatrix` is free software made available under the MIT License. For details see the the [LICENSE](https://github.com/cbirdferrer/collatrix/blob/master/LICENSE) file.

# Contributors
Clara N. Bird and KC Bierlich are the developers of this software.
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
 professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at clara.birdferrer@gmail.com. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# Contribution Guidelines

## Reporting issues

- **Search for existing issues.** Please check to see if someone else has reported the same issue.
- **Share as much information as possible.** Include operating system and version, browser and version. Also, include steps to reproduce the bug.

## Project Setup
Refer to the [README](README.md).

## Code Style

This code generally tries to adhere to [PEP-8]( <https://www.python.org/dev/peps/pep-0008/>) standards for style, howevever we emphasize the PEP-8 team's point that "A Foolish Consistency is the Hobgoblin of Little Minds". From their website...

*Some other good reasons to ignore a particular guideline:*

* When applying the guideline would make the code less readable, even for someone who is used to reading code that follows this PEP.
* To be consistent with surrounding code that also breaks it (maybe for historic reasons) -- although this is also an opportunity to clean up someone else's mess (in true XP style).
* Because the code in question predates the introduction of the guideline and there is no other reason to be modifying that code.
* When the code needs to remain compatible with older versions of Python that don't support the feature recommended by the style guide.

## Testing
Use the material provided in the [demo]( <https://github.com/cbirdferrer/collatrix/tree/master/demo>)  directory for testing

## Pull requests
- Try not to pollute your pull request with unintended changes – keep them simple and small. If possible, squash your commits.
- Try to share how your code has been tested before submitting a pull request.
- If your PR resolves an issue, include **closes #ISSUE_NUMBER** in your commit message (or a [synonym](https://help.github.com/articles/closing-issues-via-commit-messages)).
- Review
    - If your PR is ready for review, another contributor will be assigned to review your PR
    - The reviewer will accept or comment on the PR. 
    - If needed address the comments left by the reviewer. Once you're ready to continue the review, ping the reviewer in a comment.
    - Once accepted your code will be merged to `master`
    
 ## Attribution
 These guidelines are adapted from the guidelines for contributing to [MorphoMetriX]( <https://github.com/wingtorres/morphometrix>).
---
title: 'CollatriX: A GUI to collate MorphoMetriX outputs'
tags:
  - Python
  - Drones
  - Photogrammetry
  - GUI
  - Morphometry

authors:
  - name: Clara N Bird
    orcid: 0000-0001-7763-7761
    affiliation: "1, 2"
  - name: KC Bierlich
    affiliation: 2
affiliations:
 - name: Geospatial Ecology of Marine Megafauna Lab, Marine Mammal Institute, Department of Fisheries and Wildlife, Oregon State University
   index: 1
 - name: Nicholas School of the Environment, Duke University Marine Laboratory
   index: 2
date: 22 May 2020
bibliography: paper.bib

---

# Summary

CollatriX is a graphical user interface (GUI) developed using PyQt5 to collate outputs from MorphoMetriX [@Torres:2020], a photogrammetric measurement GUI designed for morphometric analysis of wild animals from aerial imagery. For each image used in MorphoMetriX, a comma-separated-values sheet (.csv) is produced containing the custom measurements (length, area, angle) and their associated labels created by the user. Hence, projects with a large number of images in their morphometric analysis have a large number of output files, creating time-intensive and tedious workflows to manually combine output files into a single data file to be used in analysis. CollatriX was designed as a user-friendly GUI overcoming this limitation by collating these measurement outputs into a single data sheet (.csv) based on the animal’s individual ID (Fig. 1). Furthermore, CollatriX has two add-on functions, one to correct for altitude error from Unoccupied Aerial Systems (UAS or drone) flights and another for calculating different animal body condition metrics, following body volume from @Christiansen:2018 and body area index (BAI) from @Burnett:2018 (Fig. 1). The framework of CollatriX was also designed to have the flexibility to accommodate and encourage other future add-on functions.

# Main Features

CollatriX will collate MorphoMetriX output files saved in various file structure formats (i.e., a single folder vs. across multiple folders). We also included an option to subset the collated data into a separate datasheet based on a specified list of Animal IDs provided by the user. A safety option was built in CollatriX to increase user efficiency in working through large image datasets while avoiding user input errors. For example, MorphoMetriX automatically scales length measurements in pixels to real world values (i.e., meters) from manually entered altitude, focal length, and pixel dimension values [@Torres:2020]. While this setup allows for each separate image to be scaled accordingly, there is potential for input errors when entering these values, especially when working through large datasets. CollatriX provides a safety option for users where the number of pixels in a length measurement is back calculated, and the measurement is recalculated using the correct values per image from a user provided csv.

# Add-on Functions

CollatriX also has two add-on functions for 1) calibrating altitude errors from a UAS flight and 2) calculating whale body condition (Fig. 1). The altitude calibration function follows the recommended Method 5 from @Burnett:2018, where measurements of a calibration object of known length are used to calculate the true altitude of the UAS to create a linear model used for correcting the altitude of images taken throughout the flight (for more detail see @Burnett:2018). The output of this add-on function can then be used as the safety for the main CollatriX function (Fig. 1).

The whale body condition add-on function calculates two commonly used metrics of cetacean body condition. If the user used MorphoMetriX to measure perpendicular widths based on a length measurement, a common method for assessing body condition in cetaceans [@Christiansen:2016; @Dawson:2017; @Burnett:2018], the function can calculate body volume of the whale following @Christiansen:2018 and Body Area Index (BAI) following @Burnett:2018. BAI is a measure of dorsal surface area normalized by length, and CollatriX can calculate the surface area using parabolas [@Burnett:2018] or trapezoids [@Christiansen:2016]. Since MorphoMetriX allows the user to specify the number of perpendicular width segments of a length measurement, CollatriX provides the flexibility for the user to specify the number of width segments to include in the body condition calculation, i.e., 20 widths (5% increments of the total length), as well as the minimum and maximum bound in which to calculate body volume or BAI (i.e., widths between 20-85% vs. 25-80% of total length).

Together, MorphoMetriX and CollatriX provide a toolkit that is flexible, easy to use, and adaptable to future projects on a variety of species and applications, as CollatriX is designed to easily incorporate other add-on functions. CollatriX has been used on MorphoMetriX outputs from several projects on a variety of cetacean species including bottlenose dolphins and Antarctic minke, dwarf minke, fin, blue, gray, and humpback whales.

# Figures

![Basic overview of CollatriX workflow using measurement outputs from MorphoMetriX (Torres and Bierlich 2020) Measurement outputs are collated into a single output file based on the ‘Image ID’. Solid arrows represent main pathway, dotted arrows represent pathway including the add-on functions (green boxes).](../images/Figure1.png)

# Acknowledgements

We acknowledge Dr. David Johnston and Dr. Leigh Torres for project support, Walter Torres for development support, and Julian Dale for beta testing. This work was supported by the Duke University Marine Laboratory and the Duke Marine Robotics and Remote Sensing Lab.


# References
# Demo

## File Table of Contents

### **Folders**

**altitude_calibration_outputs**: This folder contains example outputs from MorphoMetrix of measured calibration objects.

**measured_whale_outputs**: This folder contains example outputs from MorphoMetriX of measured whales. These measurements include perpendicular widths measured at 10% increments. Note that the manually entered altitude values are not necessarily correct and uncalibrated, so the values of the collated output will be different.

### **Altitude calibration files**
**calibration_obj_imglist.csv**: This is the list of the measured calibration object images. It includes the date, flight, altitude, focal length, and pixel dimension per image.

**image_list**: This is the list of measured whale image that includes the date, flight, focal length, pixel dimension, and altitude per image. The altitude in this datasheet is the altitude that was recorded by the drone that needs to be calibrated.

### **Main function files**
**demo_safety.csv**: This is the sheet that includes the calibrated altitudes, focal lenghts, and pixel dimensions. The output of the altitude calibration should look like this file.

**demo_animal_list.csv**: This is an example of the animal list that can be used to subset the collated datasheet by animal ID.

## Demo altitude calibration input instructions
- Use the calibration_obj_imglist.csv as the calibration object img list input
- Use the image_list.csv as the img list w/ altitudes input
- Board Length measurement name = BL
- True length of calibration object = 1
- Folder containing MorphoMetriX outputs = altitude_calibration_outputs

## Demo collatrix input instructions
- Animal ID from folder name?: either is fine
- Safety: yes, use altitude calibration file from altitude calibration function. Or provided demo_safety.csv
- List of Specific Individuals: either is fine, if yes use demo_animal_list.csv
- Output name: whatever you want
- Location of MorphoMetrix files: measured_whales_outputs folder

## Demo whale body condition input instructions
- Name of length measurement = TL
- Lower bound = 20
- Upper bound = 80
- Interval = 10
*use these settings for both body volume and BAI
