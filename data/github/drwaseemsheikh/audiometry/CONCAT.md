# audiometry

Audiometry is an open-source application framework based on WPF and .NET to create hearing test related applications. The framework is built using the Model-View-ViewModel (MVVM) software architectural pattern.

Audiometry can store, process, and visualize data corresponding to tuning fork tests including Weber, Rinne, Schwabach, absolute bone conduction, Teal, and Gelle; speech audiometry; pure-tone audiometry (PTA) including air conduction masked, air conduction unmasked, bone conduction masked, bone conduction unmasked, air conduction aided, loudness level, and sound field; impedance audiometry; bithermal caloric test; and advanced tests including alternate binaural loudness balance (ABLB), short increment sensitivity index (SISI), tone decay, and Stenger.

The application framework is extensible and can be used to develop new hearing test applications by extending the current functionality.

Audiometry is independent of specific hearing test hardware thus making it possible to be used with a wide variety of hearing test hardware.

Audiometry provides a unified and uniform interface for storing, processing, and visualizing data from a wide range of hearing tests that traditionally rely on different hardware and software to process and store data.

# Installation

``Audiometry`` can be installed on a Windows 7 or Windows 10 machine. To install the application, run the AudiometryInstaller.msi in the installer folder of the repository.

# Documentation

A description of the software can be found at https://github.com/drwaseemsheikh/audiometry/blob/master/docs/full_paper.pdf.
API docs can be found at https://drwaseemsheikh.github.io/audiometry/.
The JOSS paper is available at [![DOI](https://joss.theoj.org/papers/10.21105/joss.02016/status.svg)](https://doi.org/10.21105/joss.02016).
---
title: 'Audiometry: A model-view-viewmodel (MVVM) application framework for hearing impairment diagnosis'
tags:
  - C#
  - WPF
  - audiometry
  - hearing loss
  - model-view-viewmodel
authors:
  - name: Waseem Sheikh
    orcid: 0000-0002-4647-4565
    affiliation: 1
  - name: Nadeem Sheikh
    affiliation: 2
affiliations:
 - name: Associate Professor, Electrical and Computer Engineering, Utah Valley University, USA
   index: 1
 - name: Assistant Professor of ENT, CMH, Quetta, Pakistan
   index: 2
date: 9 Nov 2019
bibliography: paper.bib
---

# Summary

Around 466 million people worldwide (over 5% of the world's population) have disabling hearing loss, and out of these 34 million are children [@Who:2019]. Estimates suggest that by 2050, over 900 million people worldwide will have disabling hearing loss. The annual global cost of unaddressed hearing loss amounts to US$ 750 billion [@Who:2019]. Early detection of hearing loss can reduce its impact on an individual's life in addition to saving a huge cost. The existing hearing test applications are closed-source, not extensible, test for a limited number of hearing tests such as pure-tone air conduction audiometry, the audiograms generated are either incomplete or do not fully conform to the American National Standards Institute (ANSI) ANSI S3.6-1996 Specification for Audiometers [@ANSIS3.6-1996:1996], are tightly coupled with a specific vendor hardware, and do not provide an ability to implement various data analytics algorithms to draw important conclusions from the hearing test data [@ChenSmartphone:2018; @SamelliTablet:2018; @BarczikAccuracy:2018; @LivshitzApplication:2017; @AbuSmartphone:2016; @YaoBrowser:2015]. In addition, the price of proprietary hearing test software applications makes these prohibitive for underdeveloped countries which tend to have a higher prevalence of people with hearing loss. In most underdeveloped countries, hearing test data is still stored on paper and graphs such as audiograms are drawn by hand. Such a primitive system of managing hearing test data is error-prone and makes it very difficult to save, track, analyze, and reproduce hearing test data. In addition, a lack of open-source software in this domain stifles innovation.

``Audiometry`` is an open-source application framework written in C# and based on WPF and .NET to create hearing test applications. ``Audiometry`` enables accurate digital recording, search, analysis, graphical visualization, and reproduction of human audio-vestibular impairment test data to assist in hearing loss or disability diagnosis. The framework is built using the Model-View-ViewModel (MVVM) [@MvvmMicrosoft:2019; @MvvmWikipedia:2019] software architectural pattern which separates the development of graphical user interface (GUI) from the development of business and back-end logic. Some of the benefits of the MVVM pattern include reusable components, independent development of GUI and business or back-end logic, flexibility to modify GUI without having to change business or back-end logic, ease of comprehensive unit testing, faster application development time, and reduced maintenance overhead. The proposed framework makes it possible to easily extend the application functionality thus enabling other researchers and practitioners to develop their own hearing impairment diagnosis applications.

``Audiometry`` can store, search, analyze, print, and visualize data corresponding to tuning fork tests including Weber, Rinne, Schwabach, absolute bone conduction, Teal, and Gelle; speech audiometry; pure-tone audiometry (PTA); impedance audiometry; bithermal caloric test; and advanced tests including alternate binaural loudness balance (ABLB), short increment sensitivity index (SISI), tone decay, and Stenger [@KramerAudiology:2019; @GelfandEssentials:2016; @KatzHandbook:2015; @DhingraDiseases:2014; @BessAudiology:2003]. The application framework can also be used to develop new hearing test applications by extending its current functionality. ``Audiometry`` is independent of specific hearing test hardware thus making it possible to be used with a wide variety of hearing test hardware. In addition, ``Audiometry`` provides a unified and uniform interface for storing, analyzing, and visualizing data from a wide range of hearing tests which traditionally rely on different hardware and software. The software was evaluated by an otolaryngologist who found it to be very beneficial in reaching a hearing impairment diagnosis conclusion more methodically, swiftly, and accurately.

Following are examples of some of the research questions that can be investigated by the use of ``Audiometry``:

1.	The software can be used to compare the sensitivity (true-positive rate) and specificity (true-negative rate) of various hearing test equipment and methods. For example, questions like how reliable is a pure-tone audiometry test performed by a smartphone or a tablet when compared to a benchmark calibrated audiometer can be easily answered by using this software.

2.	The software can be used to determine important correlations between lifestyle, work conditions, and demographics; and the types of hearing loss.

3.	The software can be used to measure the efficacy of a certain treatment, intervention, or equipment on the progression of a hearing loss.

The current functionality of the application can be extended and enhanced in various ways. Some important future research directions include adding additional hearing impairment diagnostic intelligence into the application, using machine learning and artificial intelligence techniques to increase the accuracy of diagnosis, and a client-server based architecture of the application.

# Figures

![Pure-tone audiogram interface.](puretone1.png)

![Speech audiometry interface.](speech1.png)

![Impedance audiometry interface.](impedance1.png)

![Bithermal caloric interface.](calorigram1.png)

![Search interface.](search1.png)

![Patient interface.](patient1.png)

# Documentation
The Doxygen generated API documentation for ``Audiometry`` can be found under the docs folder. The full-length paper on ``Audiometry`` which explains its design, architecture, and implementation in detail is located in the paper folder.

# Installation

``Audiometry`` can be installed on a Windows 7 or Windows 10 machine. To install the application, run the AudiometryInstaller.msi in the installer folder of the repository. To test the application, please follow the steps listed in the test.md file under the test folder.

# References
---
name: Bug Report
about: Report a bug here
labels: bug
---

Please replace every line in curly brackets { like this } with appropriate answers, and remove this line.

## Description

{ Please provide a clear and concise description of the bug. }

## Environment

1. OS (where Audiometry runs): { e.g. Mac OS 10, Windows 10, Ubuntu 16.4, etc. }
2. Audiometry version: { e.g. Audiometry 1.0 }
3. Other environment details:

## Reproducible Steps

Steps to create the smallest reproducible scenario:
1. { e.g. Run ... }
2. { e.g. Click ... }
3. { e.g. Save ... }
4. { e.g. See error ... }

## Expected Output

{ Please describe what you expected to happen. }

## Actual Output

{ Please describe what actually happened. }
 
## Additional information

{ Any additional information, including logs or screenshots if you have any. }
---
name: Feature Request
about: Request a feature here
labels: feature
---

Please replace every line in curly brackets { like this } with appropriate answers, and remove this line.

## Problem to Solve

{ Please describe the problem you would like to solve. }

## Current Workaround

{ Please describe how you currently solve or work around this problem, given Audiometry's limitation. }

## Proposed Solution

{ Please describe the solution you would like Audiometry to provide, to solve the problem above. }

## Additional Information

{ Any additional information, including logs or screenshots if you have any. }
To test the ``Audiometry`` application, please follow the steps given below:
1. Install the application using the AudiometryInstaller.msi file under the installer folder. The application can be installed on a Windows 7 or Windows 10 machine.
2. After installation, run the application. It will prompt you to register as a user. Provide a username and password.
3. You can create a patient's profile and enter hearing test data for various hearing tests supported by the application such as pure-tone audiometry, speech audiometry, impedance audiometry, and bithermal caloric test.
4. As the hearing test data is entered, the various visualization plots are automatically generated including pure-tone audiogram, speech audiogram, tympanogram, and calorigram.
5. The application will also automatically generate various metrics for hearing disability diagnosis based on the data entered.
6. The patient data can be saved in a database by clicking on the Save menu option.
7. The patient data can be searched and opened by using the Search menu option.
8. The patient data can be modified and saved using the Save menu option.
9. The patient data can be deleted using the Delete menu option.
10. A hearing test report can be printed or a PDF file exported using the Print menu option.
OxyPlot is a cross-platform plotting library for .NET

- [Web page](http://oxyplot.org)  
- [Documentation](http://oxyplot.org/documentation)
- [Announcements](http://oxyplot.org/announcements) / [atom](http://oxyplot.org/atom.xml)
- [Discussion forum](http://discussion.oxyplot.org)
- [Source repository](http://github.com/oxyplot/oxyplot)
- [Issue tracker](http://github.com/oxyplot/oxyplot/issues)
- [NuGet packages](http://www.nuget.org/packages?q=oxyplot)
- [Stack Overflow](http://stackoverflow.com/questions/tagged/oxyplot)
- [Twitter](https://twitter.com/hashtag/oxyplot)
- [Gitter](https://gitter.im/oxyplot/oxyplot)

[![Build status](https://ci.appveyor.com/api/projects/status/mlaqnruo6ic3oe60)](https://ci.appveyor.com/project/objorke/oxyplot)

[![Gitter](https://badges.gitter.im/Join Chat.svg)](https://gitter.im/oxyplot/oxyplot?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

![Plot](http://oxyplot.org/public/images/normal-distributions.png)

#### Branches

`master` - the release branch (stable channel)  
`develop` -  the main branch with the latest development changes (pre-release channel)

See '[A successful git branching model](http://nvie.com/posts/a-successful-git-branching-model/)' for more information about the branching model in use.

#### Getting started

1. Use the NuGet package manager to add a reference to OxyPlot
2. Add a `PlotView` to your user interface
3. Create a `PlotModel` in your code
4. Bind the `PlotModel` to the `Model` property of your `PlotView`

#### Examples

You can find examples in the `/Source/Examples` folder in the code repository.

#### Contribute

See [the documentation](http://oxyplot.org/documentation/contributions) for information about how to contribute!
OxyPlot is a cross-platform plotting library for .NET

- [Web page](http://oxyplot.org)  
- [Documentation](http://oxyplot.org/documentation)
- [Announcements](http://oxyplot.org/announcements) / [atom](http://oxyplot.org/atom.xml)
- [Discussion forum](http://discussion.oxyplot.org)
- [Source repository](http://github.com/oxyplot/oxyplot)
- [Issue tracker](http://github.com/oxyplot/oxyplot/issues)
- [NuGet packages](http://www.nuget.org/packages?q=oxyplot)
- [Stack Overflow](http://stackoverflow.com/questions/tagged/oxyplot)
- [Twitter](https://twitter.com/hashtag/oxyplot)
- [Gitter](https://gitter.im/oxyplot/oxyplot)

[![Build status](https://ci.appveyor.com/api/projects/status/mlaqnruo6ic3oe60)](https://ci.appveyor.com/project/objorke/oxyplot)

[![Gitter](https://badges.gitter.im/Join Chat.svg)](https://gitter.im/oxyplot/oxyplot?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

![Plot](http://oxyplot.org/public/images/normal-distributions.png)

#### Branches

`master` - the release branch (stable channel)  
`develop` -  the main branch with the latest development changes (pre-release channel)

See '[A successful git branching model](http://nvie.com/posts/a-successful-git-branching-model/)' for more information about the branching model in use.

#### Getting started

1. Use the NuGet package manager to add a reference to OxyPlot
2. Add a `PlotView` to your user interface
3. Create a `PlotModel` in your code
4. Bind the `PlotModel` to the `Model` property of your `PlotView`

#### Examples

You can find examples in the `/Source/Examples` folder in the code repository.

#### Contribute

See [the documentation](http://oxyplot.org/documentation/contributions) for information about how to contribute!
