<p align="center"><img src="docs/Images/Logo_START_R.svg" alt="logo" height="150px"></p>

------

[![DOI](https://zenodo.org/badge/127922557.svg)](https://zenodo.org/badge/latestdoi/127922557)

Visit our website : https://thomasdenecker.github.io/START-R/

## Requirements


We use Docker to develop and manage START-R. We invite you to verify that the
following requirements are correctly satisfied before trying to bootstrap the
application:

* [Docker 1.12.6+](https://docs.docker.com/engine/installation/)

> We recommend you to follow Docker's official documentations to install
required docker tools (see links above).To help you, explanatory videos for each
operating system are available [here](https://www.bretfisher.com/installdocker/)

To help you with the installation, you will find below videos on youtube for each system:
- [Mac OS X](https://www.youtube.com/watch?v=mbSsh40_8WM)
- [Windows 10](https://www.youtube.com/watch?v=_9AWYlt86B8)
- [Linux](https://www.youtube.com/watch?v=8Iu5uqby9PY)

**Docker must be on for the duration of START-R use.**

**Important**

Note that the size of the RAM that should be allocated to the Docker depends on
the size of the studied organism genome. START-R can work with data from several
organims. For the human genome, we strongly recommend an increase in allocated memory for Docker.
Otherwise, the risk is an early termination of the analysis that will be incomplete.

A workstation or a laboratory server with 16GB of RAM is therefore well dimensioned.
To increase the allocated memory, go here for
- [Mac OS X](https://docs.docker.com/docker-for-mac/#memory)
- [Windows 10](https://docs.docker.com/docker-for-windows/#advanced)
- [Linux](https://docs.docker.com/config/containers/resource_constraints/#limit-a-containers-access-to-memory)

## Quick start

Have you read the "Requirements" section above?

### START-R project installation

Download the zip file ([here](https://github.com/thomasdenecker/START-R/archive/master.zip)), extract this file and copy the obtained folder where you want on your computer. Note that if you move the folder, the installation procedure will have to be redone.

**Reminder** : Docker must always be switched on for any installation and use of START-R !

#### Windows installation 

**IMPORTANT** : START-R needs Docker. It will only be possible to install on **Windows 10**.

In this folder, you will find a file named INSTALLATION_WINDOWS.bat. By double clicking on it, the installation will begin. This may take a little time depending on the quality of your internet connection. When the installation is completed, a new file will appear. They allow to launch the START-R applications.

#### Mac OsX installation

**In command line**

[Open a terminal](https://www.youtube.com/watch?v=QROX039ckO8) and run these commands:

```
git clone https://github.com/thomasdenecker/START-R.git
cd START-R
sudo ./INSTALLATION_MAC.sh
```

The installation will begin. This may take a little time depending on the quality of your internet connection. When the installation is completed, a new file will appear. They allow to launch the START-R applications. Once the installation is complete, use this command to launch START-R analyzer:
```
./START-R_analyzer.sh
```

and this command to launch START-R viewer
```
./START-R_viewer.sh
```

**NOTE**

You can also double click the file START-R_analyzer.sh and START-R_viewer.sh. In this situation a small manipulation is required (only once). In the Finder, right-click the file START-R_analyzer.sh (idem for START-R_viewer.sh) and select "Open with" and then "Other...".

You can select the application you want the file to be execute with. In this case it should be the Terminal. To be able to select the Terminal, you have to switch from "Recommended Applications" to "All Applications"  (the Terminal.app application can be found in the Utilities folder).

Check "Always Open With" and after clicking OK you should be able to execute you SHELL script by simply double-clicking it.

#### Linux installation

**In command line**

[Open a terminal](https://linuxconfig.org/how-to-open-a-terminal-on-ubuntu-bionic-beaver-18-04-linux) and run these commands:

```
git clone https://github.com/thomasdenecker/START-R.git
cd START-R
sudo ./INSTALLATION_LINUX.sh
```
Once the installation is complete, use this command to launch START-R analyzer:
```
sudo ./START-R_analyzer.sh
```

and this command to launch START-R viewer
```
sudo ./START-R_viewer.sh
```

### START-R application utilisation

Double click on START-R file (Windows / MacOS X) or launch the command line (Linux) and open your internet browser, typing the following url: http://localhost:3838/ for START-R analyzer and http://localhost:3839/ for START-R viewer and it should work.

**NOTE** (MAC users) : You may need to repeat the same manipulation as for the installation file (only once).

If you want to test the START-R suite on specific datasets for a differential analysis, please use the GEO datasets [GSM2111308](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2111308) for U2OS cells and [GSM2111313](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2111313) for K562 cells.

### START-R analyzer results

When an analysis is done with START-R analyzer, all results for this analysis are available in the folder START-R_analyzer/Outputs. The folder name begins with the date of analysis and successive numbers (for example, the first analysis will be named 20190611_1, the second analysis 20190611_2 ...)

## How to use 

To help you in the use of START-R, we have written a [wiki](https://github.com/thomasdenecker/START-R/wiki). If you have any questions or problems, do not hesitate to post an [issue](https://github.com/thomasdenecker/START-R/issues/new/). 

## Development

### Launch in debug mode

During development, you will probably need to get all messages (errors, warnings and notifications) in the R terminal. The following command launches the application and generates a log file in the application folder. To find the path to the application, you can look in the launch file.

START-R analyzer
```
docker run -ti --rm -p 3838:3838 -v YOUR_APPLICATION_PATH:/var/log/shiny-server -v YOUR_APPLICATION_PATH/START-R_analyzer:/srv/shiny-server tdenecker/start-r
```

START-R viewer
```
docker run -ti --rm -p 3839:3838 -v YOUR_APPLICATION_PATH:/var/log/shiny-server -v YOUR_APPLICATION_PATH/START-R_viewer:/srv/shiny-server tdenecker/start-r
```

### Connect to a R session

```
docker run -ti --rm -p 3839:3838 -v YOUR_APPLICATION_PATH:/srv/shiny-server  tdenecker/start-r R
```

**Warning**: nothing is saved in this session (package installation, ...)

### Remove folder (Only for linux user) 

To delete an analysis folder, you must use the following command :
```
sudo rm -rf dirName
```
## References
START-R use R packages. You will find below the list of packages and the installed versions in the Docker image: 

- Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan
McPherson (2017). shiny: Web Application Framework for R. R package
version 1.0.5. https://CRAN.R-project.org/package=shiny

- Thomas Lin Pedersen (2016). shinyFiles: A Server-Side File System
Viewer for Shiny. R package version 0.6.2.
https://CRAN.R-project.org/package=shinyFiles

- Dean Attali (2018). shinyjs: Easily Improve the User Experience of
Your Shiny Apps in Seconds. R package version 1.0.
https://CRAN.R-project.org/package=shinyjs

- Winston Chang (2016). shinythemes: Themes for Shiny. R package
version 1.1.1. https://CRAN.R-project.org/package=shinythemes

- Carson Sievert, Chris Parmer, Toby Hocking, Scott Chamberlain,
Karthik Ram, Marianne Corvellec and Pedro Despouy (2017). plotly:
Create Interactive Web Graphics via 'plotly.js'. R package version
4.7.1. https://CRAN.R-project.org/package=plotly

- Dean Attali (2017). colourpicker: A Colour Picker Tool for Shiny and
for Selecting Colours in Plots. R package version 1.0.
https://CRAN.R-project.org/package=colourpicker

- Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and
Smyth, G.K. (2015). limma powers differential expression analyses for
RNA-sequencing and microarray studies. Nucleic Acids Research 43(7),
e47.

- Venkatraman E. Seshan and Adam Olshen (). DNAcopy: DNA copy number
data analysis. R package version 1.44.0.

- SNPchip: R classes and methods for SNP array data R.B. Scharpf and G.
Parmigiani and J. Pevsner and I. Ruczinski 2007, Bioinformatics, Vol.
23, 627-628

- Hans W. Borchers (2018). pracma: Practical Numerical Math Functions.
R package version 2.1.4. https://CRAN.R-project.org/package=pracma

- Matt Dowle and Arun Srinivasan (2017). data.table: Extension of
'data.frame'. R package version 1.10.4-3.
https://CRAN.R-project.org/package=data.table

- RStudio and Inc. (2017). htmltools: Tools for HTML. R package version
0.3.6. https://CRAN.R-project.org/package=htmltools

- Marek Walesiak and Andrzej Dudek (2017). clusterSim: Searching for
Optimal Clustering Procedure for a Data Set. R package version
0.47-1. https://CRAN.R-project.org/package=clusterSim

- John Fox and Sanford Weisberg (2011). An {R} Companion to Applied
Regression, Second Edition. Thousand Oaks CA: Sage. URL:
http://socserv.socsci.mcmaster.ca/jfox/Books/Companion

- R Core Team (2015). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
URL https://www.R-project.org/.

Note: You will find in the application help in the choice of methods based on these packages. These helps are directly extracted from the packages. 

## Citation
If you use START-R project, please cite our paper :

**Efficient, quick and easy-to-use DNA replication timing analysis with START-R suite**

Djihad Hadjadj, Thomas Denecker, Eva Guérin, Su-Jung Kim, Fabien Fauchereau, Giuseppe Baldacci, Chrystelle Maric, Jean-Charles Cadoret

NAR Genomics and Bioinformatics, Volume 2, Issue 2, June 2020, lqaa045, https://doi.org/10.1093/nargab/lqaa045

## Contributing

Please, see the [CONTRIBUTING](CONTRIBUTING.md) file.

## Contributor Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](http://contributor-covenant.org/). By participating in this project you
agree to abide by its terms. See [CODE_OF_CONDUCT](CODE_OF_CONDUCT.md) file.

## License

START R is released under the BSD-3 License. See the bundled [LICENSE](LICENSE)
file for details.
Code of Conduct
===============

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of
experience, nationality, personal appearance, race, religion, or sexual identity
and orientation.

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
reported by contacting the project leader at thomas.denecker@u-psud.fr. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an
incident. Further details of specific enforcement policies may be posted
separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
Contributing Guidelines
=======================

First of all, **thank you** for contributing, **you are awesome**!

Before starting, you should read, agree to, and follow these three things:

* [How to contribute?](#how-to-contribute)
* [Pull Request Guidelines](#pull-request-guidelines)
* [Code of Conduct](CODE_OF_CONDUCT.md)

---

## How to contribute?

> You might be interested in [how to communicate with GitHub
> labels](https://tailordev.fr/blog/2016/09/27/communication-with-github-issues-1/).

### Report Bugs

Report bugs at: https://github.com/thomasdenecker/START-R/issues/new.

When reporting a bug, please include:

* Any details about your local setup which might be helpful in troubleshooting
* Detailed steps to reproduce the bug. Where possible, please write a test case

If you are not able to do that, that's fine! Open an issue anyway and let us
know as much information as you can. We will get back to you to determine the
problem, and (hopefully) fix it.

### Fix Bugs

Check out the [open bugs](https://github.com/thomasdenecker/START-R/issues) - anything
tagged with the **[easy pick]** label could be a good choice for newcomers (and
we are willing to help you).

We have two kind of bugs: **[critical]** and **[bug]**. We tend to fix critical
bugs as soon as possible. Feel free to come up with a patch before we do though!

### Implement Features

Look through the GitHub issues for features. Anything tagged with
**[improvement]** or **[feature]** is open to whoever wants to implement it.

If the issue is unclear or you are not sure what is expected, ask for more
information by commenting on the issue.

### Submit Feedback

Any issue with the **[question]** label is open for feedback, so feel free to
share your thoughts with us!

The best way to send feedback is to [create a new
issue](https://github.com/thomasdenecker/START-R/issues/new) on GitHub.

If you are proposing a feature:

* Explain how you envision it working. Try to be as detailed as you can
* Try to keep the scope as narrow as possible. This will help make it easier to
  implement
* Feel free to include any code you might already have, even if it is
  just a rough idea. This is a volunteer-driven project, and contributions are
  welcome :)

Your issue will be flagged as **[feature request]** first, and if we agree on
it, we will label it as **[feature]**, meaning it has been accepted and the
feature will eventually be added to the project.

## Pull Request Guidelines

Here are a few rules to follow in order to make code reviews and discussions go
more smoothly before maintainers accept and merge your work:

* you MUST test in Docker START-R
* you MUST comment your code
* you SHOULD write documentation

Please, write [commit messages that make
sense](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html),
and [rebase your branch](http://git-scm.com/book/en/Git-Branching-Rebasing)
before submitting your Pull Request.

You may be asked to [squash your
commits](http://gitready.com/advanced/2009/02/10/squashing-commits-with-rebase.html)
too. This is to "clean" your Pull Request before merging it (we don't want
commits such as `fix tests`, `fix 2`, `fix 3`, etc.).

Also, while creating your Pull Request on GitHub, you MUST write a description
which gives the context and/or explains why you are creating it.

For further information about creating a Pull Request, please read [this blog
post](http://williamdurand.fr/2013/11/20/on-creating-pull-requests/).

Thank you!
