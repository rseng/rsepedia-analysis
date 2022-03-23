# Security Policy

## Supported Versions

As we have limited resource, we only support the last minor versions
of MetaWards with security updates. For example, if the current version
is 1.1, then only versions 1.1 and 1.0 will have security updates.

| Version | Supported          |
| ------- | ------------------ |
| 1.0.x   | :white_check_mark: |
| <= 1.0.x| :x:                |

## Reporting a Vulnerability

Please report a vulnerability by emailing the Bristol RSE development
team at ask-rse@bristol.ac.uk. We will review your report as quickly
as we can, and will aim to respond to you within seven days with
the outcome. We really appreciate your help and your patience.
# Developing and packaging MetaWards

MetaWards is now fully tested and deployed using GitHub actions.

Full instructions for packaging are
[available here](https://metawards.org/packaging.html).
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
reported by contacting the project team at Christopher.Woods@bristol.ac.uk. All
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
# Installing MetaWards

MetaWards is available as a pip package or as a R module.

Full instructions for installing are
[available here](https://metawards.org/install.html).
# MetaWards

[![Build status](https://github.com/metawards/MetaWards/workflows/Build/badge.svg)](https://github.com/metawards/MetaWards/actions?query=workflow%3ABuild) [![PyPI version](https://badge.fury.io/py/metawards.svg)](https://pypi.python.org/pypi/metawards) [![Downloads](https://pepy.tech/badge/metawards)](https://pepy.tech/project/metawards) [![DOI](https://joss.theoj.org/papers/10.21105/joss.03914/status.svg)](https://doi.org/10.21105/joss.03914) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5562737.svg)](https://doi.org/10.5281/zenodo.5562737)

* For the most accurate and up to date information please [visit the project website](https://metawards.org).
* For an overview of features please [visit the features page](https://metawards.org/features).
* Please take a read of our [Journal of Open Source Software paper](https://joss.theoj.org/papers/10.21105/joss.03914)
* Please cite according to the CITATION.cff file in this repo.

## Scientific Background

MetaWards implements a stochastic metapopulation model of disease transmission. It can scale from modelling local transmission up to full national- or international-scale metapopulation models.

Please follow the [quick start guide](https://metawards.org/quickstart) to see how to quickly get up and running using MetaWards to model your own custom disease or metapopulation model.

It is was originally developed to support modelling of disease transmission in Great Britain. The complete model description and the original C code are described here;

*  *"The role of routine versus random movements on the spread of disease in Great Britain"*, Leon Danon, Thomas House, Matt J. Keeling, Epidemics, December 2009, 1 (4), 250-258; DOI:[10.1016/j.epidem.2009.11.002](https://doi.org/10.1016/j.epidem.2009.11.002)

*  *"Individual identity and movement networks for disease metapopulations"*, Matt J. Keeling, Leon Danon, Matthew C. Vernon, Thomas A. House Proceedings of the National Academy of Sciences, May 2010, 107 (19) 8866-8870; DOI:[10.1073/pnas.1000416107](https://doi.org/10.1073/pnas.1000416107)

In this model, the population is divided into electoral wards. Disease transmission between wards occurs via the daily movement of individuals. For each ward, individuals contribute to the *force of infection* (FOI) in their *home* ward during the night, and their *work* ward during the day.

This model was recently adapted to model CoVID-19 transmission in England and Wales, with result of the original C code published here;

* *"A spatial model of CoVID-19 transmission in England and Wales: early spread and peak timing"*, Leon Danon, Ellen Brooks-Pollock, Mick Bailey, Matt J Keeling, Philosophical Transactions of the Royal Society B, 376(1829); DOI:[10.1098/rstb.2020.0272](https://doi.org/10.1098/rstb.2020.0272)

This Python code is a port which can identically reproduce the outputs from the original C code as used in that work. This Python code has been optimised and parallelised, with additional testing added to ensure that development and scale-up of MetaWards has been robustly and efficiently conducted.

## Program Info

The package makes heavy use of [cython](https://cython.org) which is used with [OpenMP](https://openmp.org) to compile bottleneck parts of the code to parallelised C. This enables this Python port to run at approximately the same speed as the original C program on one core, and to run several times faster across multiple cores.

The program compiles on any system that has a working C compiler that supports OpenMP, and a working Python >= 3.7. This include X86-64 and ARM64 servers.

The software supports running over a cluster using MPI (via [mpi4py](https://mpi4py.readthedocs.io/en/stable/)) or via simple networking (via [scoop](http://scoop.readthedocs.io)).

Full instructions on how to use the program, plus example job submission scripts can be found on the [project website](https://metawards.org).

## Installation

[Full installation instructions are here](https://metawards.org/install.html).

Binary packages are uploaded to [pypi](https://pypi.python.org/pypi/metawards) for Windows, OS X and Linux (manylinux). The easiest way to install is to type in the console:

```
pip install metawards
```

(this assumes that you have pip installed and are using Python 3.7 or above - if this doesn't work please follow the [full installation instructions](https://metawards.org/install.html)).

Alternatively, you can also install from within R (or RStudio) by typing;

```
library(devtools)
install_github("metawards/rpkg")
metawards::py_install_metawards()
```

But, as you are here, I guess you want to install the latest code from GitHub ;-)

To do that, first clone and install the requirements;

```
git clone https://github.com/metawards/MetaWards
cd MetaWards
pip install -r requirements.txt
pip install -r requirements-dev.txt
```

Next, you can make using the standard Python setup.py script route.

```
CYTHONIZE=1 python setup.py build
CYTHONIZE=1 python setup.py install
```

Alternatively, you can also use the makefile, e.g.

```
make
make install
```

(assuming that `python` is version 3.7 or above)

You can run tests using pytest, e.g.

```
METAWARDSDATA="/path/to/MetaWardsData" pytest tests
```

or you can type

```
make test
```

You can generate the docs using

```
make doc
```

## Running

* [A quick start guide is here](https://metawards.org/quickstart)
* [A complete tutorial is here](https://metawards.org/tutorial)
* [Full usage instructions are here](https://metawards.org/usage.html)

You can either load and use the Python classes directly, or you can run the `metawards` front-end command line program that is automatically installed.

```
metawards --help
```

will print out all of the help for the program.

### Running an ensemble

This program supports parallel running of an ensemble of jobs using [multiprocessing](https://docs.python.org/3.7/library/multiprocessing.html) for single-node jobs, and [mpi4py](https://mpi4py.readthedocs.io/en/stable/) or [scoop](http://scoop.readthedocs.io) for multi-node cluster jobs.

Note that mpi4py and scoop are not installed by default, so you will need to install them before you run on a cluster (e.g. `pip install mpi4py` or `pip install scoop`).

[Full instructions for running on a cluster are here](https://metawards.org/cluster_usage.html)

## History

This is a Python port of the [MetaWards](https://github.com/ldanon/MetaWards) package originally written by Leon Danon. This port has been performed with Leon's support by the [Bristol Research Software Engineering Group](https://www.bristol.ac.uk/acrc/research-software-engineering).
This is the numpy.random.binomial function that I've copied out.
This reproduces the numpy random binomial numbers, and means that
I can now compile and link just this code to metawards, thereby
avoiding having to compile and link in all of numpy and/or gsl

---
name: Feature branch
about: Let developers know what you are doing in a feature branch
title: "[FEATURE BRANCH] - I am working on feature [XXX] in branch [feature_XXX]..."
labels: feature-branch

---

**What is the feature you are working on?**
A clear and concise description of what you are developing

**What is the use case for this feature?**
A clear and concise use case, ideally including an idea for the tutorial
chapter that will describe use of this feature.

**Core code or plugin?**
Will this feature change the core code? Or can it be implemented entirely
in one of the plugins, e.g. iterator, extractor etc.? Ideally most new
features should be in plugins.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions that exist
already in code, and why a new feature is needed.

**Additional information**
Add any additional information, e.g. is this a going to change the API
at all, which milestone release are you targetting, is anyone else
involved in developing this feature, does this depend on any other work,
or are there any other features depending on this work?
---
name: Bug report
about: Create a report to help us improve
title: "[BUG] - MetaWards is not working because..."
labels: bug
assignees: chryswoods

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Run the command '....'
2. Input data is located here '....'
3. See error 'Copy full error here'

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Environment (please complete the following information):**
 - OS: [e.g. Windows 10, Ubuntu 16.04, OS X Catalina)
 - Python: [e.g. 3.7. 3.8, anaconda]
 - Version: [e.g. metawards 1.0.0]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for MetaWards
title: "[FEATURE REQUEST] - I'd like MetaWards to..."
labels: enhancement
assignees: chryswoods

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Documentation or tutorial error
about: Create a report to help us improve the documentation or tutorial
title: "[DOCS] - MetaWards docs should be improved..."
labels: documentation
assignees: chryswoods

---

**Describe the documentation or tutorial problem**
A clear and concise description of the documentation of tutorial problem

**Location**
Which page (or pages) need work?

**Expected content**
A clear and concise description of what should be in the documentation or tutorial.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Environment (if applicable, please complete the following information):**
 - OS: [e.g. Windows 10, Ubuntu 16.04, OS X Catalina)
 - Python: [e.g. 3.7. 3.8, anaconda]
 - Version: [e.g. metawards 1.0.0]

**Additional context**
Add any other context about the problem here.
---
name: Something else?
about: Ask something else that doesn't fit into the other templates
title: "[MISC] - Something else about MetaWards..."
labels: miscellaneous
assignees: chryswoods

---

**Describe your request in detail below**
---
title: "Reagional Model Plots"
output: html_notebook
---
```{r}

library(tidyverse)
library(tidyverse)
library(gghighlight)
library(purrr)
library(plotly)
library(cowplot)
library(scales)

```


```{r}
wardlookup<-read.csv('~/GitHub/MetaWards/data/2011/WardsProcessing/Ward_Lookup.csv') 
lad2region<-read.csv('~/GitHub/MetaWards/data/2011/WardsProcessing/Output_Area_2011_to_Builtup_Area_Subdivision_to_Builtup_Area_to_Local_Authority_District_to_Region_December_2011_Lookup_in_England_and_Wales.csv') %>% 
  group_by(LAD11CD) %>% 
  summarise(Region=unique(RGN11NM))

wardlookupregion<-wardlookup %>% 
  inner_join(.,lad2region, by=("LAD11CD"))


allinfections=read.table(file='Testing/BrigthonNewOutput/ForMattData.dat',sep = ' ')


allinfections %>% # make long data frame
  mutate(Time=row_number()) %>% # count up time, as row number
  pivot_longer(-Time) %>% 
  mutate(Ward=as.integer(str_remove(name,'V'))) %>% 
  select(-name)->inf_long # rename name to Ward integers for easier matching


wardlookupregion %>% 
  inner_join(.,inf_long, by=c('FID'='Ward'),all.y=T,all.x=F) %>%
  group_by(Region,Time)%>% 
  summarise(Cases=sum(value)) -> region_inf

region_inf%>%
  group_by(Region) %>%
  summarise(PeakTime=which.max(Cases),PeakCases=max(Cases))


```
# Plot region epidemics in ggplot
```{r}


region_inf %>% 
  filter(Time<250) %>% 
  ggplot(aes(x=Time,y=Cases,colour=Region)) +
  geom_line(size=2,alpha=0.5)+
  theme_minimal_grid()



```
---
title: "Seasonality Model Plots"
output: html_notebook
---
```{r}

library(tidyverse)
library(tidyverse)
library(gghighlight)
library(purrr)
library(plotly)
library(cowplot)
library(scales)

```

# UV Curves



```{r UV}

pathname = 'Testing'

dir(pathname,pattern='PlayInfections.dat',recursive=TRUE)->pfilenames
dir(pathname,pattern='WorkInfections.dat',recursive=TRUE)->wfilenames


dataUVP <- tibble(filename = pfilenames) %>% # create a data frame
                                                  # holding the file names
  mutate(file_contents = map(filename,            # read files into a new data column
                             ~readIncidenceRecovered(., pathname))) %>% # use the function written above to process files
  mutate(Seasonal=str_remove(str_split(filename,'/',simplify = TRUE)[,1],'UV')) %>% 
  #mutate(Run=str_split(filename,'/',simplify = TRUE)[,2]) %>% 
  select(-filename) %>% 
  unnest(cols=file_contents)
  

# WORK infected
dataUVW <- tibble(filename = wfilenames) %>% # create a data frame
                                                  # holding the file names
  mutate(file_contents = map(filename,            # read files into a new data column
                             ~readIncidenceRecovered(., pathname))) %>% # use the function written above to process files
  mutate(Seasonal=str_remove(str_split(filename,'/',simplify = TRUE)[,1],'UV')) %>% 
  #mutate(Run=str_split(filename,'/',simplify = TRUE)[,2]) %>% 
  select(-filename) %>% 
  unnest(cols=file_contents)

bind_cols(dataUVW,dataUVP) %>% 
  mutate(Incidence=Incidence+Incidence1) %>% 
  mutate(Recovered=Recovered+Recovered1) %>% 
  select(Time,Incidence,Recovered,Seasonal)->dataUV

dataUV %>%
  filter(Time<650) %>% 
  mutate(Date=as.Date(Time,origin="2020-2-10")) %>% 
  ggplot(aes(x=Date,colour=Seasonal))+
  geom_line(aes(y=Incidence),alpha=1,size=1.5)+
#  geom_line(aes(y=Recovered),alpha=0.1)+
  theme_half_open()+
  theme(legend.position = "top")+
 #scale_y_continuous(labels = unit_format(unit = "M"))+
  xlab('Time')+
  scale_x_date(breaks=date_breaks("3 months"),
               labels = date_format("%b"))->p


```
```{r}

dataUV %>% 
  group_by(Seasonal) %>% 
  summarise(Peak=which.max(Incidence),PeakIncidence=max(Incidence),AttackRate=100*max(Recovered)/56082077)
#->p
#ggplotly(p)



```

```{r}
t=1:365
#theta=seq(0,2*pi,pi/100)

seasonality<-function(t,m){
  x=(1-m/2.0)+m*cos(2*pi*t/365)/2.0
  return(tibble(t,x))
  } 


tibble(m=c(0,0.25,0.5,0.75,1)) %>% 
  mutate(seasoncols=map(m,~seasonality(t,.))) %>% 
  unnest(cols=seasoncols) %>% 
  mutate(Seasonal=as.factor(m))->seasoncurves

seasoncurves %>% 
  mutate(Date=as.Date(t,origin="2020-1-1")) %>% 
  ggplot(aes(x=Date,y=x,colour=Seasonal))+
  geom_line(size=2,alpha=0.8)+
  xlab('')+
  ylab('Seasonal Scaling')+
  theme_minimal_grid()+
  theme(legend.position = 'none')+
  panel_border()+
  scale_y_continuous(position='right')+
  scale_x_date(breaks=date_breaks("3 months"),
               labels = date_format("%b"))->inset


inset
# put the plots together. 

ggdraw(p)+
  draw_plot(inset,.65,.5,.35,.4)

```
---
title: "OutpuMaps"
author: "Leon Danon"
date: "07/02/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Spatial Analysis

```{r}

library(tidyverse)
library(sf)
library(rgdal)
library(viridis)
library(mapview)

WardsEW <- read_sf(dsn='~/GitHub/MetaWards/data/2011/shapefile/Wards_December_2011_Boundaries_EW_BFC',layer = 'Wards_December_2011_Boundaries_EW_BFC')

st_centroid(WardsEW)->wardCentroids

wardlookup<-read.csv('~/GitHub/MetaWards/data/2011/WardsProcessing/Ward_Lookup.csv') 

wardCentroids %>% 
  inner_join(.,wardlookup,by=c('wd11cd'='WD11CD'))->wardslookup2 

#mapview(wardslookup2,zcol='FID')


```


```{r}

allinfectionsShef=read.table(file='ForMattData.dat',sep = ' ')

allinfectionsShef %>% # make long data frame
  mutate(Time=row_number()) %>% # count up time, as row number
  pivot_longer(-Time) %>% 
  mutate(Ward=as.integer(str_remove(name,'V'))) %>% 
  select(-name)->inf_longShef # rename name to Ward integers for easier matching



inf_longShef %>% 
  group_by(Ward) %>%               # by ward
  filter(cumsum(value) < 1) %>%    # cumulative sum > 1 means it's just arrived/first infection
  mutate(Arrival=max(Time)) %>%    # find the maximum before that
  filter(row_number()==1) %>%      # take the first row of group_by
  select(Ward,Arrival) %>%    # # select columns and rename and save in a data frame
  ungroup() -> arrivals
      
wardslookup2 %>% inner_join(.,arrivals, by=c('FID'='Ward')) %>% 
  mutate(AltArrival=log(1/Arrival)*100)->wardslookuparrivalShef

mapview(wardslookuparrivalShef,zcol='AltArrival',cex = 'AltArrival')


# plotme %>% 
#   ggplot(aes(x=Time,y=value,colour=name))+
#   geom_line()

```

```{r}

allinfectionsLon=read.table(file='Testing/London/0/ForMattData.dat',sep = ' ')

allinfectionsLon %>% # make long data frame
  mutate(Time=row_number()) %>% # count up time, as row number
  pivot_longer(-Time,values_to = 'Infected') %>% 
  mutate(Ward=as.integer(str_remove(name,'V'))) %>% 
  select(-name)->inf_longLon # rename name to Ward integers for easier matching

inf_longLon %>% 
  group_by(Ward) %>%  # by ward
  filter(cumsum(Infected) < 1) %>%    # cumulative sum > 1 means it's just arrived/first infection
  mutate(Arrival=max(Time)) %>%    # find the maximum before that
  filter(row_number()==1) %>%      # take the first row of group_by
  select(Ward,Arrival) %>%    # # select columns and rename and save in a data frame
  ungroup() -> arrivalsLon

inf_longLon %>% 
  group_by(Ward) %>% 
  slice(which.max(Infected))->TimeToPeak
```


```{r}

TimeToPeak %>% 
  ggplot(aes(x=Time))+
  geom_histogram()

```

```{r}
      
wardslookup2 %>% 
  inner_join(.,arrivalsLon, by=c('FID'='Ward')) %>% 
  mutate(AltArrival=log(1/Arrival)*100)->wardslookuparrivalLon

wardslookuparrivalLon %>% 
  inner_join(.,TimeToPeak, by=c('FID'='Ward'))%>% 
  mutate(AltTime=Time/100)->wardslookuparrivaltimetopeakL

mapview(wardslookuparrivaltimetopeakL,zcol='Time',cex = 'AltTime')


# TODO:
# Time to peak. DONE 
# Time to 10 cases. 
# Lookup from Postcode to Ward/OA.
```

---
title: "Plotting Model Output"
output: html_notebook
---

```{r}
library(tidyverse)
library(gghighlight)
library(purrr)
library(plotly)
library(cowplot)
library(scales)

```


# Read in model output files


```{r}

# this bit of code is now obsolete. 
#pathname='Locations'
#dir(pathname,pattern='TotalInfections.dat',recursive=TRUE)->infectionfiles

# readinfections<-function (fname,fpath){        
# #  Function to read in TotalInfections.dat' type files (one column only of infecteds, Time is row number)
#   d<-read_csv(file.path(fpath,fname),col_names=FALSE,cols(X1=col_integer())) %>% # 
#     mutate(Time=row_number()) %>% 
#     select(Time,Infecteds=X1)
#    return(d)
# }
  
# data <- tibble(filename = infectionfiles) %>% # create a data frame
                                                  # holding the file names
  # mutate(file_contents = map(filename,            # read files into a new data column
  #                            ~readinfections(., pathname))) %>% # use the function written above to process files
  # mutate(Location=str_split(filename,'/',simplify = TRUE)[,1]) %>% 
  # mutate(Run=str_split(filename,'/',simplify = TRUE)[,2]) %>% 
  # select(-filename) %>% 
  # unnest(cols=file_contents)
  


# data %>%
#   mutate(Infected=Infecteds/10^6) %>% 
#   filter(Time<250) %>% 
#   ggplot(aes(x=Time,y=Infected,colour=Run))+
#   geom_line()+
#   gghighlight(~Location)+
#   facet_wrap(~Location)+
#   theme_minimal_grid()+
#   theme(strip.text.x = element_blank())+
#   scale_y_continuous(labels = unit_format(unit = "M"))+
#   xlab('Time (days)')
# 
# 
# data %>% 
#   group_by(Run) %>% 
#   summarise(TimeToPeak=which.max(Infecteds)) %>% summary()
# 
# 
# 
#   ggplot(aes(x=TimeToPeak))+
#   geom_histogram()

```


# Number of Infected Wards

```{r}

read.csv(file='NumberWardsInfected.dat',sep=' ') %>% 
  mutate(Time=X0) %>% 
  select(-X,-X0,-X0.4) ->wardcurves


wardcurves %>% 
  pivot_longer(-Time,names_to = 'Classes',values_to = 'NumberWards')%>% 
  ggplot(aes(x=Time,y=NumberWards,colour=Classes))+
  geom_line()
```

# Incidence curves

```{r}




readincidence<-function (fname,fpath){        # read and process a single file
#  Function to read in Work/PlayInfections.dat' type files (many columns, one for each class only of infecteds, Time is row number)
  d<-read_delim(file.path(fpath,fname),col_names=FALSE,delim=' ',col_types = "iiiiiii")%>% select(Time=X1,Incidence=X4)
   return(d)
}


readIncidenceRecovered<-function (fname,fpath){        # read and process a single file
#  Function to read in Work/PlayInfections.dat' type files (many columns, one for each class only of infecteds, Time is row number)
  d<-read_delim(file.path(fpath,fname),col_names=FALSE,delim=' ',col_types = "iiiiiii")%>% select(Time=X1,Incidence=X4,Recovered=X6)
   return(d)
}
 

readPrevalence<-function (fname,fpath){        # read and process a single file
#  Function to read in Work/PlayInfections.dat' type files (many columns, one for each class only of infecteds, Time is row number)
  d<-read_delim(file.path(fpath,fname),col_names=FALSE,delim=' ',col_types = "iiiiiii")%>%  
    transmute(Time=X1,Prevalence=X4+X5)
   return(d)
}

```


```{r}

pathname='Locations'
dir(pathname,pattern='PlayInfections.dat',recursive=TRUE)->pfilenames
dir(pathname,pattern='WorkInfections.dat',recursive=TRUE)->wfilenames


#PLAY infecteds
dataIP <- tibble(filename = pfilenames) %>% # create a data frame
                                                  # holding the file names
  mutate(file_contents = map(filename,            # read files into a new data column
                             ~readIncidenceRecovered(., pathname))) %>% # use the function written above to process files
  mutate(Location=str_split(filename,'/',simplify = TRUE)[,1]) %>% 
  mutate(Run=str_split(filename,'/',simplify = TRUE)[,2]) %>% 
  select(-filename) %>% 
  unnest(cols=file_contents)
  

# WORK infected
dataIW <- tibble(filename = wfilenames) %>% # create a data frame
                                                  # holding the file names
  mutate(file_contents = map(filename,            # read files into a new data column
                             ~readIncidenceRecovered(., pathname))) %>% # use the function written above to process files
  mutate(Location=str_split(filename,'/',simplify = TRUE)[,1]) %>% 
  mutate(Run=str_split(filename,'/',simplify = TRUE)[,2]) %>% 
  select(-filename) %>% 
  unnest(cols=file_contents)

bind_cols(dataIW,dataIP) %>% 
  mutate(Incidence=Incidence+Incidence1) %>% 
  mutate(Recovered=Recovered+Recovered1) %>% 
  select(Time,Incidence,Recovered,Location,Run)->dataIR

dataIR %>%
  mutate(Date=as.Date(Time,origin='2020-02-10')) %>% 
  #mutate(Incidence=Incidence/10^6) %>% 
  #mutate(Recovered=Recovered/10^6) %>% 
  filter(Time<450) %>% 
  ggplot(aes(x=Date,group=Run))+
  geom_line(aes(y=Incidence),alpha=0.3,size=1,colour='red')+
#  geom_line(aes(y=Recovered),alpha=0.1)+
#  facet_wrap(~Location)+
  theme_minimal_grid()+
  theme(legend.position = "none")+
 #scale_y_continuous(labels = unit_format(unit = "M"))+
  xlab('Time')+
  scale_x_date(breaks=date_breaks("3 months"),
               labels = date_format("%b"))
#->p
#ggplotly(p)

dataIR %>% 
  group_by(Run) %>% 
  summarise(PeakIncidence=max(Incidence)) %>% 
  summary()
  
```

# Time To Peak Incidence

```{r}
dataIR %>% 
   group_by(Run) %>% 
   summarise(TimeToPeak=which.max(Incidence)) %>% 
   ggplot(aes(x=TimeToPeak))+
   geom_histogram(aes(y=..density..),bins=15)+
  geom_density()+
  theme_minimal_grid()+
  xlab('Time to Peak Incidence')+
  ylab('Density')

dataIR %>% 
   group_by(Run) %>% 
   summarise(TimeToPeak=which.max(Incidence)) %>% summary()
 
```






# Prevalence Curves (Correct) 

```{r}

dataPrevP <- tibble(filename = pfilenames) %>% # create a data frame
                                                  # holding the file names
  mutate(file_contents = map(filename,            # read files into a new data column
                             ~readPrevalence(., pathname))) %>% # use the function written above to process files
  mutate(Location=str_split(filename,'/',simplify = TRUE)[,1]) %>% 
  mutate(Run=str_split(filename,'/',simplify = TRUE)[,2]) %>% 
  select(-filename) %>% 
  unnest(cols=file_contents)

dataPrevW <- tibble(filename = wfilenames) %>% # create a data frame
                                                  # holding the file names
  mutate(file_contents = map(filename,            # read files into a new data column
                             ~readPrevalence(., pathname))) %>% # use the function written above to process files
  mutate(Location=str_split(filename,'/',simplify = TRUE)[,1]) %>% 
  mutate(Run=str_split(filename,'/',simplify = TRUE)[,2]) %>% 
  select(-filename) %>% 
  unnest(cols=file_contents)

bind_cols(dataPrevW,dataPrevP) %>% 
  mutate(Prevalence=Prevalence+Prevalence1) %>% 
  select(Time,Prevalence,Location,Run)->dataPrev


dataPrev %>% 
  filter(Time<450) %>% 
  ggplot(aes(x=Time,colour=Run))+
  geom_line(aes(y=Prevalence),alpha=0.2)+
  theme_minimal_grid()+
#  gghighlight(Location)+
#  facet_grid(~Location)+  
  theme(legend.position = "none")+
  xlab('Time (days)')


```

# Attack Rate

```{r}
library(kableExtra)
dataIR %>% 
  group_by(Run) %>% 
  summarise(AttackRate=max(Recovered)) %>%
  mutate(AttackRate2=100*AttackRate/56082077) %>% 
  ggplot(aes(x=AttackRate2))+
  geom_histogram(aes(y=..density..),bins=15)+
  geom_density()+
  xlab('Attack Rate (%)')+
  ylab('Density')+
  theme_minimal_grid()


dataIR %>% 
  group_by(Run) %>% 
  summarise(AttackNumbers=max(Recovered)) %>% 
  mutate(Deaths=AttackNumbers/100) %>% 
  mutate(AttackRate=100*AttackNumbers/56082077) %>%
  summary() 

```
---
title: "Meta Wards"
output: html_notebook
---

# Processing the census

The purpose of this file is to collect information, and data required for the Wards Model. We need:
1) origin destination data,


```{r libs}
library(tidyverse)
library(rgdal)
library(leaflet)
library(maptools)
library(broom)
# set factors to false
options(stringsAsFactors = FALSE)


```





# Processing Origin-Destination data


We need to generate a commuter matrix from Ward to Ward. Data are available from ONS. Below we detail the procedure if we need to repreat.



1. Download Ward lookups for Census Merged Wards to Original Wards from

https://geoportal.statistics.gov.uk/datasets/ward-to-census-merged-ward-to-local-authority-district-december-2011-lookup-in-england-and-wales

2. Then download the lookup from Output Areas to Ward level.

https://geoportal.statistics.gov.uk/datasets/output-area-to-ward-to-local-authority-district-december-2018-lookup-in-england-and-wales

3. Use left_join to combine them into a main lookup table.


```{r process }

OA2Ward = read.csv('~/GitHub/MetaWards/data/2011/Output_Area_to_Ward_to_Local_Authority_District_December_2011_Lookup_in_England_and_Wales.csv') # maps OA to Ward data

wardlookup<-read.csv('~/GitHub/MetaWards/data/2011/WardsProcessing/Ward_Lookup.csv')

OA2Ward %>% left_join(.,wardlookup,by=c('WD11CD'='WD11CD'))->OA2WardLookup

write.csv(file="WardLookupMaster.csv",OA2WardLookup)

```

4. Download the bulk data for Output Area to Output Area commuter numbers from:
https://www.nomisweb.co.uk/census/2011/bulk/rOD1

at Output Area level.

5. We need to aggregate up to Ward level, which we do below.

```{r process origin-destination}

#4

OA2OAmovements = read.csv(file='data/2011/WardsProcessing/bulk/wf01bew_oa_v1.csv', header=F)# Output area to output area commuting.

#5


OA2OAmovements %>%
  inner_join(.,OA2WardLookup, by = c('V1'='OA11CD')) %>% # match output area  to ward for first column
  inner_join(.,OA2WardLookup, by = c('V2'='OA11CD'))%>% # match output area  to ward for second column
  select(from=FID.x,to=FID.y,ObjectID=ObjectId.x,V3) %>%        # remove columns  that are not needed
  group_by(from,to) %>%
  summarize(total=sum(V3)) %>%
  ungroup->Ward2Ward



```

The Ward identifier is "FID", which seems to go from 1 to 8588, and is the index which will be used in the code.


```{r}
Ward2Ward %>% ggplot(aes(x=total)) +
  geom_freqpoly()+
  scale_x_log10()+
  scale_y_log10()+ # checking the histogram makes sense.


write.table(file='EW.dat',Ward2Ward,col.names = FALSE,row.names=F)


```

We also need a population size per ward, population working per ward and population not working per ward.



# Ward population sizes

Population sizes from here:
https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/wardlevelmidyearpopulationestimatesexperimental


```{r}

PopSizePerWard=read.csv('data/2011/PopSizePerWard.csv') # table from above

wardlookup %>%
  inner_join(.,PopSizePerWard,by=c('WD11CD'='geography.code'))->wardlookupWithPop # join with population ward lookup tables

Ward2Ward %>%
  group_by(from) %>%
  summarise(WorkSize=sum(total)) ->WorkSize # sum up the populations commuting to give us the commuting size

write.table(file='WorkSize.dat',WorkSize,col.names = FALSE,row.names=F) # write WorkPopulation per ward to file

WorkSize %>%
  inner_join(.,wardlookupWithPop,by=c("from"="FID"))->wardlookupWithPop2 # join with ward lookup file

wardlookupWithPop2 %>%
  select(Ward=from,All=Variable..All.usual.residents..measures..Value,WorkSize) %>%
  mutate(PlaySize=All-WorkSize) %>%
  mutate(PlaySize=ifelse(PlaySize<0,0,PlaySize))->wardlookupWithPop3

wardlookupWithPop3%>%
  select(from=Ward,PlaySize) %>%
  write.table(.,file="PlaySize.dat",col.name=FALSE,row.names=FALSE)

```

# Non-commuter matrix

```{r}
wardlookupWithPop3%>%
  select(from=Ward,PlaySize)->
  play


Ward2Ward %>%
  group_by(from) %>%
  mutate(rate=total/sum(total)) %>%
  ungroup() %>%
  inner_join(.,play,by='from') %>%
  mutate(PlayTotal=round(rate*PlaySize))->Ward2WardAll

Ward2WardAll %>%
  select(from,to,rate) %>%
  write.table(.,file="PlayMatrix.dat",row.names=F,col.names=F)


```

# CCG analysis

```{r}
lsoa2ccg = read.csv('~/GitHub/MetaWards/data/2011/CCG/Lower_Layer_Super_Output_Area_2011_to_Clinical_Commissioning_Group_to_Local_Authority_District_April_2017_Lookup_in_England_Version_4.csv')

lsoa2ccg %>%
  group_by(CCG17CD) %>%
  count()



```


====================
Command line options
====================

The full set of command line options for metawards-plot is below;

.. program-output:: metawards-plot --help
==================
Packaging releases
==================

MetaWards is now fully tested and deployed using GitHub actions.
The development process should be;

* New features are developed on feature branches, called ``feature-{feature}``,
  either in the `main MetaWards repository <https://github.com/metawards/MetaWards>`__
  for authorised developers, or in personal forks for
  new developers.
* Bug fixes or issue fixes are developed on fix branches, called
  ``fix-issue-{number}`` (again in either the main repository or forks).
* Pull requests are issued from these branches to ``devel``. All merge conflicts
  must be fixed in the branch and all tests must pass before the pull
  request can be merged into ``devel``. NOTE THAT ONLY AUTHORISED
  DEVELOPERS CAN ACCEPT THE PULL REQUEST. Authorised developers will
  review the pull request as quickly as they can. They will be greatly
  helped if the feature is accompanied with tests, examples and/or tutorial
  instructions, and if most changes are confined to plugins
  (e.g. :mod:`metawards.iterators`, :mod:`metawards.extractors`,
  :mod:`metawards.mixers` or :mod:`metawards.movers`).

The result of this is that "devel" should contain the fully-working and
tested, and most up-to-date version of ``metawards``. However, this
version should not be used for production runs.

.. note::

  The group of developers authorised to have access to the
  `main MetaWards repository <https://github.com/metawards/MetaWards>`__
  and to accept pull requests is not fixed,
  and will evolve over time. If you wish to join this group then
  please complete the tutorial and then demostrate your commitment
  by submitting good issues and pull requests from
  a personal fork of the repository. Please get in touch if you find
  this difficult, or follow
  `this workshop <https://chryswoods.com/beginning_git>`__ if you need
  to learn how to use Git, GitHub, feature branching, merging, pull
  requests etc.

Defining a release
------------------

We will release ``metawards`` regularly. Releases aim to be backwards
compatible and capable of being used for production runs, at least for
the functionality that is fully described in the tutorial.

We use `semantic versioning <https://semver.org>`__ and take care
not to cause breaking changes in the public API. The public API
consists of the core :mod:`metawards` module only - it does not
cover :mod:`metawards.utils` or any plugins that are in
:mod:`metawards.iterators`, :mod:`metawards.extractors`,
:mod:`metawards.mixers` or :mod:`metawards.movers` (although,
the core plugin interface API will not change). We will
endeavour to retain backwards compatibility outside of
the core :mod:`metawards` module, but cannot guarantee it.

.. note::

  It is the job of the release manager (currently
  `chryswoods <https://github.com/chryswoods>`__) to decide when it is time
  to create a new release. If you are interested in helping join the release
  management group then please feel free to get in touch.

Creating a release
------------------

To create a release first checkout the "main" branch.

.. code-block:: bash

   git checkout main
   git pull

Next, merge in all changes from the "devel" branch.

.. code-block:: bash

   git pull origin devel

Next, update the :doc:`changelog` with details about this release. This
should include the link at the top of the release that shows the commit
differences between versions. This can be easily copied from a previous
release and updated, e.g.

::

  `0.11.2 <https://github.com/metawards/MetaWards/compare/0.11.1...0.11.2>`__ - May 11th 2020


could be changed to

::

  `0.12.0 <https://github.com/metawards/MetaWards/compare/0.11.2...0.12.0>`__ - May 18th 2020

when moving from the 0.11.2 to 0.12.0 release.

Now push this change back to GitHub, using;

.. code-block:: bash

   git push

This will trigger a CI/CD run which will build and test everything on
Windows, Mac and Linux for Python 3.7 and 3.8. Everything should work,
as "devel" should have been in a release-ready state.

Testing the packages
--------------------

`GitHub actions <https://github.com/metawards/MetaWards/actions>`__ will
produce the source and binary wheels for ``metawards`` on all supported
platforms. This will be in an artifact called ``dist`` which you should
download and unpack.

.. image:: images/github_artifacts.jpg
   :alt: Image of the GitHub Actions interface showing the dist artifact

You should unpack these into the ``dist`` directory, e.g.

.. code-block:: bash

   cd dist
   unzip ~/Downloads/dist.zip

This should result in six binary wheels and once source package, e.g.

::

    metawards-0.11.1+7.g52b3671-cp37-cp37m-macosx_10_14_x86_64.whl
    metawards-0.11.1+7.g52b3671-cp37-cp37m-manylinux1_x86_64.whl
    metawards-0.11.1+7.g52b3671-cp37-cp37m-win_amd64.whl
    metawards-0.11.1+7.g52b3671-cp38-cp38-macosx_10_14_x86_64.whl
    metawards-0.11.1+7.g52b3671-cp38-cp38-manylinux1_x86_64.whl
    metawards-0.11.1+7.g52b3671-cp38-cp38-win_amd64.whl
    metawards-0.11.1+7.g52b3671.tar.gz

Try to install the package related to you machine, just to double-check
that it is working, e.g.

.. code-block:: bash

   pip install ./metawards-0.11.1+7.g52b3671-cp37-cp37m-macosx_10_14_x86_64.whl
   cd ..
   pytest tests

Once it is working, remove these temporary packages from your ``dist`` folder,

.. code-block:: bash

   rm dist/*

Tagging a new release
---------------------

Now that you are happy that the release is ready, you can tag the new
version. Do this using the ``git tag`` command, e.g.

.. code-block:: bash

   git tag -a {VERSION} -m "{VERSION} release"

replacing ``{VERSION}`` with the version number. For this 0.12.0 release
the command would be;

.. code-block:: bash

   git tag -a 0.12.0 -m "0.12.0 release"

Next, push your tag to GitHub;

.. code-block:: bash

   git push --tags

The tag will be used by automatic versioning script to generate
the version numbers of the code. Building the package
(as happens below) will automatically update the _version.py
that is included in the package to tag versions.

This will also trigger a full CI/CD to test and build the new version.
Again, it should work as this tag was taken from your fully-tested
"main" branch.

Uploading packages to pypi
--------------------------

While you are waiting for the CI/CD GitHub Actions to complete, make sure
that your version of twine is fully up to date;

.. code-block:: bash

   pip install --upgrade twine

Once GitHub actions is complete, you will see that another build artifact
is ready for download. Download this and unpack it into your ``dist``
directory as before. You should now have a ``dist`` directory that
contains six binary wheels and one source package, named according to
the release version. For example, for the 0.11.2 release we had;

.. code-block:: bash

   $ ls dist
    metawards-0.11.2-cp37-cp37m-macosx_10_14_x86_64.whl
    metawards-0.11.2-cp37-cp37m-manylinux1_x86_64.whl
    metawards-0.11.2-cp37-cp37m-win_amd64.whl
    metawards-0.11.2-cp38-cp38-macosx_10_14_x86_64.whl
    metawards-0.11.2-cp38-cp38-manylinux1_x86_64.whl
    metawards-0.11.2-cp38-cp38-win_amd64.whl
    metawards-0.11.2.tar.gz

Now you can upload to pypi using the command;

.. code-block:: bash

   python3 -m twine upload dist/*

.. note::

    You will need a username and password for pypi and to have
    permission to upload code to this project. Currently only
    the release manager has permission. If you would like
    join the release management team then please get in touch.

Testing the final release
-------------------------

Finally(!) test the release on a range of different machines by logging
in and typing;

.. code-block:: bash

   pip install metawards=={VERSION}

replacing ``{VERSION}`` with the version number, e.g. for 0.11.2

.. code-block:: bash

   pip install metawards==0.11.2

Play with the code, run the tests and run some examples. Everything should
work as you have performed lots of prior testing to get to this stage.
============
Contributing
============

We welcome all helpful contributions to improving the quality of
MetaWards. Ways you can help include;

* Finding and fixing documentation bugs and typos
* Creating more tests and adding them to the pytest library in the
  ``tests`` directory
* Taking on some of the tasks in the :doc:`snaglist`.
* Porting and testing MetaWards on different computers
* Writing new features

We accept pull requests to the devel branch and are happy to discuss
ideas via
`GitHub issues on the repo <https://github.com/metawards/MetaWards/issues>`__.

When contributing, please keep in mind our
`code of conduct <https://github.com/metawards/MetaWards/blob/devel/CODE_OF_CONDUCT.md>`__.

Before contributing we encourage everyone to
:doc:`complete the tutorial <tutorial/index>`, as this gives a good
grounding in how the model works and how the code is laid out. If you have
any problems running the tutorial then please
`raise and issue <https://github.com/metawards/MetaWards/issues>`__.

Contributing new code
---------------------

We welcome developers who want to join us to help reduce bugs and add
new functionality. Please bear in mind the following;

* We use `semantic versioning <https://semver.org>`__ and take care
  not to cause breaking changes in the public API. The public API
  consists of the core :mod:`metawards` module only - it does not
  cover :mod:`metawards.utils` or any plugins that are in
  :mod:`metawards.iterators`, :mod:`metawards.extractors`,
  :mod:`metawards.mixers` or :mod:`metawards.movers` (although,
  the core plugin interface API will not change). We will
  endeavour to retain backwards compatibility outside of
  the core :mod:`metawards` module, but cannot guarantee it.

* The core :mod:`metawards` module and data structures are now
  quite fixed, and we endeavour not to make large or breaking changes.
  This means that new functionality should ideally be added via
  one of the pluging interfaces, e.g.
  :doc:`iterators <api/index_MetaWards_iterators>`,
  :doc:`extractors <api/index_MetaWards_extractors>`,
  :doc:`mixers <api/index_MetaWards_mixers>` and
  :doc:`movers <api/index_MetaWards_movers>`. If you can't fit your code
  into a plugin then please
  `raise an issue <https://github.com/metawards/MetaWards/issues>`__
  to discuss your idea with the core developers, so that a way
  forward can be found. We really appreciate your help, and want
  to make sure that your ideas can be included in the most compatible
  way in the code.

We've :doc:`added features <devsupport>` to ``metawards`` that we hope
will make it easier to develop new code. This includes tools to make
simplify profiling and "printf" debugging, and full pytest integration.
Please feel free to
`raise an issue <https://github.com/metawards/MetaWards/issues>`__ if there
is something else you think would help make development easier.

Finally, we have a :doc:`very detailed developer guide <development>` that
we hope will help you get up to speed with development. We strive to
be a helpful, friendly and welcoming community, so if you have any
questions or anything is not clear then please get in touch with
us by `raising an issue <https://github.com/metawards/MetaWards/issues>`__.

====================
Command line options
====================

The full set of command line options for metawards is below;

.. program-output:: metawards --help
=========
Snag list
=========

Below is a series of snags with the code that we'd like to fix, but
currently lack the personpower to tackle now. If you'd like to
try and fix one of the below then get in touch with us by
`raising an issue <https://github.com/metawards/MetaWards/issues>`__,
or forking the repo and issuing a pull request when you are done.

Ideally let us know so that we can create a feature branch for your
work, so that no-one else duplicates your effort.

Snags
-----

* [feature_improve_rng] - Random number generator is potentially quite slow

We are using the numpy random number generator in single threaded,
not-vectorised mode. There is a lot of scope for speed-up, so someone
who wants to take a look would be very welcome (and appreciated).

* [feature_multinomial] - Replace ran_binomial loop in advance_foi

We think that a ran_multinomial call can replace the slow and multiple
calls to ran_binomial in advance_foi. If the maths works out, then
this should be significantly quicker.

* [feature_bigtests] - Create a larger test suite full of examples

We'd like to create a repo full of examples of MetaWards runs, which can
serve as both a learning resource and also a fuller suite of regression
and integration tests. The aim would for pushes to this repo to trigger
runs of metawards against a lot of examples. This should be a separate
repo to the main code as the examples will take up a lot of space.

====================
Running on a cluster
====================

One of the reasons for this Python port is to make it easier to run
MetaWards analyses at scale on a HPC cluster. MetaWards supports
parallelisation using MPI (via `mpi4py <https://mpi4py.readthedocs.io>`__)
or simple networking (via `scoop <https://scoop.readthedocs.io>`__).

MetaWards will automatically detect most of what it needs so that you
don't need to write a complicated HPC job script.

MetaWards will look for a ``hostfile`` via either the PBS environment
variable of ``PBS_NODEFILE``, or the slurm ``SLURM_HOSTFILE``, or
for a ``hostfile`` passed directly via the ``--hostfile`` command
line argument.

It will then use the information combined there, together with the number of
threads per model run requested by the user, and the number of
cores per compute node (set in the environment variable
``METAWARDS_CORES_PER_NODE``, or passed as the command line option
``--cores-per-node``) to work out how many parallel scoop or MPI
processes to start, and will start those in a round-robin fashion
across the cluster. Distribution of work to nodes is via the
scoop or mpi4py work pools.

What this means is that the job scripts you need to write are very simple.

Example PBS job script
======================

Here is an example job script for a PBS cluster;

::

  #!/bin/bash
  #PBS -l walltime=01:00:00
  #PBS -l select=4:ncpus=64:mem=64GB
  # The above sets 4 nodes with 64 cores each

  # source the version of metawards we want to use
  # (assumes your python environments are in $HOME/envs)
  source $HOME/envs/metawards-0.6.0/bin/activate

  # change into the directory from which this job was submitted
  cd $PBS_O_WORKDIR

  # if you need to change the path to the MetaWardsData repository,
  # then update the below line and uncomment
  #export METAWARDSDATA="$HOME/GitHub/MetaWardsData"

  metawards --additional ExtraSeedsBrighton.dat \
            --input ncovparams.csv --repeats 8 --nthreads 16

The above job script will run 8 repeats of the adjustable parameter sets
in ``ncovparams.csv``. The jobs will be run using 16 cores per model run,
over 4 nodes with 64 cores per node (so 256 cores total, running
16 model runs in parallel). The runs will take only a minute or two
to complete, hence why it is not worth requesting more than one hour
of walltime.

The above job script can be submitted to the cluster using the PBS
``qsub`` command, e.g. if the script was called ``submit.sh``, then you
could type;

.. code-block:: bash

  qsub submit.sh

You can see the status of your job using

.. code-block:: bash

  qstat -n

Example slurm job script
========================

Here is an example job script for a slurm cluster;

::

  #!/bin/bash
  #SBATCH --time=01:00:00
  #SBATCH --ntasks=4
  #SBATCH --cpus-per-task=64
  # The above sets 4 nodes with 64 cores each

  # source the version of metawards we want to use
  # (assumes your python environments are in $HOME/envs)
  source $HOME/envs/metawards-0.6.0/bin/activate

  # if you need to change the path to the MetaWardsData repository,
  # then update the below line and uncomment
  #export METAWARDSDATA="$HOME/GitHub/MetaWardsData"

  metawards --additional ExtraSeedsBrighton.dat \
            --input ncovparams.csv --repeats 8 --nthreads 16

This script does the same job as the PBS job script above. Assuming
you name this script ``submit.slm`` you can submit this job using

.. code-block:: bash

  sbatch submit.slm

You can check the status of your job using

.. code-block:: bash

  squeue -u USER_NAME

where ``USER_NAME`` is your cluster username.
=========
MetaWards
=========

.. image:: https://github.com/metawards/MetaWards/workflows/Build/badge.svg
   :target: https://github.com/metawards/MetaWards/actions?query=workflow%3ABuild
   :alt: Build Status

.. image:: https://badge.fury.io/py/metawards.svg
   :target: https://pypi.python.org/pypi/metawards
   :alt: PyPi version

.. image:: https://pepy.tech/badge/metawards
   :target: https://pepy.tech/project/metawards
   :alt: Number of downloads

.. image:: https://camo.githubusercontent.com/1f42db8f4826e53d3ffb6c3d298ced7b0ec73a2547d44c023882d47dd021b95d/68747470733a2f2f6a6f73732e7468656f6a2e6f72672f7061706572732f31302e32313130352f6a6f73732e30333931342f7374617475732e737667
   :target: https://doi.org/10.21105/joss.03914
   :alt: JOSS paper

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5562737.svg
   :target: https://doi.org/10.5281/zenodo.5562737

Please take a look at the :doc:`features <features>` to see what
MetaWards can do. Follow the :doc:`quick start guide <quickstart/index>` to
see how to quickly get up and running using MetaWards to model your own custom
disease or metapopulation model.

Citation
========

Please cite use of this software following the instructions in the 
"Cite this repository" link in the `GitHub repository <https://github.com/metawards/MetaWards>`__,
and also please cite our [Journal of Open Source Software paper](https://doi.org/10.21105/joss.03914).


Scientific Background
=====================

MetaWards implements a stochastic metapopulation model of disease
transmission. It can scale from modelling local transmission up to
full national- or international-scale metapopulation models.

This is a Python port of the
`MetaWards <https://github.com/ldanon/MetaWards>`__ package originally written
by Leon Danon. This port has been performed with Leon's support by the
`Bristol Research Software Engineering Group
<https://www.bristol.ac.uk/acrc/research-software-engineering/>`__.

It is was originally developed to support modelling of disease transmission
in Great Britain. The complete model description and the
original C code are described here;

*  *"The role of routine versus random movements on the spread of disease
   in Great Britain"*, Leon Danon, Thomas House, Matt J. Keeling,
   Epidemics, December 2009, 1 (4), 250-258; DOI:
   `10.1016/j.epidem.2009.11.002 <https://doi.org/10.1016/j.epidem.2009.11.002>`__

*  *"Individual identity and movement networks for disease metapopulations"*,
   Matt J. Keeling, Leon Danon, Matthew C. Vernon, Thomas A.
   House Proceedings of the National Academy of Sciences,
   May 2010, 107 (19) 8866-8870; DOI:
   `10.1073/pnas.1000416107 <https://doi.org/10.1073/pnas.1000416107>`__

In this model, the population is divided into electoral wards. Disease
transmission between wards occurs via the daily movement of individuals.
For each ward, individuals contribute to the *force of infection* (FOI)
in their *home* ward during the night, and their *work* ward during the
day.

This model was recently adapted to model CoVID-19 transmission in
England and Wales, with result of the original C code
published here;

* *"A spatial model of CoVID-19 transmission in England and Wales:
  early spread and peak timing"*, Leon Danon, Ellen Brooks-Pollock,
  Mick Bailey, Matt J Keeling, Philosophical Transactions of the Royal Society B, 376(1829); DOI:
  `10.1098/rstb.2020.0272 <https://doi.org/10.1098/rstb.2020.0272>`__

This Python code is a port which can identically reproduce the outputs
from the original C code as used in that work. This Python code has
been optimised and parallelised, with additional testing added to ensure
that development and scale-up of MetaWards has been robustly and
efficiently conducted.

Features
========

.. toctree::
   :maxdepth: 2

   features

Installation
============

.. toctree::
   :maxdepth: 2

   install

Model Data
==========

.. toctree::
   :maxdepth: 2

   model_data

Quick Start Guide
=================

.. toctree::
   :maxdepth: 2

   quickstart/index

Tutorial
========

.. toctree::
   :maxdepth: 2

   tutorial/index

Files
=====

.. toctree::
   :maxdepth: 2

   fileformats/index

Usage
=====

.. toctree::
   :maxdepth: 2

   usage
   cluster_usage

Getting help
============

.. toctree::
   :maxdepth: 2

   support

Contributing
============

.. toctree::
   :maxdepth: 2

   contributing
   devsupport
   roadmap
   packaging
   development
   snaglist

Documentation
=============

.. toctree::
   :maxdepth: 2

   api/index

Changelog
=========
.. toctree::
   :maxdepth: 2

   changelog

Acknowledgements
================
.. toctree::
   :maxdepth: 2

   acknowledgements

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
=========
Changelog
=========

`1.6.2 <https://github.com/metawards/MetaWards/compare/1.6.1...1.6.2>`__ - March 18th 2022

* Added links to our published JOSS paper!
* Updated links in the website
* Fixed a crash in the progress bar code if the update frequency is too small

`1.6.1 <https://github.com/metawards/MetaWards/compare/1.6.0...1.6.1>`__ - October 11th 2021
--------------------------------------------------------------------------------------------

* Added CITATION.cff file
* Fixed a bug where custom parameters were not being copied into runs correctly
* Worked on better packaging of the metawards_random library on Windows
* Upgrade to the latest version of Sphinx

`1.6.0 <https://github.com/metawards/MetaWards/compare/1.5.1...1.6.0>`__ - May 25th 2021
--------------------------------------------------------------------------------------------

* Added work_to_play iterators, that allow the weekend to be modelled
  by treating workers as players. Updated
  :doc:`the tutorial <tutorial/part03/02_weekend>` to include
  this better model.
* Added Apple M1 support. MetaWards now compiles and runs on Apple M1
  processors. OpenMP support is included, so it can also run
  in parallel. Setup scripts auto-detect if an M1 version of Python
  is used, and will compile MetaWards appropriately.
* Fixed a bug in :func:`metawards.run` whereby the MetaWards executable
  was not found when running on Windows. This stopped MetaWards from
  being run from Jupyter or RStudio on Windows. This should now
  work correctly, including with proper handling of error output
  (i.e. fixed the hanging reported in issue 169),
  catching of failures, and passing of command line arguments.
  To aid running, the API was extended with the following functions;
  :func:`metawards.find_mw_exe`, :func:`metawards.find_mw_include`
  and :func:`metawards.get_reticulate_command`.
* Added the ability to seed infections in the worker population. By default,
  infections will be seeded in the player population. However, if there
  are any remaining infections to seed (or no players in a ward) then
  the infections will be seeded in the workers.
* Added support for compiling Cython plugins that link to MetaWards
  C functions. This means that you can now write plugins that are compiled
  dynamically that use the random number generators in metawards_random.
  This is demonstrated in `tests/test_cython_iterator.py`, specifically
  in the iterator `tests/iterators/cython_iterator.pyx`.
* Fixed a bug in :func:`~metawards.VariableSet` whereby array variables
  could not be parsed if they were on the first line of the file
  (the commas confused the parser and put it into column mode).
  Added a check so that arrays can now be detected on the first line,
  and the parser put into row mode.
* Added more docs to the installation instructions to let new users
  know that they have to install MetaWardsData.

`1.5.1 <https://github.com/metawards/MetaWards/compare/1.5.0...1.5.1>`__ - February 26th 2021
--------------------------------------------------------------------------------------------

* Fixed a bug in the new :func:`~metawards.mixers.merge_matrix_multi_population`
  function that caused it to calculate incorrect values when one or more
  of the demographics had zero members (it didn't catch this divide by
  zero case). This led to too-high FOIs which produced meaningless results.
* Updated the setup.cfg to correctly state that MetaWards supports
  Python 3.9
* Updated sphinx to the latest version

`1.5.0 <https://github.com/metawards/MetaWards/compare/1.4.1...1.5.0>`__ - February 9th 2021
--------------------------------------------------------------------------------------------

* Added :func:`~metawards.mixers.merge_matrix_single_population` and
  :func:`~metawards.mixers.merge_matrix_multi_population` to account
  for different normalisation factors when merging FOIs calculated
  across multiple demographics. This is
  :doc:`described in the tutorial here <tutorial/part05/04_contacts>`.
  It allows the population dynamics to be better modelled, e.g.
  either as all demographics being part of the same population
  with equal contact probabilities, or the demographics having
  different contact probabilities, e.g. as in traditional
  age-based demographics.
* Added :func:`~metawards.mixers.mix_evenly_single_population`
  and :func:`~metawards.mixers.mix_none_single_population` as
  convenience mixers for the case where all demographics are
  part of the same population.
* Added :func:`~metawards.mixers.mix_evenly_multi_population`
  and :func:`~metawards.mixers.mix_none_multi_population` as
  convenience mixers for the case where the contacts between
  demographics depend on the number of individuals in each
  demographic.
* Fixed a small bug where MetaWards raised an exception
  in multi-demographic runs when a custom disease model
  doesn't have a `I` or `R` stage.

`1.4.1 <https://github.com/metawards/MetaWards/compare/1.4.0...1.4.1>`__ - November 18th 2020
---------------------------------------------------------------------------------------------

* Fixed a small bug in `output_core` where the Workspace of a demographic
  subnetwork (its subspace) was being zeroed incorrectly after statistics
  had been accumulated and processed, when the code was run on more than
  4 processors. This prevented extractors from getting this data for
  demographics. Simple fix :-). Bug does not impact any previous runs
  or would have been likely to have caused any issues.
* Added Python 3.9 to the GitHub Actions build matrix, and tested that
  the Python 3.9 packages are built and run correctly.
* Bumped Sphinx up to the latest version.


`1.4.0 <https://github.com/metawards/MetaWards/compare/1.3.0...1.4.0>`__ - August 14th 2020
-------------------------------------------------------------------------------------------

* MetaWards now includes ward-local parameters (e.g. cutoff and scale_uv), plus
  supports custom user ward-local parameters. This supports modelling of
  different ward-local behaviour, e.g. local control measures or
  local lockdowns. Examples of a local lockdown model is
  in :doc:`chapter 7 <tutorial/part03/07_cutoff>` and
  :doc:`chapter 8 <tutorial/part03/08_local_lockdown>` of
  :doc:`part 3 <tutorial/index_part03>` of the tutorial, plus
  :doc:`chapter 5 of part 8 <tutorial/part08/05_local_network>`.
  This can also be used to model
  :doc:`ward-local vaccination strategies <tutorial/part09/02_vaccinate>`.
* We have re-worked the go_functions. We have added
  :func:`~metawards.movers.go_ward` that can be used with
  :class:`~metawards.movers.MoveGenerator` to specify moves between
  any and all combinations of demographics, disease stages and wards
  (for players) and ward links (for workers). Individuals can move
  from worker to player, player to worker, move around the network,
  move to different networks in different demographics, move
  between different disease stages, move from susceptible to
  infected or infected to susceptible etc. etc. This enables
  some advanced modelling, e.g. of vaccinated or recovered individuals
  who gradually lose immunity and become susceptible again, movements
  associated with, e.g. start of university, and movements and
  quarantine associated with holidays.
* Added a :meth:`Ward.bg_foi <metawards.Ward.bg_foi>` per-ward
  parameter to set a background force of infection that can be used
  to drive (or suppress) infections regardless of the number of
  infecteds in each ward.
* Added a global :meth:`Parameters.bg_foi <metawards.Parameters.bg_foi>`
  to set a per-network background force of infection. This is
  useful as a way to model
  :doc:`holiday destinations as different networks <tutorial/part09/05_holidays>`.
* Added a global :meth:`Parameters.scale_uv <metawards.Parameters.scale_uv`
  make it easier to set and control the scale_uv parameter via
  a design file, parameter or via an adjustment in a demographic.
* Finalised demographic adjustment support, enabling you to create
  demographics with adjusted parameters, e.g. see the
  :doc:`holiday destinations example <tutorial/part09/05_holidays>`.
* Added :class:`~metawards.movers.MoveRecord` that can be used to
  record all moves performed by a mover or go_function. Added
  :func:`~metawards.go_record` that can move specific individuals
  according to specific moves indicated in the
  :class:`~metawards.movers.MoveRecord`. As there is a
  :func:`MoveRecord.invert <metawards.movers.MoveRecord.invert>` you can
  use this to reverse moves.
* Added :meth:`is_infected <metawards.Disease.is_infected>` parameter
  to :class:`~metawards.Disease` to mark whether a disease stage is classed
  as being an infected stage. This is useful for non-recovered,
  non-susceptible and non-infected stages, e.g. vaccinated (V) stages.
  For example see :doc:`this tutorial <tutorial/part09/02_vaccinate>`.
* MetaWards now has a `proper R package <https://github.com/metawards/rpkg>`_.
  You can now install and update
  MetaWards directly from within R. See the updated
  :doc:`installation instructions <install>` and the
  :doc:`R quickstart guide <quickstart/01_R>`.
* Added a ``--UV-max`` command line parameter so that you can specify
  the date in the year when disease transmission is highest (if UV is not
  equal to 1.0, and thus disease transmission is seasonal). This defaults
  to the first day of the outbreak.
* Optimised :func:`~metawards.iterators.advance_foi` to skip calculations
  of FOI for a stage if beta[stage] is zero. This changes the order
  of random numbers, so meaning that this version of metawards will
  give different output than older versions for the same input and
  same random number seed. We've made a similar change to the original
  C code to make sure that this has not invalidated the results.
* Added a "null" or "scratch" ward that can be used to temporarily
  store individuals during a day. This is useful when implementing more
  complex moves that involve gathering and scattering populations.
* Removed all parameters and dead code that were ported from the original
  C code but are unused.

`1.3.0 <https://github.com/metawards/MetaWards/compare/1.2.0...1.3.0>`__ - July 22nd 2020
-----------------------------------------------------------------------------------------

* Added a new :doc:`quick start guide <quickstart/index>` that quickly
  showcases the main features of MetaWards. A Python, R and command line
  version is available, so this should suit a range of audiences.
* Added support for different demographics to use different networks.
  This is described partially in the :doc:`tutorial/index_part08`,
  but mostly in the new :doc:`quick start guide <quickstart/index>`.
  This will be documented further in the tutorial in a future release
  (e.g. 1.3.1 or 1.4.0).
* Added a :func:`metawards.run` function to run MetaWards jobs from the API.
  This enables jobs to be run from within Python or R scripts, or to run
  interactively from within, e.g. RStudio or Jupyter.
* Added in R support via reticulate. You can now use the MetaWards API
  within R, plus, via the new :func:`metawards.run` function you can
  write nice tutorials or vignettes that include running the jobs.
  Aim to create a CRAN MetaWards package in a future release.
* Cleaned up the Python API so that this is as flexible as the R API.
  Made sure that key classes, like :class:`~metawards.Disease`,
  :class:`~metawards.InputFiles` and :class:`~metawards.Demographics`
  are easy to use and can serialised to/from JSON.
* New :class:`~metawards.Ward` / :class:`~metawards.Wards` API to let
  you easily create new networks in Python or R.
  You can convert :class:`~metawards.Network` to and from a
  :class:`~metawards.Wards`, and these can be saved and loaded from JSON.
  You can harmonise multiple Wards objects, which enables different
  demographics to use different networks. Also can now refer to wards
  in a network by name rather than index.
* Fixed issues with the "single" ward model. This did not assign any
  player weights, so outbreaks were incorrect. This is now fixed, and the
  single-ward model now matches a manually-created single ward model.
* Added convenience executables (metawards-python, metawards-jupyter
  and metawards-reticulate) to make it easier for users to use the
  right Python executable if many are installed on the system.
* Cleaned up the output and changed "UV" to "scale_uv" as this clashed with
  the UV command-line parameter (and confused people).
* Fixed a bug where the "population" parameter was ignored for repeated
  single-ward network runs.
* More robust reading of the traditional network file format
* Added progress bars for slow operations :-)
* Better support for sequential naming of output directories for repeated runs
* "master" branch was renamed to "main"

`1.2.0 <https://github.com/metawards/MetaWards/compare/1.1.0...1.2.0>`__ - June 26th 2020
-----------------------------------------------------------------------------------------

* Added the ability to use custom-named disease stages. You can now run any
  type of model, and are not limited to ``S``, ``E``, ``I`` and ``R``.
  Learn more in the :doc:`tutorial here <tutorial/part07/05_named_stages>`.
* Improved formatting out information output to the user regarding different
  disease stages. This includes better console output and also more
  informative output data files. Again, this is all detailed in the
  above tutorial.
* Updated all output files to support the summary data for custom
  named disease stages. Now you can collect the data you want directly
  without needing to build a custom extractor - just say which mapping
  stage you want. Again, this is described in the above tutorial.
* Added really flexible support for reading in different formats of
  additional seeds. See the :doc:`tutorial here <tutorial/part08/01_networks>`
  and the new :doc:`fileformats documentation <fileformats/index>`.
  This includes being able to read extra seeds from the command line,
  rather than needing to always write a file.
* Added in the ability to seed infections by date as well as day. Also
  seeding wards by name as well as index (e.g. ``Clifton/Bristol``).
* Added in :class:`metawards.Interpret` to consolidate all of the code
  used to interpret strings into data types. This increases the power
  and flexibility of the data parsers, and adds in new features such
  as reading in random data, or adding math functions to the
  expression support, e.g. ``pi * sqrt(3.5)`` now works.
* Added cython support for plugins. If your plugin ends with ``.pyx`` and
  you have cython installed, then it will be compiled at run time.
  This should enable you to write plugin that are both powerful and fast.
* Fixed a deadlock on Linux when using multiprocessing and OpenMP together
* Removed the unused ``.err`` file.
* Removed ``TotalInfections.dat.bz2`` file (and similar) as these were
  difficult to work with and not well understood. Replaced with
  ``total_infections.csv.bz2`` (and similar) files, which have more
  information and are easier to work with (e.g. have column names).

`1.1.0 <https://github.com/metawards/MetaWards/compare/1.0.0...1.1.0>`__ - June 11th 2020
-----------------------------------------------------------------------------------------

* Different demographics can now follow different disease pathways. This
  supports modelling of super-spreaders and hospitals, as described
  in :doc:`part 7 of the tutorial <tutorial/index_part07>`.
* Variables in demographic sub-networks can be scanned independently from
  the overal network or other sub-networks. This means you can, e.g.
  enact lock-downs in specific demographics, or scan disease parameters
  for different demographics.
* Added a :meth:`~metawards.movers.go_stage` function that moves individuals
  from and to specific disease stages in different demographics. This is
  used to support conditional branching, e.g. 20% of I2 infecteds go to
  hospital.
* Added "--star-as-E", "--star-as-R" and "--disable-star" command line
  arguments to control how the "*" state is counted in the summary outputs.
  This enables it to be counted as an extra "E" state, which makes the
  output more meaningful and more easily interpretable.
* Clarified the meaning the "day 0" and "day 1". Now "day 0" is before
  the model run starts (i.e. setup). The first iteration of the model
  run is "day 1". This is a change from previous versions, which called
  the first half of the first iteration "day 0" and the second half "day 1".
  Since seeding happens in the first half, this means that we now seed one
  day earlier than previous versions, so outbreaks are now one day ahead.
* Fixed a major bug in calculation of the demographic sub-networks
  denominators. These have not been used in production yet. If you
  are going to use demographic sub-networks then please make sure
  you use this version (1.1.0) or above.
* Added database support to :class:`~metawards.OutputFiles`, so that you
  can now write data to SQLite3 databases. This is described in a new
  part of :doc:`tutorial chapter 4 <tutorial/part04/04_rates>`.
* Added in extra output to :class:`~metawards.Workspace` so that you can
  get the populations of all disease stages for all demographics. This
  is demonstrated in a rate calculation, also in the
  :doc:`new tutorial chapter 4 <tutorial/part04/04_rates>`.
* Fixed a directory permissions bug that appeared sometimes on windows.
* Fixed an existing bug from the C code whereby user-set values of
  contrib_foi are ignored. This had no impact as these values are always 1.0.
* Fixed a bug in distribute_remainders that meant that individuals could
  sometimes still be added to a demographic even if the desired percentage
  was zero.

`1.0.0 <https://github.com/metawards/MetaWards/compare/0.12.0...1.0.0>`__ - May 23rd 2020
-----------------------------------------------------------------------------------------

* Improved "go_to" and "go_isolate" functions, which now support modelling
  self-isolation and quarantine. This is all demonstrated in a new
  part 6 of the tutorial.
* Added an InteractionMatrix class to make it easier to create more
  sophisticated interaction matricies.
* Added ability for any plugin to signal that the model run should end
  after the current iteration by raising a StopIteration exception
* Added a "--model single" mode that uses a single-ward model for
  debugging and validation purposes.
* Updated parallel runners (multiprocessing, scoop and MPI) to return
  results as they are available, so that the Console can report summaries
  and live progress.
* Added a developer's "debug" mode to the Console, complete with nice
  variable printing.
* Lots of file and text encoding fixes, particularly to fix unicode
  issues on windows.
* Finally fixed the issue on windows where the wrong plugin would
  sometimes be loaded.
* Updated all tutorial outputs to the new format.
* Fixed a runtime check exception that occurred on rare occasions on Windows.
  This didn't cause any errors in data, but did stop runs from continuing
  when the run-time test was failed.


`0.12.0 <https://github.com/metawards/MetaWards/compare/0.11.2...0.12.0>`__ - May 18th 2020
--------------------------------------------------------------------------------------------

* Switched to configargparse to have better management of command line options,
  plus adding the ability to set options using a config file. This is now
  written to the output directory of each job to support reproducibility.
* metawards-plot defaults to png output if pillow (and jpeg) are not available
* Got basic movers working and added half of the sixth part of the tutorial,
  where self-isolation is modelled.
* Added rich-console support, which has significantly altered the look and
  feel of metawards. Output is now more robust, with more info given in
  real time for parallel jobs, plus all output now also being recorded
  to output/console.txt.bz2, so that no output is lost.
* Added theming support and a "simple" theme activated using "--theme simple"
  for those that don't like colour ;-)
* Added support for setting the number of repeats for a VariableSet into
  the output file. Also can specify different number of repeats for different
  adjustable variable sets on the command line.
* Cleaned up the design file and user custom variable file parsing to use
  csv and support a wide range of formats, variable types and inputs.
  Can now directly work with dates, ints, floats, bools and strings. This
  is intelligent, and will use the best type it thinks, but it can be
  forced by the user via a d"3.4" numpy-type syntax
* Improved the robustness of the parallel runners (multiprocessing, scoop
  and mpi4py) such that errors in one job don't break all jobs. These are
  now handled individually and recorded properly. Jobs are run async so
  that results are processed and feedback is given to the user as soon
  as it is available.
* Updated all of the tutorial to use lurgy3 - accidentally had gone back
  to lurgy2 in part 5.

`0.11.2 <https://github.com/metawards/MetaWards/compare/0.11.1...0.11.2>`__ - May 11th 2020
--------------------------------------------------------------------------------------------

* Minor bugfixes
* Use last matching custom function rather than first, so
  that the examples in the tutorial work and behaviour is more natural
* Caching network builds so that they are more thoroughly tested, fixed
  bug in networks.copy that meant that independent copies weren't made.
  This bug did not impact any past results or runs.
* Added more validation tests of the mixers
* Cleaned up website typos and fixed the version switcher
* Fixed packaging problems that caused broken builds when pip installing
  from a .tgz sdist package.

`0.11.1 <https://github.com/metawards/MetaWards/compare/0.11.0...0.11.1>`__ - May 10th 2020
--------------------------------------------------------------------------------------------

* Fixed CI/CD to produce working sdist and bdist packages

`0.11.0 <https://github.com/metawards/MetaWards/compare/0.10.0...0.11.0>`__ - May 10th 2020
--------------------------------------------------------------------------------------------

* Code now fully works and has been tested on Windows :-)
* Major update of the API to support a Networks of multiple Network objects
* This has been used to support modelling multiple demographics
* Added in movers and mixers to enable a user to customise how individuals
  are moved between demographics and how the FOIs of demographics are
  merged together (e.g. via an interaction matrix). This is demonstrated
  in part 5 of the tutorial which shows how this can be used to model
  shielding
* Allow compilation using compilers that don't support OpenMP - now compiles
  even on stock OS X.
* Added more extractors and can now output files that are needed for graphics
* Added a special random number seed to support debugging
* Moved random number files to a separate library which is now properly
  compiled and linked.
* Updated CI to CI/CD and now build the OS X, Windows and ManyLinux wheels
* Updated URLs to metawards.org
* Allow multiple multi-node jobs to run from a single directory (they now
  have their own hostfiles)
* Updated metawards-plot to render multi-demographic trajectories and
  to make better animations.
* General bug fixes and speed-ups :-)

`0.10.0 <https://github.com/metawards/MetaWards/compare/0.9.0...0.10.0>`__ - April 27th 2020
--------------------------------------------------------------------------------------------

* Created all of the extract framework to support customising the output
  and analysis during a run.
* Created a better Workspace class for holding accumulated data during extract
* Completed most of the extractor tutorial
* Added in WardInfo(s) to get metadata about wards, and to support searching
  for wards via their name, code, authority and region

`0.9.0 <https://github.com/metawards/MetaWards/compare/0.8.4...0.9.0>`__ - April 24th 2020
------------------------------------------------------------------------------------------

* Merged in latest changes from the C code. Now gives complete agreement,
  including via a custom iterator that repeats the lockdown model.
* Support x/y and lat/lon coordinates and distances. Now works properly
  with the 2011UK model data
* Added an example of a lockdown parameter set scan

`0.8.5 <https://github.com/metawards/MetaWards/compare/0.8.3...0.8.5>`__ - April 22nd 2020
------------------------------------------------------------------------------------------

* Small bugfixes to support the loading of the 2011UK model data
* Cleaned up the website and added the version combo box

`0.8.3 <https://github.com/metawards/MetaWards/compare/0.8.0...0.8.3>`__ - April 21st 2020
------------------------------------------------------------------------------------------

* Fixing CI/CD so that I can build and deploy on a new tag (hopefully 0.8.2)

`0.8.0 <https://github.com/metawards/MetaWards/compare/0.7.0...0.8.0>`__ - April 21st 2020
------------------------------------------------------------------------------------------

* Automated github actions for building a versioned website plus automating
  building the packages.
* Switched default for UV parameter to 0.0, as this should not normally be 1.0
* Added custom user variables both for scanning and to act as inputs that
  may be used by custom advance and iterate functions. Detailed tutorial
  now shows how these can be used to model a lockdown.
* Improved speed of custom iterators

`0.7.1 <https://github.com/metawards/MetaWards/compare/0.6.0...0.7.1>`__ - April 17th 2020
------------------------------------------------------------------------------------------

* Small bugfixes to support all of the examples in part 3 of the tutorial

`0.7.0 <https://github.com/metawards/MetaWards/compare/0.6.0...0.7.0>`__ - April 17th 2020
------------------------------------------------------------------------------------------

* Lots of progress with the project website, including a detailed tutorial
* Support fully customisable disease models, and can adjust any disease
  parameter using a more flexible input file format
* Can record the date in a model run, plus set the starting day and date
* Broken up the iterate function into :mod:`metawards.iterators`, and
  can now have the user create their own custom iterators. Tutorial on
  how to do this will appear soon.
* Broken up the extract_data function into :mod:`metawards.extractors`,
  and will soon enable a user to create their own. Tutorial on how
  to do this will appear soon.
* Added metawards-plot to create simple plots and animations. This is
  particularly useful when working through the tutorial.
* General code cleaning, documentation improvements and nice-to-haves
  that make the code easier to use.

`0.6.0 <https://github.com/metawards/MetaWards/compare/0.5.0...0.6.0>`__ - April 9th 2020
-----------------------------------------------------------------------------------------

* Wrote an initial draft of the complete project website
* Fixed packaging problems that prevented installation of older packages
  on some systems

`0.5.0 <https://github.com/metawards/MetaWards/compare/0.4.0...0.5.0>`__ - April 8th 2020
-----------------------------------------------------------------------------------------

* Support running multiple model runs in serial or in parallel
* Support aggregation and writing of model multiple model run outputs
  to the same directory, including to a single shared CSV data file.
* Support for parallel running via multiprocessing, mpi4py or scoop

`0.4.0 <https://github.com/metawards/MetaWards/compare/0.3.1...0.4.0>`__ - April 7th 2020
-----------------------------------------------------------------------------------------

* Parallelisation of individual model runs using OpenMP
* Parallel code scales to large numbers of cores and can complete individual
  runs in 10-15 seconds.

`0.3.1 <https://github.com/metawards/MetaWards/compare/0.3.0...0.3.1>`__ - April 5th 2020
-----------------------------------------------------------------------------------------

* Minor bug fixes in packaging and misplaced commits caused by move of
  repository

`0.3.0 <https://github.com/metawards/MetaWards/compare/v0.2.0...0.3.0>`__ - April 5th 2020
------------------------------------------------------------------------------------------

* Adding in a simple profiler to support optimisation of the code
* Replaced GSL random number generator with a more liberally licensed and
  easily bundled generator extracted from numpy.
* Switched code to the https://github.com/metawards organisation
* Optimised more using cython and raw C for file reading
* Added automatic versioning of packages and files using versioneer
* Cleaned up the repository and added status badges

`0.2.0 <https://github.com/metawards/MetaWards/compare/v0.1.0...v0.2.0>`__ - March 31st 2020
--------------------------------------------------------------------------------------------

* Cythonizing the bottleneck code to bring the python code up to a comparable
  performance as the original C code.
* Added in packaging information and general repository and file cleaning.

`0.1.0 <https://github.com/metawards/MetaWards/releases/tag/v0.1.0>`__ - March 29th 2020
----------------------------------------------------------------------------------------

* Fully working Python port of the original C code that completely reproduces
  the results of the C code when given the same random number seed. However,
  it is *significantly* slower! Python port has promise, so worth exploring
  different options for speeding the code up.

`Start of the Python port <https://github.com/metawards/MetaWards/commit/ef989ece450c40fe0ddb9f22e21693c90afb432e>`__ - March 25th 2020
---------------------------------------------------------------------------------------------------------------------------------------

* Imported code from https://github.com/ldanon/metawards and began thinking
  about what the code was and trying to understand it. Decided to write
  a port as I find that if I can translate something, then I can
  understand it.
================
Acknowledgements
================

We gratefully acknowledge funding from the EPSRC (EP/N018591/1) who supported
the research software engineering time needed to create this version
of MetaWards.

We also gratefully acknowledge the many helpers and contributers to this
project.
=========================
Installation instructions
=========================

MetaWards is a Python package, but it does come with R wrappers. This
means that you can choose to install MetaWards either via Python or
via R. Note that MetaWards depends on data files that are
held in MetaWardsData. To use MetaWards, you must install MetaWardsData
following :doc:`these instructions <model_data>`.

Installation via Python
=======================

MetaWards is compiled and
`tested on Windows, Linux and OS X <https://github.com/metawards/MetaWards/actions>`__,
but should be able to compile and install on any operating system
with a working Python >= 3.7 and a C compiler.

To install MetaWards, you first need to install Python >= 3.7. To check
if you have Python 3.7 installed type;

.. code-block:: bash

    python -V

This should print out the version number of the Python you have available.
If the version number is ``2.???`` then you have a very old Python 2. If
this is the case, try;

.. code-block:: bash

    python3 -V

and see if you have a Python 3 that has a version number >= 3.7. If so,
please use ``python3`` instead of ``python``.

If you don't have Python >=3.7 installed, then you can install Python
either via your package manager, or more simply, by using
`anaconda <https://anaconda.org>`__.

Once you have a working Python >= 3.7, the easiest way to install
MetaWards is using
`pip <https://pip.pypa.io/en/stable/>`__.

.. code-block:: bash

    pip install metawards

(if this doesn't work, then you may need to use the command ``pip3``,
or you may have to `install pip <https://pip.pypa.io/en/stable/installing/>`__.
If you have trouble installing pip then we recommend that you download
and install `anaconda <https://anaconda.org>`__, which has pip included)

To install a specific version, e.g. 1.5.1, type

.. code-block:: bash

    pip install metawards==1.5.1

This will install a binary version of metawards if it is available for your
operating system / processor / version of python. If not, then
the metawards pyx code that was compiled into C using
`cython <https://cython.org>`__,
will be compiled into an executable using your system C compiler
(e.g. `gcc <https://gcc.gnu.org>`__ or `clang <https://clang.llvm.org>`__).
You can control which C compiler to use by setting the ``CC`` environment
variable, e.g.

.. code-block:: bash

    CC=gcc pip install metawards

MetaWards is written in standard C, using
`OpenMP <https://www.openmp.org>`__ for parallelisation,
and it has no external dependencies, so
it should compile without issue. Note that the system C compiler on
OS X (Mac) does not support OpenMP. In this case you the code will
compile with OpenMP disabled. If you want to use all of the cores
of your computer than you will need to install
a different compiler, e.g. installing clang via
`homebrew <https://brew.sh>`__. If you have any problems then please
`post an issue on GitHub <https://github.com/metawards/MetaWards/issues>`__.

Once installed, you can run `metawards` by typing

.. code-block:: bash

    metawards --version

This should print out some version information about MetaWards,
including whether or not it found MetaWardsData.

If this doesn't work, then it is possible that the directory into which
`metawards` has been placed is not in your PATH (this is quite
common on Windows). You can find the
location of the `metawards` executable by starting Python and
typing;

.. code-block:: python

    > import metawards
    > print(metawards.find_mw_exe())

You can then run `metawards` either by typing out the full path
to the executable that is printed, or by adding the directory
containing the executable to your PATH.

Installation via R
==================

MetaWards can be used within R via the `reticulate <https://rstudio.github.io/reticulate/>`_
package. We have built a MetaWards R package that simplifies this use.

To use this, you must first have installed metawards as above,
via the Python route. This is because the version of Python used
by default in R is too old to run MetaWards (MetaWards needs
Python 3.7 or newer, while the default in reticulate is to install
and use Python 3.6).

Once you have MetaWards (and MetaWardsData)
installed in Python, you first need to
get the reticulate command that you will need to tell R which
Python interpreter to use. We have written a function to do
this for you. Open Python and type;

.. code-block:: python

    import metawards
    print(metawards.get_reticulate_command())

You should see printed something like

.. code-block:: R

    reticulate::use_python("/path/to/python", required=TRUE)

where `/path/to/python` is the full path to the python executable
that you are using.

Next, open R or RStudio and install reticulate (if you haven't
done that already).

.. code-block:: R

   > install.packages("reticulate")

Next, you should install the MetaWards R package from GitHub.
You need to use the `devtools` library to do this.

.. code-block:: R

   > library(devtools)
   > install_github("metawards/rpkg")

Next, you need to tell reticulate to use the Python
executable you found earlier. Copy in the reticulate
command exactly as it was printed by Python, e.g.

.. code-block:: R

   > reticulate::use_python("/path/to/python", required=TRUE)

Next, load the ``metawards`` package and use the ``py_metawards_available``
to check that MetaWards can be found and loaded.

.. code-block:: R

   > metawards::py_metawards_available()
   [1] TRUE

.. note::

   In the future, once reticulate defaults to Python 3.7 or
   above, you will be able to install MetaWards directly
   by calling `metawards::py_install_metawards()`. This
   command will install metawards into the default Python
   that comes with reticulate.

Once installed, you can check if there
are any updates to MetaWards available directly in R, via;

.. code-block:: R

   > metawards::py_metawards_update_available()

and you can update to the lastest version using;

.. code-block:: R

   > metawards::py_update_metawards()

Source install
==============

You can download a source release of MetaWards from the
`project release page <https://github.com/metawards/MetaWards/releases>`__.

Once you have downloaded the file you can unpack it and change into
that directory using;

.. code-block:: bash

   tar -zxvf MetaWards-X.Y.Z.tar.gz
   cd MetaWards-X.Y.Z

where ``X.Y.Z`` is the version you downloaded. For the 1.4.0 release
this would be;

.. code-block:: bash

    tar -zxvf MetaWards-1.4.0.tar.gz
    cd MetaWards-X.Y.Z

Next you need to install the dependencies of MetaWards. Do this by typing;

.. code-block:: bash

    pip install -r requirements.txt

Now you are ready to compile and install MetaWards itself;

.. code-block:: bash

    make
    make install

You can choose the C compiler to use by setting the ``CC`` environment
variable, e.g.

.. code-block:: bash

    CC=clang make
    CC=clang make install

MetaWards is written in standard C, using
`OpenMP <https://www.openmp.org>`__ for parallelisation,
and it has no external dependencies, so
it should compile without issue. Note that the system C compiler on
OS X (Mac) does not support OpenMP. In this case you the code will
compile with OpenMP disabled. If you want to use all of the cores
of your computer than you will need to install
a different compiler, e.g. installing clang via
`homebrew <https://brew.sh>`__. If you have any problems then please
`post an issue on GitHub <https://github.com/metawards/MetaWards/issues>`__.

Just as for the normal Python install you may need to use the
`find_mw_exe()` function to find the full path to the installed
`metawards` executable if this hasn't been placed in your PATH.

For developers
==============

You can clone the MetaWards repository to your computer and install from
there;

.. code-block:: bash

    git clone https://github.com/metawards/MetaWards
    cd MetaWards
    pip install -r requirements-dev.txt

From this point you can compile as if you have downloaded from source.
As a developer you may want to run the tests and create the website.
To do this type;

.. code-block:: bash

    pytest tests
    make doc

There are shortcuts for running the quick or slow tests, e.g.

.. code-block:: bash

   make test
   make quicktest

Note that the tests assume that you have already downloaded the
model data from `MetaWardsData <https://github.com/metawards/MetaWardsData>`__
and configured this as `described here <model_data.html>`__.
=======
Roadmap
=======

MetaWards is now fully featured for everything that the development
team needs. If you would like to suggest features that should be
added to the roadmap then please
`raise a feature request here <https://github.com/metawards/MetaWards/issues/new?assignees=chryswoods&labels=enhancement&template=feature_request.md&title=%5BFEATURE+REQUEST%5D+-+I%27d+like+MetaWards+to>`_.

More details about individual feature branches, which contain the new
features being actively worked on, is available on the
`GitHub issues <https://github.com/metawards/MetaWards/issues?q=is%3Aissue+is%3Aopen+label%3Afeature-branch>`_
page.
=================
Developer Support
=================

``metawards`` comes with a set of functionality that, we hope, will make
it easier for you to develop the code.

Profiling
---------

Use the ``--profile`` command-line option to switch on profiling. This
is very detailed and is free (doesn't affect the run-time). Use this
if you think the code is slow, or to examine the scaling of parallel
sections.

This will print out a tree showing the call path and the amount of time
in each function or sub-section of the function (more on that later ;-)).

For example, running;

.. code-block:: bash

   metawards -d lurgy --profile

gives;

::

  Total time: 3548.695 ms (3548.695 ms)
    \-Network.build: 3548.695 ms (3546.863 ms)
        \-build_function: 3143.628 ms (3143.561 ms)
            \-build_wards_network: 3143.561 ms (3081.800 ms)
                \-read_work_file: 1232.249 ms
                \-fill_in_gaps: 5.017 ms
                \-build_play_matrix: 1752.456 ms (1752.430 ms)
                    \-build_play_matrix: 1752.430 ms (1751.152 ms)
                        \-allocate_memory: 353.156 ms
                        \-read_play_file: 1367.946 ms
                        \-renormalise?: 0.003 ms
                        \-renormalise_loop: 5.227 ms
                        \-fill_in_gaps: 10.288 ms
                        \-read_play_size_file: 14.532 ms
                \-resize_nodes_and_links: 92.078 ms
        \-read_and_validate: 11.067 ms
        \-add_distances: 256.749 ms
        \-add_lookup: 88.264 ms
        \-read_done_file: 3.898 ms
        \-reset_everything: 9.651 ms (9.033 ms)
            \-reset_work: 4.989 ms
            \-reset_play: 4.008 ms
            \-reset_susceptibles: 0.021 ms
            \-reset_params: 0.015 ms
        \-rescale_play_matrix: 11.787 ms
        \-move_from_play_to_work: 21.819 ms

for the timing for building the Network, and

::

  Total time: 75.417 ms (75.417 ms)
    \-timing for day 3: 75.417 ms (74.053 ms)
        \-<function advance_additional at 0x10186c0e0>: 0.014 ms (0.002 ms)
            \-additional_seeds: 0.002 ms
        \-<built-in function advance_foi>: 13.983 ms (13.942 ms)
            \-setup: 0.031 ms
            \-loop_over_classes: 13.911 ms (13.843 ms)
                \-work_0: 3.539 ms
                \-play_0: 0.014 ms
                \-work_1: 3.411 ms
                \-play_1: 0.014 ms
                \-work_2: 3.421 ms
                \-play_2: 0.013 ms
                \-work_3: 3.418 ms
                \-play_3: 0.013 ms
        \-<built-in function advance_recovery>: 6.339 ms (6.311 ms)
            \-recovery: 6.311 ms
        \-<built-in function advance_infprob>: 0.067 ms (0.042 ms)
            \-infprob: 0.042 ms
        \-<built-in function advance_fixed>: 11.027 ms (10.990 ms)
            \-fixed: 10.990 ms
        \-<built-in function advance_play>: 5.347 ms (5.297 ms)
            \-play: 5.297 ms
        \-<built-in function output_core>: 29.488 ms
        \-<function output_basic at 0x10185e8c0>: 0.710 ms
        \-<built-in function output_dispersal>: 0.428 ms
        \-<function output_incidence at 0x10185ef80>: 3.333 ms
        \-<function output_prevalence at 0x101867320>: 3.317 ms

for the timing of each iteration.

.. note::
   The timings will vary for each iteration as more individuals are
   infected. It will also vary depending on the computer and number
   of threads used.

Each level of the tree gives the total time within that tree, e.g.
``timing for day 3: 75.417 ms (74.053 ms)`` is the total time spent
in all of the sub-functions used to calculate ``day 3``. The first
number (``75.417 ms``) is the actual measured time, while the number
in brackets (``74.053 ms``) is the sum of the reported measured times for
the sub-functions. Any large different between the two indicates that
some time is not being reported.

Profiler
--------

You control what is reported using the :class:`metawards.utils.Profiler` class.
You start timing using ``profiler = profiler.start("name of section")`` and you
stop timing using ``profiler = profiler.stop()``. For example using the
iterator ``iterate_profile.py``;

.. code-block:: python

    from metawards.utils import Profiler

    import time

    def advance_profile(profiler: Profiler, **kwargs):
        p = profiler.start("timing advance_profile")

        for i in range(0, 10):
            p = p.start(f"Timing loop iteration {i}")
            time.sleep(0.02)
            p = p.stop()

        p.stop()

    def iterate_profile(**kwargs):
        return [advance_profile]

and then running using;

.. code-block:: bash

   metawards -d lurgy --iterator iterate_profile

gives this in the output;

::

      \-<function advance_profile at 0x1092bc5f0>: 228.184 ms (228.175 ms)
          \-timing advance_profile: 228.175 ms (228.038 ms)
              \-Timing loop iteration 0: 20.354 ms
              \-Timing loop iteration 1: 25.021 ms
              \-Timing loop iteration 2: 23.624 ms
              \-Timing loop iteration 3: 22.029 ms
              \-Timing loop iteration 4: 20.011 ms
              \-Timing loop iteration 5: 25.011 ms
              \-Timing loop iteration 6: 21.608 ms
              \-Timing loop iteration 7: 25.014 ms
              \-Timing loop iteration 8: 21.466 ms
              \-Timing loop iteration 9: 23.900 ms

As you can see, the ``time.sleep(0.02)`` slept for a little-over
20 milliseconds. As with all profiling, it is worth repeating runs
several times and taking an average, especially if you are interested
in plotting parallel scaling of the individual ``metawards`` functions.

.. note::
   A :class:`~metawards.utils.Profiler` is always passed as the ``profiler``
   keyword argument to all of the plugin classes (iterators, extractors etc.)

By default, all plugin functions are timed, hence why in the output
you can see all of the ``advance_functions`` and ``output_functions`` that
were called, in which order, and how long they took. This is actually
the easiest way to debug whether your plugin function has been called - just
run using ``--profile`` and see if your function is listed in the timing.

printf Debugging
----------------

``metawards`` developers are big fans of
`printf debugging <https://tedspence.com/the-art-of-printf-debugging-7d5274d6af44>`__.

All printing in ``metawards`` is handled using the
:class:`metawards.utils.Console` object, which comes with a handy
:meth:`~metawards.utils.Console.debug` function. Use this to print
debug messages, and, optionally, also the values of variables by passing
them as a list to the ``variables`` keyword, e.g. using an iterator
in the file ``iterate_debug.py``;

.. code-block:: python

    from metawards import Population
    from metawards.iterators import iterate_default
    from metawards.utils import Console

    def iterate_debug(population: Population, **kwargs):
        beta_scale = 0.5
        Console.debug("Hello!", variables=[population, beta_scale])

        return iterate_default(population=population, **kwargs)

and running ``metawards`` using the ``--debug`` keyword;

.. code-block:: bash

    metawards -d lurgy --iterator iterate_debug --debug

you will see;

::

    [12:33:05]                           Hello!                            iterate_debug.py:7

            Name │ Value
     ════════════╪═══════════════════════════════════════════════════════════════════════════
      population │ 2020-05-20: DAY: 0 S: 56082077    E: 0    I: 0    R: 0    IW: 0   UV: 1.0
                 │ TOTAL POPULATION 56082077
      beta_scale │ 0.5


printed for every iteration.

.. note::

   The time of the printout is on the top-left, and the filename and line
   for the debug statement is on the top-right. The time is only printed
   once for each second, so if you have a lot of debug statements printed
   in a single second then they won't show the time.

The debug output is only printed when you run ``metawards`` using the
``--debug`` command line argument, so it is safe to leave in
production code.

Debugging levels
----------------

You can optionally set the debugging level of your output using
the ``level`` keyword argument, e.g.

.. code-block:: python

    from metawards import Population
    from metawards.iterators import iterate_default
    from metawards.utils import Console

    def iterate_debug(population: Population, **kwargs):
        beta_scale = 0.5
        Console.debug("Hello!", variables=[population, beta_scale],
                      level=3)

        return iterate_default(population=population, **kwargs)

would only print out if the debugging level was ``3`` or above.

You can set the level using the ``--debug-level`` command line argument,
e.g.

.. code-block:: bash

   metawards -d lurgy --iterator iterate_debug --debug --debug-level 5

would set the debug level to ``5``, which is above ``3``,
and so the debug output will be printed;

::

    [12:53:02]                      Level 3: Hello!                       iterate_debug.py:10

            Name │ Value
     ════════════╪═══════════════════════════════════════════════════════════════════════════
      population │ 2020-05-22: DAY: 2 S: 56082077    E: 0    I: 0    R: 0    IW: 0   UV: 1.0
                 │ TOTAL POPULATION 56082077
      beta_scale │ 0.5

If the level was below ``3``, or was not set, then this debug output would
not be printed.

Debugging lambdas
-----------------

It is better to avoid constructing expensive debug strings if they are
not going to be printed to the screen. There are two ways to avoid this;

1. Use the :meth:`~metawards.utils.Console.debugging_enabled` function of
   :class:`~metawards.utils.Console` to see if debugging is enabled for
   the level you wish, and only call
   :meth:`Console.debug <metawards.utils.Console.debug>` if it is.

2. Put your debug string into a lambda function. The
   :meth:`Console.debug <metawards.utils.Console.debug>`
   function will call this lambda function only if the debug output
   is enabled. This is really easy to do, e.g. here we change the
   debug statement to;

.. code-block:: python

    from metawards import Population
    from metawards.iterators import iterate_default
    from metawards.utils import Console

    def iterate_debug(population: Population, **kwargs):
        beta_scale = 0.5
        Console.debug(lambda: "Hello!", variables=[population, beta_scale],
                      level=3)

        return iterate_default(population=population, **kwargs)

such that we use ``lambda: "Hello!"`` instead of
``"Hello!"``. This converts the string into a lambda function, meaning
that it should not be generated unless the debug statement is
actually printed.

Testing
-------

``metawards`` has a
`large test suite <https://github.com/metawards/MetaWards/tree/devel/tests>`__
built using `pytest <https://docs.pytest.org/en/latest/>`__.
We encourage you to
look through the tests and use these to help learn how to use the classes.
We also encourage you to write your own tests for your new code :-)
==========
Model data
==========

Getting the data
================

MetaWards uses a large amount of data to describe the network over which
the model is run. This data is stored in a separate
`MetaWardsData <https://github.com/metawards/MetaWardsData>`__ repository.

You can download this data in one of two ways;

Direct download
    Download one of the data releases from the
    `releases page <https://github.com/metawards/MetaWardsData/releases>`__.
    It is best to choose the latest release unless you are trying to
    reproduce an earler calculation.

    You can download the data as a ``zip`` or a ``tar.gz`` archive, and
    can unpack it using either;

    .. code-block:: bash

        tar -zxvf MetaWardsData-0.2.0.tar.gz

    if you downloaded the 0.2.0 release tar.gz file, or

    .. code-block:: bash

        unzip MetaWardsData-0.2.0.zip

    if you downloaded the zip file.

Git clone
    You can clone the data repository using

    .. code-block:: bash

        git clone https://github.com/metawards/MetaWardsData

Finding the data
================

Once you have downloaded the data you will need to tell metawards
where the data is located. This can be done in one of three ways;

1. Set the ``METAWARDSDATA`` environment variable equal to the full path
   to the data, e.g.

   .. code-block:: bash

        export METAWARDSDATA=$HOME/MetaWardsData

   (assuming you have placed it in your home directory)

2. Pass the full path to the data to metawards using the command line argument
   ``--repository``, e.g.

   .. code-block:: bash

        metawards --repository $HOME/MetaWardsData

    (again assumes you have placed it in your home directory)

3. Move the data to ``$HOME/GitHub/MetaWardsData`` as this is the default
   search location, e.g.

   .. code-block:: bash

        mkdir $HOME/GitHub
        mv $HOME/MetaWardsData $HOME/GitHub/

Versioning the data
===================

Provenance is very important for MetaWards. It is important to know what
inputs were used for a model run so that it is possible to recreate the
calculation and thus reproduce the result.

To aid this, MetaWardsData has a small script that you should run after
you have downloaded the data. This will version the data by finding the
git tag of the data. Run this script using;

.. code-block:: bash

    cd $METAWARDSDATA
    ./version

This will create a small file called ``version.txt`` which will look
something like this;

.. code-block:: bash

    {"version": "0.2.0","repository": "https://github.com/metawards/MetaWardsData","branch": "main"}

This file will be read by metawards during a run to embed enough information
in the outputs to enable others to download the same MetaWardsData input
as you.

Understanding the data
======================

The MetaWardData repository contains four types of data;

model_data
    The ``model_data`` directory contains the data needed to construct
    a network of wards and links. Currently the 2011 model is included,
    in the ``2011Data`` directory. The manifest of files that
    comprise this model is ``2011Data/description.json``. This json
    file is read by metawards to find and load all of the files that
    are needed for this model.

diseases
    The ``diseases`` directory contains the parameters for different
    diseases. Each disease is described in its own json file, e.g.
    the parameters for SARS-Cov-2 are in ``diseases/ncov.json``.
    Once of the main purposes of metawards is to perform parameter
    sweeps to adjust these disease parameters so that a model epidemic
    can more closely follow what is observed in reality. From there,
    more meaningful predictions should be able to be made.

parameters
    The ``parameters`` directory contains parameters used to control
    metawards model runs. These represent "good defaults" for the
    values of the parameters, and are stored in the MetaWardsData
    repository to make it easier to version control and track
    provenance of the parameters used for different jobs.
    The parameters used for jobs in March 2020 (and correct as
    of March 29th) are in ``parameters/march29.json``.

extra_seeds
    The ``extra_seeds`` directory contains files that can be used
    to seed new model disease clusters at different wards during
    a model run. Each file can contain as many additional seeds
    as needed. The format is three numbers per line;
    ``t   num   loc``, where ``t`` is the step (day) on which
    the additional infection will be seeded, ``num`` is the number
    of additional infections, and ``loc`` is the index of the
    ward in the model network in which the infection will occur.
=====
Usage
=====

metawards program
=================

MetaWards comes with a command-line program called ``metawards``.
This is installed by default into the same directory as the
``python`` executable used to build the package. To run ``metawards``
type;

.. code-block:: bash

    metawards --version

This prints out the version of ``metawards``, which should look something
like this;

.. program-output:: metawards --version

This version information gives you the provenance of this executable which
should help in reproducing output. It is written to the top of all
metawards outputs.

If you see output like this, with ``WARNING`` lines about the version
not having been committed to git...

::

    ┌──────────────────────────────────────────────────────────────────────┐
    │                                                                      │
    │                                                                      │
    │             MetaWards version 0.12.0+17.gd839bca1.dirty              │
    │                                                                      │
    │                                                                      │
    │                        https://metawards.org                         │
    │                                                                      │
    │                          Source information                          │
    │                                                                      │
    │   • repository: https://github.com/metawards/MetaWards               │
    │   • branch: feature-mover-tutorial                                   │
    │   • revision: d839bca1e32e8c8b1814d7f72667e84ead1a59d7               │
    │   • last modified: 2020-05-20T13:01:14+0100                          │
    │                                                                      │
    │  WARNING: This version has not been committed to git, so you may     │
    │  not be able to recover the original source code that was used to    │
    │  generate this run!                                                  │
    │                                                                      │
    │                      MetaWardsData information                       │
    │                                                                      │
    │   • version: 0.5.0                                                   │
    │   • repository: https://github.com/metawards/MetaWardsData           │
    │   • branch: main                                                   │
    │                                                                      │
    │                        Additional information                        │
    │                                                                      │
    │  Visit https://metawards.org for more information about metawards,   │
    │  its authors and its license                                         │
    │                                                                      │
    └──────────────────────────────────────────────────────────────────────┘

then this means that your ``metawards`` executable has been built using
source code that has not been committed to git, and is therefore not
version controlled. Do not use ``dirty`` software for production jobs
as it will not be possible to recover the software used to produce
outputs at a later date, and thus it may not be possile to reproduce
results.

Getting help
============

The ``metawards`` program has up-to-date and very comprehensive in-built
help for all of its command line options. You can print this help by
typing;

.. code-block:: bash

    metawards --help

The full help is :doc:`available here <metawards_help>`.

.. toctree::
   :hidden:

   metawards_help

Understanding the options
=========================

``metawards`` is a powerful program so it comes with a lot of options.
The most used and thus most important options are given first. These
are;

* ``--disease / -d`` : Specify the name of the disease file to load.
  These files are described in `Model data <model_data.html>`__.
  If the file exists in your path then that will be used. Otherwise the
  file will be
  searched for from the ``MetaWardsData/diseases`` directory. Note that
  you don't need to specify the file type, as this is assumed to be
  ``.json``.

* ``--input / -i`` : Specify the input file adjustable parameters that will be
  explored for the model run. You must supply an input file to be
  able to run a model. This file is described below.

* ``--line / -l`` : Specify the line number (or line numbers) of adjustable
  parameter sets from the input file. Line numbers are counted from 0, and
  multiple line numbers can be given, e.g. ``--line 3, 5, 7-10`` will
  read from lines 3, 5, 7, 8, 9, and 10 (remembering that the first line
  of the file is line 0). If multiple lines are read, then multiple model
  runs will be performed.

* ``--repeats / -r`` : specify the number of times model runs for each
  adjustable parameter set should be repeated. MetaWards model runs
  are stochastic, based on random numbers. The results for multiple
  runs must thus be processed to derive meaning.

* ``--additional / -a`` : specify the file (or files) containing additional
  seeds. These files are described in `Model data <model_data.html>`__.
  You can specify as many or few files as you wish. If the file exists
  in your path then that will be used. Otherwise the file will be
  searched for from the ``MetaWardsData/extra_seeds`` directory.
  Note that you can write the additional seeds directly, rather
  than using a file. To see how, take a look at
  :doc:`this section of the tutorial <tutorial/part08/01_networks>`.

* ``--output / -o`` : specify the location to place all output files. By default
  this will be in a new directory called ``output``. A description of the
  output files is below.

* ``--seed / -s`` : specify the random number seed to use for a run. By default
  a truly random seed will be used. This will be printed into the output,
  so that you can use it together with this option to reproduce a run.
  The same version of ``metawards`` will reproduce the same output when
  given the same input, same random number seed, and run over the
  same number of threads.

* ``--start-date`` : specify the date of *day zero* of the model outbreak.
  If this isn't specified then this defaults to the current date. This
  recognises any date that is understood by
  :func:`metawards.Interpret.date`,
  which includes dates like *today*, *tomorrow*, *Monday*, *Jan 2020* etc.
  The start date is used to trigger events based on day of week or
  date within a model outbreak (e.g. is the day a weekend)

* ``--start-day`` : specify the day of the outbreak, e.g. the default
  is ``0`` for *day zero*. This is useful if you want to start the
  model run from a later day than the start date. Note that the
  start date passed via ``--start-date`` is the date of *day zero*,
  so the first day modelled will be ``--start-day`` days
  after ``--start-date``.

* ``--nthreads`` : specify the number of threads over which to perform a
  model run. The sequence of random numbers drawn in parallel is
  deterministic for a given number of threads, but will be different
  for different numbers of threads. This is why you can only reproduce
  a model run using a combination of the same random number seed and
  same number of threads. The number of threads used for a model run
  is written into the output.

* ``--nprocs`` : specify the number of processes over which you want
  to parallelise the model runs. This is useful if you have multiple
  processors on your computer and you are running multiple model runs.
  Note that this option is set automatically for you if you are
  `running on a cluster <cluster_usage.html>`__.

Using these options, a typical ``metawards`` run can be performed
using;

.. code-block:: bash

    metawards -d ncov -a "1  5  1"

This python port of ``metawards`` was written to reproduce the output
of the original C code. This original code is bundled with this
port and is in the ``original`` directory. There are several integration
tests included in the unit testing suite that validate that the
Python code still reproduces the results generated using the C code.

Understanding the input
=======================

The input file for ``metawards`` is a
:doc:`design file <fileformats/design>` that can be as
simple as a set of lines containing
five comma-separated or space-separated values per line, e.g.

::

  0.95,0.95,0.19,0.91,0.91
  0.90,0.93,0.18,0.92,0.90

or

::

  0.90 0.93 0.18 0.92 0.90

These five values per line adjust the ``beta[2]``, ``beta[3]``,
``progress[1]``, ``progress[2]`` and ``progress[3]`` parameters of
the ``disease`` model as described in `Model Data <model_data.html>`__.

You can optionally choose which parameters will be varied by adding
a title line, e.g.

::

  beta[2]   progress[2]   progress[3]
    0.90        0.92         0.90
    0.85        0.91         0.92

specifies that you want to adjust the ``beta[2]``, ``progress[2]`` and
``progress[3]`` parameters to the specified values.

This file can adjust a lot more, include user-specified parameters,
and control the numbers of repeats and output directories. Please
see the :doc:`full file format description <fileformats/design>`
for more information.

Understanding the output
========================

The output of ``metawards`` is primarily a trajectory of the outbreak
through the model population. This is reported daily from the first
day (day 0) until the outbreak ends. For single runs this is printed
to the screen, e.g.

::

   21 58
  S: 56081959    E: 52    I: 17    R: 49    IW: 9   TOTAL POPULATION 56082025

A model run moves individuals between different states according to
whether they become infected, and then progress through the outbreak.
The codes mean;

* **S**: The number of the population who are *susceptible* to infection
* **E**: The number of the population who are *latent*, meaning they are
  infected, but not yet infectious.
* **I**: The number of the population who are *infected*, meaning they
  have symptoms and are infectious
* **R**: The number of the population who are removed from being susceptible,
  either because they have been newly infected that day, or because they
  have recovered from the infection and are no longer susceptible to infection
* **IW**: The number of electoral wards that contain at least one
  individual who was newly infected that day.

As well as being printed to the screen, this data is also written
to the CSV file ``output/results.csv.bz2`` for easy reading and analysis
using R or Python pandas.

If multiple model runs are performed, then each run is given a fingerprint
based on the adjustable parameters and repeat number. The output is written
to ``output/[fingerprint]_[repeat]/output.txt``. In addition, results
for all of the runs are combined into a single ``output/results.csv.bz2``
file for easy combined analysis.

For example, this output could be read into a pandas dataframe using

.. code-block:: python

    import pandas as pd

    df = pd.read_csv("output/results.csv.bz2")

    df # perform analysis

We run a good online workshop on
`how to use pandas for data analysis <https://milliams.com/courses/data_analysis_python/index.html>`__.

.. note::

  The ``E``, ``I`` and ``R`` stages are just the defaults, and ``metawards``
  does support custom disease stages which may not have ``E``, ``I`` or ``R``.
  To learn more, take a look at the
  :doc:`tutorial on custom named stages <tutorial/part07/05_named_stages>`.

metawards-plot
--------------

For quick and simple plots, ``metawards`` comes with the command-line
program ``metawards-plot``. This can be used to create plots of the
data in the ``output/results.csv.bz2`` file using, e.g.

.. code-block:: bash

   metawards-plot --input output/results.csv.bz2

For full help on this program using ``metawards-plot --help``.
The full help is :doc:`available here <metawards_plot_help>`.

.. toctree::
   :hidden:

   metawards_plot_help
========
Features
========

.. figure:: images/uq4map.jpg
    :width: 100%
    :figwidth: 35%
    :align: right
    :alt: Map of an outbreak modelled in England and Wales
    :figclass: align-right

    `Analysis of a simulation <https://uq4covid.github.io/vignettes/metawards_plot>`__
    run that used MetaWards to chart disease spread across England and Wales.

MetaWards is a SIR-based metapopulation disease modelling program. The
software models individuals who move as metapopulations between
home wards and work (or play) wards. This is a highly flexible program,
that can model custom geographies, networks, diseases and demographics.
You can use MetaWards directly from within Python (including Jupyter),
R (including RStudio) or from the command line (terminal / console).

**For more information, take a look at the** :doc:`quick start guide <quickstart/index>`.

The program is designed to enable researchers to model how an infection
may spread, and what impact different control measures or individual
behaviours may make.

To this end, MetaWards features;

.. figure:: images/demographic_model.jpg
    :width: 100%
    :figwidth: 35%
    :align: right
    :alt: Output during modelling different self-isolation durations
    :figclass: align-right

    Multiple connected networks enable complex scenarios to be modelled

* a flexible plugin architecture that makes it easy to implement new
  control measures. For example, the tutorial shows how this can be
  used to
  :doc:`model shielding <tutorial/part05/03_custom>`,
  :doc:`different lockdown scenarios <tutorial/part03/06_scan_lockdown>`
  and to
  :doc:`investigate necessary durations of quarantine or self-isolation <tutorial/part06/02_duration>`.

* multi-network demographic support. Multiple networks can be run as
  a single combined group, with custom plugins used to merge data
  between networks, and conditionally move individuals between
  different demographics. We've used to model shielding and self-isolation,
  hospital admissions, impact of individuals returning from holidays etc.

.. figure:: images/tutorial_5_3_1_demographics.jpg
    :width: 100%
    :figwidth: 35%
    :align: right
    :alt: Results of a shielding experiment
    :figclass: align-right

    Built in metawards-plot tool for rapid visualisation of results,
    including across multiple networks

* per-ward custom parameter support. Different wards can have different
  parameters, meaning that you can easily model local behaviour
  (e.g. local lockdowns, changes in local control measures etc.).

* complete-detail and full control over horizontal and vertical movements
  through disease stages or across demographics. We've used this to
  :doc:`model vaccination <tutorial/part09/02_vaccinate>` and also
  to :doc:`model waning immunity <tutorial/part09/03_fading_immunity>`,
  with individuals returned from the R or V stages back to S.

* flexible data output support - again handled using an array of in-built
  or user-supplied data extraction plugin functions. Output the data you
  need in the format you want to perform analysis.

* full reproducibility support baked in throughout. The code records
  enough data to make reproduction easy, with results designed to
  be the same given the same inputs, random number seed and number of threads.

* a :doc:`complete tutorial <tutorial/index>` that takes you from beginning to
  learn how to run SIR simulations, to writing powerful plugins that let you
  model complex scenarios.

* support for scanning design files for optimisation or sensitivity analysis
  of nearly all input parameters, plus any user custom parameters used
  in the main code or any plugins. These scans can use as much compute
  you have available, parallelising individual runs over multiple cores,
  and scaling multiple runs up to full supercomputers (if you are lucky
  enough to have access to one)...

.. figure:: images/rstudio.jpg
    :width: 100%
    :figwidth: 35%
    :align: right
    :alt: RStudio screenshot
    :figclass: align-right

    Run MetaWards and analyse results directly within RStudio via the
    MetaWards R interface.

* ...but - the individual code is optimised and can run happily on small
  laptops. Individual national-scale networks fit in approximately
  80 MB of memory, and model runs can take 15-90 seconds to perform.
  This scales with the number of
  demographics that are added, but high performance and low memory consumption
  are design goals. Models using only a few wards are kilobytes, and take
  less than a second.

* :doc:`flexible input files <model_data>` that would enable modelling of any
  region or country to be undertaken (subject to good input data).
  Models of the UK and England and Wales have been created, and a
  Python and R API are provided to make it easy to create custom networks.
  These can model everything from individual wards or local geographies,
  up to full national- or international-scale metapopulation models.

.. figure:: images/pandas_example.jpg
    :width: 100%
    :figwidth: 35%
    :align: right
    :alt: Analysing data in pandas
    :figclass: align-right

    Easily load compressed data files into pandas, R or Excel for analysis

* a colourful, modern and informative console output, with full unicode
  support and progress indicators. All outputs are duplicated to
  text files to ensure that no data from a run is lost.

* both command-line and API interfaces. Feel free to run MetaWards as
  a standalone program, or to use the Python or R API to embed it as part
  of a larger framework. A modular, robust design has been used, so
  feel free to take and re-use the parts of MetaWards that are most
  of use to you.

.. figure:: images/parallel_output.jpg
    :width: 100%
    :figwidth: 35%
    :align: right
    :alt: Running multiple jobs in parallel
    :figclass: align-right

    Perform multiple runs in parallel, with live summary updates as
    jobs finish

Software design
---------------

MetaWards is a modern piece of software that has been engineered following
recognised best practice, e.g. using a modular design,
lots of documentation, copious run-time and unit tests, and following
a "tutorial-driven" development philosophy.

The software is mostly Python, with C used (via cython) to accelerate
key parts. An :doc:`R interface <quickstart/01_R>`
is provided via `reticulate <https://rstudio.github.io/reticulate/>`__.
The code is parallelised using `OpenMP <https://openmp.org>`__,
with multiple model runs parallelised using multiprocessing,
`scoop <https://scoop.readthedocs.io>`__ or
`MPI, via mpi4py <https://mpi4py.readthedocs.io>`__.

We take testing very seriously, and have lots of unit, integration and
run-time tests. These are run as part of our CI/CD system deployed
via `GitHub Actions <https://github.com/metawards/MetaWards/actions>`__
The code has in-built developer
:doc:`support for debug logging and profiling <devsupport>`, with
:doc:`full API docs available <api/index>`
that we hope will help new developers get quickly up to speed.

We also take versioning and backwards compatibility seriously. We follow
the `semantic versioning <https://semver.org>`__  system for the main API,
which should give confidence
to anyone wanting to build on top of MetaWards. We also maintain
compatibility of inputs and outputs, with the strong aim that all tutorials
will be runnable, as-is, in future versions of the code.
=================
Developer's guide
=================

The source code for MetaWards is available on
`GitHub <https://github.com/metawards/MetaWards>`__.

The data needed to run a MetaWards simulation is also on GitHub
in the `MetaWardsData <https://github.com/metawards/MetaWardsData>`__
repository.

Setting up your computer
=========================

MetaWards requires Python >= 3.7, so please install this before continuing
further. `Anaconda <https://anaconda.org>`__ provides easy installers for
Python that work on a range of operating systems and that don't need
root or admin permissions.

Virtual environments
--------------------

It is recommended that you develop MetaWards in a Python
`virtual environment <https://docs.python.org/3/tutorial/venv.html>`__.
You can create a new environment in the directory ``venvs/metawards-devel``
by typing;

.. code-block:: bash

   mkdir venvs
   python -m venv venvs/metawards-devel

Feel free to place the environment in any directory you want.

Virtual environments provide sandboxes which make it easier to develop
and test code. They also allow you to install Python modules without
interfering with other Python installations.

You activate you environment by typing;

.. code-block:: bash

    source venvs/metawards-devel/bin/activate

This will update your shell so that all python commands (such as
``python``, ``pip`` etc.) will use the virtual environment. You can
deactivate the environment and return to the "standard" Python using;

.. code-block:: bash

   deactivate

If you no longer want the environment then you can remove it using

.. code-block:: bash

  rm -rf venvs/metawards-devel

Developer dependencies
----------------------

If is recommended that you have the following modules installed when
developing metawards;

* cython
* numpy
* flake8
* pytest
* sphinx (plus sphinx_issues sphinx_rtd_theme)

You can install these manually, or all at once using;

.. code-block:: bash

   pip install -r requirements-dev.txt

Coding Style
============

MetaWards is written in Python 3 (>= 3.7), with time-critical functions
written using `Cython <https://cython.org>`__ and parallelised using
`OpenMP <https://openmp.org>`__. The Cython code is written strictly
in comformant C, meaning that the package should compile and work on
any system on which Python >= 3.7 runs. We ourselves are running production
MetaWards models on ARM64 on Linux, and develop on X86-64 on Linux and
Mac laptops. The program is tested with CI/CD on Windows 10, and
we thank the windows users who've helped us make the tutorial
cross-platform.

The aim of the Python port is to provide a simple and robust API that
is a strong foundation for robust growth and scale-up of MetaWards, and
one in which unnecessary implementation details are hidden from the user.

We aim as much as possible to follow a
`PEP8 <https://www.python.org/dev/peps/pep-0008/>`__ python coding style and
recommend that developers install and use
a linter such as `flake8 <https://flake8.pycqa.org/en/latest/>`__.

For the Cython pyx code we also try to maintain a PEP8 style where possible,
and recommend using a tool such as `autopep8 <https://pypi.org/project/autopep8/>`__
to keep this style (it is used by the lead developers, so contributions
will be formatted using it eventually ;-)).

We require that non-python code is strictly C. While C++ is an excellent
language, it is too bulky for use in MetaWards and makes it more
challenging to create portable binary distributions.

For ease of installation and support, we also minimise or bundle
external dependencies (e.g. we use a
`bundled version <https://github.com/metawards/MetaWards/tree/devel/src/metawards/ran_binomial>`__
of the binomial random number generator from `numpy <https://numpy.org>`__).
This code has to run on a wide variety of architectures, operating
systems and machines - some of which don't have any graphic libraries,
so please be careful when adding a dependency.

With this in mind, we use the following coding conventions:

Naming
------

We follow a Python style naming convention.

* Packages: lowercase, singleword
* Classes: CamelCase
* Methods: snake_case
* Functions: snake_case
* Variables: snake_case
* Source Files: snake_case with a leading underscore

Functions or variables that are private should be named with a leading
underscore. This prevents them from being prominantly visible in Python's
help and tab completion.

Relative imports should be used at all times, with imports ideally delayed
until they are needed.

For example, to import the Network object into a function that is in the
utils module, you would type

.. code-block:: python

   from .._network import Network

   network = Network()
   network.run(...)

or to import the Parameters from code that is the main MetaWards package,
you would type

.. code-block:: python

   from ._parameters import Parameters

   parameters = Parameters()
   parameters.add_seeds(filename="additional_seeds.dat")

Note that many classes are Python
`dataclasses <https://realpython.com/python-data-classes/>`__, which
are really useful for quick and safe development of code.
Python dataclasses should be preferred over writing your own
data-style classes.

Modules
-------

MetaWards consists of the main module, e.g. ``metawards``, plus
a ``metawards.utils`` module that contains useful utilities.

The main module should be the focus of external developers, while
the utils module should only be needed by developers of metawards
itself.

In addition, there is a ``metawards.app`` module which contains the
code for the various command-line applications (e.g. the
metawards executable).

MetaWards uses a plugin-style interface for most new development.
The plugins are for :doc:`iterators <api/index_MetaWards_iterators>`,
:doc:`extractors <api/index_MetaWards_extractors>`,
:doc:`mixers <api/index_MetaWards_mixers>` and
:doc:`movers <api/index_MetaWards_movers>`. Ideally, most new code
should be added as one of these plugins. If you can't fit your code
into a plugin then please
`raise an issue <https://github.com/metawards/MetaWards/issues>`__
to discuss your idea with the core developers, so that a way
forward can be found. We really appreciate your help, and want
to make sure that your ideas can be included in the most compatible
way in the code.

To make MetaWards easy for new developers
to understand, we have a set of rules that will ensure that only
necessary public functions, classes and implementation details are
exposed to the Python help system.

* Module files containing implementation details are prefixed with
  an underscore, i.e. ``_parameters.py``

* Where possible, external packages are hidden inside each module,
  e.g. ``import sys as _sys``

* Each module file contains an ``__all__`` variable that lists the
  specific items that should be imported.

* The package ``__init__.py`` can be used to safely expose the required
  functionality to the user with:

.. code-block:: python

   from module import *

This results in a clean API and documentation, with all extraneous information,
e.g. external modules, hidden from the user. This is important when working
interactively, since `IPython <https://ipython.org>`__
and `Jupyter <https://jupyter.org>`__
do not respect the ``__all__`` variable when auto-completing, meaning that the
user will see a full list of the available names when hitting tab. When
following the conventions above, the user will only be able to access the
exposed names. This greatly improves the clarity of the package, allowing
a new user to quickly determine the available functionality. Any user wishing
expose further implementation detail can, of course, type an underscore to
show the hidden names when searching.

Workflow
========

Feature branches
----------------

First make sure that you are on the development branch of MetaWards:

.. code-block:: bash

   git checkout devel

Now create and switch to a feature branch. This should be prefixed with
*feature*, e.g.

.. code-block:: bash

   git checkout -b feature-process

While working on your feature branch you won't want to continually re-install
in order to make the changes active. To avoid this, you can either make use
of ``PYTHONPATH``, e.g.

.. code-block:: bash

   make
   PYTHONPATH=./build/lib.{XXX} python script.py
   PYTHONPATH=./build/lib.{XXX} pytest tests

(where ``{XXX}`` is the build directory for Cython on your computer, e.g.
``./build/lib.macosx-10.9-x86_64-3.7`` - remember that you need to
type ``make`` to rebuild any Cython code and to copy your updated
files into that directory)

or use the ``develop`` argument when running the ``setup.py`` script, i.e.

.. code-block:: bash

   python setup.py develop

(this installs your current version of metawards into your current python
environment)

Testing
=======

When working on your feature it is important to write tests to ensure that it
does what is expected and doesn't break any existing functionality. Tests
should be placed inside the ``tests`` directory, creating an appropriately
named sub-directory for any new packages.

The test suite is intended to be run using
`pytest <https://docs.pytest.org/en/latest/contents.html>`__.
When run, ``pytest`` searches for tests in all directories and files
below the current directory, collects the tests together, then runs
them. Pytest uses name matching to locate the tests. Valid names start
or end with *test*\ , e.g.:

::

   # Files:
   test_file.py       file_test.py

.. code-block:: python

   # Functions:
   def test_func():
      # code to perform tests...
      return

   def func_test():
      # code to perform tests...
      return

We use the convention of ``test_*`` when naming files and functions.

Running tests
-------------

To run the full test suite, simply type:

.. code-block:: bash

   pytest tests

To run tests for a specific sub-module, e.g.:

.. code-block:: bash

   pytest tests/utils

To only run the unit tests in a particular file, e.g.:

.. code-block:: bash

   pytest tests/test_integration.py

To run a specific unit tests in a particular file, e.g.:

.. code-block:: bash

   pytest tests/test_read_variables.py::test_parameterset

To get more detailed information about each test, run pytests using the
*verbose* flag, e.g.:

.. code-block:: bash

   pytest -v

More details regarding how to invoke ``pytest`` can be
found `here <https://docs.pytest.org/en/latest/usage.html>`__.

Writing tests
^^^^^^^^^^^^^

Basics
""""""

Try to keep individual unit tests short and clear. Aim to test one thing, and
test it well. Where possible, try to minimise the use of ``assert`` statements
within a unit test. Since the test will return on the first failed assertion,
additional contextual information may be lost.

Floating point comparisons
""""""""""""""""""""""""""

Make use of the
`approx <https://docs.pytest.org/en/latest/builtin.html#comparing-floating-point-numbers>`__
function from the ``pytest`` package for performing floating
point comparisons, e.g:

.. code-block:: python

   from pytest import approx

   assert 0.1 + 0.2 == approx(0.3)

By default, the ``approx`` function compares the result using a
relative tolerance of 1e-6. This can be changed by passing a keyword
argument to the function, e.g:

.. code-block:: python

   assert 2 + 3 == approx(7, rel=2)

Skipping tests
""""""""""""""

If you are using
`test-driven development <https://en.wikipedia.org/wiki/Test-driven_development>`__
it might be desirable to write your tests before implementing the functionality,
i.e. you are asserting what the *output* of a function should be, not how it should
be *implemented*. In this case, you can make use of
the ``pytest`` *skip* decorator
to flag that a unit test should be skipped, e.g.:

.. code-block:: python

   @pytest.mark.skip(reason="Not yet implemented.")
   def test_new_feature():
       # A unit test for an, as yet, unimplemented feature.
       ...

Parametrizing tests
"""""""""""""""""""

Often it is desirable to run a test for a range of different input parameters.
This can be achieved using the ``parametrize`` decorator, e.g.:

.. code-block:: python

   import pytest
   from operator import mul

   @pytest.mark.parametrize("x", [1, 2])
   @pytest.mark.parametrize("y", [3, 4])
   def test_mul(x, y):
       """ Test the mul function. """
       assert mul(x, y) == mul(y, x)

Here the function test_mul is parametrized with two parameters, ``x`` and ``y``.
By marking the test in this manner it will be executed using all possible
parameter pairs ``(x, y)``\ , i.e. ``(1, 3), (1, 4), (2, 3), (2, 4)``.

Alternatively:

.. code-block:: python

   import pytest
   from operator import sub
   @pytest.mark.parametrize("x, y, expected",
                           [(1, 2, -1),
                            (7, 3,  4),
                            (21, 58, -37)])
   def test_sub(x, y, expected):
       """ Test the sub function. """
       assert sub(x, y) == -sub(y, x) == expected

Here we are passing a list containing different parameter sets, with the names
of the parameters matched against the arguments of the test function.

Testing exceptions
""""""""""""""""""

Pytest provides a way of testing your code for known exceptions. For example,
suppose we had a function that raises an ``IndexError``\ :

.. code-block:: python

   def indexError():
       """ A function that raises an IndexError. """
       a = []
       a[3]

We could then write a test to validate that the error is thrown as expected:

.. code-block:: python

   def test_indexError():
       with pytest.raises(IndexError):
           indexError()

Custom attributes
"""""""""""""""""

It's possible to mark test functions with any attribute you like. For example:

.. code-block:: python

   @pytest.mark.slow
   def test_slow_function():
       """ A unit test that takes a really long time. """
       ...

Here we have marked the test function with the attribute ``slow`` in order to
indicate that it takes a while to run. From the command line it is possible
to run or skip tests with a particular mark.

.. code-block:: bash

   pytest mypkg -m "slow"        # only run the slow tests
   pytest mypkg -m "not slow"    # skip the slow tests

The custom attribute can just be a label, as in this case, or could be your
own function decorator.

Continuous integration and delivery
-----------------------------------

We use GitHub Actions to run a full continuous integration (CI)
on all pull requests to devel and
main, and all pushes to devel and main. We will not merge a pull
request until all tests pass. We only accept pull requests to devel.
We only allow pull requests from devel to main. In addition to CI,
we also perform a build of the website on pushes to devel and tags
to main. The website is versioned, so that old the docs for old
versions of the code are always available. Finally, we have set up
continuous delivery (CD) on pushes to main and devel, which build the
pypi source and binary wheels for Windows, Linux (manylinux2010)
and OS X. These are manually uploaded to pypi when we tag
releases, but we expect to automate this process soon.

Documentation
=============

MetaWards is fully documented using a combination of hand-written files
(in the ``doc`` folder) and auto-generated api documentation created from
`NumPy <https://numpy.org>`__ style docstrings.
See `here <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`__
for details. The documentation is automatically built using
`Sphinx <http://sphinx-doc.org>`__ whenever a commit is pushed to devel, which
will then update this website.

To build the documentation locally you will first need to install some
additional packages.

.. code-block:: bash

   pip install sphinx sphinx_issues sphinx_rtd_theme

Then move to the ``doc`` directory and run:

.. code-block:: bash

   make html

When finished, point your browser to ``build/html/index.html``.

Committing
==========

If you create new tests, please make sure that they pass locally before
commiting. When happy, commit your changes, e.g.

.. code-block:: bash

   git commit src/metawards/_new_feature.py tests/test_feature \
       -m "Implementation and test for new feature."

Remember that it is better to make small changes and commit frequently.

If your edits don't change the MetaWards source code, or documentation,
e.g. fixing typos, then please add ``ci skip`` to your commit message, e.g.

.. code-block:: bash

   git commit -a -m "Updating docs [ci skip]"

This will avoid unnecessarily running the
`GitHub Actions <https://github.com/metawards/MetaWards/actions>`__, e.g.
building a new MetaWards package, updating the website, etc.
(the GitHub actions are configured in the file
``.github/workflows/main.yaml``). To this end, we
have provided a git hook that will append ``[ci skip]`` if the commit only
modifies files in a blacklist that is specified in the file ``.ciignore``
(analagous to the ``.gitignore`` used to ignore untracked files). To enable
the hook, simply copy it into the ``.git/hooks`` directory:

.. code-block:: bash

    cp git_hooks/commit-msg .git/hooks

Any additional files or paths that shouldn't trigger a re-build can be added
to the ``.ciignore`` file.

Next, push your changes to the remote server, e.g.

.. code-block:: bash

   # Push to the feature branch on the main MetaWards repo, if you have access.
   git push origin feature

   # Push to the feature branch your own fork.
   git push fork feature

When the feature is complete, create a *pull request* on GitHub so that the
changes can be merged back into the development branch.
For information, see the documentation
`here <https://help.github.com/articles/about-pull-requests>`__.

Thanks
======

First, thanks to you for your interest in MetaWards and for reading this
far. We hope you enjoy having a play with the code and having a go
at adding new functionality, fixing bugs, writing docs etc.

We would also like to thank Lester Hedges and the
`BioSimSpace <https://biosimspace.org>`__ team who provided great advice
to set up the above, and from whose
`GitHub repo <https://github.com/michellab/biosimspace>`__
most of the procedures, scripts and documentation above is derived.
================
Help and support
================

We welcome new users and hope that we can make your experience of using
and working with MetaWards as easy and productive as possible.

We've followed a "tutorial-driven development" model, whereby we write
tutorials for all new functionality that is included. Hopefully that
tutorial has helped you get up and running with MetaWards.

If not, or you are having other problems, then please feel free to
get in touch with us by
`raising an issue <https://github.com/metawards/MetaWards/issues>`__
on our public GitHub repo. We will endeavour to respond as quickly
as we can.

We have provided some simple templates that should make it easier for
you to give us the information we need to answer your question efficiently.
If these don't fit, then please feel free to use the miscellaneous
template.

We strive to make the MetaWards community a welcoming, supportive and
helpful place, so please do make sure you read and follow our
`code of conduct <https://github.com/metawards/MetaWards/blob/devel/CODE_OF_CONDUCT.md>`__
when posting issues or otherwise interacting with anyone in our community.
================
User input files
================

The user input file is used to set the initial values of any parameter
that you could set in a :doc:`design file <design>`. It uses an almost
identical parsing code, and indeed, a design file with a single
row of data is also a valid user input file.

However, while you can read column-orientated data, a user input file
is best written as a row-orientated file. This file has one variable
per row, using either a space, comma, colon or equals sign to
separate the name of the variable from its initial value.

Examples
--------

::

  .isolate_ndays = 7
  .isolate_stage = 3

Set the initial value of ``isolate_ndays`` to ``7`` and
``isolate_stage`` to ``3``.

::

  # Number of days to isolate
  .isolate_ndays: 7

  # Stage at which to start isolating 
  .isolate_stage: 3

Same, except using colon as the separator, and adding comments and blank
lines to improve legibility.

::

  hospital:beta["ICU"]: 0.3
  .lockdown_start: d"next week"

Setting the initial beta parameter of the ``ICU`` stage in the ``hospital``
demographic to ``0.3``, while setting lockdown to start using the
fuzzy date ``next week``.
===================
Design (scan) files
===================

The design (or scan) file specifies which adjustable variables should be
changed during
a ``metawards`` calculation, and which value(s) they should be changed to.

This is a column-orientated flexible-format file, using one column
per adjustable variable. Like all flexible-format files you can use
commas or spaces to separate columns. The separators used in the first line
will be assumed to be used for the rest of the file.

Column headers
--------------

You should name the variables you want to adjust in the column headers.
The format of these is;

* ``variable`` : Just the variable name will adjust that single variable
  for all copies in all demographics
* ``demographic:variable`` : This will adjust the single variable in only
  the demographic named ``demographic``.
* ``variable[1]`` : This will adjust index ``1`` of the named variable
  in all demographics.
* ``variable["E"]`` : This will adjust key ``E`` of the named variable,
  e.g. where the key refers to the disease stage.
* ``demographic:variable["KEY"]`` : This will adjust key ``KEY`` of the
  named variable in the named demographic.
* ``user.myvariable`` : This will adjust a custom user variable called
  ``myvariable`` in all demographics.
* ``.myvariable`` : This will also adjust a custom user variable called
  ``myvariable`` in all demographic (we can drop the ``user`` to save
  space).
* ``demographic:.myvariable[0]`` : This will adjust index 0 of the custom
  user variable ``myvariable`` in the demographic called ``demographic``.

The full list of in-built parameters you can adjust are listed below;

.. program-output:: python get_variableset_help.py

In addition, you can create as many custom-user parameters to adjust
as you would like.

Special columns
---------------

There are two special columns;

* ``repeats`` : Specify the number of repeats of this row of values. This
  should be an integer (whole number) and gives more fine-grained control
  to specifying the number of repeats for different rows in a design.
* ``output``: Specify the output directory in which to place the output
  of the model run using this row of adjustable variables. This overrides
  the default output directory, which is named using the fingerprint
  for this set of adjustable parameters. Use this if you want to have
  control over where all of the output is written. Note that ``metawards``
  does not allow output for multiple model runs to be written to the
  same file, and will append numbers (e.g. ``x002``, ``x003``) to
  any duplicated names.

Default columns
---------------

We highly recommend that you name the columns in your design file. If you
don't, then the default columns will be used. These are;

::

   beta[2]  beta[3]  progress[1]  progress[2]  progress[3]

.. note::

   These defaults came from the original C version of MetaWards, and are
   retained so we keep backwards compatibility with the original
   input files.


Column data
-----------

Values in each column will be interpreted using :class:`metawards.Interpret`.
The code will try to guess the most appropriate data type, moving through
simple numbers,
then isoformat dates, then number expressions (including random number),
then booleans, before returning the data as a string.

You can force the type of the data using numpy-style quote characters. The
type characters recognised are;

* ``d"value"`` : forces ``value`` to be interpreted as a date (isoformat or
  fuzzy date via the `Python dateparser module <https://dateparser.readthedocs.io/en/latest/>`__).
* ``f"value"`` : forces ``value`` to be interpreted as a number (floating point
  or integer)
* ``i"value"`` : forces ``value`` to be interpreted as an integer (whole number)
* ``b"value"`` : forces ``value`` to be interpreted as a boolean (true or
  false) value. The code recognises ``true/false``, ``on/off`` and ``yes/no``,
  as well as ``1/0`` (all case-insensitive).
* ``s"value"`` : forces ``values`` to be a string. Use this if you really want
  a string or want to interpret the value yourself in your plugin code.

Examples
--------

::

  beta[1]  beta[2]  beta[3]
    0.5      0.5      0.5
    0.6      0.7      0.8

Two sets of variables, setting ``beta[1]``, ``beta[2]`` and ``beta[3]`` to
``0.5`` in the first set, and ``0.6``, ``0.7`` and ``0.8`` in the second
set.

::

   beta["I1"]  beta["I2"]  beta["I3"]
   # initial baseline
     0.5        0.5        0.5

   # increasing infectivitiy
     0.6        0.7        0.8

Same, except using the names of the disease stages and adding comments
and extra whitespace to aid legibility.

::

   .lockdown_start    .scale_rate   repeats   output
    d"March 15 2020"      0.2         5       lockdown_march
    d"April 1 2020"       0.1         3       lockdown_april

Setting the user parameters ``lockdown_start`` and ``scale_rate``
to either March 15th 2020 / 0.2 or April 1st 2020 / 0.1. Asking for
5 repeats of the March 15 data, outputting the results to directories
called ``lockdown_march``, and 3 repeats of the April lockdown,
outputting the data to ``lockdown_april``.

.. note::

   The first repeats will be written to ``lockdown_march`` and
   ``lockdown_april``, while the remainder will be written to
   numbered directories, e.g. ``lockdown_marchx002`` etc.

::

   beta[2], beta[3], progress[1], progress[2], progress[3]
    0.95,   0.95,     0.19,        0.91,         0.91
    0.90,   0.93,     0.18,        0.92,         0.90

Named columns, using comma separators, and adding extra spaces for
improved legibility.
=================
Extra seeds files
=================

This file specifies how infections should be seeded in an outbreak.
This is a column-based flexible-format file, with each line of the
file containing information about a specific seeding event.

There should be up to three pieces of information per line;

1. ``day`` : The day or date to seed the infection(s). This should be
   an integer or a date that is interpreted via
   :func:`Interpret.day_or_date <metawards.Interpret.day_or_date>` function.
   This either sets the day number on which to seed, or the exact date
   in which to seed.
2. ``ward`` : The ward (either index or name) in which the infection
   will be seeded. This should either be an integer (whole number) specifying
   the index of the ward (1-indexed), or a string that identifies the
   ward by name, e.g. ``Clifton/Bristol`` would look for the ward called
   ``Clifton`` in the authority ``Bristol``. The ward is looked up by
   name using the :func:`WardInfos.find <metawards.WardInfos.find>`
   function.
3. ``number`` : The number of infections to be seeded. This is an integer
   that is interpreted via the
   :func:`Interpret.integer <metawards.Interpret.integer>` function.
4. ``demographic`` (optional) : The demographic (either index or name)
   in which the infection will be seeded. If this is not specified,
   then the first demographic in the network is seeded. If this is
   the index it is interpreted via
   :func:`Interpret.integer <metawards.Interpret.integer>`, while if this
   is the name then this will direct look-up the demographic by name.

This is a column-formatted file. Like all flexible-format files you can use
commas or spaces to separate columns. The separators used in the first line
will be assumed to be used for the rest of the file. You can order the
columns however you wish, as long as you provide column headers
(``day``, ``ward``, ``number``, and (optionally) ``demographic``). Otherwise
the column orders are ``day``, ``number``, ``ward``,
(optional) ``demographic``.

You can add extra spaces and blank lines to make the file more readable,
and can add comments via the ``#`` character.

Examples
--------

::

   1   5   1

Seed 5 infections on day 1 (first day) of the outbreak in ward 1.

::

   1, 5, 1

Same, but using commas to separate

::

   day, ward, number
   # seed 5 infections on day 1 in ward 1
   1, 5, 1

Same, but naming the columns and adding a comment

::

   # Seeding two locations in Bristol

   day, demographic, ward, number
   tomorrow, 0, Clifton / Bristol, "rand(5, 20)"

   next week, 0, Knowle / Bristol, 5 + 3

Adding comments and changing the order of columns. Seeding by fuzzy dates,
e.g. ``tomorrow`` and ``next week`` in wards identified by
``name / authority``. Number to be seeded is a random number from 5 to 20
(inclusive) in Clifton, Bristol, and the result of the expression
``5 + 3`` (8) in Knowle.

::

   day             ward       number
   2020-12-05      1          5
   2020-12-10      1          5
   2020-12-15      1          5

Seeding using three `isoformat dates <https://pythontic.com/datetime/date/isoformat>`__
dates (5th, 10th and 15th December 2020), all in ward 1, seeding 5 infections
each day.

=======================
Input files and formats
=======================

Input data in ``metawards`` is provided via a variety of different files,
many of which are described below.

There are three types of input files:

1. Flexible-format: These are a very flexible format, and support a wide
   range of input data types and layout options. Examples include
   the :doc:`extra seeds <extraseeds>`, :doc:`design <design>`
   and :doc:`user input <userinput>` files.

2. Rigid-format: These have a rigid format, which is specific to their
   type. The main examples of this are the files used to specify
   a model (e.g. :doc:`the network, connections etc. <network>`)

3. JSON-format: These are files that are in standard
   `JSON <https://en.wikipedia.org/wiki/JSON>`__ format. Examples include
   the disease and demographics files, plus many of the files
   in MetaWardsData (e.g. the description of the model data, and the
   static parameters file).

Flexible-format files
---------------------

Flexible-format files are read using the
`Python CSV module <https://docs.python.org/3/library/csv.html>`__.
These files are either column or row based (depending on the file type),
and you are free to use
a comma or spaces as the separator (but must be consistent within a file).
Comments can be added using a ``#`` character, and blank lines are
ignored.

Data within a flexible-format file is interpreted using
:class:`metawards.Interpret`. This can interpret simple data, such as
strings, numbers, booleans (``true`` or ``false``), as well as complex
data such as dates (``next week``, ``January 10 2020``), expressions
(``10.0 / 3.0``, ``pi * 2.3**2``) and random numbers
(``rand(0,10)``, ``rand()``, ``rand(2.5, 2.6)``).

Rigid-format files
------------------

Rigid-format files are read using a custom parser for each file type. As
such, the files have a rigid format that is specified for each file type.
We plan to migrate as many rigid-format files across to either
flexible-format or JSON-format as possible.

JSON-format files
-----------------

JSON-format files are standard `JSON <https://en.wikipedia.org/wiki/JSON>`__
files that are used for small or less-structured files, e.g. specifying
the parameters for a disease, or specifying the data associated with
different demographics. These files are read using the
`Python JSON module <https://docs.python.org/3/library/json.html>`__
into dictionaries, which are interpreted by the classes associated with
each file (e.g. :class:`~metawards.Disease` in the case of disease parameters).
Many of these classes use :class:`metawards.Interpret` to interpret the
data from the JSON file, meaning that these support expressions, random
numbers etc. We plan to ensure that as much data as possible is interpreted
using :class:`metawards.Interpret`, so that there is a consistent
experience across the code.

Detailed descriptions
---------------------

.. toctree::
   :maxdepth: 2

   network
   extraseeds
   design
   userinput
=============
Network Files
=============

Network files are those that define the nodes (wards) of the network,
and the links (communiting/movement) of individuals between nodes.

The Network files are a collection of files that must all exist in the
same directory. They comprise;

* ``description.json``: This file must exist in the directory, and provides
  the metadata needed to locate all of the other necessary files.
* ``work_size / play_size``: These files list the index and population of
   workers (in ``work_size``) and players (in ``play_size``) in
   each ward in the network. The total population in each ward (node) is
   the sum of the work and play populations.
* ``work``: This file lists all of the work connections between nodes, giving
  the population of workers who commute from the source node to the
  destination.
* ``play``: This file lists all of the play connections between nodes, giving
  the **proportion** of players who make random movements from the source
  node to the destination node.
* ``position``: This file contains the coordinates of the centre of each node.
  These coordinates can be X/Y or latitude/longitude coordinates.
* ``lookup``: This file contains the names and other metadata about each
  node, e.g. the name, local authority, region etc.

These files are all described below.

description.json
----------------

This is a JSON-formatted file that provides metadata about the network model
that is contained in the directory. The file describes a simple dictionary
of the following key-value pairs;

* ``work``: Filename of the ``work`` file.
* ``work_size``: Filename of the ``work_size`` file.
* ``play``: Filename of the ``play`` file.
* ``play_size``: Filename of the ``play_size`` file.
* ``position``: Filename of the ``position`` file.
* ``coordinates``: Whether the coordinates in the ``position`` file are
  in x/y (``x/y``) or latitude/longitude (``lat/long``) format.
* ``coordinate_units``: (Optional) - the units for x/y coordinates. This should
  be either ``m`` for meters, or ``km`` for kilometers. Distances in
  ``metawards`` are always reported as kilometers.
* ``lookup``: (Optional) Filename of the ``lookup`` file.
* ``lookup_columns``: (Optional) A dictionary that gives the column numbers
  (zero-indexed)for the ``code``, ``name``, ``alternate_code``,
  ``alternate_name``, ``authority_code`` and ``authority_name`` fields for
  the ward metadata.

work
----

The ``work`` file contains the list of all work connections between nodes,
giving the number of individuals who commute from the source node to the
destination node. This is a column-based file with three columns of numbers
that can be space or comma separated. For example;

::

    1 1 290
    1 2 3
    1 5 139
    1 6 59
    1 7 17
    1 8 119
    1 9 37
    1 10 121

The first number is the (1-indexed) index of the source node, while the
second number is the (1-indexed) index of the destination node. The
third number is the number of individuals who commute from the source
to the destination node.

In this example, this lists workers who commute from ward 1 to wards
1 to 10. 290 workers commute from ward 1 to ward 1 (so work in the same
ward in which they live). 3 workers commute from ward 1 to ward 2,
139 commute from ward 1 to ward 5 etc.

.. note::

   All of the source wards have to be listed contiguously, i.e.
   you must list all of the connections where the source ward is
   equal to ``1`` before you list all of the connections where the
   source ward is equal to ``2`` etc.

work_size
---------

The ``work_size`` file contains the number of workers who reside in each ward.
This is a column-based file with two columns of numbers that can be space or
comma separated. For example;

::

    1 6800
    2 1091
    3 7148
    4 5684
    5 7226
    6 6561
    7 6904
    8 7213
    9 6715
    10 7452

The first number is the (1-indexed) index of the ward, while the second
is the number of workers in that ward. The number of workers in a ward must
equal the sum of the number of workers from the ``work`` file that say
that they commute from that ward.

In this example, ward 1 has 6800 workers, ward 2 has 1091 workers etc.

play
----

The ``play`` file contains the list of play connections between wards.
This is a column-based file with three columns of numbers that can be
space or comma separated. For example;

::

    1 1 0.0426470588235294
    1 2 0.000441176470588235
    1 5 0.0204411764705882
    1 6 0.00867647058823529
    1 7 0.0025
    1 8 0.0175
    1 9 0.00544117647058823
    1 10 0.0177941176470588

The first number is the (1-indexed) index of the source ward, while the
second is the (1-indexed) index of the destination ward. The third number
is the fraction of players who will randomly travel from the source
ward to the destination ward.

In this case, ``0.0426...`` of players will remain in ward 1, while
``0.00044...`` of players will randomly move from ward 1 to ward 2.

.. note::

   Be careful as the third number is a fraction (floating point number
   between 0 and 1) and not the number of players that move. Also note
   that the same ordering requirements as for the ``work`` file apply,
   namely that all connections for source ward 1 have to be listed before
   all connections for source ward 2 etc.

play_size
---------

The ``play_size`` file contains the number of players in each ward.
This is a column-based file with two columns of numbers that can be
space or comma separated. For example:

::

    1 8915
    2 374
    3 7012
    4 10579
    5 8703
    6 12257
    7 10533
    8 11259
    9 8592
    10 10999

The first number is the (1-indexed) index of the ward, while the second
is the number of players in that ward. In this case, there are 8915
players in ward 1, 374 players in ward 2 etc.

position
--------

The ``position`` file contains the coordinates of the center of each ward.
This is a column-based file with three columns of numbers that can be
space or comma separated. For example:

::

    1 524693.890435782 190136.324582048
    2 532169.852194767 181663.72329877
    3 522106.698233411 179737.792519091
    4 533388.404693453 193451.026467071
    5 525973.729674732 188464.548078951
    6 523953.505171282 187729.710513154
    7 520763.019103366 193360.422742775
    8 523100.570665414 189524.60815157
    9 523747.496973864 196638.656771657
    10 523019.549069799 192050.596775831

The coordinates are either x/y or latitude/longitude, depending on what
is set in the ``description.json`` file. If x/y, then the units are also
set in this file.

In this case, this is for x/y coordinates in meters. The first column
is the (1-indexed) index of the ward. The second column is the X coordinate,
while the third is the Y coordinate.

.. note::

  If lat/long was used, then the second column would be the latitude, and
  the third column the longitude.

So ward 1 has its center are (524.7km, 190.1km), while ward 2 has its center
at (532.1km, 181.7km) etc.

lookup
------

The ``lookup`` file contains metadata information about a ward that allows
it to be looked up by name, region etc. This is a column-based file with
a header row, and columns separated by spaces or commas. For example:

::

    WD11CD,WD11NM,WD11NMW,CMWD11CD,CMWD11NM,CMWD11NMW,IND,LAD11CD,LAD11NM,LAD11NMW,FID
    E05002337,Central,,E36000890,Central,,0,E06000039,Slough,,2001
    E05002338,Chalvey,,E36000891,Chalvey,,0,E06000039,Slough,,2002
    E05002339,Cippenham Green,,E36000892,Cippenham Green,,0,E06000039,Slough,,2003
    E05002340,Cippenham Meadows,,E36000893,Cippenham Meadows,,0,E06000039,Slough,,2004
    E05002341,Colnbrook with Poyle,,E36000894,Colnbrook with Poyle,,0,E06000039,Slough,,2005
    E05002342,Farnham,,E36000895,Farnham,,0,E06000039,Slough,,2006
    E05002343,Foxborough,,E36000896,Foxborough,,0,E06000039,Slough,,2007
    E05002344,Haymill,,E36000897,Haymill,,0,E06000039,Slough,,2008
    E05002345,Kedermister,,E36000898,Kedermister,,0,E06000039,Slough,,2009
    E05002346,Langley St Mary's,,E36000899,Langley St Mary's,,0,E06000039,Slough,,2010

The first row provides the names of each of the columns. These are ignored
by the code. Each row gives the information for the ward whose index is equal
to the line number (so line 1 give the information for ward at index 1).

The ``lookup_columns`` field in ``description.json`` specifies which columns
from this file to use for the name, code, authority etc. for each ward.
For example, this has set column 0 for the code, column 1 for the name,
and so the first ward has code ``E05002337`` and name ``Central``.
=============================
Extending the Model in Python
=============================

Adding a disease stage
----------------------

Continuing in ipython/jupyter from the last session, we will now extend the
disease to include an additional, less-infectious, semi-recovering stage,
which will come after I, and be called IR. We do this by inserting a new
stage, named "IR", at index 2, with ``beta`` value 0.2, and ``progress``
value 0.1

.. code-block:: python

   >>> lurgy.insert(2, name="IR", beta=0.2, progress=0.1)
   >>> print(lurgy)

   * Disease: lurgy
   * stage: ['E', 'I', 'IR', 'R']
   * mapping: ['E', 'I', 'IR', 'R']
   * beta: [0, 0.8, 0.2, 0]
   * progress: [0.25, 0.25, 0.1, 0]
   * too_ill_to_move: [0, 0, 0.0, 0]
   * start_symptom: 2

.. note::

   MetaWards is a Python program, so the index is counted from 0.
   Index 0 is E, index 1 is I and (before this call), index 2 was R.
   Inserting at index 2 will insert IR between I and R

We can now run the model using :func:`metawards.run`. This time we will
set ``silent`` to ``TRUE`` so that it doesn't print so much output
to the screen.

.. code-block:: python

   >>> results = mw.run(model=home, disease=lurgy,
                        additional=100, silent=True)

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ INFO ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Writing output to directory ./output_s839le7f

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

.. note::

   All of the output is written to the (randomly) named output directory
   indicated, e.g. for me to output_s839le7f. The full log of the run
   is recorded in the file called ``console.log.bz2`` which is in
   this directory.

We can now process and plot the results similarly to before, e.g.

.. code-block:: python

   >>> df = pd.read_csv(results)
   >>> df.plot.line(x="day", y=["S","E","I","IR","R"])

Repeating a run
---------------

MetaWards model runs are stochastic, meaning that they use random numbers.
While each individual run is reproducible (given the same random number
seed and number of processor threads), it is best to run multiple runs
so that you can look at averages.

You can perform multiple runs using the ``repeats`` argument, e.g.
to perform four runs, you should type;

.. code-block:: python

   >>> results = mw.run(model=home, disease=lurgy,
                        additional=100, silent=True, repeats=4)

If you look at the results, you will that there is a *repeat* column,
which indexes each run with a repeat number, e.g.

.. code-block:: python

   >>> df = pd.read_csv(results)
   >>> print(df)

        fingerprint  repeat  day        date      S   E   I  IR     R  IW
    0        REPEAT       1    0  2020-07-22  10000   0   0   0     0   0
    1        REPEAT       1    1  2020-07-23   9900  88  12   0     0   1
    2        REPEAT       1    2  2020-07-24   9890  74  30   6     0   1
    3        REPEAT       1    3  2020-07-25   9864  83  40  13     0   1
    4        REPEAT       1    4  2020-07-26   9835  95  48  20     2   1
    ..          ...     ...  ...         ...    ...  ..  ..  ..   ...  ..
    553      REPEAT       4  148  2020-12-17     63   0   0   1  9936   0
    554      REPEAT       4  149  2020-12-18     63   0   0   1  9936   0
    555      REPEAT       4  150  2020-12-19     63   0   0   1  9936   0
    556      REPEAT       4  151  2020-12-20     63   0   0   1  9936   0
    557      REPEAT       4  152  2020-12-21     63   0   0   0  9937   0

We can group by repeat and plot using;

.. code-block:: python

   >>> import matplotlib.pyplot as plt
   >>> df.groupby("repeat").plot.line(
        x="day", y=["S","E","I","IR","R"], ax=plt.gca())

You should get a result that looks something like this;

.. image:: ../images/py02.jpg
   :alt: Plot of the outbreak with a long recovery stage

.. note::

   With a bit more pandas you can make this plot a lot prettier ;-)

From this you can see the build-up of individuals in the yellow long
recovery (IR) stage.

Adding more wards
-----------------

Next, we will extend the model by adding more wards. We will model *home*,
*work* and *school*, so let's now add the *work* and *school* wards.

.. code-block:: python

   >>> work = mw.Ward("work")
   >>> school = mw.Ward("school")

We will now add some *workers* who will make daily, predictable movements
from *home* to *work* or *school*.

.. code-block:: python

   >>> home.add_workers(7500, destination=work)
   >>> home.add_workers(5000, destination=school)

.. note::

   The term *worker* is very broad in MetaWards. It means any individual
   that make regular, predictable movements each day. In this case, it
   refers to workers, teachers and students.

Next we need to combine these individual :class:`~metawards.Ward` objects
into a single :class:`~metawards.Wards` that represents the entire network.

.. code-block:: python

   >>> network = mw.Wards()
   >>> network.add(home)
   >>> network.add(work)
   >>> network.add(school)

Running the model
-----------------

We can now run the model. In this case, we want to seed the infection in
the *home* ward, so we need to pass this name into the ``additional``
parameter.

.. code-block:: python

   >>> results = mw.run(disease=lurgy, model=network,
                        additional="1, 100, home")

.. note::

   The format is **day number** (in this case seed on day 1), then
   **number to seed** (seeding 100 infections), then
   **ward name or number** (in this case, home)

You will see a lot of output. MetaWards does print a table to confirm
the seeding, e.g.

::

    ┏━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━┓
    ┃ Day ┃ Demographic ┃                     Ward                     ┃  Number   ┃
    ┃     ┃             ┃                                              ┃  seeded   ┃
    ┡━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━┩
    │  1  │    None     │ 1 : WardInfo(name='home', alternate_names=,  │    100    │
    │     │             │   code='', alternate_codes=, authority='',   │           │
    │     │             │        authority_code='', region='',         │           │
    │     │             │               region_code='')                │           │
    └─────┴─────────────┴──────────────────────────────────────────────┴───────────┘

The results can be processed and visualised as before, e.g.

.. code-block:: python

   >>> df = pd.read_csv(results)
   >>> df.plot.line(x="day", y=["S","E","I","IR","R"])

Complete code
-------------

The complete Python code for this part of the getting started guide is
re-copied below (this continues from the code in the last part);

.. code-block:: python

   # extend the disease model
   lurgy.insert(2, name="IR", beta=0.2, progress=0.1)

   # run the model
   results = mw.run(model=home, disease=lurgy,
                    additional=100, silent=True)

   # load and graph the results
   df = pd.read_csv(results)
   df.plot.line(x="day", y=["S","E","I","IR","R"])

   # run multiple repeats
   results = mw.run(model=home, disease=lurgy,
                    additional=100, silent=True, repeats=4)

   # load and graph the results
   df = pd.read_csv(results)
   import matplotlib.pyplot as plt
   df.groupby("repeat").plot.line(
        x="day", y=["S","E","I","IR","R"], ax=plt.gca())

   # add more wards
   work = mw.Ward("work")
   school = mw.Ward("school")
   home.add_workers(7500, destination=work)
   home.add_workers(5000, destination=school)

   # build the network
   network = mw.Wards()
   network.add(home)
   network.add(work)
   network.add(school)

   # run the model
   results = mw.run(disease=lurgy, model=network,
                    additional="1, 100, home")

   # load the graph the results
   df = pd.read_csv(results)
   df.plot.line(x="day", y=["S","E","I","IR","R"])
=============================
Multiple pathways with Python
=============================

MetaWards uses :class:`metawards.Demographics` to model different groups
of individuals in different ways. Individuals can move between
different demographics, and different demographics can experience
different disease models and move within different networks. This
is very powerful, and enables MetaWards to model multiple pathways
for different individuals.

This is explored in more depth in the :doc:`tutorial <../tutorial/index>`.
For this quick start guide, we will create three demographics;

* ``students`` : experience a mild version of the lurgy and travel to school each day
* ``teachers`` : experience the normal version of the lurgy, and also travel to school each day
* ``default`` : experience the normal version of the lurgy and either travel to work each day or stay home and play

Creating a mild disease
-----------------------

We must create a milder version of the lurgy that is written to the file
``mild_lurgy.json.bz2``;

::

  {
    "name": "mild_lurgy",
    "stage": ["E", "I", "R"],
    "beta": [0.0, 0.2, 0.0],
    "progress": [0.25, 0.5, 0.0]
  }

Creating the networks
---------------------

We now need to create the three networks for the three demographics.
We will start with the students, who will move between home and school.
This will be saved to ``students.json``.

::

  [
    {
        "id": 1,
        "info": {
        "name": "home"
        },
        "num_workers": 3000,
        "num_players": 0,
        "workers": {
        "destination": [
            2
        ],
        "population": [
            3000
        ]
        }
    },
    {
        "id": 2,
        "info": {
        "name": "school"
        },
        "num_workers": 0,
        "num_players": 0
    }
  ]

We will next do the same for the teachers, who will also move between
home and school (saving to ``teachers.json.bz2``).

::

  [
    {
        "id": 1,
        "info": {
        "name": "home"
        },
        "num_workers": 200,
        "num_players": 0,
        "workers": {
        "destination": [
            2
        ],
        "population": [
            200
        ]
        }
    },
    {
        "id": 2,
        "info": {
        "name": "school"
        },
        "num_workers": 0,
        "num_players": 0
    }
  ]

Next, we will create the default network. This will consist of some players
who stay at home, and workers who go to work.

::

  [
    {
        "id": 1,
        "info": {
        "name": "home"
        },
        "num_workers": 7000,
        "num_players": 10000,
        "workers": {
        "destination": [
            2
        ],
        "population": [
            7000
        ]
        }
    },
    {
        "id": 2,
        "info": {
        "name": "work"
        },
        "num_workers": 0,
        "num_players": 0
    }
  ]

Creating the demographics
-------------------------

Next, we create the demographics. We do this by creating
a file called ``network.json`` that contains data for each demographic that
specify the network and disease to use for each group.

::

  {
    "demographics": [
        "default",
        "teachers",
        "students"
    ],
    "diseases": [
        "lurgy.json",
        "lurgy.json",
        "mild_lurgy.json"
    ],
    "networks": [
        "default.json",
        "teachers.json",
        "students.json"
    ]
  }

.. note::

   Like before, it is easier to write these json file using the Python
   or R APIs. Any small errors in the file can cause difficult-to-debug
   errors when running metawards.

Running the model
-----------------

We can run the model by passing in the demographics. Note that we don't need
to specify the model as this is now fully specified in the demographics.

.. code-block:: bash

   metawards --disease lurgy.json --demographics demographics.json --additional "1, 100, home, default"

.. note::

   We have added ``default`` to the additional seeding to specify that the
   initial infections will be in this demographic. This is needed as a current
   limitation of MetaWards is that you can only seed infections in players,
   and only the default demographic in this example has players.

You can then process and graph the results as before;

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2

When you do this, you will notice that the number of susceptibles falls
until it reaches a number above 3200. This is because we seeded the outbreak
in the ``default`` demographic. By default, demographics do not mix with
each other, and so the outbreak does not spread to the teachers or
students.

We can control the amount of mixing of demographics using the ``mixer``
argument. This specifies a mixing function to use. We will use
:func:`~metawards.mixers.mix_evenly`, which sets that all demographics will
mix evenly with each other.

.. code-block:: bash

   metawards --disease lurgy.json --demographics demographics.json --additional "1, 100, home, default" --mixer mix_evenly
   metawards-plot -i output/results.csv.bz2

Now you should see that the outbreak spreads through the entire population.

.. note::

   The ``trajectory.csv.bz2`` file in the output directory of the run
   contains the trajectory for each of the demographics in each
   disease state. You can load this to generate demographic graphs.

What's next?
------------

This was a quick start guide to show some of the capabilities of MetaWards.
To learn more, e.g. how to create custom iterators to model lockdowns,
how to write extractors to get more detailed information output,
how to write mixers for modelling shielding etc., or how to write movers
to model conditional branching, please do now follow the
:doc:`tutorial <../tutorial/index>`.
=========================
Getting Started in Python
=========================

Prerequisites
-------------

First, you need to start an interactive Python session using a Python
interpreter with which MetaWards has been installed. To make this easier,
MetaWards comes with two applications; ``metawards-python`` and
``metawards-jupyter``, that will start an interative python or jupyter
session using the python interpreter used by MetaWards.

To start this, e.g. running a jupyter notebook, type;

.. code-block:: bash

   metawards-jupyter notebook

This should open a Jupyter notebook session in your browser. In here,
click "New" to start a new Python3 notebook, and then type the
commands for MetaWards below in the Jupyter cells.

You will also need to import pandas, as we will use this for
analysing and plotting the results.

.. code-block:: python

   >>> import pandas as pd

.. note::

   You also need to have installed MetaWardsData. If you have
   not installed MetaWardsData then you need to install it by
   following :doc:`these instructions <../model_data>`.

Importing metawards
-------------------

First we need to import the :mod:`metawards` Python module. To do this
we just need to type;

.. code-block:: python

   >>> import metawards as mw

The module comes with lots of help documentation, so feel free to use
that and tab-completion to explore the module.

Creating the disease
--------------------

You should now be in a Jupyter notebook (or ipython session) and have
imported :mod:`metawards`.

To run a simulation you need to define the :class:`~metawards.Disease`
that you want to model. MetaWards implements a SEIR-style model, but
you have complete control to define as many (or few) stages as you wish.

First, we will create a disease, which we will call ``lurgy``, that
will consist of four stages: S, E, I and R. To do this, let's create
the disease;

.. code-block:: python

   >>> lurgy = mw.Disease(name="lurgy")

Next, we will add each stage. You don't define the "S" stage, as the model
starts with a set of susceptible individuals by default. Instead, we need
to add in the E, I and R stages.

First, lets add the latent ("E") stage. Latent individuals are not
infectious, and so we will set ``beta`` (the infectivity parameter) to 0.0.
Individuals will progress quickly through this stage, so we will set
``progress`` to 0.5, meaning that 50% of individuals move to
the next stage each day.

.. code-block:: python

   >>> lurgy.add("E", beta=0.0, progress=0.5)

Next we will add the infectious ("I") stage. This will have a high ``beta``
value (0.8), but a lower progress (0.25) as we will model this as a
disease with a long symptomatic period.

.. code-block:: python

   >>> lurgy.add("I", beta=0.8, progress=0.25)

Finally, we need to add the recovered ("R") stage. We don't need to set the
``beta`` or ``progress`` values, as MetaWards will automatically recognise
this as the recovered state, and will set ``beta`` to 0 and ``progress``
to 0 automatically.

.. code-block:: python

   >>> lurgy.add("R")

You can should print this disease to the screen to confirm that everything
has been correctly set.

.. code-block:: python

   >>> print(lurgy)

   * Disease: lurgy
   * stage: ['E', 'I', 'R']
   * mapping: ['E', 'I', 'R']
   * beta: [0, 0.8, 0]
   * progress: [0.5, 0.25, 0]
   * too_ill_to_move: [0, 0, 0]
   * start_symptom: 2

.. note::

   You can save this disease to a file using
   ``lurgy.to_json("lurgy.json.bz2")``, and then load it back
   using ``lurgy = metawards.Disease.from_json("lurgy.json.bz2")``

Creating the wards (network)
----------------------------

Next, you need to define the wards (network) that will contain the individuals
who will experience the model outbreak.

We will first start with a single ward, called home.

.. code-block:: python

   >>> home = mw.Ward(name="home")

MetaWards works by assigning individuals as either `workers` or `players`.
The difference is that `workers` make fixed (predictable) movements
between different wards each day, while `players` make random movements.
Since we have just a single ward, we will start by populating it
with 10,000 players.

.. code-block:: python

   >>> home.set_num_players(10000)
   >>> print(home)

   Ward( info=home, num_workers=0, num_players=10000 )

.. note::

   You can save this Ward to a file using
   ``home.to_json("home.json.bz2")``, and then load it back
   using ``home = metawards.Ward.from_json("home.json.bz2")``

Running the model
-----------------

Now we have a disease and a network, we can now model an outbreak. To do this,
we will use the :func:`metawards.run` function.

.. code-block:: python

   >>> results = metawards.run(model=home, disease=lurgy)

This will print a lot to the screen. The key lines are these;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 5 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0
    Ending on day 5

This shows the number of people in the different stages of the outbreak.
In this case, there was no infection seeded, and so the number of infections
remained zero.

Seeding the outbreak
--------------------

We need to seed the outbreak with some additional seeds. We do this using
the ``additional`` option. This can be very powerful (e.g. adding seeds
at different days, different wards etc.), but at its simplest, it is
just the number of initial infections on the first day in the first
ward. We will start with 100 initial infections;

.. code-block:: python

   >>> results = metawards.run(model=home, disease=lurgy, additional=100)

Now you get a lot more output, e.g. for me the outbreak runs for 75 days.

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 70 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 423  E: 0  I: 1  R: 9576  IW: 0  POPULATION: 10000
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 71 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 423  E: 0  I: 1  R: 9576  IW: 0  POPULATION: 10000
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 72 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 423  E: 0  I: 1  R: 9576  IW: 0  POPULATION: 10000
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 73 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 423  E: 0  I: 1  R: 9576  IW: 0  POPULATION: 10000
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 74 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 423  E: 0  I: 1  R: 9576  IW: 0  POPULATION: 10000
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 75 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 423  E: 0  I: 0  R: 9577  IW: 0  POPULATION: 10000
    Number of infections: 0
    Ending on day 75


Visualising the results
-----------------------

The output ``results`` contains the filename of a csv file that contains
the S, E, I and R data (amongst other things). You can load and plot this
using standard R commands, e.g.

.. code-block:: python

   >>> df = pd.read_csv(results)
   >>> print(df)

        fingerprint  repeat  day        date      S    E   I     R  IW  SCALE_UV
    0        REPEAT       1    0  2020-07-21  10000    0   0     0   0       1.0
    1        REPEAT       1    1  2020-07-22   9900   76  24     0   1       1.0
    2        REPEAT       1    2  2020-07-23   9878   79  39     4   1       1.0
    3        REPEAT       1    3  2020-07-24   9840   95  49    16   1       1.0
    4        REPEAT       1    4  2020-07-25   9800  111  59    30   1       1.0
    ..          ...     ...  ...         ...    ...  ...  ..   ...  ..       ...
    103      REPEAT       1  103  2020-11-01    511    0   1  9488   0       1.0
    104      REPEAT       1  104  2020-11-02    511    0   1  9488   0       1.0
    105      REPEAT       1  105  2020-11-03    511    0   1  9488   0       1.0
    106      REPEAT       1  106  2020-11-04    511    0   1  9488   0       1.0
    107      REPEAT       1  107  2020-11-05    511    0   0  9489   0       1.0

    [108 rows x 10 columns]

We can visualise the data using;

.. code-block:: python

   >>> df.plot.line(x="day", y=["S","E","I","R"])

The result should look something like this;

.. image:: ../images/py01.jpg
   :alt: Plot of the initial outbreak

Complete code
-------------

The complete Python code for this part of the getting started guide is
re-copied below;

.. code-block:: python

   # import the required modules
   import pandas as pd
   import metawards as mw

   # create the disease
   lurgy = mw.Disease(name="lurgy")
   lurgy.add("E", beta=0.0, progress=0.25)
   lurgy.add("I", beta=0.8, progress=0.25)
   lurgy.add("R")

   # create the wards network
   home = mw.Ward(name="home")
   home.set_num_players(10000)

   # run the model
   results = metawards.run(model=home, disease=lurgy, additional=100)

   # load and graph the results
   df = pd.read_csv(results)
   df.plot.line(x="day", y=["S","E","I","R"])

=============================
Multiple pathways with Python
=============================

MetaWards uses :class:`metawards.Demographics` to model different groups
of individuals in different ways. Individuals can move between
different demographics, and different demographics can experience
different disease models and move within different networks. This
is very powerful, and enables MetaWards to model multiple pathways
for different individuals.

This is explored in more depth in the :doc:`tutorial <../tutorial/index>`.
For this quick start guide, we will create three demographics;

* ``students`` : experience a mild version of the lurgy and travel to school each day
* ``teachers`` : experience the normal version of the lurgy, and also travel to school each day
* ``default`` : experience the normal version of the lurgy and either travel to work each day or stay home and play

Creating a mild disease
-----------------------

First, we need to save the current version of the lurgy to a file called
``lurgy.json.bz2``.

.. code-block:: python

   >>> lurgy.to_json("lurgy.json.bz2")

Next, we must create a milder version of the lurgy and save this to
``mild_lurgy.json.bz2`` using;

.. code-block:: python

   >>> mild_lurgy = mw.Disease("mild_lurgy")
   >>> mild_lurgy.add("E", progress=0.25, beta=0.0)
   >>> mild_lurgy.add("I", progress=0.5, beta=0.2)
   >>> mild_lurgy.add("R")
   >>> mild_lurgy.to_json("mild_lurgy.json.bz2")

Creating the networks
---------------------

We now need to create the three networks for the three demographics.
We will start with the students, who will move between home and school.
This will be saved to ``students.json.bz2``.

.. code-block:: python

   >>> home = mw.Ward("home")
   >>> school = mw.Ward("school")
   >>> home.add_workers(3000, destination=school)
   >>> students = mw.Wards()
   >>> students.add(home)
   >>> students.add(school)
   >>> students.to_json("students.json.bz2")

We will next do the same for the teachers, who will also move between
home and school (saving to ``teachers.json.bz2``).

.. code-block:: python

   >>> home = mw.Ward("home")
   >>> school = mw.Ward("school")
   >>> home.add_workers(200, destination=school)
   >>> teachers = mw.Wards()
   >>> teachers.add(home)
   >>> teachers.add(school)
   >>> teachers.to_json("teachers.json.bz2")

Next, we will create the default network. This will consist of some players
who stay at home, and workers who go to work.

.. code-block:: python

   >>> home = mw.Ward("home")
   >>> work = mw.Ward("work")
   >>> home.set_num_players(10000)
   >>> home.add_workers(7000, destination=work)
   >>> default = mw.Wards()
   >>> default.add(home)
   >>> default.add(work)
   >>> default.to_json("default.json.bz2")

Creating the demographics
-------------------------

Next, we create the demographics. We do this by creating
:class:`~metawards.Demographic` objects for each demographic that
specify the network and disease to use for each group. These are then
combined into a single :class:`~metawards.Demographics` object.

.. code-block:: python

   >>> students = mw.Demographic("students",
                                 disease="mild_lurgy.json.bz2",
                                 network="students.json.bz2")
   >>> teachers = mw.Demographic("teachers",
                                 disease="lurgy.json.bz2",
                                 network="teachers.json.bz2")
   >>> default = mw.Demographic("default",
                                disease="lurgy.json.bz2",
                                network="default.json.bz2")
   >>> demographics = mw.Demographics()
   >>> demographics.add(default)
   >>> demographics.add(teachers)
   >>> demographics.add(students)
   >>> print(demographics)

   [
     Demographic(name='default', work_ratio=0.0, play_ratio=0.0, disease=lurgy.json.bz2, network='default.json.bz2')
     Demographic(name='teachers', work_ratio=0.0, play_ratio=0.0, disease=lurgy.json.bz2, network='teachers.json.bz2')
     Demographic(name='students', work_ratio=0.0, play_ratio=0.0, disease=mild_lurgy.json.bz2, network='students.json.bz2')
   ]

Running the model
-----------------

We can run the model by passing in the demographics. Note that we don't need
to specify the model as this is now fully specified in the demographics.

.. code-block:: python

   >>> results = mw.run(disease=lurgy, demographics=demographics,
                        additional="1, 5, home, default", silent=True)

.. note::

   We have added ``default`` to the additional seeding to specify that the
   initial infections will be in this demographic. This is needed as a current
   limitation of MetaWards is that you can only seed infections in players,
   and only the default demographic in this example has players.

You can then process and graph the results as before;

.. code-block:: python

   >>> df = pd.read_csv(results)
   >>> df.plot.line(x="day", y=["S","E","I","IR","R"])

When you do this, you will notice that the number of susceptibles falls
until it reaches a number above 3200. This is because we seeded the outbreak
in the ``default`` demographic. By default, demographics do not mix with
each other, and so the outbreak does not spread to the teachers or
students.

We can control the amount of mixing of demographics using the ``mixer``
argument. This specifies a mixing function to use. We will use
:func:`~metawards.mixers.mix_evenly`, which sets that all demographics will
mix evenly with each other.

.. code-block:: python

   >>> results = mw.run(disease=lurgy, demographics=demographics,
                        additional="1, 5, home, default",
                        mixer="mix_evenly", silent=True)
   >>> df = pd.read_csv(results)
   >>> df.plot.line(x="day", y=["S","E","I","IR","R"])

Now you should see that the outbreak spreads through the entire population.

.. note::

   The ``trajectory.csv.bz2`` file in the output directory of the run
   contains the trajectory for each of the demographics in each
   disease state. You can load this to generate demographic graphs.

Complete code
-------------

The complete Python code for this part of the getting started guide is
re-copied below (this continues from the code in the last part);

.. code-block:: python

   # save the lurgy to disk
   lurgy.to_json("lurgy.json.bz2")

   # create a milder lurgy and save to disk
   mild_lurgy = mw.Disease("mild_lurgy")
   mild_lurgy.add("E", progress=0.25, beta=0.0)
   mild_lurgy.add("I", progress=0.5, beta=0.2)
   mild_lurgy.add("R")
   mild_lurgy.to_json("mild_lurgy.json.bz2")

   # create the students network
   home = mw.Ward("home")
   school = mw.Ward("school")
   home.add_workers(3000, destination=school)
   students = mw.Wards()
   students.add(home)
   students.add(school)
   students.to_json("students.json.bz2")

   # create the teachers network
   home = mw.Ward("home")
   school = mw.Ward("school")
   home.add_workers(200, destination=school)
   teachers = mw.Wards()
   teachers.add(home)
   teachers.add(school)
   teachers.to_json("teachers.json.bz2")

   # create the default network
   home = mw.Ward("home")
   work = mw.Ward("work")
   home.set_num_players(10000)
   home.add_workers(7000, destination=work)
   default = mw.Wards()
   default.add(home)
   default.add(work)
   default.to_json("default.json.bz2")

   # now create the demographics
   students = mw.Demographic("students",
                             disease="mild_lurgy.json.bz2",
                             network="students.json.bz2")
   teachers = mw.Demographic("teachers",
                             disease="lurgy.json.bz2",
                             network="teachers.json.bz2")
   default = mw.Demographic("default",
                            disease="lurgy.json.bz2",
                            network="default.json.bz2")
   demographics = mw.Demographics()
   demographics.add(default)
   demographics.add(teachers)
   demographics.add(students)

   # run the model
   results = mw.run(disease=lurgy, demographics=demographics,
                    additional="1, 5, home, default",
                    mixer="mix_evenly", silent=True)

   # graph the results
   df = pd.read_csv(results)
   df.plot.line(x="day", y=["S","E","I","IR","R"])

What's next?
------------

This was a quick start guide to show some of the capabilities of MetaWards.
To learn more, e.g. how to create custom iterators to model lockdowns,
how to write extractors to get more detailed information output,
how to write mixers for modelling shielding etc., or how to write movers
to model conditional branching, please do now follow the
:doc:`tutorial <../tutorial/index>`.
========================
Extending the Model in R
========================

Adding a disease stage
----------------------

Continuing in R from the last session, we will now extend the disease
to include an additional, less-infectious, semi-recovering stage, which
will come after I, and be called IR. We do this by inserting a new
stage, named "IR", at index 2, with ``beta`` value 0.2, and ``progress``
value 0.1

.. code-block:: R

   > lurgy$insert(2, name="IR", beta=0.2, progress=0.1)
   > print(lurgy)

   * Disease: lurgy
   * stage: ['E', 'I', 'IR', 'R']
   * mapping: ['E', 'I', 'IR', 'R']
   * beta: [0.0, 0.8, 0.2, 0.0]
   * progress: [0.25, 0.25, 0.1, 0.0]
   * too_ill_to_move: [0.0, 0.0, 0.0, 0.0]
   * start_symptom: 2

.. note::

   MetaWards is a Python program, so the index is counted from 0.
   Index 0 is E, index 1 is I and (before this call), index 2 was R.
   Inserting at index 2 will insert IR between I and R

We can now run the model using :func:`metawards.run`. This time we will
set ``silent`` to ``TRUE`` so that it doesn't print so much output
to the screen.

.. code-block:: R

   > results <- metawards$run(model=home, disease=lurgy,
                              additional=100, silent=TRUE)

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ INFO ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Writing output to directory ./output_n81uzd7l

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

.. note::

   All of the output is written to the (randomly) named output directory
   indicated, e.g. for me to output_n81uzd7l. The full log of the run
   is recorded in the file called ``console.log.bz2`` which is in
   this directory.

We can now process and plot the results identically to before, e.g.

.. code-block:: R

   > results <- read.csv(results)
   > results <- results %>%
        pivot_longer(c("S", "E", "I", "R"),
        names_to = "stage", values_to = "count")
   > ggplot(data = results,
            mapping = aes(x=day, y=count, color=stage)) + geom_line()

Repeating a run
---------------

MetaWards model runs are stochastic, meaning that they use random numbers.
While each individual run is reproducible (given the same random number
seed and number of processor threads), it is best to run multiple runs
so that you can look at averages.

You can perform multiple runs using the ``repeats`` argument, e.g.
to perform four runs, you should type;

.. code-block:: R

   > results <- metawards$run(model=home, disease=lurgy,
                              additional=100, silent=TRUE, repeats=4)

If you look at the results, you will that there is a *repeat* column,
which indexes each run with a repeat number, e.g.

.. code-block:: R

   > results <- read.csv(results)
   > print(results)

       fingerprint repeat. day       date     S    E    I   IR    R IW SCALE_UV
    1       REPEAT       1   0 2020-07-21 10000    0    0    0    0  0        1
    2       REPEAT       1   1 2020-07-22  9900   82   18    0    0  1        1
    3       REPEAT       1   2 2020-07-23  9887   81   27    5    0  1        1
    4       REPEAT       1   3 2020-07-24  9869   77   44    9    1  1        1
    5       REPEAT       1   4 2020-07-25  9826  102   50   20    2  1        1
    6       REPEAT       1   5 2020-07-26  9783  113   67   34    3  1        1
    7       REPEAT       1   6 2020-07-27  9724  149   73   48    6  1        1
    8       REPEAT       1   7 2020-07-28  9653  174   96   64   13  1        1
    9       REPEAT       1   8 2020-07-29  9573  209  118   80   20  1        1
    10      REPEAT       1   9 2020-07-30  9472  254  145   99   30  1        1

.. note::

   Because ``repeat`` is a keyword in R, the column is automatically renamed
   as ``repeat.``

We can pivot and graph these runs using;

.. code-block:: R

   > results <- results %>%
        pivot_longer(c("S", "E", "I", "IR", "R"),
        names_to = "stage", values_to = "count")
   > ggplot(data = results,
            mapping = aes(x=day, y=count, color=stage)) + geom_point()

.. note::

   We have used ``geom_point()`` rather than ``geom_line()`` as this better
   shows the different runs. With a bit more R you could adjust the
   point shape to match the repeat number.

You should get a result that looks something like this;

.. image:: ../images/r02.jpg
   :alt: Plot of the outbreak with a long recovery stage

From this you can see the build-up of individuals in the green long
recovery (IR) stage.

Adding more wards
-----------------

Next, we will extend the model by adding more wards. We will model *home*,
*work* and *school*, so let's now add the *work* and *school* wards.

.. code-block:: R

   > work <- metawards$Ward("work")
   > school <- metawards$Ward("school")

We will now add some *workers* who will make daily, predictable movements
from *home* to *work* or *school*.

.. code-block:: R

   > home$add_workers(7500, destination=work)
   > home$add_workers(5000, destination=school)

.. note::

   The term *worker* is very broad in MetaWards. It means any individual
   that make regular, predictable movements each day. In this case, it
   refers to workers, teachers and students.

Next we need to combine these individual :class:`~metawards.Ward` objects
into a single :class:`~metawards.Wards` that represents the entire network.

.. code-block:: R

   > network <- metawards$Wards()
   > network$add(home)
   > network$add(work)
   > network$add(school)

Running the model
-----------------

We can now run the model. In this case, we want to seed the infection in
the *home* ward, so we need to pass this name into the ``additional``
parameter.

.. code-block:: R

   > results <- metawards$run(disease=lurgy, model=network,
                              additional="1, 100, home")

.. note::

   The format is **day number** (in this case seed on day 1), then
   **number to seed** (seeding 100 infections), then
   **ward name or number** (in this case, home)

You will see a lot of output. MetaWards does print a table to confirm
the seeding, e.g.

::

    ┏━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━┓
    ┃ Day ┃ Demographic ┃                     Ward                     ┃  Number   ┃
    ┃     ┃             ┃                                              ┃  seeded   ┃
    ┡━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━┩
    │  1  │    None     │ 1 : WardInfo(name='home', alternate_names=,  │    100    │
    │     │             │   code='', alternate_codes=, authority='',   │           │
    │     │             │        authority_code='', region='',         │           │
    │     │             │               region_code='')                │           │
    └─────┴─────────────┴──────────────────────────────────────────────┴───────────┘

The results can be processed and visualised as before, e.g.

.. code-block:: R

   > results <- read.csv(results)
   > results <- results %>%
        pivot_longer(c("S", "E", "I", "IR", "R"),
        names_to = "stage", values_to = "count")
   > ggplot(data = results,
            mapping = aes(x=day, y=count, color=stage)) + geom_point()

Complete code
-------------

The complete R code for this part of the getting started guide is
re-copied below (this continues from the code in the last part);

.. code-block:: R

   # add the IR stage between the I and R stages
   lurgy$insert(2, name="IR", beta=0.2, progress=0.1)

   # create the network of home, work and school wards
   work <- metawards$Ward("work")
   school <- metawards$Ward("school")
   network <- metawards$Wards()

   home$add_workers(7500, destination=work)
   home$add_workers(5000, destination=school)

   network$add(home)
   network$add(work)
   network$add(school)

   # run the model using the updated disease and network
   results <- metawards$run(disease=lurgy, model=network,
                            additional="1, 100, home")

   # plot the resulting trajectory
   results <- read.csv(results)
   results <- results %>%
        pivot_longer(c("S", "E", "I", "IR", "R"),
        names_to = "stage", values_to = "count")
   ggplot(data = results,
          mapping = aes(x=day, y=count, color=stage)) + geom_point()
=================
Quick Start Guide
=================

.. note::

   Please make sure you have installed ``metawards``, e.g. by following these
   :doc:`installation instructions <../install>`. You can test this by
   typing ``metawards --version`` on the terminal/console. Note that this
   should show that MetaWardsData has been found. If you have
   not installed MetaWardsData then you need to install it by
   following :doc:`these instructions <../model_data>`.

You have three choices for how you use MetaWards, and thus three choices
for how you will follow this quick start guide;

1. :doc:`R interface <01_R>`. You can choose to use MetaWards from directly
   within R (including interactively, e.g. within RStudio).

2. :doc:`Python interface <01_python>`. You can choose to use MetaWards from
   directly within Python (including interactively, e.g. within ipython or
   Jupyter notebook or lab).

3. :doc:`Command-line interface <01_console>`. You can choose to use MetaWards
   from the command line (terminal/console/command prompt).

There is no one right or wrong interface to use, so please feel free to use
the one (or many) that fit your needs.

This quick start should take about 30 minutes to follow. It is a precursor
to the :doc:`detailed tutorial <../tutorial/index>`. If you want to learn more,
please follow that tutorial, or read the
:doc:`detailed API documentation <../api/index>`. Please choose to follow
either the R, Python and console versions.

R / RStudio
-----------
.. toctree::
   :maxdepth: 2

   01_R
   02_R
   03_R

Python / Jupyter
----------------
.. toctree::
   :maxdepth: 2

   01_python
   02_python
   03_python

Console / Command line
----------------------
.. toctree::
   :maxdepth: 2

   01_console
   02_console
   03_console

Next steps
----------

Go to the :doc:`tutorial <../tutorial/index>` to learn more.
========================
Multiple pathways with R
========================

MetaWards uses :class:`metawards.Demographics` to model different groups
of individuals in different ways. Individuals can move between
different demographics, and different demographics can experience
different disease models and move within different networks. This
is very powerful, and enables MetaWards to model multiple pathways
for different individuals.

This is explored in more depth in the :doc:`tutorial <../tutorial/index>`.
For this quick start guide, we will create three demographics;

* ``students`` : experience a mild version of the lurgy and travel to school each day
* ``teachers`` : experience the normal version of the lurgy, and also travel to school each day
* ``default`` : experience the normal version of the lurgy and either travel to work each day or stay home and play

Creating a mild disease
-----------------------

First, we need to save the current version of the lurgy to a file called
``lurgy.json.bz2``.

.. code-block:: R

   > lurgy$to_json("lurgy.json.bz2")

Next, we must create a milder version of the lurgy and save this to
``mild_lurgy.json.bz2`` using;

.. code-block:: R

   > mild_lurgy <- metawards$Disease("mild_lurgy")
   > mild_lurgy$add("E", progress=0.25, beta=0.0)
   > mild_lurgy$add("I", progress=0.5, beta=0.2)
   > mild_lurgy$add("R")
   > mild_lurgy$to_json("mild_lurgy.json.bz2")

Creating the networks
---------------------

We now need to create the three networks for the three demographics.
We will start with the students, who will move between home and school.
This will be saved to ``students.json.bz2``.

.. code-block:: R

   > home <- metawards$Ward("home")
   > school <- metawards$Ward("school")
   > home$add_workers(3000, destination=school)
   > students <- metawards$Wards()
   > students$add(home)
   > students$add(school)
   > students$to_json("students.json.bz2")

We will next do the same for the teachers, who will also move between
home and school (saving to ``teachers.json.bz2``).

.. code-block:: R

   > home <- metawards$Ward("home")
   > school <- metawards$Ward("school")
   > home$add_workers(200, destination=school)
   > teachers <- metawards$Wards()
   > teachers$add(home)
   > teachers$add(school)
   > teachers$to_json("teachers.json.bz2")

Next, we will create the default network. This will consist of some players
who stay at home, and workers who go to work.

.. code-block:: R

   > home <- metawards$Ward("home")
   > work <- metawards$Ward("work")
   > home$set_num_players(10000)
   > home$add_workers(7000, destination=work)
   > default <- metawards$Wards()
   > default$add(home)
   > default$add(work)
   > default$to_json("default.json.bz2")

Creating the demographics
-------------------------

Next, we create the demographics. We do this by creating
:class:`~metawards.Demographic` objects for each demographic that
specify the network and disease to use for each group. These are then
combined into a single :class:`~metawards.Demographics` object.

.. code-block:: R

   > students <- metawards$Demographic("students",
                                       disease="mild_lurgy.json.bz2",
                                       network="students.json.bz2")
   > teachers <- metawards$Demographic("teachers",
                                       disease="lurgy.json.bz2",
                                       network="teachers.json.bz2")
   > default <- metawards$Demographic("default",
                                      disease="lurgy.json.bz2",
                                      network="default.json.bz2")
   > demographics <- metawards$Demographics()
   > demographics$add(default)
   > demographics$add(teachers)
   > demographics$add(students)
   > print(demographics)

   [
     Demographic(name='default', work_ratio=0.0, play_ratio=0.0, disease=lurgy.json.bz2, network='default.json.bz2')
     Demographic(name='teachers', work_ratio=0.0, play_ratio=0.0, disease=lurgy.json.bz2, network='teachers.json.bz2')
     Demographic(name='students', work_ratio=0.0, play_ratio=0.0, disease=mild_lurgy.json.bz2, network='students.json.bz2')
   ]

Running the model
-----------------

We can run the model by passing in the demographics. Note that we don't need
to specify the model as this is now fully specified in the demographics.

.. code-block:: R

   > results <- metawards$run(disease=lurgy, demographics=demographics,
                              additional="1, 5, home, default", silent=TRUE)

.. note::

   We have added ``default`` to the additional seeding to specify that the
   initial infections will be in this demographic. This is needed as a current
   limitation of MetaWards is that you can only seed infections in players,
   and only the default demographic in this example has players.

You can then process and graph the results as before;

.. code-block:: R

   > results <- read.csv(results)
   > results <- results %>%
        pivot_longer(c("S", "E", "I", "IR", "R"),
        names_to = "stage", values_to = "count")
   > ggplot(data = results,
            mapping = aes(x=day, y=count, color=stage)) + geom_point()

When you do this, you will notice that the number of susceptibles falls
until it reaches a number above 3200. This is because we seeded the outbreak
in the ``default`` demographic. By default, demographics do not mix with
each other, and so the outbreak does not spread to the teachers or
students.

We can control the amount of mixing of demographics using the ``mixer``
argument. This specifies a mixing function to use. We will use
:func:`~metawards.mixers.mix_evenly`, which sets that all demographics will
mix evenly with each other.

.. code-block:: R

   > results = metawards$run(disease=lurgy, demographics=demographics,
                             additional="1, 5, home, default",
                             mixer="mix_evenly", silent=TRUE)
   > results <- read.csv(results)
   > results <- results %>%
        pivot_longer(c("S", "E", "I", "IR", "R"),
        names_to = "stage", values_to = "count")
   > ggplot(data = results,
            mapping = aes(x=day, y=count, color=stage)) + geom_point()

Now you should see that the outbreak spreads through the entire population.

.. note::

   The ``trajectory.csv.bz2`` file in the output directory of the run
   contains the trajectory for each of the demographics in each
   disease state. You can load this to generate demographic graphs.

Complete code
-------------

The complete R code for this part of the getting started guide is
re-copied below (this continues from the code in the last part);

.. code-block:: R

   # save the lurgy to disk
   lurgy$to_json("lurgy.json.bz2")

   # create a milder lurgy and save to disk
   mild_lurgy <- metawards$Disease("mild_lurgy")
   mild_lurgy$add("E", progress=0.25, beta=0.0)
   mild_lurgy$add("I", progress=0.5, beta=0.2)
   mild_lurgy$add("R")
   mild_lurgy$to_json("mild_lurgy.json.bz2")

   # create the students network
   home <- metawards$Ward("home")
   school <- metawards$Ward("school")
   home$add_workers(3000, destination=school)
   students <- metawards$Wards()
   students$add(home)
   students$add(school)
   students$to_json("students.json.bz2")

   # create the teachers network
   home <- metawards$Ward("home")
   school <- metawards$Ward("school")
   home$add_workers(200, destination=school)
   teachers <- metawards$Wards()
   teachers$add(home)
   teachers$add(school)
   teachers$to_json("teachers.json.bz2")

   # create the default network
   home <- metawards$Ward("home")
   work <- metawards$Ward("work")
   home$set_num_players(10000)
   home$add_workers(7000, destination=work)
   default <- metawards$Wards()
   default$add(home)
   default$add(work)
   default$to_json("default.json.bz2")

   # now create the demographics
   students <- metawards$Demographic("students",
                                     disease="mild_lurgy.json.bz2",
                                     network="students.json.bz2")
   teachers <- metawards$Demographic("teachers",
                                     disease="lurgy.json.bz2",
                                     network="teachers.json.bz2")
   default <- metawards$Demographic("default",
                                    disease="lurgy.json.bz2",
                                    network="default.json.bz2")

   demographics <- metawards$Demographics()
   demographics$add(default)
   demographics$add(teachers)
   demographics$add(students)

   # run the model
   results = metawards$run(disease=lurgy, demographics=demographics,
                           additional="1, 5, home, default",
                           mixer="mix_evenly", silent=TRUE)

   # graph the results
   results <- read.csv(results)
   results <- results %>%
        pivot_longer(c("S", "E", "I", "IR", "R"),
        names_to = "stage", values_to = "count")
   ggplot(data = results,
          mapping = aes(x=day, y=count, color=stage)) + geom_point()

What's next?
------------

This was a quick start guide to show some of the capabilities of MetaWards.
To learn more, e.g. how to create custom iterators to model lockdowns,
how to write extractors to get more detailed information output,
how to write mixers for modelling shielding etc., or how to write movers
to model conditional branching, please do now follow the
:doc:`tutorial <../tutorial/index>`.
====================
Getting Started in R
====================

Prerequisites
-------------

Open R or RStudio.

Make sure you have installed MetaWards according to the
:doc:`R installation instructions <../install>`. Check this is
the case by typing;

.. code-block:: R

   > metawards::py_metawards_available()
   [1] TRUE

If you don't see ``TRUE`` returned, then double-check your installation.

.. note::

   You also need to have installed MetaWardsData. If you have
   not installed MetaWardsData then you need to install it by
   following :doc:`these instructions <../model_data>`.

Now make sure that you have installed the
`tidyverse <https://www.tidyverse.org>`__, e.g. via;

.. code-block:: R

   > install.packages("tidyverse")
   > library(tidyverse)

Then, finally, load the metawards R module;

.. code-block:: R

   > library(metawards)

Creating the disease
--------------------

You should now be in R (or RStudio) and have imported :mod:`metawards`.
To run a simulation you need to define the :class:`~metawards.Disease`
that you want to model. MetaWards implements a SEIR-style model, but
you have complete control to define as many (or few) stages as you wish.

First, we will create a disease, which we will call ``lurgy``, that
will consist of four stages: S, E, I and R. To do this, let's create
the disease;

.. code-block:: R

   > lurgy <- metawards$Disease(name="lurgy")

Next, we will add each stage. You don't define the "S" stage, as the model
starts with a set of susceptible individuals by default. Instead, we need
to add in the E, I and R stages.

First, lets add the latent ("E") stage. Latent individuals are not
infectious, and so we will set ``beta`` (the infectivity parameter) to 0.0.
Individuals will progress quickly through this stage, so we will set
``progress`` to 0.5, meaning that 50% of individuals move to
the next stage each day.

.. code-block:: R

   > lurgy$add("E", beta=0.0, progress=0.5)

Next we will add the infectious ("I") stage. This will have a high ``beta``
value (0.8), but a lower progress (0.25) as we will model this as a
disease with a long symptomatic period.

.. code-block:: R

   > lurgy$add("I", beta=0.8, progress=0.25)

Finally, we need to add the recovered ("R") stage. We don't need to set the
``beta`` or ``progress`` values, as MetaWards will automatically recognise
this as the recovered state, and will set ``beta`` to 0 and ``progress``
to 0 automatically.

.. code-block:: R

   > lurgy$add("R")

You can should print this disease to the screen to confirm that everything
has been correctly set.

.. code-block:: R

   > print(lurgy)

   * Disease: lurgy
   * stage: ['E', 'I', 'R']
   * mapping: ['E', 'I', 'R']
   * beta: [0, 0.8, 0]
   * progress: [0.5, 0.25, 0]
   * too_ill_to_move: [0, 0, 0]
   * start_symptom: 2

.. note::

   You can save this disease to a file using
   ``lurgy$to_json("lurgy.json.bz2")``, and then load it back
   using ``lurgy = metawards$Disease$from_json("lurgy.json.bz2")``

Creating the wards (network)
----------------------------

Next, you need to define the wards (network) that will contain the individuals
who will experience the model outbreak.

We will first start with a single ward, called home.

.. code-block:: R

   > home <- metawards$Ward(name="home")

MetaWards works by assigning individuals as either `workers` or `players`.
The difference is that `workers` make fixed (predictable) movements
between different wards each day, while `players` make random movements.
Since we have just a single ward, we will start by populating it
with 10,000 players.

.. code-block:: R

   > home$set_num_players(10000)
   > print(home)

   Ward( info=home, num_workers=0, num_players=10000 )

.. note::

   You can save this Ward to a file using
   ``home$to_json("home.json.bz2")``, and then load it back
   using ``home = metawards$Ward$from_json("home.json.bz2")``

Running the model
-----------------

Now we have a disease and a network, we can now model an outbreak. To do this,
we will use the :func:`metawards.run` function.

.. code-block:: R

   > results <- metawards$run(model=home, disease=lurgy)

This will print a lot to the screen. The key lines are these;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 5 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0
    Ending on day 5

This shows the number of people in the different stages of the outbreak.
In this case, there was no infection seeded, and so the number of infections
remained zero.

Seeding the outbreak
--------------------

We need to seed the outbreak with some additional seeds. We do this using
the ``additional`` option. This can be very powerful (e.g. adding seeds
at different days, different wards etc.), but at its simplest, it is
just the number of initial infections on the first day in the first
ward. We will start with 100 initial infections;

.. code-block:: R

   > results <- metawards$run(model=home, disease=lurgy, additional=100)

Now you get a lot more output, e.g. for me the outbreak runs for 75 days.

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 70 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 423  E: 0  I: 1  R: 9576  IW: 0  POPULATION: 10000
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 71 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 423  E: 0  I: 1  R: 9576  IW: 0  POPULATION: 10000
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 72 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 423  E: 0  I: 1  R: 9576  IW: 0  POPULATION: 10000
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 73 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 423  E: 0  I: 1  R: 9576  IW: 0  POPULATION: 10000
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 74 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 423  E: 0  I: 1  R: 9576  IW: 0  POPULATION: 10000
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 75 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 423  E: 0  I: 0  R: 9577  IW: 0  POPULATION: 10000
    Number of infections: 0
    Ending on day 75


Visualising the results
-----------------------

The output ``results`` contains the filename of a csv file that contains
the S, E, I and R data (amongst other things). You can load and plot this
using standard R commands, e.g.

.. code-block:: R

   > results <- read.csv(results)
   > print(results)
       fingerprint repeat. day       date     S    E    I    R IW SCALE_UV
    1       REPEAT       1   0 2020-07-20 10000    0    0    0  0        1
    2       REPEAT       1   1 2020-07-21  9900   57   43    0  1        1
    3       REPEAT       1   2 2020-07-22  9859   66   66    9  1        1
    4       REPEAT       1   3 2020-07-23  9807   86   82   25  1        1
    5       REPEAT       1   4 2020-07-24  9747  101  112   40  1        1
    6       REPEAT       1   5 2020-07-25  9654  140  130   76  1        1
    7       REPEAT       1   6 2020-07-26  9548  183  165  104  1        1
    8       REPEAT       1   7 2020-07-27  9433  215  203  149  1        1
    9       REPEAT       1   8 2020-07-28  9280  252  269  199  1        1
    10      REPEAT       1   9 2020-07-29  9082  318  341  259  1        1
    ...

To visualise the data we need to tidy it up so that we can group by S, E, I and R.

.. code-block:: R

   > results <- results %>%
        pivot_longer(c("S", "E", "I", "R"),
        names_to = "stage", values_to = "count")
   > print(results)
   # A tibble: 304 x 8
      fingerprint repeat.   day date          IW SCALE_UV stage count
      <chr>         <int> <int> <chr>      <int>    <dbl> <chr> <int>
    1 REPEAT            1     0 2020-07-20     0        1 S     10000
    2 REPEAT            1     0 2020-07-20     0        1 E         0
    3 REPEAT            1     0 2020-07-20     0        1 I         0
    4 REPEAT            1     0 2020-07-20     0        1 R         0
    5 REPEAT            1     1 2020-07-21     1        1 S      9900
    6 REPEAT            1     1 2020-07-21     1        1 E        57
    7 REPEAT            1     1 2020-07-21     1        1 I        43
    8 REPEAT            1     1 2020-07-21     1        1 R         0
    9 REPEAT            1     2 2020-07-22     1        1 S      9859
   10 REPEAT            1     2 2020-07-22     1        1 E        66
   # … with 294 more rows

You can graph S, E, I and R against day using;

.. code-block:: R

   > ggplot(data = results,
            mapping = aes(x=day, y=count, color=stage)) + geom_line()

The result should look something like this;

.. image:: ../images/r01.jpg
   :alt: Plot of the initial outbreak

Complete code
-------------

The complete R code for this part of the getting started guide is
re-copied below;

.. code-block:: R

   # Load the dependencies / libraries
   library(tidyverse)
   library(metawards)

   # Create the disease
   lurgy <- metawards$Disease(name="lurgy")
   lurgy$add("E", beta=0.0, progress=0.25)
   lurgy$add("I", beta=0.8, progress=0.25)
   lurgy$add("R")

   # Create the model network
   home <- metawards$Ward(name="home")
   home$set_num_players(10000)

   # Run the model
   results <- metawards$run(model=home, disease=lurgy, additional=100)

   # Read the tidy the results
   results <- read.csv(results)
   results <- results %>%
        pivot_longer(c("S", "E", "I", "R"),
        names_to = "stage", values_to = "count")

   # Graph the results
   ggplot(data = results,
          mapping = aes(x=day, y=count, color=stage)) + geom_line()
==================================
Extending the Model in the Console
==================================

Adding a disease stage
----------------------

Continuing in the terminal/console from the last session, we will now extend
the disease to include an additional, less-infectious, semi-recovering stage,
which will come after I, and be called IR. We do this by inserting a new
stage, named "IR", at index 2, with ``beta`` value 0.2, and ``progress``
value 0.1. Edit your ``lurgy.json`` file to read;

::

    {
      "name": "lurgy",
      "stage": ["E", "I", "IR", "R"],
      "beta": [0.0, 0.8, 0.2, 0.0],
      "progress": [0.5, 0.25, 0.1, 0.0]
    }

We can now run the model again;

.. code-block:: bash

   metawards --disease lurgy.json --model network.json --additional 100

We can now process and plot the results similarly to before, e.g.

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2

Repeating a run
---------------

MetaWards model runs are stochastic, meaning that they use random numbers.
While each individual run is reproducible (given the same random number
seed and number of processor threads), it is best to run multiple runs
so that you can look at averages.

You can perform multiple runs using the ``repeats`` argument, e.g.
to perform four runs, you should type;

.. code-block:: bash

   metawards --disease lurgy.json --model network.json --additional 100 --repeats 4

If you look at the results, you will that there is a *repeat* column,
which indexes each run with a repeat number, e.g.

::

    fingerprint,repeat,day,date,S,E,I,IR,R,IW,SCALE_UV
    REPEAT,1,0,2020-07-22,10000,0,0,0,0,0,1.0
    REPEAT,1,1,2020-07-23,9900,48,52,0,0,1,1.0
    REPEAT,1,2,2020-07-24,9863,60,70,7,0,1,1.0
    REPEAT,1,3,2020-07-25,9802,88,79,31,0,1,1.0
    REPEAT,1,4,2020-07-26,9727,111,110,49,3,1,1.0
    REPEAT,1,5,2020-07-27,9637,142,135,78,8,1,1.0

The ``metawards-plot`` command will automatically graph these repeats, e.g.

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2

You should get a resulting image (``output/overview.png``)
that looks something like this;

.. image:: ../images/console02.jpg
   :alt: Plot of the outbreak with a long recovery stage

.. note::

   You will need to load the data into Python pandas, R or Excel to
   visualise more data, as metawards-plot is limited in its
   visualisation capabilities.

Adding more wards
-----------------

Next, we will extend the model by adding more wards. We will model *home*,
*work* and *school*, so let's now add the *work* and *school* wards.
Edit your ``network.json`` file to read;

::

  [
    {
      "id": 1,
      "info": {
          "name": "home"
      },
      "num_workers": 12500,
      "num_players": 10000,
      "workers": {
        "destination": [2, 3],
      "population": [7500, 5000]
      }
    },
    {
      "id": 2,
      "info": {
          "name": "work"
      }
    },
    {
      "id": 3,
      "info": {
          "name": "school"
      }
    }
  ]

.. warning::

   Writing these network files by hand is quite difficult, and any small
   errors will cause issues. It is much better to construct this file
   using the Python or R API, and then save the file using the
   :func:`metawards.Wards.to_json` function.

.. note::

   The term *worker* is very broad in MetaWards. It means any individual
   that make regular, predictable movements each day. In this case, it
   refers to workers, teachers and students.

Running the model
-----------------

We can now run the model. In this case, we want to seed the infection in
the *home* ward, so we need to pass this name into the ``additional``
parameter.

.. code-block:: bash

   metawards --disease lurgy.json --model network.json --additional "1, 100, home"

.. note::

   The format is **day number** (in this case seed on day 1), then
   **number to seed** (seeding 100 infections), then
   **ward name or number** (in this case, home)

You will see a lot of output. MetaWards does print a table to confirm
the seeding, e.g.

::

    ┏━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━┓
    ┃ Day ┃ Demographic ┃                     Ward                     ┃  Number   ┃
    ┃     ┃             ┃                                              ┃  seeded   ┃
    ┡━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━┩
    │  1  │    None     │ 1 : WardInfo(name='home', alternate_names=,  │    100    │
    │     │             │   code='', alternate_codes=, authority='',   │           │
    │     │             │        authority_code='', region='',         │           │
    │     │             │               region_code='')                │           │
    └─────┴─────────────┴──────────────────────────────────────────────┴───────────┘

The results can be processed and visualised as before, e.g.

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2
==============================
Getting Started in the Console
==============================

Prerequisites
-------------

First, you need to open the terminal (or console or command prompt) and
type;

.. code-block:: bash

   metawards --version

You should see that MetaWards version information is printed to the screen.
If not, then you need to :doc:`install MetaWards <../install>`.

.. note::

   The output should show that MetaWardsData has been found. If you have
   not installed MetaWardsData then you need to install it by
   following :doc:`these instructions <../model_data>`.

You also need to be able to use a text editor, e.g. notepad, vim, emacs,
nano or pico. This quick start will use nano. Please use the editor
that you prefer.

Creating the disease
--------------------

You should now be at an open terminal.

To run a simulation you need to define the :class:`~metawards.Disease`
that you want to model. MetaWards implements a SEIR-style model, but
you have complete control to define as many (or few) stages as you wish.

First, we will create a disease, which we will call ``lurgy``, that
will consist of four stages: S, E, I and R. While you can use the
Python or R API to create this model in Python or R, you can also write
the required JSON data file directly at the console. Create a new
text file called ``lurgy.json``, e.g. by typing;

.. code-block:: bash

   nano lurgy.json

and type in the below;

::

  {
    "name": "lurgy",
    "stage": ["E", "I", "R"],
    "beta": [0.0, 0.8, 0.0],
    "progress": [0.5, 0.25, 0.0]
  }

This defines the three stages, "E", "I" and "R". You don't define the "S"
stage, as the model starts with a set of susceptible individuals by default.

The ``beta`` (infectivity) parameters are set such that individuals
are not infectious during the latent ("E") stage or recovered ("R") stage
(``beta`` equals 0), but are quite infectious in the "I" stage
(``beta`` equals 0.8).

The ``progress`` parameter is set so that individuals progress quickly
through the "E" stage (``progress`` equals 0.5, meaning that 50% of
individuals move to the next stage each day), while progress through
the "I" stage is slower (``progress`` equals 0.25). The ``progress``
value for the "R" stage must be 0, as once recovered, the individual
no longer moves through the model.

Creating the wards (network)
----------------------------

Next, you need to define the wards (network) that will contain the individuals
who will experience the model outbreak.

We will first start with a single ward, called home, in a file called
``network.json``.

::

  [
    {
      "info": {
          "name": "home"
      },
      "num_players": 10000
    }
  ]

MetaWards works by assigning individuals as either `workers` or `players`.
The difference is that `workers` make fixed (predictable) movements
between different wards each day, while `players` make random movements.
Since we have just a single ward, called "home", and we start by populating it
with 10,000 players.

Running the model
-----------------

Now we have a disease and a network, we can now model an outbreak. To do this,
we will run the ``metawards`` program directly.

.. code-block:: bash

   metawards --disease lurgy.json --model network.json

This will print a lot to the screen. The key lines are these;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 5 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 10000  E: 0  I: 0  R: 0  IW: 0  POPULATION: 10000
    Number of infections: 0
    Ending on day 5

This shows the number of people in the different stages of the outbreak.
In this case, there was no infection seeded, and so the number of infections
remained zero.

Seeding the outbreak
--------------------

We need to seed the outbreak with some additional seeds. We do this using
the ``additional`` option. This can be very powerful (e.g. adding seeds
at different days, different wards etc.), but at its simplest, it is
just the number of initial infections on the first day in the first
ward. We will start with 100 initial infections;

.. code-block:: bash

   metawards --disease lurgy.json --model network.json --additonal 100

.. note::

   MetaWards writes its output to a directory called ``output``. You can
   change this using the ``--output`` argument. By default, MetaWards will
   check before overwriting output. To remove this check, pass in the
   ``--force-overwrite-output`` option.

Now you get a lot more output, e.g. for me the outbreak runs for 71 days.

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 67 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 520  E: 1  I: 1  R: 9478  IW: 1  POPULATION: 10000
    Number of infections: 2

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 68 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 520  E: 0  I: 1  R: 9479  IW: 0  POPULATION: 10000
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 69 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 520  E: 0  I: 1  R: 9479  IW: 0  POPULATION: 10000
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 70 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 520  E: 0  I: 1  R: 9479  IW: 0  POPULATION: 10000
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 71 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 520  E: 0  I: 0  R: 9480  IW: 0  POPULATION: 10000
    Number of infections: 0
    Ending on day 71


Visualising the results
-----------------------

The output is written to the ``output`` directory. In here is a
comma-separated file called
``results.csv.bz2`` (MetaWards automatically bzip2 compresses most files to
save space). This contains the full trajectory, e.g. reading this
via ``bunzip2 -kc output/results.csv.bz2`` should show something that
starts with;

::

    fingerprint,repeat,day,date,S,E,I,R,IW,SCALE_UV
    REPEAT,1,0,2020-07-22,10000,0,0,0,0,1.0
    REPEAT,1,1,2020-07-23,9900,45,55,0,1,1.0
    REPEAT,1,2,2020-07-24,9867,48,66,19,1,1.0
    REPEAT,1,3,2020-07-25,9818,71,76,35,1,1.0
    REPEAT,1,4,2020-07-26,9755,98,99,48,1,1.0
    REPEAT,1,5,2020-07-27,9685,112,129,74,1,1.0
    REPEAT,1,6,2020-07-28,9587,151,158,104,1,1.0
    REPEAT,1,7,2020-07-29,9461,213,185,141,1,1.0
    REPEAT,1,8,2020-07-30,9317,260,235,188,1,1.0
    REPEAT,1,9,2020-07-31,9130,300,326,244,1,1.0
    REPEAT,1,10,2020-08-01,8869,406,399,326,1,1.0

We can visualise the data by loading into Python (pandas), R or Excel.
MetaWards also comes with a quick plotting program called ``metawards-plot``.
Use this to visualise the results using;

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2

.. note::

   This program may prompt you to install additional Python modules, e.g.
   pandas and matplotlib

This should produce a resulting image (``output/overview.png``)
that looks something like this;

.. image:: ../images/console01.jpg
   :alt: Plot of the initial outbreak
===============================
Part 4 - Customising the output
===============================

In this fourth part of the tutorial you will:

* Understand more about the other files that are output by ``metawards``.
* Learn about `metawards.extractors` and
  `~metawards.extractors.extract_default`.
* Write a custom extractor that can output different information to a file
  during a *model run*.
* Learn how to extract geographical information and include this in
  the outputs.


.. toctree::
   :maxdepth: 2

   part04/01_extractors
   part04/02_default
   part04/03_regional
   part04/04_rates
======================================================
Part 7 - Different pathways for different demographics
======================================================

In this seventh part of the tutorial you will;

* Learn how to give different disease pathways to different demographics.
* Explore how this can be used to model different disease behaviours
  and progress for different groups
* See how this can be combined with move functions to model conditional
  pathways, e.g. if individuals are sent to hospital
* See how you can name the different stages of the disease, e.g. to
  support collecting statistics such as the number of hospital patients

.. toctree::
   :maxdepth: 2

   part07/01_pathways
   part07/02_scanning
   part07/03_conditional
   part07/04_icu
   part07/05_named_stages
   part07/06_hospital
========================================
Part 5 - Modelling multiple demographics
========================================

In this fifth part of the tutorial you will:

* Understand how different demographics can be modelled using ``metawards``
* Learn about `metawards.Demographics` and how to distribute a population
  between multiple demographics
* Learn about `metawards.mixers`
* Learn how to write mixers to model shielding

.. toctree::
   :maxdepth: 2

   part05/01_demographics
   part05/02_mixers
   part05/03_custom
   part05/04_contacts
=====================================
Part 9 - Advanced Moves and Scenarios
=====================================

In this ninth part of the tutorial you will;

* Discover :class:`~metawards.movers.MoveGenerator` and
  :func:`~metawards.movers.go_ward`, and learn how this can be
  used to get full control to define very intricate movements
  of individuals.
* Use :meth:`Disease.is_infected <metawards.Disease.is_infected>` to
  create the non-infected Vaccinated (V) stage, and use
  :func:`~metawards.movers.go_ward` to model vaccination, both nationally
  and on a ward-by-ward basis.
* Use :class:`~metawards.movers.MoveGenerator` to create advanced moves
  that model individuals losing immunity due to vaccination or
  recovery after a period of time.
* Use :class:`~metawards.movers.MoveGenerator` to model movement of
  individuals from their home to their university ward at the start
  of the university year. This can include movements between
  workers and players.
* Use :meth:`Ward.bg_foi <metawards.Ward.bg_foi>` to set a background
  force of infection in wards, and use this to model hyperthetical
  foreign holiday destinations with different background rates of
  infection. Use :func:`~metawards.movers.go_ward` to model individuals
  going on holiday to these destinations, and the effect of different
  numbers of days of quarantine on return.

.. toctree::
   :maxdepth: 2

   part09/01_move_ward
   part09/02_vaccinate
   part09/03_fading_immunity
   part09/04_university
   part09/05_holidays
============================================
Part 8 - Creating your own models / networks
============================================

In this eighth part of the tutorial you will;

* Learn how the model network of wards and metapopulation movements is
  defined.
* Explore how to create your own model networks using the Python or R
  interface.
* Learn how to run MetaWards directly from within a Python or R script.
* Build more complex networks with ward-local parameters

.. toctree::
   :maxdepth: 2

   part08/01_networks
   part08/02_custom_network
   part08/03_reticulate
   part08/04_using_custom
   part08/05_local_network
   part08/06_network_demographics
========
Tutorial
========

This tutorial will take you through most aspects of running, understanding,
analysing and customising metawards. The tutorial is designed to be run
through in order.

.. warning::

    This tutorial models a **completely fictional** disease called the
    `lurgy <https://en.wiktionary.org/wiki/lurgy>`__.
    Any similarity to any real diseases is unintended, and the results or
    outputs from this tutorial **should not be used** to infer or
    imply anything about any real outbreak.

.. toctree::
   :maxdepth: 2

   index_part01
   index_part02
   index_part03
   index_part04
   index_part05
   index_part06
   index_part07
   index_part08
   index_part09
=================================
Part 3 - Customising the outbreak
=================================

In this third part of the tutorial you will:

* Understand the purpose of `metawards.iterators.iterate_default`
* Learn how to use different iterators on different days, so that
  the program can model weekends
* Create a custom iterator to represent a lockdown
* Learn how to trigger a custom iterator to start at a specific date
* Create an adjustable variable so that you can model lockdowns with
  different levels of effectiveness
* Setup models in which lockdowns are ended, and then explore what
  happens next
* Create models in which lockdowns can be stopped and started in
  response to changing conditions

.. toctree::
   :maxdepth: 2

   part03/01_iterators
   part03/02_weekend
   part03/03_lockdown
   part03/04_dynamic
   part03/05_scanning
   part03/06_scan_lockdown
   part03/07_cutoff
   part03/08_local_lockdown
============================
Part 1 - Modelling outbreaks
============================

In this part of the tutorial you will:

* Install ``metawards`` and perform the first *model runs* that simulate
  the outbreak of a fictional disease known as the lurgy.
* Learn how to seed an outbreak in specific city.
* Learn how to run a reproducible calculation.
* Learn how to perform multiple *model runs* and analyse the output
  that is collated in the ``results.csv.bz2`` file.
* Discover ``metawards-plot`` and use that to make simple plots.

.. toctree::
   :maxdepth: 2

   part01/01_getting_started
   part01/02_repeating
   part01/03_plotting
===========================
Part 2 - Refining the model
===========================

In this second part of the tutorial you will:

* Learn how to change a disease model by adding in an extra "infectious
  but asymptomatic" stage to the lurgy.
* Learn how to perform parameter sweeps to investigate how disease
  parameters change the population trajectories during a *model run*.
* Understand random error by running a parameter sweep with many
  repeats.
* Learn how to run ``metawards`` on a HPC cluster
* Learn how to examine a hypothesis and conduct a large model run on a
  cluster to gather evidence.
* Learn how to change the underlying geographic model data

.. toctree::
   :maxdepth: 2

   part02/01_disease
   part02/02_adjustable
   part02/03_analysis
   part02/04_cluster
   part02/05_refining
   part02/06_model_data
====================================
Part 6 - Moving between demographics
====================================

In this sixth part of the tutorial you will:

* Understand how you can move individuals between demographics using
  :doc:`metawards.movers <../api/index_MetaWards_movers>`.
* Use ``go functions`` to model the movement of individuals as they
  go into isolation and then are released back to their normal life
* Learn how to write custom movers to investigate when individuals
  should isolate, and how many days they should stay quarantined.
* Learn how to conditionally send a fraction of individuals to
  another demographic, so you can investigate different levels of
  compliance with self-isolation orders, and different compliance
  rates for quarantine

.. toctree::
   :maxdepth: 2

   part06/01_movers
   part06/02_duration
   part06/03_scan_duration
   part06/04_compliance
   part06/05_quarantine
==========================================
How many days of self-isolation is needed?
==========================================

In the last page you saw that about 1000 individuals were still
infected after seven days of self-isolation, and were released
back into the community. This was just 1-1.5% of the total number of
infections, so it was unlikely to have a big impact. We can
investigate this impact by scanning the number of days of
self-isolation required. Edit your ``move_isolate.py`` to read;

.. code-block:: python

    from metawards import Population
    from metawards import Networks
    from metawards.movers import go_isolate, go_to
    from metawards.utils import Console

    def move_isolate(network: Networks, population: Population, **kwargs):
        user_params = network.params.user_params

        ndays = user_params["isolate_ndays"]
        isolate_stage = user_params["isolate_stage"]

        if ndays > 7:
            Console.error(f"move_isolate supports a maximum of 7 days of "
                        f"isolation, so {ndays} is too many!")
            raise ValueError("Too many days of isolation requested")
        elif ndays <= 0:
            # just send infected individuals straight to "released"
            # (this is the control)
            func = lambda **kwargs: go_isolate(
                                        go_from="home",
                                        go_to="released",
                                        self_isolate_stage=isolate_stage,
                                        **kwargs)
            return [func]

        day = population.day % ndays
        isolate = f"isolate_{day}"

        go_isolate_day = lambda **kwargs: go_isolate(
                                            go_from="home",
                                            go_to=isolate,
                                            self_isolate_stage=isolate_stage,
                                            **kwargs)

        go_released = lambda **kwargs: go_to(go_from=isolate,
                                            go_to="released",
                                            **kwargs)

        return [go_released, go_isolate_day]

This move function will now read the number of days to isolate from
the ``isolate_ndays`` user parameter. It will also move an individual
into self-isolation when they reach disease stage ``isolate_stage``.

There is a little error-catching, e.g. as we only have seven available
``isolate_N`` demographics, we can only model up to seven days of
self-isolation. If more than seven days is requested this calls
:meth:`Console.error <metawards.utils.Console.error>` to write an
error to the output, and raises a Python ``ValueError``.

To add a control, we've set that if ``isolate_ndays`` is zero, then
infected individuals are sent straight from ``home`` to ``released``.
This way we can compare runs that use self-isolation against a run
where self-isolation is not performed.

.. note::

    We could of course increase the number of days of self-isolation that
    could be modelled by editing ``demographics.json`` and
    ``mix_isolate.py`` and just increasing the number of ``isolate_N``
    demographics.

Next create a file called ``scan_isolate.dat`` and copy in;

::

    # Scan through self-isolation from 0 to 7 days,
    # with self-isolation starting from disease
    # stages 2-4

    .isolate_stage   .isolate_ndays
            2               0
            2               1
            2               2
            2               3
            2               4
            2               5
            2               6
            2               7

            3               0
            3               1
            3               2
            3               3
            3               4
            3               5
            3               6
            3               7

            4               0
            4               1
            4               2
            4               3
            4               4
            4               5
            4               6
            4               7

This file will scan ``isolate_stage`` from 2 to 4, while scanning
``isolate_ndays`` from 0 to 7.

You can run this job locally using the command;

.. code-block:: bash

   metawards -d lurgy4 -D demographics.json -a ExtraSeedsLondon.dat --mixer mix_isolate --mover move_isolate --extractor extract_none --nsteps 365 -i scan_isolate.dat

Given that it is good to repeat the runs several times, and there are a lot
of jobs, you may want to run this on a cluster. To do this, the PBS
job script could look like;

.. code-block:: bash

    #!/bin/bash
    #PBS -l walltime=12:00:00
    #PBS -l select=4:ncpus=64:mem=64GB
    # The above sets 4 nodes with 64 cores each

    source $HOME/envs/metawards/bin/activate

    # change into the directory from which this job was submitted
    cd $PBS_O_WORKDIR

    metawards -d lurgy4 -D demographics.json -a ExtraSeedsLondon.dat \
              --mixer mix_isolate --mover move_isolate --extractor extract_none \
              --nsteps 365 -i scan_isolate.dat --repeats 8 \
              --nthreads 16 --force-overwrite-output --no-spinner --theme simple

while the slurm job script would be;

.. code-block:: bash

    #!/bin/bash
    #SBATCH --time=01:00:00
    #SBATCH --ntasks=4
    #SBATCH --cpus-per-task=64
    # The above sets 4 nodes with 64 cores each

    source $HOME/envs/metawards/bin/activate

    metawards -d lurgy4 -D demographics.json -a ExtraSeedsLondon.dat \
              --mixer mix_isolate --mover move_isolate --extractor extract_none \
              --nsteps 365 -i scan_isolate.dat --repeats 8 \
              --nthreads 16 --force-overwrite-output --no-spinner --theme simple

.. note::
   Notice that we've set the output theme to ``simple`` using
   ``--theme simple`` and have switched off the progress spinner
   using ``--no-spinner`` as these are unnecessary when run in batch
   mode on a cluster. All of the output will be written to the
   ``output/console.log.bz2`` file.

Analysing the results
---------------------

Once the job has completed you can generate and animate the overview plots
using ``metawards-plot``, e.g. via

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2
   metawards-plot --animate -i output/overview_*.jpg

The resulting animation should look something like this;

.. image:: ../../images/tutorial_6_3_1.gif
   :alt: Animated overview of different self-isolation durations

It is clear from these plots that self-isolation can work *if* it is
started at an early stage, e.g. stage 2 or 3. Self-isolating for 5-7 days
significantly reduces the numbers infected by the lurgy after a year,
if started as soon as possible.

However, what is also clear is that self-isolation has little impact
if it starts too late, e.g. by stage 4 of the lurgy. In this case,
even using the full seven days of self-isolation, more than 40 M
individuals are infected with six months.

The reason is clear when looking at the demographic plots. These can
be generated using;

.. code-block:: bash

   metawards-plot -i output/*/trajectory.csv.bz2
   metawards-plot --animate -i output/*x001/demographics.jpg -o demographics.gif

Mine are shown here;

.. image:: ../../images/tutorial_6_3_2.gif
   :alt: Animated demographic plots for different self-isolation durations

When caught early, most of the infected individuals are in one of the
``isolate_N`` demographics. However, when starting isolation from stage 4,
the vast majority of infected individuals are in the ``home`` demographic.
The epidemic experienced exponential growth within the ``home`` demographic
before individuals self-isolated, meaning the impact of self-isolation
on the course of the outbreak was minimal.

Relections on results
---------------------

For the lurgy, these results indicate that there is a very short window
of time between an individual showing symptoms and electively going
into self-isolation, for which self-isolation is effective.

Self-isolation (and by extensions, testing, contact-tracing etc.) would
not work for the lurgy in this model if individuals waited until
they'd reached stage 4 (showing larger symptoms and being more infectious),
or if the time it took to confirm an infection or trace infected
individuals was longer than the time it took to progress from stage
3 to stage 4.

It has to be done immediately from stage 3, when individuals start to
notice symptoms (even though, at this stage, the symptoms are not
limiting their movement).

It should also be noted that, in this model, ``beta`` for the asymptomatic
infectious stage had already been much reduced to account for society
adopting some control measures.
=======================
Movers and go functions
=======================

Demographics in ``metawards`` are a powerful concept that enables the
modelling of a wide variety of different scenarios. Just as *work*
and *play* have very general meanings in ``metawards``, so to do
*demographics*. We use it to mean any group of individuals. It is fluid,
in the sense that an individual can move between different demographics
during a *model run*, with the constraint that they can only belong
to one demographic at a time.

Demographic for self-isolation
------------------------------

Individuals are moved between demographics during a model run using
:doc:`mover functions <../../api/index_MetaWards_movers>`. These are
plugins that return the ``go functions`` that are used to make individuals
go from one demographic to another.

This is best demonstrated by example. In this example we will use
demographics to model the effect of self-isolation or
quarantine during an outbreak.

First, create a new ``demographics.json`` file that contains;

::

    {
        "demographics" : ["home", "isolate"],
        "work_ratios"  : [ 1.0,      0.0   ],
        "play_ratios"  : [ 1.0,      0.0   ]
    }

This specifies two demographics:

1. ``home`` - this holds the entire population and represents individuals
   behaving "normally", e.g. continuing to *work* and *play*.
2. ``isolate`` - this currently has no members. We will use this demographic
   to represent individuals who are self-isolating or in quarantine, e.g.
   they will not contribute to the force of infection of any ward.

Moving individuals to isolation
-------------------------------

Next, create a custom :doc:`move function <../../api/index_MetaWards_movers>`
called ``move_isolate`` by creating a file called ``move_isolate.py``
and copying in the below;

.. code-block:: python

    from metawards.movers import go_isolate

    def move_isolate(**kwargs):
        func = lambda **kwargs: go_isolate(go_from="home",
                                           go_to="isolate",
                                           self_isolate_stage=2,
                                            **kwargs)

        return [func]

This defines a custom :doc:`move function <../../api/index_MetaWards_movers>`
called ``move_isolate``. This returns the
``go function`` :meth:`~metawards.movers.go_isolate` that is
provided in :mod:`metawards.movers`. This
:meth:`~metawards.movers.go_isolate` function scans through the
demographics idenfied by ``go_from`` to search for individuals who
are showing signs of infection, i.e. individuals in a disease stage
that is greater or equal to ``self_isolate_stage``.

:meth:`~metawards.movers.go_isolate` moves these infected individuals
from their existing demographic into the new demographic identified
by ``go_to``.

This go function has several parameters that must be set before it
can be returned by ``move_isolate``. We set these parameters by using
`lambda <https://chryswoods.com/parallel_python/lambda.html>`__ to create
a new anonymous go function where those arguments are bound to fixed
values.

.. note::
   Here is a `good explanation of lambda and argument binding <https://chryswoods.com/parallel_python/lambda.html>`__
   if you've never seen this before. In this case we have bound
   ``go_from`` to equal ``"home"``, ``go_to`` to equal ``"isolate"``,
   and ``self_isolate_stage`` to equal ``2`` . This means that these values
   will be used every time the ``go_isolate`` function returned
   from ``move_isolate`` is called.

Mixing without infection
------------------------

Next, create a :doc:`mixer <../../api/index_MetaWards_mixers>` in
``mix_isolate.py`` and copy in the below;

.. code-block:: python

    from metawards.mixers import merge_using_matrix

    def mix_isolate(network, **kwargs):

        matrix = [ [1.0, 0.0],
                   [0.0, 0.0] ]

        network.demographics.interaction_matrix = matrix

        return [merge_using_matrix]

This mixer specifies an interaction matrix where the only contribution
to the FOIs comes from the ``home`` demographic (``matrix[0][0] == 1``).
The ``isolate`` demographic makes no contribution to the FOI
(``matrix[0][1]``, ``matrix[1][0]`` and ``matrix[1][1]`` are all zero).

Running the model
-----------------

You can run the simulation by passing in your custom mover using the
``--mover`` command line argument, and your custom mixer using the
``--mixer`` command line argument. We will seed the infection using
``ExtraSeedsBrighton.dat`` and will use the parameters from ``lurgy3.json``
which you should copy into this directory. Run the job using;

.. code-block:: bash

   metawards -d lurgy3 -D demographics.json -a ExtraSeedsBrighton.dat --mover move_isolate --mixer mix_isolate

.. note::
   Note that we are using the ``lurgy3`` parameters that were
   :doc:`optimised earlier <../part02/05_refining>`. These include the
   long-lived asymptomatic but infectious stage 3 of the disease.

You should see a trajectory that looks something like this;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
          home  S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
       isolate  S:        0  E: 0  I: 0  R: 0  IW: 0  POPULATION:        0
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    seeding demographic 0 play_infections[0][2124] += 5
    S: 56082072  E: 5  I: 0  R: 0  IW: 0  POPULATION: 56082077
          home  S: 56082072  E: 5  I: 0  R: 0  IW: 0  POPULATION: 56082077
       isolate  S:        0  E: 0  I: 0  R: 0  IW: 0  POPULATION:        0
    Number of infections: 5

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082072  E: 0  I: 5  R: 0  IW: 0  POPULATION: 56082077
          home  S: 56082072  E: 0  I: 5  R: 0  IW: 0  POPULATION: 56082077
       isolate  S:        0  E: 0  I: 0  R: 0  IW: 0  POPULATION:        0
    Number of infections: 5

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082072  E: 0  I: 5  R: 0  IW: 0  POPULATION: 56082077
          home  S: 56082072  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082072
       isolate  S:        0  E: 0  I: 5  R: 0  IW: 0  POPULATION:        5
    Number of infections: 5

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082072  E: 0  I: 5  R: 0  IW: 0  POPULATION: 56082077
          home  S: 56082072  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082072
       isolate  S:        0  E: 0  I: 5  R: 0  IW: 0  POPULATION:        5
    Number of infections: 5
    ...
    ...
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 20 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082072  E: 0  I: 1  R: 4  IW: 0  POPULATION: 56082077
          home  S: 56082072  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082072
       isolate  S:        0  E: 0  I: 1  R: 4  IW: 0  POPULATION:        5
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 21 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082072  E: 0  I: 1  R: 4  IW: 0  POPULATION: 56082077
          home  S: 56082072  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082072
       isolate  S:        0  E: 0  I: 1  R: 4  IW: 0  POPULATION:        5
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 22 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082072  E: 0  I: 0  R: 5  IW: 0  POPULATION: 56082077
          home  S: 56082072  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082072
       isolate  S:        0  E: 0  I: 0  R: 5  IW: 0  POPULATION:        5
    Number of infections: 0
    Infection died ... Ending on day 22

The infection was seeded with five individuals on day 1. They had a
latent infection for a day (``E == 5``), before developing symptoms
on day 2 (``I == 5``). At the beginning of day 3 then were moved
into the ``isolate`` demographic, in which they were unable to
infect others, and so progressed through the disease until they had
all recovered by day 22.

The asymptomatic stage
----------------------

Self-isolation appeared to have worked well. However, we neglected
to account for the asymptomatic ``stage 3`` of the lurgy. We modelled
the disease such that symptoms only appeared in stage 3, but individuals
were infectious from stage 2. We need to update our ``move_isolate`` function
so that individuals only self-isolate at stage 3, when they realise
that they have symptoms. Edit ``move_isolate.py`` and change it to read;

.. code-block:: python

    from metawards.movers import go_isolate

    def move_isolate(**kwargs):
        func = lambda **kwargs: go_isolate(go_from="home",
                                           go_to="isolate",
                                           self_isolate_stage=3,
                                            **kwargs)

        return [func]

(we have just changed ``self_isolate_stage`` from ``2`` to ``3``).

Now run ``metawards`` again using;

.. code-block:: bash

   metawards -d lurgy3 -D demographics.json -a ExtraSeedsBrighton.dat --mover move_isolate --mixer mix_isolate

We now have a completely different outbreak. Asymptomatic (and thus not
self-isolating) individuals were able to spread the infection to others
before going into isolation. This spread was exponential, and so
the epidemic lasted for a long time, with the vast majority of Individuals
in the model being infected. For example;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 370 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 11394138  E: 0  I: 1  R: 44687938  IW: 0  POPULATION: 56082077
          home  S: 11394138  E: 0  I: 0  R:        0  IW: 0  POPULATION: 11394138
       isolate  S:        0  E: 0  I: 1  R: 44687938  IW: 0  POPULATION: 44687939
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 371 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 11394138  E: 0  I: 1  R: 44687938  IW: 0  POPULATION: 56082077
          home  S: 11394138  E: 0  I: 0  R:        0  IW: 0  POPULATION: 11394138
       isolate  S:        0  E: 0  I: 1  R: 44687938  IW: 0  POPULATION: 44687939
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 372 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 11394138  E: 0  I: 0  R: 44687939  IW: 0  POPULATION: 56082077
          home  S: 11394138  E: 0  I: 0  R:        0  IW: 0  POPULATION: 11394138
       isolate  S:        0  E: 0  I: 0  R: 44687939  IW: 0  POPULATION: 44687939
    Number of infections: 0
    Infection died ... Ending on day 373

Here, the outbreak lasted for 372 days, with ~45M infections.

Adjusting progress
------------------

The issue here is that the amount of time spent in the asymptomatic but
infectious stage was very long (``progress[3] == 0.2``) and the
infectiousness of asymptomatic individuals was very high
(``beta[3] == 0.4``). During a real outbreak it is likely that individuals
would take actions that would reduce the chance of infection even from
asymptomatic carriers, e.g. by generally being more wary of one another,
washing hands, wearing masks etc. To account for this, we should reduce
the value of ``beta[3]`` to a lower value, e.g. to ``0.2``. Copy
``lurgy3.json`` to ``lurgy4.json`` and update that to read;

::

    { "name"             : "The Lurgy",
      "version"          : "May 18th 2020",
      "author(s)"        : "Christopher Woods",
      "contact(s)"       : "christopher.woods@bristol.ac.uk",
      "reference(s)"     : "Completely ficticious disease - no references",
      "beta"             : [0.0, 0.0, 0.2, 0.5, 0.5, 0.0],
      "progress"         : [1.0, 1.0, 0.2, 0.5, 0.5, 0.0],
      "too_ill_to_move"  : [0.0, 0.0, 0.0, 0.5, 0.8, 1.0],
      "contrib_foi"      : [1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
    }

(we have just changed ``beta[3]`` from ``0.4`` to ``0.2``)

.. note::

    In a real outbreak you should scan the value of ``beta[3]`` to match
    against observations. We are not doing this now as the lurgy is a
    fictional disease.

Now run the model using the command;

.. code-block:: bash

    metawards -d lurgy4 -D demographics.json -a ExtraSeedsBrighton.dat --mover move_isolate --mixer mix_isolate --nsteps 365

.. note::
   We've switched to using ``lurgy4`` and have limited the run to modelling
   just a single year (365 days)

You should see output similar to;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 362 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 55673047  E: 4015  I: 35748  R: 369267  IW: 3073  POPULATION: 56082077
          home  S: 55673047  E: 4015  I: 23938  R:   4094  IW: 3073  POPULATION: 55705094
       isolate  S:        0  E:    0  I: 11810  R: 365173  IW:    0  POPULATION:   376983
    Number of infections: 39763

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 363 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 55668968  E: 4094  I: 35800  R: 373215  IW: 2990  POPULATION: 56082077
          home  S: 55668968  E: 4094  I: 23987  R:   4079  IW: 2990  POPULATION: 55701128
       isolate  S:        0  E:    0  I: 11813  R: 369136  IW:    0  POPULATION:   380949
    Number of infections: 39894

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 364 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 55664887  E: 4079  I: 35966  R: 377145  IW: 3028  POPULATION: 56082077
          home  S: 55664887  E: 4079  I: 24148  R:   4081  IW: 3028  POPULATION: 55697195
       isolate  S:        0  E:    0  I: 11818  R: 373064  IW:    0  POPULATION:   384882
    Number of infections: 40045
    Exiting model run early
    Infection died ... Ending on day 365

.. note::

  You may find that the outbreak dies out quite quickly. The number of
  infections is low at the start, and the action of low ``beta`` and
  the move to self-isolating does quench the outbreak during some runs.
  However, once it catches light, the outbreak will spread to approximately
  40,000 individuals within one year.

The spread of the infection was significantly reduced by the reduction
in ``beta[3]`` for the asymptomatic stage. This demonstrates how small changes
in ``beta``, e.g. caused by increased hand-washing, masks etc., can have
a big impact on the spread of the disease in the model.
=======================
Self-isolation duration
=======================

In the last section we saw how self-isolation and a population that
took steps to reduce transmissability of the virus could dramatically
reduce the spread of the disease. However, in that model infected
individuals were moved into self-isolation for the entire duration
of the outbreak. This is clearly unrealistic.

Using demographics to represent days
------------------------------------

Typical advice to someone who is self-isolating is that they should
self-isolate for a set number of days. We can model this by using
different self-isolation demographics to represent the different
days that individuals start their self-isolation. For example,
if self-isolation was for seven days, then we could have a
self-isolation demographic for each day of the week. Once a week
is up, then the individuals who are self-isolating in that
day-demographic are released and moved to the "released" demographic.
Newly infected individuals for that day are then moved into
the now-empty day-demographic.

To do this, create a new demographics file called ``demographics.json``
and copy in the below;

::

    {
      "demographics" : ["home", "released",
                        "isolate_0", "isolate_1", "isolate_2",
                        "isolate_3", "isolate_4", "isolate_5",
                        "isolate_6" ],
      "work_ratios"  : [ 1.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ],
      "play_ratios"  : [ 1.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    }

This creates the ``home`` demographic, plus one ``isolate`` demographic
for each day of the week. There is also a ``released`` demographic
that will be used to release individuals from self-isolation.

Moving daily
------------

We start with all individuals placed into the ``home`` demographic. We will
now write a custom move function that will move individuals into the
assigned ``isolate_N`` demographic for the day in which they develop
symptoms. This move function will move them into the ``released`` demographic
once they have spent seven days in self-isolation. To do this,
create a move function called ``move_isolate.py`` and copy in the below;

.. code-block:: python

    from metawards import Population
    from metawards.movers import go_isolate, go_to

    def move_isolate(population: Population, **kwargs):
        day = population.day % 7
        isolate = f"isolate_{day}"

        go_isolate_day = lambda **kwargs: go_isolate(
                                            go_from="home",
                                            go_to=isolate,
                                            self_isolate_stage=3,
                                            **kwargs)

        go_released = lambda **kwargs: go_to(go_from=isolate,
                                             go_to="released",
                                             **kwargs)

        return [go_released, go_isolate_day]

This function works out which ``isolate_N`` demographic to use based
on the day of the week (``population.day % 7`` returns a number from ``0``
to ``6``).

It then creates two ``go functions``. The first, ``go_isolate_day``
is a :meth:`~metawards.movers.go_isolate` that moves infected
individuals from ``home`` into the ``isolate_N`` demographic of that day.

The second, ``go_released`` calls :meth:`~metawards.movers.go_to` to
send all individuals who are in that ``isolate_N`` demographic
to the ``released`` demographic.

The ``move_isolate`` function returns ``go_released`` first, so that
everyone who ends their self-isolation leaves before ``go_isolate_day``
then sends in the new cohort of infected individuals.

Mixing home and released
------------------------

Next, create a ``mixing function`` that merges the FOIs of the ``home``
and ``released`` demographics evenly, while making sure that everyone
in the ``isolate_N`` demographics is isolated and does not contribute
to any FOI.

Do this by creating a mixing function called ``mix_isolate.py`` and
copying in the below;

.. code-block:: python

    from metawards import Networks
    from metawards.mixers import merge_using_matrix, InteractionMatrix

    def mix_isolate(network: Networks, **kwargs):
        matrix = InteractionMatrix.ones(n=2)
        matrix.resize(2 + 7, value=0.0)

        network.demographics.interaction_matrix = matrix

        return [merge_using_matrix]

.. note::

   Note that we are using :func:`~metawards.mixers.merge_using_matrix`.
   This may not be the right choice depending on how we want the
   population dynamics to mix, e.g.
   :func:`~metawards.mixers.merge_matrix_single_population` or
   :func:`~metawards.mixers.merge_matrix_multi_population` may
   be a better choice. :doc:`See here for more information <../part05/04_contacts>`.

Here we use :class:`metawards.mixers.InteractionMatrix` to simplify the
creation of the interation matrix. We first create a 2x2 matrix;

::

  [ [1, 1],
    [1, 1] ]

using :meth:`InteractionMatrix.ones(n=2) <metawards.mixers.InteractionMatrix.ones>`.
We then resize this to be a 9x9 matrix using
:meth:`InteractionMatrix.resize(2 + 7, value=0.0) <metawards.mixers.InteractionMatrix.resize>`,
where the new values are equal to
zero. You can double-check that this matrix is correct using, e.g.
ipython or a jupyter notebook;

.. code-block:: python

    >>> from metawards.mixers import InteractionMatrix
    >>> matrix = InteractionMatrix.ones(n=2)
    >>> print(matrix)
    | 1.000, 1.000 |
    | 1.000, 1.000 |
    >>> matrix.resize(2 + 7, value=0.0)
    >>> print(matrix)
    | 1.000, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
    | 1.000, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
    | 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
    | 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
    | 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
    | 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
    | 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
    | 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
    | 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |

.. note::
   You could have created this matrix manually, but that is error-prone.
   The :class:`~metawards.mixers.InteractionMatrix` class has lots
   of helper functions that are useful for setting interactions between
   different demographics.

With this mixer created, you can now run ``metawards`` using;

.. code-block:: bash

   metawards -d lurgy4 -D demographics.json -a ExtraSeedsLondon.dat --mixer mix_isolate --mover move_isolate --extractor extract_none --nsteps 365

.. note::
   We've limited the number of days to model to 365 (one year), as
   self-isolation significantly slows down the spread of the disease,
   and modelling more than a year is unhelpful. We've also here used
   the :func:`~metawards.extractors.extract_none` extractor to limit
   the amount of output. Outputting data can be a little slow when
   there are a large number of demographics.

What do you see? In some cases self-isolation will cause the outbreak to
quickly die out. However, for most runs, we see that the infectious
asymptomatic allows new infections to be seeded before the individual
develops symptoms and moves into self-isolation.

Of more interest, we also see that as the outbreak grows, about 1% of
infected individuals have not recovered after 7 days. They leave
self-isolation and are able to contribute to the *force of infection*
in the ``home`` demographic.

::

    ...
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 11 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082064  E: 0  I: 6  R: 7  IW: 2  POPULATION: 56082077
           home  S: 56082064  E: 0  I: 6  R: 2  IW: 2  POPULATION: 56082072
       released  S:        0  E: 0  I: 0  R: 0  IW: 0  POPULATION:        0
      isolate_0  S:        0  E: 0  I: 0  R: 1  IW: 0  POPULATION:        1
      isolate_1  S:        0  E: 0  I: 0  R: 1  IW: 0  POPULATION:        1
      isolate_2  S:        0  E: 0  I: 0  R: 0  IW: 0  POPULATION:        0
      isolate_3  S:        0  E: 0  I: 0  R: 1  IW: 0  POPULATION:        1
      isolate_4  S:        0  E: 0  I: 0  R: 0  IW: 0  POPULATION:        0
      isolate_5  S:        0  E: 0  I: 0  R: 2  IW: 0  POPULATION:        2
      isolate_6  S:        0  E: 0  I: 0  R: 0  IW: 0  POPULATION:        0
    Number of infections: 6
    ...
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 39 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56081874  E: 19  I: 95  R: 89  IW: 13  POPULATION: 56082077
           home  S: 56081874  E: 19  I: 71  R: 13  IW: 13  POPULATION: 56081977
       released  S:        0  E:  0  I:  0  R: 48  IW:  0  POPULATION:       48
      isolate_0  S:        0  E:  0  I:  2  R:  5  IW:  0  POPULATION:        7
      isolate_1  S:        0  E:  0  I:  0  R:  8  IW:  0  POPULATION:        8
      isolate_2  S:        0  E:  0  I:  3  R:  4  IW:  0  POPULATION:        7
      isolate_3  S:        0  E:  0  I:  5  R:  2  IW:  0  POPULATION:        7
      isolate_4  S:        0  E:  0  I: 13  R:  0  IW:  0  POPULATION:       13
      isolate_5  S:        0  E:  0  I:  0  R:  3  IW:  0  POPULATION:        3
      isolate_6  S:        0  E:  0  I:  1  R:  6  IW:  0  POPULATION:        7
    Number of infections: 114
    ...
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 139 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56064409  E: 566  I: 4253  R: 12849  IW: 524  POPULATION: 56082077
           home  S: 56064409  E: 566  I: 2938  R:   585  IW: 524  POPULATION: 56068498
       released  S:        0  E:   0  I:   28  R: 10553  IW:   0  POPULATION:    10581
      isolate_0  S:        0  E:   0  I:   23  R:   374  IW:   0  POPULATION:      397
      isolate_1  S:        0  E:   0  I:   42  R:   377  IW:   0  POPULATION:      419
      isolate_2  S:        0  E:   0  I:   85  R:   341  IW:   0  POPULATION:      426
      isolate_3  S:        0  E:   0  I:  121  R:   272  IW:   0  POPULATION:      393
      isolate_4  S:        0  E:   0  I:  218  R:   232  IW:   0  POPULATION:      450
      isolate_5  S:        0  E:   0  I:  364  R:   115  IW:   0  POPULATION:      479
      isolate_6  S:        0  E:   0  I:  434  R:     0  IW:   0  POPULATION:      434
    Number of infections: 4819
    ...
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 364 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 54653198  E: 14724  I: 126080  R: 1288075  IW: 6273  POPULATION: 56082077
           home  S: 54653198  E: 14724  I: 84809  R:   14613  IW: 6273  POPULATION: 54767344
       released  S:        0  E:     0  I:  1014  R: 1218000  IW:    0  POPULATION:  1219014
      isolate_0  S:        0  E:     0  I: 14010  R:       0  IW:    0  POPULATION:    14010
      isolate_1  S:        0  E:     0  I:   834  R:   12575  IW:    0  POPULATION:    13409
      isolate_2  S:        0  E:     0  I:  1490  R:   11998  IW:    0  POPULATION:    13488
      isolate_3  S:        0  E:     0  I:  2482  R:   11212  IW:    0  POPULATION:    13694
      isolate_4  S:        0  E:     0  I:  4150  R:    9327  IW:    0  POPULATION:    13477
      isolate_5  S:        0  E:     0  I:  6926  R:    6897  IW:    0  POPULATION:    13823
      isolate_6  S:        0  E:     0  I: 10365  R:    3453  IW:    0  POPULATION:    13818
    Number of infections: 140804

As you can see above, by day 364, there were 1014 infected individuals in
the ``released`` demographic, indicating that they have left self-isolation
too early.
========================
Staying under quarantine
========================

From the previous page, we saw that at least 90% of individuals needed to
comply with self-isolation to limit the outbreak.
However, that run assumed that everyone who entered self-isolation remained
quarantined for the full seven days. It is important that we also model
the impact of individuals breaking quarantine early.

go_early functions
------------------

We can represent a fraction of individuals
leaving quarantine early each day by passing the keyword
argument ``fraction`` to
:meth:`~metawards.movers.go_to`. Do this by updating your
``move_isolate.py`` to read;

.. code-block:: python

    from metawards import Population
    from metawards import Networks
    from metawards.movers import go_isolate, go_to

    def move_isolate(network: Networks, population: Population, **kwargs):
        user_params = network.params.user_params

        ndays = 7
        isolate_stage = 3
        compliance_fraction = 0.9

        # fraction who remain in isolation, counting from the longest to
        # shortest stay in isolation.
        remain = [0.9, 0.9, 0.95, 0.95, 1.00, 1.00]

        day = population.day % 7
        isolate = f"isolate_{day}"

        go_early = []

        # have to define this functions one-by-one and not in a loop
        # otherwise python will bind all functions to the value of i
        # of the last iteration of the loop
        go_early.append(lambda **kwargs: go_to(
                                go_from=f"isolate_{(day + 1) % 7}",
                                go_to="released",
                                fraction=(1.0 - remain[0]),
                                **kwargs))
        go_early.append(lambda **kwargs: go_to(
                                go_from=f"isolate_{(day + 2) % 7}",
                                go_to="released",
                                fraction=(1.0 - remain[1]),
                                **kwargs))
        go_early.append(lambda **kwargs: go_to(
                                go_from=f"isolate_{(day + 3) % 7}",
                                go_to="released",
                                fraction=(1.0 - remain[2]),
                                **kwargs))
        go_early.append(lambda **kwargs: go_to(
                                go_from=f"isolate_{(day + 4) % 7}",
                                go_to="released",
                                fraction=(1.0 - remain[3]),
                                **kwargs))
        go_early.append(lambda **kwargs: go_to(
                                go_from=f"isolate_{(day + 5) % 7}",
                                go_to="released",
                                fraction=(1.0 - remain[4]),
                                **kwargs))
        go_early.append(lambda **kwargs: go_to(
                                go_from=f"isolate_{(day + 6) % 7}",
                                go_to="released",
                                fraction=(1.0 - remain[5]),
                                **kwargs))

        go_isolate_day = lambda **kwargs: go_isolate(
                                            go_from="home",
                                            go_to=isolate,
                                            self_isolate_stage=isolate_stage,
                                            fraction=compliance_fraction,
                                            **kwargs)

        go_released = lambda **kwargs: go_to(go_from=isolate,
                                             go_to="released",
                                             **kwargs)

        return go_early + [go_released, go_isolate_day]


.. note::
   It would be nicer and less error-prone if we could create the
   ``go_early`` functions in a loop. However, this would not work
   because of the way that Python lambda functions bind their arguments.
   If we did this, the arguments from the last iteration of the loop
   would be used for all of the ``go_early`` functions, i.e. we would
   try to move individuals out of the same ``isolate_N`` demographic
   six times.

.. note::
   Note that we've set ``compliance`` to 0.9 based on the results of the
   last scan.

Here, we've created a new set of go functions called ``go_release_early``.
There is one for each ``isolate_N`` demographic *except* for the
demographic to which individuals will be moved on each day.

This ``go_release_early`` function moves a fraction of individuals from
the ``isolate_N`` demographic to ``released``, representing that fraction
breaking their quarantine early. This fraction is taken from the list
``remain``, which counts up from ``0.90`` to ``1.00``. The first value
(``0.90``) is the fraction for individuals that have been isolating the longest
(six days), and that will remain in isolation that day (i.e. 90% will remain,
while 10% will break quarantine early). The last value (``1.00``) is the
fraction for the individuals who only entered isolation the previous day,
i.e. everyone remains in isolation for at least one day. These ``go_early``
functions are then added before ``go_released`` and ``go_isolate_day``.

Now run ``metawards`` using your ``move_isolate.py`` via;

.. code-block:: bash

   metawards -d lurgy4 -D demographics.json -a ExtraSeedsLondon.dat --mixer mix_isolate --mover move_isolate --nsteps 365

You should see that the disease spreads, now both from individuals who
choose not to self-isolate, and now also from individuals who break
their quarantine early. Graphing the output via;

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2

gives an overview plot that should look something like;

.. image:: ../../images/tutorial_6_5_1.jpg
   :alt: Effect of compliance of remaining in quarantine

It is clear that individuals who break quarantine early contribute
significantly to further spread of the outbreak (> 14 million have
experienced infection after one year in this *model run* compared
to ~6 million if all individuals strictly observed the full seven
days of quarantine).

Scanning quarantine
-------------------

The next step is to scan through different compliance levels for
different numbers of days. To do this, update your
``move_isolate.py`` to read;

.. code-block:: python

    from metawards import Population
    from metawards import Networks
    from metawards.movers import go_isolate, go_to

    def move_isolate(network: Networks, population: Population, **kwargs):
        user_params = network.params.user_params

        ndays = 7
        isolate_stage = 3
        compliance_fraction = 0.9

        # fraction who remain in isolation, counting from the longest to
        # shortest stay in isolation.
        remain = user_params["remain"]

        day = population.day % 7
        isolate = f"isolate_{day}"

        go_early = []

        # have to define this functions one-by-one and not in a loop
        # otherwise python will bind all functions to the value of i
        # of the last iteration of the loop
        go_early.append(lambda **kwargs: go_to(
                                go_from=f"isolate_{(day + 1) % 7}",
                                go_to="released",
                                fraction=(1.0 - remain[0]),
                                **kwargs))
        go_early.append(lambda **kwargs: go_to(
                                go_from=f"isolate_{(day + 2) % 7}",
                                go_to="released",
                                fraction=(1.0 - remain[1]),
                                **kwargs))
        go_early.append(lambda **kwargs: go_to(
                                go_from=f"isolate_{(day + 3) % 7}",
                                go_to="released",
                                fraction=(1.0 - remain[2]),
                                **kwargs))
        go_early.append(lambda **kwargs: go_to(
                                go_from=f"isolate_{(day + 4) % 7}",
                                go_to="released",
                                fraction=(1.0 - remain[3]),
                                **kwargs))
        go_early.append(lambda **kwargs: go_to(
                                go_from=f"isolate_{(day + 5) % 7}",
                                go_to="released",
                                fraction=(1.0 - remain[4]),
                                **kwargs))
        go_early.append(lambda **kwargs: go_to(
                                go_from=f"isolate_{(day + 6) % 7}",
                                go_to="released",
                                fraction=(1.0 - remain[5]),
                                **kwargs))

        go_isolate_day = lambda **kwargs: go_isolate(
                                            go_from="home",
                                            go_to=isolate,
                                            self_isolate_stage=isolate_stage,
                                            fraction=compliance_fraction,
                                            **kwargs)

        go_released = lambda **kwargs: go_to(go_from=isolate,
                                            go_to="released",
                                            **kwargs)

        return go_early + [go_released, go_isolate_day]

The only change is that we now read the ``remain`` list from the
custom user variable called ``remain``.

Create a scan file called ``scan_remain.dat`` and copy in the below;

::

 .remain[0]  .remain[1]  .remain[2]  .remain[3]  .remain[4]  .remain[5]
    1.00        1.00        1.00        1.00        1.00        1.00

    0.95        1.00        1.00        1.00        1.00        1.00
    0.95        0.95        1.00        1.00        1.00        1.00
    0.95        0.95        0.95        1.00        1.00        1.00
    0.95        0.95        0.95        0.95        1.00        1.00
    0.95        0.95        0.95        0.95        0.95        1.00
    0.95        0.95        0.95        0.95        0.95        0.95

    0.90        0.95        0.95        0.95        0.95        0.95
    0.90        0.90        0.95        0.95        0.95        0.95
    0.90        0.90        0.90        0.95        0.95        0.95
    0.90        0.90        0.90        0.90        0.95        0.95
    0.90        0.90        0.90        0.90        0.90        0.95
    0.90        0.90        0.90        0.90        0.90        0.90

    0.85        0.90        0.90        0.90        0.90        0.90
    0.85        0.85        0.90        0.90        0.90        0.90
    0.85        0.85        0.85        0.90        0.90        0.90
    0.85        0.85        0.85        0.85        0.90        0.90
    0.85        0.85        0.85        0.85        0.85        0.90
    0.85        0.85        0.85        0.85        0.85        0.85

    0.80        0.85        0.85        0.85        0.85        0.85
    0.80        0.80        0.85        0.85        0.85        0.85
    0.80        0.80        0.80        0.85        0.85        0.85
    0.80        0.80        0.80        0.80        0.85        0.85
    0.80        0.80        0.80        0.80        0.80        0.85
    0.80        0.80        0.80        0.80        0.80        0.80

This will scan the percentage of individuals who should remain in quarantine
each day from ``1.00`` (100%) to ``0.80`` (80%). It does this in increments
of 0.05, starting from the longest period of isolation (7 days) and moving
that to the shortest period (less than 1 day).

There are a large number of jobs, and repeats are needed to properly
sample the outbreak. Here are job scripts to run this job on either a
PBS or Slurm cluster;

.. code-block:: bash

    #!/bin/bash
    #PBS -l walltime=12:00:00
    #PBS -l select=4:ncpus=64:mem=64GB
    # The above sets 4 nodes with 64 cores each

    source $HOME/envs/metawards/bin/activate

    # change into the directory from which this job was submitted
    cd $PBS_O_WORKDIR

    metawards -d lurgy4 -D demographics.json -a ExtraSeedsLondon.dat \
            --mixer mix_isolate --mover move_isolate --extractor extract_none \
            --nsteps 365 -i scan_remain.dat --repeats 8 \
            --nthreads 16 --force-overwrite-output --no-spinner --theme simple

.. code-block:: bash

    #!/bin/bash
    #SBATCH --time=01:00:00
    #SBATCH --ntasks=4
    #SBATCH --cpus-per-task=64
    # The above sets 4 nodes with 64 cores each

    source $HOME/envs/metawards/bin/activate

    metawards -d lurgy4 -D demographics.json -a ExtraSeedsLondon.dat \
            --mixer mix_isolate --mover move_isolate --extractor extract_none \
            --nsteps 365 -i scan_remain.dat --repeats 8 \
            --nthreads 16 --force-overwrite-output --no-spinner --theme simple

Analysis
--------

Once you have run the jobs, you can generate the animation of the overview
plots using;

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2
   metawards-plot --animate output/overview*

You should get an animation that looks something like this;


.. image:: ../../images/tutorial_6_5_2.gif
   :alt: Scanning compliance of remaining in quarantine

Reflections on results
----------------------

Again, it is clear that self-isolation only works well if there is very
high compliance with remaining in quarantine for the full seven days.
As 5% of people break quarantine after 2-3 days, the size of the outbreak
grows from ~6 million to ~15 million. This grows to >20 million if 10%
break quarantine.

All of these results suggest that, for the lurgy, that self-isolation and
quarantine are only effective strategies if;

1. they start as quickly as possible once an individual is infectious
   (which may be difficult due to the early stage of the lurgy being largely
   asymptomatic - track and trace systems would likely need to be used
   and be extremely fast and effective),
2. are observed by the vast majority (>90%) of infected individuals, and
3. >95% of individuals remain in quarantine for the full seven days.

As noted at the start of this tutorial, the lurgy is not a serious disease,
so would not warrant the level of social control needed to achieve these
three requirements (and, indeed, a population that does not perceive the
lurgy to be serious would be unlikely to comply to the level needed).

This should not be surprising given where each *model run* starts - just
five infected individuals in one ward. Any disease that is so contagious that
an infection in such a small group of individuals grows quickly into a
widespread outbreak will always be very difficult to control.
======================
Conditional Compliance
======================

In the last section we saw that self-isolation will only work if
individuals self-isolate at an early stage of the lurgy, and for 5-7 days.

However, that model and conclusion is flawed. The model assumed that
everyone who had reached the self-isolation stage would elect to
self-isolate, and would remain their for the full required duration.

Unfortunately, this is not realistic. Many individuals would not
self-isolate, e.g. because they don't feel that ill, don't recognise
the symptoms, have financial or social pressures to keep going out
etc. Equally, many may break quarantine early, if they feel better or
more mobile.

To model this, we need to conditionally send only a percentage of
individuals into self-isolation. And then, for each day quarantine,
we need to release an increasing percentage back to their normal
behaviour.

Moving a fraction of individuals
--------------------------------

Both the :meth:`~metawards.movers.go_isolate` and
:meth:`~metawards.movers.go_to` go functions accept the ``fraction``
keyword argument. This argument sets the percentage (or fraction)
of the population who should move. By default this is ``1.0`` (representing
100%). This fraction is used to sample a random number of individuals
according to a binomial distribution.

We can adjust this value to examine impact of reduced self-isolation
compliance on the outbreak.

To do this, update your ``move_isolate.py`` move function to;

.. code-block:: python

    from metawards import Population
    from metawards import Networks
    from metawards.movers import go_isolate, go_to

    def move_isolate(network: Networks, population: Population, **kwargs):
        ndays = 7
        isolate_stage = 3
        compliance_fraction = 0.5

        day = population.day % ndays
        isolate = f"isolate_{day}"

        go_isolate_day = lambda **kwargs: go_isolate(
                                            go_from="home",
                                            go_to=isolate,
                                            self_isolate_stage=isolate_stage,
                                            fraction=compliance_fraction,
                                            **kwargs)

        go_released = lambda **kwargs: go_to(go_from=isolate,
                                            go_to="released",
                                            **kwargs)

        return [go_released, go_isolate_day]

Here we have added ``compliance_fraction``, set to ``0.5`` to represent,
on average, to 50% of individuals complying with the need to go into
self-isolation. This fraction is passed as ``fraction`` to the
``go_isolate`` function.

Run ``metawards`` using;

.. code-block:: bash

   metawards -d lurgy4 -D demographics.json -a ExtraSeedsLondon.dat --mixer mix_isolate --mover move_isolate --nsteps 365

You should see that infection spreads quickly, as the 50% of individuals who
each day decide not to self-isolate from stage 2 of the lurgy infect the
susceptible population. The outbreak should end with approximately
two thirds of those infected recovering in ``released`` (and so had
isolated), and one third recovering in ``home`` (having never isolated).
A plot of the the demographics shows this more clearly;

.. code-block:: bash

   metawards-plot -i output/trajectory.csv.bz2

.. image:: ../../images/tutorial_6_4_1.jpg
   :alt: Demographics for 50% self-isolation compliance

Scanning compliance
-------------------

We can investigate the effect of different levels of compliance by
scanning through ``compliance_fraction``. Update your ``move_isolate.py``
to read;

.. code-block:: python

    from metawards import Population
    from metawards import Networks
    from metawards.movers import go_isolate, go_to

    def move_isolate(network: Networks, population: Population, **kwargs):
        user_params = network.params.user_params

        ndays = 7
        isolate_stage = 3
        compliance_fraction = user_params["compliance"]

        day = population.day % 7
        isolate = f"isolate_{day}"

        go_isolate_day = lambda **kwargs: go_isolate(
                                            go_from="home",
                                            go_to=isolate,
                                            self_isolate_stage=isolate_stage,
                                            fraction=compliance_fraction,
                                            **kwargs)

        go_released = lambda **kwargs: go_to(go_from=isolate,
                                            go_to="released",
                                            **kwargs)

        return [go_released, go_isolate_day]

The only change here is that we now set ``compliance_fraction`` from the
``compliance`` user defined parameter.

Now create a scan file called ``scan_compliance.dat`` that contains;

::

    .compliance
        1.00
        0.95
        0.90
        0.85
        0.80
        0.75
        0.70
        0.65
        0.60
        0.55
        0.50

This scans ``compliance`` from ``1.00`` to ``0.50`` in increments of ``0.05``.

This is the command to run a single scan locally;

.. code-block:: bash

   metawards -d lurgy4 -D demographics.json -a ExtraSeedsLondon.dat --mixer mix_isolate --mover move_isolate --extractor extract_none --nsteps 365 -i scan_compliance.dat

You will likely need a cluster to perform repeats, so here is a suitable
PBS and slurm job file;

::

    #!/bin/bash
    #PBS -l walltime=12:00:00
    #PBS -l select=4:ncpus=64:mem=64GB
    # The above sets 4 nodes with 64 cores each

    source $HOME/envs/metawards/bin/activate

    # change into the directory from which this job was submitted
    cd $PBS_O_WORKDIR

    metawards -d lurgy4 -D demographics.json -a ExtraSeedsLondon.dat \
            --mixer mix_isolate --mover move_isolate --extractor extract_none \
            --nsteps 365 -i scan_compliance.dat --repeats 8 \
            --nthreads 16 --force-overwrite-output --no-spinner --theme simple

::

    #!/bin/bash
    #SBATCH --time=01:00:00
    #SBATCH --ntasks=4
    #SBATCH --cpus-per-task=64
    # The above sets 4 nodes with 64 cores each

    source $HOME/envs/metawards/bin/activate

    metawards -d lurgy4 -D demographics.json -a ExtraSeedsLondon.dat \
            --mixer mix_isolate --mover move_isolate --extractor extract_none \
            --nsteps 365 -i scan_compliance.dat --repeats 8 \
            --nthreads 16 --force-overwrite-output --no-spinner --theme simple

You can generate an overview animation using;

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2
   metawards-plot --animate output/overview*

The resulting animation should look like this;

.. image:: ../../images/tutorial_6_4_2.gif
   :alt: Animation showing effect on the outbreak of different levels
         of compliance with self-isolation

From this it is clear the high compliance with self-isolation is needed
to prevent a large-scale outbreak. You can see exactly the effect
by plotting the average number of recovereds as a function of the level
of compliance at the end of the year (365 days). This can be done
in Python pandas, R or Excel. For example, here is how you would make
this plot in pandas;

.. code-block:: python

   >>> import pandas as pd
   >>> import matplotlib as mpl
   >>> df = pd.read_csv("output/results.csv.bz2")
   >>> day_365 = df[df["day"] == 365]
   >>> day_365_mean = day_365.groupby("fingerprint").mean()
   >>> ax = day_365_mean.plot.line(x=".compliance", y="R")
   >>> ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
   >>> fig = ax.get_figure()
   >>> fig.tight_layout()
   >>> fig.savefig("compliance.jpg")

The figure I get is shown here;

.. image:: ../../images/tutorial_6_4_3.jpg
   :alt: Effect of compliance on the number of infections (recovereds)

For the lurgy, compliance of 90% still leads to ~6 million infections over
the year, while compliance of 80% leads to ~20 million. If compliance drops
to 50%, then ~39 million individuals are infected.
========================
Custom advance functions
========================

So far we have just been creating our own custom iterator functions.
You can also create your own custom advance functions.

Create a new python file called ``advance.py`` and copy in the
below text;

.. code-block:: python

  from metawards.utils import Console

  def advance_function(**kwargs):
      Console.debug("Hello advance_function")

  def iterate_advance(**kwargs):
      Console.debug("Hello iterate_advance")
      return [advance_function]

This defines two functions. One is the advance function we will use.
This takes ``**kwargs`` as arguments, but returns nothing.

The ``iterate_advance`` function is our iterator, which returns
a list of just our ``advance_function`` to run.

Use this iterator via the ``metawards`` command;

.. code-block:: bash

  metawards -d lurgy3 --additional ExtraSeedsLondon.dat --iterator advance --debug

You should see that your "Hellos" from these functions are printed.

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    [15:35:18]                       Hello iterate_advance                       advance.py:9
                                    Hello advance_function                       advance.py:5
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

Creating advance_lockdown
-------------------------

We can write a custom advance function that represents a lockdown. To do
this we will assume that there are no new work infections during a lockdown.
We want to create an iterator that advances normally before the lockdown
is triggered, and then advances using only the play function afterwards.

To do this, create a new python file called ``lockdown.py``, and copy in;

.. code-block:: python

  from metawards.iterators import iterate_working_week, \
                                  advance_infprob, \
                                  advance_fixed, \
                                  advance_play
  from metawards.utils import Console

  def advance_lockdown(**kwargs):
      Console.debug("We are on lockdown")
      advance_infprob(**kwargs)
      advance_play(**kwargs)

  def iterate_lockdown(population, **kwargs):
      if population.day > 20:
          return [advance_lockdown]
      else:
          return iterate_working_week(population=population,
                                      **kwargs)

Here, the ``iterate_lockdown`` function behaves identically
to ``iterate_working_week`` for the first 20 days of the
outbreak. However, after day 20, ``advance_lockdown`` is called.
This function calls :meth:`~metawards.iterators.advance_infprob`
to advance the infection probabilities before calling
:meth:`~metawards.iterators.advance_play` to advance the
"play" infections only. This is the equivalent of making
every day into a weekend, e.g. there is no fixed or predictable
travelling from home to work (or school).

Run this iterator and draw a graph using;

.. code-block:: bash

  metawards -d lurgy3 --additional ExtraSeedsLondon.dat --iterator lockdown
  metawards-plot -i output/results.csv.bz2 --format jpg --dpi 150

You should see that the overview plot is very similar to your
"weekend" only graphs that you created :doc:`on the last page <02_weekend>`.
This is unsurprising, as this "lockdown" has not reduced the random
infections.

Scaling the infection rate
--------------------------

To model a reduction in the rate of new infections, the
:meth:`~metawards.iterators.advance_infprob` function can take
and extra ``scale_rate`` argument that is used to scale all of the
infection rates that are calculated. This is a blunt tool, but it can
be used to model the reduced infection rates that a lockdown aim to
achieve.

To set ``scale_rate``, edit your ``lockdown.py`` file to contain;

.. code-block:: python

  from metawards.iterators import iterate_working_week, \
                                  advance_infprob, \
                                  advance_fixed, \
                                  advance_play

  from metawards.utils import Console

  def advance_lockdown(**kwargs):
      Console.debug("We are on lockdown")
      advance_infprob(scale_rate=0.25, **kwargs)
      advance_play(**kwargs)

  def iterate_lockdown(population, **kwargs):
      if population.day > 20:
          return [advance_lockdown]
      else:
          return iterate_working_week(population=population,
                                      **kwargs)

All we have done is set ``scale_rate=0.25`` in
:meth:`~metawards.iterators.advance_infprob`. This represents
a four-fold reduction in the infection rate.

Run the model and generate graphs again using;

.. code-block:: bash

  metawards -d lurgy3 --additional ExtraSeedsLondon.dat --iterator lockdown
  metawards-plot -i output/results.csv.bz2 --format jpg --dpi 150

You should now see a dramatic reduction in the infection, e.g.
my overview graph looks like this;

.. image:: ../../images/tutorial_3_3_overview.jpg
   :alt: Overview image of a quick lockdown
===============
Local lockdowns
===============

In the :doc:`last chapter <07_cutoff>` it was implied that the parameters
for equation 1, and the action of the ``cutoff`` parameter, were global,
and applied to all wards equally. If this was true, then it would make
it difficult for ``metawards`` to model ward-specific behaviour, such
as local lockdowns, or differences in local behaviour.
In reality, ``metawards`` provides support for ward-specific values of
the scaling and cutoff parameters.

Limiting movement at a local level
----------------------------------

Every ward can have its own value
of the cutoff parameter. Movement between two wards is only permitted
if the distance between wards is less than the minimum of the
cutoff distance of each ward, and the value of
:meth:`Parameters.dyn_dist_cutoff <metawards.Parameters.dyn_dist_cutoff>`,
e.g. if this condition is true;

3. :math:`D_\text{ij} < \text{min}( C_i, C_j, C_\text{global} )`

where;

* :math:`D_\text{ij}` is the distance between the centres of wards
  :math:`i` and :math:`j`,
* :math:`C_i` is the local cutoff parameter for ward :math:`i`,
* :math:`C_j` is the local cutoff parameter for ward :math:`j`, and
* :math:`C_\text{global}` is the global cutoff distance set via
  :meth:`Parameters.dyn_dist_cutoff <metawards.Parameters.dyn_dist_cutoff>`

Scaling FOI at a local level
----------------------------

Every ward can have its own value of the FOI scaling parameter. In reality,
equation 1 is actually;

4. :math:`F(w) = S \times S_l(w) \times U(t) \times \sum_s [ C_s \beta_s N_s(w) ]`

where;

* :math:`S_l(w)` is the local scaling factor for ward :math:`w`. This acts
  together with :math:`S`, :math:`U(t)` and :math:`C_s` to give you a lot
  of control over how infectious individuals at each disease stage in each
  ward contribute to the ward's FOI.

Reading and writing local parameters
------------------------------------

You can read and write local ward parameters from within a custom
iterator. The local ward parameters are stored in the
:class:`network.nodes <metawards.Nodes>` object. This provides arrays
which are indexed by ward ID (i.e. counting up from ``1`` to ``nnodes + 1``).

:meth:`network.nodes.scale_uv <metawards.Nodes.scale_uv>` contains the
:math:`S_l(w)` values for each node, while
:meth:`network.nodes.cutoff <metawards.Nodes.cutoff>` contains the
cutoff values for each node (in kilometers).

The scaling factors for each node default to ``1.0`` (meaning no local scaling),
while the cutoff defaults to a large distance that is greater than the
distance between any two points on Earth (meaning thus no local cutoff).

As these values are used as part of the FOI calculation, they must be
read and set during the ``foi`` stage of the model day. For example,
create a new iterator called ``lockdown.py`` and copy in the below;

.. code-block:: python

    from metawards.iterators import iterate_default

    def advance_lockdown(network, **kwargs):
        scale_uv = network.nodes.scale_uv
        cutoff = network.nodes.cutoff

        # set a lockdown in even-numbered wards
        for i in range(2, network.nnodes, 2):
            scale_uv[i] = 0.0
            cutoff[i] = 0.0


    def iterate_lockdown(stage, **kwargs):
        # get the default functions for this stage
        funcs = iterate_default(stage, **kwargs)

        if stage == "foi":
            # make sure that advance_lockdown is called
            # first in the foi stage
            return [advance_lockdown] + funcs
        else:
            return funcs

This iterator defines the ``advance_lockdown`` advance function. This
simply loops over all even-numbered wards and sets
:meth:`network.nodes.scale_uv <metawards.Nodes.scale_uv>` and
:meth:`network.nodes.cutoff <metawards.Nodes.cutoff>` to zero.
In effect, this places half of the country into extreme lockdown,
where the disease is unable to spread.

The ``iterate_lockdown`` function takes the ``stage`` parameter. This
tells ``metawards`` that this iterator wants to specify the advance
functions to call at different stages. By default, this returns the
default advance functions for that stage (as returned by
:func:`~metawards.iterators.iterate_default`). For the ``foi`` stage,
this return ``advance_lockdown`` before the default functions, thereby
ensuring that ``advance_lockdown`` changes the
:meth:`network.nodes.scale_uv <metawards.Nodes.scale_uv>` and
:meth:`network.nodes.cutoff <metawards.Nodes.cutoff>` ward-local parameters
before they are used to calculate the force of infection of each ward.

You can run this iterator using;

.. code-block:: bash

   metawards -d lurgy3 -a ExtraSeedsLondon.dat --iterator lockdown

You should see that the infection spreads to only half of the country,
as the lurgy can only infect the half of the population that are
resident of visiting wards that are not in complete lockdown. For example,
I see;

::

    ...

    ─────────────────────────────────────────────── Day 67 ───────────────────────────────────────────────
    S: 45526353  E: 871045  I: 5537877  R: 4146802  IW: 4293  POPULATION: 56082077
    Number of infections: 6408922

    ─────────────────────────────────────────────── Day 68 ───────────────────────────────────────────────
    S: 44593470  E: 901872  I: 5911611  R: 4675124  IW: 4294  POPULATION: 56082077
    Number of infections: 6813483

    ─────────────────────────────────────────────── Day 69 ───────────────────────────────────────────────
    S: 43632470  E: 932883  I: 6269338  R: 5247386  IW: 4294  POPULATION: 56082077
    Number of infections: 7202221

    ...

    ────────────────────────────────────────────── Day 184 ───────────────────────────────────────────────
    S: 28691853  E: 0  I: 2  R: 27390222  IW: 0  POPULATION: 56082077
    Number of infections: 2

    ────────────────────────────────────────────── Day 185 ───────────────────────────────────────────────
    S: 28691853  E: 0  I: 2  R: 27390222  IW: 0  POPULATION: 56082077
    Number of infections: 2

    ────────────────────────────────────────────── Day 186 ───────────────────────────────────────────────
    S: 28691853  E: 0  I: 0  R: 27390224  IW: 0  POPULATION: 56082077
    Number of infections: 0
    Ending on day 186

Custom ward-local parameters
----------------------------

You can also read and write your own custom ward-local parameters.
You do this by calling
:meth:`network.nodes.get_custom <metawards.Nodes.get_custom>`. For example;

.. code-block:: python

   my_params = network.nodes.get_custom("my_params", default=0.0)

will return the custom ward-local parameters called ``my_params``. If these
don't exist, then they are created, with each ward given a default
starting value specified by ``default`` (here ``0.0``). The return value
is the array indexed by ward ID. This can be read and written in an identical
way to :meth:`network.nodes.scale_uv <metawards.Nodes.scale_uv>` and
:meth:`network.nodes.cutoff <metawards.Nodes.cutoff>`.

.. note::

   Custom ward-level parameters are always stored as an array of
   floating point numbers.

You can use custom parameters to store or manipulate extra ward-level data.
For example, edit your ``lockdown.py`` iterator to read;

.. code-block:: python

    from metawards.iterators import iterate_default
    from metawards.utils import Console

    def advance_lockdown(network, workspace, **kwargs):
        # get the ward-specific scaling and cutoff parameters
        scale_uv = network.nodes.scale_uv
        cutoff = network.nodes.cutoff

        # get the custom parameter 'in_lockdown' which we will
        # initialise to 0 (meaning false)
        in_lockdown = network.nodes.get_custom("in_lockdown", default=0)

        # count of number of case-free days per ward
        case_free_days = network.nodes.get_custom("case_free_days", default=0)

        # get the total number of infections from the workspace
        I_in_wards = workspace.I_in_wards

        # loop over all wards
        for i in range(1, network.nnodes + 1):
            # is this ward in lockdown?
            if in_lockdown[i]:
                # has the number of infections dropped to zero? If so,
                # then leave lockdown
                if I_in_wards[i] == 0:
                    # we need 28 case-free days before releasing lockdown
                    if case_free_days[i] > 28:
                        Console.debug(f"Ward {i} leaving lockdown")
                        # stay on high vigilence, so keep actions that
                        # reduce beta to 20% of normal
                        scale_uv[i] = 0.2
                        cutoff[i] = 99999.99
                        in_lockdown[i] = 0
                        case_free_days[i] = 0
                    else:
                        case_free_days[i] += 1
                        Console.debug(f"Ward {i} case_free_days equals {case_free_days[i]}")
                else:
                    case_free_days[i] = 0

            # if not, then enter lockdown if the number of infections
            # goes above 5
            elif I_in_wards[i] > 5:
                Console.debug(f"Ward {i} entering lockdown")
                in_lockdown[i] = 1
                case_free_days[i] = 0

                # stop all travel and enact measures that
                # will scale down beta to 1% of normal
                cutoff[i] = 0.0
                scale_uv[i] = 0.01

        # get the number of wards in lockdown
        num_lockdown = int(sum(in_lockdown))

        if num_lockdown > 0:
            Console.print(f"Number of wards in lockdown equals {num_lockdown}")


    def iterate_lockdown(stage, **kwargs):
        # get the default functions for this stage
        funcs = iterate_default(stage=stage, **kwargs)

        if stage == "foi":
            return [advance_lockdown] + funcs
        else:
            return funcs

In this case ``advance_lockdown`` will move individual wards in and out
of local lockdowns depending on the number of infections in that ward
(read from the :meth:`workspace.I_in_wards <metawards.Workspace.I_in_wards>`
array from the passed :class:`~metawards.Workspace` object).

Two ward-local custom parameters are used to record whether or not
a ward is in a local lockdown;

* ``in_lockdown`` is 1 if the ward is in lockdown, and 0 if it is not
* ``case_free_days`` is the count of the number of consecutive days in
  a ward without an infection (really detectable infection, e.g. an
  individual in the ``I`` state).

The ``advance_lockdown`` function works by looping over all wards and
seeing if the ward is in lockdown by checking the ``in_lockdown`` value
for that ward. If it is, and if the number of infections is zero, then
it checks if more than 28 case-free days have passed. If they have, then
the local lockdown is relaxed, travel is allowed (cutoff is set to
a large value) and the scaling factor is increased to ``0.2`` (implying
that measures such as mask wearing, physical distancing etc. are still
followed).

If 28 days have not passed, then the number of case-free days is incremented.

If the ward is not in local lockdown, then if the number of local detected
infections goes above 5 then a local lockdown is initiated. Travel is
halted (cutoff is set to 0) and stringent measures are taken such that
the scaling factor is ``0.01`` (implying that ``beta`` is scaled down by 99%).

Finally, the number of wards in lockdown is calculated as the sum of
the ``in_lockdown`` custom parameter, and is printed to the screen.

.. note::

   Note that there are some additional :meth:`Console.debug <metawards.utils.Console.debug>`
   statements in the function that print out debug lines when wards move
   in and out of lockdown.

You can run this iterator using;

.. code-block:: bash

   metawards -d lurgy3 -a ExtraSeedsLondon.dat --iterator lockdown

You should see that the number of wards in lockdown increases as the
disease spreads. The spread is slowed down, but as wards come out of
lockdown they are sometimes re-infected, and have to re-lockdown. You
may see wave like behaviour as the disease is slowly brought under
control. For example, for me, the plot of the outbreak, produced via;

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2

shows the following ``output/overview.jpg`` plot;

.. image:: ../../images/tutorial_3_8.jpg
   :alt: Outbreak controlled using local lockdowns

This wave behaviour is more clear if we make the entering of exiting of
local lockdown more extreme. For example, update your ``lockdown.py``
to read;

.. code-block:: python

    from metawards.iterators import iterate_default
    from metawards.utils import Console

    def advance_lockdown(network, workspace, **kwargs):
        # get the ward-specific scaling and cutoff parameters
        scale_uv = network.nodes.scale_uv
        cutoff = network.nodes.cutoff

        # get the custom parameter 'in_lockdown' which we will
        # initialise to 0 (meaning false)
        in_lockdown = network.nodes.get_custom("in_lockdown", default=0)

        # count of number of case-free days per ward
        case_free_days = network.nodes.get_custom("case_free_days", default=0)

        # get the total number of infections from the workspace
        I_in_wards = workspace.I_in_wards

        # loop over all wards
        for i in range(1, network.nnodes + 1):
            # is this ward in lockdown?
            if in_lockdown[i]:
                # has the number of infections dropped to zero? If so,
                # then leave lockdown
                if I_in_wards[i] == 0:
                    # we need 28 case-free days before releasing lockdown
                    if case_free_days[i] > 28:
                        Console.debug(f"Ward {i} leaving lockdown")
                        # completely relax the lockdown
                        scale_uv[i] = 1.0
                        cutoff[i] = 99999.99
                        in_lockdown[i] = 0
                        case_free_days[i] = 0
                    else:
                        case_free_days[i] += 1
                        Console.debug(f"Ward {i} case_free_days equals {case_free_days[i]}")
                else:
                    case_free_days[i] = 0

            # if not, then enter lockdown if the number of infections
            # goes above 5
            elif I_in_wards[i] > 5:
                Console.debug(f"Ward {i} entering lockdown")
                in_lockdown[i] = 1
                case_free_days[i] = 0

                # stop all travel and enact measures that
                # stop all local transmission (beta is 0)
                cutoff[i] = 0.0
                scale_uv[i] = 0.0

        # get the number of wards in lockdown
        num_lockdown = int(sum(in_lockdown))

        if num_lockdown > 0:
            Console.print(f"Number of wards in lockdown equals {num_lockdown}")


    def iterate_lockdown(stage, **kwargs):
        # get the default functions for this stage
        funcs = iterate_default(stage=stage, **kwargs)

        if stage == "foi":
            return [advance_lockdown] + funcs
        else:
            return funcs

The only change is that ``scale_uv[i]`` is set to ``0.0`` for wards that
are in lockdown (i.e. there is no more spread), while ``scale_uv[i]``
is returned to ``1.0`` for wards that leave lockdown. This extreme switching
when entering and leaving lockdown causes waves of infection that
spread across wards, e.g. when I run this model I see;

.. image:: ../../images/tutorial_3_8_2.jpg
   :alt: Outbreak with local lockdowns resulting in waves of disease
=====================
Where is the weekend?
=====================

It may not have escaped your attention that every day is a work day
in this model. While this may seem unrealistic, we must remember that
these are random, imperfect models, based on very noisy data.
Adding more "realism" may be counter-productive, especially as
modern working patterns mean that there is blurring of the line between
work days and weekends.

This doesn't mean that we can't model a weekend. Indeed, ``metawards``
is really flexible and you can customise exactly what is performed
for each model day.

Creating the weekend
--------------------

Create a new directory called ``weekend`` and copy into it your
``lurgy3.json`` disease parameters. Change into this directory and
create a new file called ``weekend.py``, and copy into it the below
code.

.. code-block:: python

    from metawards.utils import Console

    def iterate_weekend(**kwargs):
        Console.print("Hello iterate_weekend")

        return []

This is a simple function called ``iterate_weekend``. It takes an
unspecified number of
`keyword arguments (**kwargs) <https://book.pythontips.com/en/latest/args_and_kwargs.html>`__
(more about these later). It returns an empty list (``[]``). All it does
is print ``Hello iterate_weekend`` to the screen.

.. note::
   Notice that you *must* print using the
   :meth:`Console.print <metawards.utils.Console.print>` function of
   :meth:`~metawards.utils.Console`. This ensures that all printing
   goes to the right place and stays sane when multiple processes
   and threads all try to print at the same time. It also ensures
   that everything that is printed to the screen also gets printed
   to a file for safekeeping (``output/console.log.bz2``). It is
   very important that information is not lost when running a job.

You can run this function by starting ``ipython`` in this directory
and typing;

.. code-block:: ipython

   In [1]: import weekend
   In [2]: weekend.iterate_weekend()
   Hello iterate_weekend
   Out[2]: []

You can tell ``metawards`` to call this function every iteration
using the ``--iterator`` command-line argument. Type;

.. code-block:: bash

   metawards metawards -d lurgy3 --additional ExtraSeedsLondon.dat --iterator weekend

You should see a very different outbreak to what you have before, e.g.

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello iterate_weekend
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello iterate_weekend
    seeding play_infections[0][255] += 5
    S: 56082072  E: 5  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 5

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello iterate_weekend
    S: 56082072  E: 4  I: 1  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 5

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello iterate_weekend
    S: 56082072  E: 2  I: 3  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 5

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello iterate_weekend
    S: 56082072  E: 2  I: 2  R: 1  IW: 0  POPULATION: 56082077
    Number of infections: 4

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 5 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello iterate_weekend
    S: 56082072  E: 1  I: 2  R: 2  IW: 0  POPULATION: 56082077
    Number of infections: 3

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 6 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello iterate_weekend
    S: 56082072  E: 0  I: 2  R: 3  IW: 0  POPULATION: 56082077
    Number of infections: 2

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 7 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello iterate_weekend
    S: 56082072  E: 0  I: 1  R: 4  IW: 0  POPULATION: 56082077
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 8 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello iterate_weekend
    S: 56082072  E: 0  I: 0  R: 5  IW: 0  POPULATION: 56082077
    Number of infections: 0
    Infection died ... Ending on day 9

What happened here? Well, just as you imported ``weekend`` into ``ipython``
and called the ``iterate_weekend`` function, so too has ``metawards``.
The ``--integrator`` option tells ``metawards`` to import the ``weekend``
module. ``metawards`` then automatically found the first function in that
module whose name started with ``iterate``, in this case ``iterate_weekend``.

Then, ``metawards`` called this function for every iteration of the
**model run**.

You can name your function whatever you want, e.g. edit ``weekend.py``
to read;

.. code-block:: python

  from metawards.utils import Console

  def another_function(**kwargs):
      Console.print("Hello another_function")

      return []


  def iterate_weekend(**kwargs):
      Console.print("Hello iterate_weekend")

      return []

This has added another function called ``another_function``. You can tell
``metawards`` to use this function using
``--iterator weekend::another_function``. Try running this using the
command below;

.. code-block:: bash

  metawards -d lurgy3 --additional ExtraSeedsLondon.dat --iterator weekend::another_function

You should see ``Hello another_function`` is now printed for
every iteration.

.. warning::
   Sometimes you may see ``metawards`` exit with a warning that it can't
   find your iterator function. This is likely because there is a typo
   or syntax error in your iterator. ``metawards`` does its best to
   detect these and report them to you, so check above the error in the
   output to see if there is anything helpful. If not, then run your
   iterator in python to see if you get any errors, e.g. if your iterator
   is in a file called ``iterator.py`` then type ``python iterator.py``.
   If there is an error, then that will be printed to the screen.

Printing debug output
---------------------

In general, you should only print things to the screen if they will be useful
for the user of the program. Sometimes when developing you want to print
some debugging output that can verify that everything is working. To do this,
using :meth:`Console.debug <metawards.utils.Console.debug>`. For example,
change your iterator to;

.. code-block:: python

    from metawards.utils import Console

    def iterate_weekend(**kwargs):
        Console.debug("Hello iterate_weekend")

        return []

Now, you will only see this print output if the ``--debug`` option is passed
to ``metawards``, e.g.

.. code-block:: bash

    metawards -d lurgy2 --iterator weekend --debug

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    [15:23:08]                       Hello iterate_weekend                       weekend.py:5
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

Note that the time of the debug string, and the line and file of the debug
statement are included. You can also easily print the values of variables
using the ``variables`` keyword argument to
:meth:`~metawards.utils.Console.debug`, e.g.

.. code-block:: python

    from metawards.utils import Console

    def iterate_weekend(**kwargs):
        a = 42
        b = "This is a string"

        Console.debug("Hello iterate_weekend", variables=[a, b])

        return []

.. code-block:: bash

    metawards -d lurgy2 --iterator weekend --debug

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    [15:25:38]                       Hello iterate_weekend                       weekend.py:8

     Name │ Value
    ══════╪══════════════════
        a │ 42
        b │ This is a string

More information about debug strings, debugging levels, and how you can leave
these debug strings in production code is
:doc:`available here <../../devsupport>`.

Advancing the outbreak
----------------------

You may have noticed that the disease outbreak was not advancing during
any of the runs using your custom weekend iterator. The output showed
that five initial infections were seeded. These progressed through
the disease stages until all five individuals moved into the **R**
state.

The reason the disease hasn't advanced is because you haven't supplied
any functions that are used to advance the outbreak. The job of
the iterator function is to return the functions that are needed to
advance an outbreak (so-called ``advance functions``).

You can write an advance function by editing ``weekend.py`` to contain;

.. code-block:: python

  from metawards.iterators import advance_infprob, advance_play
  from metawards.utils import Console

  def iterate_weekend(**kwargs):
      Console.debug("Hello iterate_weekend")

      return [advance_infprob, advance_play]

In this code you have imported the :meth:`~metawards.iterators.advance_infprob`
and :meth:`~metawards.iterators.advance_play` advance functions.
These were described on the :doc:`last page <01_iterators>`. By returning
them from ``iterate_weekend`` you have told ``metawards`` to call them,
one after another, to advance the outbreak. If you now run
``metawards`` using this new ``weekend.py`` via;

.. code-block:: bash

   metawards -d lurgy3 --additional ExtraSeedsLondon.dat --iterator weekend

you will see that the outbreak now advances throughout the population.
However, each day now only progresses new infections using the "play" mode
:meth:`~metawards.iterators.advance_play`. The "work" mode
:meth:`~metawards.iterators.advance_fixed`, is not used, meaning
that every day is now modelled as like a weekend.

Create an overview graph of your "weekend only" run and compare it to
the results from the "weekday only" runs in
:doc:`part 2 <../part02/05_refining>`. Do you see a difference?

My graph is shown below;

.. image:: ../../images/tutorial_3_2_1_overview.jpg
   :alt: Overview image of a weekend only run

It is clear that the outbreak is now much smaller, peaking at 7 million
as opposed to nearly 20 million. The peak is also broadened
out, with the outbreak lasting months rather than weeks.

Changing iterators with time
----------------------------

A week of only weekends is also not realistic. We can however create
a function that can choose which advance functions to return based
on the day of the outbreak.

To do this, create a new python file called ``week.py`` and copy into
it the code below;

.. code-block:: python

  from metawards.iterators import advance_infprob, \
                                  advance_fixed, \
                                  advance_play
  from metawards.utils import Console


  def iterate_week(population, **kwargs):
      date = population.date

      Console.debug(f"Creating functions for {date}")

      if date.weekday() < 5:
          Console.debug("This is a weekday")
          return [advance_infprob,
                  advance_fixed,
                  advance_play]
      else:
          Console.debug("This is a weekend")
          return [advance_infprob,
                  advance_play]

This has created an ``iterate_week`` function. This has a slightly
different signature to ``iterate_weekend``, in that it accepts
the ``population`` argument. Every iterator is passed a lot of
arguments, most of which are ignored by the ``**kwarg`` variables.

When you need an argument you name it in the function. In this case,
we need the ``population`` argument. This is a
:class:`~metawards.Population` object, which contains the distribution
of the population across the different **S**, **E**, **I** states,
plus the current date of the outbreak (
:meth:`Population.date <~metawards.Population.date>`).

The ``date`` is a standard `Python date object <https://docs.python.org/3/library/datetime.html>`__.
The ``.weekday()`` function returns a number from 0-6 to correspond
with Monday to Sunday (0 is Monday, 6 is Sunday).

If the weekday is less than 5, then the day must be a weekday. Hence
the ``iterate_week`` function returns the infprob, fixed and play
advance functions. Otherwise, the day must be a weekend, and so
only the infprob and play advance functions are returned.

Run ``metawards`` using this new iterator and see what happens;

.. code-block:: bash

  metawards -d lurgy3 --additional ExtraSeedsLondon.dat --iterator week
  metawards-plot -i output/results.csv.bz2 --format jpg --dpi 150

You should see something similar to this;

.. image:: ../../images/tutorial_3_2_2_overview.jpg
   :alt: Overview image of a weekend only run

There is a significant spread in the infection during weekdays,
but then this growth falls back at weekends.

.. note::
  This "week" iterator is so important that it is supplied
  as the :meth:`metawards.iterators.iterate_working_week`
  iterator. You can use this via the command line option
  ``--iterator iterate_working_week``. Similarly there
  is :meth:`metawards.iterators.iterate_weekday` function
  to iterate as a weekday only, and
  :meth:`metawards.iterators.iterate_weekend` to iterate
  as weekends only.

.. note::
  By default the outbreak is modelled to start from today.
  You can control the start date using the ``--start-date``
  command line option.

Changing iterators with stage
-----------------------------

This method of modelling the weekend, while simple, is not correct.
We have modelled the weekend as a time when only the players move
and can be infected. The workers are ignored, meaning that this model
assumes that the workers spend their weekends at home, not interacting
with anyone else, and thus have zero risk of contracting an infection.

To model the weekend properly, we have to do something to advance the
workers. One option is to advance workers at the weekend in the same way
that we advance players. Because the data structure used to model
workers is different, we can't just call the
:meth:`~metawards.iterators.advance_play` function on the workers.
Instead, metawards comes with
:meth:`~metawards.iterators.advance_work_to_play`. This advances the
workers as if they were players.

Using this, the `iterate_week` iterator could look like;

.. code-block:: python

  from metawards.iterators import advance_infprob, \
                                  advance_work_to_play, \
                                  advance_play
  from metawards.utils import Console


  def iterate_week(population, **kwargs):
      date = population.date

      Console.debug(f"Creating functions for {date}")

      if date.weekday() < 5:
          Console.debug("This is a weekday")
          return [advance_infprob,
                  advance_fixed,
                  advance_play]
      else:
          Console.debug("This is a weekend")
          return [advance_infprob,
                  advance_work_to_play,
                  advance_play]

However, this would also not be correct. The issue is that, as well
as advancing the infection, the iterator is also responsible for
calculating the force of infection (FOI) resulting from infected
individuals. This is calculated during the "foi" stage of the day.

Each model day is divided into a series of stages;

1. `initialise` : used to initialise any variables that day
2. `setup` : used to set up any additional variables or infections
3. `foi` : used to calculate the FOI resulting from existing infections
4. `infect` : used to calculate and enact new infections
5. `analyse` : used to calculate / analyse the data from infections

In addition, for a model run, there is a `finalise` stage,
which is called once at the end of a model run, at which
all outputs from a model run are written to disk or a database.
Then, for a collection of model runs, there is a
`summary` stage which is called once at the end of all model runs,
to calculate and write out summary statistics from all runs.

By default, a custom iterator will only specify the advance functions
that are used for the `infect` stage. The advance functions used for
other stages are supplied by
:meth:`~metawards.iterators.iterate_default`. By default, the
:meth:`~metawards.iterators.advance_foi` function is used to
calculate the FOI. This assumes that all workers behave like
workers, and all players behave like players. It does not work
when we use :meth:`~metawards.iterators.advance_foi_work_to_play`,
because this makes the workers behave like players.

We thus need to use the :meth:`~metawards.iterators.advance_foi_work_to_play`
function instead. To do this, we need to tell metawards to use that
function in the `foi` stage. We can control which stage our iterator
operates by passing in the `stage` argument. For example;

.. code-block:: python

  from metawards.iterators import advance_infprob, \
                                  advance_work_to_play, \
                                  advance_play, \
                                  advance_foi, \
                                  advance_foi_work_to_play, \
                                  advance_recovery
  from metawards.utils import Console


  def iterate_week(stage, population, **kwargs):
      date = population.date

      Console.debug(f"Creating functions for {date}")

      is_weekend = date.weekday() < 5


      if is_weekend:
          Console.debug("This is a weekend")
          if stage == "foi":
              return [advance_foi_work_to_play,
                      advance_recovery]
          elif stage == "infect":
              return [advance_infprob,
                      advance_work_to_play,
                      advance_play]
          else:
              return iterator_default(stage=stage, **kwargs)

       else:
          Console.debug("This is a weekday")
          if stage == "foi":
              return [advance_foi,
                      advance_recovery]
          elif stage == "infect:
              return [advance_infprob,
                      advance_fixed,
                      advance_play]
          else:
              return iterator_default(stage=stage, **kwargs)

.. note::

   :meth:`~metawards.iterators.advance_recovery` is another advance
   function that must be returned at the `foi` stage. This advance
   function is used to advance infected individuals along
   the stages of the disease (e.g. E to I to R).

.. note::

   :meth:`~metawards.iterators.iterate_default` returns the default
   set of advance functions for each stage. It makes sense
   to call this for the stages that you are not explicitly
   handling in your iterator, e.g.
   `return iterate_default(stage=stage, **kwargs)` will return
   the default advance functions for the specified stage.

With those changes, the `iterate_week` iterator would work as a model
where all workers become players at the weekends.

To make things easier, metawards provides a built-in
:meth:`~metawards.iterators.iterate_weekend` iterator for
iterating a weekend day, plus a
:meth:`~metawards.iterators.iterate_working_week` iterator that
will use :meth:`~metawards.iterators.iterate_default` for weekdays,
and :meth:`~metawards.iterators.iterate_weekend` for weekends.
Use these if you want to model a working week where workers
behave like players at the weekend.
=====================
How a day is modelled
=====================

You have now used multiple *model runs* to refine the disease parameters
for the lurgy. However, these parameters have been refined for the
model that ``metawards`` uses to represent interactions of individuals.

Work and play
-------------

``metawards`` models disease transmission via two main modes;

1. Infections at "work"

  These are random infections that individuals acquire in the electoral
  ward in which they "work", based on random interactions with others
  who are also in that "work" ward. In this nomenclature, "work" means
  any regular (daily), predictable travel that an individual makes.

1. Infections at "play"

  These are random infections that individuals acquire in the
  electoral ward in which they live, based on random interactions
  with others who are also in their ward. In this nomenclature, "play"
  means random, unpredictable interactions that are not regular.

.. warning::
  Do not read too much into the definitions of "work" and "play". They
  have very broad meanings in ``metawards``, and, in essence, capture
  the difference between predicatable daily travel and interactions
  (commuting, colleagues in an office, students in a school), and
  the random interations (partying, shopping, playing in the park).
  Most of the hard data science behind ``metawards`` is constructing
  the information in
  `MetaWardsData <https://github.com/metawards/MetaWardsData>`__
  to gain this very broad overview of how individuals move and interact.

What is a normal day?
---------------------

``metawards`` uses an ``iterator`` to iterate the model outbreak forward
day by day. All of the iterators are in the
:doc:`metawards.iterators <../../api/index_MetaWards_iterators>` module.
The default iterator is :class:`~metawards.iterators.iterate_default`.

An iterator applies a sequence of functions that advance the disease step
by step through each day. These ``advance_functions`` control exactly
what happens in each electoral ward on each day.

The :func:`~metawards.iterators.iterate_default` iterator applies the
following ``advance_functions`` in sequence;

1. :func:`~metawards.iterators.advance_additional` is applied to
   add any additional seeds of the disease in the ward,
   which leads to new infections. These additional seeds represent, e.g.
   new sources of infection arriving in the ward via outside travel.

2. :func:`~metawards.iterators.advance_foi` is applied to advance the
   calculation of the *force of infection (foi)* for each ward. This must
   be called at the beginning of the day after
   :func:`~metawards.iterators.advance_additional`, as the *foi* parameters
   are used to guide the path of the outbreak in each ward for the
   rest of the day.

3. :func:`~metawards.iterators.advance_recovery` is applied to all
   individuals in a ward who are infected. This will see whether an
   individual progresses from one stage of the disease to the next.
   This decision is based on the **progress** disease parameter for the stage
   that the individual is at.

4. :func:`~metawards.iterators.advance_infprob` is applied to recalculate
   the infection probabilities needed to guide new infections. These are
   based on the *foi* parameters for each ward and the number of
   individuals who are at each stage of the disease (based on the
   **contrib_foi** disease parameter).

5. :func:`~metawards.iterators.advance_fixed` is applied to advance
   all infections that are based on fixed (predictable) movements
   of individuals (the so-called "work" mode). Infected individuals
   continue to "work" unless they become too symptomatic
   (**too_ill_to_move**, based on the parameter for the stage of the
   disease at which the individual is at).

6. Finally :func:`~metawards.iterators.advance_play` is applied to
   advance all infections that are based on random movements of
   individuals (the so-called "play" mode). Infected individuals
   continue to "play" unless they become too symptomatic
   (again controlled by the **too_ill_to_move** disease parameters).

Once all of these functions have been applied, then a day is considered
complete. The statistics for the day, e.g. numbers of individuals
who are in the **S**, **E**, **I**, **IW**, and **R** states are
collected and printed, and then the next day begins and all of
these functions are applied again.
=================
Scanning lockdown
=================

Now that we can create and scan custom variables, we can write a
proper lockdown iterator that enables us to explore different
scenarios.

Create a file called ``lockdown.inp`` and copy in the below;

::

    # Full lockdown (red)
    .scale_rate[0] = 0.05
    .can_work[0]  = False

    # Relaxed lockdown (yellow)
    .scale_rate[1] = 0.1
    .can_work[1]  = False

    # More relaxed lockdown (green)
    .scale_rate[2] = 0.1
    .can_work[2]  = True

This has defined three lockdown states, ranging from "red" (full lockdown
with strong reduction in transmission rate and working) to
"green" (relaxed lockdown with weaker reduction in transmission rate
and work allowed).

To use this data create an iterator in a file called ``lockdown.py`` and
copy in the below;

.. code-block:: python

    from metawards.iterators import advance_infprob, advance_fixed, \
                                    advance_play, iterate_working_week
    from metawards.utils import Console

    def get_lockdown_state(population):
        if not hasattr(population, "lockdown_state"):
            population.lockdown_state = -1
            population.is_locked_down = False

        if population.total > 5000:
            if population.lockdown_state == -1:
                Console.print(f"Lockdown started on {population.date}")
                population.lockdown_state = 0
                population.is_locked_down = True

            elif population.lockdown_state > 0:
                Console.print(f"Restarting lockdown on {population.date}")
                population.lockdown_state = 0
                population.is_locked_down = True

        elif population.total > 3000:
            if population.lockdown_state == 2:
                Console.print(f"Re-entering relaxed (yellow) on {population.date}")
                population.lockdown_state = 1

        elif population.total < 2000:
            if population.lockdown_state == 0:
                Console.print(f"Entering relaxed (yellow) on {population.date}")
                population.lockdown_state = 1

            elif population.total < 1000:
                if population.lockdown_state == 1:
                    Console.print(f"Entering relaxed (green) on {population.date}")
                    population.lockdown_state = 2

        return population.lockdown_state

    def advance_lockdown(network, population, **kwargs):
        params = network.params
        state = get_lockdown_state(population)
        scale_rate = params.user_params["scale_rate"][state]
        can_work = params.user_params["can_work"][state]
        Console.debug("State", variables=[scale_rate, can_work])

        advance_infprob(scale_rate=scale_rate,
                        network=network, population=population,
                        **kwargs)
        advance_play(network=network, population=population,
                    **kwargs)

        if can_work:
            advance_fixed(network=network, population=population,
                        **kwargs)

    def iterate_custom(network, population, **kwargs):
        params = network.params
        state = get_lockdown_state(population)

        if population.is_locked_down:
            Console.debug("Locked down")
            return [advance_lockdown]
        else:
            Console.debug("Normal working week day")
            return iterate_working_week(network=network,
                                        population=population,
                                        **kwargs)

The ``get_lockdown_state`` function is the most complex and different.
It uses the number of infecteds (``population.total``) to decide which
``lockdown_state`` should be used. This is an integer, with ``-1``
meaning no lockdown, ``0`` being "red", ``1`` "yellow" and ``2`` "green".

Whether or not the population is locked down is stored in the
``population.is_locked_down`` variable. If this is "False" then
``iterate_lockdown`` simply returns the result of
:meth:`~metawards.iterators.iterate_working_week`. Otherwise,
it returns the ``advance_lockdown`` function that we've defined.

This ``advance_lockdown`` function obtains the ``scale_rate`` and
``can_work`` custom user parameters from the
:class:`~metawards.Parameters` objects in the model
:class:`~metawards.Network`.

It calls :meth:`~metawards.iterators.advance_infprob` with
the set ``scale_rate`` scaling factor, before calling
:meth:`~metawards.iterators.advance_play`, and then, if
``can_work`` is "True", :meth:`~metawards.iterators.advance_fixed`.

Run ``metawards`` using the below commands and see what you get;

.. code-block:: bash

   metawards -d lurgy3 -a ExtraSeedsLondon.dat  -u lockdown.inp --iterator lockdown --debug
   metawards-plot -i output/results.csv.bz2

I see;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 36 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                                Normal working week day                     lockdown.py:63
    S: 56070689  E: 1663  I: 5889  R: 3836  IW: 1352  POPULATION: 56082077
    Number of infections: 7552

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 37 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Lockdown started on 2020-06-26
                                        Locked down                           lockdown.py:60
                                            State                              lockdown.py:43

            Name │ Value
     ════════════╪═══════
      scale_rate │ 0.05
        can_work │ False

    S: 56070608  E: 2118  I: 7192  R: 2159  IW: 80  POPULATION: 56082077
    Number of infections: 9310
    ...
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 51 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Entering relaxed (yellow) on 2020-07-10
                                        Locked down                           lockdown.py:60
                                            State                              lockdown.py:43

            Name │ Value
     ════════════╪═══════
      scale_rate │ 0.1
        can_work │ False

    S: 56069562  E: 36  I: 1518  R: 10961  IW: 55  POPULATION: 56082077
    Number of infections: 1554
    ...
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 55 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Entering relaxed (green) on 2020-07-14
                                        Locked down                           lockdown.py:60
                                            State                              lockdown.py:43

            Name │ Value
     ════════════╪═══════
      scale_rate │ 0.1
        can_work │ True

    S: 56069369  E: 46  I: 852  R: 11810  IW: 59  POPULATION: 56082077
    Number of infections: 898
    ...
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 187 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                                        Locked down                           lockdown.py:60
                                            State                              lockdown.py:43

            Name │ Value
     ════════════╪═══════
      scale_rate │ 0.1
        can_work │ True

    S: 56068649  E: 0  I: 0  R: 13428  IW: 0  POPULATION: 56082077
    Number of infections: 0
    Infection died ... Ending on day 188

with the overview graph as here;

.. image:: ../../images/tutorial_3_6_1_overview.jpg
   :alt: Overview image of a lockdown with custom parameters

Running on a cluster
--------------------

Now that this is working, we can scan through lots of different lockdown
scenarios by creating an input file that varies the ``scale_rate`` and
``can_work`` parameters. Create an input file called ``scan.csv`` and
copy in the following;

::

    # Adjust "red" state from 0.05 to 0.20
    # while adjusting "yellow" from "green" + 0.05 to 0.25
    # while adjusting "green" from "yellow" if working, or
    #                              "yellow" + 0.05 if not

    .scale_rate[0]  .scale_rate[1]  .scale_rate[2]  .can_work[2]
    # first set allow working in "green"
        0.05           0.10            0.10           True
        0.05           0.15            0.15           True
        0.05           0.20            0.20           True
        0.05           0.25            0.25           True
        0.10           0.15            0.15           True
        0.10           0.20            0.20           True
        0.10           0.25            0.25           True
        0.15           0.20            0.20           True
        0.15           0.25            0.25           True
        0.20           0.25            0.25           True

    # second set prevent working in "green"
        0.05           0.10            0.15           False
        0.05           0.15            0.20           False
        0.05           0.20            0.25           False
        0.05           0.25            0.30           False
        0.10           0.15            0.20           False
        0.10           0.20            0.25           False
        0.10           0.25            0.30           False
        0.15           0.20            0.25           False
        0.15           0.25            0.30           False
        0.20           0.25            0.30           False

.. note::
  Note that we have added comments to this file using '#' - these
  are useful to help your future self understand what you were doing

Copy all of the files onto a cluster and submit the job where you
repeat each adjustable variable set 16 times. I used the PBS
job script;

.. code-block:: bash

    #!/bin/bash
    #PBS -l walltime=12:00:00
    #PBS -l select=4:ncpus=64:mem=64GB
    # The above sets 4 nodes with 64 cores each

    # source the version of metawards we want to use
    source $HOME/envs/metawards-0.8.0/bin/activate

    # change into the directory from which this job was submitted
    cd $PBS_O_WORKDIR

    metawards --additional ExtraSeedsLondon.dat \
            --disease lurgy3 -u lockdown.inp \
            --iterator lockdown \
            --input scan.csv --repeats 16 --nthreads 8 \
            --force-overwrite-output

Submit your job (e.g. ``qsub jobscript.sh``) and then wait for it to
finish. Once it has completed, generate the overview and average
graphs via;

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2
   metawards-plot --animate output/overview*.jpg
   metawards-plot --animate output/average*.jpg

What do you see?

I get a range of scenarios, from outbreaks that are controlled until
they die out, through oscillating outbreaks where the population is
forever moved between the "green" and "yellow" lockdown states,
through to outbreaks that grow despite lockdown. These can all
be seen here;

.. image:: ../../images/tutorial_3_6.gif
   :alt: Overview image of a lockdown with custom parameters
=================================
Responding to changing conditions
=================================

We've now created a lock-down advance function, but are currently
triggering this function in an iterator based on a fixed number
of days since the outbreak started.

A better approach would be to trigger the lock-down based on the
number of individuals who are detected as infected in the model.

To do this, edit your ``lockdown.py`` script and copy in the following;

.. code-block:: python

  from metawards.iterators import iterate_working_week, \
                                  advance_infprob, \
                                  advance_fixed, \
                                  advance_play
  from metawards.utils import Console

  def advance_lockdown(**kwargs):
      Console.debug("We are on lockdown")
      advance_infprob(scale_rate=0.25, **kwargs)
      advance_play(**kwargs)

  def iterate_lockdown(population, **kwargs):
      if not hasattr(population, "lockdown_state"):
          population.lockdown_state = "before"
          population.is_locked_down = False

      if population.lockdown_state == "before":
          if population.total > 5000:
              population.lockdown_state = "lockdown"
              population.lockdown_started = population.day
              population.is_locked_down = True

      if population.is_locked_down:
          return [advance_lockdown]
      else:
          return iterate_working_week(population=population,
                                      **kwargs)

The first thing we do here is see if the ``population`` has a
``lockdown_state`` variable using the standard Python
`hasattr function <https://docs.python.org/3/library/functions.html#hasattr>`__.
This variable won't exist on the first call to
to ``iterate_lockdown``, and so here we
set the ``lockdown_state`` to ``before``
and set the flag ``population.is_locked_down`` to ``False``.

Next, we check if the lockdown is in the ``before`` state. If it is,
then if the total infected population is greater than 5000 we change
the ``lockdown_state`` to ``lockdown``, save the day the lockdown
started to ``population.lockdown_started``, and set the flag
``population.is_locked_down`` to ``True``.

Finally, we either return our ``advance_lockdown`` advance function,
or the standard advance functions for a working week depending
on the value of the ``population.is_locked_down`` flag.

Run the model and draw the overview graph using;

.. code-block:: bash

  metawards -d lurgy3 --additional ExtraSeedsLondon.dat --iterator lockdown
  metawards-plot -i output/results.csv.bz2 --format jpg --dpi 150

You should now see that the lockdown takes effect some time after the
infected population grows above 5000. This tips the curve and reduces
the spread of the disease. You can see what my graphs looked like here;

.. image:: ../../images/tutorial_3_4_1_overview.jpg
   :alt: Overview image of a automatic lockdown

Releasing lockdown
------------------

We can use the data in ``population`` to decide when to release the
lockdown as well. For example, we could release when the size of the
infected population drops below 2000. To do this, edit your ``lockdown.py``
file to read;

.. code-block:: python

  from metawards.iterators import iterate_working_week, \
                                  advance_infprob, \
                                  advance_fixed, \
                                  advance_play
  from metawards.utils import Console

  def advance_lockdown(**kwargs):
      Console.debug("We are on lockdown")
      advance_infprob(scale_rate=0.25, **kwargs)
      advance_play(**kwargs)

  def iterate_lockdown(population, **kwargs):
      if not hasattr(population, "lockdown_state"):
          population.lockdown_state = "before"
          population.is_locked_down = False

      if population.lockdown_state == "before":
          if population.total > 5000:
              population.lockdown_state = "lockdown"
              population.lockdown_started = population.day
              population.is_locked_down = True

      elif population.lockdown_state == "lockdown":
          if population.total < 2000:
              population.lockdown_state = "after"
              population.lockdown_ended = population.day
              population.is_locked_down = False

      if population.is_locked_down:
          return [advance_lockdown]
      else:
          return iterate_working_week(population=population,
                                      **kwargs)

Run the model as before and see what happens...

To start, the lockdown has worked and the number of infections has fallen,
with the number falling below 2000 on day 78. However, releasing the
lockdown completely leads to a rapid growth in the infection,
with over 85,000 infected three weeks after lockdown ended. This
is unsurprising, as there was still a lot of infected individuals
remaining once lockdown ended, and a large population that was
still susceptible to infection (as you can see from the print
of my run and the overview graph below).

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 76 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                                    We are on lockdown                        lockdown.py:9
    S: 56055907  E: 201  I: 2192  R: 23777  IW: 215  POPULATION: 56082077
    Number of infections: 2393

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 77 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                                    We are on lockdown                        lockdown.py:9
    S: 56055697  E: 223  I: 2144  R: 24013  IW: 202  POPULATION: 56082077
    Number of infections: 2367

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 78 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                                    We are on lockdown                        lockdown.py:9
    S: 56055515  E: 210  I: 2109  R: 24243  IW: 174  POPULATION: 56082077
    Number of infections: 2319

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 79 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    [15:41:25]                        We are on lockdown                        lockdown.py:9
    S: 56055304  E: 182  I: 2099  R: 24492  IW: 205  POPULATION: 56082077
    Number of infections: 2281

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 80 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                                    We are on lockdown                        lockdown.py:9
    S: 56055121  E: 211  I: 2055  R: 24690  IW: 175  POPULATION: 56082077
    Number of infections: 2266

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 81 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                                    We are on lockdown                        lockdown.py:9
    S: 56054939  E: 183  I: 2004  R: 24951  IW: 176  POPULATION: 56082077
    Number of infections: 2187

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 82 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                                    We are on lockdown                        lockdown.py:9
    S: 56054770  E: 182  I: 1964  R: 25161  IW: 166  POPULATION: 56082077
    Number of infections: 2146

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 83 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56053695  E: 169  I: 1923  R: 26290  IW: 916  POPULATION: 56082077
    Number of infections: 2092

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 84 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56052615  E: 1075  I: 1879  R: 26508  IW: 920  POPULATION: 56082077
    Number of infections: 2954

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 85 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56051552  E: 1080  I: 2730  R: 26715  IW: 908  POPULATION: 56082077
    Number of infections: 3810

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 86 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56050161  E: 1063  I: 3598  R: 27255  IW: 1187  POPULATION: 56082077
    Number of infections: 4661

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 87 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56049050  E: 1391  I: 4452  R: 27184  IW: 966  POPULATION: 56082077
    Number of infections: 5843

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 88 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56047600  E: 1111  I: 5589  R: 27777  IW: 1236  POPULATION: 56082077
    Number of infections: 6700

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 89 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56044865  E: 1450  I: 6352  R: 29410  IW: 2072  POPULATION: 56082077
    Number of infections: 7802

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 90 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56041834  E: 2735  I: 7386  R: 30122  IW: 2246  POPULATION: 56082077
    Number of infections: 10121

.. image:: ../../images/tutorial_3_4_4_overview.jpg
   :alt: Overview image of a automatically released lockdown

Note how the second wave of infection makes the initial wave almost invisible
in this graph. The only visible evidence is the small peak in the number
of infected wards (**IW**) plot.

Relaxing, not removing lockdown
-------------------------------

The problem is that we treated lockdown like a binary switch, and
immediately went back to normal once it was lifted.

Instead, we need to release the lockdown in stages. To model this,
edit your ``lockdown.py`` to contain the following.

.. code-block:: python

    from metawards.iterators import iterate_working_week, \
                                    advance_infprob, \
                                    advance_fixed, \
                                    advance_play
    from metawards.utils import Console

    def advance_lockdown(population, **kwargs):
        Console.debug("We are in lockdown",
                      variables=[population.lockdown_scale_rate])

        advance_infprob(population=population,
                        scale_rate=population.lockdown_scale_rate,
                        **kwargs)
        advance_play(population=population, **kwargs)

    def iterate_lockdown(population, **kwargs):
        try:
            population.lockdown_state
        except Exception:
            population.lockdown_state = "before"
            population.is_locked_down = False
            population.lockdown_scale_rate = 0.05

        if population.lockdown_state == "before":
            if population.total > 5000:
                population.lockdown_state = "lockdown"
                population.lockdown_started = population.day
                population.is_locked_down = True

        elif population.lockdown_state == "lockdown":
            if population.total < 2000:
                population.lockdown_state = "relaxed_lockdown"
                population.lockdown_ended = population.day
                population.lockdown_scale_rate = 0.10
                population.is_locked_down = True

        elif population.lockdown_state == "relaxed_lockdown":
            if population.total < 1000:
                population.lockdown_scale_rate = 0.20
            else:
                population.lockdown_scale_rate = 0.10

        if population.is_locked_down:
            return [advance_lockdown]
        else:
            return iterate_working_week(population=population,
                                        **kwargs)

In this code we have created a new lockdown state that we've called
``relaxed_lockdown``. This is entered when the number of infections
drops below 2000. In this state controls can be released that
correspond to now only halving the infection rate (``scale_rate``
is increased to 0.10 from 0.05 during the strong lockdown).
In the ``relaxed_lockdown`` state the infected population
is always checked. If it is below 1000 then the lockdown can be
relaxed even more, with the ``scale_rate`` increasing from 0.10
to 0.20. However, if the infected population rises above 1000,
then the lockdown is tightened and the ``scale_rate`` is lowered
again to 0.10.

Have a go at running using this iterator. What do you see? In my
case I see the model moving from lockdown (``scale_factor==0.05``),
through relaxed lockdown (``scale_factor==0.1``) to light
lockdown (``scale_factor==0.2``) during the outbreak, which
is brought under control. The overview plots are here;

.. image:: ../../images/tutorial_3_4_2_overview.jpg
   :alt: Overview image of a automatically relaxing lockdown

There is a small second peak as the lockdown is relaxed, but
this seems to be under control.

.. warning::
  Remember, we cannot read too much into single **model runs**
  as these are very stochastic simulations. We would need to
  run models many times and average before we could gain real
  insight.

Returning to work
-----------------

Because Python is dynamically typed, we can set whatever flags
or add whatever data we want to the ``population`` object that
we need (or indeed to any Python object).

Let's now add an extra flag that will be used by
``advance_lockdown`` to call ``advance_fixed`` if the lockdown
has been lifted sufficiently for people to return to work.
Copy the below into your ``lockdown.py`` file;

.. code-block:: python

    from metawards.iterators import iterate_working_week, \
                                    advance_infprob, \
                                    advance_fixed, \
                                    advance_play
    from metawards.utils import Console


    def advance_lockdown(population, **kwargs):
        Console.debug("We are in lockdown",
                      variables=[population.lockdown_scale_rate,
                                 population.is_work_locked_down])

        advance_infprob(population=population,
                        scale_rate=population.lockdown_scale_rate,
                        **kwargs)

        advance_play(population=population, **kwargs)

        if not population.is_work_locked_down:
            advance_fixed(population=population, **kwargs)


    def iterate_lockdown(population, **kwargs):
        if not hasattr(population, "lockdown_state"):
            population.lockdown_state = "before"
            population.is_locked_down = False
            population.lockdown_scale_rate = 0.25
            population.is_work_locked_down = False

        if population.lockdown_state == "before":
            if population.total > 5000:
                population.lockdown_state = "lockdown"
                population.lockdown_started = population.day
                population.is_locked_down = True
                population.is_work_locked_down = True

        elif population.lockdown_state == "lockdown":
            if population.total < 2000:
                population.lockdown_state = "relaxed_lockdown"
                population.lockdown_ended = population.day
                population.lockdown_scale_rate = 0.5
                population.is_locked_down = True
                population.is_work_locked_down = False

        elif population.lockdown_state == "relaxed_lockdown":
            population.is_work_locked_down = False

            if population.total < 1000:
                population.lockdown_scale_rate = 0.75
            else:
                population.lockdown_scale_rate = 0.5

        if population.is_locked_down:
            return [advance_lockdown]
        else:
            return iterate_working_week(population=population,
                                        **kwargs)

This is getting longer, but I hope you can see that all we have
added is a ``population.is_work_locked_down`` flag, plus some
extra code to flip this between ``True`` and ``False``. This flag
is read by ``advance_lockdown``, which calls ``advance_fixed``
if the flag is ``False``. We've also added a check to see if the
infected population rises above 5000 while in "relaxed lockdown",
and if it does, to re-enter full lockdown.

Run the model and plot the graphs. What do you see? Do you get
a graph similar to below?

.. image:: ../../images/tutorial_3_4_3_overview.jpg
   :alt: Overview image of a automatically relaxing lockdown
=====================
Alternative lockdowns
=====================

The flexibility of ``metawards`` means that there are multiple different
ways you could choose to model a lockdown.

To understand how, we must look at how ``metawards`` calculates the
*force of infection* (FOI) of each ward. The FOI of a ward is used
to calculate the infection rate, which determines the rate by which
individuals in a ward become infected.

Force of infection
-------------------

The *force of infection* (FOI) is calculated for each ward individually,
using the equation;

1. :math:`F(w) = S \times U(t) \times \sum_s [ C_s \beta_s N_s(w) ]`

where;

* :math:`F(w)` is the *force of infection*, :math:`F`, calculated for a specific ward, :math:`w`,
* :math:`S` is a constant scaling parameter, set via :meth:`Population.scale_uv <metawards.Population.scale_uv>`,
* :math:`U(t)` is a seasonal scaling function (UV) calculated for the specified day :math:`t`,
* :math:`\sum_s` is the sum over all disease stages, :math:`s`,
* :math:`C_s` is the ``contrib_foi`` disease parameter for stage :math:`s`,
* :math:`\beta_s` is the ``beta`` disease parameter for stage :math:`s`, and
* :math:`N_s(w)` is the number of infected individuals in ward :math:`w` in disease stage :math:`s`.

The seasonal scaling function, :math:`U(t)` is calculated via;

2. :math:`U(t) = 1.0 - \frac{U_v}{2.0} + U_v \frac{\text{cos}(2 \pi t / 365)}{2.0}`

where;

* :math:`U_v` is the ``UV`` parameter supplied by the user ``--UV`` command line argument, and
* :math:`t` is the number of days between the current date of the model run and the
  date at which seasonal spread is at a maximum, ``UV_max``, set via the user
  ``--UV-max`` command line argument.

The seasonal scaling function acts to scale down the force of infection
from a seasonal maximum (typically 1st January to model winter diseases in
the northern hemisphere) to a seasonal minimum 6 months later. :math:`U_v`
should be set to between 0 and 1, and is the scaling factor at the
six month minimum. The :math:`U(t)` functions for different values of
:math:`U_v` are shown below;

.. image:: ../../images/uv.jpg
   :alt: Plots of U(t) for different values of U_v

Metapopulation movements
------------------------

The force of infection calculation is based on the number of infected individuals
at each disease stage in each ward, :math:`N_s(w)`. The metapopulation
part of ``metawards`` is because individuals move between wards, and thus
will contribute to :math:`N_s(w)` differently depending on when and where
they are.

The force of infection is calculated both for *daytime* and *nighttime*.
Individuals contribute to the force of infection calculation for wards
they visit during the *daytime*, and to their home ward at *nighttime*.

Workers visit their fixed work ward during the *daytime*, as long as they
are not too sick to go to work (controlled by ``too_sick_to_move``, and
in which case they stay at home). Players
will visit a randomly chosen play ward during the *daytime*.

However, these movements are controlled by a cutoff parameter,
:meth:`Parameters.dyn_dist_cutoff <metawards.Parameters.dyn_dist_cutoff>`.
If the distance between the home and work ward, or home and play ward, is
greater than the cutoff (measured in kilometers), then the worker or player
will stay at home, and contribute to the force of infection of their home
ward.

.. note::

   This is all calculated in the :func:`~metawards.iterators.advance_foi`
   advance function. If you want to change this, you can write your
   own version of :func:`~metawards.iterators.advance_foi` and set that
   to be used instead via a custom iterator.

Rates of infection
------------------

The daytime and nighttime forces of infection for each ward are converted
into daytime and nighttime infection probabilities by the
:func:`~metawards.iterators.advance_infprob` advance function. These give
the probability that an individual staying in the ward that day or night
would be infected. These probablities are used by the
:func:`~metawards.iterators.advance_fixed` and
:func:`~metawards.iterators.advance_play` advance functions to determine
whether workers or players, as they move between wards each day, would
be infected.

Enacting lockdown
-----------------

Based on this knowledge, we could enact a lockdown by adjusting the
parameters used to calculate the force of infection, and the parameters
used to control movement of individuals between wards. For example,
we could reduce :meth:`Population.scale_uv <metawards.Population.scale_uv>`,
thereby reducing :math:`S` in equation 1. This would have the same effect
on scaling down the FOI as scaling down :math:`\beta`. This would mean that
you could relate different values of :math:`S` to different lockdown
control measures, e.g. closing schools, wearing masks, limiting number
of contacts etc.

We could also reduce
:meth:`Parameters.dyn_dist_cutoff <metawards.Parameters.dyn_dist_cutoff>`
to, e.g. 5 km, to prevent most work and play movements. Indeed, we could
even reduce this to 0 km to stop all movement between wards.

A good example of an
`alternative lockdown model is here <https://github.com/metawards/MetaWards/tree/devel/examples/lockdown>`__.
This is provided as an example in the MetaWards GitHub repository, and
enacts lockdown by directly changing these two parameters.
This has the effect of reducing the contribution from each infected
individual to the overall *force of infection* of each ward, and reducing
the movement of individuals between wards.

There are many parameters to adjust. You can also add these
to your scan to investigate their impact.
The full list of built-in adjustable parameters is below;

.. program-output:: python get_variableset_help.py
================
Custom variables
================

In the last example we got ``metawards`` to model a lockdown. However,
we put a lot of hard-coded parameters into ``iterate_lockdown`` function
we wrote. One of the main purposes of ``metawards`` is to enable us
to scan through multiple variables to see their affect on the
population trajectories.

User variables
--------------

``metawards`` supports the use of any custom user variables. For example,
create a new file called ``custom.inp`` and copy in the below;

::

    # first state
    user.scale[0] = 0.2
    user.flag[0]  = False

    # second state
    user.scale[1] = 0.5
    user.flag[1]  = True

    # third state
    user.scale[2] = 0.7
    user.flag[2]  = True

.. note::
  Comments start with a '#'. The file format is very flexible. You can
  use an '=' sign or a ':' or even nothing, e.g. "user.scale[1] 0.5".
  You can abbreviate "user.variable" to ".variable" too if you want.

This file is going to be used to define three states. The first (at index 0)
has a scale factor of 0.2, and a flag set to ``False``. The second has
a scale factor of 0.5 and a flag set to ``True``. The third has a
scale factor of 0.7 and a flag set also to ``True``.

We can pass this input file into ``metawards`` using the ``--user-variables``
(or ``-u``) parameter, e.g. try typing;

.. code-block:: bash

   metawards -d lurgy3 --user-variables custom.inp

In the output you should see a line that reads;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Custom parameters and seeds ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Adjusting variables to (.scale[0]=0.2, .flag[0]=False, .scale[1]=0.5, .flag[1]=True,
    .scale[2]=0.7, .flag[2]=True)[repeat 1]

.. note::
  See how "user.variable" has been abbreviated to ".variable" in this output

This shows that ``metawards`` has read your custom parameters into the
program. They are stored as ``user_params`` in a
:class:`~metawards.Parameters` object that is stored in the
:class:`~metawards.Network` that holds the network of wards being modelled.

Reading user variables in custom iterators
------------------------------------------

You can access these parameters in an iterator or advance function using
``network.params.user_params[X]`` when ``X`` is the name of the parameter
you want. For example, ``network.params.user_params["scale[0]"] returns
``0.2``. For example, copy the below into a file called ``custom.py``;

.. code-block:: python

    from metawards.utils import Console

    def get_state(population):
        if population.day < 2:
            return 0
        elif population.day < 4:
            return 1
        else:
            return 2

    def advance_lockdown(network, population, **kwargs):
        params = network.params
        state = get_state(population)
        scale = params.user_params["scale"][state]
        Console.debug("Hello advance_lockdown", variables=[scale])

    def advance_relaxed(network, population, **kwargs):
        params = network.params
        state = get_state(population)
        scale = params.user_params["scale"][state]
        Console.debug("Hello advance_relaxed", variables=[scale])

    def iterate_custom(network, population, **kwargs):
        params = network.params
        state = get_state(population)
        flag = params.user_params["flag"][state]
        Console.debug("Hello iterate_custom", variables=[scale, flag])

        if flag:
            return [advance_lockdown]
        else:
            return [advance_relaxed]

This code defines four functions:

1. ``get_state`` - this returns the state that the population should be
   set to. In this case, the state is 0, 1 or 2 depending on the day
   of the outbreak.

2. ``iterate_custom`` - this gets the state using ``get_state``. It then
   looks up the ``flag`` custom parameter at index ``state``.
   If the flag is True, then it returns ``advance_lockdown``. Otherwise
   it returns ``advance_relaxed``.

3. ``advance_lockdown`` - this gets the state using ``get_state``. It then
   looks up the ``scale`` custom parameter at index ``state``.
   It prints this to the screen.

4. ``advance_relaxed`` - this does the same as ``advance_lockdown``, but
   prints a different message to the screen.

Use this iterator by running ``metawards`` via;

.. code-block:: bash

   metawards -d lurgy3 -u custom.inp --iterator custom --debug

You should now see printed to the screen something very similar to the below;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    [15:56:50]                       Hello iterate_custom                        custom.py:31

      Name │ Value
    ═══════╪═══════
     state │ 0
      flag │ False

                                    Hello advance_relaxed                       custom.py:24

      Name │ Value
    ═══════╪═══════
     scale │ 0.2

    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                                    Hello iterate_custom                        custom.py:31

      Name │ Value
    ═══════╪═══════
     state │ 0
      flag │ False

                                    Hello advance_relaxed                       custom.py:24

      Name │ Value
    ═══════╪═══════
     scale │ 0.2

    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                                    Hello iterate_custom                        custom.py:31

      Name │ Value
    ═══════╪═══════
     state │ 1
      flag │ True

                                    Hello advance_lockdown                       custom.py:17

      Name │ Value
    ═══════╪═══════
     scale │ 0.5

    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                                    Hello iterate_custom                        custom.py:31

      Name │ Value
    ═══════╪═══════
     state │ 1
      flag │ True

                                    Hello advance_lockdown                       custom.py:17

      Name │ Value
    ═══════╪═══════
     scale │ 0.5

    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                                    Hello iterate_custom                        custom.py:31

      Name │ Value
    ═══════╪═══════
     state │ 2
      flag │ True

    [15:56:51]                      Hello advance_lockdown                       custom.py:17

      Name │ Value
    ═══════╪═══════
     scale │ 0.7

    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0
    Infection died ... Ending on day 5


Hopefully the change between state, functions called and values of the
scale factor printed makes sense and follows what you expected.

Scanning custom variables
-------------------------

You can also adjust your custom variables by scanning in the same way
that we adjusted in-built variables like **beta** and **progress**
in :doc:`an earlier part of this tutorial <../part02/02_adjustable>`.

In this case, we use the full (``user.variable``) or abbreviated
(``.variable``) names as titles in the ``metawards`` input file.

For example, create a new file called ``scan.csv`` and copy in the below;

::

    .scale[0]  .flag[0]
      0.1        False
      0.2        False
      0.3        False
      0.1         True
      0.2         True
      0.3         True

This tells ``metawards`` to perform six *model runs*, with ``user.scale[0]``
varied from 0.1-0.3, and ``user.flag[0]`` varied from "True" to "False".

Perform these model runs using;

.. code-block:: bash

   metawards -d lurgy3  -u custom.inp --iterator custom -i scan.csv

You should get output that is very similar to this;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ MULTIPROCESSING ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 1 of 6                                                                  │
    │  (.scale[0]=0.1, .flag[0]=False)[repeat 1]                                             │
    │  2020-05-25: DAY: 5 S: 56082077    E: 0    I: 0    R: 0    IW: 0   UV: 1.0   TOTAL     │
    │  POPULATION 56082077                                                                   │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 2 of 6                                                                  │
    │  (.scale[0]=0.2, .flag[0]=False)[repeat 1]                                             │
    │  2020-05-25: DAY: 5 S: 56082077    E: 0    I: 0    R: 0    IW: 0   UV: 1.0   TOTAL     │
    │  POPULATION 56082077                                                                   │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 3 of 6                                                                  │
    │  (.scale[0]=0.3, .flag[0]=False)[repeat 1]                                             │
    │  2020-05-25: DAY: 5 S: 56082077    E: 0    I: 0    R: 0    IW: 0   UV: 1.0   TOTAL     │
    │  POPULATION 56082077                                                                   │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 4 of 6                                                                  │
    │  (.scale[0]=0.1, .flag[0]=True)[repeat 1]                                              │
    │  2020-05-25: DAY: 5 S: 56082077    E: 0    I: 0    R: 0    IW: 0   UV: 1.0   TOTAL     │
    │  POPULATION 56082077                                                                   │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 5 of 6                                                                  │
    │  (.scale[0]=0.2, .flag[0]=True)[repeat 1]                                              │
    │  2020-05-25: DAY: 5 S: 56082077    E: 0    I: 0    R: 0    IW: 0   UV: 1.0   TOTAL     │
    │  POPULATION 56082077                                                                   │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 6 of 6                                                                  │
    │  (.scale[0]=0.3, .flag[0]=True)[repeat 1]                                              │
    │  2020-05-25: DAY: 5 S: 56082077    E: 0    I: 0    R: 0    IW: 0   UV: 1.0   TOTAL     │
    │  POPULATION 56082077                                                                   │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Writing a summary of all results into the csv file                                    │
    │  /Users/chris/GitHub/tutorial/weekend/output/results.csv.bz2. You can use this to      │
    │  quickly look at statistics across all runs using e.g. R or pandas                     │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘

A quick look in the ``output/results.csv.bz2`` file, e.g. using pandas,
shows that the fingerprint and columns for custom variables are
constructed identially to in-built variables, e.g.

.. code-block:: python

    >>> import pandas as pd
    >>> df = pd.read_csv("output/results.csv.bz2")
    >>> print(df)
        ingerprint  repeat  .scale[0]  .flag[0]  day        date         S  E  I  R  IW   UV
    0        0i1vF       1        0.1     False    0  2020-05-20  56082077  0  0  0   0  1.0
    1        0i1vF       1        0.1     False    1  2020-05-21  56082077  0  0  0   0  1.0
    2        0i1vF       1        0.1     False    2  2020-05-22  56082077  0  0  0   0  1.0
    3        0i1vF       1        0.1     False    3  2020-05-23  56082077  0  0  0   0  1.0
    4        0i1vF       1        0.1     False    4  2020-05-24  56082077  0  0  0   0  1.0
    5        0i1vF       1        0.1     False    5  2020-05-25  56082077  0  0  0   0  1.0
    6        0i2vF       1        0.2     False    0  2020-05-20  56082077  0  0  0   0  1.0
    7        0i2vF       1        0.2     False    1  2020-05-21  56082077  0  0  0   0  1.0
    8        0i2vF       1        0.2     False    2  2020-05-22  56082077  0  0  0   0  1.0
    9        0i2vF       1        0.2     False    3  2020-05-23  56082077  0  0  0   0  1.0
    10       0i2vF       1        0.2     False    4  2020-05-24  56082077  0  0  0   0  1.0
    11       0i2vF       1        0.2     False    5  2020-05-25  56082077  0  0  0   0  1.0
    12       0i3vF       1        0.3     False    0  2020-05-20  56082077  0  0  0   0  1.0
    13       0i3vF       1        0.3     False    1  2020-05-21  56082077  0  0  0   0  1.0
    14       0i3vF       1        0.3     False    2  2020-05-22  56082077  0  0  0   0  1.0
    15       0i3vF       1        0.3     False    3  2020-05-23  56082077  0  0  0   0  1.0
    16       0i3vF       1        0.3     False    4  2020-05-24  56082077  0  0  0   0  1.0
    17       0i3vF       1        0.3     False    5  2020-05-25  56082077  0  0  0   0  1.0
    18       0i1vT       1        0.1      True    0  2020-05-20  56082077  0  0  0   0  1.0
    19       0i1vT       1        0.1      True    1  2020-05-21  56082077  0  0  0   0  1.0
    20       0i1vT       1        0.1      True    2  2020-05-22  56082077  0  0  0   0  1.0
    21       0i1vT       1        0.1      True    3  2020-05-23  56082077  0  0  0   0  1.0
    22       0i1vT       1        0.1      True    4  2020-05-24  56082077  0  0  0   0  1.0
    23       0i1vT       1        0.1      True    5  2020-05-25  56082077  0  0  0   0  1.0
    24       0i2vT       1        0.2      True    0  2020-05-20  56082077  0  0  0   0  1.0
    25       0i2vT       1        0.2      True    1  2020-05-21  56082077  0  0  0   0  1.0
    26       0i2vT       1        0.2      True    2  2020-05-22  56082077  0  0  0   0  1.0
    27       0i2vT       1        0.2      True    3  2020-05-23  56082077  0  0  0   0  1.0
    28       0i2vT       1        0.2      True    4  2020-05-24  56082077  0  0  0   0  1.0
    29       0i2vT       1        0.2      True    5  2020-05-25  56082077  0  0  0   0  1.0
    30       0i3vT       1        0.3      True    0  2020-05-20  56082077  0  0  0   0  1.0
    31       0i3vT       1        0.3      True    1  2020-05-21  56082077  0  0  0   0  1.0
    32       0i3vT       1        0.3      True    2  2020-05-22  56082077  0  0  0   0  1.0
    33       0i3vT       1        0.3      True    3  2020-05-23  56082077  0  0  0   0  1.0
    34       0i3vT       1        0.3      True    4  2020-05-24  56082077  0  0  0   0  1.0
    35       0i3vT       1        0.3      True    5  2020-05-25  56082077  0  0  0   0  1.0================
Plotting outputs
================

You should now have a ``results.csv.bz2`` file in the ``output`` directory,
which contains the results of four *model runs* of the outbreak of the
lurgy that was seeded in London.

You can plot graphs of the result using the
:doc:`metawards-plot <../../metawards_plot_help>` command.
To run this, type;

.. code-block:: bash

   metawards-plot --input output/results.csv.bz2

.. note::
   ``metawards-plot`` uses `pandas <https://pandas.pydata.org>`__ and
   `matplotlib <https://matplotlib.org>`__ for plotting. If you don't have
   these on your computer then you will see an error message giving
   instructions on how to install the packages. Note that the default
   format of the output is jpeg. You can change the format using the
   ``--format`` option, e.g. ``--format png`` or ``--format pdf``.
   You may need to install
   `Pillow <https://pillow.readthedocs.io/en/stable/>`__
   to support output in some file formats.

Understanding the graphs
------------------------

This will create two sets of graphs;

* ``output/overview.pdf``
    This is an overview of the **E**, **I**, **IW** and **R** values
    from each day of the model outbreak for each of the four *model runs*.
    Your graph should look something like this;

.. image:: ../../images/tutorial_1_3_overview.jpg
   :alt: Overview image of the outbreak of the lurgy

* ``output/average.pdf``
    This shows the average of the **E**, **I**, **IW** and **R** values,
    with the standard deviation shown as the error bars. Your graph should
    look something like this;

.. image:: ../../images/tutorial_1_3_average.jpg
   :alt: Overview image of the outbreak of the lurgy

The **E** graph shows the total number of latent infections. It should be
slightly ahead of and of a similar shape to the **I** graph, which shows
the total number of infections.

The **IW** graph shows the number of wards with at least one infection.
This cannot grow to more than the total number of wards, hence why
you see this graph topping out at ~8588, as this is the number of wards.

The **R** graph shows the number of individuals removed from the epidemic
(e.g. as they may have recovered). This should have an "S" shape,
showing exponential growth in the initial stage of the epidemic that
tails off as the number of individuals susceptible to infection is
reduced (e.g. as immunity in the population is built up).

.. note::

    Your graphs may look a little different in the exact numbers, but should
    be similar in shape. The purpose of this type of modelling is not to
    make exact numerical predictions, but to instead understand trends
    and timelines.

Jupyter notebooks
-----------------

The ``metawards-plot`` command can be used for quick-and-simple plots.
If you want something more complex then you can take a look at the functions
in the `metawards.analysis` module. In addition, we have a
:download:`Jupyter notebook <../../notebooks/1_3_plotting.ipynb>`
which you can look at which breaks down exactly how ``metawards-plot``
uses pandas and matplotlib to render these two graphs.
===============
Getting Started
===============

This tutorial assumes that you have installed ``metawards``. To check
that you have, type the following into your console;

.. code-block:: bash

   metawards --version

and then press return. You should see something similar to the below
printed to your screen.

.. command-output:: metawards --version

If you don't see this, or the output includes a warning about not being
about to find `MetaWardsData`, then please try
:doc:`installing MetaWards <../../install>` or
:doc:`installing and configuring MetaWardsData <../../model_data>` again.

.. warning::

  This tutorial is written for ``metawards`` version |MetaWardsVersion| or
  higher. If you are using an older version then please upgrade.

Introducing the Lurgy
---------------------

`The Lurgy <https://en.wiktionary.org/wiki/lurgy>`__ is a
**completely ficticious** disease that we will use throughout this
tutorial. We can run a simulation of an outbreak of the lurgy using
the ``--disease`` command line argument. Type the following;

.. code-block:: bash

   metawards --disease lurgy

Press return and you should see a lot of output printed. Near the end
of the output you will see these lines;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0
    Infection died ... Ending on day 5

.. note::
   Do not worry if you don't see exactly this output. You may be using
   a different version of ``metawards`` compared to the one used to write
   this tutorial. The main thing to look for is the line
   ``Infection died ... Ending on day 5``

The ``--disease`` option also has the shorthand ``-d``. You can get the same
output as above by typing;

.. code-block:: bash

   metawards -d lurgy

.. note::
   This time when you ran ``metawards`` it stopped to say that the output
   directory already exists, and if you want to remove it.

The ``metawards`` program takes care not to overwrite any of your output.
By default a lot of output files from this run have been written to a
directory called ``output`` (we will take a look at these files later).
``metawards`` will ask you if you want to remove any existing output.
Press ``y`` and hit return to do so. If you want to automatically
remove existing output then use the ``--force-overwrite-output`` option,
e.g.

.. code-block:: bash

   metawards -d lurgy --force-overwrite-output

You can also set the output directory using the ``--output`` or ``-o``
options, e.g.

.. code-block:: bash

   metawards -d lurgy -o output2

Seeding an outbreak
-------------------

The key output from ``metawards`` are the lines which read;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0
    Infection died ... Ending on day 5

These tell you how long the outbreak lasted (in this case, 5 days),
together with how many people were infected. These are the numbers next
to the codes;

* **S**: The number of the population who are *susceptible* to infection
* **E**: The number of the population who are *latent*, meaning they are
  infected, but not yet infectious.
* **I**: The number of the population who are *infected*, meaning they
  have symptoms and are infectious
* **R**: The number of the population who are removed from being susceptible,
  either because they have been newly infected that day, or because they
  have recovered from the infection and are no longer susceptible to infection
* **IW**: The number of electoral wards that contain at least one
  individual who was newly infected that day.

For more information about these values, please
`read <https://doi.org/10.1016/j.epidem.2009.11.002>`__
`the <https://doi.org/10.1073/pnas.1000416107>`__
`papers <https://doi.org/10.1101/2020.02.12.20022566>`__ detailed
in the :doc:`scientific background <../../index>`.

From this output it is clear that no-one has been infected by the lurgy.
This is because we haven't yet seeded any outbreaks. We can seed an
outbreak in a specific electoral ward by using an additional seeds file.
In this case, we will seed an infection of the lurgy in London using
the `ExtraSeedsLondon.dat <https://github.com/metawards/MetaWardsData/blob/main/extra_seeds/ExtraSeedsLondon.dat>`__
file that comes in ``MetaWardsData``. You specify the additional seeds
file to use via the ``--additional`` or ``-a`` options.

Try typing the below into your console and press return;

.. code-block:: bash

   metawards -d lurgy -a ExtraSeedsLondon.dat

Now the program will run for a long time (minutes), and you will see
an outbreak move through the population. The final lines of your output
may look something like this;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 214 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 11780863  E: 0  I: 1  R: 44301213  IW: 1  POPULATION: 56082077
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 215 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 11780863  E: 1  I: 0  R: 44301213  IW: 0  POPULATION: 56082077
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 216 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 11780863  E: 0  I: 1  R: 44301213  IW: 0  POPULATION: 56082077
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 217 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 11780863  E: 0  I: 1  R: 44301213  IW: 0  POPULATION: 56082077
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 218 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 11780863  E: 0  I: 0  R: 44301214  IW: 0  POPULATION: 56082077
    Number of infections: 0
    Infection died ... Ending on day 219

.. note::

   Do not worry if your numbers are different. All will be explained :-)

Repeating a calculation
-----------------------

``metawards`` runs a stochastic simulation. This means that random numbers
are used in the decisions on how individuals in the model are infected,
and how quickly they progress through the infection. This means that
every ``metawards`` run is different.

Fortunately, ``metawards`` prints enough information in the output
to enable a job to be repeated. Look the for line the reads;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Repeating this run ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    To repeat this job use the command;
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │metawards –repeats 1 –seed 85564100 –additional ExtraSeedsLondon.dat –output output     │
    │–disease lurgy –start-date 2020-05-20 –start-day 0 –parameters march29 –repository      │
    │/Users/chris/GitHub/MetaWardsData –population 57104043 –nsteps 730 –UV 0.0 –nthreads 4  │
    │–nprocs 1 –max-nodes 16384 –max-links 4194304                                           │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Or alternatively use the config.yaml file that will be written to the output directory and
    use the command;
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │metawards -c config.yaml                                                                │
    └────────────────────────────────────────────────────────────────────────────────────────┘

This is the command line that you can use to repeat a job (note that
the command line you see will be different). We have been careful to
write ``metawards`` so that it gives the same output when you use
the same inputs, using the same version of ``metawards`` and same version
of data in ``MetaWardsData``, for the same random number seed and running
the calculation over the same number of threads. We consider it a bug
if ``metawards`` is not reproducible, and ask that you
`submit an issue <https://github.com/metawards/MetaWards/issues>`__ if
you find you cannot repeat a run.

As the command line can be quite long, ``metawards`` will also print out
a config file in the ``output`` directory called ``output/config.yaml``.
This file contains everything needed to reproduce the calculation, which
can be re-run using the command;

.. code-block:: bash

   metawards -c config.yaml

(assuming you have copied the ``config.yaml`` file into your current
directory)

Using config files for inputs
-----------------------------

You can use this config file directly to run a job using the ``--config``
or ``-c`` options, e.g.

.. code-block:: bash

   metawards --config config.yaml

This should repeat the calculation that generated this config. You can
also edit this file and use it to store commonly used options, e.g.
if you always want to model the lurgy, then the config file would read;

::

  disease: lurgy

and you could use this via

.. code-block:: bash

   metawards -c config.yaml

.. note::
  ``metawards`` uses the `ConfigArgParse <https://pypi.org/project/ConfigArgParse/>`__
  python module for parsing command line arguments. Options can be passed
  on the command line, in a yaml or ini format config file, or in
  some identified cases as an environment variable. If an arg is
  specified in more than one place, then commandline values override
  environment variables which override config file values which
  override defaults.
===================
Multiple model runs
===================

In the last page you successfully performed a single run modelling
an outbreak of the lurgy that started in London. This run
(which we call a *model run*) is stochastic, meaning that the results
will be slightly different every time it is performed.

To gain confidence in any predictions, we need to perform a *model run*
multiple times, and average over the results.

Performing multiple model runs
------------------------------

``metawards`` has the command line option ``--repeats`` (or ``-r``) to
set the number of times a ``model run`` should be repeated. For example,
run the below command to repeat the *model run* four times;

.. code-block:: bash

  metawards -d lurgy -a ExtraSeedsLondon.dat --repeats 4

``metawards`` will automatically use as many of the cores on your computer
as it can to parallelise the jobs. On my computer, the output shows;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Running the model ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Using random number seed: 87340504
    Running 4 jobs using 4 process(es)

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ MULTIPROCESSING ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

I have four processor cores on my laptop, so I see the four repeats run
in parallel using four processes, with each *model run* performed
using 1 thread. You will see a different distribution of threads
and processes if you have a different number of cores on your computer.
You can set the number of processes that ``metawards`` should use via
the ``--nprocs`` command line option. You can set the number of threads
that ``metawards`` should use via the ``--nthreads`` command line option.

This calculation may take some time (2-3 minutes). This time, instead
of seeing a summary of the outbreak, ``metawards`` will show a summary
of the different *model run* jobs. Something similar to this should
be printed;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ MULTIPROCESSING ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 1 of 4                                                                  │
    │  (NO_CHANGE)[repeat 1]                                                                 │
    │  2021-01-13: DAY: 238 S: 11784852    E: 0    I: 0    R: 44297225    IW: 0   UV: 1.0    │
    │  TOTAL POPULATION 56082077                                                             │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 2 of 4                                                                  │
    │  (NO_CHANGE)[repeat 2]                                                                 │
    │  2021-01-05: DAY: 230 S: 11770162    E: 0    I: 0    R: 44311915    IW: 1   UV: 1.0    │
    │  TOTAL POPULATION 56082077                                                             │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 3 of 4                                                                  │
    │  (NO_CHANGE)[repeat 3]                                                                 │
    │  2021-02-05: DAY: 261 S: 11789449    E: 0    I: 0    R: 44292628    IW: 1   UV: 1.0    │
    │  TOTAL POPULATION 56082077                                                             │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 4 of 4                                                                  │
    │  (NO_CHANGE)[repeat 4]                                                                 │
    │  2021-01-04: DAY: 229 S: 11779688    E: 0    I: 0    R: 44302389    IW: 0   UV: 1.0    │
    │  TOTAL POPULATION 56082077                                                             │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Writing a summary of all results into the csv file                                    │
    │  /Users/chris/GitHub/tutorial/test/output/results.csv.bz2. You can use this to         │
    │  quickly look at statistics across all runs using e.g. R or pandas                     │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘


.. note::
   ``metawards`` prints a progress spinner to the screen while the jobs
   are running, so show you that it hasn't crashed. You can switch off
   the spinner using the ``--no-spinner`` option if it annoys you.
   Similarly, ``metawards`` by default writes output to the screen using
   a very colourful theme. You can change this to a more simple and
   less-colourful theme by passing in the option ``--theme simple``.

In this case, all four outbreaks completed within 229-261 days, while the
number of the population who progressed to the '**R**' state were all
around 44.3 million.

The results.csv.bz2 file
------------------------

The day-by-day progress of each the outbreak for each *model run* is
recorded in the output file ``results.csv.bz2``. This is a comma-separated
file that has been compressed using `bzip2 <https://en.wikipedia.org/wiki/Bzip2>`__.

You can read this file easily using
`Python Pandas <https://pandas.pydata.org>`__  or with
`R <https://www.r-project.org>`__. You can even import this into Excel
(although you may need to uncompress this file first using ``bunzip2``).

For example, if you have Pandas installed, then you can read this
file via an `ipython <https://ipython.org>`__ or
`Jupyter notebook <https://jupyter.org>`__ session via;

.. code-block:: python

  >>> import pandas as pd
  >>> df = pd.read_csv("output/results.csv.bz2")
  >>> df
      fingerprint  repeat  day        date         S  E  I         R  IW
  0     NO_CHANGE       1    0  2020-04-20  56082077  0  0         0   0
  1     NO_CHANGE       1    1  2020-04-21  56082077  0  0         0   0
  2     NO_CHANGE       1    2  2020-04-22  56082072  5  0         0   0
  3     NO_CHANGE       1    3  2020-04-23  56082072  0  5         0   0
  4     NO_CHANGE       1    4  2020-04-24  56082068  0  5         4   4
  ..          ...     ...  ...         ...       ... .. ..       ...  ..
  929   NO_CHANGE       4  224  2020-11-30  11782419  0  4  44299654   0
  930   NO_CHANGE       4  225  2020-12-01  11782419  0  3  44299655   0
  931   NO_CHANGE       4  226  2020-12-02  11782419  0  1  44299657   0
  932   NO_CHANGE       4  227  2020-12-03  11782419  0  1  44299657   0
  933   NO_CHANGE       4  228  2020-12-04  11782418  0  0  44299659   1

  [934 rows x 9 columns]

Each repeat is given its own number, which is in the ``repeat`` column.
The day of the outbreak is given in the ``day`` column. This counts up
from *day zero* when the outbreak started, to the last day when the
outbreak was over. You can control the start day of the outbreak using
the ``--start-day`` command line option.

The ``date`` column contains the date of each day in the outbreak. By
default, ``metawards`` assumes that *day zero* is today. You can set the
date of *day zero* using the ``--start-date`` command line option, e.g.
``--start-date tomorrow`` would start tomorrow, while
``--start-date Jan 1`` would start on January 1st this year.

The values of **S**, **E**, **I**, **R** and **IW** for each repeat for
each day are then given in their correspondingly named columns.

The *fingerprint* column not used for this calculation - we will see what it is
later.
======================
Using a Custom Network
======================

You can run a simulation using a custom network by passing filename of
the JSON file that contains the network to ``metawards`` via the
``--model`` or ``-m`` parameter.

For example, to use the ``custom_network.json.bz2`` file from the last section,
together with the ``lurgy4.json`` disease model from previous chapters,
and seed the outbreak with 5 infections in London on day 1
you would run;

.. code-block:: bash

   metawards -d lurgy4 -m custom_network.json.bz2 -a "1 5 London"

You should see that the model runs very quickly, producing output similar
to;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Loading additional seeds from the command line
    ┏━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
    ┃ Day ┃ Demographic ┃                                   Ward                                    ┃   Number   ┃
    ┃     ┃             ┃                                                                           ┃   seeded   ┃
    ┡━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
    │  1  │    None     │ 2 : WardInfo(name='London', alternate_names=, code='', alternate_codes=,  │     5      │
    │     │             │        authority='', authority_code='', region='', region_code='')        │            │
    └─────┴─────────────┴───────────────────────────────────────────────────────────────────────────┴────────────┘
    seeding play_infections[0][2] += 5
    S: 20345  E: 5  I: 0  R: 0  IW: 0  POPULATION: 20350
    Number of infections: 5

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 20345  E: 0  I: 5  R: 0  IW: 0  POPULATION: 20350
    Number of infections: 5

    ...
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 129 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 2895  E: 0  I: 1  R: 17454  IW: 0  POPULATION: 20350
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 130 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 2895  E: 0  I: 1  R: 17454  IW: 0  POPULATION: 20350
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 131 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 2895  E: 0  I: 1  R: 17454  IW: 0  POPULATION: 20350
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 132 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 2895  E: 0  I: 0  R: 17455  IW: 0  POPULATION: 20350
    Number of infections: 0
    Ending on day 132

Running MetaWards from within Python or R
-----------------------------------------

It is also possible to run your custom network by passing it in directly
to the :func:`metawards.run` function in Python or R. For example,
in Python;

.. code-block:: python

   >>> from metawards import Ward, run
   >>> bristol = Ward(name="Bristol")
   >>> bristol.add_workers(500)
   >>> bristol.set_num_players(750)
   >>> london = Ward(name="London")
   >>> london.add_workers(8500)
   >>> london.set_num_players(10000)
   >>> bristol.add_workers(500, destination=london)
   >>> london.add_workers(100, destination=bristol)
   >>> wards = bristol + london
   >>> print(wards)
   [ Ward( info=Bristol, id=1, num_workers=1000, num_players=750 ), Ward( info=London, id=2, num_workers=8600, num_players=10000 ) ]
   >>> results = run(model=wards, additional=5)

This would create the wards, and then run the model. This will run in a
new directory called ``output_XXXX`` (where XXXX is replaced by a random
string). The ``results`` variable holds the full path to the resulting
``results.csv.bz2`` file for this run. The arguments to
:func:`metawards.run` match those of the command line program. Any Python
objects (e.g. the wards, disease, demographics) can be passed in as
Python objects. They will be automatically converted to JSON files and
passed to the ``metawards`` processed in the background.

.. note::

   You can use the ``+`` operator to add multiple individual ward objects
   together to create the wards, e.g. ``wards = bristol + london``.

.. note::

   You can force :func:`metawards.run` to use a specified output directory
   by passing in the ``output`` argument. You will need to set
   ``force_overwrite_output`` to True to overwrite any existing output.
   You can silence the printing to the screen by passing in
   ``silent = True``.

You can achieve the same in R by typing;

.. code-block:: R

   > library(metawards)
   > bristol <- metawards$Ward(name="Bristol")
   > bristol$add_workers(500)
   > bristol$set_num_players(750)
   > london <- metawards$Ward(name="London")
   > london$add_workers(8500)
   > london$set_num_players(10000)
   > bristol$add_workers(500, destination=london)
   > london$add_workers(100, destination=bristol)
   > wards = metawards$Wards()
   > wards$add(bristol)
   > wards$add(london)
   > print(wards)
   [ Ward( info=Bristol, id=1, num_workers=1000, num_players=750 ), Ward( info=London, id=2, num_workers=8600, num_players=10000 ) ]
   > results <- metawards$run(model=wards, additional=5)

.. note::

   R does not support adding individual ward objects together to get Wards

======================
Creating Networks in R
======================

While ``metawards`` is a Python module, you can use the :mod:`metawards`
module directly in R.

This is because the `reticulate project <https://rstudio.github.io/reticulate/>`__
lets you embed and use Python directly in your R scripts.

Installing MetaWards in R
-------------------------

You can install MetaWards by starting R and typing;

.. code-block:: R

    > library(devtools)
    > install_github("metawards/rpkg")

This will install the MetaWards R package.

Next, you need to install MetaWards itself. The R package provides
a convenient function to support this. Type;

.. code-block:: R

   > metawards::py_install_metawards()

This will download and install the :mod:`metawards` module into the
Python interpreter associated with reticulate. If you want to specify
the Python interpreter manually, you would need to type;

.. code-block:: R

   > reticulate::use_python("/path/to/python", required = TRUE)

before calling ``py_install_metawards()``. Here ``/path/to/python``
is the path to the Python interpreter you want to use.

You can double-check that MetaWards is available and working by typing;

.. code-block:: R

   > metawards::py_metawards_available()
   [1] TRUE

You can get the version of MetaWards Python installed using;

.. code-block:: R

   > metawards::py_version_metawards()
   [1] "1.3.0"

You can check if updates to MetaWards are available using;

.. code-block:: R

   > metawards::py_metawards_update_available()

and can update MetaWards in Python to the latest version using;

.. code-block:: R

   > metawards::py_update_metawards()

Using metawards in R
--------------------

To load the :mod:`metawards` module type;

.. code-block:: R

   > library(metawards)

This loads all of the :mod:`metawards` Python objects into the
``metawards`` namespace in R. You can then call those objects directly
as you would in Python.

For example, we can create the same custom network containing Bristol
and London in R as we did in Python via;

.. code-block:: R

   > wards <- metawards$Wards()
   > bristol <- metawards$Ward(name="Bristol")
   > bristol$add_workers(500, destination=bristol)
   > bristol$set_num_players(750)
   > print(bristol)
   Ward( id=1, name=Bristol, num_workers=500, num_players=750 )
   > london <- metawards$Ward(name="London")
   > london$add_workers(8500, destination=london)
   > london$set_num_players(10000)
   > print(london)
   Ward( id=2, name=London, num_workers=8500, num_players=10000 )
   > bristol$add_workers(500, destination=london)
   > london$add_workers(100, destination=bristol)
   > wards$add(bristol)
   > wards$add(london)
   > print(wards)
   [ Ward( id=1, name=Bristol, num_workers=1000, num_players=750 ), Ward( id=2, name=London, num_workers=8600, num_players=10000 ) ]
   > wards$to_json("custom_network.json", indent=2)
   [1] "/path/to/custom_network.json.bz2"

This should result in a (compressed) file called ``custom_network.json.bz2``,
which should have identical contents as if you have run these commands
in Python, e.g.

::

  [
    {
      "id": 1,
      "info": {
        "name": "Bristol"
      },
      "num_workers": 1000,
      "num_players": 750,
      "workers": {
        "destination": [
          1,
          2
        ],
        "population": [
          500,
          500
        ]
      }
    },
    {
      "id": 2,
      "info": {
        "name": "London"
      },
      "num_workers": 8600,
      "num_players": 10000,
      "workers": {
        "destination": [
          1,
          2
        ],
        "population": [
          100,
          8500
        ]
      }
    }
  ]

Going further
-------------

This was a simple example. However, I hope this is enough to show you how
you can begin to use the Python :mod:`metawards` module within R using
`reticulate <https://rstudio.github.io/reticulate/>`__. More details about
reticulate, more advanced ways of calling Python, plus how to set up
code completion and inline help are also available at the
`reticulate project webpage <https://rstudio.github.io/reticulate/>`__.
=============================================
Different networks for different demographics
=============================================

Custom networks significantly increase the flexibility of ``metawards``.
You can use them to create models for different countries. And you can
use them to use the concept of metapopulations to model different
environments.

For example, we can use a custom network to model ``home``, ``school``
and ``work``, and then use different demographics to model
students, teachers and everyone else. This way
you could explore how closing schools could impact disease spread.

This concept is explored in detail in the
:doc:`quick start guide <../../quickstart/index>`. The key takeaway you should
have now, having worked through this
:doc:`tutorial <../index>` as well as the
:doc:`quick start <../../quickstart/index>`, is that you can create;

* different :class:`~metawards.Disease` objects that contain different
  disease stages and parameters,
* different :class:`~metawards.Demographic` objects that represent
  different groups of individuals, each of which can experience
  different :class:`~metawards.Disease` objects (and thus advance
  along different pathways,
* who can each experience metapopulation movements on different
  :class:`~metawards.Network` networks (built flexibly via the
  :class:`~metawards.Ward` / :class:`~metawards.Wards` objects),
* where each ward in each network can have its own parameters, thereby
  enabling modelling of ward-local behaviour for different demographics at
  different disease stages in different networks.

This flexibility for the model is then enhanced by building custom
:mod:`~metawards.iterators` that enable you to write your own code
to advance individuals from one disease stage to the next, and control
how susceptible individuals are infected.

You can build custom :mod:`~metawards.mixers` to control how the force
of infections calculated for different demographics are mixed and merged
together, thereby controlling how they interact (from not interacting,
to evenly interacting, via custom interactions for the interaction matrix).

You can build custom :mod:`~metawards.movers` to control vertical
movements of individuals between demographics, that complement the
horizontal movement of individuals along disease stages.

And you can write custom :mod:`~metawards.extractors` to analyse the
data from the simulation on the fly, and write whatever data you wish
out to disk or database.

And, from here you can creata custom user parameters that control all
of this model, and write scan/design files that can run multiple
simulations to scan through these custom parameters, as well as
(nearly) all in-built and disease parameters, optionally repeating
runs to check for statistical variance.

And (finally!) you can do all of this from within a single Python or R
script using the Python or R APIs, or you can run from the command line,
or from within a batch script that will automatically parallelise
the runs over a supercomputing cluster.
================================
Specifying ward-local parameters
================================

In :doc:`part 3 <../part03/08_local_lockdown>` you learned how to read
and write ward-local parameters, e.g.
:meth:`network.nodes.scale_uv <metawards.Nodes.scale_uv>` and
:meth:`network.nodes.cutoff <metawards.Nodes.cutoff>`.

You can also set these parameters (and custom parameters) using
the :class:`~metawards.Ward` object.

For example, in Python;

.. code-block:: python

   >>> from metawards import Ward
   >>> bristol = Ward(name="bristol")
   >>> bristol.set_position(lat=51.4545, long=2.5879)
   >>> bristol.set_cutoff(15.0)
   >>> bristol.set_scale_uv(0.5)
   >>> bristol.set_custom("in_lockdown", 0.0)
   >>> bristol.set_custom("case_free_days", 21)
   >>> print(bristol.to_json(indent=2))
   {
     "id": null,
     "position": {
       "lat": 51.4545,
       "long": 2.5879
     },
     "info": {
       "name": "Bristol"
     },
     "num_workers": 0,
     "num_players": 0,
     "scale_uv": 0.5,
     "cutoff": 15.0,
     "custom": {
       "in_lockdown": 0.0,
       "case_free_days": 21.0
     }
   }

or in R

.. code-block:: R

   > library(metawards)
   > bristol <- metawards$Ward(name="bristol")
   > bristol$set_position(lat=51.4545, long=2.5879)
   > bristol$set_cutoff(15.0)
   > bristol$set_scale_uv(0.5)
   > bristol$set_custom("in_lockdown", 0.0)
   > bristol$set_custom("case_free_days", 21)
   > print(bristol$to_json(indent=2))
   {
     "id": null,
     "position": {
       "lat": 51.4545,
       "long": 2.5879
     },
     "info": {
       "name": "Bristol"
     },
     "num_workers": 0,
     "num_players": 0,
     "scale_uv": 0.5,
     "cutoff": 15.0,
     "custom": {
       "in_lockdown": 0.0,
       "case_free_days": 21.0
     }
   }

would create a ward called ``bristol``. The position of the centre of the
ward is set to a specified latitude and longitude. The ``cutoff`` distance
is set to 15 kilometers, and the ``scale_uv`` parameter is set to 0.5.
Two custom parameters are added; ``in_lockdown`` which is set to 0 and
``case_free_days`` which is set to 21.

.. note::

   You can set the position of a ward using either X/Y coordinates
   (which should be in kilometers) or latitude / longitude. All wards
   in a single network should use the same coordinates scheme.

.. note::

   You can add as many custom parameters as you like. They will all be
   stored as (double precision) floating point numbers. This means
   that ``True`` will be converted to ``1.0`` and ``42`` will be
   converted to ``42.0``.

Networks with ward-local parameters
-----------------------------------

You can combine individual ward objects together into a single Wards
object, even if they have different ward-local parameters. For example,
continue the above scripts to add a ward called ``london``, e.g. in
Python;

.. code-block:: python

   >>> london = Ward(name="london")
   >>> london.set_position(lat=51.5074, long= 0.1278)
   >>> london.set_custom("in_lockdown", 1)
   >>> bristol.add_workers(50, destination=london)
   >>> wards = bristol + london

You can now convert this into a :class:`metawards.Network`. We can then
read the parameters using the :class:`network.nodes <metawards.Nodes>`
object, as we did in :doc:`part 3 <../part03/08_local_lockdown>`.

.. code-block:: python

   >>> from metawards import Network
   >>> network = Network.from_wards(wards)
   Calculating distances...
   Total links distance equals 273.9213284848716
   Total play distance equals 0.0
   Total distance equals 273.9213284848716
   Network loaded. Population: 50, Workers: 50, Players: 0
   >>> print(network.nodes.x)
   array('d', [0.0, 51.4545, 51.5074])
   >>> print(network.nodes.y)
   array('d', [0.0, 2.5879, 0.1278])
   >>> print(network.nodes.scale_uv)
   array('d', [1.0, 0.5, 1.0])
   >>> print(network.nodes.cutoff)
   array('d', [99999.99, 15.0, 99999.99])
   >>> print(network.nodes.get_custom("in_lockdown"))
   array('d', [0.0, 0.0, 1.0])
   >>> print(network.nodes.get_custom("case_free_days"))
   array('d', [0.0, 21.0, 0.0])

The data has been correctly set (remembering that there is no ward at
index 0). Note that if ``scale_uv`` is not set, then it defaults to ``1.0``.
Similarly, if ``cutoff`` is not set then it defaults to a large distance
that is greater than the distance between two points on Earth (``99999.99``).
Custom parameters are defaulted to ``0.0`` if they are not set, e.g.
note how ``case_free_days`` for ``london`` is ``0.0``.

The same code in R would read;

.. code-block:: R

   > london <- metawards$Ward(name="london")
   > london$set_position(lat=51.5074, long= 0.1278)
   > london$set_custom("in_lockdown", 1)
   > bristol$add_workers(50, destination=london)
   > wards <- metawards$Wards()
   > wards$add(bristol)
   > wards$add(london)
   > network <- metawards$Network$from_wards(wards)
   Calculating distances...
   Total links distance equals 273.9213284848716
   Total play distance equals 0.0
   Total distance equals 273.9213284848716
   Network loaded. Population: 50, Workers: 50, Players: 0
   > print(network$nodes$x)
   array('d', [0.0, 51.4545, 51.5074])
   > print(network$nodes$y)
   array('d', [0.0, 2.5879, 0.1278])
   > print(network$nodes$scale_uv)
   array('d', [1.0, 0.5, 1.0])
   > print(network$nodes$cutoff)
   array('d', [99999.99, 15.0, 99999.99])
   > print(network$nodes$get_custom("in_lockdown"))
   array('d', [0.0, 0.0, 1.0])
   > print(network$nodes$get_custom("case_free_days"))
   array('d', [0.0, 21.0, 0.0])
===============================
Understanding the Model Network
===============================

At its core, ``metawards`` implements a meta-population model. Individuals
are grouped into workers and players, and distributed across wards. The
force of infection of a ward caused by infection of an individual, based on
their movements and behaviour is calculated, and used to decide whether
other individuals who reside or visit that ward should be infected.

The :class:`~metawards.Network` is the collection of
:class:`~metawards.Nodes` and :class:`~metawards.Links` that describe
the individual wards, and the movements between wards of their
residents.

Up to now, you have used a single underlying :class:`~metawards.Network`
for every demographic and model run in this tutorial. This does not need
to be the case, and you can choose to assign different
:class:`~metawards.Network` objects to different demographics in
a simulation. You would do this, for example, to model different connections
or movements between wards for different demographics, e.g. workers
or school students. Or to use different networks to represent movements
related to holidays or foreign travel.

The single-node network
-----------------------

You can set the overall network used by default by all demographics via
the ``--model`` or ``-m`` command-line argument. The default value
is ``2011Data``, which is based on 2011 census data, and models every
electoral ward in England and Wales.

There are two other models supplied;

* ``2011UK`` : Again, this is based on the 2011 census, but includes Scotland
  and Northern Ireland. This model is still a work in progress.
* ``single`` : This is a single-ward network that is used for testing, or
  when you don't want to model geographic behaviour.

You can run a single-ward simulation using the command line;

.. code-block:: bash

   metawards -d lurgy_home -m single -P 100

.. note::

   We are using the ``lurgy_home.json`` file from the last chapter of this
   tutorial. Note that we also need to set the population of this single
   ward using the ``-P`` or ``--population`` command line argument. In this
   case, we are setting the population to 100 individuals.

When you run, you should see that ``metawards`` is siginficantly faster.
Also, as the outbreak is not seeded, nothing much happens.

Seeding the outbreak
--------------------

Up to now, you have been using the ``ExtraSeedsLondon.dat`` file to seed
every outbreak. This file contains the single line;

::

    1       5   255

The three numbers instruct ``metawards`` to seed the outbreak on day 1
(the first number), infecting 5 individuals (the second number) in
ward 255 (the third number).

.. warning::

   The name of this file is misleading, as it is really seeding
   the ward with index 255, not seeding a ward in London. This will
   only be seeding London if the network this is used with has a
   ward in London at index 255.

If you try to seed using this file you will get an error, e.g.

.. code-block:: bash

   metawards -d lurgy_home -m single -P 100 -a ExtraSeedsLondon.dat

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ ERROR ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Unable to seed the infection using (1, 255, 5, None). The error was <class 'IndexError'>: Invalid node index
    255. Number of nodes in this container equals 2. Please double-check that you are trying to seed a node that
    exists in this network.

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Traceback (most recent call last):

    [lots of output]

    IndexError: Invalid node index 255. Number of nodes in this container equals 2

.. note::

   The traceback can be long and complex, and is really only of use for ``metawards``
   developers. You can normally work out what has gone wrong by scrolling up
   to before the traceback, and seeing if there is a ``ERROR`` printed immediately
   before.

We have provided an ``ExtraSeedsOne.dat`` file, which seeds 5 infections
in ward 1 on the first day of the outbreak. This contains the line;

::

  1       5        1

which says to seed 5 individuals on the first day in ward 1.

.. note::

   Nodes are indexed from 1 rather than 0. This means that ``ward[1]``
   is the first node, and ``ward[nnodes]`` is the last node in the
   network.

You can use this with single-ward networks, e.g.

.. code-block:: bash

   metawards -d lurgy_home -m single -P 100 -a ExtraSeedsOne.dat

However, creating a file to seed an outbreak is inconvenient, particularly
when you only want to seed a single ward. You can, optionally, pass
the seeding information as the argument, instead of the filename. Thus
this will work;

.. code-block:: bash

   metawards -d lurgy_home -m single -P 100 -a "1 5 1"

Equally, you can seed on multiple days, e.g. seeding 5 individuals on day 1,
and then 10 individuals on day 2, via;

.. code-block:: bash

   metawards -d lurgy_home -m single -P 100 -a "1 5 1\n2 10 1"

The contents of the string is interpreted identically to if it had been
read from a file, with ``\n`` representing a newline character.

.. note::

   The extra seeds file has a flexible and powerful format, e.g. supporting
   seeding by date, seeding random wards or by random amounts etc.
   More information on the format of this file can
   be :doc:`found here <../../fileformats/extraseeds>`.
=========================
Creating a Custom Network
=========================

While the default networks supplied in ``MetaWardsData`` model the UK, there
is nothing to stop you creating your own network to model any country or
region. Indeed, you can use the concepts of ``wards``, ``workers`` and
``players`` in a more generic way to model lots of different environments,
e.g.

* Using wards to represent different university buildings, and then track
  disease spread between buildings as staff and students move around.
* Using wards to represent care homes, hospitals and homes in a single
  region, and then model the motion of staff, patients and the general
  population between these different environments.

Network file formats
--------------------

A lot of data needs to be loaded to define the network. There are two
file formats for specifying this data;

1. :doc:`A set of fixed-format files <../../fileformats/network>` that
   contain the data in a set of files that are contained in a single
   directory. This is an older format that is used predominantly for
   older model networks.

2. A JSON-format file that contains all of the data needed to describe the
   network in a single file. This file should only be manipulated or
   edited using the Python API described below.

Creating and editing Networks in Python
---------------------------------------

The best way to create a new network is to use the Python API. A
:class:`~metawards.Network` is edited or created via the
:class:`~metawards.Wards` class. This represents an editable collection
of individual :class:`~metawards.Ward` objects, each of which
represents a ward. The :class:`~metawards.Ward` provides functions
for setting the name and metadata for a ward, plus adding work and
play connections to other wards.

For example, we can interactively create a new :class:`~metawards.Network`
using the :class:`~metawards.Ward` and :class:`~metawards.Wards` classes
in, e.g. ipython or a jupyter notebook;

First, we will import the necessary classes and create our
:class:`~metawards.Wards` object, which we will call ``wards``;

.. code-block:: python

   >>> from metawards import Ward, Wards, Network
   >>> wards = Wards()

Next, we will create a :class:`~metawards.Ward` object to represent
Bristol (which we will call ``bristol``). We will add `500`
workers who will work in Bristol, and ``750`` players.

.. code-block:: python

   >>> bristol = Ward(name="Bristol")
   >>> bristol.add_workers(500)
   >>> bristol.set_num_players(750)

Next, we will create a :class:`~metawards.Ward` object to represent
London (which we will call ``london``). We will add ``8600``
workers and ``10000`` players.

   >>> london = Ward(name="London")
   >>> london.add_workers(8500)
   >>> london.set_num_players(10000)

Now, we will add some commuters. We will have ``500`` Bristolians
commute each day to London, while ``100`` Londoners
will commute each day to Bristol.

.. code-block:: python

   >>> bristol.add_workers(500, destination=london)
   >>> london.add_workers(100, destination=bristol)

We can confirm that the information is correct by printing, e.g.

.. code-block:: python

   >>> print(bristol)
   Ward( name=Bristol, num_workers=1000, num_players=750 )

   >>> print(london)
   Ward( name=London, num_workers=8600, num_players=10000 )

Next, we add the two :class:`~metawards.Ward` objects to our
:class:`~metawards.Wards` object that represents the entire model.

.. code-block:: python

   >>> wards.add(bristol)
   >>> wards.add(london)
   >>> print(wards)
   [ Ward( info=Bristol, id=1, num_workers=1000, num_players=750 ), Ward( info=London, id=2, num_workers=8600, num_players=10000 ) ]

.. note::

   Note that each :class:`~metawards.Ward` in the :class:`~metawards.Wards`
   collection has been automatically assigned an ID number (1 for Bristol
   and 2 for London). You can refer to the wards by their ID number, but
   should, in general, leave ``metawards`` to automatically generate
   and manage these IDs.

We can now save this set of :class:`~metawards.Wards` to a file by converting
this to a data dictionary, and then serialising that dictionary to JSON,
which we will stream to the file ``custom_network.json.bz2``.

.. code-block:: python

   >>> wards.to_json("custom_network.json", indent=2)

The resulting JSON file will look something like this;

::

  [
    {
      "id": 1,
      "info": {
        "name": "Bristol"
      },
      "num_workers": 1000,
      "num_players": 750,
      "workers": {
        "destination": [
          1,
          2
        ],
        "population": [
          500,
          500
        ]
      }
    },
    {
      "id": 2,
      "info": {
        "name": "London"
      },
      "num_workers": 8600,
      "num_players": 10000,
      "workers": {
        "destination": [
          1,
          2
        ],
        "population": [
          100,
          8500
        ]
      }
    }
  ]

.. note::

   Note that the exact format of the JSON will change as ``metawards``
   evolves. We will retain backwards compatibility, meaning that newer
   versions of ``metawards`` will be able to read old files, but older
   versions may not be able to read new files.

   Note that the file will be automatically compressed using bzip2. You
   can disable this by setting ``auto_bzip=False``.

   Note also that ``indent=2`` just sets the indentation used for printing.
   You can set whatever indentation you want, including not setting any.
   It won't affect the included information - just its human-readability.

You can load this JSON file into a Wards object using;

.. code-block:: python

   >>> wards = Wards.from_json("custom_network.json.bz2")
   >>> print(wards)
   [ Ward( id=1, name=Bristol, num_workers=1000, num_players=750 ), Ward( id=2, name=London, num_workers=8600, num_players=10000 ) ]
==============================
Multi-demographic sub-networks
==============================

So far every member of the population during a ``metawards`` *model run* has
been treated equally. There was no concept of different groups of individuals,
e.g. school children, hospital patients, holiday-makers etc. being
modelled any differently. The only distinction was between *workers*,
who made fixed movements (potentially between wards)
during the day, and random interactions within their home ward in
the evening, and *players*, who just made random interactions within
their home ward.

Understanding Network
---------------------

The :class:`~metawards.Network` class models this network of *workers* and
*players*. It holds a set of :class:`~metawards.Nodes` that represent each
electoral ward, and a set of :class:`~metawards.Links` that represent the
fixed motions between wards made by the "workers".

The :class:`~metawards.Node` class holds the number of *players* who
are susceptible to infection (S). This is held in the variable
:data:`metawards.Node.play_suscept`.

The :class:`~metawards.Link` class contains the number of *workers* who are
susceptible to infection (S) who travel regularly between their
"home" ward and their "work" ward. This is held in the
variable :data:`metawards.Link.suscept`.

During a *model run* these values will be reduced as individuals are
infected and then progressed through the different disease stages.

.. note::

  These variables are stored as double precision floating point numbers,
  despite representing integer values. This is because it
  is more efficient to use double precision numbers when they are used
  for calculations during a *model run*.

Understanding Networks
----------------------

To model multiple demographics, we need to create space for multiple different
sets of *workers* and *players*. This is achieved in ``metawards`` by
giving each demographic its own :class:`~metawards.Network`. All of these
demographic sub-networks are merged and managed via the
:class:`~metawards.Networks` class.

We create a :class:`~metawards.Networks` object by first specifying the way
we would like to distribute the population between different demographics.
We do this by writing a *demographics* file, which is a simple
`JSON-format <https://guide.couchdb.org/draft/json.html>`__
file. First, create a file called ``demographics.json`` and copy in the below;

::

    {
      "demographics" : ["red", "blue"],
      "work_ratios"  : [ 0.0,   1.0  ],
      "play_ratios"  : [ 0.5,   0.5  ]
    }

This file specifies that the population will be distributed between
two demographics, *red* and *blue*, named in the ``demographics``
field. You can have as many or as few demographics as you wish, and
can name them in any way you want.

The ``work_ratios`` lists the ratio (or percentage) of the *worker* population
that should belong to each demographic. For example, here all of the
*workers* are in the *blue* demographic, and none of the workers are
in the *red* demographic.

The ``play_ratios`` lists the ratio (or percentage) of the *player* population
that should belong to each demographic. For example, here 50% of the
*players* are in the *red* demographic, and 50% of the *players* are in the
*blue* demographic.

.. note::
  You can specify the work and play ratios using either numbers between
  0.0 and 1.0, or you can pass in strings that are interpreted using
  :func:`metawards.Interpret.number`, e.g. "50%", "1/4" or
  "(10+15)%". The only requirement is that the sum of ratios must
  equal 1.0 (or 100%), as every individual must be assigned to one
  of the demographics.

Now that you have created the ``demographics.json`` file, you can tell
``metawards`` to use it via the ``--demographics``, or ``-D``,
command line argument. Run ``metawards`` using;

.. code-block:: bash

   metawards -d lurgy3 -D demographics.json

In the output you should see lines such as;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Specialising into demographics ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Demographics demographics.json
    loaded from demographics.json
    version: None
    author(s): unknown
    contact(s): unknown
    references(s): none
    repository: None
    repository_branch: None
    repository_version: None
    demographics = [
        Demographic(name='red', work_ratio=0.0, play_ratio=0.5, adjustment=None)
        Demographic(name='blue', work_ratio=1.0, play_ratio=0.5, adjustment=None)
    ]
    Seeding generator used for demographics with seed 4751828
    Using random number seed: 4751828
    Number of differences is 0.0 + 0.0
    Specialising network - population: 56082077
    red - population: 16806530
    blue - population: 39275547

These show that your demographics file was read correctly. In this case,
this has specialised the :class:`~metawards.Network` which modelled a
population of 56082077 individuals into a :class:`~metawards.Networks`
which has a *red* population of 16806530 and a *blue* population of
39275547.

.. warning::
  A fixed random number seed is used to assign left-over individuals
  to a random demographic. For example, 10 individuals cannot be divided
  equally between 3 demographics, so one randomly chosen demographic
  will have 4 individuals, while the other two will have 3. This
  division is performed by ``metawards`` in every single
  :class:`~metawards.Node` and every single :class:`~metawards.Link`,
  to ensure that every individual is allocated. This random seed is
  hard-coded to ``4751828``. Or, you can set it for a demographic
  by adding ``"random_seed" = number`` to the *demographics* file,
  e.g. ``"random_seed" = 10859403``.

Once the :class:`~metawards.Networks` had been specialised, the *model run*
was performed as before. Now, the output shows the S, E, I, R values
for both the overall total population, and also for the demographic
sub-network populations, e.g.

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
       red  S: 16806530  E: 0  I: 0  R: 0  IW: 0  POPULATION: 16806530
      blue  S: 39275547  E: 0  I: 0  R: 0  IW: 0  POPULATION: 39275547
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
       red  S: 16806530  E: 0  I: 0  R: 0  IW: 0  POPULATION: 16806530
      blue  S: 39275547  E: 0  I: 0  R: 0  IW: 0  POPULATION: 39275547
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
       red  S: 16806530  E: 0  I: 0  R: 0  IW: 0  POPULATION: 16806530
      blue  S: 39275547  E: 0  I: 0  R: 0  IW: 0  POPULATION: 39275547
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
       red  S: 16806530  E: 0  I: 0  R: 0  IW: 0  POPULATION: 16806530
      blue  S: 39275547  E: 0  I: 0  R: 0  IW: 0  POPULATION: 39275547
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
       red  S: 16806530  E: 0  I: 0  R: 0  IW: 0  POPULATION: 16806530
      blue  S: 39275547  E: 0  I: 0  R: 0  IW: 0  POPULATION: 39275547
    Number of infections: 0
    Infection died ... Ending on day 5

In this case no infection was seeded, so nothing appears to happen.

We can seed an infection just as before, by using the ``--additional``
(or ``-a``) option, e.g. now run;

.. code-block:: bash

   metawards -d lurgy3 -D demographics.json -a ExtraSeedsLondon.dat

You should see output similar (but not identical) to;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 249 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56079607  E: 0  I: 2  R: 2468  IW: 0  POPULATION: 56082077
       red  S: 16804060  E: 0  I: 2  R: 2468  IW: 0  POPULATION: 16806530
      blue  S: 39275547  E: 0  I: 0  R:    0  IW: 0  POPULATION: 39275547
    Number of infections: 2

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 250 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56079607  E: 0  I: 2  R: 2468  IW: 0  POPULATION: 56082077
       red  S: 16804060  E: 0  I: 2  R: 2468  IW: 0  POPULATION: 16806530
      blue  S: 39275547  E: 0  I: 0  R:    0  IW: 0  POPULATION: 39275547
    Number of infections: 2

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 251 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56079607  E: 0  I: 1  R: 2469  IW: 0  POPULATION: 56082077
       red  S: 16804060  E: 0  I: 1  R: 2469  IW: 0  POPULATION: 16806530
      blue  S: 39275547  E: 0  I: 0  R:    0  IW: 0  POPULATION: 39275547
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 252 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56079607  E: 0  I: 1  R: 2469  IW: 0  POPULATION: 56082077
       red  S: 16804060  E: 0  I: 1  R: 2469  IW: 0  POPULATION: 16806530
      blue  S: 39275547  E: 0  I: 0  R:    0  IW: 0  POPULATION: 39275547
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 253 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56079607  E: 0  I: 0  R: 2470  IW: 0  POPULATION: 56082077
       red  S: 16804060  E: 0  I: 0  R: 2470  IW: 0  POPULATION: 16806530
      blue  S: 39275547  E: 0  I: 0  R:    0  IW: 0  POPULATION: 39275547
    Number of infections: 0
    Infection died ... Ending on day 254

By default, infections are seeded into the first demographic (in this case
*red*). This demographic are *players*, so only interact in their home
ward via random interactions. As such, the infection did not spread
beyond that home ward and so it died out quite quickly.

Seeding different demographics
------------------------------

You can seed different demographics by specifying the demographic in
the additional seeding file. Create a new seeding file called
``ExtraSeedsLondonBlue.dat`` and copy in the below;

::

  1  5  255  blue

The format of this file is a list of lines that say which wards should
be seeded. In this case, there is just one line containing four values.

* The first value (``1``) is the day of seeding, in this case on day 1.
* The second value  (``5``) is the number of individuals to infect, in
  this case 5.
* The third value (``255``) is the index of the ward to infect. You can find
  the index of the ward you want using the :class:`~metawards.WardInfos`
  object, e.g. via ``network.info.find("...")``.
* The fourth value (``blue``) is the name or index of the demographic
  you want to seed.

.. note::

  If you want, you could have specified the demographic in this file by
  its index (``1``) rather than by its name (``blue``). It is up to you.

The *blue* demographic contains all of the *workers*, so we would expect
to see a different outbreak. Perform a *model run* using;

.. code-block:: bash

  metawards -d lurgy3 -D demographics.json -a ExtraSeedsLondonBlue.dat

You should see a more sustained outbreak, ending in a similar way to this;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 147 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 16864032  E: 0  I: 1  R: 39218044  IW: 0  POPULATION: 56082077
       red  S: 16806530  E: 0  I: 0  R:        0  IW: 0  POPULATION: 16806530
      blue  S:    57502  E: 0  I: 1  R: 39218044  IW: 0  POPULATION: 39275547
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 148 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 16864032  E: 0  I: 1  R: 39218044  IW: 0  POPULATION: 56082077
       red  S: 16806530  E: 0  I: 0  R:        0  IW: 0  POPULATION: 16806530
      blue  S:    57502  E: 0  I: 1  R: 39218044  IW: 0  POPULATION: 39275547
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 149 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 16864032  E: 0  I: 0  R: 39218045  IW: 0  POPULATION: 56082077
       red  S: 16806530  E: 0  I: 0  R:        0  IW: 0  POPULATION: 16806530
      blue  S:    57502  E: 0  I: 0  R: 39218045  IW: 0  POPULATION: 39275547
    Number of infections: 0
    Infection died ... Ending on day 150

Because the *blue workers* could move between wards, they were able to carry
the infection across the country, meaning that most members of the *blue*
demographic were infected.
==========================
Mixers and Merge-functions
==========================

In the last page you modelled outbreaks that were seeded either in
the *red* population, who were all *players*, and the *blue* population,
which contained both *players* and all of the *workers*.

This resulted in two different outbreaks, which was understandable given
that only *workers* move between wards. One issue was that the demographics
were completely separate - infected individuals in the *red* population
couldn't infect individuals in the *blue* population and vice-versa.

This is because the :class:`~metawards.Network` for each demographic was
advanced independently, and no mixing took place.

Force of infection (FOI)
------------------------

Infected individuals in a ward contribute to the *force of infection* (FOI)
of that ward. The higher the FOI, the more likely it is that other
individuals in a ward will become infected.

The FOIs for each demographic sub-network are calculated independently,
based only on the infected individuals from that demographic. We can thus
use the FOI calculation as a way of enabling different demographics to
mix. We do this in ``metawards`` using custom functions, called
:doc:`../../api/index_MetaWards_mixers`. These choose different
``merge functions`` that are
used to mix and merge calculations of FOIs across different demographics.

The default mixer is :func:`~metawards.mixers.mix_none`, which, as the
name suggests, performs no mixing between demographics.

Mix evenly
----------

An alternative is to use :func:`~metawards.mixers.mix_evenly`. This evenly
mixes all demographics. It sets the FOI in each ward in each demographic
equal to the sum of the FOIs in that ward across all demographics.

You can choose this mixer using the ``--mixer`` argument. Run ``metawards``
using;

.. code-block:: bash

  metawards -d lurgy3 -D demographics.json -a ExtraSeedsLondon.dat --mixer mix_evenly

You should now see that the outbreak spreads from the initial infection in
the *red* demographic to the *blue workers*, who then spread it around
the country to both the *red* and *blue* groups.

The trajectories for each of the demographics are written into a single
csv file called ``output/trajectory.csv.bz2``. This can be loaded
into tools such as Excel, R or Pandas for detailed analysis. For example,
in Pandas you can load the file using;

.. code-block:: python

   >>> import pandas as pd
   >>> df = pd.read_csv("output/trajectory.csv.bz2")
   >>> print(df)
         day        date demographic         S  E  I         R  IW
    0      0  2020-05-07     overall  56082077  0  0         0   0
    1      0  2020-05-07         red  16806501  0  0         0   0
    2      0  2020-05-07        blue  39275576  0  0         0   0
    3      1  2020-05-08     overall  56082077  0  0         0   0
    4      0  2020-05-08         red  16806501  0  0         0   0
    ..   ...         ...         ...       ... .. ..       ...  ..
    388    0  2020-09-13         red   1609032  0  0  15197469   0
    389    0  2020-09-13        blue   1101473  0  1  38174102   0
    390  130  2020-09-14     overall   2710505  0  0  53371572   0
    391    0  2020-09-14         red   1609032  0  0  15197469   0
    392    0  2020-09-14        blue   1101473  0  0  38174103   0

    [393 rows x 8 columns]

You can see that the values for each day for each demographic are printed
in turn. The column called ``demographic`` holds the name of the demographic,
with the overall total network called ``overall``.

You can plot this data easily using ``metawards-plot``, e.g.

.. code-block:: bash

   metawards-plot -i output/trajectory.csv.bz2

This will produce a graph called ``output/demographics.jpg``.
This should look similar to this;

.. image:: ../../images/tutorial_5_2_demographics.jpg
   :alt: Disease trajectory across the red and blue demographics

Notice how infections in the *red* demographic lag behind those of the
*blue* demographic. This makes sense, as it is only the workers from the
*blue* demographic that can be infected outside their home ward. The
*blue workers* are thus the source of the infection for the *red* demographic.
======================================
Contact probabilities and Demographics
======================================

MetaWards implements the mixing of demographics by merging together
the *force of infections* (FOIs) calculated for each network. By default
this mixing takes no account of the number of individuals (N) in each
ward in each demographic.

The FOI, :math:`F(w)` calculated in each ward, :math:`w` is calculated
separately for each demographic via the equations
:doc:`described here <../part03/07_cutoff>`.

The merged FOI for ward :math:`w` for demographic :math:`i`,
:math:`F_i'(w)`, is calculated as
the sum of the FOIs calculated for ward :math:`w` across
all demographics, :math:`i, j, k...`
scaling by row :math:`i` of the interaction matrix :math:`M` via;

:math:`F_i'(w) = M_{ii} F_i(w) + M_{ij} F_j(w) + M_{ik} F_k(w) + ...`

This is then divided by the number of individuals in ward :math:`w` in
demographic :math:`i`, :math:`N_i(w)` before it is converted into a probability
of infection.

This normalisation by :math:`N_i(w)` means that it is up to you to
account for the contact probability in the interation matrix. The individual
terms should account for both the different infectivity of the different
demographics, as well as the probability that individuals in one
demographic will come into contact with individuals of other demographics.

However, accounting for these different contact rates using just
the interaction matrix is difficult (or impossible) for cases where
you have multiple wards that have different numbers of individuals. This
is because different normalisation factors would be needed for the
FOIs calculated for different wards in different demographics.

Modelling demographics as part of the same population
-----------------------------------------------------

If all demographics are part of the same population then we assume
that any member in a ward of any demographic is equally likely to contact
any other member of the same ward in any other demographic.

In this case, the normalisation factor should be the total number of
individuals in that ward summed over all demographics,
:math:`N(w)`. This is calculated via;

:math:`N(w) = N_i(w) + N_j(w) + N_k(w) + ...`

The merge function :func:`~metawards.mixers.merge_matrix_single_population`
calculates the merged FOI via the equation;

:math:`F_i'(w) = M_{ii} F_i(w) \frac{N_i(w)}{N(w)} + M_{ij} \frac{N_i(w)}{N(w)} F_j(w) + M_{ik} \frac{N_i(w)}{N(w)} F_k(w) + ...`

Because this FOI is divided by :math:`N_i(w)` before it is converted into
an infection probability, this is equivalent to using :math:`N(w)` as the
normalisation factor for ward :math:`w` in all demographics.

This therefore models the demographics as a single population where individuals
in a ward across all demographics have an equal probability of contact with
each other. The interaction matrix has the effect of scaling :math:`\beta`,
i.e. the probability of infection at each contact between members of
different demographics. You can pass in a custom interaction matrix
if you want control of this, or can use the convenience mixer functions
:func:`~metawards.mixers.mix_evenly_single_population`,
and :func:`~metawards.mixers.mix_none_single_population` if you want
to use an interactions matrix of all ones (:math:`M = |1|`) that
represents even infections, or a identity interaction matrix with
a diagonal of 1 and off-diagonal of 0, that represents no infections between
members of different demographics, respectively.

Modelling demographics as multiple populations
----------------------------------------------

It may be that the demographics represent different populations that
have a differnet probability of meeting each other. This would be
the case for the classic use of demographics, e.g. using demographics
to represent different age groups, and having a different contact
probability between individuals of different age groups.

In this case, the FOI calculated from each demographic must be normalised
by the number of individuals in that demographic before it is combined
into the merged demographic. To achieve this,
the merge function :func:`~metawards.mixers.merge_matrix_multi_population`
calculates the merged FOI via the equation;

:math:`F_i'(w) = M_{ii} F_i(w) \frac{N_i(w)}{N_i(w)} + M_{ij} \frac{N_i(w)}{N_j(w)} F_j(w) + M_{ik} \frac{N_i(w)}{N_k(w)} F_k(w) + ...`

Because this FOI is divided by :math:`N_i(w)` before it is converted into
an infection probability, this is equivalent to using :math:`N_x(w)` as the
normalisation factor for ward :math:`w` for demographic :math:`x`.

This therefore models the demographics as being part of separate populations
where the probability of contact between individuals in the same ward
in different demographics depends on the number of individuals in each
demographic.

The interaction matrix has the effect of scaling :math:`\beta`,
i.e. the probability of infection at each contact between members of
different demographics. You can pass in a custom interaction matrix
if you want control of this, or can use the convenience mixer functions
:func:`~metawards.mixers.mix_evenly_multi_population`,
and :func:`~metawards.mixers.mix_none_multi_population` if you want
to use an interactions matrix of all ones (:math:`M = |1|`) that
represents even infections, or a identity interaction matrix with
a diagonal of 1 and off-diagonal of 0, that represents no infections between
members of different demographics, respectively.

Modelling more complex population dynamics
------------------------------------------

For more complex population dynamics you will need to write your own
custom mixer function that accounts for the number of individuals
in each demographic in each ward, :math:`N_x(w)`. Because this will
be potentially an expensive function, you should write it as a
cython plugin using the `.pyx` extension. MetaWards will automatically
compile and load this plugin at runtime.

You should take inspiration from :func:`~metawards.mixers.merge_matrix_single_population`
and :func:`~metawards.mixers.merge_matrix_multi_population`. These show
how the FOIs calculated for each ward are combined with the interaction
matrix and the number of individuals in each ward in each demographic.
Note that the code is slightly more complex than described here, as the
FOIs and number of individuals are calculated separately both for
daytime and nighttime, and the number of individuals is the sum
of the number of players (who make random movements during the day)
and the number of workers (who make predictable movements during the day).

If you have any questions or would like help, please feel free to get
in touch by `raising an issue on the GitHub repository <https://github.com/metawards/MetaWards/issues>`_.
===============================================
Using custom merge functions to model shielding
===============================================

So far you have modelled multiple demographics, but they either
interacted fully with each other during an outbreak, or have
not interacted at all.

One of the benefits of using demographics is that we can model
different levels of interaction between different groups. To do
this, we will need to set up an interaction matrix between
demographics, and to call this via a custom mixer.

Interaction matrix
------------------

An interaction matrix specifies how much the *force of infection* (FOI)
calculated from one demographic influences another. It is an asymmetric
matrix, e.g. for our *red* and *blue* demographics it could be;

.. list-table::
   :widths: 30 35 35
   :header-rows: 1
   :stub-columns: 1

   * -
     - red
     - blue
   * - red
     - 0.2
     - 0.1
   * - blue
     - 0.8
     - 1.0

The diagonals of this matrix show how much infections in each demographic
affect that demographic, i.e. how easily the infection spreads within the
demographic. It has a similar affect to ``scale_uv``, though is not
equivalent.

.. note::

  ``scale_uv`` scales the precursor to the calculation of each ward's day
  and night FOIs, while the *interaction matrix* mixes the day and night
  FOIs for each ward together.

The top-left value (0.2) says that *red* demographic is adopting lockdown
measures that scale down the FOI between *red* individuals by 0.2.

In contrast, the bottom-right value (1.0) says that no measures are being
taken within the *blue* demographic, and the FOI calculated from infections
within that demographic are not scaled.

The top-right value (0.1) gives the amount by which the *red* demographic is
affected by infections within the *blue* demographic.
This describes how easily infections from the *blue*
demographic can be transmitted to the *red* demographic. This low value of
0.1 indicates that care is taken by members of the *blue* demographic to
not infect members of the *red* demographic.

The bottom-left value (0.8) is the reverse, namely how much the *blue*
demographic is affected by infections within the *red* demographic,
i.e. how easily infections
from the *red* demographic can be transmitted to the *blue* demographic.
The larger value of 0.8 indicates that the *blue* demographic are still
being infected by members of the *red* demographic,
e.g. perhaps because they are treating or caring for them and are thus
exposed to large numbers of members.

The values in this matrix correspond to a potential shielding scenario,
whereby measures are taken by the *blue* demographic to shield the
*red* demographic from infection.

Custom mixers
-------------

We can use this interaction matrix by creating a custom mixer.
Create a file called ``shield.py`` and add the following;

.. code-block:: python

    from metawards.mixers import merge_using_matrix

    def mix_shield(network, **kwargs):
        matrix = [ [0.2, 0.1],
                   [0.8, 1.0] ]

        network.demographics.interaction_matrix = matrix

        return [merge_using_matrix]

Here we have created a new mixer called ``mix_shield``. This function
has two purposes:

1. It sets :data:`~metawards.Demographics.interaction_matrix` equal to
   the desired interaction matrix in the networks'
   :class:`~metawards.Demographics` object.
   The matrix is as described above, and represents shielding of
   the *red* demographic by the *blue* demographic.

2. It returns the :meth:`~metawards.mixers.merge_using_matrix` merge
   function as the function to use to merge the FOIs.
   :meth:`~metawards.mixers.merge_using_matrix`
   reads the interaction matrix from
   ``network.demographics.interaction_matrix``, meaning that the
   value set in step 1 will be used for the merge.

You set the mixer to use using the ``--mixer`` flag, e.g. run ``metawards``
using;

.. code-block:: bash

   metawards -d lurgy3 -D demographics.json -a ExtraSeedsLondon.dat --mixer shield

You should see that, while the infection moves through most of the *blue*
demographic, it is relatively contained within the *red* demographic.

You can plot the trajectory using;

.. code-block:: bash

   metawards-plot -i output/trajectory.csv.bz2

You should see a plot similar to this;

.. image:: ../../images/tutorial_5_3_1_demographics.jpg
   :alt: Disease trajectory for a shielding scenario for the red demographic

Adjusting shielding parameters
------------------------------

This has worked well, in that the shielded *red* demographic has been
protected from the disease. However, using scaling factors of 0.2 and
0.1 is quite extreme, especially over the four months of the model
outbreak.

We can use adjustable parameters to investigate how much shielding is
needed to protect the *red* demographic. To do this, update your
``shield.py`` file to contain;

.. code-block:: python

    from metawards.mixers import merge_using_matrix

    def mix_shield(network, **kwargs):
        params = network.params

        red_red = params.user_params["red_red"]
        red_blue = params.user_params["red_blue"]
        blue_red = params.user_params["blue_red"]
        blue_blue = params.user_params["blue_blue"]

        matrix = [ [red_red , red_blue ],
                   [blue_red, blue_blue] ]

        network.demographics.interaction_matrix = matrix

        return [merge_using_matrix]

Here we have adapted our ``mix_shield`` function to get the values for
the interaction matrix from adjustable user parameters that have been
set using :doc:`same mechanism as before <../part03/05_scanning>`.
In this case we have called the parameters ``red_red``, for the impact
of *red* on *red*, ``red_blue`` for the impact of *blue* on *red* etc.

We then need to create an input file to set the initial values of these
parameters. Create such a file called "shield.inp" and copy in;

::

    .red_red   = 0.2
    .red_blue  = 0.1
    .blue_red  = 0.8
    .blue_blue = 1.0

Finally, we would like to scan through the different value of
``red_red`` and ``red_blue`` to see how much the *red* demographic
needs to be shielded. Create a scan file called ``scan.dat`` and copy in;

::

  .red_red  .red_blue
     0.2       0.1
     0.2       0.2
     0.2       0.3
     0.2       0.4
     0.2       0.5

     0.3       0.1
     0.3       0.2
     0.3       0.3
     0.3       0.4
     0.3       0.5

     0.4       0.1
     0.4       0.2
     0.4       0.3
     0.4       0.4
     0.4       0.5

     0.5       0.1
     0.5       0.2
     0.5       0.3
     0.5       0.4
     0.5       0.5

This scans ``red_red`` between 0.2 and 0.5 while scanning ``red_blue``
from 0.1 to 0.5

You can run these jobs using this command;

.. code-block:: bash

   metawards -d lurgy3 -D demographics.json -a ExtraSeedsLondonBlue.dat --mixer shield --user-variables shield.inp -i scan.dat

or, alternatively if you have a cluster you could use a job script such
as this to run multiple repeats (always a good idea for a stochastic
simulation).

.. code-block:: bash

    #!/bin/bash
    #PBS -l walltime=12:00:00
    #PBS -l select=4:ncpus=64:mem=64GB
    # The above sets 4 nodes with 64 cores each

    # Assume you have metawards in $HOME/envs/metawards
    source $HOME/metawards/bin/activate

    # change into the directory from which this job was submitted
    cd $PBS_O_WORKDIR

    metawards -u shield.inp -i scan.dat -d lurgy3 \
              -D demographics.json -a ExtraSeedsLondonBlue.dat \
              --mixer shield \
              --repeats 8 --nthreads 16 --force-overwrite-output

if you are using PBS, or

::

    #!/bin/bash
    #SBATCH --time=01:00:00
    #SBATCH --ntasks=4
    #SBATCH --cpus-per-task=64
    # The above sets 4 nodes with 64 cores each

    # Assume you have metawards in $HOME/envs/metawards
    source $HOME/metawards/bin/activate

    metawards -u shield.inp -i scan.dat -d lurgy3 \
              -D demographics.json -a ExtraSeedsLondonBlue.dat \
              --mixer shield \
              --repeats 8 --nthreads 16 --force-overwrite-output

This job may take a while (likely 1-2 minutes per *model run*, and then
scaled by number of jobs divided by number of cores). In my case,
this took about 16 minutes on 256 cores of
`Catalyst <https://www.bristol.ac.uk/news/2018/april/supercomputer-collaboration.html>`__.

Once you have performed this calculation you can generate an animation
of the overview graphs using;

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2
   metawards-plot --animate output/overview*.jpg -o output/overview.gif

(assuming all of your output is in the ``output`` directory)

Your animation should look something like this;

.. image:: ../../images/tutorial_5_3_2.gif
   :alt: Overview image of shielding with custom parameters

You can also generate the individual demographic trajectory overview plots
and animate those using;

.. code-block:: bash

   metawards-plot -i output/*x001/trajectory.csv.bz2
   metawards-plot --animate output/*x001/demographics.jpg -o output/demographics.gif

.. note::
   We've only generated and animated the first repeat here (directories
   are all named "SOMETHINGx001"). This makes processing quicker and the
   resulting animation smaller, as each repeat has almost the same plot.

Your animation should look something like this;

.. image:: ../../images/tutorial_5_3_3.gif
   :alt: Demographic trajectories when shielding with custom parameters

From this scan it is clear that the ``red-blue`` scale has a much bigger
impact on the success of shielding than ``red-red``.
This suggests, at least in this model,
that it is more important for the *blue* demographic to take care when
interacting with the *red* demographic than it is to control the level
of lockdown of the *red* demographic.
==========================
Refining the disease model
==========================

In the last part of the tutorial you plotted graphs that showed the progress
of the lurgy across four *model runs* of the outbreak which were seeded
in London.

While this was good, the *model runs* were performed using an initial,
possibly poor model of the lurgy. Our next step is to customise the
``disease file`` to build a more representative model.

Understanding the disease file
------------------------------

The model of the lurgy is based on the
`disease file <https://github.com/metawards/MetaWardsData/blob/main/diseases/lurgy.json>`__
contained in the
`MetaWardsData <https://github.com/metawards/MetaWardsData>`__ repository.

The first step to customise this file is to copy it into our current
directory and rename it to ``lurgy2.json``.
Assuming you have set the ``METAWARDSDATA`` environment
variable equal to the path to this directory, type;

.. code-block:: bash

  cp $METAWARDSDATA/diseases/lurgy.json ./lurgy2.json

on Linux or Mac OS X, or type;

.. code-block:: bash

  copy $Env:METAWARDSDATA\diseases\lurgy.json lurgy2.json

on Windows.

.. note::
   The above command work on Linux, Mac (OS X) or Windows. If you are following
   this tutorial using a different operating system and know the correct
   copy command to use, then please send it to us by sending us a
   `pull request <https://github.com/metawards/MetaWards>`__
   or `posting an issue <https://github.com/metawards/MetaWards/issues>`__.

The disease file for the lurgy looks like the following;

::

  { "name"             : "The Lurgy",
    "version"          : "April 16th 2020",
    "author(s)"        : "Christopher Woods",
    "contact(s)"       : "christopher.woods@bristol.ac.uk",
    "reference(s)"     : "Completely ficticious disease - no references",
    "beta"             : [0.0, 0.0, 0.5, 0.5, 0.0],
    "progress"         : [1.0, 1.0, 0.5, 0.5, 0.0],
    "too_ill_to_move"  : [0.0, 0.0, 0.5, 0.8, 1.0],
    "contrib_foi"      : [1.0, 1.0, 1.0, 1.0, 0.0]
  }

The file is in `JSON <https://en.wikipedia.org/wiki/JSON>`__ format. This
is easy for computers to read, and relatively easy for us to read too ;-)

Notice that the file contains metadata, e.g. **version**, **author(s)**,
**references** that embed some extra context into where the data for
this file comes from. Feel free to edit the file and add your metadata.

The key parameters for the disease are the set of numbers associated
with **beta**, **progress**, **too_ill_to_move** and **contrib_foi**.

These lists provide the values for those parameters in the model for the
different stages of the disease. In this case, the lurgy has five
stages (also called classes). An individual who is infected progresses
along each stage (or class) in turn according to the value of the parameters;

* **beta** - the rate of infection of the disease. For the first two stages
  the value of **beta** is zero, meaning that infected individuals are
  not infectious. However, individuals are infectious in stages three
  and four, when **beta** is 0.5. The final fifth stage is when the individual
  is removed from the outbreak, when the **beta** value is zero again.

* **progess** - rate of progression of an individual through the stages
  of the disease. Individuals progress quickly through the first two
  stages (where **progress** is 1.0), and then progress slows down
  in the third and fourth stages (**progress** is 0.5). The fifth stage
  is the final stage, and so the **progress** value is zero.

* **too_ill_to_move** - the proportion of infected individuals who are
  unable to move. This starts low in the initial stages, e.g. at
  stage three individuals can move (**too_ill_to_move** is 0.5) and
  are quite infectious (**beta** is 0.5). This progresses up to 1.0
  for the final stage, when the individual is removed from the oubreak.

* **contrib_foi** - this is how much the infected individuals contribute
  to the force of infection in their home electoral ward or during
  their work commute. In this case, individuals infected with
  this model lurgy contribute fully (**contrib_foi** is 1.0) for all
  but the final stage, when they are removed from the outbreak.

Modifying the disease file
--------------------------

With the current model, the lurgy only becomes infectious when
the infected individual are already showing symptoms and have
reduced mobility. For the lurgy that we want to model we need to
represent the lurgy as a disease that is infectious for a period
of time before the individual is showing symptoms. To this end,
we will an extra stage in the middle for which **beta** is high
but **too_ill_to_move** is still zero. Edit your copy of
``lurgy2.json`` to read;

::

  { "name"             : "The Lurgy",
    "version"          : "April 16th 2020",
    "author(s)"        : "Christopher Woods",
    "contact(s)"       : "christopher.woods@bristol.ac.uk",
    "reference(s)"     : "Completely ficticious disease - no references",
    "beta"             : [0.0, 0.0, 0.5, 0.5, 0.5, 0.0],
    "progress"         : [1.0, 1.0, 1.0, 0.5, 0.5, 0.0],
    "too_ill_to_move"  : [0.0, 0.0, 0.0, 0.5, 0.8, 1.0],
    "contrib_foi"      : [1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
  }

Once you have saved the file you can perform a single *model run*
using this file via;

.. code-block:: bash

   metawards -d lurgy2 -a ExtraSeedsLondon.dat

.. note::

  ``metawards`` automatically used your ``lurgy2.json`` file as it
  found the file in your current directory. You can pass a full path to
  your file, with or without the ``.json`` extension

This will run the ``metawards`` *model run* for your new version of the
lurgy. Notice that near the top of the output you have your parameters
and the metadata printed to the screen, e.g.

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Disease ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    • Disease: lurgy2
    • loaded from: lurgy2.json
    • repository: /Users/chris/GitHub/MetaWardsData
    • repository_branch: None
    • repository_version: None
    • beta: [0.0, 0.0, 0.5, 0.5, 0.5, 0.0]
    • progress: [1.0, 1.0, 1.0, 0.5, 0.5, 0.0]
    • too_ill_to_move: [0.0, 0.0, 0.0, 0.5, 0.8, 1.0]
    • contrib_foi: [1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
    • start_symptom: 3

Again, this helps someone reproduce this output in the future.

This model run may take longer, as, intuitively, you would expect that
the changes we have made mean that more individuals are likely to be
infected. Indeed, for the run I performed, copied below,
the outbreak lasted for 196 days, involving nearly 50m individuals.

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 189 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 6225870  E: 1  I: 3  R: 49856203  IW: 0  POPULATION: 56082077
    Number of infections: 4

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 190 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 6225870  E: 0  I: 2  R: 49856205  IW: 0  POPULATION: 56082077
    Number of infections: 2

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 191 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 6225870  E: 0  I: 2  R: 49856205  IW: 0  POPULATION: 56082077
    Number of infections: 2

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 192 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 6225870  E: 0  I: 2  R: 49856205  IW: 0  POPULATION: 56082077
    Number of infections: 2

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 193 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 6225870  E: 0  I: 2  R: 49856205  IW: 0  POPULATION: 56082077
    Number of infections: 2

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 194 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 6225870  E: 0  I: 2  R: 49856205  IW: 0  POPULATION: 56082077
    Number of infections: 2

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 195 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 6225870  E: 0  I: 1  R: 49856206  IW: 0  POPULATION: 56082077
    Number of infections: 1

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 196 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 6225870  E: 0  I: 0  R: 49856207  IW: 0  POPULATION: 56082077
    Number of infections: 0
    Infection died ... Ending on day 197


An overview plot of the outbreak, created using

.. code-block:: bash

  metawards-plot -i output/results.csv.bz2

shows the much higher peak of the outbreak compared to our original
model of the lurgy.

.. image:: ../../images/tutorial_2_1_overview.jpg
   :alt: Overview image of the outbreak of second version of the lurgy
====================
Running on a cluster
====================

You have now used ``metawards`` to perform one *model run* for several
different combinations of the **beta** and **too_ill_to_move**
disease parameters for your model of the lurgy.

It is important to repeat each *model run* several times in order to
reduce random error and interpret the change in the population trajectory
in response to changes in disease parameters.

This is supported using the ``--repeats`` command line argument. You
can use this to specify the number of times you want to repeat each
*model run*, e.g. ``--repeats 8`` would repeat each **model run**
eight times.

.. note::
   You can pass multiple values to ``--repeats`` if you want different
   adjustable variable sets to be repeated different numbers of times.
   Equally, you can add a ``repeats`` column to the ``lurgyparams.csv``
   file that gives the number of times each line should be repeated.
   If you use both options, then the repeats are multipled, e.g. if
   ``repeats`` is ``2`` in the model file, and ``--repeats 3`` is passed
   as a command-line argument, then the adjustable variables on that
   line of the model file will be repeated 2 x 3 == 6 times.

To this end, ``metawards`` natively supports running across multiple
compute nodes of a
`slurm <https://slurm.schedmd.com>`__ or
`PBS-style <https://en.wikipedia.org/wiki/Portable_Batch_System>`__
High Performance Computing (HPC) cluster.

Installing on a cluster
-----------------------

The first thing to do is to get ``metawards`` installed on your cluster.
A good option is to use a Python environment, as this should help make
it easier to return to a previous ``metawards`` installation if you
need to repeat a job.
:doc:`Take a look here <../../development>` to learn
how to install ``metawards`` into a Python environment.

Setting up the job
------------------

Next, create a directory for your cluster job, and into this copy your
``lurgy2.json`` file. You could also copy your ``lurgyparams.csv`` file,
but we will take the opportunity of running on a cluster to run a more
fine-grained parameter sweep. Rather than write the ``lurgyparams.csv``
file by hand, we will write now a simple script that can generate it
for us.

Create a file called ``create_params.py`` and copy in the below;

.. code-block:: python

  import sys

  b0 = float(sys.argv[1])
  b1 = float(sys.argv[2])
  bdel = float(sys.argv[3])

  i0 = float(sys.argv[4])
  i1 = float(sys.argv[5])
  idel =float(sys.argv[6])

  print("beta[2]  too_ill_to_move[2]")

  b = b0
  while b <= b1:
      i = i0
      while i <= i1:
          print("  %.2f     %.2f" % (b, i))
          i += idel
      b += bdel

Run this script using;

.. code-block:: bash

  python create_params.py 0.3 0.71 0.05 0.0 0.5 0.05 > lurgyparams.csv

This will create ``lurgyparams.csv`` that describes 99 model runs. These
will vary **beta** between 0.3 to 0.7 inclusive, in steps of 0.05, while
also varying **too_ill_to_move** between 0.0 to 0.5 inclusive, also
in steps of 0.05. The first few lines of this file are shown below;

::

  beta[2]  too_ill_to_move[2]
    0.30     0.00
    0.30     0.05
    0.30     0.10
    0.30     0.15
    0.30     0.20
    0.30     0.25
    0.30     0.30
    0.30     0.35
    0.30     0.40

Writing a job script
--------------------

We now need to write a job script that will submit a run the job to the
cluster queueing system.
:doc:`Example job scripts for SLURM and PBS are here <../../cluster_usage>`.

I am running on the `Catalyst ARM64 cluster <https://www.bristol.ac.uk/news/2018/april/supercomputer-collaboration.html>`__,
which uses PBS. The ``metawards`` command I need is very similar to before,
but now I am going to run 16 repeats, use 8 cores per *model run*, and
will force the overwriting of output to make sure that my jobs don't
hang on a prompt. The job-script I used, which I called ``jobscript.sh``,
is copied here;

.. code-block:: bash

  #!/bin/bash
  #PBS -l walltime=12:00:00
  #PBS -l select=4:ncpus=64:mem=64GB

  # source the version of metawards we want to use
  source $HOME/envs/metawards-devel/bin/activate

  # change into the directory from which this job was submitted
  cd $PBS_O_WORKDIR

  export METAWARDS_CORES_PER_NODE="64"
  export METAWARDSDATA="$HOME/GitHub/MetaWardsData"

  metawards --additional ExtraSeedsLondon.dat \
            --disease lurgy2.json \
            --input lurgyparams.csv --repeats 16 --nthreads 8 \
            --force-overwrite-output

The ``PBS`` commands at the top tell the queueing system that I want to run
for a maximum of 12 hours using four 64-core nodes (256 cores in total).

I've then activated my ``metawards-devel`` python environment that was in
``$HOME/envs/metawards-devel``.

To help distribute work, ``metawards`` needs to know how many cores there
are on each compute nodes. This is set using the
``METAWARDS_CORES_PER_NODE`` environment variable (or alternatively could
be passed using the ``--cores-per-node`` command-line argument).
I've also used the ``METAWARDSDATA`` environment variable to locate
the MetaWardsData repository data.

You may have to modify this script for your cluster and queueing system.

Running the HPC job
-------------------

Once you have written the job script, you should submit it using your
job submission command. As I used a PBS cluster, I used;

.. code-block:: bash

   qsub jobscript.sh

I could then check the status of the job using

.. code-block:: bash

   qstat -n

Processing the output
---------------------

The job will take a while. 99 *model runs* with 16 repeats each is
1584 total runs, so you may want to go to lunch or leave this running
overnight.

In my case, the job took 2 hours in total to run. Once complete, the
``results.csv.bz2`` file contains all of the population trajectories
and can be analysed in an identical way as before. If you want, you can
:download:`my results.csv.bz2 file here <output1/results.csv.bz2>`.

You can then produce graphs and animations using;

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2 --format jpg --dpi 150
   metawards-plot --animate output/overview*.jpg

The resulting animation of the overview plots is shown below.

.. image:: ../../images/tutorial_2_4.gif
   :alt: Overview animation of the outbreak of the lurgy

==========
Model data
==========

``metawards`` is a geographical SIR model. This means that a separate
progression of states from **S** through to **R** is performed in every
*ward* that is modelled. A *ward* is just a geographic area or cell.
It originally represented a single electoral ward in the UK (hence the name),
but there is no technical reason why this has to be the case.

Infected individuals who are resident in a *ward* contribute to the
*force of infection* (FOI) of that ward. The higher the FOI, the more
likely it is that another resident will be infected.

In addition, some residents (called *workers*) make regular movements
each day between their home *ward* and a work *ward*. They contribute
to and are exposed to the FOI of their *work* ward, and so can
spread the infection between wards.

The number of wards, geographic data, populations and population flows
between them are altogether called the *model data*. ``metawards`` was
first developed to model electoral wards in England and Wales, with
data based on the 2011 census. The default model used in ``metawards``
is thus this ``2011Data`` model.

You can set the model used by ``metawards`` using the ``--model`` or ``-m``
command line argument, e.g.

.. code-block:: bash

   metawards -d lurgy3 -m 2011Data

would perform a *model run* using your new lurgy3 parameters on the default
``2011Data`` model data. We have been, and will continue to use
this ``2011Data`` default model for
all of the examples in this tutorial. However, please remember that
you can create your own *model data* to represent any location, region,
country or interacting sets of groups as you wish.

Alternative models
------------------

Constructing this *model data* involves a lot of data science and careful
consideration. The models are held in the
`MetaWardsData <https://github.com/metawards/MetaWardsData>`__ repository,
and the current models include;

* ``2011Data`` : default data for England and Wales from the 2011 census

* ``2011UK`` : work-in-progress model for the whole of the UK including Northern Ireland.

In addition, there is a special model called ``single``, that is used to
perform validation *model runs* where the entire population is resident
in only a single ward. This is useful if you want to validate the underlying
models used in ``metawards``, or want to see if results obtained are
influenced by geography.

You should specify the size of the population when you use the ``single``
model via the ``--population`` or ``-P`` keyword, e.g.

.. code-block:: bash

   metawards -d lurgy3 -a ExtraSeedsOne.dat -m single -P 50000000 --nsteps 2000

would perform a run with 50 million individuals in a single ward.

.. note::

    We have used ``ExtraSeedsOne.dat`` to seed this run, as this file
    seeds 5 infections into the first ward of a model on day 1. We've
    also set ``nsteps`` to 2000 to stop the model timing out after
    two years. This is because the lack of geographic mixing significantly
    slows the progression of the disease.

After running this command, you can create an overview plot using;

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2

The output from your run should look something like this;

.. image:: ../../images/tutorial_2_6.jpg
   :alt: Overview from a run using a single ward

.. note::

    You may find that the disease dies out quickly in a small number of runs.
    This is due to the random nature of the model, and random chance means
    that the disease in those runs could not catch light. If it does catch,
    then the epidemic should progress for ~1000 days.
====================
Adjustable Variables
====================

We have now adjusted the **disease file** for the lurgy, which has
resulted in a *model run* that has involved a larger section of the
population in the disease outbreak. However, the values we chose
for the parameters in the **disease file** were made up. They had
no basis in reality, which means that it is difficult to use this
model to make good predictions.

``metawards`` is designed to support parameter searches that seek to
adjust the disease parameters so that predictions from *model runs*
can be matched to observed reality. The idea is that variables are
adjusted until the *model runs* can reliably capture the observed
behaviour of an outbreak, at which point the model may then act
as a good prediction of the future.

Choosing what to adjust
-----------------------

In the last section we introduced a new stage to the lurgy in which
infected individuals were very mobile (**too_ill_to_move** was 0.0),
but also very infectious (**beta** was 0.5).

This was introduced as stage 3 of the lurgy. We need therefore to tune
the **beta[2]** and **too_ill_to_move[2]** parameters (remembering that
we count from zero). What we want to do is vary these two parameters
and see how they affect the model outbreak.

To do this, in your current directory create a file called ``lurgyparams.csv``
and copy in the text below;

::

  beta[2]  too_ill_to_move[2]
    0.3         0.00
    0.4         0.00
    0.5         0.00
    0.3         0.25
    0.4         0.25
    0.5         0.25
    0.3         0.50
    0.4         0.50
    0.5         0.50

This file has two columns of numbers; ``beta[2]`` which shows how to
adjust the **beta[2]** value, and ``too_ill_to_move[2]`` which shows
how to adjust the **too_ill_to_move[2]** value. Each row of this file
represents a single model run using this pair of values. For example,
the line ``0.3       0.00`` will set **beta[2]** to 0.3 and
**too_ill_to_move[2]** to 0.0.

.. note::
   The format of this file is very flexible. Columns can be space-separated
   or comma-separated, as long as this is consistent in the entire file.
   You can choose to adjust any of the parameters for any stage
   of a disease. Simply title the column with the parameter name, e.g.
   ``beta``, ``progress``, ``too_ill_to_move`` or ``contrib_foi``, and
   add the stage number in square brackets, e.g. ``beta[0]``. Remember
   that we count from zero, so stage one has index zero.

You can pass ``metawards`` this input file using the ``--input``
(or ``-i``) command-line option. This input file has nine pairs
of values of **beta[2]** and **too_ill_to_move[2]**, which means that
nine *model runs* need to be performed. To run them all, use this command;

.. code-block:: bash

   metawards -d lurgy2 -a ExtraSeedsLondon.dat --input lurgyparams.csv

This will run all nine jobs. If you have nine or more processor cores
on your computer then all of them will be run in parallel (with
individual *model runs* then using any other cores you have). If, like me,
you are running this on your laptop, then this may take 10+ minutes
to complete.

.. note::
  If you want to distribute this work over a set of disconnected
  computers, you can tell ``metawards`` to only adjust parameters
  from a subset of the lines of the input file. To do this,
  use the ``--line`` (or ``-l``) command line argument to specify
  the line or lines to process. Lines are counted from 0 being the
  top line (containing the header), and multiple lines or ranges
  can be specified, e.g. ``-l 1 2 3`` will use lines one to three,
  while ``-l 4-6`` would use lines four to six (inclusive).

.. warning::
  We are only going to use one repeat of each pair of values, which means
  that our results will suffer from a lot of random error. Ideally you
  should ask for multiple repeats using the ``--repeats`` command-line
  argument. A good value would be at least eight repeats. In this case,
  eight repeats would require 72 *model runs*. If you parallelise
  each run over 16 cores then this needs 1152 cores. Fortunately, for
  these jobs, ``metawards`` runs well in parallel
  :doc:`across a slurm or PBS cluster <../../cluster_usage>`.

Once the jobs are complete (which took 15 minutes on my laptop), you
should have output that looks similar to this;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ MULTIPROCESSING ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 1 of 9                                                                  │
    │  (beta[2]=0.3, too_ill_to_move[2]=0)[repeat 1]                                         │
    │  2020-12-15: DAY: 209 S: 7978312    E: 0    I: 0    R: 48103765    IW: 1   UV: 1.0     │
    │  TOTAL POPULATION 56082077                                                             │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 2 of 9                                                                  │
    │  (beta[2]=0.4, too_ill_to_move[2]=0)[repeat 1]                                         │
    │  2020-11-28: DAY: 192 S: 7046459    E: 0    I: 0    R: 49035618    IW: 0   UV: 1.0     │
    │  TOTAL POPULATION 56082077                                                             │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 3 of 9                                                                  │
    │  (beta[2]=0.5, too_ill_to_move[2]=0)[repeat 1]                                         │
    │  2020-11-17: DAY: 181 S: 6221586    E: 0    I: 0    R: 49860491    IW: 1   UV: 1.0     │
    │  TOTAL POPULATION 56082077                                                             │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 4 of 9                                                                  │
    │  (beta[2]=0.3, too_ill_to_move[2]=0.25)[repeat 1]                                      │
    │  2020-12-13: DAY: 207 S: 8010218    E: 0    I: 0    R: 48071859    IW: 0   UV: 1.0     │
    │  TOTAL POPULATION 56082077                                                             │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 5 of 9                                                                  │
    │  (beta[2]=0.4, too_ill_to_move[2]=0.25)[repeat 1]                                      │
    │  2020-12-02: DAY: 196 S: 7071208    E: 0    I: 0    R: 49010869    IW: 1   UV: 1.0     │
    │  TOTAL POPULATION 56082077                                                             │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 6 of 9                                                                  │
    │  (beta[2]=0.5, too_ill_to_move[2]=0.25)[repeat 1]                                      │
    │  2020-12-01: DAY: 195 S: 6260263    E: 0    I: 0    R: 49821814    IW: 0   UV: 1.0     │
    │  TOTAL POPULATION 56082077                                                             │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 7 of 9                                                                  │
    │  (beta[2]=0.3, too_ill_to_move[2]=0.5)[repeat 1]                                       │
    │  2021-01-01: DAY: 226 S: 8031161    E: 0    I: 0    R: 48050916    IW: 0   UV: 1.0     │
    │  TOTAL POPULATION 56082077                                                             │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 8 of 9                                                                  │
    │  (beta[2]=0.4, too_ill_to_move[2]=0.5)[repeat 1]                                       │
    │  2020-11-27: DAY: 191 S: 7103861    E: 0    I: 0    R: 48978216    IW: 0   UV: 1.0     │
    │  TOTAL POPULATION 56082077                                                             │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Completed job 9 of 9                                                                  │
    │  (beta[2]=0.5, too_ill_to_move[2]=0.5)[repeat 1]                                       │
    │  2020-11-16: DAY: 180 S: 6293958    E: 0    I: 0    R: 49788119    IW: 0   UV: 1.0     │
    │  TOTAL POPULATION 56082077                                                             │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘
    ┌────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                        │
    │  Writing a summary of all results into the csv file                                    │
    │  /Users/chris/GitHub/tutorial/test/output/results.csv.bz2. You can use this to         │
    │  quickly look at statistics across all runs using e.g. R or pandas                     │
    │                                                                                        │
    └────────────────────────────────────────────────────────────────────────────────────────┘


Output directories
------------------

The output files for each repeat are placed into subdirectories of the
main output directory. By default, these subdirectories are named
according to the fingerprint of the adjustable variables used for each
run, e.g. listing the contents of the output directory using;

.. code-block:: bash

   $ ls output
   0i3v0i0x001     0i3v0i5x001     0i4v0i25x001    0i5v0i0x001     0i5v0i5x001     console.log.bz2
   0i3v0i25x001    0i4v0i0x001     0i4v0i5x001     0i5v0i25x001    config.yaml     results.csv.bz2

The fingerprint is a unique key used for each run, e.g.
``0i3v0i0x001`` refers to the run using values ``0.3 0.0``, and the
first repeat. The ``i`` represents a decimal point, ``v`` is used to
separate values, and ``x001`` means the first repeat.

Similarly, ``0i4v0i25x001`` refers to the first repeat of the values
``0.4 0.25``.

Sometimes you may want to specify the names of the output directories
yourself. You can do this by adding a ``output`` column to your scan file,
e.g.

::

  beta[2]  too_ill_to_move[2]      output
    0.3         0.00           beta_0i3_ill_0i00
    0.4         0.00           beta_0i4_ill_0i00
    0.5         0.00           beta_0i5_ill_0i00
    0.3         0.25           beta_0i3_ill_0i25
    0.4         0.25           beta_0i4_ill_0i25
    0.5         0.25           beta_0i5_ill_0i25
    0.3         0.50           beta_0i3_ill_0i50
    0.4         0.50           beta_0i4_ill_0i50
    0.5         0.50           beta_0i5_ill_0i50

Running ``metawards`` with this file would place output in the following
directories;

.. code-block:: bash

   $ ls output
   beta_0i3_ill_0i00 beta_0i4_ill_0i00 beta_0i5_ill_0i00 config.yaml
   beta_0i3_ill_0i25 beta_0i4_ill_0i25 beta_0i5_ill_0i25 console.log.bz2
   beta_0i3_ill_0i50 beta_0i4_ill_0i50 beta_0i5_ill_0i50 results.csv.bz2

If you run multiple repeats of these jobs, e.g. using the ``--repeats``
keyword via;

.. code-block:: bash

   metawards -d lurgy2 -a ExtraSeedsLondon.dat --input lurgyparams.csv --repeats 2

then the repeat number will be automatically added to the directory names,
e.g.

.. code-block:: bash

   $ ls output
   beta_0i3_ill_0i00     beta_0i4_ill_0i00     beta_0i5_ill_0i00     config.yaml
   beta_0i3_ill_0i00x002 beta_0i4_ill_0i00x002 beta_0i5_ill_0i00x002 console.log.bz2
   beta_0i3_ill_0i25     beta_0i4_ill_0i25     beta_0i5_ill_0i25     results.csv.bz2
   beta_0i3_ill_0i25x002 beta_0i4_ill_0i25x002 beta_0i5_ill_0i25x002
   beta_0i3_ill_0i50     beta_0i4_ill_0i50     beta_0i5_ill_0i50
   beta_0i3_ill_0i50x002 beta_0i4_ill_0i50x002 beta_0i5_ill_0i50x002
====================
Analysing the output
====================

The ``results.csv.bz2`` file contains all of the population trajectories
from the nine model runs. You can explore this using Python pandas, R,
or Excel :doc:`as you did before <../part01/02_repeating>`. Using
ipython or Jupyter notebooks with pandas, we can load up the file;

.. code-block:: python

   >>> import pandas as pd
   >>> df = pd.read_csv("output/results.csv.bz2")
   >>> print(df)
		fingerprint  repeat  beta[2]  too_ill_to_move[2]  day        date         S  E  I         R  IW   UV
	0        0i3v0i0       1      0.3                 0.0    0  2020-05-12  56082077  0  0         0   0  1.0
	1        0i3v0i0       1      0.3                 0.0    1  2020-05-13  56082077  0  0         0   0  1.0
	2        0i3v0i0       1      0.3                 0.0    2  2020-05-14  56082072  5  0         0   0  1.0
	3        0i3v0i0       1      0.3                 0.0    3  2020-05-15  56082072  0  5         0   0  1.0
	4        0i3v0i0       1      0.3                 0.0    4  2020-05-16  56082066  0  5         6   5  1.0
	...          ...     ...      ...                 ...  ...         ...       ... .. ..       ...  ..  ...
	1769     0i5v0i5       1      0.5                 0.5  169  2020-10-28   6302422  1  0  49779654   0  1.0
	1770     0i5v0i5       1      0.5                 0.5  170  2020-10-29   6302422  0  1  49779654   0  1.0
	1771     0i5v0i5       1      0.5                 0.5  171  2020-10-30   6302422  0  1  49779654   0  1.0
	1772     0i5v0i5       1      0.5                 0.5  172  2020-10-31   6302422  0  1  49779654   0  1.0
	1773     0i5v0i5       1      0.5                 0.5  173  2020-11-01   6302422  0  0  49779655   0  1.0

	[1774 rows x 12 columns]

This is very similar to before, except now we have extra columns giving
the values of the variables that are being adjusted (columns
``beta[2]`` and ``too_ill_to_move[2]``. We also now have a use for the
``fingerprint`` column, which contains a unique identifier for each
pair of adjustable variables.

.. note::
   In the fingerprint the ``i`` character represents a decimal point
   and ``v`` separates variables. For example ``0.3  0.0`` becomes
   ``0i3v0i0``, while ``0.5 0.5`` becomes ``0i5v0i5``.

Finding peaks
-------------

We can use ``.groupby`` to group the results with the same fingerprint
together. Then the ``.max`` function can be used to show the maximum
values of selected columns from each group, e.g.

.. code-block:: python

  >>> df.groupby("fingerprint")[["day", "E","I", "IW", "R"]].max()
				day        E         I    IW         R
	fingerprint
	0i3v0i0      223  1867602   9124999  8588  48107545
	0i3v0i25     213  1806030   8849753  8588  48080295
	0i3v0i5      203  1925482   9373635  8588  48050991
	0i4v0i0      191  1926941   9441320  8588  49044169
	0i4v0i25     209  2013614   9815751  8588  49007451
	0i4v0i5      196  2049927   9994566  8588  48979256
	0i5v0i0      177  2108287  10278260  8588  49861436
	0i5v0i25     180  2093600  10215078  8588  49814876
	0i5v0i5      173  2070325  10128212  8588  49779655

From this, we can see that higher peaks occured for higher values
of **beta**, which is expected. However, different values of
**too_ill_to_move** had little impact on the peaks.

.. warning::

  Do not over-interpret the results of single runs, such as the above.
  There is a lot of random
  error in these calculations and multiple *model runs* must be
  averaged over to gain a good understanding.

Plotting the output
-------------------

There are lots of plots you would likely want to draw, so it is recommended
that you use a tool such as R, Pandas or Excel to create the plots that
will let you explore the data in full. For a quick set of plots, you
can again use ``metawards-plot`` to generate some overview plots. To
do this type;

.. code-block:: bash

  metawards-plot -i output/results.csv.bz2 --format jpg --dpi 150

.. note::
   We have used the 'jpg' image format here are we want to create animations.
   You can choose from many different formats, e.g. 'pdf' for publication
   quality graphs, 'png' etc. Use the ``--dpi`` option to set the
   resolution when creating bitmap (png, jpg) images.

As there are multiple fingerprints, this will produce multiple overview
graphs (one overview per fingerprint, and if you have run multiple
repeats, then one average per fingerprint too).

The fingerprint value is included in the graph name, and they will
all be plotted on the same axes. This means that they could be joined
together into an animation. As well as plotting, ``metawards-plot`` has
an animation mode that can be used to join images together. To run this,
use;

.. code-block:: bash

  metawards-plot --animate output/overview*.jpg

.. note::
   You can only animate image files (e.g. jpeg, png). You can't animate
   pdfs (yet - although
   `pull requests welcome <https://github.com/metawards/MetaWards/pulls>`__).
   Also, animation relies on you installing
   `Pillow <https://pillow.readthedocs.io/en/stable/>`__ to create the gifs
   and (optional but recommended)
   `gifsicle <https://www.lcdf.org/gifsicle/>`__ and
   `pygifsicle <https://pypi.org/project/pygifsicle/>`__ to optimise the gifs
   (this reduces their size by 5-10 times)

Here is the animation.

.. image:: ../../images/tutorial_2_3.gif
   :alt: Animated overview graphs from the parameter sweep

Jupyter notebook
----------------

In addition, to the ``metawards-plot`` command, we also have a
:download:`Jupyter notebook <../../notebooks/2_3_analysis.ipynb>`
which you can look at which breaks down exactly how ``metawards-plot``
uses pandas and matplotlib to render multi-fingerprint graphs.
==================
Refining the model
==================

You have now run a large number of repeats of lots of different combinations
of the **beta** and **too_ill_to_move** disease parameters of your model
of the lurgy.

Looking at the results, it is clear the greater the value of **beta**,
the larger the outbreak. This makes sense, as **beta** controls the
the rate of infection.

The interpretation for **too_ill_to_move** is less clear. It doesn't
immediately seem to have a big impact on the size or duration of the
infection.

Testing a hypothesis
--------------------

One hypothesis why this was the case is because the **progress** parameter,
which represents the rate at which individuals move through our new
stage of the lurgy, was 1.0. This meant that individuals in the model
did not spend long in the infectious (**beta** == 1.0) but mobile and
symptom-free (**too_ill_to_move** == 0.0) state.

The hypothesis is that decreasing the rate at which individuals move
through this state should increase the outbreak as they will have more
time to unknowingly spread the infection.

To test this hypothesis, we can create a new set of adjustable variables
in which we vary **progress** as well as **beta** and **too_ill_to_move**.

Setting up the job
------------------

Create a new directory on your cluster and copy in your ``lurgy2.json``
disease file. Then create a new script called ``create_params.py`` and
copy into this;

.. code-block:: python

  import sys

  b0 = float(sys.argv[1])
  b1 = float(sys.argv[2])
  bdel = float(sys.argv[3])

  i0 = float(sys.argv[4])
  i1 = float(sys.argv[5])
  idel =float(sys.argv[6])

  p0 = float(sys.argv[7])
  p1 = float(sys.argv[8])
  pdel = float(sys.argv[9])

  print("beta[2]  too_ill_to_move[2]  progress[2]")

  b = b0
  while b <= b1:
      i = i0
      while i <= i1:
          p = p0
          while p <= p1:
              print("  %.2f     %.2f    %.2f" % (b, i, p))
              p += pdel

          i += idel
      b += bdel

Run this script using;

.. code-block:: bash

  python create_params.py 0.4 0.71 0.1 0.0 0.5 0.1 0.2 0.81 0.2 > lurgyparams.csv

This will create a set of adjustable variable sets for *model runs* that
vary **beta** between 0.4 and 0.7 in steps of 0.1,
**too_ill_to_move** between 0.0 and 0.5 in steps of 0.1, and
**progress** between 0.2 and 0.8 in steps of 0.2.

The first few lines of this file are shown below;

::

  beta[2]  too_ill_to_move[2]  progress[2]
    0.40     0.00    0.20
    0.40     0.00    0.40
    0.40     0.00    0.60
    0.40     0.00    0.80
    0.40     0.10    0.20
    0.40     0.10    0.40
    0.40     0.10    0.60
    0.40     0.10    0.80
    0.40     0.20    0.20

.. note::
  The file is space-separated, so it is not a problem that the columns
  don't line up with their titles. With a little more python you could
  make them line up. What approach would you take?

Running the job
---------------

You can now submit this job using a copy of the job script that you used
before. For example, I used;

.. code-block:: bash

  #!/bin/bash
  #PBS -l walltime=12:00:00
  #PBS -l select=4:ncpus=64:mem=64GB

  # source the version of metawards we want to use
  source $HOME/envs/metawards-devel/bin/activate

  # change into the directory from which this job was submitted
  cd $PBS_O_WORKDIR

  export METAWARDS_CORES_PER_NODE="64"
  export METAWARDSDATA="$HOME/GitHub/MetaWardsData"

  metawards --additional ExtraSeedsLondon.dat \
            --disease lurgy2.json \
            --input lurgyparams.csv --repeats 16 --nthreads 8 \
            --force-overwrite-output

I submitted usig the command;

.. code-block:: bash

  qsub jobscript.sh

Processing the output
---------------------

The job will take a while. In my case, the job took 2 hours to run.
Once complete, the ``results.csv.bz2`` file contains all of the
population trajectories and can be analysed in an identical way
as before. If you want, you can
:download:`my results.csv.bz2 file here <output2/results.csv.bz2>`.

You can then produce graphs and animations using;

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2 --format jpg --dpi 150
   metawards-plot --animate output/overview*.jpg

The resulting animation of the overview plots is shown below.

.. image:: ../../images/tutorial_2_5.gif
   :alt: Overview animation of the outbreak of the lurgy

Conclusion from experiment
--------------------------

It is clear from these graphs that the rate of progress through the new
stage of the lurgy we added has an impact on the form of the
outbreak. The lower value of **progress** in the new infectious and
asymptomatic stage leads to a wider outbreak that moves more quickly
through the population.

The (fictional) lurgy is a mild and mildly infectious disease that doesn't
cause too many symptoms. From the graphs, it is clear that this is best
modelled using a low value of **beta** and a low value of **progress**
for the new stage we added. We will leave the value of
**too_ill_to_move** at 0.0 to capture asymptomatic nature of this
stage. With this choice made,
please create a new version of the lurgy disease parameters called
``lurgy3.json`` and copy in the below;

::

  { "name"             : "The Lurgy",
    "version"          : "April 17th 2020",
    "author(s)"        : "Christopher Woods",
    "contact(s)"       : "christopher.woods@bristol.ac.uk",
    "reference(s)"     : "Completely ficticious disease - no references",
    "beta"             : [0.0, 0.0, 0.4, 0.5, 0.5, 0.0],
    "progress"         : [1.0, 1.0, 0.2, 0.5, 0.5, 0.0],
    "too_ill_to_move"  : [0.0, 0.0, 0.0, 0.5, 0.8, 1.0],
    "contrib_foi"      : [1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
  }

Once you have this, verify that you can run a ``metawards`` simulation,
and then plot the graphs using;

.. code-block:: bash

   metawards -d lurgy3.json --additional ExtraSeedsLondon.dat
   metawards-plot -i output/results.csv --format jpg --dpi 150

If it works, then you should obtain a graph that looks something like this;

.. image:: ../../images/tutorial_2_5_2.jpg
   :alt: Overview from a run using the refined parameters

.. warning::

  Remember that this is a completely fictional disease and is not
  related to any real outbreak. The graphs above are purely
  illustrative and chosen so that people following this tutorial
  can quickly see the impact of their work.
====================
Named Disease Stages
====================

In the last section you learned how to use demographics with different
disease stages to model hospital and ICU admissions. While this worked,
the calculation of statistics from the simulation was slightly hacky,
as the disease stages were still labelled ``E``, ``I`` and ``R``, when
we really wanted to refer to them as ``H1``, ``H2`` etc.

Fortunately, ``metawards`` supports custom naming of disease stages.
You can do this by adding a ``stage`` field to the disease file.

Simple example
--------------

For example, here is a simple disease file that uses stages ``A``,
``B`` and ``C``. Please create the file ``named.json`` and copy in the
below;

::

    {
        "stage"            : [ "A", "B", "C" ],
        "beta"             : [ 0.0, 0.5, 0.0 ],
        "progress"         : [ 1.0, 1.0, 0.0 ],
        "too_ill_to_move"  : [ 0.0, 0.0, 0.0 ],
        "contrib_foi"      : [ 1.0, 1.0, 0.0 ],
        "start_symptom"    : 1
    }

.. note::

   Note that we've not included the ``name``, ``author`` or other metadata
   fields as these are not needed for this simple example. These are optional
   fields. We recommend you include them when you want to publish a disease
   file.

This file defines three disease stages, called ``A``, ``B`` and ``C``.
The first stage (``A``) is not infectious, as ``beta["A"]`` is ``0.0``.
The infectious stage is ``B``, as ``beta["B"]`` is ``0.5``. The final
stage is ``C``, which is not infectious, and is where the disease ends
(``progress["C"]`` is ``0.0``).

Run ``metawards`` using this disease file via;

.. code-block:: bash

   metawards -a ExtraSeedsLondon.dat -d named.json --nsteps 20

You should see output similar to;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Loading additional seeds from /Users/chris/GitHub/MetaWardsData/extra_seeds/ExtraSeedsLondon.dat
    (1, 255, 5, None)
    S: 56082077  A: 0  B: 0  C: 0  IW: 0  POPULATION: 56082077

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    seeding play_infections[0][255] += 5
    S: 56082072  A: 0  B: 5  C: 0  IW: 0  POPULATION: 56082077
    Number of infections: 5

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082067  A: 5  B: 0  C: 5  IW: 2  POPULATION: 56082077
    Number of infections: 10

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082067  A: 0  B: 5  C: 5  IW: 0  POPULATION: 56082077
    Number of infections: 10

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082063  A: 4  B: 0  C: 10  IW: 3  POPULATION: 56082077
    Number of infections: 14

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 5 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082063  A: 0  B: 4  C: 10  IW: 0  POPULATION: 56082077
    Number of infections: 14

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 6 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082062  A: 1  B: 0  C: 14  IW: 1  POPULATION: 56082077
    Number of infections: 15

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 7 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082062  A: 0  B: 1  C: 14  IW: 0  POPULATION: 56082077
    Number of infections: 15

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 8 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082062  A: 0  B: 0  C: 15  IW: 0  POPULATION: 56082077
    Number of infections: 15

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 9 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082062  A: 0  B: 0  C: 15  IW: 0  POPULATION: 56082077
    Number of infections: 15

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 10 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082062  A: 0  B: 0  C: 15  IW: 0  POPULATION: 56082077
    Number of infections: 15

.. note::

   Note that the simulation gets stuck in the ``C`` state. This is because
   any individual who is not in ``S`` or ``R`` is counted as an infection,
   and so the 15 individuals in ``C`` are counted as infecteds. To prevent
   the model running forever we set the maximum number of days to
   20 via ``--nsteps 20``.

As you can see, the output now records movement from ``S`` to ``A``, ``B``
and then ``C``. This data is also recorded in the output files, e.g.

.. code-block:: python

   >> import pandas as pd
   >> df = pd.read_csv("output/results.csv.bz2")
   >> df.head()
      fingerprint  repeat  day        date         S  E  I  A  B   C  R  IW   UV
    0      REPEAT       1    0  2020-06-23  56082077  0  0  0  0   0  0   0  1.0
    1      REPEAT       1    1  2020-06-24  56082072  0  0  0  5   0  0   0  1.0
    2      REPEAT       1    2  2020-06-25  56082067  0  0  5  0   5  0   2  1.0
    3      REPEAT       1    3  2020-06-26  56082067  0  0  0  5   5  0   0  1.0
    4      REPEAT       1    4  2020-06-27  56082063  0  0  4  0  10  0   3  1.0
   >> df = pd.read_csv("output/trajectory.csv.bz2")
   >> df.head()
       day        date demographic         S  E  I  A  B   C  R  IW
    0    0  2020-06-23     overall  56082077  0  0  0  0   0  0   0
    1    1  2020-06-24     overall  56082072  0  0  0  5   0  0   0
    2    2  2020-06-25     overall  56082067  0  0  5  0   5  0   2
    3    3  2020-06-26     overall  56082067  0  0  0  5   5  0   0
    4    4  2020-06-27     overall  56082063  0  0  4  0  10  0   3

Additional columns have been added to the tables in these files for the
``A``, ``B`` and ``C`` states.

Sub-stages example
------------------

You can have multiple named sub-stages of each stage, e.g. instead of
having a single infectious ``B`` stage, you can have ``B1``, ``B2`` and
``B3``. The totals reported for a the ``B`` stage will be the sum of
the number of individuals in each sub-stage. For example, edit
``named.json`` to read;

::

    {
        "stage"            : [ "A", "B1", "B2", "B3", "C" ],
        "beta"             : [ 0.0, 0.2,  0.8,  0.1,  0.0 ],
        "progress"         : [ 1.0, 1.0,  1.0,  1.0,  0.0 ],
        "too_ill_to_move"  : [ 0.0, 0.0,  0.2,  0.8,  0.0 ],
        "contrib_foi"      : [ 1.0, 1.0,  1.0,  1.0,  0.0 ],
        "start_symptom"    : 1
    }

Here we've expanded the ``B`` stage into three infectious sub-stages
(``B1``, ``B2`` and ``B3``), similar to the three stages of the lurgy.

Run ``metawards`` using this disease file via;

.. code-block:: bash

   metawards -a ExtraSeedsLondon.dat -d named.json --nsteps 20

You should see in the output that the population of ``A``, ``B`` and ``C``
are summarised, e.g.

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Loading additional seeds from /Users/chris/GitHub/MetaWardsData/extra_seeds/ExtraSeedsLondon.dat
    (1, 255, 5, None)
    S: 56082077  A: 0  B: 0  C: 0  IW: 0  POPULATION: 56082077

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    seeding play_infections[0][255] += 5
    S: 56082072  A: 0  B: 5  C: 0  IW: 0  POPULATION: 56082077
    Number of infections: 5

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082071  A: 1  B: 5  C: 0  IW: 1  POPULATION: 56082077
    Number of infections: 6

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082067  A: 4  B: 6  C: 0  IW: 4  POPULATION: 56082077
    Number of infections: 10

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082066  A: 1  B: 5  C: 5  IW: 1  POPULATION: 56082077
    Number of infections: 11

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 5 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082064  A: 2  B: 6  C: 5  IW: 2  POPULATION: 56082077
    Number of infections: 13

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 6 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082060  A: 4  B: 7  C: 6  IW: 4  POPULATION: 56082077
    Number of infections: 17

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 7 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082058  A: 2  B: 7  C: 10  IW: 2  POPULATION: 56082077
    Number of infections: 19

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 8 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082053  A: 5  B: 8  C: 11  IW: 4  POPULATION: 56082077
    Number of infections: 24

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 9 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082053  A: 0  B: 11  C: 13  IW: 0  POPULATION: 56082077
    Number of infections: 24

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 10 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082049  A: 4  B: 7  C: 17  IW: 4  POPULATION: 56082077
    Number of infections: 28

These are also summarised in the ``output/results.csv.bz2`` and
``output/trajectory.csv.bz2`` files.

However, the actual populations in each individual stage are given in the
``play_infections.csv.bz2`` (play infections), ``work_infections.csv.bz2``
(work infections) and ``number_infected_wards.csv.bz2`` (number of infected
wards) files, e.g.

.. code-block:: python

   >>> import pandas as pd
   >>> df = pd.read_csv("output/total_infections.csv.bz2")
   >>> df.head()
       day  A  B1  B2  B3  C
    0    1  0   5   0   0  0
    1    2  1   0   5   0  0
    2    3  4   1   0   5  0
    3    4  1   4   1   0  5
    4    5  2   1   4   1  5
   >>> df = pd.read_csv("output/number_infected_wards.csv.bz2")
   >>> df.head()
       day  A  B1  B2  B3  C
    0    1  0   1   0   0  0
    1    2  1   0   1   0  0
    2    3  4   1   0   1  0
    3    4  1   4   1   0  1
    4    5  2   1   4   1  1

These files are very useful if you want to see, e.g. how many workers
are infected at each different stage on each day, or how many wards
have a population infected in the ``B1`` state on each day.

Scanning named stage parameters
-------------------------------

You can also use the name of a stage when scanning disease parameters.
For example, create a file called ``scan.dat`` and copy in the below;

::

    beta["B1"]  beta["B2"]
      0.2         0.7
      0.3         0.8

Hopefully you can see that this will adjust the ``beta`` parameters for
the ``B1`` and ``B2`` stages. You can run this file using;

.. code-block:: bash

    metawards -a ExtraSeedsLondon.dat -d named.json --nsteps 20 -i scan.dat

and should see that the specified variables are indeed scanned, e.g.

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Adjustable parameters to scan ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    • (beta["B1"]=0.2, beta["B2"]=0.7)[repeat 1]
    • (beta["B1"]=0.3, beta["B2"]=0.8)[repeat 1]

    [...]

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ MULTIPROCESSING ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                                            │
    │  Completed job 1 of 2                                                                                      │
    │  (beta["B1"]=0.2, beta["B2"]=0.7)[repeat 1]                                                                │
    │  2020-07-13: DAY: 20  S: 56081987  A: 6  B: 21  C: 63  IW: 6  UV: 1.0  TOTAL POPULATION 56082077           │
    │                                                                                                            │
    └────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
    Computing model run ✔
    ┌────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
    │                                                                                                            │
    │  Completed job 2 of 2                                                                                      │
    │  (beta["B1"]=0.3, beta["B2"]=0.8)[repeat 1]                                                                │
    │  2020-07-13: DAY: 20  S: 56081794  A: 47  B: 93  C: 143  IW: 43  UV: 1.0  TOTAL POPULATION 56082077        │
    │                                                                                                            │
    └────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

Mapping stages to summaries
---------------------------

By default, the population of a disease sub-stage is summed into a summary
value that has the same name (but missing the sub-stage number). So ``B1``,
``B2`` and ``B3`` sub-stages will accumulate into the ``B`` stage.

You can control this mapping via the ``mapping`` value in the disease file.
You can set a disease stage to map to any individual stage (e.g. you could
map ``B1`` to be ``B1`` only), to any grouped stage (e.g. you could map
``C`` to map to the grouped ``B`` stage), or to any of the standard mapped
stages (``E``, ``I``, ``R`` or ``*``).

For example, you could output every stage to the summary via;

::

    {
        "stage"            : [ "A", "B1", "B2", "B3", "C" ],
        "mapping"          : [ "A", "B1", "B2", "B3", "C" ],
        "beta"             : [ 0.0, 0.2,  0.8,  0.1,  0.0 ],
        "progress"         : [ 1.0, 1.0,  1.0,  1.0,  0.0 ],
        "too_ill_to_move"  : [ 0.0, 0.0,  0.2,  0.8,  0.0 ],
        "contrib_foi"      : [ 1.0, 1.0,  1.0,  1.0,  0.0 ],
        "start_symptom"    : 1
    }

Running ``metawards`` using this file will tell it to output every stage,
e.g.

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Loading additional seeds from /Users/chris/GitHub/MetaWardsData/extra_seeds/ExtraSeedsLondon.dat
    (1, 255, 5, None)
    S: 56082077  A: 0  B1: 0  B2: 0  B3: 0  C: 0  IW: 0  POPULATION: 56082077

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    seeding play_infections[0][255] += 5
    S: 56082072  A: 0  B1: 5  B2: 0  B3: 0  C: 0  IW: 0  POPULATION: 56082077
    Number of infections: 5

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082068  A: 4  B1: 0  B2: 5  B3: 0  C: 0  IW: 4  POPULATION: 56082077
    Number of infections: 9

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082064  A: 4  B1: 4  B2: 0  B3: 5  C: 0  IW: 3  POPULATION: 56082077
    Number of infections: 13

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082061  A: 3  B1: 4  B2: 4  B3: 0  C: 5  IW: 3  POPULATION: 56082077
    Number of infections: 16

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 5 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082056  A: 5  B1: 3  B2: 4  B3: 4  C: 5  IW: 4  POPULATION: 56082077
    Number of infections: 21

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 6 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082046  A: 10  B1: 5  B2: 3  B3: 4  C: 9  IW: 9  POPULATION: 56082077
    Number of infections: 31

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 7 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082044  A: 2  B1: 10  B2: 5  B3: 3  C: 13  IW: 2  POPULATION: 56082077
    Number of infections: 33

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 8 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082036  A: 8  B1: 2  B2: 10  B3: 5  C: 16  IW: 8  POPULATION: 56082077
    Number of infections: 41

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 9 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082018  A: 18  B1: 8  B2: 2  B3: 10  C: 21  IW: 14  POPULATION: 56082077
    Number of infections: 59

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 10 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082010  A: 8  B1: 18  B2: 8  B3: 2  C: 31  IW: 8  POPULATION: 56082077
    Number of infections: 67

Alternatively, you can map your named stages to standard named accumulators,
e.g.

::

    {
        "stage"            : [ "A", "B1", "B2", "B3", "C" ],
        "mapping"          : [ "E", "I",  "I",  "I",  "R" ],
        "beta"             : [ 0.0, 0.2,  0.8,  0.1,  0.0 ],
        "progress"         : [ 1.0, 1.0,  1.0,  1.0,  0.0 ],
        "too_ill_to_move"  : [ 0.0, 0.0,  0.2,  0.8,  0.0 ],
        "contrib_foi"      : [ 1.0, 1.0,  1.0,  1.0,  0.0 ],
        "start_symptom"    : 1
    }

would count ``A`` as a latent ``E`` stage, ``B1`` to ``B3`` would be
infected ``I`` stages, and ``C`` would be accumulated as a ``R`` stage.

Running with this file would give;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Loading additional seeds from /Users/chris/GitHub/MetaWardsData/extra_seeds/ExtraSeedsLondon.dat
    (1, 255, 5, None)
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    seeding play_infections[0][255] += 5
    S: 56082072  E: 0  I: 5  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 5

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082072  E: 0  I: 5  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 5

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082068  E: 4  I: 5  R: 0  IW: 4  POPULATION: 56082077
    Number of infections: 9

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082068  E: 0  I: 4  R: 5  IW: 0  POPULATION: 56082077
    Number of infections: 4

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 5 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082068  E: 0  I: 4  R: 5  IW: 0  POPULATION: 56082077
    Number of infections: 4

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 6 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082064  E: 4  I: 4  R: 5  IW: 4  POPULATION: 56082077
    Number of infections: 8

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 7 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082064  E: 0  I: 4  R: 9  IW: 0  POPULATION: 56082077
    Number of infections: 4

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 8 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082064  E: 0  I: 4  R: 9  IW: 0  POPULATION: 56082077
    Number of infections: 4

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 9 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082062  E: 2  I: 4  R: 9  IW: 2  POPULATION: 56082077
    Number of infections: 6

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 10 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56082061  E: 1  I: 2  R: 13  IW: 1  POPULATION: 56082077
    Number of infections: 3

while the original stage names are still accessible in the
``output/total_infections.csv.bz2``, ``output/number_wards_infected.csv.bz2``
files etc.

One advantage of doing this is that now, ``C`` is correctly interpreted
as an ``R`` state, and so ``metawards`` will exit correctly once the
outbreak has died out and all individuals are left in ``S`` or ``C``.
=================
Modelling the ICU
=================

We have modelled the hospital as single patient population. However,
for a small proportion of cases, we need to model an extra stage
in the hospital for patients that need to go to the
intensive care unit (ICU).

We will use the ``H2`` state to model this. We will model that 10%
of patients in the ``H2`` state go to the ``H2`` state of an ``ICU``
demographic, while the remainder are immediately moved back home
to the ``R`` state, as they will have recovered from the disease.

To do this, modify your ``demographics.json`` file to read;

::

    {
        "demographics" : ["home", "staff", "patients", "icu"],
        "work_ratios"  : [ 0.99,   0.01,     0.00,     0.00 ],
        "play_ratios"  : [ 1.00,   0.00,     0.00,     0.00 ],
        "diseases"     : [ null,   null,   "lurgy_hospital", "lurgy_hospital" ]
    }

.. note::

   Now you can see the reason for the ``H2`` state - it is being used to
   model the ICU. Using two states like this enables us to use the
   same ``lurgy_hospital`` disease file to model the hospital.

Next, we need to modify the ``mix_hospital.py`` to read;

.. code-block:: python

    from metawards.mixers import merge_using_matrix

    def mix_shield(network, **kwargs):
        matrix = [ [1.0, 1.0, 0.0, 0.0],
                   [0.0, 0.1, 0.1, 0.1],
                   [0.0, 0.1, 0.0, 0.0],
                   [0.0, 0.1, 0.0, 0.0] ]

        network.demographics.interaction_matrix = matrix

        return [merge_using_matrix]

.. note::

  The only change here is adding in the fourth row and column for the ICU
  population. They don't contribute to the FOI of each other, other patients
  or home, but can infect and be infected by the hospital staff.

.. note::

   Note that we are using :func:`~metawards.mixers.merge_using_matrix`.
   This may not be the right choice depending on how we want the
   population dynamics to mix, e.g.
   :func:`~metawards.mixers.merge_matrix_single_population` or
   :func:`~metawards.mixers.merge_matrix_multi_population` may
   be a better choice. :doc:`See here for more information <../part05/04_contacts>`.

Next we need to modify ``extract_hospital.py`` to obtain the ICU
statistics. Edit the file and copy in the below;

.. code-block:: python

    from metawards.extractors import extract_default


    def output_patients(network, population, workspace, output_dir, **kwargs):
        # Open the file "patients.csv" in the output directory,
        # using the supplied headers for the columns
        FILE = output_dir.open("patients.csv",
                               headers=["day", "H1", "H2", "ICU"],
                               sep=",")

        # Now get the workspace for the "patients" demographic
        index = network.demographics.get_index("patients")
        subspace = workspace.subspaces[index]

        # The total population at each infection stage is the sum
        # of the work and play infections
        inf_tot = [inf + pinf for inf, pinf in
                   zip(subspace.inf_tot, subspace.pinf_tot)]

        H1 = inf_tot[2]
        H2 = inf_tot[3]

        # Now get the ICU demographic
        index = network.demographics.get_index("icu")
        subspace = workspace.subspaces[index]

        inf_tot = [inf + pinf for inf, pinf in
                   zip(subspace.inf_tot, subspace.pinf_tot)]

        ICU = inf_tot[3]

        FILE.write(str(population.day) + ",")
        FILE.write(",".join([str(x) for x in [H1, H2, ICU]]) + "\n")


    def extract_patients(**kwargs):
        # return all of the functions from "extract_default"
        # plus our new "output_i1"
        funcs = extract_default(**kwargs)
        funcs.append(output_patients)
        return funcs

.. note::

   The change here is that we extract only ``H1`` and ``H2`` from the
   ``patients`` demographic, before getting what we will call ``ICU``
   from the ``icu`` demographic.

Multiple go functions go home
-----------------------------

Finally, we will now update the ``move_hospital.py`` file so that
we will have four "go functions":

* First we start with the function that moves 20% of the ```home``
  and ``staff`` ``I2`` population to ``H1`` patients.

* Next, we move 10% of the ``H2`` ``patients`` to the same stage
  in the ``icu`` demographic. We will refer to this as ``ICU``.

* Next, we move the remainder of ``H2`` ``patients`` to ``R`` in ``home``,
  as these patients have now fully recovered and can go home.

* Finally, we move all ``R`` ``patients`` and ``icu`` members to ``R``
  in ``home`` as they have fully recovered and can go home.

You can implement this by editing your ``move_hospital.py`` file and
copying in the below;

.. code-block:: python

    from metawards.movers import go_stage


    def move_hospital(**kwargs):
        # move 20% of I2 home/staff population to H1 patients
        func1 = lambda **kwargs: go_stage(go_from=["home", "staff"],
                                          go_to="patients",
                                          from_stage=4,
                                          to_stage=2,
                                          fraction=0.2,
                                          **kwargs)

        # move 10% of H2 patients to H2 ICU
        func2 = lambda **kwargs: go_stage(go_from="patients",
                                          go_to="icu",
                                          from_stage=3,
                                          to_stage=3,
                                          fraction=0.1,
                                          **kwargs)

        # move the remainder of H2 patients to home R
        func3 = lambda **kwargs: go_stage(go_from="patients",
                                          go_to="home",
                                          from_stage=3,
                                          to_stage=-1,
                                          fraction=1.0,
                                          **kwargs)

        # move R ICU and patients to home R
        func4 = lambda **kwargs: go_stage(go_from=["patients", "icu"],
                                          go_to="home",
                                          from_stage=-1,
                                          to_stage=-1,
                                          fraction=1.0,
                                          **kwargs)

        return [func1, func2, func3, func4]

You can then run ``metawards`` using the command;

.. code-block:: bash

   metawards -D demographics.json -d lurgy4 --mixer mix_hospital --mover move_hospital --extract extract_hospital -a ExtraSeedsLondon.dat

You should see patients arriving in hospital, with some moving to the ICU.
By the end of the outbreak everyone has recovered and has returned home.

You can plot the demographics trajectory using;

.. code-block:: bash

   metawards-plot -i output/trajectory.csv.bz2

You should see a plot similar to this;

.. image:: ../../images/tutorial_7_4_1.jpg
   :alt: Demographic trajectories for the simple hospital plus ICU model

The ICU population is just visible on this plot, and is seen to lag behind
the patient population. You can see this more clearly by plotting the data
that was output to the ``output/patients.csv.bz2`` file, e.g. using
pandas;

.. code-block:: python

   >>> import pandas as pd
   >>> df = pd.read_csv("output/patients.csv.bz2")
   >>> df.plot(x="day")
   >>> import matplotlib.pyplot as plt
   >>> plt.savefig("hospital.jpg")

You should see output something similar to this;

.. image:: ../../images/tutorial_7_4_2.jpg
   :alt: Populations of the H1, H2 and ICU states

Similarly, we can extract the peak patient and ICU populations, via;

.. code-block:: python

    >>> import pandas as pd
    >>> df = pd.read_csv("output/patients.csv.bz2")
    >>> df[ df["H1"] == df["H1"].max() ]
         day       H1      H2     ICU
    127  127  1890553  472957  172004
    >>> df[ df["ICU"] == df["ICU"].max() ]
         day       H1      H2     ICU
    132  132  1779346  445783  180845

This again shows that the time around 130 days since the start of the
outbreak would be most challenging, with a peak of nearly 1.9 million
normal patients, and over 180 thousand ICU patients.

Note that this is a very simplified model and data fitting would be
needed to optimise the various parameters (e.g. the interaction matrix
or the percentages of population who move from, e.g. ``H2`` to ``ICU``).
Also this is missing lots of other movements.

However, we hope that this gives you a good idea of how you can use
demographics, mixing functions / interaction matrices, plus
move functions to conditionally move to different disease stages in
different demographics, to model a wide range of different scenarios.
=====================
Scanning demographics
=====================

Now that we can model super-spreaders, it is important to be able to
scan through their parameters to so that we can investigate their impact.

There are three parameters that we want to scan;

* **beta[2]** for ``lurgy_super``, which controls how infectious super-spreaders
  are in the asymptomatic (or pre-symptomatic) phase.
* **progress[2]** for ``lurgy_super``, which controls how long the super-spreaders
  remain in this phase.
* The percentage of super-spreaders in the normal population.

Scanning disease parameters
---------------------------

We can scan the ``lurgy_super`` parameters using a scan file (also called
a design file), as you have done :doc:`previously <../part02/02_adjustable>`.
The difference now is that instead of changing the disease parameter for
all demographics, we want to change it only for the ``super`` demographic.

To do this, create a design file called ``scan.dat`` and copy in the below;

::

    super:beta[2]  super:progress[2]
        0.5             0.5
        0.5             0.6
        0.5             0.7
        0.5             0.8
        0.5             0.9
        0.5             1.0

        0.6             0.5
        0.6             0.6
        0.6             0.7
        0.6             0.8
        0.6             0.9
        0.6             1.0

        0.7             0.5
        0.7             0.6
        0.7             0.7
        0.7             0.8
        0.7             0.9
        0.7             1.0

        0.8             0.5
        0.8             0.6
        0.8             0.7
        0.8             0.8
        0.8             0.9
        0.8             1.0

        0.9             0.5
        0.9             0.6
        0.9             0.7
        0.9             0.8
        0.9             0.9
        0.9             1.0

        1.0             0.5
        1.0             0.6
        1.0             0.7
        1.0             0.8
        1.0             0.9
        1.0             1.0

This scans through ``beta[2]`` from 0.5 to 1.0, while scanning ``progress[2]``
also from 0.5 to 1.0. We specify that we are only changing these parameters
for the ``super`` demographic by prefixing ``beta[2]`` and ``super[2]`` with
``super:``.

.. note::

   Adding ``demographic_name:`` before a parameter means that you will only
   change that parameter for the demographic called ``demographic_name``.
   Note that the name ``overall`` is reserved, and refers to the overall
   set of networks. Setting a parameter using ``overall:`` will only change
   the parameter in the overall :class:`~metawards.Networks` object, and
   not the demographic sub-networks. To change in all networks you should
   not specify the demographic name.

Scanning demographic percentage
-------------------------------

As of the current version of ``metawards`` it is not possible to use a
scan/design file to change the percentage of individuals in different
demographics. Instead, we need to create several ``demographics.json``
files, and then run ``metawards`` independently for each file.
To do this, create the following four demographics files, which we
will use to scan the super-spreader percentage from 5% to 20%.

``demographics_05.json``

::

    {
        "demographics" : ["home", "super"],
        "work_ratios"  : [ 0.95, 0.05 ],
        "play_ratios"  : [ 0.95, 0.05 ],
        "diseases"     : [ null, "lurgy_super" ]
    }

``demographics_10.json``

::

    {
        "demographics" : ["home", "super"],
        "work_ratios"  : [ 0.90, 0.10 ],
        "play_ratios"  : [ 0.90, 0.10 ],
        "diseases"     : [ null, "lurgy_super" ]
    }

``demographics_15.json``

::

    {
        "demographics" : ["home", "super"],
        "work_ratios"  : [ 0.85, 0.15 ],
        "play_ratios"  : [ 0.85, 0.15 ],
        "diseases"     : [ null, "lurgy_super" ]
    }

``demographics_20.json``

::

    {
        "demographics" : ["home", "super"],
        "work_ratios"  : [ 0.80, 0.20 ],
        "play_ratios"  : [ 0.80, 0.20 ],
        "diseases"     : [ null, "lurgy_super" ]
    }

Running the models
------------------

We now have four demographic files to run, each of which have 36 parameter
combinations to scan, which we would like to repeat 8 times. This will be
1152 individual *model runs*, so we need to use a cluster. Here are example
slurm and PBS job submission scripts for these runs;

::

    #!/bin/bash
    #PBS -l walltime=12:00:00
    #PBS -l select=4:ncpus=64:mem=64GB
    # The above sets 4 nodes with 64 cores each

    source $HOME/envs/metawards/bin/activate

    # change into the directory from which this job was submitted
    cd $PBS_O_WORKDIR

    metawards -d lurgy_home -D demographics_05.json -a ExtraSeedsLondon.dat \
            --extractor extract_none -i scan.dat --repeats 8 \
            --nthreads 16 --force-overwrite-output --mixer mix_evenly \
            --no-spinner --theme simple \
            --output output_05

    metawards -d lurgy_home -D demographics_10.json -a ExtraSeedsLondon.dat \
            --extractor extract_none -i scan.dat --repeats 8 \
            --nthreads 16 --force-overwrite-output --mixer mix_evenly \
            --no-spinner --theme simple \
            --output output_10

    metawards -d lurgy_home -D demographics_15.json -a ExtraSeedsLondon.dat \
            --extractor extract_none -i scan.dat --repeats 8 \
            --nthreads 16 --force-overwrite-output --mixer mix_evenly \
            --no-spinner --theme simple \
            --output output_15

    metawards -d lurgy_home -D demographics_20.json -a ExtraSeedsLondon.dat \
            --extractor extract_none -i scan.dat --repeats 8 \
            --nthreads 16 --force-overwrite-output --mixer mix_evenly \
            --no-spinner --theme simple \
            --output output_20

::

    #!/bin/bash
    #SBATCH --time=12:00:00
    #SBATCH --ntasks=4
    #SBATCH --cpus-per-task=64
    # The above sets 4 nodes with 64 cores each

    source $HOME/envs/metawards/bin/activate

    metawards -d lurgy_home -D demographics_05.json -a ExtraSeedsLondon.dat \
            --extractor extract_none -i scan.dat --repeats 8 \
            --nthreads 16 --force-overwrite-output --mixer mix_evenly \
            --no-spinner --theme simple \
            --output output_05

    metawards -d lurgy_home -D demographics_10.json -a ExtraSeedsLondon.dat \
            --extractor extract_none -i scan.dat --repeats 8 \
            --nthreads 16 --force-overwrite-output --mixer mix_evenly \
            --no-spinner --theme simple \
            --output output_10

    metawards -d lurgy_home -D demographics_15.json -a ExtraSeedsLondon.dat \
            --extractor extract_none -i scan.dat --repeats 8 \
            --nthreads 16 --force-overwrite-output --mixer mix_evenly \
            --no-spinner --theme simple \
            --output output_15

    metawards -d lurgy_home -D demographics_20.json -a ExtraSeedsLondon.dat \
            --extractor extract_none -i scan.dat --repeats 8 \
            --nthreads 16 --force-overwrite-output --mixer mix_evenly \
            --no-spinner --theme simple \
            --output output_20

.. note::

   Notice how the ``--output`` command line option is used to direct the
   output from each different demographic file to a different output
   directory. Also notice how we run the four ``metawards`` calculations
   one after another. There is little need to run them in parallel as
   each calculation is already parallelising its 36 x 8 = 288 *model runs*,
   each of which is running over 16 cores. This means that you would need
   to be running over more than 4608 cores before it is worth parallelising
   the four individual ``metawards`` calculations. In the above job scripts,
   we've only asked for 256 cores. You should adjust the request depending
   on how many cores are available on your cluster.

The jobs will take a while to run, e.g. in my case it took 90 minutes using
256 cores, with each individual *model run* taking about 3 minutes.

Analysing the results
---------------------

The first stage in the analysis is to look at the four overview plots.
These can be generated via;

.. code-block:: bash

   metawards-plot -i output_05/results.csv.bz2
   metawards-plot --animate output_05/overview_*.jpg
   metawards-plot -i output_10/results.csv.bz2
   metawards-plot --animate output_10/overview_*.jpg
   metawards-plot -i output_15/results.csv.bz2
   metawards-plot --animate output_15/overview_*.jpg
   metawards-plot -i output_20/results.csv.bz2
   metawards-plot --animate output_20/overview_*.jpg

You will see plots similar to this (which is for the ``output_20`` 20%
super-spreader demographic)

.. image:: ../../images/tutorial_7_2_1.gif
   :alt: Overview of the scan of the 20% super-spreader demographic

You should see in these that, regardless of the percentage of super-spreaders,
the larger the value of ``super:beta[2]``, the more intense the outbreak,
and the smaller the value of ``super:progress[2]``, the more intense
the outbreak. This makes sense, as you would expect a stronger outbreak
the more infectious the super-spreaders are, and the longer they spend in
the infective state.

By eye you can see that this effect is greatest for the run with the
largest percentage of super-spreaders. We can plot the demographics
for the first run of the ``super:beta[2]==1.0`` and ``super:progress[2]==1.0``
parameters, and animate using the commands;

.. code-block:: bash

   metawards-plot -i output_*/1i0v0i5x001/trajectory.csv.bz2
   metawards-plot --animate output_*/1i0v0i5x001/demographics.jpg -o demographics_1i0v0i5.gif --order filename

This should result in an animation that looks something like this;

.. image:: ../../images/tutorial_7_2_2.gif
   :alt: Overview of the scan of the 20% super-spreader demographic

You can see in this plot that the greater the percentage of super-spreaders,
the faster the outbreak and the more individuals who are infected. Again,
this is what you would expect. As we are modelling the lurgy, we don't have
real data to compare against. For a real outbreak, you would fit the parameters
for the super-spreaders to match observed data.
===========================
Disease pathways and stages
===========================

Up to this point, every individual modelled in ``metawards`` progresses
through the same stages of the same disease pathway. This pathway
is defined in the disease file (e.g. ``lurgy4.json``) and progresses
an individual from being susceptible to infection (``S``) to the
final stage of the disease (normally called ``R`` to represent individuals
who are removed from the outbreak).

Stages of the lurgy
-------------------

For the lurgy, the disease stages are;

* ``S`` - this is the starting point of all individuals

* ``Stage 0`` - this is a holding state. Individuals are moved into this
  state as soon as they are infected. This is used to record internally
  in ``metawards`` to record new infections each day. The ``progress``
  value for this state (``progress[0]``) should be ``1.0``, to show
  that individuals will immediately progress to ``Stage 1`` (the latent, ``E``
  state) the next day. Equally, the ``beta`` value for this state
  (``beta[0]``) should be ``0.0`` as individuals in this state should
  not contribute to the force of infection.

* ``Stage 1`` - this is the latent or ``E`` state. Individuals are held
  in this state with a duration defined by ``progress[1]``. The ``beta[1]``
  value should be ``0.0`` as latent individuals are not infectious and
  should not thus contribute to the force of infection.

* ``Stages 2-4`` - these are the infected or ``I`` states. The lurgy has
  three such states, with ``Stage 2`` (or ``I1``) representing the
  asymptomatic infectious state (low ``beta`` but zero ``too_ill_to_move``),
  then ``Stage 3`` (or ``I2``) representing the initial symptomatic
  state (medium ``beta`` and medium ``too_ill_to_move``), and ``Stage 4``
  (or ``I3``) representing the highly symptomatic state (medium ``beta``
  and high ``too_ill_to_move``).

* ``Stage 5`` - this is the removed or ``R`` state. Individuals in this
  state cannot infect others, and so ``beta`` is ``0.0`` and ``too_ill_to_move``
  is ``1.0`` (no need to model movements of non-infectious or infectable
  individuals). The ``progress`` value is ``0.0`` as, once removed,
  individuals remain in this stage for the remainder of the *model run*.

Modelling a fraction of asymptomatic super-spreaders
----------------------------------------------------

It is really useful to have the ability to model different demographics
following different disease pathways. For example, up to now we have
modelled the lurgy as having an asymptomatic infectious phase which
everyone will move through. However, in reality, only a small proportion
of those infected by the lurgy will become these asymptomatic
"super-spreaders". We would like to investigate how the percentage
of these super-spreaders, and their mobility and infectivity affects the
progression of the outbreak. To do this, we will create two demographics;

1. **home**, which will contain a population who do not progress through
   the asymptomatic infectious stage, and

2. **super**, which will contain a population of super-spreaders who do
   move through the asymptomatic infectious phase.

To start, we need to create disease files for the **home** and **super**
demographics. Do this by creating two files, first ``lurgy_home.json`` that
should contain;

::

  { "name"             : "The Lurgy",
    "version"          : "June 2nd 2020",
    "author(s)"        : "Christopher Woods",
    "contact(s)"       : "christopher.woods@bristol.ac.uk",
    "reference(s)"     : "Completely ficticious disease - no references",
    "beta"             : [0.0, 0.0, 0.5, 0.5, 0.0],
    "progress"         : [1.0, 1.0, 0.5, 0.5, 0.0],
    "too_ill_to_move"  : [0.0, 0.0, 0.5, 0.8, 1.0],
    "contrib_foi"      : [1.0, 1.0, 1.0, 1.0, 0.0]
  }

and ``lurgy_super.json`` that should contain;

::

  { "name"             : "The Lurgy",
    "version"          : "June 2nd 2020",
    "author(s)"        : "Christopher Woods",
    "contact(s)"       : "christopher.woods@bristol.ac.uk",
    "reference(s)"     : "Completely ficticious disease - no references",
    "beta"             : [0.0, 0.0, 0.8, 0.2, 0.1, 0.0],
    "progress"         : [1.0, 1.0, 0.5, 0.5, 0.5, 0.0],
    "too_ill_to_move"  : [0.0, 0.0, 0.0, 0.1, 0.0, 1.0],
    "contrib_foi"      : [1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
  }

In these files we have removed the asymptomatic infectious phase from the
``home`` demographic, and then changed the ``super`` demographic to
experience a highly infectious asymptomatic stage, followed by a
very mild disease of decreasing infectiousness.

Next, we need to create the ``demographics.json`` file that should contain;

::

    {
        "demographics" : ["home", "super"],
        "work_ratios"  : [ 0.9, 0.1 ],
        "play_ratios"  : [ 0.9, 0.1 ],
        "diseases"     : [ null, "lurgy_super" ]
    }

This describes the two demographics, with 90% of individuals in the
``home`` demographic, and 10% in the ``super`` demographic. The new
line here, ``disease_stages``, specifies the file for the disease stages
for each demographic, e.g. ``home`` will follow the default disease, while
``super`` will follow ``lurgy_super``.

.. note::

   ``null`` in a json file means "nothing". In this case, "nothing" means
   that the ``home`` demographic should use the disease parameters from
   the global disease file set by the user.

Run ``metawards`` using the command;

.. code-block:: bash

  metawards -d lurgy_home -D demographics.json --mixer mix_evenly -a ExtraSeedsLondon.dat

.. note::

   Here we set the global disease file, via ``-d lurgy_home`` to ``lurgy_home``.
   This will be used by the ``home`` demographic. In theory we could have
   specified ``lurgy_home`` directly in the ``demographics.json`` file,
   but this would mean we would have to specify it twice, once there and once
   globally. It is better to define things once only, as this leads to fewer
   bugs.

.. note::

   Note that we are using :func:`~metawards.mixers.mix_evenly`.
   This may not be the right choice depending on how we want the
   population dynamics to mix, e.g.
   :func:`~metawards.mixers.mix_evenly_single_population` or
   :func:`~metawards.mixers.mix_evenly_multi_population` may
   be a better choice. :doc:`See here for more information <../part05/04_contacts>`.

You should see that the disease spreads quickly via the super-spreaders, with
all, or almost all becoming infected. You can see this in the demographics
plot, produced via;

.. code-block:: bash

    metawards-plot -i output/trajectory.csv.bz2

You should see something similar to this, which shows that the infection
burns quickly through the super-spreader demographic, moving through that
entire demographic in just a couple of months.

.. image:: ../../images/tutorial_7_1.jpg
   :alt: Demographic trajectory including the super-spreaders

.. note::

   It is counter-intuitive that the super-spreaders are all infected quickly,
   and complete the outbreak long before the general population. This is
   due to the model setup. The super-spreaders all interact with each other
   within their demographic, and so, as their ``beta`` value is high,
   there is a high probability that they will infect each other quickly.
   Their limited population naturally means that it takes less time until
   all members of the this population are infected.

   While not wholly realistic, this does make some practical sense,
   as real super-spreaders will
   go about their normal day during an outbreak as they do not noticeably
   become ill. As they will continue normally, it could make sense that
   they outbreak could burn through that population quickly,
   until the point where there are few remaining. Meanwhile, the larger,
   more general population experiences a slower outbreak.

   One way to counter this effect would be to use an interaction matrix
   to slow down the rate of infection in the super-spreader demographic.
   This would reflect the reality that super-spreaders are more dispersed,
   and so less likely to interact with one another than with members
   of the general population. An interaction matrix of, e.g.
   (1, 1, 1, 0.5) may thus be appropriate, although data fitting would
   be needed to find the exact values.
====================
Modelling a Hospital
====================

Now you've seen how to use and name different disease stages, we can
put all of this together to fully model a hospital.

Disease files
-------------

First, we will create a set of disease files for the different demographics;

First, ``lurgy_home.json``

::

    {
        "stage"            : [ "E1", "E2", "I1", "I2", "I3", "R" ],
        "beta"             : [  0.0,  0.0,  0.2,  0.5,  0.5, 0.0 ],
        "progress"         : [  1.0,  1.0,  0.2,  0.5,  0.5, 0.0 ],
        "too_ill_to_move"  : [  0.0,  0.0,  0.0,  0.5,  0.8, 1.0 ],
        "contrib_foi"      : [  1.0,  1.0,  1.0,  1.0,  1.0, 0.0 ],
        "start_symptom"    : 3
    }

This is the same as ``lurgy4.json`` except we have now named the individual
disease stages. This will be used to model the disease in the general
population and in hospital staff.

Next, ``lurgy_hospital.json``, used to model hospital patients;

::

    {
        "stage"            : [ "H1", "H2" ],
        "beta"             : [  0.2, 0.2  ],
        "progress"         : [  0.2, 0.2  ],
        "too_ill_to_move"  : [  1.0, 1.0  ],
        "contrib_foi"      : [  1.0, 1.0  ],
        "start_symptom"    : 1
    }

Next, ``lurgy_icu.json``, used to model intensive care patients;

::

    {
        "stage"            : [ "ICU", "R" ],
        "beta"             : [  0.2,  0.0 ],
        "progress"         : [  0.2,  0.0 ],
        "too_ill_to_move"  : [  1.0,  1.0 ],
        "contrib_foi"      : [  1.0,  0.0 ],
        "start_symptom"    : 1
    }

And, finally, we need to have an overall disease file that defines the
stages that we want to use for mapping. This just names the mapping
stages, in the order we want to see them reported. This is in
``lurgy_overall.json``;

::

   { "stage" : ["E", "I", "H", "ICU", "R"] }

.. note::

   We don't need to set any disease parameters as these are set by the
   disease files used by the different demographics.

Demographics
------------

Next, we need to update the ``demographics.json`` file to read;

::

    {
        "demographics" : ["home", "staff", "patients", "icu"],
        "work_ratios"  : [ 0.99,   0.01,     0.00,     0.00 ],
        "play_ratios"  : [ 1.00,   0.00,     0.00,     0.00 ],
        "diseases"     : [ "lurgy_home", "lurgy_home",
                           "lurgy_hospital", "lurgy_icu" ]
    }

Here, we've set the ``home`` and ``staff`` demographics to use ``lurgy_home``,
while ``patients`` will use ``lurgy_hospital`` and ``icu`` will use
``lurgy_icu``.

Mixers, movers and extractors
-----------------------------

We can use the same mixer as before, e.g. ``mix_hospital.py`` should be;

.. code-block:: python

    from metawards.mixers import merge_using_matrix

    def mix_shield(network, **kwargs):
        matrix = [ [1.0, 1.0, 0.0, 0.0],
                   [0.0, 0.1, 0.1, 0.1],
                   [0.0, 0.1, 0.0, 0.0],
                   [0.0, 0.1, 0.0, 0.0] ]

        network.demographics.interaction_matrix = matrix

        return [merge_using_matrix]

We can use the same mover as before, except we can now take advantage
of the named disease stages to specify the moves using demographic and
disease name. This is significantly less error-prone than using indexes,
e.g. update ``move_hospital.py`` to read;

.. code-block:: python

    from metawards.movers import go_stage


    def move_hospital(**kwargs):
        # move 20% of I2 home/staff population to H1 patients
        func1 = lambda **kwargs: go_stage(go_from=["home", "staff"],
                                          go_to="patients",
                                          from_stage="I2",
                                          to_stage="H1",
                                          fraction=0.2,
                                          **kwargs)

        # move 10% of H2 patients to ICU1 ICU
        func2 = lambda **kwargs: go_stage(go_from="patients",
                                          go_to="icu",
                                          from_stage="H2",
                                          to_stage="ICU",
                                          fraction=0.1,
                                           **kwargs)

        # move the remainder of H2 patients to home R
        func3 = lambda **kwargs: go_stage(go_from="patients",
                                          go_to="home",
                                          from_stage="H2",
                                          to_stage="R",
                                          fraction=1.0,
                                          **kwargs)

        # move R ICU and H2 patients to home R
        func4 = lambda **kwargs: go_stage(go_from=["patients", "icu"],
                                          go_to="home",
                                          from_stage=["H2", "R"],
                                          to_stage="R",
                                          fraction=1.0,
                                          **kwargs)

        return [func1, func2, func3, func4]

As for the extractor, well we don't really need to use it now as the
data for each stage will be written already into the summary and
individual output files.

Running the job
---------------

Run ``metawards`` using;

.. code-block:: bash

   metawards -D demographics.json -d lurgy_overall --mixer mix_hospital --mover move_hospital -a ExtraSeedsLondon.dat --nsteps 40

.. note::

  We've limited here to running just 40 steps to enable a quick demonstration.
  Feel free to run the full model if you have time.

You should see that statistics for the mapped disease stages for all of the
demographics are written nicely to the screen, e.g.

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 36 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56079510  E: 632  I: 1246  H: 117  ICU: 7  R: 565  IW: 509  POPULATION: 56082077
    ┏━━━━━━━━━━┳━━━━━━━━━━┳━━━━━┳━━━━━━┳━━━━━┳━━━━━┳━━━━━┳━━━━━┳━━━━━━━━━━━━┓
    ┃          ┃    S     ┃  E  ┃  I   ┃  H  ┃ ICU ┃  R  ┃ IW  ┃ POPULATION ┃
    ┡━━━━━━━━━━╇━━━━━━━━━━╇━━━━━╇━━━━━━╇━━━━━╇━━━━━╇━━━━━╇━━━━━╇━━━━━━━━━━━━┩
    │   home   │ 55975294 │ 626 │ 1240 │     │     │ 560 │ 298 │  55977720  │
    │  staff   │  104216  │  6  │  6   │     │     │  5  │  1  │   104233   │
    │ patients │    0     │     │      │ 117 │     │     │ 86  │    117     │
    │   icu    │    0     │     │      │     │  7  │  0  │  7  │     7      │
    ├──────────┼──────────┼─────┼──────┼─────┼─────┼─────┼─────┼────────────┤
    │  total   │ 56079510 │ 632 │ 1246 │ 117 │  7  │ 565 │ 509 │  56082077  │
    └──────────┴──────────┴─────┴──────┴─────┴─────┴─────┴─────┴────────────┘
    Number of infections: 2002

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 37 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56079060  E: 795  I: 1415  H: 148  ICU: 6  R: 653  IW: 622  POPULATION: 56082077
    ┏━━━━━━━━━━┳━━━━━━━━━━┳━━━━━┳━━━━━━┳━━━━━┳━━━━━┳━━━━━┳━━━━━┳━━━━━━━━━━━━┓
    ┃          ┃    S     ┃  E  ┃  I   ┃  H  ┃ ICU ┃  R  ┃ IW  ┃ POPULATION ┃
    ┡━━━━━━━━━━╇━━━━━━━━━━╇━━━━━╇━━━━━━╇━━━━━╇━━━━━╇━━━━━╇━━━━━╇━━━━━━━━━━━━┩
    │   home   │ 55974849 │ 789 │ 1405 │     │     │ 645 │ 381 │  55977688  │
    │  staff   │  104211  │  6  │  10  │     │     │  6  │  5  │   104233   │
    │ patients │    0     │     │      │ 148 │     │     │ 103 │    148     │
    │   icu    │    0     │     │      │     │  6  │  2  │  6  │     8      │
    ├──────────┼──────────┼─────┼──────┼─────┼─────┼─────┼─────┼────────────┤
    │  total   │ 56079060 │ 795 │ 1415 │ 148 │  6  │ 653 │ 622 │  56082077  │
    └──────────┴──────────┴─────┴──────┴─────┴─────┴─────┴─────┴────────────┘
    Number of infections: 2364

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 38 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56078586  E: 924  I: 1625  H: 168  ICU: 9  R: 765  IW: 696  POPULATION: 56082077
    ┏━━━━━━━━━━┳━━━━━━━━━━┳━━━━━┳━━━━━━┳━━━━━┳━━━━━┳━━━━━┳━━━━━┳━━━━━━━━━━━━┓
    ┃          ┃    S     ┃  E  ┃  I   ┃  H  ┃ ICU ┃  R  ┃ IW  ┃ POPULATION ┃
    ┡━━━━━━━━━━╇━━━━━━━━━━╇━━━━━╇━━━━━━╇━━━━━╇━━━━━╇━━━━━╇━━━━━╇━━━━━━━━━━━━┩
    │   home   │ 55974378 │ 916 │ 1614 │     │     │ 758 │ 392 │  55977666  │
    │  staff   │  104208  │  8  │  11  │     │     │  6  │  3  │   104233   │
    │ patients │    0     │     │      │ 168 │     │     │ 115 │    168     │
    │   icu    │    0     │     │      │     │  9  │  1  │  9  │     10     │
    ├──────────┼──────────┼─────┼──────┼─────┼─────┼─────┼─────┼────────────┤
    │  total   │ 56078586 │ 924 │ 1625 │ 168 │  9  │ 765 │ 696 │  56082077  │
    └──────────┴──────────┴─────┴──────┴─────┴─────┴─────┴─────┴────────────┘
    Number of infections: 2726

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 39 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56078059  E: 1001  I: 1914  H: 189  ICU: 11  R: 903  IW: 763  POPULATION: 56082077
    ┏━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━┳━━━━━━┳━━━━━┳━━━━━┳━━━━━┳━━━━━┳━━━━━━━━━━━━┓
    ┃          ┃    S     ┃  E   ┃  I   ┃  H  ┃ ICU ┃  R  ┃ IW  ┃ POPULATION ┃
    ┡━━━━━━━━━━╇━━━━━━━━━━╇━━━━━━╇━━━━━━╇━━━━━╇━━━━━╇━━━━━╇━━━━━╇━━━━━━━━━━━━┩
    │   home   │ 55973854 │ 995  │ 1899 │     │     │ 894 │ 440 │  55977642  │
    │  staff   │  104205  │  6   │  15  │     │     │  7  │  3  │   104233   │
    │ patients │    0     │      │      │ 189 │     │     │ 135 │    189     │
    │   icu    │    0     │      │      │     │ 11  │  2  │ 11  │     13     │
    ├──────────┼──────────┼──────┼──────┼─────┼─────┼─────┼─────┼────────────┤
    │  total   │ 56078059 │ 1001 │ 1914 │ 189 │ 11  │ 903 │ 763 │  56082077  │
    └──────────┴──────────┴──────┴──────┴─────┴─────┴─────┴─────┴────────────┘
    Number of infections: 3115

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 40 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56077432  E: 1154  I: 2225  H: 210  ICU: 12  R: 1044  IW: 874  POPULATION: 56082077
    ┏━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━┳━━━━━━┳━━━━━┳━━━━━┳━━━━━━┳━━━━━┳━━━━━━━━━━━━┓
    ┃          ┃    S     ┃  E   ┃  I   ┃  H  ┃ ICU ┃  R   ┃ IW  ┃ POPULATION ┃
    ┡━━━━━━━━━━╇━━━━━━━━━━╇━━━━━━╇━━━━━━╇━━━━━╇━━━━━╇━━━━━━╇━━━━━╇━━━━━━━━━━━━┩
    │   home   │ 55973231 │ 1147 │ 2208 │     │     │ 1033 │ 525 │  55977619  │
    │  staff   │  104201  │  7   │  17  │     │     │  7   │  4  │   104232   │
    │ patients │    0     │      │      │ 210 │     │      │ 147 │    210     │
    │   icu    │    0     │      │      │     │ 12  │  4   │ 12  │     16     │
    ├──────────┼──────────┼──────┼──────┼─────┼─────┼──────┼─────┼────────────┤
    │  total   │ 56077432 │ 1154 │ 2225 │ 210 │ 12  │ 1044 │ 874 │  56082077  │
    └──────────┴──────────┴──────┴──────┴─────┴─────┴──────┴─────┴────────────┘
    Number of infections: 3601

Equally, you will find the individual ``ICU`` and ``R`` populations for
the ``icu`` demographic is ``output/total_infections_icu.csv.bz2``,
and the individual populations for ``H1`` and ``H2`` in the
``patients`` demographic in ``output/total_infections_patients.csv.bz2``.

In addition, the ``output/results.csv.bz2`` and ``output/trajectory.csv.bz2``
files record the totals in the ``H`` and ``ICU`` stages, as well as the
traditional ``S``, ``E``, ``I`` and ``R``, e.g.

.. code-block:: python

   >>> import pandas as pd
   >>> df = pd.read_csv("output/results.csv.bz2")
   >>> df.head()
      fingerprint  repeat  day        date         S  E  I  H  ICU  R  IW   UV
    0      REPEAT       1    0  2020-06-23  56082077  0  0  0    0  0   0  1.0
    1      REPEAT       1    1  2020-06-24  56082072  5  0  0    0  0   1  1.0
    2      REPEAT       1    2  2020-06-25  56082072  0  5  0    0  0   0  1.0
    3      REPEAT       1    3  2020-06-26  56082071  1  5  0    0  0   1  1.0
    4      REPEAT       1    4  2020-06-27  56082069  3  5  0    0  0   3  1.0
====================
Conditional Pathways
====================

In the last example we saw how each demographic can be given a different
disease pathway. In this section, we will see how to use move functions
to move individuals between pathways. We did this before
:doc:`when we modelled quarantine. <../part06/05_quarantine>`
Now, we will see how, by using different disease stages, we can create
a range of different pathways through which different individuals
can conditionally progress.

Modelling a hospital
--------------------

In this example we will create a model of a hospital. The demographics
involved will be;

* ``patients`` : Patients at the hospital who are infected with the lurgy
* ``staff`` : Hospital staff who care for the patients
* ``home`` : The general population who are neither patients or staff

Members of the ``staff`` and ``home`` demographics will move along the
stages of the lurgy according to ``lurgy4.json``. Hospital ``patients``
will move through two hospital states, which we will refer to as
``H1`` and ``H2``. These are defined in the file ``lurgy_hospital.json``,
which you should create and copy in the below;

::

  { "name"             : "The Lurgy (hospital patients)",
    "version"          : "June 10th 2020",
    "author(s)"        : "Christopher Woods",
    "contact(s)"       : "christopher.woods@bristol.ac.uk",
    "reference(s)"     : "Completely ficticious disease - no references",
    "beta"             : [0.0, 0.0, 0.2, 0.2, 0.0],
    "progress"         : [1.0, 1.0, 0.2, 0.2, 0.0],
    "too_ill_to_move"  : [0.0, 0.0, 1.0, 1.0, 1.0],
    "contrib_foi"      : [0.0, 0.0, 1.0, 1.0, 0.0]
  }

.. note::

   Every disease must include the ``*`` and ``E`` stages, which are
   stages 0 and 1, and the ``R`` stage which is the last stage. This
   means that ``H1`` and ``H2`` are stages 2 and 3 in this file, both
   of which have ``beta == 0.2``, ``progress == 0.2`` and
   ``too_ill_to_move == 1.0``. While the ``H1`` and ``H2`` stages
   are identical now, we will look to change them later in this
   tutorial.

The aim of the model will be that 20% of those suffering in the ``I2``
stage of ``lurgy4.json`` will be moved to the ``H1`` stage of
``lurgy_hospital.json``, and will then progress to ``H2`` and then
``R`` in the hospital.

We will next set that 1% of the worker population are hospital staff, while,
initially, nobody has the lurgy, and so there are no hospital patients.
We can set this using the file ``demographics.json``, which you should
create and copy in the below;

::

    {
        "demographics" : ["home", "staff", "patients"],
        "work_ratios"  : [ 0.99,   0.01,     0.00 ],
        "play_ratios"  : [ 1.00,   0.00,     0.00 ],
        "diseases"     : [ null,   null,   "lurgy_hospital" ]
    }

Next, we need a mixing function that will model the interactions between
hospital staff, patients and the general population. We will use the
following interaction matrix;

.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * -
     - home
     - staff
     - patients
   * - home
     - 1.0
     - 1.0
     - 0.0
   * - staff
     - 0.0
     - 0.1
     - 0.1
   * - patients
     - 0.0
     - 0.1
     - 0.0

This matrix sets that the ``home`` demographic is taking no precautions,
and are fully exposed to each other, and also fully exposed to members
of the ``staff`` demographic.

The ``staff`` demographic are taking lots of precautions, and are not
exposed to the ``home`` demographic, and only lightly exposed to other
``staff`` or members of the ``patient`` demographic (matrix values
are ``0.1``).

Finally, the ``patients`` group are isolated from one another, and so
do not infect one another. At the moment this is moot, as all patients
in this model are already infected. However, you could add a susceptible
patient population who do not have the lurgy, and then use this final
element in the matrix to control the force of infection between patients.

.. note::

   This is not an entirely realistic interaction matrix, but is set up
   to demonstrate how the only route of infection of ``staff`` is from
   ``patients``.

To use this interaction matrix, create a mixer in ``mix_hospital.py``
and copy in the below.

.. code-block:: python

    from metawards.mixers import merge_using_matrix

    def mix_shield(network, **kwargs):
        matrix = [ [1.0, 1.0, 0.0],
                   [0.0, 0.1, 0.1],
                   [0.0, 0.1, 0.0] ]

        network.demographics.interaction_matrix = matrix

        return [merge_using_matrix]

.. note::

   This is essentially an identical mixing function as that used
   in the :doc:`shielding example <../part05/03_custom>`.

.. note::

   Note that we are using :func:`~metawards.mixers.merge_using_matrix`.
   This may not be the right choice depending on how we want the
   population dynamics to mix, e.g.
   :func:`~metawards.mixers.merge_matrix_single_population` or
   :func:`~metawards.mixers.merge_matrix_multi_population` may
   be a better choice. :doc:`See here for more information <../part05/04_contacts>`.

Next, we want to make sure that we have output the populations in the
``H1`` and ``H2`` states. We can do this by creating a custom extractor
that writes the populations of all of the disease stages for just
the ``patients`` demographic to a file called ``patients.csv``. To do this,
create the file ``extract_patients.py`` and copy in the below;

.. code-block:: python

    from metawards.extractors import extract_default


    def output_patients(network, population, workspace, output_dir, **kwargs):
        # Open the file "patients.csv" in the output directory,
        # using the supplied headers for the columns
        FILE = output_dir.open("patients.csv",
                               headers=["day", "*", "E", "H1", "H2", "R"],
                               sep=",")

        # Now get the workspace for the "patients" demographic
        index = network.demographics.get_index("patients")
        subspace = workspace.subspaces[index]

        # The total population at each infection stage is the sum
        # of the work and play infections
        inf_tot = [inf + pinf for inf, pinf in
                   zip(subspace.inf_tot, subspace.pinf_tot)]

        FILE.write(str(population.day) + ",")
        FILE.write(",".join([str(x) for x in inf_tot]) + "\n")


    def extract_patients(**kwargs):
        # return all of the functions from "extract_default"
        # plus our new "output_patients"
        funcs = extract_default(**kwargs)
        funcs.append(output_patients)
        return funcs

Conditionally moving to hospital
--------------------------------

The new part of this example is that we need to add a move function that
will move 20% of individuals who are in the ``I2`` stage to
hospital, in the ``H1`` stage. We can do this by writing a move function
into the file ``move_hospital.py``, which you should create and
copy in the below;

.. code-block:: python

    from metawards.movers import go_stage


    def move_hospital(**kwargs):
        func = lambda **kwargs: go_stage(go_from=["home", "staff"],
                                         go_to="patients",
                                         from_stage=4,
                                         to_stage=2,
                                         fraction=0.2,
                                         **kwargs)

        return [func]

This move function returns :meth:`~metawards.movers.go_stage`. This is
very similar to :meth:`~metawards.movers.go_to`, except you also specify
the ``from_stage`` and ``to_stage``, which are the stage(s) to move from,
and the stage to move to. In this case, we will move 20% of individuals
from the ``I2``
stage from the ``home`` and ``staff`` demographics, which is stage 4
of ``lurgy4.json``. We will move these individuals to stage 2, which is
``H2``, in the ``patients`` demographic.

Now this is set, we can run the model using;

.. code-block:: bash

   metawards -D demographics.json -d lurgy4 --mixer mix_hospital --mover move_hospital --extract extract_hospital -a ExtraSeedsLondon.dat

You should see that the epidemic starts in the ``home`` demographic, with
no infections in ``staff`` until after a number of the more ill ``home``
individuals are moved to hospital as ``patients``. In this model, the
~70,000 staff are overwhelmed with millions of patients. This means that,
despite their care, and the lack of interaction with ``home``,
all of the members of the ``staff`` demographic become infected within
the first ~130 days of the epidemic.

You can plot this using the command;

.. code-block:: bash

   metawards-plot -i output/trajectory.csv.bz2

The resulting plot should look something like this;

.. image:: ../../images/tutorial_7_3.jpg
   :alt: Demographic trajectories for the simple hospital model

You can see that the infections for the ``patients`` lag behind that of the
general population, as individuals move to hospital at a late stage in the
disease, and then remain in hospital for a long time. You can also see
in the ``IW`` plot how the ``staff`` are quickly overwhelmed and are all
infected within ~4 months.

In addition to this plot, you will also see in the output directory
the ``patients.csv.bz2`` file that was written by ``extract_patients``.
We can load this into R, pandas or Excel to find the maximum occupancy
in the hospital for the two different stages.

.. code-block:: python

   >>> import pandas as pd
   >>> df = pd.read_csv("output/patients.csv.bz2")
   >>> df["H1"].max()
   1799688
   >>> df["H2"].max()
   2166228
   >>> df[ df["H1"] == df["H1"].max() ]
        day  *  E       H1       H2        R
   126  126  0  0  1799688  2090894  6647271
   >>> df[ df["H2"] == df["H2"].max() ]
        day  *  E       H1       H2        R
   130  130  0  0  1742752  2166228  8348305

This shows that days 126-130 were particularly difficult, with nearly
4 million patients in hospital, while nearly all staff were infected.
=====================
Modelling vaccination
=====================

Let's now use :class:`~metawards.movers.MoveGenerator` and
:func:`~metawards.movers.go_ward` to model vaccination.

First, we will create a new version of the lurgy that includes
a vaccination stage. To do this in Python open ipython or
jupyter and type;

.. code-block:: python

   >>> from metawards import Disease
   >>> lurgy = Disease("lurgy5")
   >>> lurgy.add("E", beta=0.0, progress=1.0)
   >>> lurgy.add("I1", beta=0.4, progress=0.2)
   >>> lurgy.add("I2", beta=0.5, progress=0.5, too_ill_to_move=0.5)
   >>> lurgy.add("I3", beta=0.5, progress=0.8, too_ill_to_move=0.8)
   >>> lurgy.add("R")
   >>> lurgy.add("V", beta=0.0, progress=0.0, is_infected=False)
   >>> lurgy.to_json("lurgy5.json", indent=2, auto_bzip=False)

or, in R/RStudio type;

.. code-block:: R

   > library(metawards)
   > lurgy <- metawards$Disease("lurgy5")
   > lurgy$add("E", beta=0.0, progress=1.0)
   > lurgy$add("I1", beta=0.4, progress=0.2)
   > lurgy$add("I2", beta=0.5, progress=0.5, too_ill_to_move=0.5)
   > lurgy$add("I3", beta=0.5, progress=0.8, too_ill_to_move=0.8)
   > lurgy$add("R")
   > lurgy$add("V", beta=0.0, progress=0.0, is_infected=FALSE)
   > lurgy$to_json("lurgy5.json", indent=2, auto_bzip=FALSE)

or, simply copy the below into the file ``lurgy5.json``;

::

    {
        "name": "lurgy5",
        "stage":           ["E", "I1", "I2", "I3", "R", "V"],
        "mapping":         ["E", "I",  "I",  "I",  "R", "V"],
        "beta":            [0.0, 0.4,  0.5,  0.5,  0.0, 0.0],
        "progress":        [1.0, 0.2,  0.5,  0.8,  0.0, 0.0],
        "too_ill_to_move": [0.0, 0.0,  0.5,  0.8,  0.0, 0.0],
        "contrib_foi":     [1.0, 1.0,  1.0,  1.0,  1.0, 1.0],
        "is_infected": [true, true, true, true, false, false],
        "start_symptom": 2
    }

This has added an extra stage, called V, which is used to hold
vaccinated individuals. Because vaccinated individuals are
not infected, we have to set ``is_infected`` to ``false`` for
this stage (note how this is set automatically for the R stage).

Also, once vaccinated, individuals will remain in the V stage. Thus
``progress`` for this stage is set to ``0.0``, just like it is
automatically set for the R stage.

We can run a quick test of this model using;

.. code-block::

   metawards -d lurgy5.json -m single -a 5

This uses the single-ward model and infects 5 individuals on day 1 with
the disease described in ``lurgy5.json``.

You should see that the outbreak progresses as before, with many of the
1000-strong population of the ward infected, and, currently, no-one
vaccinated.

Vaccinating using go_ward
-------------------------

We can model vaccination by writing a mover that uses
:class:`~metawards.movers.go_ward` to move up to 50 individuals per
day from S to V (e.g. representing a vaccination capacity of
50 vaccinations per day).

Do this by creating a file called ``move_vaccinate.py`` and copying
in the below;

.. code-block:: python

    from metawards.movers import MoveGenerator, go_ward


    def move_vaccinate(**kwargs):
        capacity = 50  # number of vaccinations per day

        def go_vaccinate(**kwargs):
            gen = MoveGenerator(from_stage="S", to_stage="V", number=capacity)
            go_ward(generator=gen, **kwargs)

        return [go_vaccinate]

Run the model using;

.. code-block:: bash

   metawards -d lurgy5.json -m single -a 5 --mover move_vaccinate.py

You should see that the number of vaccinated individuals climbs by 50 each
day. This reduces the number of susceptibles, meaning that the infection
is dampened before it can get going.

Demand-driven vaccination
-------------------------

This model started vaccination on the first day of the outbreak. We
can instead trigger vaccination based on a threshold of infections,
e.g. because maybe the vaccine has a limited supply. Edit
``move_vaccinate.py`` to add a trigger where vaccination starts
only when the number of infections grows above 100, e.g.

.. code-block:: python

    from metawards.movers import MoveGenerator, go_ward
    from metawards.utils import Console


    def move_vaccinate(**kwargs):
        capacity = 50  # number of vaccinations per day
        trigger = 100  # number of infections to trigger vaccination

        def go_vaccinate(network, population, **kwargs):
            params = network.params

            # Are we vaccinating? If this has not been set then we are not
            is_vaccinating = params.user_params.get("is_vaccinating", False)

            if not is_vaccinating and population.total >= trigger:
                # trigger vaccination
                params.user_params["is_vaccinating"] = True
                is_vaccinating = True
                Console.info("Starting vaccination")

            if is_vaccinating:
                gen = MoveGenerator(from_stage="S", to_stage="V", number=capacity)
                go_ward(generator=gen, network=network,
                        population=population, **kwargs)

        return [go_vaccinate]

Run the model using;

.. code-block:: bash

   metawards -d lurgy5.json -m single -a 5 --mover move_vaccinate.py

You should see now in the output that vaccination starts when the number
of infections rises about 100, which for me happened on day 17.
Vaccination of the remaining susceptibles completed by day 27, with
the outbreak ending with 528 vaccinated individuals and 472 who
contracted the lurgy and recovered.

National vaccination
--------------------

We can extend the model to vaccinate on a per-ward basis if the number
of infections climbs above a threshold in an individual ward. To do this,
modify ``move_vaccinate.py`` to read;

.. code-block:: python

    from metawards import WardID
    from metawards.movers import MoveGenerator, go_ward
    from metawards.utils import Console


    def move_vaccinate(**kwargs):
        capacity = 50  # number of vaccinations per day
        trigger = 20  # number of infections to trigger vaccination
        stop_trigger = 5 # number of infections to stop vaccination

        def go_vaccinate(network, workspace, **kwargs):
            # create ward-local is_vaccinating parameters, that defaults to 0.0
            is_vaccinating = network.nodes.get_custom("is_vaccinating", 0.0)

            vaccinate = []

            I_in_wards = workspace.I_in_wards

            for i in range(1, network.nnodes + 1):
                if not is_vaccinating[i]:
                    if I_in_wards[i] >= trigger:
                        is_vaccinating[i] = 1.0

                if is_vaccinating[i]:
                    if I_in_wards[i] < stop_trigger:
                        is_vaccinating[i] = 0.0
                    else:
                        vaccinate.append(WardID(i, all_commute=True))
                        vaccinate.append(WardID(i))

            nv = int(sum(is_vaccinating))

            if nv > 0:
                print(f"Vaccinating in wards: number = {nv}")
                gen = MoveGenerator(from_ward=vaccinate,
                                    from_stage="S", to_stage="V", number=capacity)
                go_ward(generator=gen, network=network,
                        workspace=workspace, **kwargs)

        return [go_vaccinate]

You can run this model using;

.. code-block:: bash

   metawards -d lurgy5.json -m 2011Data --mover move_vaccinate.py -a ExtraSeedsLondon.dat

In this case, vaccination starts in a ward if the number of infections
grows to 50 or above. Then, there is capacity for 20 vaccinations per day
per ward and per ward-link. If the number of infections in the ward drops
below 5 then vaccination is stopped.

For my run, I see that all wards enter the vaccination program, with
eventually ~44 M vaccinations and ~10 M infections. The progress
plot produced via ``metawards-plot`` is shown below;

.. image:: ../../images/tutorial_9_2.jpg
   :alt: Demographic trajectories for the ward-local vaccinations
==============================
Modelling population movements
==============================

We can use :func:`~metawards.movers.go_ward` to model population movements.
There are many mass population movements that occur through the year,
e.g. people travelling to holiday destinations in the summer, people
travelling home to family at Christmas, and students going to university
in the autumn. Each of these have the potential to spread a local
disease outbreak and seed it in many new locations.

Ward 0, the null ward
---------------------

To perform this move, we first need to gather together all of the students,
into a single "scratch" or "null" ward, from where they can then
be dispersed to the various university wards. We want to do this because
it is much more efficient than trying to sample the double loop of
from each ward to each university ward.

To help with this, ``metawards`` provides a single null (or scratch) ward,
which is not used in any calculation. This is ward 0. It is at index 0
in the internal arrays used by ``metawards``. You can use this ward
as a temporary space to hold individuals during the day. However, you
must ensure that there are no individuals left in this ward at the
end of the day, when :func:`~metawards.extractors.output_core`
is called.

Moving to university
--------------------

We will model the first day of university, when (in this ideal model),
all students will travel from their home ward to live in their
university ward. Students will be modelled as workers in the university
ward, who commute to the university ward each day.

Create a mover called ``move_university.py`` and copy in the below;

.. code-block:: python

    from metawards import WardID
    from metawards.movers import MoveGenerator, MoveRecord, go_ward
    from metawards.utils import Console, ran_int


    def move_university(population, **kwargs):

        def go_null(network, **kwargs):
            record = MoveRecord()
            gen = MoveGenerator(to_ward=WardID(0,0), fraction=0.02)

            Console.print("Moving students to the null ward...")
            go_ward(generator=gen, network=network, record=record, **kwargs)

            nstudents = 0

            for r in record:
                nstudents += r[-1]

            Console.print(f"Number of students: {nstudents}")

            ninfected = network.links.weight[0] - network.links.suscept[0]

            Console.print(f"Number of infected students: {int(ninfected)}")

        def go_university(network, rngs, **kwargs):
            # get the first random number generator from the
            # list of thread-local generators
            rng = rngs[0]

            # how many wards are there?
            nwards = network.nnodes

            # generate a list of 50 random wards that
            # represent 50 university towns
            uni_wards = [ran_int(rng, 1, nwards) for _ in range(0, 50)]

            for uni_ward in uni_wards[0:-1]:
                nstudents = 0

                record = MoveRecord()
                gen = MoveGenerator(from_ward=WardID(0,0),
                                    to_ward=WardID(uni_ward, uni_ward),
                                    fraction=0.075)
                go_ward(generator=gen, network=network,
                        rngs=rngs, record=record, **kwargs)

                for r in record:
                    nstudents += r[-1]

                Console.print(f"{nstudents} went to university ward {uni_ward}")

            # move the remainder of the students to the last university
            record = MoveRecord()
            gen = MoveGenerator(from_ward=WardID(0,0),
                                to_ward=WardID(uni_wards[-1], uni_wards[-1]))

            go_ward(generator=gen, network=network,
                    rngs=rngs, record=record, **kwargs)

            nstudents = 0

            for r in record:
                nstudents += r[-1]

            Console.print(f"{nstudents} went to university ward {uni_wards[-1]}")


        if population.day == 20:
            return [go_null, go_university]
        else:
            return []

There is a lot to discuss with this code, so we will explore each bit
separately.

This code defines a mover (``move_university``), which, on the 20th
day of the outbreak, returns two go functions (``go_null`` and
``go_university``). This is controlled by the code at the end of the
function;

.. code-block:: python

        if population.day == 20:
            return [go_null, go_university]
        else:
            return []

go_null
-------

This is the ``go_null`` function again;

.. code-block:: python

        def go_null(network, **kwargs):
            record = MoveRecord()
            gen = MoveGenerator(to_ward=WardID(0,0), fraction=0.02)

            Console.print("Moving students to the null ward...")
            go_ward(generator=gen, network=network, record=record, **kwargs)

            nstudents = 0

            for r in record:
                nstudents += r[-1]

            Console.print(f"Number of students: {nstudents}")

            ninfected = network.links.weight[0] - network.links.suscept[0]

            Console.print(f"Number of infected students: {int(ninfected)}")

This creates a :class:`~metawards.movers.MoveGenerator` that moves
2% of individuals from every ward and disease stage to become workers
in the null ward.

.. note::

   Workers in the null ward have a :class:`~metawards.WardID` of
   ``WardID(0,0)``, while players have ``WardID(0)``, or, simply ``0``.

We use a :class:`~metawards.movers.MoveRecord` to record the moves,
and count up the total number of individuals moved, which is the
total number of students.

The total number of workers in a ward is held in
:meth:`network.links.weight <metawards.Links.weight>`. The current number
of susceptible workers in a ward is held in
:meth:`network.links.suscept <metawards.Links.suscept>`. The null ward
is at index 0, and so;

.. code-block:: python

   ninfected = network.links.weight[0] - network.links.suscept[0]

gives the number of infected workers in the null ward.

.. note::

   Similarly, :meth:`network.nodes.save_play_suscept <metawards.Nodes.save_play_suscept>`
   is the number of players in a ward, while
   :meth:`network.nodes.play_suscept <metawards.Nodes.play_suscept>` is the
   number of susceptible players in a ward. The difference between these
   is the number of infected individuals in a ward.

go_university
-------------

This is the ``go_university`` function again;

.. code-block:: python

        def go_university(network, rngs, **kwargs):
            # get the first random number generator from the
            # list of thread-local generators
            rng = rngs[0]

            # how many wards are there?
            nwards = network.nnodes

            # generate a list of 50 random wards that
            # represent 50 university towns
            uni_wards = [ran_int(rng, 1, nwards) for _ in range(0, 50)]

            for uni_ward in uni_wards[0:-1]:
                nstudents = 0

                record = MoveRecord()
                gen = MoveGenerator(from_ward=WardID(0,0), to_ward=uni_ward,
                                    fraction=0.075)
                go_ward(generator=gen, network=network,
                        rngs=rngs, record=record, **kwargs)

                for r in record:
                    nstudents += r[-1]

                Console.print(f"{nstudents} went to university ward {uni_ward}")

            # move the remainder of the students to the last university
            record = MoveRecord()
            gen = MoveGenerator(from_ward=WardID(0,0), to_ward=uni_wards[-1])

            go_ward(generator=gen, network=network,
                    rngs=rngs, record=record, **kwargs)

            nstudents = 0

            for r in record:
                nstudents += r[-1]

            Console.print(f"{nstudents} went to university ward {uni_wards[-1]}")

Because this is an illustrative model, we are not using any real world
data. As such, we first need to randomly generate the IDs of 50 wards
that will be designated as the locations of the 50 universities in this
model.

To do this, we use ``rngs``, which is a list of random number generators,
one for each running thread of ``metawards``. As this is a single-threaded
function, we will take the first random number generator, and
assign that to the variable ``rng``, via,

.. code-block:: python

           rng = rngs[0]

The 50 random ward IDs are 50 integers randomly picked between 1 and
the number of wards in the network, using the :func:`~metawards.utils.ran_int`
function;

.. code-block:: python

            # how many wards are there?
            nwards = network.nnodes

            # generate a list of 50 random wards that
            # represent 50 university towns
            uni_wards = [ran_int(rng, 1, nwards) for _ in range(0, 50)]

.. note::

   Note that we haven't checked for duplicated random ward IDs because
   this is just meant to be a simple, illustrative model.

Next, we loop over the first 49 of those random wards, and create
a :class:`~metawards.movers.MoveGenerator` that moves 7.5% of the
students (on average) from being workers in the null ward
(``WardID(0,0)``), to being workers in the university ward
(``WardID(uni_ward, uni_ward)``).

.. code-block:: python

            for uni_ward in uni_wards[0:-1]:
                nstudents = 0

                record = MoveRecord()
                gen = MoveGenerator(from_ward=WardID(0,0),
                                    to_ward=WardID(uni_ward, uni_ward),
                                    fraction=0.075)
                go_ward(generator=gen, network=network,
                        rngs=rngs, record=record, **kwargs)

                for r in record:
                    nstudents += r[-1]

                Console.print(f"{nstudents} went to university ward {uni_ward}")

Finally, we move the remaining students from the null ward to the 50th
university;

.. code-block:: python

            # move the remainder of the students to the last university
            record = MoveRecord()
            gen = MoveGenerator(from_ward=WardID(0,0), to_ward=uni_wards[-1])

            go_ward(generator=gen, network=network,
                    rngs=rngs, record=record, **kwargs)

            nstudents = 0

            for r in record:
                nstudents += r[-1]

            Console.print(f"{nstudents} went to university ward {uni_wards[-1]}")

50 universities should have 2% of the students each, so why did we use 7.5%?
The reason is because ``fraction`` is used with a random binomial distribution
to sample a random number of individuals with that fraction. This
underestimates the number as we draw more students. Thus the number
of students sampled reduces as we loop through the 49 universities,
to the point that far too many are remaining to go into the 50th university.
Through trial and error, we found that 7.5% gave a reasonable number remaining
for the 50th university.

Performing the run
------------------

We will run this model using ``lurgy5.json`` developed in the last chapter.
We can run our mover using;

.. code-block:: bash

   metawards -d lurgy5.json -m 2011Data --mover move_university.py -a ExtraSeedsLondon.dat --nsteps 40

We are using the 2011Data England and Wales model, and have seeded five infections
in London. We've limited the run to 40 steps as we are only interested in the
change in the outbreak caused by the movement of students to universities
which, in this model, occurs on day 20.

You can plot the results using;

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2

The resulting graph shows that the movement to university had no
perceptible impact on the outbreak, e.g.

.. image:: ../../images/tutorial_9_4.jpg
   :alt: Outbreak during which students went to university

This is unsurprising when you look at the number of students going to
university, as printed in the output. For example, I got;

::

    Moving students to the null ward...
    Number of students: 1122060
    Number of infected students: 25
    84257 went to university ward 8253
    77500 went to university ward 7782
    72120 went to university ward 5122
    66461 went to university ward 2257
    61307 went to university ward 7124
    56873 went to university ward 2241
    53190 went to university ward 1776
    49065 went to university ward 4896
    44821 went to university ward 7752
    42100 went to university ward 4987
    38481 went to university ward 6981
    35484 went to university ward 4100
    32644 went to university ward 5628
    30619 went to university ward 6048
    28363 went to university ward 4480
    26001 went to university ward 5377
    24077 went to university ward 1395
    22415 went to university ward 3738
    20795 went to university ward 3999
    19507 went to university ward 5433
    17863 went to university ward 8282
    16316 went to university ward 4203
    15206 went to university ward 6343
    13909 went to university ward 5494
    12823 went to university ward 6732
    12102 went to university ward 1645
    11238 went to university ward 8275
    10269 went to university ward 4975
    9476 went to university ward 8339
    8869 went to university ward 5837
    7943 went to university ward 2057
    7668 went to university ward 4077
    6915 went to university ward 1506
    6469 went to university ward 1405
    5863 went to university ward 2197
    5563 went to university ward 6986
    5137 went to university ward 6863
    4851 went to university ward 6086
    4229 went to university ward 2494
    3999 went to university ward 4986
    3747 went to university ward 7534
    3447 went to university ward 1583
    3108 went to university ward 2696
    2966 went to university ward 2020
    2597 went to university ward 398
    2473 went to university ward 2202
    2329 went to university ward 5530
    2089 went to university ward 8114
    1978 went to university ward 3328
    24538 went to university ward 2293

Of the 1.1M students who moved, only 25 were infected. These 25 were
dispersed over 50 wards, and so had little impact on an outbreak
that was already spreading rapidly due to normal movements.
========================
Creating intricate moves
========================

Up to now you have used the top-level go functions (e.g.
:func:`~metawards.movers.go_isolate`, :func:`~metawards.movers.go_stage`
and :func:`~metawards.movers.go_to`) to perform more complex
moves of individuals between demographics and/or disease stages.

These go functions are built on top of lower-level function
:func:`~metawards.movers.go_ward`. This is a generic go function that
can move individuals between any and all combinations of
demographics, disease stages, wards (for players)
(for workers) and :class:`~metawards.PersonType` (worker or player).

It uses :class:`~metawards.movers.MoveGenerator` to define the move.
This is a generator that generates all of the moves to perform
based on the arguments to its constructor.

Using MoveGenerator
-------------------

You construct a :class:`~metawards.movers.MoveGenerator` object by specifying
the ``from`` and ``to`` states of individuals for the moves to be
generated. You can specify one or more of;

* ``from_demographic`` / ``to_demographic`` : Moving individuals between
  different demographics (and thus different networks)
* ``from_stage`` / ``to_stage`` : Moving individuals between different
  disease stages
* ``from_ward`` / ``to_ward`` : Moving individuals between different
  wards.

You can use the index or the name of the demographic, stage or ward when
specifying the move.

For example;

.. code-block:: python

   from metawards.movers import MoveGenerator

   # Move from disease stage E to I
   gen = MoveGenerator(from_stage="E", to_stage="I")

   # Move from disease stage 0 to 1
   gen = MoveGenerator(from_stage=0, to_stage=1)

   # Move from demographic "red" to "blue"
   gen = MoveGenerator(from_stage="red", to_stage="blue")

   # Move from ward Bristol to London
   gen = MoveGenerator(from_ward="Bristol", to_ward="London")

   # Move from ward with ID 5 to ID 7
   gen = MoveGenerator(from_ward=5, to_ward=7)

You can also specify multiple from and to stages, e.g.

.. code-block:: python

   from metawards.movers import MoveGenerator

   # Move from both E and I to R
   gen = MoveGenerator(from_stage=["E", "I"], to_stage="R")

   # Move from S to E, and also move from I to R
   gen = MoveGenerator(from_stage=["S", "I"], to_stage=["E", "R"])

   # Move from red and blue to green
   gen = MoveGenerator(from_demographic=["red", "blue"], to_demographic="green")

   # Move from Bristol and Oxford to London
   gen = MoveGenerator(from_ward=["Bristol", "Oxford"], to_ward="London")

If a state isn't specified, then all matching individuals will move,
keeping the states the same as much as possible. So a move from
disease stage E to I that does not specify a move of demographic or ward
would move all individuals in every demographic and every ward from
E to I.

You can limit this by specifying multiple states per move, e.g.

.. code-block:: python

   from metawards.movers import MoveGenerator

   # Move everyone in the blue demographic from E to I
   gen = MoveGenerator(from_demographic="blue", from_stage="E", to_stage="I")

   # Move everyone in the I stage from Bristol to Oxford
   gen = MoveGenerator(from_stage="I", from_ward="Bristol", to_ward="Oxford")

   # Move everyone in the red demographic and E stage to the blue
   # demographic and I stage
   gen = MoveGenerator(from_demographic="red", from_stage="E",
                       to_demographic="blue", to_stage="I")

   # Move everyone in Bristol from the red demographic and E stage to the blue
   # demographic and I stage
   gen = MoveGenerator(from_ward="Bristol",
                       from_demographic="red", from_stage="E",
                       to_demographic="blue", to_stage="I")

   # Move everyone from Bristol who is in the red demographic and E stage
   # to Oxford to the blue demographic and I stage
   gen = MoveGenerator(from_ward="Bristol", to_ward="Oxford",
                       from_demographic="red", from_stage="E",
                       to_demographic="blue", to_stage="I")

   # Move everyone who is in the red demographic to Oxford from wherever
   # they are now to Oxford
   gen = MoveGenerator(from_demographic="red", to_ward="Oxford")

.. note::

   Note how you can specify the ``from`` without a ``to``. This would mean
   move everyone who matches the ``from``. Note also how you can specify
   the ``to`` without the ``from``. In this case, it moves everyone
   to the ``to``.

WardID and commuters
--------------------

``from_ward`` and ``to_ward`` are a little more complex because individuals
are either workers or players.

Players are modelled in a single ward, and so can be identified just by
the ward ID or name. Thus ``from_ward="Bristol"`` means only the players
who reside in Bristol.

Workers are modelled in a ward-link (or ward-connection). This is a link
from the home ward of the worker to the commute ward where they work each
day. We need to specify both the home and commute ward to identify a worker.
To do this, the :class:`metawards.WardID` class is used, e.g.

.. code-block:: python

   from metawards import WardID
   from metawards.movers import MoveGenerator

   # Move all of the workers who live in Bristol and commute to Oxford
   # to become players who live in London
   gen = MoveGenerator(from_ward=WardID("Bristol", "Oxford"), to_ward="London")

Often you want to identify all workers who reside in a ward, not just those
that commute between two wards. To do this, you need to set
``all_commute`` to true, e.g.

.. code-block:: python

   from metawards import WardID
   from metawards.movers import MoveGenerator

   # Move all of the workers who live in Bristol
   # to become players who live in London
   gen = MoveGenerator(from_ward=WardID("Bristol", all_commute=True),
                       to_ward="London")

This enables you to specify moves with a lot of detail, e.g.

.. code-block:: python

   from metawards import Network, Ward, Network
   from metawards.movers import MoveGenerator

   # create the Bristol, Oxford and London wards
   bristol = Ward("Bristol")
   london = Ward("London")
   oxford = Ward("Oxford")

   # Add the connections for commuters
   bristol.add_workers(0, destination=london)
   bristol.add_workers(0, destination=oxford)
   london.add_workers(0, destination=bristol)
   london.add_workers(0, destination=oxford)

   # build a network from these wards
   network = Network.from_wards(bristol+london+oxford)

   # Move everyone who lives in Bristol from the red to blue demographic
   gen = MoveGenerator(from_wards=["bristol",
                                    WardID("bristol", all_commute=True)],
                       from_demographic="red", to_demographic="blue")

   # Move only workers who live in Bristol from the red to blue demographic
   gen = MoveGenerator(from_wards=WardID("bristol", all_commute=True),
                       from_demographic="red", to_demographic="blue")

   # Move only workers who live in Bristol and commute to London
   # from the red to blue demographic
   gen = MoveGenerator(from_wards=WardID("bristol", "london"),
                       from_demographic="red", to_demographic="blue")

   # Move only workers who commute to London from red to blue
   gen = MoveGenerator(from_wards=[WardID("bristol", "london"),
                                   WardID("oxford", "london")],
                       from_demographic="red", to_demographic="blue")

   # Move all workers who commute to Oxford to become players in Oxford
   gen = MoveGenerator(from_wards=[WardID("bristol", "oxford"),
                                   WardID("london", "oxford")],
                       to_wards="Oxford")

   # Move all players in Bristol to become workers who commute to London
   gen = MoveGenerator(from_wards="Bristol",
                       to_wards=WardID("Bristol", "London"))

Moving a number or fraction
---------------------------

You can also control how many individuals will be moved using either
the ``number`` or ``fraction`` arguments. ``number`` specifies the
maximum number of individuals in an individual ward or ward-link who can move.
``fraction`` specifies the fraction (percentage) of individuals from an
individual ward or ward-link who can move. ``fraction`` should be
between 0 and 1. If it is not 1, then the fraction of individuals
are sampled according to the random binomial distribution. For example;

.. code-block:: python

   from metawards.movers import MoveGenerator

   # move 50% of individuals from the red to blue demographics
   gen = MoveGenerator(from_demographic="red", to_demographic="blue",
                       fraction=0.5)

   # move up to 10 individuals from each ward or ward-link from the
   # S stage to the E stage. If there are more than or equal to
   # 10 matching individuals, then all 10 will be moved. Else, only
   # the number who match will be moved.
   gen = MoveGenerator(from_stage="S", to_stage="E", number=10)

   # move 25% of the maximum number of 10 players in Bristol to play in Oxford.
   # In this case, 25% of the up-to 10 individuals in Bristol will be sampled
   # using the binomial distribution.
   gen = MoveGenerator(from_ward="Bristol", to_ward="Oxford",
                       number=10, fraction=0.25)

Using go_ward
-------------

:class:`~metawards.movers.MoveGenerator` is used to generate the moves that
are made by :func:`~metawards.movers.go_ward`. This is a go_function that
you can use in a move function. For example, let's create now a
custom go_function that will use :class:`~metawards.movers.MoveGenerator`
and :func:`~metawards.movers.go_ward` to move individuals from the
R stage back to S. This would imply that as soon as they have recovered,
they are not immune and can be infected again.

To do this, create a move function in a file called ``move_cycle.py``
and copy in the below;

.. code-block:: python

   from metawards.movers import MoveGenerator, go_ward


   def move_cycle(**kwargs):

       # Create the go-function
       def go_cycle(**kwargs):
           gen = MoveGenerator(from_stage="R", to_stage="S")
           go_ward(generator=gen, **kwargs)

       # Return this function to be called
       return [go_cycle]

.. note::

   We've put go_cycle inside move_cycle as this is cleaner than
   having it as a function defined in global scope. This style will
   also be used in later pages in this tutorial as it will enable
   information to be passed between multiple go functions.

You can run this model using;

.. code-block:: bash

   metawards --mover move_cycle -m single -d lurgy -a 5

.. note::

   Here we are using the original lurgy disease model and the single
   ward network for speed. We have seeded the outbreak with 5 infections.

You should see that the outbreak cycles forever (cutting off at the
automatic 2-year - 720 day mark). The plot of ``results.csv.bz2``
(e.g. produced using ``metawards-plot``) shows the outbreak becoming
random noise once R individuals are moved back into S, e.g.

.. image:: ../../images/tutorial_9_1.jpg
   :alt: Demographic trajectories for a cyclic model
=======================
Holidays and quarantine
=======================

Another use for :func:`~metawards.movers.go_ward` is to provide a good
model of the impact of individuals taking holidays to infected
destinations.

Background force of infection
-----------------------------

We can model this by making use of the background force of infection
(``bg_foi``) parameter. This can be set on a per-ward basis via
:meth:`Ward.bg_foi <metawards.Ward.bg_foi>`, or globally for a
network via the :meth:`Parameters.bg_foi <metawards.Parameters.bg_foi>`
parameters.

The background FOI is the starting value for the FOI calculation in each
ward. If this is positive, then this implies that the ward has a background
outbreak that will drive infections in addition to any impact on the
FOI from infected individuals. If this is negative, then this implies
that the ward has a background mitigation that reduces the FOI by a
fixed amount regardless of the number of infected individuals (subject
to the FOI not dropping below zero).

It is a constant that is added to the FOI for each ward. The global
:meth:`Parameters.bg_foi <metawards.Parameters.bg_foi>` is added first
to all wards, and then the per-ward
:meth:`Ward.bg_foi <metawards.Ward.bg_foi>` is added for each ward.

Holiday to a different demographic
----------------------------------

We will model the holiday destination as being a different demographic,
where the global background FOI is set to a high positive value. This means
that individuals who are moved to this demographic will be exposed to
infection. When they then return from the holiday demographic back to
their home demographic, they will carry the infection with them and
then drive the epidemic back at home.

We will start by creating two demographic, ``home`` and ``holiday``.
You can do this using Python, R, or by copying the file from below,
e.g. in Python;

.. code-block:: python

    >>> from metawards import Demographic, VariableSet
    >>> home = Demographic("home")
    >>> home.work_ratio = 1.0
    >>> home.play_ratio = 1.0
    >>> holiday = Demographic("holiday")
    >>> holiday.work_ratio = 0.0
    >>> holiday.play_ratio = 0.0
    >>> adjustment = VariableSet()
    >>> adjustment["scale_uv"] = 0.0
    >>> adjustment["dyn_dist_cutoff"] = 0.0
    >>> adjustment["bg_foi"] = 0.05
    >>> holiday.adjustment = adjustment
    >>> demographics = home + holiday
    >>> demographics.to_json("demographics.json", indent=2, auto_bzip=False)

or in R;

.. code-block:: R

    > library(metawards)
    > home <- metawards$Demographic("home")
    > home$work_ratio <- 1.0
    > home$play_ratio <- 1.0
    > holiday <- metawards$Demographic("holiday")
    > holiday$work_ratio <- 0.0
    > holiday$play_ratio <- 0.0
    > adjustment <- metawards$VariableSet()
    > adjustment$set_value("scale_uv", 0.0)
    > adjustment$set_value("dyn_dist_cutoff", 0.0)
    > adjustment$set_value("bg_foi", 0.05)
    > holiday$adjustment <- adjustment
    > demographics <- metawards$Demographics()
    > demographics$add(home)
    > demographics$add(holiday)
    > demographics$to_json("demographics.json", indent=2, auto_bzip=False)

or copy this text to the file ``demographics.json``.

::

  {
    "demographics": ["home", "holiday"],
    "work_ratios":  [  1.0,     0.0   ],
    "play_ratios":  [  1.0,     0.0   ],
    "adjustments": [
        null,
        {
            "variables": {
                "scale_uv": 0,
                "dyn_dist_cutoff": 0,
                "bg_foi": 0.05
            }
        }
    ]
  }

This create ``home``, which will start with all of the population
of the network. ``holiday`` is created as an empty network.

``holiday`` has its parameters adjusted according to;

* :meth:`Parameters.scale_uv <metawards.Parameters.scale_uv>` is set to
  zero. This stops members of the ``holiday`` demographic from infecting
  each other. They will only experince the force of infection as set
  by the background FOI. This is realistic, as holidaymakers are unlikely
  to travel together in large groups, and so are unlikely to infect
  one another. This of course could be changed by setting scale_uv
  to a value above zero.
* :meth:`Parameters.dyn_dist_cutoff <metawards.Parameters.dyn_dist_cutoff>` is
  also set to zero. This stops holidaymakers from travelling outside
  their ward. This is because this is because they are not really in
  their home ward - it is used just as a marker so we can remember which
  ward they came from so that they can return their at the end of their
  holiday.
* :meth:`Parameters.bg_foi <metawards.Parameters.bg_foi>` is set to 0.05.
  This is a reasonable starting value for the FOI, which will drive
  a small outbreak infecting ~0.2% of susceptibles per day.

go_holiday and go_home
----------------------

Next we will create ``move_holiday``. Open a file called ``move_holiday.py``
and copy in the below;

.. code-block:: python

    from metawards.movers import MoveGenerator, go_ward


    def move_holiday(population, **kwargs):

        def go_holiday(**kwargs):
            gen = MoveGenerator(from_demographic="home",
                                to_demographic="holiday",
                                fraction=0.005)
            go_ward(generator=gen, **kwargs)

        def go_home(**kwargs):
            gen = MoveGenerator(from_demographic="holiday",
                                to_demographic="home")
            go_ward(generator=gen, **kwargs)

        if population.day == 1:
            return [go_holiday]
        elif population.day == 14:
            return [go_home]
        else:
            return []

This defines ``move_holiday``. This move function returns ``go_holiday``
on day 1 of the model, which will send individuals on holiday. It then
returns ``go_home`` on day 14 of the outbreak, which will bring
the holidaymakers back home.

Both ``go_holiday`` and ``go_home`` are simple go functions that
move individuals between the ``home`` and ``holiday`` demographics.
``go_holiday`` moves 0.5% of the population to ``holiday``, while
``go_home`` brings them all back home.

We can run this model using;

.. code-block:: bash

   metawards -d lurgy5.json -D demographics.json -m 2011Data --mover move_holiday.py --nsteps 28

.. note::

   We are still using lurgy5.json from the last chapter, and the 2011 England
   and Wales model. We have set ``nsteps`` to 28 as we are only interested
   in how the epidemic spreads in the first two weeks after the
   holidaymakers return.

.. note::

   Note also that we haven't seeded the infection. Holidaymakers picked
   up the infection based only on the ``bg_foi`` of the holiday demographic.

You will see that the ~280k individuals went on holiday, with about
800 per day becoming infected;

::

    ─────────────────────────────────────────────── Day 0 ────────────────────────────────────────────────
    S: 56082077  E: 0  I: 0  V: 0  R: 0  IW: 0  POPULATION: 56082077
    ┏━━━━━━━━━┳━━━━━━━━━━┳━━━┳━━━┳━━━┳━━━┳━━━━┳━━━━━━━━━━━━┓
    ┃         ┃    S     ┃ E ┃ I ┃ V ┃ R ┃ IW ┃ POPULATION ┃
    ┡━━━━━━━━━╇━━━━━━━━━━╇━━━╇━━━╇━━━╇━━━╇━━━━╇━━━━━━━━━━━━┩
    │  home   │ 56082077 │ 0 │ 0 │ 0 │ 0 │ 0  │  56082077  │
    │ holiday │    0     │ 0 │ 0 │ 0 │ 0 │ 0  │     0      │
    ├─────────┼──────────┼───┼───┼───┼───┼────┼────────────┤
    │  total  │ 56082077 │ 0 │ 0 │ 0 │ 0 │ 0  │  56082077  │
    └─────────┴──────────┴───┴───┴───┴───┴────┴────────────┘


    ─────────────────────────────────────────────── Day 1 ────────────────────────────────────────────────
    S: 56081266  E: 811  I: 0  V: 0  R: 0  IW: 742  POPULATION: 56082077
    ┏━━━━━━━━━┳━━━━━━━━━━┳━━━━━┳━━━┳━━━┳━━━┳━━━━━┳━━━━━━━━━━━━┓
    ┃         ┃    S     ┃  E  ┃ I ┃ V ┃ R ┃ IW  ┃ POPULATION ┃
    ┡━━━━━━━━━╇━━━━━━━━━━╇━━━━━╇━━━╇━━━╇━━━╇━━━━━╇━━━━━━━━━━━━┩
    │  home   │ 55801092 │  0  │ 0 │ 0 │ 0 │  0  │  55801092  │
    │ holiday │  280174  │ 811 │ 0 │ 0 │ 0 │ 742 │   280985   │
    ├─────────┼──────────┼─────┼───┼───┼───┼─────┼────────────┤
    │  total  │ 56081266 │ 811 │ 0 │ 0 │ 0 │ 742 │  56082077  │
    └─────────┴──────────┴─────┴───┴───┴───┴─────┴────────────┘

    Number of infections: 811

    ─────────────────────────────────────────────── Day 2 ────────────────────────────────────────────────
    S: 56080493  E: 773  I: 811  V: 0  R: 0  IW: 721  POPULATION: 56082077
    ┏━━━━━━━━━┳━━━━━━━━━━┳━━━━━┳━━━━━┳━━━┳━━━┳━━━━━┳━━━━━━━━━━━━┓
    ┃         ┃    S     ┃  E  ┃  I  ┃ V ┃ R ┃ IW  ┃ POPULATION ┃
    ┡━━━━━━━━━╇━━━━━━━━━━╇━━━━━╇━━━━━╇━━━╇━━━╇━━━━━╇━━━━━━━━━━━━┩
    │  home   │ 55801092 │  0  │  0  │ 0 │ 0 │  0  │  55801092  │
    │ holiday │  279401  │ 773 │ 811 │ 0 │ 0 │ 721 │   280985   │
    ├─────────┼──────────┼─────┼─────┼───┼───┼─────┼────────────┤
    │  total  │ 56080493 │ 773 │ 811 │ 0 │ 0 │ 721 │  56082077  │
    └─────────┴──────────┴─────┴─────┴───┴───┴─────┴────────────┘

    Number of infections: 1584

By the time they return, there are over 8000 infected holidaymakers,
who spread the virus throughout the country. There is a rapid
increase in the number of infected wards as the infected holidaymakers
make their random (player) or fixed (worker) movements between wards;

::

    ─────────────────────────────────────────────── Day 13 ───────────────────────────────────────────────
    S: 56072036  E: 718  I: 5657  V: 0  R: 3666  IW: 662  POPULATION: 56082077
    ┏━━━━━━━━━┳━━━━━━━━━━┳━━━━━┳━━━━━━┳━━━┳━━━━━━┳━━━━━┳━━━━━━━━━━━━┓
    ┃         ┃    S     ┃  E  ┃  I   ┃ V ┃  R   ┃ IW  ┃ POPULATION ┃
    ┡━━━━━━━━━╇━━━━━━━━━━╇━━━━━╇━━━━━━╇━━━╇━━━━━━╇━━━━━╇━━━━━━━━━━━━┩
    │  home   │ 55801092 │  0  │  0   │ 0 │  0   │  0  │  55801092  │
    │ holiday │  270944  │ 718 │ 5657 │ 0 │ 3666 │ 662 │   280985   │
    ├─────────┼──────────┼─────┼──────┼───┼──────┼─────┼────────────┤
    │  total  │ 56072036 │ 718 │ 5657 │ 0 │ 3666 │ 662 │  56082077  │
    └─────────┴──────────┴─────┴──────┴───┴──────┴─────┴────────────┘

    Number of infections: 6375

    ─────────────────────────────────────────────── Day 14 ───────────────────────────────────────────────
    S: 56069022  E: 3014  I: 5727  V: 0  R: 4314  IW: 2396  POPULATION: 56082077
    ┏━━━━━━━━━┳━━━━━━━━━━┳━━━━━━┳━━━━━━┳━━━┳━━━━━━┳━━━━━━┳━━━━━━━━━━━━┓
    ┃         ┃    S     ┃  E   ┃  I   ┃ V ┃  R   ┃  IW  ┃ POPULATION ┃
    ┡━━━━━━━━━╇━━━━━━━━━━╇━━━━━━╇━━━━━━╇━━━╇━━━━━━╇━━━━━━╇━━━━━━━━━━━━┩
    │  home   │ 56069022 │ 3014 │ 5727 │ 0 │ 4314 │ 2396 │  56082077  │
    │ holiday │    0     │  0   │  0   │ 0 │  0   │  0   │     0      │
    ├─────────┼──────────┼──────┼──────┼───┼──────┼──────┼────────────┤
    │  total  │ 56069022 │ 3014 │ 5727 │ 0 │ 4314 │ 2396 │  56082077  │
    └─────────┴──────────┴──────┴──────┴───┴──────┴──────┴────────────┘

    Number of infections: 8741

    ─────────────────────────────────────────────── Day 15 ───────────────────────────────────────────────
    S: 56065953  E: 3069  I: 8059  V: 0  R: 4996  IW: 2389  POPULATION: 56082077
    ┏━━━━━━━━━┳━━━━━━━━━━┳━━━━━━┳━━━━━━┳━━━┳━━━━━━┳━━━━━━┳━━━━━━━━━━━━┓
    ┃         ┃    S     ┃  E   ┃  I   ┃ V ┃  R   ┃  IW  ┃ POPULATION ┃
    ┡━━━━━━━━━╇━━━━━━━━━━╇━━━━━━╇━━━━━━╇━━━╇━━━━━━╇━━━━━━╇━━━━━━━━━━━━┩
    │  home   │ 56065953 │ 3069 │ 8059 │ 0 │ 4996 │ 2389 │  56082077  │
    │ holiday │    0     │  0   │  0   │ 0 │  0   │  0   │     0      │
    ├─────────┼──────────┼──────┼──────┼───┼──────┼──────┼────────────┤
    │  total  │ 56065953 │ 3069 │ 8059 │ 0 │ 4996 │ 2389 │  56082077  │
    └─────────┴──────────┴──────┴──────┴───┴──────┴──────┴────────────┘

    Number of infections: 11128

The outbreak grows rapidly, until by day 28 there have been over 300k
infections;

::

    ─────────────────────────────────────────────── Day 28 ───────────────────────────────────────────────
    S: 55768352  E: 65734  I: 194526  V: 0  R: 53465  IW: 8454  POPULATION: 56082077
    ┏━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━┳━━━━━━━━┳━━━┳━━━━━━━┳━━━━━━┳━━━━━━━━━━━━┓
    ┃         ┃    S     ┃   E   ┃   I    ┃ V ┃   R   ┃  IW  ┃ POPULATION ┃
    ┡━━━━━━━━━╇━━━━━━━━━━╇━━━━━━━╇━━━━━━━━╇━━━╇━━━━━━━╇━━━━━━╇━━━━━━━━━━━━┩
    │  home   │ 55768352 │ 65734 │ 194526 │ 0 │ 53465 │ 8454 │  56082077  │
    │ holiday │    0     │   0   │   0    │ 0 │   0   │  0   │     0      │
    ├─────────┼──────────┼───────┼────────┼───┼───────┼──────┼────────────┤
    │  total  │ 55768352 │ 65734 │ 194526 │ 0 │ 53465 │ 8454 │  56082077  │
    └─────────┴──────────┴───────┴────────┴───┴───────┴──────┴────────────┘

    Number of infections: 260260

Modelling quarantine
--------------------

One method of preventing holidaymakers from spreading the infection when
they return home is to require them all to enter quarantine. We can
model this by creating a third, ``quarantine`` demographic. You can
do this in python by typing;

.. code-block:: python

    >>> from metawards import Demographic, VariableSet
    >>> home = Demographic("home")
    >>> home.work_ratio = 1.0
    >>> home.play_ratio = 1.0
    >>> holiday = Demographic("holiday")
    >>> holiday.work_ratio = 0.0
    >>> holiday.play_ratio = 0.0
    >>> adjustment = VariableSet()
    >>> adjustment["scale_uv"] = 0.0
    >>> adjustment["dyn_dist_cutoff"] = 0.0
    >>> adjustment["bg_foi"] = 0.05
    >>> holiday.adjustment = adjustment
    >>> quarantine = Demographic("quarantine")
    >>> quarantine.work_ratio = 0.0
    >>> quarantine.play_ratio = 0.0
    >>> adjustment = VariableSet()
    >>> adjustment["scale_uv"] = 0.0
    >>> adjustment["dyn_dist_cutoff"] = 0.0
    >>> quarantine.adjustment = adjustment
    >>> demographics = home + holiday + quarantine
    >>> demographics.to_json("demographics.json", indent=2, auto_bzip=False)

or in R;

.. code-block:: R

    > library(metawards)
    > home <- metawards$Demographic("home")
    > home$work_ratio <- 1.0
    > home$play_ratio <- 1.0
    > holiday <- metawards$Demographic("holiday")
    > holiday$work_ratio <- 0.0
    > holiday$play_ratio <- 0.0
    > adjustment <- metawards$VariableSet()
    > adjustment$set_value("scale_uv", 0.0)
    > adjustment$set_value("dyn_dist_cutoff", 0.0)
    > adjustment$set_value("bg_foi", 0.05)
    > holiday$adjustment <- adjustment
    > quarantine <- metawards$Demographic("quarantine")
    > quarantine$work_ratio <- 0.0
    > quarantine$play_ratio <- 0.0
    > adjustment <- metawards$VariableSet()
    > adjustment$set_value("scale_uv", 0.0)
    > adjustment$set_value("dyn_dist_cutoff", 0.0)
    > quarantine$adjustment <- adjustment
    > demographics <- metawards$Demographics()
    > demographics$add(home)
    > demographics$add(holiday)
    > demographics$add(quarantine)
    > demographics$to_json("demographics.json", indent=2, auto_bzip=False)

or copy this text to the file ``demographics.json``.

::

  {
    "demographics": ["home", "holiday", "quarantine"],
    "work_ratios": [  1.0,     0.0,      0.0 ],
    "play_ratios": [  1.0,     0.0,      0.0 ],
    "adjustments": [
        null,
        {
        "variables": {
            "scale_uv": 0,
            "dyn_dist_cutoff": 0,
            "bg_foi": 0.05
            }
        },
        {
        "variables": {
            "scale_uv": 0,
            "dyn_dist_cutoff": 0
            }
        }
    ]
  }

The ``qurantine`` demographic has both
:meth:`Parameters.scale_uv <metawards.Parameters.scale_uv>` and
:meth:`Parameters.dyn_dist_cutoff <metawards.Parameters.dyn_dist_cutoff>`
as zero, just like the holiday demographic, because holidaymakers
should be self-isolating and thus not infecting one another. The
:meth:`Parameters.bg_foi <metawards.Parameters.bg_foi>` parameter is
not set, meaning that it keeps its default value of zero, meaning that
there is no background driver for more infections.

go_quarantine
-------------

We now need to write a ``go_quarantine`` function and add it to
``move_holiday``. Edit ``move_holiday.py`` and copy in the below;

.. code-block:: python

    from metawards.movers import MoveGenerator, go_ward


    def move_holiday(population, **kwargs):

        def go_holiday(**kwargs):
            gen = MoveGenerator(from_demographic="home",
                                to_demographic="holiday",
                                fraction=0.005)
            go_ward(generator=gen, **kwargs)

        def go_quarantine(**kwargs):
            gen = MoveGenerator(from_demographic="holiday",
                                to_demographic="quarantine")
            go_ward(generator=gen, **kwargs)

        def go_home(**kwargs):
            gen = MoveGenerator(from_demographic="quarantine",
                                to_demographic="home")
            go_ward(generator=gen, **kwargs)

        holiday_length = 14
        quarantine_length = 14

        if population.day == 1:
            return [go_holiday]
        elif population.day == holiday_length:
            return [go_quarantine]
        elif population.day == holiday_length+quarantine_length:
            return [go_home]
        else:
            return []

Here we have changed ``go_holiday`` so that individuals go to
``quarantine`` instead of ``home``. We have then added ``go_quarantine``
that moves all individuals from ``quarantine`` to ``home``.
We've then set ``move_holiday`` to return ``go_quarantine`` at the
end of their 14-day holiday, and then return ``go_home`` at the end
of the 14-day quarantine.

We can run this using;

.. code-block:: bash

   metawards -d lurgy5.json -D demographics.json -m 2011Data --mover move_holiday.py --nsteps 48

and then plot using

.. code-block:: bash

   metawards-plot -i output/results.csv.bz2

What I see is that the 14 days of quarantine is not quite long enough for the
lurgy. There are still ~400 infections, which drive a further 300 infections
when the holidaymakers leave quarantine. This can be see in the infection
graphs, e.g.

.. image:: ../../images/tutorial_9_5.jpg
   :alt: Outbreak amongst holidaymakers, spreading after quarantine

What next?
----------

You could extend the above model in many ways;

* You could model different lengths of quarantine by setting quarantine_length
  via a custom user parameter, and then scanning it to see the impact
  on the outbreak.
* You could model different holiday destinations with different background
  FOIs, to see how you could classify destinations according to the risk
  they pose to home.
* You could model different levels of compliance with quarantine by having
  a fraction of holidaymakers go straight home (set ``fraction`` to a value
  less than one in ``go_quarantine`` and have the remainder ``go_home``).
* You could model quarantine fatigue, similarly to how you modelled
  self-isolation fatigue in :doc:`part 6 <../part06/04_compliance>`.
* You could model the impact of holidaymakers during different stages
  in an outbreak at home, by seeding the outbreak at home, and then
  having the holiday start later in the year.
* You could model the impact of different lockdown measures at home,
  and see how they are disrupted by holidaymakers returning from holiday.

.. note::

   Congratulations on reaching the end of the tutorial. Hopefully this
   has given you a good insight into what you could do with ``metawards``.
   If you have any questions then please
   `post an issue to our GitHub repository <https://github.com/metawards/MetaWards/issues>`_.
   Please feel free to post issues to request more tutorials, or to
   ask for more documentation.

.. note::

   If you are moving on to developing ``metawards``, or writing more
   complex code on top of ``metawards``, then you should read the
   developer documentation for the :mod:`metawards` module, which can be
   :doc:`found here <../../api/index>`.

===============
Fading immunity
===============

Immunity to a disease is not always permanent. For some diseases, immunity
fades over time, meaning that an individual can be infected multiple times.

We can model this by using :func:`~metawards.movers.go_ward` to move
individuals from the R stage back to the S stage after a period of time.

Slowing progress
----------------

The :meth:`Disease.progress <metawards.Disease.progress>` parameter controls
the rate at which an individual will move through a disease stage.
Advancement through the disease stages is controlled by
:func:`~metawards.iterators.advance_recovery`. This simply samples
the fraction, :meth:`Disease.progress <metawards.Disease.progress>`,
from the number of individuals who are at that disease stage via
a random binomial distribution, and advances them to the next stage.

The key code that does this is summarised here;

.. code-block:: python

    # loop from the penultimate disease stage back to the first stage
    for i in range(N_INF_CLASSES-2, -1, -1):
        ...

        # get the progress parameter for this disease stage
        disease_progress = params.disease_params.progress[i]

        # loop over all ward-links for this disease stage
        for j in range(1, nlinks_plus_one):
            # get the number of workers in this link at this stage
            inf_ij = infections_i[j]

            if inf_ij > 0:
                # sample l workers from this stage based on disease_progress
                l = _ran_binomial(rng, disease_progress, inf_ij)

                if l > 0:
                    # move l workers from this stage to the next stage
                    infections_i_plus_one[j] += l
                    infections_i[j] -= l

        # loop over all nodes / wards for this disease stage
        for j in range(1, nnodes_plus_one):
            # get the number of players in this ward at this stage
            inf_ij = play_infections_i[j]

            if inf_ij > 0:
                # sample l players from this stage based on disease_progress
                l = _ran_binomial(rng, disease_progress, inf_ij)

                if l > 0:
                    # move l players from this stage to the next stage
                    play_infections_i_plus_one[j] += l
                    play_infections_i[j] -= l

Time since recovery
-------------------

To measure time since recovery, we can add extra "post-recovery" stages.
Individuals will be set to move slowly through those "post-recovery"
stages, until, when a particular post-recovery stage is reached, a fraction
of individuals are deemed to have lost immunity to the disease, and
are moved back to the S stage.

To do this, we will create a new version of the lurgy with these extra
post-recovery stages, which we will call ``R1`` to ``R10``. To do this
in Python, open ipython or jupyter and type;

.. code-block:: python

    >>> from metawards import Disease
    >>> lurgy = Disease("lurgy6")
    >>> lurgy.add("E", beta=0.0, progress=1.0)
    >>> lurgy.add("I1", beta=0.4, progress=0.2)
    >>> lurgy.add("I2", beta=0.5, progress=0.5, too_ill_to_move=0.5)
    >>> lurgy.add("I3", beta=0.5, progress=0.8, too_ill_to_move=0.8)
    >>> R_progress = 0.5
    >>> lurgy.add("R", progress=R_progress)
    >>> for i in range(1, 11):
    ...    lurgy.add(f"R{i}", beta=0.0, progress=R_progress)
    >>> lurgy.to_json("lurgy6.json", indent=2, auto_bzip=False)

or, in R/RStudio you could type;

.. code-block:: R

    > library(metawards)
    > lurgy <- metawards$Disease("lurgy6")
    > lurgy$add("E", beta=0.0, progress=1.0)
    > lurgy$add("I1", beta=0.4, progress=0.2)
    > lurgy$add("I2", beta=0.5, progress=0.5, too_ill_to_move=0.5)
    > lurgy$add("I3", beta=0.5, progress=0.8, too_ill_to_move=0.8)
    > R_progress <- 0.5
    > lurgy$add("R", progress=R_progress)
    > for(i in 1:10) {
          stage <- sprintf("R%d", i)
          lurgy$add(stage, beta=0.0, progress=R_progress)
      }
    > lurgy$to_json("lurgy6.json", indent=2, auto_bzip=False)

or simply copy the below into ``lurgy6.json``;

::

  {
    "name": "lurgy6",
    "stage": ["E", "I1", "I2", "I3", "R", "R1", "R2",
                "R3", "R4", "R5", "R6", "R7", "R8",
                "R9", "R10"],
    "mapping": ["E", "I", "I", "I", "R", "R", "R",
                "R", "R", "R", "R", "R", "R", "R",
                "R"],
    "beta": [0.0, 0.4, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "progress": [1.0, 0.2, 0.5, 0.8, 0.5, 0.5, 0.5,
                0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                0.5],
    "too_ill_to_move": [0.0, 0.0, 0.5, 0.8, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0],
    "contrib_foi": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    "start_symptom": 2,
    "is_infected": [true, true, true, true, false, false,
                    false, false, false, false, false, false,
                    false, false, false]
  }

This creates 10 post-recovery stages, with a progress parameter for these
stages of 0.5. This means, at quickest, an individual would take
10 days to progress from R to R10, but on average, this will take much
longer (as 50% move from one stage to the next each day).

Moving from R to S
------------------

Next, we need to add a move function that will move a fraction of individuals
from R10 back to S, to represent that fraction losing immunity.

Create a mover called ``move_immunity.py`` and copy in the below;

.. code-block:: python

    from metawards.movers import MoveGenerator, go_ward, MoveRecord
    from metawards.utils import Console


    def move_immunity(**kwargs):
        def go_immunity(**kwargs):
            record = MoveRecord()
            gen = MoveGenerator(from_stage="R10", to_stage="S",
                                fraction=0.5)

            go_ward(generator=gen, record=record, **kwargs)

            if len(record) > 0:
                nlost = record[0][-1]
                Console.print(f"{nlost} individual(s) lost immunity today")

        return [go_immunity]

This will move 50% of R10 individuals back to S each day.

.. note::

   We have used :class:`~metawards.movers.MoveRecord` to record the moves
   performed by :func:`~metawards.movers.go_ward`. This keeps a complete
   record of exactly how many individuals were moved, and the full
   details of that move. In this case, there will only be a single
   move (``record[0]``), and the number of individuals who were moved
   is the last value in the record (``record[0][-1]``).

You can run the model using;

.. code-block:: bash

   metawards -d lurgy6.json -m single -a 5 --move move_immunity.py

(using the single-ward model, seeding with 5 initial infection).

You should see that the outbreak oscillates as individuals who have
lost immunity are re-infected. For example, the graph I get
(from ``metawards-plot``) are shown below;

.. image:: ../../images/tutorial_9_3_1.jpg
   :alt: Outbreak trajectory when individuals can lose immunity

.. note::

   This is just an illustrative example. Individuals lose immunity in
   this model far more quickly than would be expected for a real disease.
   You could modify this example to use
   :doc:`custom user variables <../part02/02_adjustable>` to
   scan through different values of ``progress`` for each
   of the post-recovery stages, to better model a more realistic disease.

Vaccination and boosters
------------------------

You can apply the same method to model fading immunity after a vaccination.
This could be used to best plan how often booster doses should be deployed.

To do this, we will modify our lurgy to model to include vaccination
and post-vaccination stages. For example, in Python (in ipython/Jupyter);

.. code-block:: python

    >>> from metawards import Disease
    >>> lurgy = Disease("lurgy7")
    >>> lurgy.add("E", beta=0.0, progress=1.0)
    >>> lurgy.add("I1", beta=0.4, progress=0.2)
    >>> lurgy.add("I2", beta=0.5, progress=0.5, too_ill_to_move=0.5)
    >>> lurgy.add("I3", beta=0.5, progress=0.8, too_ill_to_move=0.8)
    >>> R_progress = 0.5
    >>> V_progress = 0.5
    >>> lurgy.add("R", progress=R_progress)
    >>> for i in range(1, 10):
    ...    lurgy.add(f"R{i}", beta=0.0, progress=R_progress)
    >>> lurgy.add("R10", beta=0.0, progress=0.0)
    >>> lurgy.add("V", progress=V_progress, is_infected=False)
    >>> for i in range(1, 11):
    ...     lurgy.add(f"V{i}", beta=0.0, progress=V_progress,
    ...               is_infected=False)
    >>> lurgy.to_json("lurgy7.json", auto_bzip=False)

.. code-block:: R

    > library(metawards)
    > lurgy <- metawards$Disease("lurgy7")
    > lurgy$add("E", beta=0.0, progress=1.0)
    > lurgy$add("I1", beta=0.4, progress=0.2)
    > lurgy$add("I2", beta=0.5, progress=0.5, too_ill_to_move=0.5)
    > lurgy$add("I3", beta=0.5, progress=0.8, too_ill_to_move=0.8)
    > R_progress <- 0.5
    > V_progress <- 0.5
    > lurgy$add("R", progress=R_progress)
    > for(i in 1:9) {
          stage <- sprintf("R%d", i)
          lurgy$add(stage, beta=0.0, progress=R_progress)
      }
    > lurgy.add("R10", beta=0.0, progress=0.0)
    > lurgy.add("V", progress=V_progress, is_infected=False)
    > for(i in 1:10) {
          stage <- sprintf("V%d", i)
          lurgy$add(stage, beta=0.0, progress=V_progress)
      }
    > lurgy.to_json("lurgy7.json", auto_bzip=False)

or copy the below into ``lurgy7.json``

::

  {
    "name": "lurgy7",
    "stage": ["E", "I1", "I2", "I3", "R", "R1",
              "R2", "R3", "R4", "R5", "R6", "R7",
              "R8", "R9", "R10", "V", "V1", "V2",
              "V3", "V4", "V5", "V6", "V7", "V8",
              "V9", "V10"],
    "mapping": ["E", "I", "I", "I", "R", "R", "R",
                "R", "R", "R", "R", "R", "R", "R",
                "R", "V", "V", "V", "V", "V", "V",
                "V", "V", "V", "V", "V"],
    "beta": [0.0, 0.4, 0.5, 0.5, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0],
    "progress": [1.0, 0.2, 0.5, 0.8, 0.5, 0.5,
                 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                 0.5, 0.5, 0.0, 0.5, 0.5, 0.5,
                 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                 0.5, 0.5],
    "too_ill_to_move": [0.0, 0.0, 0.5, 0.8, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0],
    "contrib_foi": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0, 1.0, 1.0],
    "start_symptom": 2,
    "is_infected": [true, true, true, true, false, false,
                    false, false, false, false, false, false,
                    false, false, false, false, false, false,
                    false, false, false, false, false, false,
                    false, false]
  }

.. note::

   Note how ``progress`` for the ``R10`` stage is set to 0 to prevent
   ``R10`` automatic progression from ``R10`` to ``V``.

.. note::

   Note also that we have to manually set ``is_infected`` to false for
   the V stages. This is set automatically to false only for R stages.

Next, modify ``move_immunity.py`` to read;

.. code-block:: python

    from metawards.movers import MoveGenerator, go_ward, MoveRecord
    from metawards.utils import Console


    def move_immunity(**kwargs):
        def go_vaccinate(population, **kwargs):
            if population.day <= 10:
                gen = MoveGenerator(from_stage="S", to_stage="V",
                                    number=100)
                go_ward(generator=gen, population=population, **kwargs)

        def go_immunity(**kwargs):
            record = MoveRecord()
            gen = MoveGenerator(from_stage="R10", to_stage="S",
                                fraction=0.5)

            go_ward(generator=gen, record=record, **kwargs)

            if len(record) > 0:
                Console.print(f"{record[0][-1]} individual(s) lost immunity today")

        def go_booster(**kwargs):
            gen = MoveGenerator(from_stage="V10", to_stage="V",
                                fraction=0.2)
            go_ward(generator=gen, **kwargs)

            record = MoveRecord()
            gen = MoveGenerator(from_stage="V10", to_stage="S")
            go_ward(generator=gen, record=record, **kwargs)

            if len(record) > 0:
                Console.print(f"{record[0][-1]} individual(s) didn't get their booster")

        return [go_vaccinate, go_booster, go_immunity]


Here, we've added a ``go_vaccinate`` function that, for the first
10 days of the outbreak, moves up to 100 individuals
per day from S to V. This will, in effect, vaccinate all of the
population of 1000 individuals who are not infected.

Next, we've added a ``go_booster`` function that samples 20% of the ``V10``
stage to give them a booster vaccine that returns them to ``V``. The
remaining individuals in ``V10`` miss their booster dose, and are
returned to ``S``.

You can run the model using;

.. code-block:: bash

   metawards -d lurgy7.json -m single -a 5 --move move_immunity.py

You should see that the infection nearly dies out, as nearly everyone
is vaccinated. However, a small number of lingering infections spark
a second outbreak amongst individuals who miss their booster shot,
leading then to a cycle of infection and losing immunity, e.g.

::

    ─────────────────────────────────────────────── Day 9 ────────────────────────────────────────────────
    S: 89  E: 0  I: 6  V: 900  R: 5  IW: 0  POPULATION: 1000
    Number of infections: 6

    ─────────────────────────────────────────────── Day 10 ───────────────────────────────────────────────
    S: 0  E: 0  I: 6  V: 989  R: 5  IW: 0  POPULATION: 1000
    Number of infections: 6

    ─────────────────────────────────────────────── Day 11 ───────────────────────────────────────────────
    S: 0  E: 0  I: 6  V: 989  R: 5  IW: 0  POPULATION: 1000
    Number of infections: 6

    ─────────────────────────────────────────────── Day 12 ───────────────────────────────────────────────
    1 individual(s) didn't get their booster
    S: 1  E: 0  I: 6  V: 988  R: 5  IW: 0  POPULATION: 1000
    Number of infections: 6

    ...

    ─────────────────────────────────────────────── Day 23 ───────────────────────────────────────────────
    54 individual(s) didn't get their booster
    S: 309  E: 0  I: 2  V: 680  R: 9  IW: 0  POPULATION: 1000
    Number of infections: 2

    ─────────────────────────────────────────────── Day 24 ───────────────────────────────────────────────
    72 individual(s) didn't get their booster
    2 individual(s) lost immunity today
    S: 383  E: 0  I: 2  V: 608  R: 7  IW: 0  POPULATION: 1000
    Number of infections: 2

    ─────────────────────────────────────────────── Day 25 ───────────────────────────────────────────────
    60 individual(s) didn't get their booster
    S: 443  E: 0  I: 2  V: 548  R: 7  IW: 0  POPULATION: 1000
    Number of infections: 2

    ─────────────────────────────────────────────── Day 26 ───────────────────────────────────────────────
    64 individual(s) didn't get their booster
    1 individual(s) lost immunity today
    S: 507  E: 1  I: 2  V: 484  R: 6  IW: 1  POPULATION: 1000
    Number of infections: 3

    ─────────────────────────────────────────────── Day 27 ───────────────────────────────────────────────
    49 individual(s) didn't get their booster
    1 individual(s) lost immunity today
    S: 557  E: 0  I: 3  V: 435  R: 5  IW: 0  POPULATION: 1000
    Number of infections: 3

    ─────────────────────────────────────────────── Day 28 ───────────────────────────────────────────────
    44 individual(s) didn't get their booster
    S: 599  E: 2  I: 3  V: 391  R: 5  IW: 1  POPULATION: 1000
    Number of infections: 5

    ...

    ─────────────────────────────────────────────── Day 40 ───────────────────────────────────────────────
    8 individual(s) didn't get their booster
    S: 776  E: 7  I: 39  V: 162  R: 16  IW: 1  POPULATION: 1000
    Number of infections: 46

    ─────────────────────────────────────────────── Day 41 ───────────────────────────────────────────────
    8 individual(s) didn't get their booster
    S: 776  E: 8  I: 45  V: 154  R: 17  IW: 1  POPULATION: 1000
    Number of infections: 53

    ...

    ────────────────────────────────────────────── Day 103 ───────────────────────────────────────────────
    11 individual(s) lost immunity today
    S: 289  E: 32  I: 329  V: 2  R: 348  IW: 1  POPULATION: 1000
    Number of infections: 361

    ────────────────────────────────────────────── Day 104 ───────────────────────────────────────────────
    12 individual(s) lost immunity today
    S: 258  E: 43  I: 319  V: 2  R: 378  IW: 1  POPULATION: 1000
    Number of infections: 362

    ────────────────────────────────────────────── Day 105 ───────────────────────────────────────────────
    1 individual(s) didn't get their booster
    11 individual(s) lost immunity today
    S: 233  E: 37  I: 325  V: 1  R: 404  IW: 1  POPULATION: 1000
    Number of infections: 362

.. note::

   Again, this is just an illustrative example. Immunity from vaccination would
   be expected to last for much longer than a couple of weeks.
   You could use adjustable variables (and custom user-adjustable
   variables) to scan through the progress along the post-vaccination
   stages, the numbers vaccinated each day, and different percentages
   of individuals who take a booster, to better model a real situation.

.. note::

   We have used :meth:`Disease.progress <metawards.Disease.progress>`
   to slow movement along the post-recovery and post-vaccinated stages.
   An alternative method would be to write a custom iterator to
   replace :func:`~metawards.iterators.advance_recovery`. This could
   slow down movement programmatically, e.g. by only testing for
   advancement along the ``R`` and ``V`` stages every 10 days, as
   opposed to every day. We've designed ``metawards`` to be very
   flexible, so that you have many choices for how you want to
   model different scenarios.
======================
Extracting by location
======================

Up to this point, we have only looked at the summary data from each
*model run*. There is a lot more that can be explored, as ``metawards``
models the outbreak in every ward in the model (e.g. every electoral
ward in the UK).

Metadata about those wards is loaded into
:meth:`~metawards.Network.info` object, which is of class
:class:`~metawards.WardInfos`. This class can be queried to get
the index of wards according to their name, their official ID code,
or the local authority or region in which it belongs.

Searching using WardInfo
------------------------

For example, we can find all of the wards that match the name
"Clifton" using the below code. Open up ``ipython`` or a Jupyter
notebook and type;

.. code-block:: python

   >>> from metawards import Network, Parameters
   >>> params = Parameters.load()
   >>> params.set_input_files("2011Data")
   >>> params.set_disease("lurgy")
   >>> network = Network.build(params)

This has now built a network object that you can query, e.g.

.. code-block:: python

   >>> clifton = network.info.find("Clifton")
   >>> print(clifton)
   [154, 403, 829, 3612, 3662, 3670, 3703, 3766, 3974, 3975, 8134, 8327, 8328]
   >>> for ward in clifton:
   ...     print(network.info[ward])
   WardInfo(name='Clifton-with-Maidenway', alternate_names=['Clifton-with-Maidenway'], code='E05002101', alternate_codes=['E36000654'], authority='Torbay', authority_code='E06000027', region='', region_code='')
   WardInfo(name='Clifton', alternate_names=['Clifton'], code='E05003119', alternate_codes=['E36001870'], authority='Allerdale', authority_code='E07000026', region='', region_code='')
   WardInfo(name='Clifton and Bradley', alternate_names=['Clifton and Bradley'], code='E05003350', alternate_codes=['E36002100'], authority='Derbyshire Dales', authority_code='E07000035', region='', region_code='')
   WardInfo(name='Skelton, Rawcliffe and Clifton Without', alternate_names=['Skelton, Rawcliffe and Clifton Without'], code='E05001763', alternate_codes=['E36000299'], authority='York', authority_code='E06000014', region='', region_code='')
   WardInfo(name='Clifton', alternate_names=['Clifton'], code='E05001980', alternate_codes=['E36000533'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Clifton East', alternate_names=['Clifton East'], code='E05001981', alternate_codes=['E36000534'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Clifton', alternate_names=['Clifton'], code='E05001747', alternate_codes=['E36000283'], authority='York', authority_code='E06000014', region='', region_code='')
   WardInfo(name='Clifton', alternate_names=['Clifton'], code='E05001648', alternate_codes=['E36000184'], authority='Blackpool', authority_code='E06000009', region='', region_code='')
   WardInfo(name='Clifton North', alternate_names=['Clifton North'], code='E05001831', alternate_codes=['E36000367'], authority='Nottingham', authority_code='E06000018', region='', region_code='')
   WardInfo(name='Clifton South', alternate_names=['Clifton South'], code='E05001832', alternate_codes=['E36000368'], authority='Nottingham', authority_code='E06000018', region='', region_code='')
   WardInfo(name='Clifton', alternate_names=['Clifton'], code='E05005188', alternate_codes=['E36003801'], authority='Fylde', authority_code='E07000119', region='', region_code='')
   WardInfo(name='Cliftonville East', alternate_names=['Cliftonville East'], code='E05005087', alternate_codes=['E36003700'], authority='Thanet', authority_code='E07000114', region='', region_code='')
   WardInfo(name='Cliftonville West', alternate_names=['Cliftonville West'], code='E05005088', alternate_codes=['E36003701'], authority='Thanet', authority_code='E07000114', region='', region_code='')

This has returned all wards that have "Clifton" in the name. The search is
actually performed as a `regular expression <https://chryswoods.com/intermediate_python/regexp.html>`__,
and is case-insensitive. You can pass a regular expression string directly,
e.g. ```r"^(Clifton)$"``` would match "Clifton" at the beginning (```^```) and
end (```$```) of the string, i.e. it only matches wards that exactly
match "Clifton". Try this by typing;

.. code-block:: python

   >>> clifton = network.info.find(r"^(Clifton)$")
   >>> for ward in clifton:
   ...    print(network.info[ward])
   WardInfo(name='Clifton', alternate_names=['Clifton'], code='E05003119', alternate_codes=['E36001870'], authority='Allerdale', authority_code='E07000026', region='', region_code='')
   WardInfo(name='Clifton', alternate_names=['Clifton'], code='E05001980', alternate_codes=['E36000533'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Clifton', alternate_names=['Clifton'], code='E05001747', alternate_codes=['E36000283'], authority='York', authority_code='E06000014', region='', region_code='')
   WardInfo(name='Clifton', alternate_names=['Clifton'], code='E05001648', alternate_codes=['E36000184'], authority='Blackpool', authority_code='E06000009', region='', region_code='')
   WardInfo(name='Clifton', alternate_names=['Clifton'], code='E05005188', alternate_codes=['E36003801'], authority='Fylde', authority_code='E07000119', region='', region_code='')

Searching by local authority
----------------------------

There are many "Clifton"s in the UK. We can limit to the "Clifton" in Bristol
by specifying the local authority. Again, this can be a regular expression,
although in this case just searching for "Bristol" will be enough.

.. code-block:: python

   >>> clifton = network.info.find(r"^(Clifton)$", authority="Bristol")
   >>> for ward in clifton:
   ...     print(ward)
   WardInfo(name='Clifton', alternate_names=['Clifton'], code='E05001980', alternate_codes=['E36000533'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')

.. note::
   In this case the dataset does not include regional data (the ```region```
   is empty). If regional data was available you could search by region
   using ```network.info.find("Clifton", region="South West")``.

This searching is very powerful. For example, we can now search for all
wards that are in the same local authority as "Clifton, Bristol", e.g.

.. code-block:: python

   >>> clifton = network.info.find(r"^(Clifton)$", authority="Bristol")[0]
   >>> clifton = newwork.info[clifton]
   >>> authority_code = clifton.authority_code
   >>> wards = network.info.find(authority=authority_code)
   >>> for ward in wards:
   ...     print(ward)
   WardInfo(name='Brislington West', alternate_names=['Brislington West'], code='E05001978', alternate_codes=['E36000531'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Cabot', alternate_names=['Cabot'], code='E05001979', alternate_codes=['E36000532'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Clifton', alternate_names=['Clifton'], code='E05001980', alternate_codes=['E36000533'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Clifton East', alternate_names=['Clifton East'], code='E05001981', alternate_codes=['E36000534'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Cotham', alternate_names=['Cotham'], code='E05001982', alternate_codes=['E36000535'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Easton', alternate_names=['Easton'], code='E05001983', alternate_codes=['E36000536'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Eastville', alternate_names=['Eastville'], code='E05001984', alternate_codes=['E36000537'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Filwood', alternate_names=['Filwood'], code='E05001985', alternate_codes=['E36000538'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Frome Vale', alternate_names=['Frome Vale'], code='E05001986', alternate_codes=['E36000539'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Hartcliffe', alternate_names=['Hartcliffe'], code='E05001987', alternate_codes=['E36000540'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Henbury', alternate_names=['Henbury'], code='E05001988', alternate_codes=['E36000541'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Hengrove', alternate_names=['Hengrove'], code='E05001989', alternate_codes=['E36000542'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Henleaze', alternate_names=['Henleaze'], code='E05001990', alternate_codes=['E36000543'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Hillfields', alternate_names=['Hillfields'], code='E05001991', alternate_codes=['E36000544'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Horfield', alternate_names=['Horfield'], code='E05001992', alternate_codes=['E36000545'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Kingsweston', alternate_names=['Kingsweston'], code='E05001993', alternate_codes=['E36000546'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Ashley', alternate_names=['Ashley'], code='E05001972', alternate_codes=['E36000525'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Avonmouth', alternate_names=['Avonmouth'], code='E05001973', alternate_codes=['E36000526'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Bedminster', alternate_names=['Bedminster'], code='E05001974', alternate_codes=['E36000527'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Bishopston', alternate_names=['Bishopston'], code='E05001975', alternate_codes=['E36000528'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Bishopsworth', alternate_names=['Bishopsworth'], code='E05001976', alternate_codes=['E36000529'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Brislington East', alternate_names=['Brislington East'], code='E05001977', alternate_codes=['E36000530'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Knowle', alternate_names=['Knowle'], code='E05001994', alternate_codes=['E36000547'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Lawrence Hill', alternate_names=['Lawrence Hill'], code='E05001995', alternate_codes=['E36000548'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Lockleaze', alternate_names=['Lockleaze'], code='E05001996', alternate_codes=['E36000549'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Redland', alternate_names=['Redland'], code='E05001997', alternate_codes=['E36000550'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='St George East', alternate_names=['St George East'], code='E05001998', alternate_codes=['E36000553'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='St George West', alternate_names=['St George West'], code='E05001999', alternate_codes=['E36000554'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Southmead', alternate_names=['Southmead'], code='E05002000', alternate_codes=['E36000551'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Southville', alternate_names=['Southville'], code='E05002001', alternate_codes=['E36000552'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Stockwood', alternate_names=['Stockwood'], code='E05002002', alternate_codes=['E36000555'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Stoke Bishop', alternate_names=['Stoke Bishop'], code='E05002003', alternate_codes=['E36000556'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Westbury-on-Trym', alternate_names=['Westbury-on-Trym'], code='E05002004', alternate_codes=['E36000557'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Whitchurch Park', alternate_names=['Whitchurch Park'], code='E05002005', alternate_codes=['E36000558'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')
   WardInfo(name='Windmill Hill', alternate_names=['Windmill Hill'], code='E05002006', alternate_codes=['E36000559'], authority='Bristol, City of', authority_code='E06000023', region='', region_code='')

.. note::

   It is true that we could have achieved this by just searching for
   Bristol alone. However, this method of searching for ward+authority
   is more robust against multiple authorities having similar names.
   For example, searching for the authority "Newcastle" returns
   both "Newcastle upon Tyne" and "Newcastle-under-Lyme".

Using location in an extractor
------------------------------

We can use the above search to track the total number of infections in each
of the wards in Bristol. Create a new python file called ``location.py``
and copy in the below;

.. code-block:: python

    matched_wards = None
    headers = []

    def output_location(network, population, workspace, output_dir, **kwargs):
        ward = "clifton"
        authority = "bristol"

        global matched_wards, headers

        if matched_wards is None:
            # This is performed only once, when this function is first called
            ward = network.info.find(name=ward, authority=authority)[0]
            ward = network.info[ward]
            authority_code = ward.authority_code
            matched_wards = network.info.find(authority=authority_code)

            headers = []
            headers.append("day")

            for ward in matched_wards:
                headers.append(f"'{network.info[ward].name}'")

        # open the file called "authority.dat", e.g. "bristol.dat"
        # Note we are using comma separators and have put the ward
        # names in single quotes to make the output easier to parse
        locfile = output_dir.open(f"{authority}.dat", headers=headers, sep=",")

        locfile.write(str(population.day))

        for ward in matched_wards:
            total = workspace.total_inf_ward[ward]
            locfile.write("," + str(total))

        locfile.write("\n")

    def extract_location(**kwargs):
        from metawards.extractors import extract_default

        return extract_default(**kwargs) + [output_location]

Save the file and run ``metawards`` using this extractor via

.. code-block:: bash

   metawards --extractor location

You should see that a new output file called ``bristol.dat.bz2`` was
created. Loading this up into pandas should show;

.. code-block:: python

   >>> import pandas as pd
   >>> df = pd.read_csv("output/bristol.dat.bz2", index_col="day")
   >>> print(df)
         'Brislington West'  'Cabot'  ...  'Whitchurch Park'  'Windmill Hill'
    day                               ...
    0                     0        0  ...                  0                0
    1                     0        0  ...                  0                0
    2                     0        0  ...                  0                0
    3                     0        0  ...                  0                0
    4                     0        0  ...                  0                0

    [5 rows x 35 columns]

=================
Custom extractors
=================

You have now learned how to use custom iterators to customise the
advancement of the outbreak during a model run.

In a similar way, ``metawards`` provides custom extractors that
enable you to customise the output that is produced and written
to a file (or files).

Hello extractors
----------------

You create an extractor in an almost identical manner as an iterator.
Start by creating a python file called ``hello.py`` and copy in the
below;

.. code-block:: python

  from metawards.utils import Console

  def extract_hello(**kwargs):
      Console.print("Hello extract_hello")

      return []

The extractor is passed using the ``--extractor`` command-line argument.
Run ``metawards`` using;

.. code-block:: bash

   metawards --extractor hello

You should see output something similar to this;

::

    Importing a custom extractor from hello
    Loaded hello from hello.py
    <function extract_hello at 0x1068599e0>
    Building a custom extractor for <function extract_hello at 0x1068599e0>
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 0 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello extract_hello
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 1 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello extract_hello
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello extract_hello
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 3 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello extract_hello
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 4 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Hello extract_hello
    S: 56082077  E: 0  I: 0  R: 0  IW: 0  POPULATION: 56082077
    Number of infections: 0
    Infection died ... Ending on day 5

extract_XXX and output_XXX
--------------------------

At the end of each model day, ``metawards`` calls the
:meth:`~metawards.utils.extract` function. This calls your ``extract_XXX``
function. The signature is very similar to the custom iterator functions,
namely it should take ``**kwargs``, and then return a list of functions
that :meth:`~metawards.utils.extract` will then call to output data
(what we term ``output_XXX`` functions).

At the moment, nothing is being written to the output directory. We
can change this by adding an ``output_XXX`` function. For example,
create a new python file called ``population.py`` and copy in
the below;

.. code-block:: python

    from metawards.utils import Console

    def output_population(population, output_dir, **kwargs):
        Console.debug("Hello output_population")

        # create an output file called 'population.dat'
        popfile = output_dir.open("population.dat")

        # write the population to this file
        popfile.write(f"{population.day} {population.date.isoformat()} "
                      f"{population.susceptibles} {population.latent} "
                      f"{population.total} {population.recovereds}\n")

    def extract_population(**kwargs):
        Console.debug("hello extract_population")

        return [output_population]

This defines two functions;

* ``extract_population``, which tells ``metawards`` to use your
  ``output_population`` function,

* and ``output_population`` that uses the passed
  :class:`population <metawards.Population>` and
  :class:`output_dir <metawards.OutputFiles>` objects to write
  the population of the different disease states to a file
  in the output directory called ``population.dat``.

Use this extractor using the command;

.. code-block:: bash

   metawards --extractor population

If you take a look in the ``output`` directory you should see that a file
called ``population.dat.bz2`` has been created. You can take a look at
this in R, Python pandas or excel. For example, we can load this in
pandas using;

.. code-block:: python

   >>> import pandas as pd
   >>> df = pd.read_csv("output/population.dat.bz2", sep=" ", header=None)
   >>> print(df)
         0           1         2  3  4  5
      0  0  2020-04-26  56082077  0  0  0
      1  1  2020-04-27  56082077  0  0  0
      2  2  2020-04-28  56082077  0  0  0
      3  3  2020-04-29  56082077  0  0  0
      4  4  2020-04-30  56082077  0  0  0

.. note::
   ``metawards`` will auto-compress all files written into the output
   directory. If you don't want this, then use the command-line argument
   ``--no-auto-bzip``.

Notice that there are no headers to the columns. We can add a header
by passing in the headers to the
:meth:`~metawards.OutputFiles.open` function, e.g. change ``population.py``
to read;

.. code-block:: python

    from metawards.utils import Console

    def output_population(population, output_dir, **kwargs):
        Console.debug("Hello output_population")

        # create an output file called 'population.dat'
        popfile = output_dir.open("population.dat",
                                  headers=["day", "date", "S", "E",
                                           "I", "R"])

        # write the population to this file
        popfile.write(f"{population.day} {population.date.isoformat()} "
                      f"{population.susceptibles} {population.latent} "
                      f"{population.total} {population.recovereds}\n")

    def extract_population(**kwargs):
        Console.debug("hello extract_population")

        return [output_population]

Run ``metawards`` again, and now if you load the ``population.dat.bz2``
file into pandas (or R or Excel) you will see something similar to;

.. code-block:: python

  >>> import pandas as pd
  >>> df = pd.read_csv("output/population.dat.bz2", sep=" ", index_col="day")
  >>> print(df)
               date         S  E  I  R
    day
    0    2020-04-26  56082077  0  0  0
    1    2020-04-27  56082077  0  0  0
    2    2020-04-28  56082077  0  0  0
    3    2020-04-29  56082077  0  0  0
    4    2020-04-30  56082077  0  0  0

.. note::
  Note how I have used ``index_col`` to set the ``day`` as the index
  in pandas

Occasional functions
--------------------

Just as with iterators, we can choose to only call the output function
on specific days. For example, to only output the population to the
file on even days, change ``population.py`` to read;

.. code-block:: python

    from metawards.utils import Console

    def output_population(population, output_dir, **kwargs):
        Console.debug("Hello output_population")

        # create an output file called 'population.dat'
        popfile = output_dir.open("population.dat",
                                headers=["day", "date", "S", "E",
                                        "I", "R"])

        # write the population to this file
        popfile.write(f"{population.day} {population.date.isoformat()} "
                    f"{population.susceptibles} {population.latent} "
                    f"{population.total} {population.recovereds}\n")


    def extract_population(population, **kwargs):
        Console.debug("hello extract_population")

        if population.day % 2 == 0:
            return [output_population]
        else:
            return []

Run ``metawards`` using this extractor and you should see that the
``population.dat.bz2`` file contains output only for days 0, 2, and 4.

.. note::
   The line ``population.day % 2 == 0`` takes the remainder division
   of ``population.day`` with 2. Any day that is divisible by 2 will
   return 0. You can output every ``N`` days using
   ``population.day % N == 0``.

.. note::
   You are also able to only print out on other conditions, e.g.
   when the *model run* reaches a certain date, or when the
   infected population grows above a certain size.

Exiting early
-------------

Sometimes you may want to exit a *model run* early if a condition
is reached. The best way to do this is to raise a Python
`StopIteration <https://docs.python.org/3/library/exceptions.html>`__
exception. This will signal to ``metawards`` that the *model run* should
stop at the end of the current iteration (other functions that are
part of that iteration can still complete, and any output written
for that iteration will still be recorded).

For example, you could use this output function to stop the *model run*
once the number of infections reaches 2000. Copy the below into
``extract_stop.py``;

.. code-block:: python

    from metawards.extractors import extract_default

    def output_stop(population, **kwargs):
        if population.infecteds > 2000:
            raise StopIteration


    def extract_stop(**kwargs):
        output_funcs = extract_default(**kwargs)

        output_funcs.append(output_stop)

        return output_funcs

This extractor uses all of the functions of
:meth:`~metawards.extractors.extract_default`, plus a new custom
output function called ``output_stop``. This compares the number
of infections (:data:`population.infecteds <metawards.Population.infecteds>`),
and if this is more than 2000, then it raises a Python
`StopIteration <https://docs.python.org/3/library/exceptions.html>`__.

Run ``metawards`` using;

.. code-block:: bash

    metawards -d lurgy3 -a ExtraSeedsLondon.dat --extractor extract_stop

You should see that the *model run* is stopped once the number of infections
is greater than 2000, e.g.

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 29 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56078417  E: 566  I: 1275  R: 1819  IW: 501  POPULATION: 56082077
    Number of infections: 1841

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 30 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56077705  E: 650  I: 1555  R: 2167  IW: 543  POPULATION: 56082077
    <function output_stop at 0x105412e60> has indicated that the model run should stop early. Will finish the
    run at the end of this iteration
    Number of infections: 2205
    Exiting model run early due to function request
    Infection died ... Ending on day 31

You can use this to stop a *model run* for any reason you want, e.g.
a calculated condition has been reached, the model is unstable or
uses parameters that are uninteresting. Another option is to use this to
stop ``metawards`` from running for more than a specified amount of time.

To do this, create an extractor called ``extract_stop_time.py`` and
copy in;

.. code-block:: python

    from metawards.extractors import extract_default
    from metawards.utils import Console

    from datetime import datetime


    def output_stop_time(network, **kwargs):
        if not hasattr(network.params, "_start_model_time"):
            network.params._start_model_time = datetime.now()
            return

        runtime = datetime.now() - network.params._start_model_time

        Console.print(f"Runtime is {runtime.total_seconds()} seconds")

        if runtime.total_seconds() > 5:
            Console.warning(f"Runtime exceeded 5 seconds!")
            raise StopIteration


    def extract_stop_time(**kwargs):
        output_funcs = extract_default(**kwargs)

        output_funcs.append(output_stop_time)

        return output_funcs


This uses the Python
`datetime <https://docs.python.org/3/library/datetime.html>`__ module to
calculate the time since ``output_stop_time`` was first called.

.. note::

    We've recorded this start time by adding an attribute to ``network.params``
    called ``_start_model_time``. Adding attributes like this to the
    ``network.params`` object is a good way to store parameters between
    model runs, or to initialise values at the start of a model run.
    Any parameters are guaranteed to be cleared between runs, and
    the threading model means that anything you read/write is thread
    safe and will not interfere with other runs.

Run this extractor using;

.. code-block:: bash

    metawards -d lurgy3 -a ExtraSeedsLondon.dat --extractor extract_stop_time

You should see that the run ends after five seconds, e.g.;

::

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 38 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56064800  E: 2313  I: 5934  R: 9030  IW: 1784  POPULATION: 56082077
    Runtime is 4.538544 seconds
    Number of infections: 8247

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 39 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56061567  E: 2816  I: 7023  R: 10671  IW: 2026  POPULATION: 56082077
    Runtime is 4.831688 seconds
    Number of infections: 9839

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Day 40 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    S: 56057698  E: 3233  I: 8359  R: 12787  IW: 2306  POPULATION: 56082077
    Runtime is 5.156103 seconds

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ WARNING ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Runtime exceeded 5 seconds!

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    <function output_stop_time at 0x10aa3ec20> has indicated that the model run should stop early. Will
    finish the run at the end of this iteration
    Number of infections: 11592
    Exiting model run early due to function request
    Infection died ... Ending on day 41
=============================
Accumulating into a Workspace
=============================

The default function that is always called by
:meth:`~metawards.utils.extract` is
:meth:`metawards.extractors.output_core`.

This core extractor performs the bulk of the work of accumulating all of
the infection data into a single :class:`metawards.Workspace` object.

This :class:`metawards.Workspace` object contains;

* :meth:`~metawards.Workspace.inf_tot`: The total population size at each
  of the different disease stages from the ``work`` infections.

* :meth:`~metawards.Workspace.pinf_tot`: The total population size at each
  of the different disease stages from the ``play`` infections.

* :meth:`~metawards.Workspace.n_inf_wards`: The number of wards with at
  least one member for each disease stage.

* :meth:`~metawards.Workspace.total_inf_ward`: The size of the infected
  population in each ward (the prevalence).

* :meth:`~metawards.Workspace.total_new_inf_ward`: The number of new infections
  on this day in each ward.

* :meth:`~metawards.Workspace.incidence`: The incidence of infection in each
  ward.

Output incidence
----------------

This :class:`~metawards.Workspace` contains data that can be easily output,
e.g. the :meth:`metawards.extractors.output_incidence` supplied extractor
writes the incidence to a file called ``incidence.dat.bz2``. For example,
you can call this from your ``population.py`` extractor by changing it
to read;

.. code-block:: python

    from metawards.utils import Console

    def output_population(population, output_dir, **kwargs):
        Console.debug("Hello output_population")

        # create an output file called 'population.dat'
        popfile = output_dir.open("population.dat",
                                  headers=["day", "date", "S", "E",
                                           "I", "R"])

        # write the population to this file
        popfile.write(f"{population.day} {population.date.isoformat()} "
                    f"{population.susceptibles} {population.latent} "
                    f"{population.total} {population.recovereds}\n")


    def extract_population(population, **kwargs):
        Console.debug("hello extract_population")

        from metawards.extractors import output_incidence

        if population.day % 2 == 0:
            return [output_population, output_incidence]
        else:
            return [output_incidence]

.. note::
   See how we are calling ``output_incidence`` on every day, but
   ``output_population`` only on even days.

If you run this extractor using

.. code-block:: bash

   metawards --extractor population

You will now see that you get a file called ``incidence.dat.bz2``
in the output directory. This will be a big matrix of mostly zeroes,
as no infection has been seeded.

Default outputs
---------------

The default extractor is :meth:`~metawards.extractors.extract_default`.
This returns;

* :meth:`~metawards.extractors.output_basic`: Writes out basic information
  to the files ``NumberWardsInfected.dat``, ``TotalInfections.dat`` etc.

* :meth:`~metawards.extractors.output_dispersal`: Calculates and writes out
  the geographic disperal of the outbreak to ``MeanXY.dat``, ``VarXY.day``
  and ``Dispersal.dat``

* :meth:`~metawards.extractors.output_prevalence`: Writes the (large)
  prevalence matrix to ``prevalence.dat``.

* :meth:`~metawards.extractors.output_incidence`: Writes the (large) incidence
  matrix to ``incidence.dat``.

You can use :meth:`~metawards.extractors.extract_default`
either by not supplying an extractor
via the ``--extractor`` command line argument, or by specifying
``--extract-default``.

Have a go using;

.. code-block:: bash

   metawards --extractor extract_default

As well as :meth:`~metawards.extractors.extract_default`, there is also
:meth:`~metawards.extractors.extract_small` (only extracting the
"small" files),
"meth:`~metawards.extractors.extract_large` (extract everything, including
producing large trajectory files) and
:meth:`~metawards.extractors.extract_none` (extract nothing - useful if
you want to restrict output only to ``results.csv.bz2``).
========================
Extracting rates by ward
========================

One of the major applications of a custom extractor is to enable you
to add your own live-analysis that is performed while the simulation
is running. This is much faster than having ``metawards`` write all
data to disk, and then running analysis in a post-processing step.

To analyse data, you need to understand the
:class:`~metawards.Workspace` class.

Workspace
---------

The :class:`~metawards.Workspace` class provides a workspace in which data
is accumulated and analysis performed. It is cleared between iterations,
and holds the following information;

Data indexed by disease stage:
    * ``inf_tot``: The total number of workers in each disease stage, summed over all wards.

    * ``pinf_tot``: The total number of players in each disease stage, summed over all wards.

    * ``n_inf_wards``: The number of wards which have at least one individual in the disease stage,

Data indexed by ward:
    * ``total_inf_ward``: The total number of infections in each ward

    * ``total_new_inf_ward``: The total number of new infections in each ward

    * ``incidence``: The incidence of infection in each ward

    * ``S_in_wards``: The total population in the ``S`` state in each ward

    * ``E_in_wards``: The total population in the ``E`` state in each ward

    * ``I_in_wards``: The total population in the ``I`` state in each ward

    * ``R_in_wards``: The total population in the ``R`` state in each ward

Data indexed by disease stage and then by ward
   * ``ward_inf_tot``: The total population in each disease stage in each ward

Extracting the population in I1
-------------------------------

We can use a custom extractor to report the total number of individuals
who are in the first ``I`` stage (I1) for each ward for each day.

To do this, create a new extractor called ``extract_i1.py`` and copy
in the below;

.. code-block:: python

    from metawards.extractors import extract_default


    def output_i1(population, workspace, output_dir, **kwargs):
        # Open the file "total_i1.csv" in the output directory
        FILE = output_dir.open("total_i1.csv")

        ward_inf_tot = workspace.ward_inf_tot

        # The I1 state is stage 2
        I1_inf_tot = ward_inf_tot[2]

        FILE.write(str(population.day) + ",")
        FILE.write(",".join([str(x) for x in I1_inf_tot]) + "\n")


    def extract_i1(**kwargs):
        # return all of the functions from "extract_default"
        # plus our new "output_i1"
        funcs = extract_default(**kwargs)
        funcs.append(output_i1)
        return funcs

This defines a new output function called ``output_i1``. This calls the
``open`` function of the :class:`~metawards.OutputFiles` object held
in ``output_dir`` to open the file ``total_i1.csv`` in the output
directory.

.. note::

   The :class:`~metawards.OutputFiles` class is clever enough to only
   open the file once, and will automatically close it when needed. It
   will also ensure that the file is opened in the correct output
   directory for a *model run* and will compress the file using bz2
   if the ``--auto-bzip`` command-line option has been passed (the default),
   or will disable automatic compression if ``--no-auto-bzip`` is set.

The function then gets ``I1_inf_tot`` from the third disease stage data
held in ``workspace.ward_inf_tot``. The third stage is the first ``I``
stage as the first stage (``ward_inf_tot[0]``) is a special ``*`` stage,
used for advanced bookkeeping, while the second stage (``ward_inf_tot[1]``)
is the latent, ``E`` stage.

.. note::

   The ``*`` stage is used to help evaluate how many wards see new infections.
   Individuals are moved into the ``*`` stage at the point of infection,
   and are moved into the ``E`` stage on the day after infection. By default,
   individuals in the ``*`` stage are counted into ``R``, which is why this
   does not appear to rise continuously. You can control how individuals
   in the ``*`` stage are counted using either ``--star-is-E`` to count them
   as ``E`` (additional latent stage), ``--star-is-R`` to count them as
   ``R`` (default behaviour) or ``--disable-star`` to remove the ``*``
   stage and instead treat this first stage as the one and only ``E`` stage.
   Note that if you do this, you will need to update the disease files to
   remove the ``*`` stage, and to update your above extractor as now
   stage 1 will be the first ``I`` stage.

Now that we have the population in the ``I1`` stage in each ward in
``I1_inf_tot``, we write this as a comma-separated line to the file,
starting each line with the day number obtained from the passed
:class:`~metawards.Population` typed ``population`` object.

To use your extractor run ``metawards`` using;

.. code-block:: bash

   metawards -d lurgy3 --extract extract_i1 -a ExtraSeedsLondon.dat --nsteps 30

.. note::

   Note that we've set ``nsteps`` to 30 for illustration only,
   just to limit the runtime and the size of the file. In a real production
   run you wouldn't need to set the number of steps.

You should see that your file called ``total_i1.csv.bz2`` has been created
in the output directory, and that this contains the populations of the ``I1``
state for each ward.

Calculating rates
-----------------

As well as outputting raw data, you can also perform some simple analysis
that is run live during the *model run*. For example, you may want to record
the number of individuals entering each state, so that you can calculate
the rate of progress across states.

To do this, you will need to save the ``ward_inf_tot`` data from the previous
day's state. You can do this by adding it as a custom attribute to the
workspace.

Create a new extractor by creating ``extract_rate.py`` and copying in the below;

.. code-block:: python

    from metawards.extractors import extract_default
    from copy import deepcopy


    def output_rate(population, workspace, output_dir, **kwargs):
        if not hasattr(workspace, "output_rate_previous"):
            # This is the first day, so we cannot calculate the rate.
            # Instead, just save today's data so that it can be
            # be used tomorrow
            workspace.output_rate_previous = deepcopy(workspace.ward_inf_tot)
            return

        # get yesterday's data
        ward_inf_previous = workspace.output_rate_previous

        # get today's data
        ward_inf_tot = workspace.ward_inf_tot

        # calculate and write the difference between the two to files for
        # each disease stage...
        for i in range(0, workspace.n_inf_classes):
            FILE = output_dir.open(f"rate_{i}.csv")

            FILE.write(str(population.day))

            # loop over the data for each ward and write the
            # difference to the file
            for old, new in zip(ward_inf_previous[i], ward_inf_tot[i]):
                FILE.write("," + str(new - old))

            FILE.write("\n")

        # save today's data so that it can be used tomorrow
        workspace.output_rate_previous = deepcopy(ward_inf_tot)


    def extract_rate(**kwargs):
        funcs = extract_default(**kwargs)
        funcs.append(output_rate)
        return funcs

This extractor looks a little more complex, but it builds on what you have
seen before. It defines ``output_rate``, which if the function that will
output the rates, and ``extract_rate`` which returns all of the functions
from ``extract_default``, plus your new ``output_rate`` function.

The first job of ``output_rate`` is to determine if it has been called on
the first day of the model. If it has, then there is no previous data
from "yesterday" that can be used to calculate the rate. The function
detects if this is the case by checking for a new custom attribute
that will be under the control of this function. We will call this
attribute ``output_rate_previous``, so to minimise the risk of a name
collision. If this attribute doesn't exist, then we must be on the first
day. We this save today's data so that it can be used tomorrow.

If the attribute does exist, then we can calculate a rate. We do that by
getting yesterday's data from ``output_rate_previous`` and todays data
from ``workspace.ward_inf_tot``. We then loop over all of the disease
stages, and open an output file for each stage (called ``rate_{i}.csv``).
We then write into this file the day, then the difference between today's
and yesterday's population in ward, for this ``ith`` disease stage.

Finally, we save today's data into ``workspace.output_rate_previous``, so
that it can be used tomorrow as yesterday's data.

Run this extractor in ``metawards`` using;

.. code-block:: bash

   metawards -d lurgy3 --extract extract_rate -a ExtraSeedsLondon.dat --nsteps 30

(again, we are limiting this to 30 steps just for demonstration reasons)

You should see that you have files ``rate_0.csv.bz2``, ``rate_1.csv.bz2`` etc.
now created in the output directory. If you look at these files you should
see that they contain the differences between the populations in each ward
for each disease stage between each day.

Writing to a database
---------------------

In the last example we wrote rates to a large number of files. While this has
worked, the data is beginning to get so large and multi-dimensional that
we are reaching the limits of what a CSV or other data file can
reasonably support. As data sizes get larger, it is better to start
writing data to a database.

The :class:`~metawards.OutputFiles` class has in-built support for opening
connections to `SQLite3 <https://docs.python.org/3/library/sqlite3.html>`__
databases. To use this, we call the function
:meth:`~metawards.OutputFiles.open_db`. For example, let's now create a new
extractor that will output the size of the population at each disease stage
for each day, and the change compared to the previous day. To do this,
open a file called ``extract_db.py`` and copy in the below;

.. code-block:: python

    from metawards.extractors import extract_default


    def create_tables(N_INF_CLASSES):
        # return a function that creates the tables
        # for the specified number of disease classes

        def initialise(conn):
            # create a table for the values...
            values = ",".join([f"stage_{i} int" for i in range(0, N_INF_CLASSES)])
            c = conn.cursor()
            c.execute(f"create table totals(day int, S int, {values})")

            # create a table for the rates...
            c.execute(f"create table deltas(day int, S int, {values})")

            conn.commit()

        return initialise


    def output_db(population, workspace, output_dir, **kwargs):
        have_yesterday = hasattr(workspace, "output_rate_previous")

        # get today's data
        inf_tot = workspace.inf_tot
        pinf_tot = workspace.pinf_tot
        S = population.susceptibles

        N_INF_CLASSES = workspace.n_inf_classes

        # open a database to hold the data - call the 'create_tables'
        # function on this database when it is first opened
        conn = output_dir.open_db("stages.db",
                                initialise=create_tables(N_INF_CLASSES))

        c = conn.cursor()

        # get the values for today
        today = [population.day, S] + [inf+pinf for inf, pinf in zip(inf_tot, pinf_tot)]

        # convert this to a string
        today_str = ",".join([str(t) for t in today])

        # write these to the database
        c.execute(f"insert into totals VALUES ({today_str})")

        if hasattr(workspace, "output_rate_db"):
            yesterday = workspace.output_rate_db

            # calculate the difference in all columns of today and yesterday
            deltas = [t - y for t, y in zip(today, yesterday)]
            # (except for the day, which should be today)
            deltas[0] = today[0]

            delta_str = ",".join([str(d) for d in deltas])

            # write this to the database
            c.execute(f"insert into deltas values ({delta_str})")

        conn.commit()

        # save today's data so that it can be used tomorrow
        workspace.output_rate_db = today


    def extract_db(**kwargs):
        funcs = extract_default(**kwargs)
        funcs.append(output_db)
        return funcs

Here, we have created a new function called ``create_tables`` that is called
to create a function that is returned and passed to
:meth:`~metawards.OutputFiles.open_db`. This function creates two tables
in the database; ``totals`` which contains the total population at each
disease stage, and ``deltas``, which contains the difference from the
previous day.

Next, we have ``output_db``. This function calls
:meth:`~metawards.OutputFiles.open_db` to create the connection, ``conn``
to the SQLite3 database. This connection is a standard Python
`sqlite3 connection object <https://docs.python.org/3/library/sqlite3.html>`__.

We calculate the total population in each stage as the sum of
``inf_tot`` (the workers at each stage) and ``pinf_tot`` (the players
at each stage). We prepend the number of susceptibles and also the
day number.

We then write this, as today's data, to the database via a cursor.

Next, we check if there is any data from yesterday by looking for the
custom attribute ``workspace.output_rate_db``. If there is, then we
get this data, and then calculate the difference from the previous day.
This is then written to the ``deltas`` table via the cursor.

Finally, we commit the changes to the database, and then save today's
data to ``workspace.output_rate_db`` so that it can be used tomorrow.

Run this extractor by typing;

.. code-block:: bash

   metawards -d lurgy3 --extract extract_db -a ExtraSeedsLondon.dat --nsteps 30

(again, we limit to 30 days just for illustration purposes)

Once this has finished, you should see a file called ``output/stages.db.bz2``.

Uncompress this file and then examine it using any SQLite3 database viewer,
e.g.

.. code-block:: bash

    # sqlite3 output/stages.db
    SQLite version 3.31.1 2020-01-27 19:55:54
    Enter ".help" for usage hints.
    sqlite> .dump
    PRAGMA foreign_keys=OFF;
    BEGIN TRANSACTION;
    CREATE TABLE totals(day int, S int, stage_0 int,stage_1 int,stage_2 int,stage_3 int,stage_4 int,stage_5 int);
    INSERT INTO totals VALUES(0,56082077,0,0,0,0,0,0);
    INSERT INTO totals VALUES(1,56082072,0,5,0,0,0,0);
    INSERT INTO totals VALUES(2,56082072,0,0,5,0,0,0);
    INSERT INTO totals VALUES(3,56082067,5,0,3,2,0,0);
    INSERT INTO totals VALUES(4,56082064,3,5,2,1,2,0);
    INSERT INTO totals VALUES(5,56082061,3,3,6,1,2,1);
    INSERT INTO totals VALUES(6,56082051,10,3,8,1,3,1);
    INSERT INTO totals VALUES(7,56082039,12,10,9,3,2,2);
    INSERT INTO totals VALUES(8,56082026,13,12,17,4,2,3);
    INSERT INTO totals VALUES(9,56082002,24,13,26,4,4,4);
    INSERT INTO totals VALUES(10,56081982,20,24,34,8,3,6);
    INSERT INTO totals VALUES(11,56081950,32,20,55,5,8,7);
    INSERT INTO totals VALUES(12,56081887,63,32,66,10,6,13);
    INSERT INTO totals VALUES(13,56081827,60,63,76,30,6,15);
    INSERT INTO totals VALUES(14,56081733,94,60,128,25,20,17);
    INSERT INTO totals VALUES(15,56081588,145,94,169,31,21,29);
    INSERT INTO totals VALUES(16,56081405,183,145,222,58,22,42);
    INSERT INTO totals VALUES(17,56081131,274,183,317,78,42,52);
    INSERT INTO totals VALUES(18,56080808,323,274,434,105,57,76);
    INSERT INTO totals VALUES(19,56080318,490,323,618,144,82,102);
    INSERT INTO totals VALUES(20,56079687,631,490,798,229,96,146);
    INSERT INTO totals VALUES(21,56078886,801,631,1134,273,159,193);
    INSERT INTO totals VALUES(22,56077748,1138,801,1539,356,224,271);
    INSERT INTO totals VALUES(23,56076190,1558,1138,2028,488,301,374);
    INSERT INTO totals VALUES(24,56074172,2018,1558,2759,652,388,530);
    INSERT INTO totals VALUES(25,56071530,2642,2018,3760,883,525,719);
    INSERT INTO totals VALUES(26,56067981,3549,2642,4997,1226,714,968);
    INSERT INTO totals VALUES(27,56063321,4660,3549,6723,1522,961,1341);
    INSERT INTO totals VALUES(28,56057020,6301,4660,8943,2095,1221,1837);
    INSERT INTO totals VALUES(29,56048552,8468,6301,11808,2834,1631,2483);
    CREATE TABLE deltas(day int, S int, stage_0 int,stage_1 int,stage_2 int,stage_3 int,stage_4 int,stage_5 int);
    INSERT INTO deltas VALUES(1,-5,0,5,0,0,0,0);
    INSERT INTO deltas VALUES(2,0,0,-5,5,0,0,0);
    INSERT INTO deltas VALUES(3,-5,5,0,-2,2,0,0);
    INSERT INTO deltas VALUES(4,-3,-2,5,-1,-1,2,0);
    INSERT INTO deltas VALUES(5,-3,0,-2,4,0,0,1);
    INSERT INTO deltas VALUES(6,-10,7,0,2,0,1,0);
    INSERT INTO deltas VALUES(7,-12,2,7,1,2,-1,1);
    INSERT INTO deltas VALUES(8,-13,1,2,8,1,0,1);
    INSERT INTO deltas VALUES(9,-24,11,1,9,0,2,1);
    INSERT INTO deltas VALUES(10,-20,-4,11,8,4,-1,2);
    INSERT INTO deltas VALUES(11,-32,12,-4,21,-3,5,1);
    INSERT INTO deltas VALUES(12,-63,31,12,11,5,-2,6);
    INSERT INTO deltas VALUES(13,-60,-3,31,10,20,0,2);
    INSERT INTO deltas VALUES(14,-94,34,-3,52,-5,14,2);
    INSERT INTO deltas VALUES(15,-145,51,34,41,6,1,12);
    INSERT INTO deltas VALUES(16,-183,38,51,53,27,1,13);
    INSERT INTO deltas VALUES(17,-274,91,38,95,20,20,10);
    INSERT INTO deltas VALUES(18,-323,49,91,117,27,15,24);
    INSERT INTO deltas VALUES(19,-490,167,49,184,39,25,26);
    INSERT INTO deltas VALUES(20,-631,141,167,180,85,14,44);
    INSERT INTO deltas VALUES(21,-801,170,141,336,44,63,47);
    INSERT INTO deltas VALUES(22,-1138,337,170,405,83,65,78);
    INSERT INTO deltas VALUES(23,-1558,420,337,489,132,77,103);
    INSERT INTO deltas VALUES(24,-2018,460,420,731,164,87,156);
    INSERT INTO deltas VALUES(25,-2642,624,460,1001,231,137,189);
    INSERT INTO deltas VALUES(26,-3549,907,624,1237,343,189,249);
    INSERT INTO deltas VALUES(27,-4660,1111,907,1726,296,247,373);
    INSERT INTO deltas VALUES(28,-6301,1641,1111,2220,573,260,496);
    INSERT INTO deltas VALUES(29,-8468,2167,1641,2865,739,410,646);
    COMMIT;
==========================
MetaWards.utils API Detail
==========================

.. automodule:: metawards.utils

.. toctree::
   :maxdepth: 1
================
MetaWards.themes
================

These are the classes that implement the different themes used
to colour the console output (and style the spinner).

.. toctree::
   :maxdepth: 1

   index_api_MetaWards_themes
===================
MetaWards.iterators
===================

These are the iterator functions that are used by the model to advance
the outbreak from day to day. The functions divide into two main types,
the most important of which are described below;

``iterate_functions`` which control which functions are called to advance the outbreak
    * :meth:`~metawards.iterators.iterate_default`
    * :meth:`~metawards.iterators.iterate_custom`
    * :meth:`~metawards.iterators.iterate_weekday`
    * :meth:`~metawards.iterators.iterate_weekend`
    * :meth:`~metawards.iterators.iterate_working_week`

``advance_functions`` which perform the actual work of advancing the outbreak
    * :meth:`~metawards.iterators.advance_foi`
    * :meth:`~metawards.iterators.advance_additional`
    * :meth:`~metawards.iterators.advance_infprob`
    * :meth:`~metawards.iterators.advance_fixed`
    * :meth:`~metawards.iterators.advance_play`
    * :meth:`~metawards.iterators.advance_work`

All of the above functions (and the many others in metawards.iterators) are
described :doc:`in more detail here <./index_api_MetaWards_iterators>`;

.. toctree::
   :maxdepth: 1

   index_api_MetaWards_iterators
===========================
MetaWards.movers API Detail
===========================

.. automodule:: metawards.movers

.. toctree::
   :maxdepth: 1
===============
MetaWards.utils
===============

These are utility functions that are used by the top-level MetaWards
package to run a model. The functions are heavily inspired by the
C functions from the original program.

The functions divide into four main types, the key functions of which
are described below;

Setting up the network
    * :meth:`~metawards.utils.build_wards_network`
    * :meth:`~metawards.utils.add_wards_network_distance`
    * :meth:`~metawards.utils.initialise_infections`
    * :meth:`~metawards.utils.initialise_play_infections`
    * :meth:`~metawards.utils.move_population_from_play_to_work`
    * :meth:`~metawards.utils.move_population_from_work_to_play`

Performing a model run
    * :meth:`~metawards.utils.run_model`
    * :meth:`~metawards.utils.iterate`
    * :meth:`~metawards.utils.iterate_weekend`

Extracting data from each iteration
    * :meth:`~metawards.utils.extract`

Performing multiple model runs in parallel
    * :meth:`~metawards.utils.run_models`
    * :meth:`~metawards.utils.run_worker`

All of the above functions (and the many others in metawards.utils) are
described :doc:`in more detail here <./index_api_MetaWards_utils>`;

.. toctree::
   :maxdepth: 1

   index_api_MetaWards_utils
===============================
MetaWards.extractors API Detail
===============================

.. automodule:: metawards.extractors

.. toctree::
   :maxdepth: 1
====================
MetaWards API Detail
====================

.. automodule:: metawards

.. toctree::
   :maxdepth: 1
=============================
MetaWards.analysis API Detail
=============================

.. automodule:: metawards.analysis

.. toctree::
   :maxdepth: 1
==================
MetaWards.analysis
==================

These are analysis functions that make it easy to produce graphs
or perform analyses. These use `pandas <https://pandas.pydata.org>`__
and `matplotlib <https://matplotlib.org>`__, which must both be
installed in able to use the `MetaWards.analysis` functions.

These functions are designed to either be called from within a
`Jupyter notebook <https://jupyter.org>`__, or are used by some of
the ``metawards`` command line tools to create quick outputs. An
example of such a tool is :doc:`metawards-plot <../metawards_plot_help>`.

Core functions include;

* :meth:`~metawards.analysis.save_summary_plots`
    Quick function that creates and saves summary plots from the
    ``results.csv.bz2`` files that are produced by ``metawards``

* :meth:`~metawards.analysis.import_graphics_modules`
    Safely import the graphics modules, and give a useful error message
    if this doesn't work. `metawards.analysis` needs Pandas and matplotlib,
    but these are often not available or are broken on clusters. We safely
    import these modules in this function to prevent any other part
    of ``metawards`` from being affected by a bad pandas or matplotlib
    install.

* :meth:`~metawards.analysis.create_average_plot`
    Create the average plot and return the resulting matplotlib figure.
    Use this to get and view figures in, e.g. a Jupyter notebook

* :meth:`~metawards.analysis.create_overview_plot`
    Create the overview plot and return the resulting matplotlib figure.
    Use this to get and view figures in, e.g. a Jupyter notebook

.. toctree::
   :maxdepth: 1

   index_api_MetaWards_analysis
=============
Documentation
=============

The metawards program and associated Python module is composed
of six main units;

* The :doc:`metawards <index_MetaWards_app>` command line program, which is
  what most people would use to run metawards. This program comes with
  full help that is provided by typing ``metawards --help``.

* The :doc:`metawards <index_MetaWards>` top-level Python module, which is
  what you should use if you want to write new Python programs that
  integrate metawards functionality. These could include plugging
  metawards as a model into a higher-level statictical
  analysis package.

* The :doc:`metawards.iterators <index_MetaWards_iterators>` module,
  which contains
  all of the package-supplied iterators that are used to advance the
  model infection from day to day. These are used to customise exactly
  how an infection progresses, and how control measures are applied.
  This is described in the :doc:`tutorial <../tutorial/index_part03>`.

* The :doc:`metawards.extractors <index_MetaWards_extractors>` module,
  which contains
  all of the package-supplied extractors that are used to extract and
  output data gathered live during a model run and write them to files.
  These are used to customise exactly what information is gathered,
  how it is processed, and where it is written as the infection progresses.
  You can see how to use these to customise output in
  the :doc:`tutorial <../tutorial/index_part04>`.

* The :doc:`metawards.mixers <index_MetaWards_mixers>` module, which
  contains all of
  the package-supplied mixers that are used to merge and mix values
  calculated across multiple demographic sub-networks. These are used
  to customise exactly how different demographics interact and how
  the disease should move between demographics. You can see how to use
  mixers in the :doc:`tutorial <../tutorial/index_part05>`.

* The :doc:`metawards.movers <index_MetaWards_movers>` module, which
  contains all of package-supplied movers that are used to move
  individuals between different demographics during an outbreak
  (e.g. move between home and hospital, work and holiday etc.).
  You can see how to use movers in the
  :doc:`tutorial <../tutorial/index_part06>`.

* The :doc:`metawards.analysis <index_MetaWards_analysis>` contains
  analysis functions that can be used to process the results of a
  ``metawards`` run to generate insights and useful graphics.

* The :doc:`metawards.utils <index_MetaWards_utils>` utility Python module
  which contains lots of functions that are used by metawards to build and
  perform model runs. These functions are internal to metawards and are
  not designed to be used outside of this program.

* The :doc:`metawards.themes <index_MetaWards_themes>` module contains
  the themes that are used to style the console output and spinners
  used by MetaWards. Choose the theme using the ``--theme`` option,
  e.g. ``--theme default`` or ``--theme simple``.

.. toctree::
   :maxdepth: 1

   index_MetaWards_app
   index_MetaWards
   index_MetaWards_iterators
   index_MetaWards_extractors
   index_MetaWards_mixers
   index_MetaWards_movers
   index_MetaWards_analysis
   index_MetaWards_utils
   index_MetaWards_themes
================
MetaWards.mixers
================

These are the mixer functions that are used to mix and merge data
between demographics during a multi-demographic *model run*.
The functions divide into two main types;


``mix_functions``, which control which functions will be called to merge data
    * :meth:`~metawards.mixers.mix_default`
    * :meth:`~metawards.mixers.mix_custom`
    * :meth:`~metawards.mixers.mix_evenly`
    * :meth:`~metawards.mixers.mix_none`

``merge_functions``, which perform the actual work of merging data across demographics
    * :meth:`~metawards.mixers.merge_evenly`

All of the above functions (and the many others in
:mod:`metawards.mixers`) are
described :doc:`in more detail here <index_api_MetaWards_mixers>`;

.. toctree::
   :maxdepth: 1

   index_api_MetaWards_mixers
===========================
MetaWards.mixers API Detail
===========================

.. automodule:: metawards.mixers

.. toctree::
   :maxdepth: 1
===========================
MetaWards.themes API Detail
===========================

.. automodule:: metawards.themes

.. toctree::
   :maxdepth: 1
==========================
MetaWards.app API Detail
==========================

.. automodule:: metawards.app

.. toctree::
   :maxdepth: 1
==============================
MetaWards.iterators API Detail
==============================

.. automodule:: metawards.iterators

.. toctree::
   :maxdepth: 1
=============
MetaWards.app
=============

This module contains all of the functions needed to write the
``metawards`` command line application.

The ``metawards`` command line application will automatically run itself
in one of three modes;

main
    This is the default mode, and is used when you normally run metawards
    to perform calculations using a single computer. Key functions that are
    used in ``main`` mode are;

    * :meth:`~metawards.app.cli`: The main command line interface
    * :meth:`~metawards.app.parse_args`: Parse command line arguments

    This mode is used to run either single model runs, or to run multiple
    runs in parallel on the same computer via a
    `multiprocessing pool <https://docs.python.org/3.4/library/multiprocessing.html?highlight=process#using-a-pool-of-workers>`__.

worker
    This is the mode used when ``metawards`` detects that it should be
    running as a worker process as part of a simulation run on a large
    cluster. If this is the case, then the application will go to sleep
    until it is called upon as part of a parallel run on multiple
    computers via either a
    `scoop mapping pool <https://scoop.readthedocs.io/en/0.7/usage.html#mapping-api>`__,
    or an
    `mpi4py MPI pool <https://mpi4py.readthedocs.io/en/stable/mpi4py.futures.html#mpipoolexecutor>`__.

    The worker will automatically detect whether scoop or mpi4py is being
    used based on which of these modules have been loaded via the
    ``-m scoop`` or ``-m mpi4py`` command line injector.

    The key function used in ``worker`` mode is;

    * :meth:`~metawards.app.get_parallel_scheme` : Auto-detect scoop or mpi4py

supervisor
    This is the mode used when ``metawards`` detects that it is being
    run across multiple computer nodes in a cluster. This is detected
    automatically via environment variables from common queueing systems
    that supply a ``hostfile``, or by manually specifying a
    ``hostfile`` via the ``--hostfile`` command line argument.

    The ``metawards`` application will choose to parallelise over the
    cluster using `mpi4py <https://mpi4py.readthedocs.io>`__ or
    `scoop <https://scoop.readthedocs.io>`__, depending on what the
    user specifies, which modules are installed, or automatically
    choosing in the order scoop > mpi4py (this order is simply because
    MPI is more error-prone to use automatically).

    Key functions used in ``supervisor`` mode are;

    * :meth:`~metawards.app.get_hostfile`: Get the hostfile
    * :meth:`~metawards.app.get_cores_per_node`: Work out the number of cores per compute node
    * :meth:`~metawards.app.get_threads_per_task`: Work out the number of threads per model run
    * :meth:`~metawards.app.scoop_supervisor`: Run a scoop supervisor to start and manage jobs
    * :meth:`~metawards.app.mpi_supervisor`: Run a mpi4py supervisor to start and manage jobs

All of the above functions (and others in the metawards.app package) are
described :doc:`in more detail here <index_api_MetaWards_app>`;

.. toctree::
   :maxdepth: 1

   index_api_MetaWards_app
====================
MetaWards.extractors
====================

These are the extractor functions that are used during a model run
to extract data from the outbreak to process and write live to files.
The functions divide into two main types,
the most important of which are described below;

``extract_functions`` which determine which functions will be called to output data
    * :meth:`~metawards.extractors.extract_default`
    * :meth:`~metawards.extractors.extract_custom`
    * :meth:`~metawards.extractors.extract_large`
    * :meth:`~metawards.extractors.extract_small`
    * :meth:`~metawards.extractors.extract_none`

``output_functions`` which perform the actual work of outputting data
    * :meth:`~metawards.extractors.output_basic`
    * :meth:`~metawards.extractors.output_core`
    * :meth:`~metawards.extractors.output_dispersal`
    * :meth:`~metawards.extractors.output_incidence`
    * :meth:`~metawards.extractors.output_prevalence`

All of the above functions (and the many others in metawards.extractors) are
described :doc:`in more detail here <./index_api_MetaWards_extractors>`;

.. toctree::
   :maxdepth: 1

   index_api_MetaWards_extractors
=========
MetaWards
=========

This is the top-level Python package that provides the core objects to
build and run MetaWards models.

The package centers around a few core objects:

:class:`~metawards.Nodes` / :class:`~metawards.Node`
    These represent the individual electoral wards in
    which progress of the disease is tracked.

:class:`~metawards.Links` / :class:`~metawards.Link`
    These represent the connections between electoral wards,
    including how people commute between wards for work. These
    therefore provide the routes via which the disease
    can spread.

:class:`~metawards.Parameters`
    This is the holder of all model parameters that can describe,
    e.g. the proportion of day versus night, cutoff distance
    for transmission between wards etc.
    Example parameter sets are held in the
    `MetaWardsData <https://github.com/metawards/MetaWardsData>`__
    repository.

:class:`~metawards.Disease`
    This is the holder of disease parameters, used to change the model
    to represent different types of disease outbreaks. Different
    disease models are held in the
    `MetaWardsData <https://github.com/metawards/MetaWardsData>`__
    repository.

:class:`~metawards.InputFiles`
    This holds information about all of the input files that are
    used to build the network of wards and links. All of the
    input data is held in the
    `MetaWardsData <https://github.com/metawards/MetaWardsData>`__
    repository.

:class:`~metawards.Network`
    This is the complete network of :class:`~metawards.Nodes` and
    :class:`~metawards.Links` built using
    the data loaded from the :class:`~metawards.InputFiles`,
    set up to model the
    disease whose parameters are in a
    :class:`~metawards.Disease`, using model run
    parameters held in a :class:`~metawards.Parameters` object.
    A Network is self-contained, containing everything needed for a model run.
    The :meth:`~metawards.Network.run` function runs the model.

:class:`~metawards.VariableSet` / :class:`~metawards.VariableSets`
    The model contains adjustable variables,
    which must be adjusted so that the model can match observed
    real-life data. A :class:`~metawards.VariableSet` contains a set
    of variables to
    be adjusted, while :class:`~metawards.VariableSets` is a collection of such
    changes that should be explored over multiple model runs.

:class:`~metawards.Population` / :class:`~metawards.Populations`
    A model run will result in a trajectory of
    changes in populations of people who have progressed along
    different stages of the disease (e.g. from being
    *susceptible to infection* (S), to being
    *removed from the outbreak* (R) - either because they
    have recovered or have sadly perished). A
    :class:`~metawards.Population` holds
    the numbers in each state for a single day in the outbreak,
    while :class:`~metawards.Populations` holds the full trajectory.

:class:`~metawards.OutputFiles`
    Manages all of the output files that are produced by the program.
    This can be told to auto-compress (bzip2) all files as they
    are being written.

All of the above classes (and others in the top-level package) are
described :doc:`in more detail here <index_api_MetaWards>`;

.. toctree::
   :maxdepth: 1

   index_api_MetaWards
================
MetaWards.movers
================

These are the mover functions that are used to move individuals between
different demographics during a multi-demographics *model run*.

The functions divide into two main types;

``move_functions`` that determine which functions are called to move individuals
    * :meth:`~metawards.movers.move_default`
    * :meth:`~metawards.movers.move_custom`

``go_functions`` that perform the work of moving individuals between demographics
    * :meth:`~metawards.movers.go_isolate`

All of the above functions (and the many others in
:mod:`metawards.movers`) are
described :doc:`in more detail here <./index_api_MetaWards_movers>`;

.. toctree::
   :maxdepth: 1

   index_api_MetaWards_movers
