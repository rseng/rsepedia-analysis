# Welcome to fishRman [![DOI](https://joss.theoj.org/papers/10.21105/joss.03467/status.svg)](https://doi.org/10.21105/joss.03467)
`fishRman` is an Open-Source Shiny R Dashboard to easily query, download, analyse and visualise Global Fishing Watch data on fishing effort. 

## Useful links
You can use the dashboard and learn more about the code, the software, the data, and the people behind them via:

- The [**Handbook**](https://raw.githubusercontent.com/Shyentist/fish-r-man/main/www/doc/Handbook.pdf) constitutes `fishRman`'s official instructions for use. Regardless of the prior knowledge of the user, reading the **Handbook** is key to the correct usage of the software.
- Website of [**Open-Source for Marine and Ocean Sciences (OSMOS)**](https://osmos.xyz/), which is our research group. If you like our projects, and would like a more interactive role, consider joining our [**Discord server**](https://discord.com/invite/W2unKxKbp7), following us on [**Twitter**](https://twitter.com/osmos_xyz), or [**donating**](https://www.buymeacoffee.com/osmos).
- The **Dashboard** can be accessed here: **https://shyentist.shinyapps.io/fish-r-man/**

## How to quote
- Software: Buonomo, P., (2021). fishRman: A Shiny R Dashboard improving Global Fishing Watch data availability. *Journal of Open Source Software*, 6(66), 3467, **https://doi.org/10.21105/joss.03467**
- Data: Global Fishing Watch, (2021). **https://globalfishingwatch.org/**

## How to become a contributor
First things first, you might want to take a look at the Issues page. There, I listed some tasks that I reckon the project needs. If you don't find anything that you can do, you can still run the code and make up **your own** opinion on what is needed. From there, you can either file this new Issue and forget about it, or file it and handle it yourself. The pull requests should be made onto the `pull-requests` branch.

## Get in touch
If you are interested in contributing, planning, chatting, or if you have any issues, you can contact me at **pasqualebuonomo@hotmail.it**---
title: 'fishRman: A Shiny R Dashboard improving Global Fishing Watch data availability'

tags:
  - R
  - fisheries
  - marine biology
  - global fishing watch
  - AIS data
  - dashboard
  - shiny
  - spatial analysis
authors:
  - name: Pasquale Buonomo
    orcid: 0000-0002-1848-9313
    affiliation: 1
affiliations:
 - name: Open-Source for Marine and Ocean Sciences (OSMOS)
   index: 1
date: 08 October 2021
bibliography: paper.bib
---

# Summary
One of the burdens of fisheries scientists is the scarcity or lack of consistent, 
extensive data on the subject. When such data do exist, they are often only available:

- To universities or other research institutions;
- Through bureaucratic ordeals;
- For a fee.

This issue has been tackled by Global Fishing Watch, an independent, international, 
non-profit organization promoting ocean sustainability through greater transparency, 
visualizing, tracking and sharing data about global fishing activity for free [@GFW].

While the datasets are indeed publicly available, they are also rather large and quite 
difficult to manage, since they require proficiency in coding. In fact, at present, the 
most notable reading material instructing on the use of the datasets targets an audience 
who is proficient in the languages R [@ClavelleR], Python, JavaScript [@ClavelleGEE], or 
SQL [@Mayorga] to download, filter, summarise, and visualise the data.


# Statement of need
Life sciences will soon need a widespread integration of computational approaches to store, 
manage, analyse, and visualise datasets that are quickly growing in size and complexity [@carey]. 
This is rather concerning, given how, although the number of published papers reporting the 
use of the R statistical language [@R] increased fivefold from 2007 to 2018 in the field of ecology 
[@lai],  most life science majors do not offer basic programming courses [@mariano].

Designed with ease of use in mind, `fishRman` is intended for a public of researchers,
students, managers, and stakeholders in the fields of fisheries science, life sciences, 
and economics, with little to no proficiency in programming, data analysis, or both, who
intend to query, download, filter, analyse, and visualise Global Fishing Watch data. 

Users who can program in R may also benefit from the software to avoid writing lines of code 
for  what has already been implemented in the dashboard, in order to focus on other aspects 
of their research, or even customize the source code to better meet their specific needs.

Users with a deeper understanding of statistics and fisheries science, and with prior knowledge 
of the datasets, only need to get acquainted with the software, while users that are new to the 
field can easily learn what they need to know via `fishRman`’s official instructions for use, the 
Handbook. Regardless of the prior knowledge of the user, reading the Handbook, which is available 
in the software itself and in the [GitHub repository](https://github.com/Shyentist/fish-r-man), is 
key to the correct usage of the software.

The user-friendly interface [@shiny; @shinyBS; @shinyjs; @shinyWidgets] allows users to 
easily interact with the SQL query constructor, seamlessly building [@glue; @stringi; @countrycode]
and running queries [@bigrquery; @DBI]. In a few clicks, users are able to analyse retrieved 
data in several different ways, such as:

- visualising the top n-th percentile of the dataframe for any percentage [@viridis; @sf; @maps; @ggplot], 
a key passage in assessing how fishing effort overlaps fishing stocks, protected or restricted areas, or 
another country's jurisdiction. 
- calculating the fishing effort exerted by specific countries via certain geartypes, in precise areas. This
is vital in assessing who is fishing where, when, and how they are doing it, so that fisheries management plans can
address the right issues even at an international level with a clear understanding of each country's responsibilities.
- producing time series of fishing effort with a daily, monthly, or yearly frequency [@dplyr; @tidyverse], which
is indispensable when searching for patterns to compare to species' life-cycles, seafood prices over the years,
compliance to maritime and market laws, and overall consistency with data from third-parties.

# Acknowledgements
Part of this work was carried out with the financial support from [Open-Source for Marine and Ocean Sciences (OSMOS)](https://osmos.xyz/)


# References# Contributor Covenant Code of Conduct

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
reported by contacting the project team at pasqualebuonomo@hotmail.it. All
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
---
name: Bug report
about: Create a report to help us improve
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
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
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
