---
title: 'DearScholar: A mobile application to conduct qualitative and quantitative diary research'
tags:
  - data collection
  - diary studies
  - longitudinal research
authors:
  - name: Peter M. Kruyen
    orcid: 0000-0003-0109-1744
    affiliation: "1"
affiliations:
 - name: Assistant Professor, Institute for Management Research, Radboud University, the Netherlands
   index: 1
date: 6 November 2020
bibliography: paper.bib
---

# Background 
Increasingly, scholars in the social and behavioral sciences prefer longitudinal data over cross-sectional data to explore new research phenomena, test effects, and build theories. Indeed, top-tier journals in, among others, management, psychology, and organization science nowadays tend to reject studies in which no longitudinal data has been used to investigate causal relations [e.g., @Bono; @Jonge; @Rico] in particular when testing mediation effects [cf. @Kline]. 

Diary studies are one particular class of research methods in which “self-report instruments [are] used repeatedly …  to investigate social, psychological, and physiological processes, within everyday situations…” [@Bolger,578]. While there is a rich tradition of diary studies [cf. @Iida; @Ohly], these methods are used relatively infrequently to collect data compared to other methods to collect longitudinal data such as experiments, panel studies, and archival research. Given the administrative burden of diary studies for both respondents (i.e., research participants) and researchers, their relative unfamiliarity (or unpopularity) is understandable.

However, because of their potential to obtain a better understanding of both *between* and, especially *within* individual differences over time, scholars across disciplines call for more research in which diary methods are applied, such as in public administration [e.g.,@Bakker; @Grimmelikhuijsen], marketing [e.g., @Elliott], and health research [e.g., @Jones].

# DearScholar and other diary research apps
DearScholar is a hybrid, open-source smartphone application (app) that can be used on iOS devices (iPhones and iPads) and Android devices (basically all other smartphones and tablets) to conduct diary studies and, obviously, for other types of longitudinal research such as repeated-survey designs and log studies. DearScholar's aim is to facilitate the research process for both respondents and researchers. Researchers can specify the number of measurement occasions (i.e., measurement schedule), survey layout, and question format. Respondents only have to download the app from the App Store[^1] or Google Play[^2], fill out their assigned credentials, and start participating.

A limited number of alternative (commercial) diary research apps have been developed, including @Indeemo, @LifeData, @Open, @PIEL, @RedCap, and @Teamscope. Acknowledging their value, some of these apps target the researcher as primary respondent instead of research participants; some are rather expensive for (large-scale) projects and make it difficult to change diary tasks during the study period; some are especially designed to collect qualitative data or quantitative data for one particular type of platform only (often Android); some store data outside the European Union, which is problematic for European researchers; and—last but not least—most alternative apps are closed-source projects.

# Use cases
Currently, DearScholar is applied in a study by Glenn Houtgraaf MSc, Dr. Peter M. Kruyen and Prof. Dr. Sandra van Thiel to investigate work-related creativity in government organizations. The app is used to follow about 100 participants over a period of six months, asking them closed- and open-ended questions at bi-weekly measurement occasions to investigate creative processes. In 2021, Liesbeth Faas MSc, Dr. Peter M. Kruyen, and Prof. Dr. Sandra van Thiel will replicate this study in local care teams.

# Acknowledgement
DearScholar is developed within the context of the research program "The creative public servant: Observations, explanations and consequences" with project number 406.18.R8.028, financed by the Dutch Research Council (NWO). The author wants to express his gratitude to Prof. Dr. Sandra van Thiel for her encouragements; to Glenn Houtgraaf MSc, Liesbeth Faas MSc, the ICT Services (Radboud University, Nijmegen, the Netherlands), and both reviewers and the  editor at the Journal of Open Source Software for their advice, testing, and feedback; and last, but not least, all (pilot) respondents for their effort and feedback during the developmental process.

# Resources
Visit the project page on GitHub[^3] for all resources (e.g., the manual, source code, and guides).

[^1]: Link to App Store: https://apps.apple.com/us/app/dearscholar/id1483121589?ls=1
[^2]: Link to Google Play: https://play.google.com/store/apps/details?id=net.peterkruyen.dearscholar
[^3]: Link to GitHub: https://github.com/pmkruyen/dearscholar

# References
[![status](https://joss.theoj.org/papers/1896b88f26b987b9c7a07035751afd7b/status.svg)](https://joss.theoj.org/papers/1896b88f26b987b9c7a07035751afd7b)

[![DOI](https://zenodo.org/badge/263641327.svg)](https://zenodo.org/badge/latestdoi/263641327)

[![Build Status](https://travis-ci.com/pmkruyen/dearscholar.svg?branch=master)](https://travis-ci.com/pmkruyen/dearscholar)

# Overview
Because collecting longitudinal data becomes more important in academic research nowadays and because limitations with available tools, DearScholar has been developed.

DearScholar allows researchers to easily and orderly collect rich and diverse qualitative and quantitative data over short and long periods of time to answer research questions about inter- and intra-individual changes, developments, and processes.

Developed as hybrid app in Cordova (using html, javascript, and css), DearScholar can be used to collect data on iOS devices (iPhones and iPads), Android devices, and--in the future--in web browsers too.

This page provides general information about the app, including an overview of current features. Background information, a summary of the current features, and current research projects can be found on the [project page](https://peterkruyen.net/dearscholar.html). 

The [Wiki section](https://github.com/pmkruyen/dearscholar/wiki) provides installation instructions, an overview of the settings and options, (automatic) testing options, and other details.

## Interested in using DearScholar in your academic research project? 
* Contact the main author (p.m.kruyen@fm.ru.nl) for credentials (username and password); download the app on the App Store or Google Play; take your time to test the app; send us feedback, and discuss with us how to implement DearScholar in your project :rocket:.

* For iOS devices (iPhones and iPads), download the app on [the App Store](https://apps.apple.com/us/app/dearscholar/id1577072187);
* For Android devices, download the app on [Google Play](https://play.google.com/store/apps/details?id=net.peterkruyen.dearscholar).

* Currently, the official deployed version of DearScholar stores research data on a secure server in the Netherlands only. If you want to use your own server, contact the main author to discuss how to make this possible.

## Current research projects
DearScholar is used in the following projects:
* Diary study on public servants' creativity (2020 - 2021) by Glenn Houtgraaf MSc, dr. Peter M. Kruyen, and prof. dr. Sandra van Thiel.
* Diary study on creativity in local-care teams (2021 - 2022) by Liesbeth Faas MSc, dr. Peter M. Kruyen, and prof. dr. Sandra van Thiel.
* Diary study on nurses' work engagement (2022) by Renée Vermeulen MSc and Evelien van Leeuwen MSc.

## Interested in contributing to DearScholar?
Super cool. Head to the [Wiki section](https://github.com/pmkruyen/dearscholar/wiki) for the installation instructions, settings and options, (automatic) testing, and other details.  Clone the project, post issues, and commit updates :icecream:.

## Acknowledgement
DearScholar is developed within the context of the research program "The creative public servant: Observations, explanations and consequences" with project number 406.18.R8.028, financed by the Dutch Research Council (NWO). The author wants to express his gratitude to Prof. Dr. Sandra van Thiel for her encouragements; to Glenn Houtgraaf MSc, Liesbeth Faas MSc, the ICT Services (Radboud University, Nijmegen, the Netherlands), and both reviewers and the editor at the Journal of Open Source Software for their advice, testing, and feedback; and last, but not least, all (pilot) respondents for their effort and feedback during the developmental process.

## License
Copyright (c) 2022 P.M. Kruyen, Institute for Management Research, Radboud University, the Netherlands. This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License (GPL3) as published by the Free Software Foundation. Radboud University, hereby disclaims all copyright interest in the program “DearScholar” (which offers an app to conduct longitudinal qualitative and quantitative diary, log 
and survey research) written by P.M. Kruyen. Commercial license available, please contact Radboud Innovation, Radboud University, the Netherlands. Radboud Innovation – Technology Transfer Office.

## Citation
Kruyen, P. M., (2020). DearScholar: A mobile application to conduct qualitative and quantitative diary research. Journal of Open Source Software, 5(55), 2506, https://doi.org/10.21105/joss.02506

# Current features
## Let respondents install DearScholar on their mobile device
1) A respondent has to download DearScholar from the App Store (iPhones and iPads) or Google Play (Android devices);
2) When opening DearScholar for the first time (left screenshot below), a respondent is asked to: 
    * allow push notifications; 
    * fill out the username and password that (s)he has received from the researcher; 
    * choose a 4-digit PIN if the device does not support Touch ID or Face Recognition; and 
    * agree to the informed consent form.
3) When everything goes well, DearScholar pulls the required survey tables and settings from the server, and the respondent is directed to the measurement schedule (homepage, right screenshot below).

<p align="center">
  <kbd><img src=https://github.com/pmkruyen/dearscholar/blob/master/screenshots/1.png width="350"/></kbd>
  <kbd><img src=https://github.com/pmkruyen/dearscholar/blob/master/screenshots/2.png width="350"/></kbd>
</p>

## Let respondents answer questions
When logging in to DearScholar, a respondent is directed to the measurement schedule (homepage) with all measurement occasions (dates). Future measurement occasions—measurement occasions beyond the current date—are locked and marked with a 'closed lock' icon. 

When a respondent clicks on a particular measurement occasion in the measurement schedule, (s)he is directed to a survey screen that displays all survey modules for that measurement date. Each module can be opened by clicking on the designated icon (see the screenshot below).

<p align="center">
  <kbd><img src=https://github.com/pmkruyen/dearscholar/blob/master/screenshots/3.png width="350"/></kbd>
</p>

DearScholar supports an unlimited number of survey modules, unlimited number of question pages in each module, a specification of which questions to appear on which page, simple branching and skipping logic, and different types of questions (binary questions, open questions, multiple-choice items, and rating sliders) to collect both quantitative and qualitative data (see the screenshots below for two examples). Researchers can make certain questions mandatory. Respondents can type their answers in the answer boxes... or--because DearScholar works with speech recognition--dictate their answers.

<p align="center">
  <kbd><img src=https://github.com/pmkruyen/dearscholar/blob/master/screenshots/4.png width="350"/></kbd>
  <kbd><img src=https://github.com/pmkruyen/dearscholar/blob/master/screenshots/5.png width="350"/></kbd>
</p>

When a respondent has completed a module (i.e., for that module, all questions have been answered and the data has been successfully uploaded to the server) the module icon turns green on the survey page. If all mandatory modules have been completed, the link to the measurement occasion turns grey and is marked with a 'sun' icon on the homepage.

DearScholar also includes links to additional, optional survey modules which can be found on the homepage's menu. Respondents can start these additional modules in between measurement occasions to report their thoughts once they occur (cf. event-sampling).

Respondents’ answers are not only send to a server, but also saved in the DearScholar. Respondents can access their previous answers by clicking on completed measurement modules, facilitating respondents to keep track of, and reread their own responses.

## Push notifications
DearScholar can be used with Google's Firebase to get the push notifications working *or* alternatively, a private push notification server can be used (using Node.JS for example). Contact the main author for advice on setting up such a server.

## Additional features
DearScholar has the capacity to send short in-app messages to specific respondents. For example, respondents can be sent thank-you messages to show engagement, small encouragements when respondents have missed a measurement occasions, ask follow-up questions, or invited respondents to elaborate on particular answers over the phone or through email. These in-app messages are displayed on a separate message screen. An envelope icon appears in the app’s status bar when new messages have been sent. When respondents have read the message, the researcher is noted.
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
reported by contacting the project team at p.m.kruyen@fm.ru.nl. All
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
<a href="https://www.patreon.com/vladimirkharlampidi"><img src="https://framework7.io/i/support-badge.png" height="20"></a>

# Framework7

Full Featured Mobile HTML Framework For Building iOS & Android Apps

## Supporting Framework7

Framework7 is an MIT-licensed open source project with its ongoing development made possible entirely by the support of these awesome [backers](https://github.com/framework7io/framework7/blob/master/BACKERS.md). If you'd like to join them, please consider [becoming a backer or sponsor on Patreon.](https://www.patreon.com/vladimirkharlampidi)

<table>
  <tr>
    <td align="center" valign="middle">
      <a href="https://privicy.com/" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/privicy.png">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://codersrank.io/?utm_source=partner&utm_medium=referral&utm_campaign=vladimir" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/codersrank.svg">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://www.sparheld.de/" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/sparheld.jpg">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://www.thoriumbuilder.com/" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/thorium.png">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://app-valley.vip/" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/appvalley.jpg">
      </a>
    </td>
  </tr>
  <tr>
    <td align="center" valign="middle">
      <a href="http://mytommy.com" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/tommy.png">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://www.securcom.me/" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/securcom.png">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="http://ananyamultitech.com/" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/ananya.png">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://kqapi.us" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/kqapius.png">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://writersperhour.com/write-my-paper" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/writersperhour.jpg">
      </a>
    </td>
  </tr>
  <tr>
    <td align="center" valign="middle">
      <a href="https://monovm.com/linux-vps/" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/monovm.jpg">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://unicorn.io" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/unicorn.svg">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://www.colognewebdesign.de/" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/colognewebdesign.png">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://rise.co" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/rise.png">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://usave.co.uk/utilities/broadband" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/usave.png">
      </a>
    </td>
  </tr>
  <tr>
    <td align="center" valign="middle">
      <a href="https://bid4papers.com/write-my-essay.html" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/bid4papers.png">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://kidoverse.app" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/kidoverse.png">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://www.cyberbrain.nl/" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/cyberbrain.png">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://hicapps.cl/web/" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/hicapps.png">
      </a>
    </td>
    <td align="center" valign="middle">
      <a href="https://blokt.com/" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/blokt.png">
      </a>
    </td>
  </tr>
  <tr>
    <td align="center" valign="middle">
      <a href="https://wappler.io/" target="_blank">
        <img width="160" src="https://framework7.io/i/sponsors/wappler.png">
      </a>
    </td>
    <td></td>
    <td></td>
    <td></td>
  </tr>
</table>

## Getting Started
  * [Getting Started Guide](https://framework7.io/docs/introduction.html)
  * [Installation Guide](https://framework7.io/docs/installation.html)
  * [App Layout](https://framework7.io/docs/app-layout.html)
  * [Initialize App](https://framework7.io/docs/init-app.html)

## Forum

If you have questions about Framework7 or want to help others you are welcome to special forum at https://forum.framework7.io/

## Docs

Documentation available at https://framework7.io/docs/

## Tutorials

Tutorials available at https://framework7.io/tutorials/

## Showcase

Appstore apps made with Framework7: https://framework7.io/showcase/
