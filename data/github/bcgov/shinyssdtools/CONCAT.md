
<!-- README.md is generated from README.Rmd. Please edit that file -->

# shinyssdtools

<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/Lifecycle-Maturing-007EC6)](https://github.com/bcgov/repomountie/blob/master/doc/lifecycle-badges.md)
[![R-CMD-check](https://github.com/bcgov/shinyssdtools/workflows/R-CMD-check/badge.svg)](https://github.com/bcgov/shinyssdtools/actions)
[![Apache
license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02848/status.svg)](https://doi.org/10.21105/joss.02848)
<!-- badges: end -->

`shinyssdtools` is a shiny application for fitting Species Sensitivity
Distributions (SSDs) using
[`ssdtools`](https://github.com/bcgov/ssdtools).

## Utilization

The most recent version is available online at
<https://bcgov-env.shinyapps.io/ssdtools/>.

To install the development version from
[GitHub](https://github.com/bcgov/shinyssdtools) and deploy locally use

``` r
# install.packages("devtools")
devtools::install_github("bcgov/shinyssdtools")
library(shinyssdtools)
shinyssdtools::run_app()
```

## Features

In addition to being a Graphical User Interface to the core
functionality in the [`ssdtools`](https://github.com/bcgov/ssdtools)
package, `shinyssdtools` also provides

-   a bilingual (English/French) interface;
-   generation of R scripts for reproducibility;
-   customization and downloads of plots and tables

## Information

For more information including how to cite `shinyssdtools` see [Dalgarno
(2021)](https://doi.org/10.21105/joss.02848).

For a review of `ssdtools` and `shinysddtools` in the context of other
SSD software packages see [Fox et
al. (2021)](https://onlinelibrary.wiley.com/doi/10.1002/etc.4925).

## Assistance

To report bugs/issues/feature requests, please file an
[issue](https://github.com/bcgov/shinyssdtools/issues/).

## Contribution

If you would like to contribute, please see our
[CONTRIBUTING](CONTRIBUTING.md) guidelines.

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/bcgov/shinyssdtools/blob/master/.github/CONTRIBUTING.md).
By participating in this project you agree to abide by its terms.

## Deploying to shinyapps.io

### Manually

Run the `deploy-app.R` script in the scripts directory (after setting
the account argument to be your `shinyapps.io` account name).

### Automatically

If your `shinyapps.io` account name is the same as your GitHub account
name simply make a commit in the master or dev branch and include
`deploy app` in the message (after setting `SHINYAPPS_TOKEN` and
`SHINYAPPS_SECRET` in your repository GitHub secrets). This triggers the
`deploy-app.yml` GitHub action.

## License

The code is released under the Apache License 2.0

Copyright 2021 Province of British Columbia

Licensed under the Apache License, Version 2.0 (the “LICENSE”); you may
not use this file except in compliance with the License. You may obtain
a copy of the License at

<http://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an “AS IS” BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

------------------------------------------------------------------------

<a rel="LICENSE" href="https://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence"
style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a><br /><span
xmlns:dct="http://purl.org/dc/terms/"
property="dct:title">ssdtools</span> by <span
xmlns:cc="http://creativecommons.org/ns#"
property="cc:attributionName">the Province of British Columbia </span>
is licensed under a
<a rel="LICENSE" href="https://creativecommons.org/licenses/by/4.0/">
Creative Commons Attribution 4.0 International License</a>.
# shinyssdtools 0.0.3.9001 (2021-04-03)

- Add `run_app()` alternative to `run_ssdtools_app()`.


# shinyssdtools 0.0.3.9000 (2021-04-02)

- Internal changes only.


---
title: 'shinyssdtools: A web application for fitting Species Sensitivity Distributions (SSDs)'
date: '2020-11-18'
authors:
- affiliation: 1
  name: Seb Dalgarno
  orcid: 0000-0002-3658-4517
affiliations:
- name: Poisson Consulting Ltd., Nelson, British Columbia
  index: 1
bibliography: paper.bib
tags:
   - R
   - Shiny
   - ssdtools
   - species sensitivity distributions
---

# Summary
The species sensitivity distribution (SSD) is the most widely used method for getting water quality benchmarks to characterize effects of chemical contaminants for water quality or ecological risk assessment [@fox_2020].
This typically involves estimating the concentration of a chemical that affects 5% of the species considered [@posthuma_species_2001].
The `ssdtools` R package [@ssdtools] has recently advanced SSD methods by providing model averaging using information-theoretic criteria and the construction of confidence intervals using bootstrapping [@fox_2020]. 

`shinyssdtools` is a web-based graphical user interface (GUI) to the `ssdtools` R package.
`shinyssdtools`, which was developed in the Shiny web framework [@shiny], is an R package in its own right.
As well as providing access to the core functionality in `ssdtools`, it also offers the following value-added features: a bilingual (English/French) interface, customization and downloads of plots and tables and generation of R scripts for reproducibility. 
`shinyssdtools` can be accessed from the web or it can be run locally by installing the R package.

### Graphical User Interface

The `shinyssdtools` web application has six navigational tabs: 

1. Data
   - Upload a dataset or enter data manually.
1. Fit
   - Select distributions.
   - Calculate information-theoretic criteria.
1. Predict
   - Estimate the concentration that affects a specific percentage of the species.
   - Calculate confidence limits using bootstrapping.
1. R code
   - Copy the R code required to reproduce the results.
1. About
   - version information, explanation of abbreviations and references
1. User guide
   - Step-by-step guide to proper use of the application.

![shinyssdtools user interface](shinyssdtools_ui.png)

# Installation

The `shinyssdtools` application is available at https://bcgov-env.shinyapps.io/ssdtools/.
`shinyssdtools` is bundled as an R package [@r] to allow the user to install and run locally using just three lines of R code:

```r
install.packages('remotes')
remotes::install_github('bcgov/shinyssdtools')
shinyssdtools::run_ssdtools_app()
```

# Statement of Need
`ssdtools` is an R package that has provided recent advances in SSD methodology [@fox_2020]. 
`shinyssdtools` provides access to this functionality via a simple, modern GUI without requiring users to be familiar with the R programming language. 
Data can be easily uploaded to the application, the interface can be viewed in multiple languages (French and English) and R code output is provided so that analyses can be reproduced or shared. 
`shinyssdtools` was initially developed for the Province of British Columbia with input from the governments of Canada and Australia and has been used by the governments of British Columbia and Canada to derive water quality benchmarks.

# Contribution
 
All historical and existing SSD software was recently reviewed by @fox_2020 who considered `ssdtools`, `shinyssdtools` and `SSD Toolbox` to be the most important contributions because they provide model averaging.

`SSD Toolbox` is standalone software performing similar functionality to `shinyssdtools` developed by the US Environmental Protection Agency [@ssdtoolbox_2020]. 
It can be downloaded as a Windows executable file and requires installation of version 9.5 of the MATLAB® Runtime Compiler (MCR) from Mathworks, which requires 3.75 GB of hard disk space. 
`shinyssdtools` provides a more appealing user interface, is open-source and does not require local installation of bulky software. 

The similarly named `shinyssd` is an alternative open source Shiny web application to fit SSDs that is also bundled as an R package [@dandrea_shinyssd_2019]. 
`shinyssdtools` contributes by being bilingual; providing additional distributions including the gamma, Gompertz and log-Gumbel; by allowing the user to model average and by providing R scripts to replicate the analysis.

# Acknowledgements

We acknowledge contributions from Angeline Tillmanns, Marianne Métivier, Andy Teucher, David Fox, Carl Schwarz and Joe Thorley.
Development of `shinyssdtools` was initially funded by the Ministry of Environment and Climate Change Strategy, British Columbia. The governments of British Columbia, Canada and Australia have also contributed to its development.

## References
# Getting help with shinyssdtools

Thanks for using shinyssdtools!
Before filing an issue, there are a few places to explore and pieces to put together to make the process as smooth as possible.

## Make sure its new

Before opening a new issue, be sure to [search issues and pull requests](https://github.com//issues) to make sure the bug hasn't been reported and/or already fixed in the development version. 
By default, the search will be pre-populated with `is:issue is:open`. 
You can [edit the qualifiers](https://help.github.com/articles/searching-issues-and-pull-requests/)  (e.g. `is:pr`, `is:closed`) as needed. 
For example, you'd simply remove `is:open` to search _all_ issues in the repo, open or closed.

## Make a reprex

Start by making a minimal **repr**oducible **ex**ample using the  [reprex](https://reprex.tidyverse.org/) package. 
If you haven't heard of or used reprex before, you're in for a treat! 
Seriously, reprex will make all of your R-question-asking endeavors easier (which is a pretty insane ROI for the five to ten minutes it'll take you to learn what it's all about). 
For additional reprex pointers, check out the [Get help!](https://www.tidyverse.org/help/) section of the tidyverse site.
Thank you for taking the time to submit a pull request!

To maximize the chances of acceptance:

* The title of your PR should briefly describe the change.

* The body of your PR should contain `Fixes #issue-number` (if relevant).

* Commit/merge messages to be included in NEWS.md should begin with `-`.

* Code should follow the tidyverse [style guide](https://style.tidyverse.org).

* Documentation should use roxygen2, with Markdown syntax.

* Contributions should include unit tests (using `testthat`).

For more information see [Contributing](/.github/CONTRIBUTING.md).
# Contributor Code of Conduct

As contributors and maintainers of this project, and in the interest of
fostering an open and welcoming community, we pledge to respect all people who
contribute through reporting issues, posting feature requests, updating
documentation, submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free
experience for everyone, regardless of level of experience, gender, gender
identity and expression, sexual orientation, disability, personal appearance,
body size, race, ethnicity, age, religion, or nationality.

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery
* Personal attacks
* Trolling or insulting/derogatory comments
* Public or private harassment
* Publishing other's private information, such as physical or electronic
  addresses, without explicit permission
* Other unethical or unprofessional conduct

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

By adopting this Code of Conduct, project maintainers commit themselves to
fairly and consistently applying these principles to every aspect of managing
this project. Project maintainers who do not follow or enforce the Code of
Conduct may be permanently removed from the project team.

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community.

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting a project maintainer at ali.azizishirazi@gov.bc.ca. All complaints will be reviewed and investigated 
and will result in a response that is deemed necessary and appropriate to the 
circumstances. Maintainers are obligated to maintain confidentiality with regard 
to the reporter of an incident.


This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 1.3.0, available at
[http://contributor-covenant.org/version/1/3/0/][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/3/0/

---
*This project was created using the [bcgovr](https://github.com/bcgov/bcgovr) package.*
## How to contribute
Government employees, public and members of the private sector are encouraged to contribute to the repository by **forking and submitting a pull request**. 

(If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git) and  check out a more detailed guide to [pull requests](https://help.github.com/articles/using-pull-requests/).)

Pull requests will be evaluated by the repository guardians on a schedule and if deemed beneficial will be committed to the master.

All contributors retain the original copyright to their stuff, but by contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users **under the terms of the license under which this project is distributed.**

---
*This project was created using the [bcgovr](https://github.com/bcgov/bcgovr) package.*
---
name: Feature request
about: Suggest an idea for this project
---
  
## Your Idea

Please briefly describe your idea.
---
name: Bug report
about: Describe a bug you're experiencing
---

## The Bug

Please briefly describe your problem and what output you expect.

## A Reprex

Please include a minimal reproducible example (AKA a reprex). 
If you've never heard of a [reprex](http://reprex.tidyverse.org/) before, start by reading <https://www.tidyverse.org/help/#reprex>.

```r
# insert reprex here
```
Cette application web ajuste les fonctions de distribution de la sensibilité des espèces aux données de concentration. L'application est basée sur le progiciel R ssdtools, et partage les mêmes fonctions. Les changements et les mises à niveau apportés à ssdtools vont résulter en modifications de cette application. Il est recommandé, lorsque l'estimation du HC5 généré par cette application est rapportée, d'inclure aussi la version de ssdtools et le nom des fonctions de distribution ajustées aux données.

Les colonnes du tableau de l'évaluation de la qualité de l’ajustement des courbes de distributiont sont la distribution (dist), la statistique d’Anderson-Darling (ad), la statistique de Kolmogorov-Smirnov (ks), la statistique de Cramer-von-Mises (cmv), le critère d’information Akaike (aic), le critère d’information Akaike corrigé pour la taille de l’échantillon (aicc), le critère d’information Bayésien (bic), la différence entre AICc (delta) et la pondération des critères d'information AICc (coefficient de pondération). L’estimation de la fonction de distribution finale est basée sur l’inférence multimodèle (à partir de l’AICc). La concentration présentant un risque est la concentration estimée d’une substance affectant un centile (seuil) sélectionné de l’ensemble des espèces.

Pour citer l’application R ‘ssdtools’:
Thorley, J. and Schwarz C., (2018). ssdtools: An R package to fit Species Sensitivity Distributions. Journal of Open Source Software, 3(31), 1082. https://doi.org/10.21105/joss.01082

Pour citer l’application web :
Seb Dalgarno (2018) ssdtools: A shiny web app to analyse species sensitivity distributions. Prepared by Poisson Consulting for the Ministry of the Environment, British Columbia. https://bcgov-env.shinyapps.io/ssdtools/

Pour plus d'information sur l'utilisation de l'inférence multimodèle afin d'obtenir des estimations de HC5, veuillez consulter:
[Schwarz, C.J. and A.R. Tillmanns. 2019. Improving statistical methods to derive species sensitivity distributions. Water Science Series, WSS2019-07, Province of British Columbia, Victoria.](http://a100.gov.bc.ca/appsdata/acat/documents/r57400/2_1568399094009_8398900200.pdf)This app **fits species sensitivity distributions to concentration data**. The app is built from the R package [ssdtools](https://github.com/bcgov/ssdtools), and shares the same functionality.


*Hint: Find and click the info icons  throughout the app to find more information on a particular input.*  

### Step 1: Provide data 

* Data should be provided for **only one chemical** at a time. 
* Each species should have only one concentration value. 
* Data must have **at least one column** containing **at least 8 distinct, positive, non-missing concentration values**. 
* Optionally, **species and group** columns can be included, which are used to label and color plot output, respectively.  
* Any additional columns are accepted but are not used by any functions.


<center>

Concentration&nbsp;&nbsp; | Species&nbsp;&nbsp; | Group &nbsp;
--- | --- | ---
2.1 | Oncorhynchus mykiss &nbsp; | Fish
2.4 | Ictalurus punctatus &nbsp;| Fish  
4.1 | Micropterus salmoides &nbsp;| Fish
10  | Brachydanio rerio &nbsp;| Fish
15.6 | Carassius auratus &nbsp;| Fish
18.3 | Pimephales promelas &nbsp;| Fish 
6 | Daphnia magna &nbsp;| Invertebrate
10 | Opercularia bimarginata &nbsp;| Invertebrate

</center>

There are three options to provide data to the app:  

1. **Use the demo Boron dataset**. 
    - Quickly preview the app functionality on a dataset that 'works'. 
    - Citation: [Canadian Council of Ministers of the Environment. 2009. Canadian water quality guidelines for the protection of aquatic life: Boron. In: Canadian  environmental  quality guidelines, 2009, Canadian Council of  Ministers of the Environment, Winnipeg.](http://ceqg-rcqe.ccme.ca/download/en/324/)
2. **Upload a csv file**. 
    - Excel file formats are not accepted. If you have an excel file, try exporting a worksheet to csv. 
3. **Fill out the interactive table**. 
    - Species and Group columns are optional. Click on a cell to begin entering data. Right-click on the table to delete/insert rows or columns. Column names cannot be changed. 
    
Finally, preview the data provided in the table on the right hand side of the tab.  

<center>
![Provide data in tab `1. Data`](https://media.giphy.com/media/fjyG7Wbgp1RKV8sS9s/giphy.gif)
</center>

### Step 2: Fit distributions 

1. Specify **which column contains concentration values**. The app attempts to guess which column contains concentration values based on data column names. This may need to be corrected.
2. **Select or deselect distributions to fit the data**.  Note that if two or more models have overlapping fits then support for this model shape will be over inflated in the model averaged parameters.  Please see the article [here](https://bcgov.github.io/ssdtools/articles/distributions.html) for more information.  The outputs may take a moment to update.
3. Format the plot using inputs in the sidebar and **download plot and goodness of fit table** as png and csv files, respectively.

<center>
![Fit distributions in tab `2. Fit`](https://media.giphy.com/media/yIjxGcFKt0zK2vj38j/giphy.gif)
</center>

Additional information about the **goodness of fit table**:
The columns in the goodness of fit table are the distribution (dist), the Anderson-Darling statistic (ad), the Kolmogorov-Smirnov statistic (ks), the Cramer-von Mises statistic (cvm), Akaike’s Information Criterion (aic), Akaike’s Information Criterion corrected for sample size (aicc), Bayesian Information Criterion (bic), the AICc difference (delta) and the AICc based Akaike weight (weight). The prediction is the model averaged (using aicc) estimate of the fit. The percent hazard concentration is the concentration of the chemical which is predicted to affect that percent of the species tested.

### Step 3: Predict hazard concentration or percent of species effected
1. Select the **threshold % species affected** to calculate **estimated hazard concentration** OR select **concentration** to calculate the percentage of species affected by a specified concentration. This affects the plot (dotted line), text displayed below the plot and calculations of confidence limits.  
2. Select the number of **bootstrap samples used to calculate confidence limits**. The recommended number of samples is 10,000, although this can take around 3 minutes to process. Select lower number of bootstrap samples to reduce processing time.  

<center>

Bootstrap Samples &nbsp;&nbsp; | Estimated Processing Time
--- | ---
10,000 &nbsp; | 45 seconds
5,000 &nbsp;| 20 seconds 
1,000 &nbsp;| 10 seconds
500 &nbsp;| 5 seconds

</center>

3. Since confidence limits take time to calculate, they are not calculated automatically; you must press the `Get CL` button.
4. **Format plot** using various inputs in sidebar and **download plot and table** as png and csv file, respectively.

<center>
![Get hazard concentration estimates and confidence limits in tab `3. Predict`](https://media.giphy.com/media/xKb9nQsPFlqTCGzUgF/giphy.gif)
</center>

### Step 4: Get R code

Copy R code to reproduce outputs programmatically. Code is dynamically generated based on user inputs and functions executed within the app (e.g., code for generating confidence limits will appear after 'Get CL' button is clicked). 

<center>
![Get R code to reproduce results programmatically in tab `R Code`](https://media.giphy.com/media/XIgsL03rRnEfn8nNas/giphy.gif)
</center>

To generate a graph with confidence bands, copy the R code and paste in R.  Then set ci = TRUE in the predict and ssd_plot functions.

 

This app fits species sensitivity distributions to concentration data. The app is built from the R package ssdtools, and shares the same functionality. Changes and upgrades in ssdtools will result in changes to this app. It is recommended that when reporting the HC5 estimates generated using this app that the version of ssdtools and the name of the distributions fit to the dataset are listed.

The columns in the goodness of fit table are the distribution (dist), the Anderson-Darling statistic (ad), the Kolmogorov-Smirnov statistic (ks), the Cramer-von Mises statistic (cvm), Akaike's Information Criterion (aic), Akaike's Information Criterion corrected for sample size (aicc), Bayesian Information Criterion (bic), the AICc difference (delta) and the AICc based Akaike weight (weight). The prediction is the model averaged (using aicc) estimate of the fit. The percent hazard concentration is the concentration of the chemical which is predicted to affect that percent of the species tested.

To cite package ssdtools in publications use:
Thorley, J. and Schwarz C., (2018). ssdtools: An R package to fit Species Sensitivity Distributions. Journal of Open Source Software, 3(31), 1082. https://doi.org/10.21105/joss.01082

To cite the web app use:
Seb Dalgarno (2018) ssdtools: A shiny web app to analyse species sensitivity distributions. Prepared by Poisson Consulting for the Ministry of the Environment, British Columbia. https://bcgov-env.shinyapps.io/ssdtools/

For more information on using model averaging to generate HC5 estimates, please see:
[Schwarz, C.J. and A.R. Tillmanns. 2019. Improving statistical methods to derive species sensitivity distributions. Water Science Series, WSS2019-07, Province of British Columbia, Victoria.](http://a100.gov.bc.ca/appsdata/acat/documents/r57400/2_1568399094009_8398900200.pdf)## Data Sources

The `CCME data.csv` data file is provided&mdash;with permission to use and redistribute&mdash;by the [Canadian Council of the Ministers of the Environment (CCME)](http://ceqg-rcqe.ccme.ca/en/index.html).

The citations and data sources are as follows:

- Boron: [Canadian Council of Ministers of the Environment. 2009. Canadian water quality guidelines for the protection of aquatic life: Boron. In: Canadian  environmental  quality guidelines, 2009, Canadian Council of  Ministers of the Environment, Winnipeg.](http://ceqg-rcqe.ccme.ca/download/en/324/)
- Cadmium: [Canadian  Council of Ministers  of  the  Environment.  2014.  Canadian  water  quality  guidelines  for  the  protection  of  aquatic  life: Cadmium. In: Canadian environmental quality guidelines, 1999, Canadian Council of Ministers of the Environment, Winnipeg.](http://ceqg-rcqe.ccme.ca/download/en/148/)
- Chloride: [Canadian Council of Ministers of the Environment. 2011. Canadian water quality guidelines for the protection of aquatic life: Chloride. In: Canadian environmental quality guidelines, 1999, Canadian Council of Ministers of the Environment, Winnipeg. ](http://ceqg-rcqe.ccme.ca/download/en/337/)
- Endosulfan: [Canadian Council of Ministers of the Environment. 2010. Canadian water quality guidelines for the protection of aquatic life:   Endosulfan.   In:   Canadian environmental quality guidelines, 1999, Canadian Council of Ministers of the Environment, Winnipeg.](http://ceqg-rcqe.ccme.ca/download/en/327/)
- Glyphosate: [Canadian Council of Ministers of the Environment. 2012. Canadian water quality guidelines for the protection of aquatic life:  Glyphosate.  In: Canadian  environmental quality guidelines, Canadian Council of Ministers of the Environment,   Winnipeg.](http://ceqg-rcqe.ccme.ca/download/en/182/)
- Uranium: [Canadian Council of Ministers of the Environment. 2011. Canadian water quality guidelines for the protection of aquatic life: Uranium. In: Canadian environmental quality guidelines, 1999, Canadian Council of Ministers of the Environment, Winnipeg.](http://ceqg-rcqe.ccme.ca/download/en/328/)
- Silver: [Canadian   Council   of   Ministers   of   the   Environment.   2015.   Canadian water quality guidelines for the protection of aquatic life: Silver. In: Canadian environmental quality guidelines, 1999, Canadian Council of Ministers of the Environment, Winnipeg.](http://ceqg-rcqe.ccme.ca/download/en/355/)
Cette application web **estime les fonctions de distribution de la sensibilité des espèces grâce aux données de concentration.** L’application est faite à partir de l’application R [ssdtools](https://github.com/bcgov/ssdtools), et partage les mêmes fonctions. 

*Astuce: Clic sur les icônes info pour obtenir plus d’informations.*

### Étape 1: Fournir les données

* Les données doivent être fournies pour **une substance chimique** à la fois.
* Une seule valeur doit être entrée par espèce. 
* L’ensemble des données doit contenir **au moins une colonne** comportant au **minimum 8 concentrations distinctes,positives et sans valeur manquante**.
* Des colonnes pour **les espèces et les groupes** taxonomiques peuvent également être ajoutées. Des étiquettes et des couleurs sont alors disponibles pour permettre leur identification dans le graphique. 
* Des colonnes additionnelles sont possibles mais elles ne sont reliées à aucune fonction. 

<center>

Concentration&nbsp;&nbsp; | Espèce&nbsp;&nbsp; | Groupe &nbsp;
--- | --- | ---
2.1 | Oncorhynchus mykiss &nbsp; | Poisson
2.4 | Ictalurus punctatus &nbsp;| Poisson 
4.1 | Micropterus salmoides &nbsp;| Poisson
10  | Brachydanio rerio &nbsp;| Poisson
15.6 | Carassius auratus &nbsp;| Poisson
18.3 | Pimephales promelas &nbsp;| Poisson
6 | Daphnia magna &nbsp;| Invertébré
10 | Opercularia bimarginata &nbsp;| Invertébré

</center>

Trois options sont disponibles pour déposer les données dans l’application:

1. **Utiliser les données démo du Bore.**
    - Permet de connaître rapidement les fonctions de l’application avec un ensemble «fonctionnel» de données.
    - Citation: [Canadian Council of Ministers of the Environment. 2009. Canadian water quality guidelines for the protection of aquatic life: Boron. In: Canadian  environmental  quality guidelines, 2009, Canadian Council of  Ministers of the Environment, Winnipeg.](http://ceqg-rcqe.ccme.ca/download/en/324/)
2. **Télécharger un fichier .csv**.
    - Le format Excel n’est pas accepté. Si vous avez un fichier Excel, il faut l’exporter dans une feuille de calcul .csv.
3. **Remplir le tableau interactif.**
    - Cliquer sur une cellule pour commencer à entrer les données. Faire un clic-droit pour ajouter ou supprimer des colonnes ou des rangées. Le nom des colonnes ne peut pas être modifiés. 

Enfin, visualiser les données fournies dans le tableau à droite de l'onglet.

<center>
![Provide data in tab `1. Data`](https://media.giphy.com/media/fjyG7Wbgp1RKV8sS9s/giphy.gif)
</center>

### Étape 2: Ajustement des distributions 

1. Spécifier **quelle est la colonne qui contient les valeurs de concentrations**. L’application tente de deviner quelle est la colonne contenant les valeurs de concentrations à l’aide des noms des colonnes. Cela peut toutefois nécessiter une correction. 
2. **Sélectionner (ou désélectionner) les distributions à ajuster aux données.** À noter que s’il y a un chevauchement dans l’ajustement de deux ou plusieurs fonctions de distribution, il y a aura alors une exagération de la forme de cet ajustement dans l’inférence multimodèle. Consultez cet article pour plus d’information.  La fonction peut prendre quelques secondes pour se mettre à jour. 
3. Mettre en forme le graphique à l’aide des entrées de la barre latérale et **télécharger le graphique et le tableau de l'évaluation de la qualité de l’ajustement des courbes de distribution** respectivement dans des fichiers .png et .csv. 

<center>
![Fit distributions in tab `2. Fit`](https://media.giphy.com/media/yIjxGcFKt0zK2vj38j/giphy.gif)
</center>

Information additionnelle sur **le tableau de l'évaluation de la qualité de l’ajustement des courbes de distribution**:
Les colonnes du tableau sont la distribution (dist), la statistique d’Anderson-Darling (ad), la statistique de Kolmogorov-Smirnov (ks), la statistique de Cramer-von-Mises (cmv), le critère d’information Akaike (aic), le critère d’information Akaike corrigé pour la taille de l’échantillon (aicc), le critère d’information Bayésien (bic), la différence AICc (delta) et le coefficient AICc basé sur le poids Akaike (weight). L’estimation de la fonction de distribution finale  est basée sur l’inférence multimodèle (à partir de l’AICc). La concentration présentant un risque est la concentration estimée d’une substance affectant le pourcentage sélectionné de l’ensemble des espèces.

### Étape 3: Estimation de la concentration présentant un risque ou du pourcentage d’espèces affectées
1. Sélectionner le **seuil (%) des espèces affectées** pour calculer **l’estimation de la concentration présentant un risque** OU  sélectionner **concentration** pour calculer le pourcentage d’espèces affecté par la concentration sélectionnée d’un substance. Cela affecte le graphique (ligne pointillée), le texte sous le graphique et le calcul des bornes de l’intervalle de confiance.
2. Sélectionner le nombre de **simulations bootstrap pour le calcul des bornes de l’intervalle de confiance.** Il est recommandé d’utiliser 10 000 simulations. Le calcul prendra environ 3 minutes.  Sélectionner un nombre moindre de simulations bootstrap pour réduire le temps de calcul. 

<center>

Nombre de simulation Bootstrap &nbsp;&nbsp; | Temps de calcul estimé
--- | ---
10 000 &nbsp; | 5 secondes
5 000 &nbsp;| 10 secondes
1 000 &nbsp;| 20 secondes
500 &nbsp;| 45 secondes

</center>

3. Les intervalles de confiance ne sont pas calculés automatiquement car cela prend un certain temps. Il faut cliquer sur le bouton `Obtenir bornes`.
4. **Le graphique peut être mis en forme** à l’aide des entrées disponibles dans la barre latérale et **le graphique et le tableau peuvent être téléchargés** respectivement dans des fichiers .png et .csv.
    - Si l’étiquette des espèces n’est pas entièrement visible sur le graphique, ajuster l’axe des X et la position des étiquettes en utilisant la fonction étiquette 'X-axis max' et 'Adjust label'.
    - L’éventail de couleurs provient de [ColorBrewer](http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3).
    - Si vous utilisez la même variable pour la couleur et la forme, fournir le même titre de légende permettra de combiner les deux en une seule légende. 
    
<center>
![Get hazard concentration estimates and confidence limits in tab `3. Predict`](https://media.giphy.com/media/xKb9nQsPFlqTCGzUgF/giphy.gif)
</center>

### Étape 4: obtenir le code R

Copier le code R pour reproduire la programmation. Le code est ajouté après chaque exécution dans l’application (par exemple : le code qui génère l’estimation des bornes de l’intervalle de confiance apparaitra après que `Obtenir bornes` aura été cliqué.  

<center>
![Get R code to reproduce results programmatically in tab `R Code`](https://media.giphy.com/media/XIgsL03rRnEfn8nNas/giphy.gif)
</center>

Pour générer un graphique avec l’intervalle de confiance, copier le code R et le coller dans R. Par la suite, régler ci = TRUE dans les fonctions d’estimation et de ssd_plot (predict and ssd_plot functions). 
 


---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# shinyssdtools

<!-- badges: start -->
[![Lifecycle: maturing](https://img.shields.io/badge/Lifecycle-Maturing-007EC6)](https://github.com/bcgov/repomountie/blob/master/doc/lifecycle-badges.md)
[![R-CMD-check](https://github.com/bcgov/shinyssdtools/workflows/R-CMD-check/badge.svg)](https://github.com/bcgov/shinyssdtools/actions)
[![Apache license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02848/status.svg)](https://doi.org/10.21105/joss.02848)
<!-- badges: end -->

`shinyssdtools` is a shiny application for fitting Species Sensitivity Distributions (SSDs) using [`ssdtools`](https://github.com/bcgov/ssdtools).

## Utilization

The most recent version is available online at <https://bcgov-env.shinyapps.io/ssdtools/>.

To install the development version from [GitHub](https://github.com/bcgov/shinyssdtools) and deploy locally use

``` r
# install.packages("devtools")
devtools::install_github("bcgov/shinyssdtools")
library(shinyssdtools)
shinyssdtools::run_app()
```

## Features

In addition to being a Graphical User Interface to the core functionality in the [`ssdtools`](https://github.com/bcgov/ssdtools) package, `shinyssdtools` also provides 

- a bilingual (English/French) interface;
- generation of R scripts for reproducibility;
- customization and downloads of plots and tables

## Information

For more information including how to cite `shinyssdtools` see [Dalgarno (2021)](https://doi.org/10.21105/joss.02848).

For a review of `ssdtools` and `shinysddtools` in the context of other SSD software packages see [Fox et al. (2021)](https://onlinelibrary.wiley.com/doi/10.1002/etc.4925).

## Assistance

To report bugs/issues/feature requests, please file an [issue](https://github.com/bcgov/shinyssdtools/issues/).

## Contribution

If you would like to contribute, please see our [CONTRIBUTING](CONTRIBUTING.md) guidelines.

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/bcgov/shinyssdtools/blob/master/.github/CONTRIBUTING.md). By participating in this project you agree to abide by its terms.

## Deploying to shinyapps.io

### Manually

Run the `deploy-app.R` script in the scripts directory (after setting the account argument to be your `shinyapps.io` account name).

### Automatically

If your `shinyapps.io` account name is the same as your GitHub account name simply make a commit in the master or dev branch and include `deploy app` in the message 
(after setting `SHINYAPPS_TOKEN` and `SHINYAPPS_SECRET` in your repository GitHub secrets). 
This triggers the `deploy-app.yml` GitHub action.

## License

The code is released under the Apache License 2.0

Copyright 2021 Province of British Columbia

Licensed under the Apache License, Version 2.0 (the "LICENSE");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at 

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

-----

<a rel="LICENSE" href="https://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence"
style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a><br /><span
xmlns:dct="http://purl.org/dc/terms/" property="dct:title">ssdtools</span> by <span
xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">the Province of British Columbia
</span> is licensed under a <a rel="LICENSE" href="https://creativecommons.org/licenses/by/4.0/">
Creative Commons Attribution 4.0 International License</a>.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run-app.R
\name{run_app}
\alias{run_app}
\alias{run_ssdtools_app}
\title{Run Shiny ssdtools Application}
\usage{
run_app()

run_ssdtools_app()
}
\description{
Run Shiny ssdtools Application
}
\section{Functions}{
\itemize{
\item \code{run_ssdtools_app}: 
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkgtemplate-package.R
\docType{package}
\name{shinyssdtools-package}
\alias{shinyssdtools}
\alias{shinyssdtools-package}
\title{shinyssdtools: ssdtools Shiny App}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

Shiny application for fitting species
    sensitivity distributions using ssdtools.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/bcgov/shinyssdtools}
  \item Report bugs at \url{https://github.com/bcgov/shinyssdtools/issues}
}

}
\author{
\strong{Maintainer}: Joe Thorley \email{joe@poissonconsulting.ca} (\href{https://orcid.org/0000-0002-7683-4592}{ORCID}) [contributor]

Authors:
\itemize{
  \item Seb Dalgarno \email{seb@northbeachconsulting.ca} (\href{https://orcid.org/0000-0002-3658-4517}{ORCID})
}

Other contributors:
\itemize{
  \item Ayla Pearson \email{ayla@poissonconsulting.ca} (\href{https://orcid.org/0000-0001-7388-1222}{ORCID}) [contributor]
}

}
\keyword{internal}
