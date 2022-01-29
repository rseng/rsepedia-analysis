# PyVOLCANS

> PyVOLCANS: A Python package to flexibly explore similarities and differences between volcanic systems

The main goal of PyVOLCANS is to help alleviate data-scarcity issues in volcanology, and contribute to developments in a range of topics, including (but not limited to): quantitative volcanic hazard assessment at local to global scales, investigation of magmatic and volcanic processes, and even teaching and scientific outreach. We hope that future users of PyVOLCANS will include any volcano scientist or enthusiast with an interest in exploring the similarities and differences between volcanic systems worldwide. Please visit our [wiki pages](https://github.com/BritishGeologicalSurvey/pyvolcans/wiki) for more information.

A new article on PyVOLCANS has been published in the Journal of Open Source Software. Please check its content here [![DOI](https://joss.theoj.org/papers/10.21105/joss.03649/status.svg)](https://doi.org/10.21105/joss.03649), and refer to it when using PyVOLCANS. Very many thanks!

## Installation instructions

PyVOLCANS can be installed from PyPI as follows:

```
pip install pyvolcans
```

This method adds `pyvolcans` to the virtual environment PATH so that it can be
used from any directory.


## API documentation and example usage

Users interact with PyVOLCANS via the command-line tool. The `--help` command describes the possible options.

```
$ pyvolcans --help

usage: pyvolcans [-h] [--apriori [APRIORI [APRIORI ...]]]
                 [-Ts TECTONIC_SETTING] [-G ROCK_GEOCHEMISTRY] [-M MORPHOLOGY]
                 [-Sz ERUPTION_SIZE] [-St ERUPTION_STYLE] [--count COUNT] [-w]
                 [-ovd] [-oad] [-W] [-v] [-pa] [-S] [-V]
                 volcano

positional arguments:
  volcano               Set target volcano name or Smithsonian ID (VNUM)

optional arguments:
  -h, --help            show this help message and exit
  --apriori [APRIORI [APRIORI ...]]
                        Provide one or more a priori analogue volcanoes
  -Ts TECTONIC_SETTING, --tectonic_setting TECTONIC_SETTING
                        Set tectonic setting weight (e.g. '0.2' or '1/5')
  -G ROCK_GEOCHEMISTRY, --rock_geochemistry ROCK_GEOCHEMISTRY
                        Set rock geochemistry weight (e.g. '0.2' or '1/5')
  -M MORPHOLOGY, --morphology MORPHOLOGY
                        Set volcano morphology weight (e.g. '0.2' or '1/5')
  -Sz ERUPTION_SIZE, --eruption_size ERUPTION_SIZE
                        Set eruption size weight (e.g. '0.2' or '1/5')
  -St ERUPTION_STYLE, --eruption_style ERUPTION_STYLE
                        Set eruption style weight (e.g. '0.2' or '1/5')
  --count COUNT         Set the number of top analogue volcanoes
  -w, --write_csv_file  Write list of top analogue volcanoes as .csv file
  -ovd, --output_volcano_data
                        Output volcano data (ID profile) for the selected
                        target volcano in a json-format file. NB. Verbose mode
                        needs to be activated to be able to use this feature.
  -oad, --output_analogues_data
                        Output volcano data (ID profile) for all the top
                        analogue volcanoes, for the selected target volcano
                        and weighting scheme, in a json-format file. NB.
                        Verbose mode needs to be activated to be able to use
                        this feature.
  -W, --website         Open GVP website for top analogue volcano
  -v, --verbose         Print debug-level logging output, ID profile for the
                        selected target volcano, and include single-criterion
                        analogy values, besides the total analogy values, in
                        the PyVOLCANS results.
  -pa, --plot_apriori   Generate bar plots displaying: (1) values of single-
                        criterion and total analogy between the target volcano
                        and any 'a priori' analogues chosen by the user; and
                        (2) percentages of 'better analogues' (for the target
                        volcano) than each of the 'a priori' analogues,
                        considering all volcanoes in the GVP database.
  -S, --save_figures    Save all generated figures
  -V, --version         Print PyVOLCANS package version and exit
```

Calling PyVOLCANS with a volcano name returns the top 10 analogue volcanoes:

```
$ pyvolcans Hekla

Top 10 analogue volcanoes for Hekla, Iceland (372070):
              name       country  smithsonian_id  total_analogy
       Torfajokull       Iceland          372050       0.941676
       Bardarbunga       Iceland          373030       0.921407
      Prestahnukur       Iceland          371070       0.919877
        Langjokull       Iceland          371080       0.915929
           Hengill       Iceland          371050       0.911855
 Brennisteinsfjoll       Iceland          371040       0.907751
        Kverkfjoll       Iceland          373050       0.906833
       Fremrinamar       Iceland          373070       0.905074
           Ecuador       Ecuador          353011       0.901611
     Marion Island  South Africa          234070       0.892960
```

In verbose mode (e.g. `$ pyvolcans Hekla --verbose`), PyVOLCANS provides further information regarding the weights selected for each volcanological criteria, as well as the values for each single-criterion analogy, which combined make up the `total_analogy` values (NB. The total analogy is a weighted average of the single-criterion analogies. Please see equation 1 in [Tierz et al., 2019](https://doi.org/10.1007/s00445-019-1336-3) for more details).

As of [PyVOLCANS v1.3.0](https://github.com/BritishGeologicalSurvey/pyvolcans/releases/tag/v1.3.0), the verbose mode also provides the _ID profile_ (i.e. summary of volcanological data available for VOLCANS calculations) for the specific volcano of interest. Please see [Tierz et al., 2019](https://doi.org/10.1007/s00445-019-1336-3) for further details on how these data are used to compute single-criterion and total-analogy values.

For example:

```
$ pyvolcans Hekla --verbose

PyVOLCANS: Supplied weights: {'tectonic_setting': 0.2, 'geochemistry': 0.2, 'morphology': 0.2, 'eruption_size': 0.2, 'eruption_style': 0.2}

Top 10 analogue volcanoes for Hekla, Iceland (372070):
              name       country  smithsonian_id  total_analogy  ATs        AG        AM       ASz       ASt
       Torfajokull       Iceland          372050       0.941676  0.2  0.188235  0.187584  0.180280  0.185577
       Bardarbunga       Iceland          373030       0.921407  0.2  0.188235  0.187584  0.172727  0.172861
      Prestahnukur       Iceland          371070       0.919877  0.2  0.188235  0.189474  0.169091  0.173077
        Langjokull       Iceland          371080       0.915929  0.2  0.188235  0.177193  0.169091  0.181410
           Hengill       Iceland          371050       0.911855  0.2  0.192157  0.173684  0.169091  0.176923
 Brennisteinsfjoll       Iceland          371040       0.907751  0.2  0.164706  0.184211  0.169091  0.189744
        Kverkfjoll       Iceland          373050       0.906833  0.2  0.188235  0.187584  0.169091  0.161923
       Fremrinamar       Iceland          373070       0.905074  0.2  0.188235  0.168421  0.169091  0.179327
           Ecuador       Ecuador          353011       0.901611  0.2  0.164706  0.194737  0.169091  0.173077
     Marion Island  South Africa          234070       0.892960  0.2  0.141176  0.200000  0.169091  0.182692

ID profile for Hekla, Iceland (372070):
{
  "name": "Hekla",
  "country": "Iceland",
  "smithsonian_id": 372070,
  "tectonic_setting": {
    "0.0": "Rift Oceanic Crust"
  },
  "geochemistry": {
    "Foidite": 0.0,
    "Phonolite": 0.0,
    "Trachyte": 0.0,
    "Trachyandesite/Basaltic trachyandesite": 0.0,
    "Phono-tephrite/Tephri-phonolite": 0.0,
    "Tephrite/Basanite/Trachybasalt": 0.0,
    "Basalt": 0.25,
    "Andesite": 0.25,
    "Dacite": 0.25,
    "Rhyolite": 0.25
  },
  "morphology": 0.39473684210526316,
  "eruption_size": {
    "VEI leq 2": 0.26666666666666666,
    "VEI 3": 0.35,
    "VEI 4": 0.2833333333333333,
    "VEI 5": 0.1,
    "VEI 6": 0.0,
    "VEI 7": 0.0,
    "VEI 8": 0.0
  },
  "eruption_style": {
    "Lava flow and/or fountaining": 0.8769230769230769,
    "Ballistics and tephra": 0.6923076923076923,
    "Phreatic and phreatomagmatic activity": 0.0,
    "Water-sediment flows": 0.15384615384615385,
    "Tsunamis": 0.015384615384615385,
    "Pyroclastic density currents": 0.09230769230769231,
    "Edifice collapse/destruction": 0.0,
    "Caldera formation": 0.0
  }
}
```

### Volcano naming conventions

Please note that some volcano names are composed of more than one word, such as Rincón de la Vieja (Costa Rica) or St. Helens (USA). In these cases, please wrap around the volcano name using quotation marks. For example:

```
$ pyvolcans "Rincon de la Vieja"

Top 10 analogue volcanoes for Rincon de la Vieja, Costa Rica (345020):
           name     country  smithsonian_id  total_analogy
    Sorikmarapi   Indonesia          261120       0.965415
         Zaozan       Japan          283190       0.963941
         Mahawu   Indonesia          266110       0.959122
 Akita-Yakeyama       Japan          283260       0.958577
         Lascar       Chile          355100       0.956498
     Miravalles  Costa Rica          345030       0.955624
     Hakkodasan       Japan          283280       0.955043
           Poas  Costa Rica          345040       0.954254
     Midagahara       Japan          283080       0.952637
       Maruyama       Japan          285061       0.951902
```

Please also note that some naming conventions used in the [Holocene Volcano List](https://volcano.si.edu/list_volcano_holocene.cfm) of the Global Volcanism Program (GVP) include the sorting of some words in the volcano name, separated by commas. For example: "Fournaise, Piton de la" (France), "Sawad, Harra Es-" (Yemen) or "Bravo, Cerro" (Colombia). One of the functionalities of PyVOLCANS is to provide a list of volcano-name suggestions, if the volcano name introduced by the user contains a typo and/or is arranged in a different word order. Please see the following command example:

```
$ pyvolcans "Nevados de Chillan"

PyVOLCANS: Nevados de Chillan not found! Did you mean:
                 name          country  smithsonian_id
  Chillan, Nevados de            Chile          357070
     Chachani, Nevado             Peru          354007
    Huila, Nevado del         Colombia          351050
      Casiri, Nevados             Peru          354060
    Cuernos de Negros      Philippines          272010
   Carrán-Los Venados            Chile          357140
      Chaine des Puys           France          210020
     Ruiz, Nevado del         Colombia          351020
             Red Hill    United States          327812
 Incahuasi, Nevado de  Chile-Argentina          355125
```

```
$ pyvolcans "Chillan, Nevados de"

Top 10 analogue volcanoes for Chillan, Nevados de, Chile (357070):
                name          country  smithsonian_id  total_analogy
          Guallatiri            Chile          355020       0.972271
 San Pedro-San Pablo            Chile          355070       0.970770
            San Jose  Chile-Argentina          357020       0.968865
             Galeras         Colombia          351080       0.966499
         Peuet Sague        Indonesia          261030       0.965005
         Zhupanovsky           Russia          300120       0.964308
             Bulusan      Philippines          273010       0.963809
          Chiginagak    United States          312110       0.963711
          Villarrica            Chile          357120       0.961412
            Cotopaxi          Ecuador          352050       0.960387
```

Finally, please be aware of possible synonyms and subfeatures of the volcanic systems listed in the [Holocene Volcano List](https://volcano.si.edu/list_volcano_holocene.cfm) of GVP. For example, "[Fagradalsfjall](https://volcano.si.edu/volcano.cfm?vn=371030&vtab=Subfeatures)" for "[Krýsuvík-Trölladyngja](https://volcano.si.edu/volcano.cfm?vn=371030&vtab=GeneralInfo)" (Iceland) or "[Sakurajima](https://volcano.si.edu/volcano.cfm?vn=282080&vtab=Subfeatures)" for "[Aira](https://volcano.si.edu/volcano.cfm?vn=282080&vtab=GeneralInfo)" (Japan). In these situations, we recommend PyVOLCANS users perform a simple Google Search using: "Fagradalsfjall GVP" or "Sakurajima GVP". In general, the first search result should point to the GVP website of the volcanic system of interest. Users can then use that volcano name to run their volcano analogues searches via `pyvolcans`.

For a comprehensive description of the purpose, input arguments and output variables for each of the functions and methods used by `pyvolcans`, please follow the link to the [source scripts](https://github.com/BritishGeologicalSurvey/pyvolcans/tree/main/pyvolcans).

Please also visit our [wiki pages](https://github.com/BritishGeologicalSurvey/pyvolcans/wiki) to find out more details on the usage of PyVOLCANS, as well as several example outputs for different commands.


## Community

### For users

The best way to get in touch to ask questions, submit bug reports or contribute code is by submitting an [issue](https://github.com/BritishGeologicalSurvey/pyvolcans/issues) in this repository.

### For developers

To modify the code, first clone the PyVOLCANS repository into your local machine:
```bash
git clone https://github.com/BritishGeologicalSurvey/pyvolcans
```

Then move to the root directory of the project (`pyvolcans`, which contains the `setup.py` file - 
please also check the `requirements.txt` file for a full list of dependencies in the code),
and run the following command to install PyVOLCANS in development mode, preferably within a
clean virtual environment:

```bash
python -m pip install -e .[dev]
```

The `-e` flag makes the files in the current working directory available
throughout the virtual environment and, therefore, changes are reflected straight away.
With this installation, it is no longer required to set the PYTHONPATH.

The `[dev]` part installs packages required for development e.g. `pytest`.

Run tests with:

```bash
pytest -v test
```

### Maintainers

`PyVOLCANS` was created by and is maintained by British Geological Survey
Volcanology and Digital Capabilities.

+ Pablo Tierz ([PTierz](https://github.com/PTierz))
+ Vyron Christodoulou ([mobiuscreek](https://github.com/mobiuscreek))
+ John A Stevenson ([volcan01010](https://github.com/volcan01010))

## Licence

`PyVOLCANS` and the associated `VOLCANS` Matlab scripts are distributed under the [LGPL v3.0 licence](LICENSE).
Copyright: © BGS / UKRI 2021
---
title: 'PyVOLCANS: A Python package to flexibly explore similarities and differences between volcanic systems'
tags:
  - Python
  - volcanology
  - volcanic hazard assessment
  - data scarcity
  - analogue volcanoes
authors:
  - name: Pablo Tierz^[Corresponding author]
    orcid: 0000-0001-8889-9900
    affiliation: 1
  - name: Vyron Christodoulou^[Now at The Data Lab, The Bayes Centre, Edinburgh, UK]
    orcid: 0000-0003-3835-3891
    affiliation: 1
  - name: John A. Stevenson
    orcid: 0000-0002-2245-1334
    affiliation: 1
  - name: Susan C. Loughlin
    affiliation: 1
affiliations:
 - name: British Geological Survey, The Lyell Centre, Edinburgh, UK.
   index: 1
date: 13 May 2021
bibliography: paper.bib
---

# Summary

There are over 1,400 volcanoes on Earth that have either erupted or shown signs of volcanic activity (e.g. fumaroles or hot springs) in, approximately, the last 12,000 years.
Of these, around 40-50 are erupting at any given time [@Siebert:2010; @GVP:2013].
Volcanoes provide a range of economic benefits, such as fertile soils, geothermal energy or valuable mineralisations, create a strong sense of belonging among local populations, and fascinate visitors.
However, volcanic systems can also generate hazardous phenomena, which may threaten local inhabitants, tourists and infrastructure at distances of up to tens or hundreds of kilometres.

In order to understand and quantify volcanic hazard, volcano scientists are faced with many questions.
How often do eruptions occur?
How big are they?
What style of eruption is possible (e.g. mainly explosive or effusive)?
From where on the volcano is eruptive activity sourced?
What areas around the volcanic system may be impacted?
Will there be any early warning signals?

Quantitative data to address these questions are scarce [@Loughlin:2015].
While a handful of volcanoes (e.g. Etna, Italy; Kīlauea, USA; Merapi, Indonesia) have been extensively studied, hundreds of volcanic systems around the world remain poorly-understood.
One possible mitigation to the issue of data scarcity in volcanology and volcanic hazard assessment is the use of _analogue volcanoes_ [@Newhall:2002; @Newhall:2017].
These are volcanoes with similar characteristics to a data-scarce volcano of interest.
Data and insights from the well-studied volcano(es) can be used to provide estimates for important variables, such as the number of eruptions during specific time windows or the size of those eruptions. Such methods have been used for many years, here we present the first tool enabling a structured and harmonised approach that can be applied worldwide.


# Statement of need

`PyVOLCANS` (Python VOLCano ANalogues Search) is an open-source tool that addresses the need for an objective, data-driven method for selection of analogue volcanoes.
It is based on the results of VOLCANS [@Tierz:2019], a first-of-its-kind method to quantify the analogy (or similarity) between volcanic systems, based on a structured combination of five volcanological criteria: tectonic setting, rock geochemistry, volcano morphology, eruption size, and eruption style.
`PyVOLCANS` provides a command-line interface to make the results from the VOLCANS study easily accessible to a wide audience.
`PyVOLCANS` is a versatile tool for volcano scientists, with potential applications ranging from investigating commonalities between volcanic systems [@Cashman:2014] to supporting probabilistic volcanic hazard assessment at local, regional and global scales. Exploring similarities and differences between volcanic systems using `PyVOLCANS` can also be useful for teaching and scientific outreach purposes.

Users can easily derive data-driven sets of _top_ analogue volcanoes (i.e. those with highest analogy) to any volcanic system listed in the reference database for recent global volcanism: the Volcanoes of the World Database, hosted by the Global Volcanism Program of the Smithsonian Institution [@GVP:2013].
Users can also choose the number of _top_ analogue volcanoes to investigate and can customise the importance (i.e. weight) that is given to each of the five aforementioned volcanological criteria.
Additionally, users can select a number of _a priori_ analogue volcanoes (i.e. volcanoes deemed as analogues by other means, such as expert knowledge) and assess their values of analogy with the target volcano to see how well they match on different criteria and if other volcanoes could be a better choice (\autoref{fig:figure1}).

![Values of single-criterion (colours) and total analogy (bar heights) between an example target volcano, Fuego (Guatemala)\*, and five _a priori_ analogues [please see @Tierz:2019, for more details].
ATs: Analogy in Tectonic setting; AG: Analogy in rock Geochemistry; AM: Analogy in volcano Morphology; ASz: Analogy in eruption Size; ASt: Analogy in eruption Style.
\*Number between brackets denotes the unique volcano identifier used by the GVP database.
\label{fig:figure1}](figure.png)

The results from the VOLCANS study have been used in research to explore the volcanological factors that influence the development of particular volcano morphologies (Philippa White, unpublished thesis); to constrain potential hazardous phenomena and hazard scenarios at a given target volcano, based on its analogue volcanoes [@Simmons:2021]; to quantify probability distributions of eruption sizes and probabilities of occurrence of diverse hazardous phenomena [@Tierz:2020]; or even to explore volcano analogies at regional scales, by generating sets of analogue volcanoes for tens of volcanic systems. The last two example applications have played a key role in developing quantitative hazard analyses for Ethiopian volcanoes, within the RiftVolc project (please see `PyVOLCANS` documentation for further details on these analyses).
Moreover, the future potential of VOLCANS/`PyVOLCANS`, particularly in the field of volcanic hazard assessment, has also been recognised in recent relevant publications in the area [@Papale:2021; @Marzocchi:2021]

We hope that the release of `PyVOLCANS` will encourage further studies based on data-driven selection of analogue volcanoes and that such analyses will continue to grow in number and diversity of their scientific purposes.


# Acknowledgements

The research leading to these results has been mainly supported by the UK National Capability Funding (Innovation Flexible Fund programme).
We most sincerely thank reviewers Jamie Farquharson and Meghan Jones, and Associate Editor Jay Hariharan, for the invaluable feedback and suggestions that they have provided on `PyVOLCANS`. This has certainly resulted in significant improvements to the code, which we hope future users will benefit from.
We would like to warmly thank Eliza Calder for all her work during the development of the VOLCANS method, and Sarah Ogburn for being one of the first people who _convinced_ us that we should develop an open-access application of VOLCANS, sooner rather than later.
Declan Valters is thanked for support with Python programming, and Fabio Dioguardi for his internal review.
Moreover, we would like to sincerely thank a number of colleagues with whom we shared very insightful conversations about analogue volcanoes and/or `PyVOLCANS`: Chris Newhall, Isla Simmons, Adriano Pimentel, Julia Crummy, Gezahegn Yirgu, Charlotte Vye-Brown, Lara Smale, Karen Fontijn, Ben Clarke, Susanna Jenkins, Elly Tennant, Pierre Barbillon, Elaine Spiller, Philippa White, Teresa Ubide, Sebastián García, Victoria Olivera, Jeremy Pesicek, Vanesa Burgos Delgado, Einat Lev, Jonty Rougier, Willy Aspinall, Paolo Papale, Monse Cascante and Thomas Giachetti.

# References
# VOLCano ANalogues Search (VOLCANS)

This is a README file that describes the structure and content of the Matlab scripts (Mathworks, 2012) that integrate
the VOLCANS method presented by [Tierz, Loughlin and Calder (2019)](https://doi.org/10.1007/s00445-019-1336-3).

This set of scripts was used to derive all the [analogy matrices](https://github.com/BritishGeologicalSurvey/pyvolcans/tree/main/pyvolcans/VOLCANS_mat_files/analogy_mats)
used by `PyVOLCANS`.

Please note that the VOLCANS Matlab code is provided _as-is_: in other words, without any cleaning for public use.
The eventual aim of us, creators of `PyVOLCANS` ([@PTierz](https://github.com/PTierz), [@mobiuscreek](https://github.com/mobiuscreek),
[@volcan01010](https://github.com/volcan01010)), is to port the entire Matlab code into Python, including it as part of future releases
of the `PyVOLCANS` code.

In the following, we provide: (1) a brief summary of the general process developed to generate the analogy matrices (for more details,
please see [Tierz et al., 2019](https://doi.org/10.1007/s00445-019-1336-3)); (2) overall descriptions of the main purpose of each of the
Matlab (`.m`) files available in the [current folder](https://github.com/BritishGeologicalSurvey/pyvolcans/tree/matlab-scripts/pyvolcans/VOLCANS_matlab_scripts); and (3) overall descriptions of the data contained in each of the csv files available
in folder [../VOLCANS_csv_files](https://github.com/BritishGeologicalSurvey/pyvolcans/tree/matlab-scripts/pyvolcans/VOLCANS_csv_files). Detailed specifications of the numbering adopted for categorical variables, as well as the units
in which numerical variables are displayed (unless this is described in the header row of the csv files) are provided in
[../VOLCANS_csv_files/VOLCANS_code_numbers](https://github.com/BritishGeologicalSurvey/pyvolcans/tree/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOLCANS_code_numbers).


## Volcano analogy

VOLCANS interprets volcano analogy as a quantitative measure of the similarity between any two volcanic systems listed in the Holocene
[Volcanoes Of The World database](https://volcano.si.edu/list_volcano_holocene.cfm) (v. 4.6.7), hosted by the Global Volcanism Program
(GVP) of the Smithsonian Institution. A volcano analogy of 0 denotes _no analogy_ between the two volcanoes, while a volcano analogy of
1 implies that the volcanoes are _perfect analogues_, considering the data available.

Five different volcanological criteria are considered in VOLCANS: tectonic setting, rock geochemistry, volcano morphology, eruption size
and eruption style. Volcano analogies based on only one criterion are termed: ''single-criterion analogies''. Volcano analogies based on
a structured combination (a weighted average) of two or more volcanological criteria are termed: ''multi-criteria or total analogies''.
Provided that each single-criterion analogy is a number between 0-1, and that the sum of weights defining the weighted average is equal
to 1, then each multi-criteria, or total analogy, is also a number between 0-1 (please see [Methods](https://link.springer.com/article/10.1007/s00445-019-1336-3#Sec6)
in Tierz et al., 2019, for more details). Please note that if there is no data available for a given volcano and a particular volcanological
criterion, then the corresponding single-criterion analogy values between that volcano and any other volcano in the GVP database are set to 0.

VOLCANS uses distance metrics to calculate each set of single-criterion analogies between any two volcanic systems in the GVP database.
These distance metrics are either based on: (a) linear distances between volcano-specific values, when the criterion is described by a
variable with single values for each volcano (tectonic setting, volcano morphology); (b) areal differences between volcano-specific Empirical
Cumulative Distribution Functions (ECDFs), when the criterion is described by a variable with a probability distribution (or histogram) for
each volcano (rock geochemistry, eruption size); or (c) normalised sum of differences between the frequency of occurrence of different groups
of hazardous events during eruptions at a given volcano (eruption style). Please see [Methods](https://link.springer.com/article/10.1007/s00445-019-1336-3#Sec6)
and [Table 2](https://link.springer.com/article/10.1007/s00445-019-1336-3/tables/2) in Tierz et al. (2019) for more details.

## Matlab scripts

The following [Matlab scripts](https://github.com/BritishGeologicalSurvey/pyvolcans/tree/matlab-scripts/pyvolcans/VOLCANS_matlab_scripts) have the purpose of calculating volcano analogy matrices that contain the single-criterion volcano analogy values
for each of the possible combinations of two volcanoes in the GVP database. All the analogy matrices are N x N, where N (= 1439) is the total number of
volcanic systems listed in the GVP v4.6.7 database used in [Tierz et al. (2019)](https://doi.org/10.1007/s00445-019-1336-3).
Please also note that any single-criterion analogy between volcanoes X and Y, e.g. rock geochemistry (AG<sub>XY</sub>), is symmetric.
That is: AG<sub>XY</sub> = AG<sub>YX</sub>.

### get_final_AT_allcross.m

It calculates the final single-criterion analogy matrix for tectonic setting. It requires to load the file [ATmatrices.mat](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_mat_files/data_mats/ATmatrices.mat), which is derived
from the script [votw_analysis.m](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_matlab_scripts/votw_analysis.m).

### get_final_AG_allcross.m

It calculates the final single-criterion analogy matrix for rock geochemistry. It requires to load the file [AGmatrices_ALU_QUET.mat](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_mat_files/data_mats/AGmatrices_ALU_QUET.mat), which is
derived from the script [votw_analysis.m](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_matlab_scripts/votw_analysis.m).

### get_final_AM_allcross.m

It calculates the final single-criterion analogy matrix for volcano morphology. It requires to load the file [AMmatrices_QUET.mat](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_mat_files/data_mats/AMmatrices_QUET.mat), which is
derived from the script [morphology_processing.m](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_matlab_scripts/morphology_processing.m).

### get_final_ASz_allcross.m

It calculates the final single-criterion analogy matrix for eruption size. It requires to load the file [ASzmatrices_SINA.mat](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_mat_files/data_mats/ASzmatrices_SINA.mat), which is
derived from the script [eruption_size_processing.m](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_matlab_scripts/eruption_size_processing.m).

### get_final_ASt_allcross.m

It calculates the final single-criterion analogy matrix for eruption style. It requires to load the file [AStmatrices_SINA.mat](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_mat_files/data_mats/AStmatrices_SINA.mat), which is
derived from the script [eruption_style_processing.m](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_matlab_scripts/eruption_style_processing.m).

### votw_analysis.m

It performs a preliminary analysis of the volcano data in the GVP database (v4.6.7), and derives the data required to calculate the single-criterion analogy matrices for tectonic setting and rock geochemistry. It requires to import the data in the file: [VOTW467_8May18_volcano_data.csv](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOTW467_8May18_volcano_data.csv) to carry out the aforementioned tasks. Please see below for a brief description of the csv files.

### morphology_processing.m

It performs a preliminary analysis of the morphology databases used by VOLCANS:
[Pike and Clow (1981)](https://www.researchgate.net/publication/259487495_Revised_classification_of_terrestrial_volcanoes_and_catalog_of_topographic_dimensions_with_new_results_on_edifice_volume)
and [Grosse et al. (2014)](https://doi.org/10.1007/s00445-013-0784-4), and derives the data required to calculate the single-criterion analogy matrix for volcano morphology.
It requires to import the data in the files: [PC81_GR2014_data.csv](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/PC81_GR2014_data.csv), and [VOTW467_8May18_volcano_data.csv](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOTW467_8May18_volcano_data.csv) to carry out the aforementioned tasks.
Please see below for a brief description of the csv files.

### eruption_size_processing.m

It derives the data required to calculate the single-criterion analogy matrix for eruption size. It requires to import the data in the files:
[VOTW467_8May18_eruption_data.csv](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOTW467_8May18_eruption_data.csv), [VOTW467_8May18_volcano_data.csv](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOTW467_8May18_volcano_data.csv), and [MeadMagill2014_June2018.csv](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/MeadMagill2014_June2018.csv) to carry out the aforementioned task.
Please see below for a brief description of the csv files, and [Methods](https://link.springer.com/article/10.1007/s00445-019-1336-3#Sec6) in Tierz et al. (2019)
for more details on the procedure applied in `VOLCANS` to account for under-recording of eruptions.

### eruption_style_processing.m

It derives the data required to calculate the single-criterion analogy matrix for eruption style. It requires to import the data in the files:
[VOTW467_8May18_event_data.csv](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOTW467_8May18_event_data.csv), [VOTW467_8May18_eruption_data.csv](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOTW467_8May18_eruption_data.csv), and [VOTW467_8May18_volcano_data.csv](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOTW467_8May18_volcano_data.csv) to carry out the aforementioned task.
Please see below for a brief description of the csv files.

## Data (csv) files

### VOTW467_8May18_volcano_data.csv

It contains the list of Holocene volcanoes in the GVP database (v4.6.7), together with data on their tectonic setting, required to calculate the analogy
with the same name, as well as data on their major and minor rock types: the latter required to calculate the analogy in rock geochemistry.
NB. It corresponds with [ESM2](https://static-content.springer.com/esm/art%3A10.1007%2Fs00445-019-1336-3/MediaObjects/445_2019_1336_MOESM2_ESM.csv) in Tierz et al. (2019).

### VOTW467_8May18_eruption_data.csv

It contains the list of Holocene eruptions from volcanoes in the GVP database (v4.6.7), together with other data required to calculate the analogy in
eruption size (e.g. Volcanic Explosivity Index or eruption date).
NB. It corresponds with [ESM3](https://static-content.springer.com/esm/art%3A10.1007%2Fs00445-019-1336-3/MediaObjects/445_2019_1336_MOESM3_ESM.csv) in Tierz et al. (2019).

### VOTW467_8May18_event_data.csv

It contains the list of ''eruptive events'' during Holocene eruptions at volcanoes in the GVP database (v4.6.7), together with their grouping according to different hazardous phenomena: the latter required to calculate the analogy in eruption style (please see [Methods](https://link.springer.com/article/10.1007/s00445-019-1336-3#Sec6) in Tierz et al., 2019,
for more details).
NB. It corresponds with [ESM4](https://static-content.springer.com/esm/art%3A10.1007%2Fs00445-019-1336-3/MediaObjects/445_2019_1336_MOESM4_ESM.csv) in Tierz et al. (2019). 

### MeadMagill2014_June2018.csv

It contains the median values for the date of completeness (i.e. date after which the eruptive record for a given volcano is considered complete) extracted
from [Table 1](https://link.springer.com/article/10.1007/s00445-014-0874-y/tables/1) in Mead and Magill (2014).

### PC81_GR2014_data.csv

It contains the list of Holocene volcanoes in the GVP database (v4.6.7) for which there is morphological data available in the
[Pike and Clow (1981)](https://www.researchgate.net/publication/259487495_Revised_classification_of_terrestrial_volcanoes_and_catalog_of_topographic_dimensions_with_new_results_on_edifice_volume)
and [Grosse et al. (2014)](https://doi.org/10.1007/s00445-013-0784-4) volcano morphology databases. These data are used to calculate the analogy in volcano morphology.
NB. It corresponds with [ESM5](https://static-content.springer.com/esm/art%3A10.1007%2Fs00445-019-1336-3/MediaObjects/445_2019_1336_MOESM5_ESM.csv) in Tierz et al. (2019).

## Authorship and Licence

`VOLCANS` was devised by [Pablo Tierz](https://www.bgs.ac.uk/people/tierz-lopez-pablo/) ([@PTierz](https://github.com/PTierz)),
[Susan C. Loughlin](https://www.bgs.ac.uk/people/loughlin-susan/) (British Geological Survey),
and [Eliza S. Calder](https://www.research.ed.ac.uk/en/persons/eliza-calder) (University of Edinburgh).
All Matlab scripts were written by [Pablo Tierz](https://www.bgs.ac.uk/people/tierz-lopez-pablo/) ([@PTierz](https://github.com/PTierz)).

`PyVOLCANS` and the associated VOLCANS Matlab scripts are distributed under the [LGPL v3.0 licence](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/main/LICENSE).
Copyright: © BGS / UKRI 2021.

# VOLCANS: code numbering

This README file details the code numbering adopted in the different csv files inside the folder:
[../../VOLCANS_csv_files](https://github.com/BritishGeologicalSurvey/pyvolcans/tree/matlab-scripts/pyvolcans/VOLCANS_csv_files).
This numbering was adopted to ease the calculations performed by the Matlab scripts that compose the VOLCANS method.

Please note that, in some cases, the code numbering can be directly checked on `.xls` files located inside
the [current folder](https://github.com/BritishGeologicalSurvey/pyvolcans/tree/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOLCANS_code_numbers) (please see descriptions below). When this is not available, a legend is provided in this file, showing the correspondence between categorical variables and code numbering (please see next sub-sections).

Please also note that the `NO DATA` value used across all data files is `-9999`.

## VOTW467_8May18_volcano_data.csv

Code numbering for all the categorical variables can be found in the file [VOTW467_8May18_Holocene_list_textdata_portable.xls](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOLCANS_code_numbers/VOTW467_8May18_Holocene_list_textdata_portable.xls)
(the variable displaying the code number in the spreadsheet is given between brackets):

Country (country#), Region (reg#), Subregion (subreg#), Primary Volcano Type (volctype#), Activity Evidence (procdate1),
Dominant Rock Type (domrock#), Tectonic Setting (tectset#).

Please also note the following four aspects:

1. The variable ''Volcano Number'' (also expressed as ''VNUM'' in other data files) denotes a unique volcano identifier that
the GVP database assigns to each volcanic system listed in the database.

2. The variables ''Major Rock 1-5'' and ''Minor Rock 1-5'' use the same code numbering as the variable ''Dominant Rock Type''.

3. The variable ''Last Known Eruption'' is given in years from current era, with negative values expressing dates before
current era (BCE). The same convention is applicable to variables ''Start Year'' in file [VOTW467_8May18_eruption_data.csv](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOTW467_8May18_eruption_data.csv),
''Eruption Start Year'' in file [VOTW467_8May18_event_data.csv](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOTW467_8May18_event_data.csv), and ''Date of completeness'' in [MeadMagill2014_June2018.csv](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/MeadMagill2014_June2018.csv).

4. The variables ''Latitude'' and ''Longitude'' (also included in the file [VOTW467_8May18_Eruption_list_textdata_portable.xls](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOLCANS_code_numbers/VOTW467_8May18_Eruption_list_textdata_portable.xls))
correspond to the approximate spatial location of the volcanic system, and are expressed in decimal degrees.

## VOTW467_8May18_eruption_data.csv

Code numbering for the categorical variable ''VEI Modifier'' can be found in the file [VOTW467_8May18_Eruption_list_textdata_portable.xls](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOLCANS_code_numbers/VOTW467_8May18_Eruption_list_textdata_portable.xls)
(please see column ''VEImod#'').

For the variable ''Eruption Category'', the following code numbering is applied:

1 = Confirmed Eruption

3 = Uncertain Eruption 

## VOTW467_8May18_event_data.csv

Code numbering for the categorical variable ''Event Type'' can be found in the file [VOTW467_8May18_Event_list_textdata_portable.xls](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOLCANS_code_numbers/VOTW467_8May18_Event_list_textdata_portable.xls)
(please see column ''event#'').

For the variable ''group#'', the following code numbering is applied (please also see [Tierz et al., 2019](https://doi.org/10.1007/s00445-019-1336-3)):

1 = Lava flow and/or fountaining

2 = Ballistics and Tephra

3 = Phreatic and phreatomagmatic activity 

4 = Water-sediment flows

5 = Tsunamis

6 = Pyroclastic density currents

7 = Edifice collapse/destruction

8 = Caldera formation

## MeadMagill2014_June2018.csv

Code numbering for the categorical variable ''Region/Country code'' corresponds with the values in column ''Name''.
These values are the same as those found in file [VOTW467_8May18_Holocene_list_textdata_portable.xls](https://github.com/BritishGeologicalSurvey/pyvolcans/blob/matlab-scripts/pyvolcans/VOLCANS_csv_files/VOLCANS_code_numbers/VOTW467_8May18_Holocene_list_textdata_portable.xls), for variables
''Region (reg#)'' and ''Country (country#)''. The variable ''Region/Country'' denotes whether a given entry corresponds
with a region (0) or a country (1).

## PC81_GR2014_data.csv

The following legend and code numbering is applicable to the merged database of volcano morphology:

- `VNUM`: Volcano Number (i.e. GVP unique volcano identifier)

- `Sub-feature`: The entry corresponds with a sub-feature (1) of a unique VNUM entry in GVP v4.6.7, or does not (0). In the latter case, the entry represents the _main_ volcanic feature of the corresponding VNUM.

- `small crater?`: The entry is (1) or is not (0) associated with a ''small crater'' flag according to the [Grosse et al. (2014)](https://doi.org/10.1007/s00445-013-0784-4) database.

- `NB`: It indicates whether there are any relevant remarks about the calculation of the morphology features displayed in the entry (these remarks coming from [Pike and Clow (1981)](https://www.researchgate.net/publication/259487495_Revised_classification_of_terrestrial_volcanoes_and_catalog_of_topographic_dimensions_with_new_results_on_edifice_volume), hereinafter PC81, and [Grosse et al. (2014)](https://doi.org/10.1007/s00445-013-0784-4), hereinafter GR2014):
    
    * 1	= One or more dimensions uncertain and subject to revision.
    * 2 = Height and width calculated from volumetric information only.
    * 3 = Estimate of lake depth included in total crater depth.
    * 4 = Island volcano: height and width down to sea level only.
    * 5 = Height and width down to sea level only, plus one or more dimensions uncertain and subject to revision.

- `W*_T from`: It indicates from which database the value of `W*` (volcano's edifice half-width, please see below) was taken from:

    * `W*_T from` = 0 &rarr; [PC81](https://www.researchgate.net/publication/259487495_Revised_classification_of_terrestrial_volcanoes_and_catalog_of_topographic_dimensions_with_new_results_on_edifice_volume) database.
    * `W*_T from` = 1 &rarr; [GR2014](https://doi.org/10.1007/s00445-013-0784-4) database.

- `d`: Average diameter of crater rim-crest ([PC81](https://www.researchgate.net/publication/259487495_Revised_classification_of_terrestrial_volcanoes_and_catalog_of_topographic_dimensions_with_new_results_on_edifice_volume)) or ''Crater Width'' variable in [GR2014](https://doi.org/10.1007/s00445-013-0784-4).

- `h`: Average depth of crater floor below rim crest ([PC81](https://www.researchgate.net/publication/259487495_Revised_classification_of_terrestrial_volcanoes_and_catalog_of_topographic_dimensions_with_new_results_on_edifice_volume)) or ''Crater Depth'' variable in [GR2014](https://doi.org/10.1007/s00445-013-0784-4).

- `H`: Average height of rim crest above pre-volcano topographic datum ([PC81](https://www.researchgate.net/publication/259487495_Revised_classification_of_terrestrial_volcanoes_and_catalog_of_topographic_dimensions_with_new_results_on_edifice_volume)) or ''Height Max'' variable in [GR2014](https://doi.org/10.1007/s00445-013-0784-4).

- `W*`:	W\* = W + (1/2 · d) in [PC81](https://www.researchgate.net/publication/259487495_Revised_classification_of_terrestrial_volcanoes_and_catalog_of_topographic_dimensions_with_new_results_on_edifice_volume), where W denotes the average half-width of the volcano, i.e. flank rim crest to edge of edifice; or W\* = 1/2 · W<sub>basal</sub> in [GR2014](https://doi.org/10.1007/s00445-013-0784-4), where W<sub>basal</sub> denotes the ''Basal Width'' variable.

- `C`: Circularity of rim crest, i.e. area inscribed circle/area circumscribed circle (only available in [PC81](https://www.researchgate.net/publication/259487495_Revised_classification_of_terrestrial_volcanoes_and_catalog_of_topographic_dimensions_with_new_results_on_edifice_volume)).

- `ave. ei`: Ellipticity index [average contours] variable in [GR2014](https://doi.org/10.1007/s00445-013-0784-4).

- `T`: T = d/(2W + d) if taken from [PC81](https://www.researchgate.net/publication/259487495_Revised_classification_of_terrestrial_volcanoes_and_catalog_of_topographic_dimensions_with_new_results_on_edifice_volume); or T = W<sub>summit</sub>/W<sub>basal</sub> (''Summit Width/Basal Width'' variable) if taken from [GR2014](https://doi.org/10.1007/s00445-013-0784-4).

- `sv`: Secondary vents, i.e. ''Sec. Peaks [total]'' variable in [GR2014](https://doi.org/10.1007/s00445-013-0784-4).
