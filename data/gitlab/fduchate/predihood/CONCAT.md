---
title: "Predihood: an open-source tool for predicting neighbourhoods' information"
tags:
  - Python
  - MongoDB
  - data management
  - neighbourhood
  - prediction
  - machine learning
authors:
  - name: Nelly Barret
    orcid: 0000-0002-3469-4149
    affiliation: 1
  - name: Fabien Duchateau
    orcid: 0000-0001-6803-917X
    affiliation: 1
  - name: Franck Favetta
    orcid: 0000-0003-2039-3481
    affiliation: 1
affiliations:
 - name: LIRIS UMR5205, Université Claude Bernard Lyon 1, Lyon, France
   index: 1
date: 16 September 2020
bibliography: paper.bib
---
# Introduction

Neighbourhoods are a widespread concept in studies from diverse domains such as health, social sciences, or biology. For instance, Japanese researchers investigated the relationships between social factors and health by taking into account behavioral risks and housing and neighbourhood environments [@takada2014japanese]. In a British study, authors describe how living areas impact physical activities, from which they determine a walkability index at the neighbourhood level for improving future urban planning [@frank2010development]. Another survey describes the luxury effect, i.e., the impact of wealthy neighbourhoods on the surrounding biodiversity [@leong2018biodiversity]. Several works focus on qualifying neighbourhoods using social networks. For instance, the Livehoods project defines and computes neighbourhoods dynamics [@cranshaw2012livehoods], while the Hoodsquare project detects similar areas based on Foursquare check-ins [@zhang2013hoodsquare]. Crowd-based systems are interesting but may sometimes be biased, and they require technical skills to extract relevant data. DataFrance is an interface that integrates data from several sources, such as indicators provided by the National Institute of Statistics (INSEE), geographical information from the National Geographic Institute (IGN), and surveys from newspapers for prices (L'Express). DataFrance enables the visualization of hundreds of indicators but makes it difficult to judge the main characteristics of a neighbourhood and is limited to France. Despite all of these existing applications, there is no simple tool to visualize and predict insights about neighbourhoods.

The Predihood tool fills this gap by defining neighbourhoods, their characteristics, and variables to be predicted. It also includes a cartographic interface for searching and displaying information about neighbourhoods. Domain experts can provide a few examples of annotated neighbourhoods, and Predihood provides a configuration interface for using popular machine-learning algorithms in order to predict variables for the remaining neighbourhoods. One of the most recent applications of Predihood was measuring the impact and influence of a neighbourhood's environment on the decision-making process when people move to another city [@data2020]. 

Predihood mainly targets non-programmers users (e.g., researchers in social sciences or history) due to its simplicity for running and configuring predictive algorithms. It can be extended to other application domains: measuring the pollution degree in neighbourhoods, determining whether a certain neighbourhood is suitable as a stopover for migratory birds, predicting neighbourhood evolution based on historical data, etc. This paper describes the main features of Predihood and how to extend them.

# Methodology

Predihood provides the following functionalities:

- adding new neighbourhoods and indicators to describe them;
- predicting variables of a neighbourhood by configuring and using predefined algorithms;
- adding new predictive algorithms.

To facilitate understanding, we describe and illustrate these functionalities based on a simple example that aims to evaluate which neighbourhood is preferable for migratory birds to make a temporary stop. We only include three indicators per neighbourhood: the percent of greens, the percent of buildings and the degree of human pressure. A single variable `migration zone` qualifies a neighbourhood from _favorable_ to _unfavorable_.

## Adding new datasets

As Predihood is a general-purpose application, it enables contributors to add their own datasets. The key concept of a dataset is the neighbourhood, which is represented as a [GeoJSON object](https://geojson.org/) including:

- a geometry (multi-polygons), which describes the shape of the neighbourhood. This crucial data is not only useful for cartographic visualization but also enable automatic calculations such as area;
- properties, which are divided into two categories. Descriptive information (e.g., name, city postcode) mainly aims to improve display while a set of quantitative indicators is used to predict the values of the variable.

Besides, some neighbourhoods have to be manually annotated, a task typically performed by domain experts. To add a new dataset, it is necessary to store them as GeoJSON and make them accessible by Predihood, for instance, in a document-oriented database. 

To build a dataset for the _bird migration_ example, it is necessary to collect and integrate data sources about neighbourhoods located in the studied geographic area. Values for the three indicators should be provided, and a value to the `migration zone` variable should be assigned to a few neighbourhoods.

## Predicting

Machine learning algorithms require data preparation by grouping relevant properties and variables. We illustrate this step on the _bird migration_ dataset, as shown in Figure \ref{fig:dataset}. Predihood produces a table composed of the identifier and the name of the neighbourhood (grey columns), its indicators (yellow columns) that could be normalized by factors such as density of population (green columns), and optionally the assessment of researchers for the `migration zone` variable (blue column). The objective of Predihood is to automatically fill question marks for neighbourhoods that are not yet assessed.

![A subset of the _bird-migration_ dataset.\label{fig:dataset}](doc/predihood-indicators-bird.png)

To perform prediction, a selection process first selects subsets of relevant indicators, since too many indicators may degrade performance. The tool automatically reduces the number of indicators, e.g., by removing indicators with a unique value or those highly correlated using Spearman coefficient. Then, Predihood generates 7 lists of indicators, containing from 10 to 100 indicators. The current version of Predihood includes 8 predictive algorithms from [scikit-learn](https://scikit-learn.org/) (e.g., Random Forest, KNeighbours) [@scikit-learn], depending on the algorithm, small or large lists of indicators may be more effective. 

Predihood provides a cartographic web interface based on [Leaflet](https://leafletjs.com/) and [Open Street Map](https://www.openstreetmap.org/), as shown in Figure \ref{fig:cartographic-interface}. For this example we use the search query "_lyon_" (left panel) and all neighbourhoods containing this query in their name or their city names are shown in blue on the map. We select the neighbourhood [_Le parc_](https://en.wikipedia.org/wiki/Parc_de_la_Tête_d%27or) and run the Random Forest classifier: migration is considered _very favorable_ in this area (for the seven lists of indicators). Indeed, this park seems relevant for bird migration as it has nice green areas for birds and it is a healthy environment for them.

![Screenshot of the cartographic interface of Predihood.\label{fig:cartographic-interface}](doc/predihood-predictions-bird.png)

## Adding new algorithms

Because prediction is a complex task, testing specific algorithms tuned with different parameters and comparing their results may help increase the overall quality. In order to facilitate this task, Predihood proposes a generic and easy-to-use programming structure for machine learning algorithms, based on Scikit-learn algorithms. Thus, experts can implement hand-made algorithms and run experiments in Predihood. Adding a new algorithm only requires 4 steps:

1. Create a new class that represents your algorithm, e.g. `MyOwnClassifier`, which inherits from `Classifier`;
2. Implement the core of your algorithm by coding the `fit()` and `predict()` functions. The `fit` function aims at fitting your classifier on assessed neighbourhoods while the `predict` function aims at predicting variable(s) for a given neighbourhood;
3. Add the `get_params()` function, which returns a dictionary of parameters along with their value, in order to be compatible with the Scikit-learn framework;
4. Comment your classifier with the [Numpy style](https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard) [@harris2020array] so that Predihood automatically extracts its parameters and enables their tuning in the dedicated interface.

Below is a simple example to illustrate the aforementioned steps. Note that new algorithms are directly loaded into Predihood.

```python
# file ./algorithms/MyOwnClassifier.py
from predihood.classes.Classifier import Classifier


class MyOwnClassifier(Classifier):
  """Description of the classifier.
  Parameters
  ------------
  a : float, default=0.01
    Description of a.
  b : int, default=10
    Description of b.
  """

  def __init__(self, a=0.01, b=10):
    self.a = a
    self.b = b

  def fit(self, X, y):
    # do stuff here

  def predict(self, df):
    # do stuff here

  def get_params(self, deep=True):
    # suppose this estimator has parameters "a" and "b"
    return {"a": self.a, "b": self.b}
```

To facilitate experiments, Predihood provides an interface for easily tuning and testing algorithms on a dataset, as shown in Figure \ref{fig:tuning-interface}. The left panel stands for the selection of an algorithm and the tuning of its parameters and hyper parameters, such as training and test sizes. Note that two options enable to remove outliers and under-represented neighbourhoods (for a given variable) without directly modifying the dataset. On the right, the table illustrates the accuracies obtained for each list of indicators (generated during the selection process) and each variable. Results can be exported in XLS with the blue download icon. Here, we notice that the new algorithm _MyOwnClassifier_ has been chosen, and its parameters (_a_ and _b_) can be configured. We have performed 2 runs, the former with the Random Forest classifier and the latter with MyOwnClassifier. Best predictive results are achieved with all indicators (green cells).

![Screenshot of algorithmic interface of Predihood.\label{fig:tuning-interface}](doc/predihood-accuracies-bird.png)

# Current applications of Predihood

Our Predihood tool was presented during the DATA conference [@data2020] posing the evaluation of whether people choose a similar neighbourhood environment when moving elsewhere as the main research challenge. The tool is bundled with data from France using the [mongiris](https://gitlab.liris.cnrs.fr/fduchate/mongiris) project (in which unit divisions named [IRIS](https://www.insee.fr/en/metadonnees/definition/c1523/) stand for neighbourhoods), this dataset contains about 50,000 neighbourhoods with 640 indicators (about population, shops, buildings, etc.). Six environment variables have been defined (_building type_, _building usage_, _landscape_, _social class_, _morphological position_ and _geographical position_), and 270 neighbourhoods were annotated by social science researchers (one to two hours per neighbourhood to investigate building and streets pictures, parked cars, facilities and green areas from services such as Google Street View). Prediction results achieved by Predihood using 6 algorithms from Scikit-learn range from 30% to 65% accuracy depending on the environment variable, and designing new algorithms could help improving these results.

The open-source project is available here: [https://gitlab.com/fduchate/predihood](https://gitlab.com/fduchate/predihood).

# Acknowledgements

This work has been partially funded by LABEX IMU (ANR-10-LABX-0088) from Université de Lyon, in the context of the program "Investissements d'Avenir" (ANR-11-IDEX-0007) from the French Research Agency (ANR). In addition to Scikit-learn and Numpy, Predihood relies on other dependencies, namely Pandas [@jeff_reback_2020_4161697], seaborn [@Waskom2021] and matplotlib [@Hunter:2007].

# References
# Predihood

Predihood is an application for predicting information about neighbourhoods (e.g., environment characteristics, bird migration possibilities, health issues). It makes it very easy, even for non programmers, to use predictive algorithms. Besides, the tool is extensible: new datasets and predictive algorithms can be added into the system. 
A cartographic interface enables the visualization of neighbourhoods along with their indicators (which describe them, such as the number of bakeries, the average income, or the number of houses over 250m^2) and the prediction results for selected neighbourhoods. A tuning interface enables to configure and test different machine learning algorithms on a given dataset.

Predihood includes a dataset of 50,000 French neighbourhoods (`hil`) used in a [research project](https://imu.universite-lyon.fr/appels-en-cours-et-bilans/2017-en-cours/hil-artificial-intelligence-to-facilitate-property-searches-system-of-recommendations-with-spatial-and-non-spatial-visualisation-for-property-search-2017/) to predict the environment of a neighbourhood (e.g., social class, type of landscape) based on hundreds of indicators (about population, shops, buildings, etc.).
The tool also includes a small test dataset (`bird-migration`) to demonstrate how to add new datasets. 

Predihood is provided under a [GNU General Public License v3.0](https://gitlab.com/fduchate/predihood/-/blob/master/LICENSE).

Contributions are welcome, following the [community guidelines](https://gitlab.com/fduchate/predihood/-/blob/master/CONTRIBUTING.md).

Predihood has been published in Journal of Open Source Software (JOSS). Please cite our work if you use it in a scientific publication.

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02805/status.svg)](https://doi.org/10.21105/joss.02805)

## Installation

Predihood includes two components: the Python application ([predihood repository](https://gitlab.com/fduchate/predihood)) and the data management library ([mongiris repository](https://gitlab.liris.cnrs.fr/fduchate/mongiris)).

Note that the two repositories takes about 2.8GB (including datasets) on the disk.

### Installation using Docker (recommended)

This method requires the use of [Docker](https://www.docker.com/). This method builds two Docker images (for a total size of 5.7GB on the disk, including all libraries and datasets loads).

First, clone the two repositories with the following commands. Note that both cloned repositories should be placed into a new (empty) directory.

```
git clone https://gitlab.liris.cnrs.fr/fduchate/mongiris
git clone https://gitlab.com/fduchate/predihood.git
```

Go into the downloaded `predihood/` directory and run in a terminal:

```
docker-compose up
```

This command deploys two containers, one for the application (`predihood`) and the other for the database (`db-predihood`). On the first run, the database container imports two datasets (`hil` and `bird-migration`), which may take a few minutes. 
Note that the application container may generate a database connection error (if the database container is not ready before timeout, which could occur at the first run on low machines), but it automatically restarts.

After some logging information, go to [http://127.0.0.1:8081/](http://127.0.0.1:8081/) in your browser (preferably Firefox or Chrome) to use the application.


### Manual installation

Although the Docker method is highly recommended, it is also possible to manually install Predihood and its dependencies. Note that some issues may occur due to version or package conflicts.

Requirements:

- Python, version >=3.8
- [MongoDB](https://www.mongodb.com/), version >=4 for importing the database about neighbourhoods.


First, clone the two repositories with the following commands:

```
git clone https://gitlab.liris.cnrs.fr/fduchate/mongiris
git clone https://gitlab.com/fduchate/predihood.git
```

Next, go in the `mongiris` directory and install the mongiris application:

```
python3 -m pip install -e .
```

Note that the download time may be quite long, as the mongiris API includes two datasets (760 MB).

Then import datasets into the MongoDB database: run the MongoDB server (`mongod`) and execute the following commands (from the MongoDB's executable directory if needed):

```
# import dataset 'hil' as a MongoDB dump
./mongorestore --archive=/path/to/dump-dbinsee.bin
# import dataset 'bird-migration' as a collection of JSON documents
./mongoimport --db=dbmigration -c=collmigration --file=/path/to/dump-bird-neighbourhoods.json	
./mongoimport --db=dbmigration -c=collindic --file=/path/to/dump-bird-indicators.json
```

where `/path/to/` is the path to the dataset files (provided with the package mongiris in `mongiris/data/dumps/`). A tip is to move the dataset files into the MongoDB binary (`PATH/TO/MONGODB/bin`). You may have to create these folders for Mongodb: `data/db` under `PATH/TO/MONGODB/bin` and run `./mongod --dbpath=./data/db`. 

Finally, go in the `predihood` directory and install the predihood application:

```
python3 -m pip install -e . -r requirements.txt
```

For running Predihood, go in the `predihood/predihood/` directory (which contains `main.py`) and run in a terminal:

```
python3 main.py [path/to/config.json]
```

The application may have an argument which is the path to the configuration file of the dataset to be loaded. By default, the dataset _hil_ is loaded (see the _Datasets_ section for more information).

After some logging information, go to [http://localhost:8081/](http://localhost:8081/) in your browser (preferably Firefox or Chrome) to use the application.


## Example usage

For the cartographic interface, an example would be:

1. Type a query in the panel on the left, e.g. "Lyon". This will display all neighbourhoods that contain "Lyon" in their name or their township.
2. Click on a neighbourhood (which are the small areas in blue). A tooltip will appear with some information about the neighbourhood. There are more information (list of all indicators) when clicking on the "More details" link.
3. In order to predict variables of the neighbourhood, you have to choose the classifier. The "Random Forest" classifier is recommended by default. After some seconds, predictions will appear in the tooltip. Prediction results can be exported as tablesheets (XLS) by clicking on the download button (in the popup)
4. Now, we want a prediction for several neighbourhoods. Select them on the map using a right-click (the list of selected ones is updated in the left panel). When all relevant neighbourhoods have been selected, select a classifier in the list and click on the button "predict selected neighbourhoods". Prediction results can be exported as tablesheets (XLS) by clicking on the download button (right of the button)

![Screenshot of the cartographic interface of Predihood](doc/predihood-predictions.png)

For the algorithmic interface, an example would be:

1. Choose an algorithm 
2. Tune it as desired
3. Click on "Train, test and evaluate" button. When computing accuracies is done, a table shows results for each environment variable and each list of indicators. 

![Screenshot of algorithmic interface of Predihood](doc/predihood-accuracies.png)


## Tests

Tests are in `predihood/predihood/tests.py` file.

Within a Docker installation, tests can be run as follows:
```
# docker-compose up is running
# docker ps lists running containers to obtain ID_CONTAINER for predihood
docker exec -it <ID_CONTAINER> bash
cd predihood/predihood/
python3 tests.py
```

With a local installation, run the tests using:
```
cd predihood/predihood/
python3 tests.py
```

## Documentation

The documentation of the code is in `predihood/doc/`. It is also available online at [https://nellybarret.gitlab.io/documentation-for-predihood](https://nellybarret.gitlab.io/documentation-for-predihood).

## Datasets

Dataset configuration files are stored in the `predihood/predihood/datasets/` directory. Predihood currently includes two datasets, _hil_ (50,000 neighbourhoods, 550 indicators and 6 environment variables to predict) and _bird-migration_ (769 neighbourhoods, 3 indicators, 1 variable to predict).

Datasets are stored into MongoDB according to the [GeoJSON format](https://geojson.org/).

### Using another dataset

To use an existing dataset (i.e., data already loaded into MongoDB), it is necessary to specify the path to the configuration file for this dataset.

- Using Docker, edit the `predihood/docker-compose.yml` to change the `CONFIG` environment option:

```
CONFIG=datasets/hil/config.json  # to use dataset hil
CONFIG=datasets/bird-migration/config.json  # to use dataset bird-migration
```
And run `docker-compose up` to use the mentioned dataset.

- With the manual installation, run predihood with an argument referring to the configuration file of the desired dataset:

```
python3 main.py datasets/<name-of-dataset>/config.json
python3 main.py datasets/bird-migration/config.json  # to use dataset bird-migration
```

### Importing a new dataset

To import another dataset, it is necessary to follow these instructions: 

1. Create a MongoDB database (called `DATABASE_NAME`), and add the following collections:

        - a collection called COLLECTION_NAME which contains data about neighbourhoods. Each document describes a single neighbourhood under the [GeoJSON format](https://geojson.org/) and includes indicators (short name and value).
        - a collection called "collindic" which contains the indicators used in your neighbourhoods (both short name and full name).

2. Create a CSV file containing human expertise (stored as `predihood/datasets/<name-of-dataset>/expertise.csv`). Each line entails a neighbourhoods' identifier and the expertized value for each variable (see `VARIABLES_VALUES` below).

3. Create a configuration file for your dataset (in `predihood/datasets/<name-of-dataset>/config.json`), in the JSON format, which contains this information:

        - DATABASE_NAME which is the name of your database in MongoDB;
        - COLLECTION_NAME which is the name of the collection in the database that contains information about neighbourhoods;
        - VARIABLES_VALUES is a dictionary which contains (at least in English) the set of variables to predict. For each variable, there are the 'label' (the description of the variable), the 'values' (the values that the variable can handle), the 'low_influence_value' (the value from values which has the least impact on the dataset while filling its missing values) and the 'median_value' (the median value from values, used for filling the missing values of character strings in the dataset);
        - NORMALIZATION corresponds to the unity with which the dataset will be normalised. It can be "None", "population" or "density". We recommend to use density (if possible) to have better results;
		- VARIABLE_REMOVE_LOW_REPRESENTATIVITY corresponds to the variable name used for removing neighbourhoods with lowest representativity when predicting in the cartographic interface.

Examples of these required files are presented in the next part for the _bird migration_ dataset.

To load a new dataset into a Docker image, check the file [mongiris/import-data.sh](https://gitlab.liris.cnrs.fr/fduchate/mongiris/blob/master/import-data.sh) which is automatically run when deploying the container. Data can be loaded either as a MongoDB dump (command `mongorestore`, as shown for database _hil_) or as a sequence of JSON documents (command `mongoimport`, as shown for database _bird-migration_).

### Example of the _bird-migration_ dataset

A fake dataset about bird migration is provided (769 neighbourhoods, 3 indicators, 1 variable). Its configuration file is in  `predihood/predihood/datasets/bird-migration` and its dump files are in `mongiris/mongiris/data/dumps/`. The objective is to predict whether a neighbourhood is suitable for migrating birds to stop by. The single variable accepts 4 values, from _favorable_ to _unfavorable_. The 3 indicators represent the _percent of greens_, the _percent of buildings_ and the _degree of human pressure_ in a neighbourhood.

The following command enables the creation of a MongoDB database with the two required collections. They are already loaded when using the Docker installation.

```
./mongoimport --db=dbmigration -c=collmigration --file=mongiris/mongiris/data/dump-bird-neighbourhoods  # neighbourhoods' collection
./mongoimport --db=dbmigration -c=collindic --file=mongiris/mongiris/data/dump-bird-indicators.json   # indicators' collection
```

The file `predihood/datasets/bird-migration/example-neighbourhood.json` shows an example of neighbourhood (including the value for each of the three indicators). Here is an simplified extract from this file:

```
{
	"_id": "5be32b9df3f0b960b1f8afb2",
	"geometry": {
		"type": "Polygon",
		"coordinates": [
			[
				[
					4.8261667,
					45.7619681
				],
				...
				[
					4.8261667,
					45.7619681
				]
			]
		]
	},
	"type": "Feature",
	"properties": {
		"NAME": "Saint-Georges",
		"CITY_NAME": "Lyon 5e Arrondissement",
		"ID": "693850103",
		"raw_indicators": {
			"percent_greens": 1,
			"human_pressure": 110,
			"percent_built": 99
		}
	}
}
```

The file `predihood/datasets/bird-migration/example-collindic.txt` shows the content of the `collindic` collection (3 documents, one for each indicator):

```
{ "_id" : "5ff1e6f3b0c86c7361a2637e", "short_label" : "human_pressure", "full_label" : "Human pressure on the area" }
{ "_id" : "5ff1e6eeb0c86c7361a2637d", "short_label" : "percent_built", "full_label" : "Percentage of building areas" }
{ "_id" : "5ff1e6e6b0c86c7361a2637c", "short_label" : "percent_greens", "full_label" : "Percentage of greens areas" }
```


Following is an example of `expertise.csv` file (simplified from `predihood/datasets/bird-migration/expertise.csv`), which contains manually expertized neighbourhoods:

```
id_neighbourhood;variable1
693860101;Favorable
693860102;Very favorable
690340801;Favorable
692660201;Favorable
693860104;Not much favorable
693860103;Not much favorable
690340402;Favorable
690340602;Not much favorable
692560101;Very favorable
693860302;Unfavorable
693860303;Unfavorable
```

Here is a commented example of `config.json` file (simplified version of `bird-migration/config.json`):

```
{
  "DATABASE_NAME": "dbmigration",                   # name of the database to connect to
  "COLLECTION_NAME": "collmigration",               # name of the collection containing neighbourhoods
  "VARIABLES_VALUES": {
    "fr": {                                         # variables to predict (in the main language)
      "variable1": {                                # first variable
        "label": "Zone de migration",               # label of the first variable
        "values": ["Favorable", "Défavorable"],     # possible values for the first variable
        "low_influence_value": "Défavorable",       # optional parameter, used for missing values in expertise that are filled in with this value 
        "median_value": "Favorable"                 # optional parameter, used for missing values in expertise (only for character strings variables)
      },
      ...                                           # next variables (if any)
    },
    "en": {                                         # variables to predict (in another language)
      "variable1": {                                # variables should be in the same order as in the main language
        "label": "Migration zone",
        "values": ["Favorable","Unfavorable"],      # values should be in the same order as in the main language
        "low_influence_value": "Unfavorable",
        "median_value": "Favorable"
        },
        ...                                         # next variables (if any)
    },
    ...                                             # next languages (if any)
  },
  "NORMALIZATION": "None",                          # optional parameter, used for normalizing all indicators using the provided indicator (e.g., density)
  "VARIABLE_REMOVE_LOW_REPRESENTATIVITY": "None"    # optional parameter, used for removing neighbourhoods with the lowest representative value for the mentioned variable (e.g., variable1)
}
```


## Note for configuring Predihood in PyCharm

Instead of running Predihood in a console, you can configure your IDE, here for PyCharm. Create a new configuration and set the following parameters:

- Script path: `path/to/predihood/predihood/main.py`
- Python interpreter: add the path to your current Python interpreter
- Working directory: `path/to/predihood/predihood`

You can also run the tests by creating a second configuration and set:

- Script path: `path/to/predihood/predihood/tests.py`
- Python interpreter: add the path to your current Python interpreter
- Working directory: `path/to/predihood/predihood`

# Contributing to Predihood

The following is a set of guidelines for contributing to this repository.

The Code of Conduct is adapted from the Contributor Covenant [Code of Conduct v2.0](https://www.contributor-covenant.org/version/2/0/code_of_conduct.html).

The Contributing sections are inspired from the [Eiffel community guidelines](https://github.com/eiffel-community/.github).

## Code of Conduct 

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

Report to one of the project's maintainer any abuse.

## How to propose changes

Anyone is welcome to propose changes to the contents of this repository by creating a _new issue ticket_ in GitHub. These requests may concern anything contained in the repository: changes to documentation, bug fixes, requests for additional information, additional features, etc.

When posting a new issue, try to be as precise as possible and formulate your arguments for your request carefully. Keep in mind that collaborating on software development is often an exercise in finding workable compromises between multiple and often conflicting needs. In particular, pay attention to the following:
1. What type of change is requested?
2. Why would you like to see this change?
3. Can you provide any concrete examples?
4. Which arguments in favor can you think of?
5. Which arguments against can you think of, and why are they outweighed by the arguments in favor?

Also, keep in mind that just as anyone is welcome to propose a change, anyone is welcome to disagree with and criticize that proposal.

### Closing issues

An issue can be closed by any member of the repository maintainers' team. This can happen in various ways, for varying reasons:
1. Issues without conclusion and no activity for at least 1 month may be closed, as a mere housekeeping activity. For instance, an issue met with requests for further clarification, but left unanswered by the original author, may simply be removed.
2. Issues may simply be rejected if found unfeasible or undesirable. In such cases they shall also be responded to, providing a polite and concise explanation as to why the proposal is rejected.
3. Issues may be closed because they are implemented. Following the successful merging of a pull request addressing an issue, it will be closed.

## How to contribute

While we welcome requests for changes (in the form of issues), we absolutely love ready solutions (in the form of Pull Requests). The best requests are the ones with Pull Requests to go along with them.

Contributions can be made by anyone using the standard [GitHub Fork and Pull model](https://help.github.com/articles/about-pull-requests). When making a pull request, keep a few things in mind.
1. Always explicitly connect a pull request to an issue (as indicated by the issue template).
2. Make sure you target the correct branch. If you are unsure which branch is appropriate, ask in the issue thread.
3. Pull Requests will be publicly reviewed, criticized, and potentially rejected. Don't take it personally.

### Reviewing and merging pull requests

Pull requests can be merged by of the repository maintainers' team.
1. A pull request should be approved by at least two maintainers (including the one doing the merging).
3. When merging, ensure that the description reflects the change. That description should include an issue reference, and should focus on *why* the change was made, to provide the reader with context.

### Closing pull requests

If the author of a pull request has not made any changes or status updates in one month time, any member of the repository maintainers' team should notify the author that they will close the pull request (see below for an example). The member should wait a week after putting the notice before closing the pull request.

  > No one has made an update to this pull request for one month. Please add a status update or we will close this pull request.


# GitHub Octicons

[![npm version](https://img.shields.io/npm/v/octicons.svg)](https://www.npmjs.org/package/octicons)
[![Build Status](https://travis-ci.org/primer/octicons.svg?branch=master)](https://travis-ci.org/primer/octicons)

> Octicons are a scalable set of icons handcrafted with <3 by GitHub.

## Install

**NOTE:** The compiled files are located in `build/`. This directory is located in the published npm package. Which means you can access it when you `npm install octicons`. You can also build this directory by following the [building octicons directions](#building-octicons). The files in the `lib/` directory are the raw source files and are not compiled or optimized.

#### npm

This repository is distributed with [npm][npm]. After [installing npm][install-npm], you can install `octicons` with this command.

```
$ npm install octicons --save
```

## Usage

For all the usages, we recommend using the CSS located in `build/build.css`. This is some simple CSS to normalize the icons and inherit colors.

### Node

After installing `npm install octicons` you can access the icons like this.

```js
var octicons = require("octicons")
octicons.alert
// { keywords: [ 'warning', 'triangle', 'exclamation', 'point' ],
//   path: '<path d="M8.865 1.52c-.18-.31-.51-.5-.87-.5s-.69.19-.87.5L.275 13.5c-.18.31-.18.69 0 1 .19.31.52.5.87.5h13.7c.36 0 .69-.19.86-.5.17-.31.18-.69.01-1L8.865 1.52zM8.995 13h-2v-2h2v2zm0-3h-2V6h2v4z"/>',
//   height: '16',
//   width: '16',
//   symbol: 'alert',
//   options:
//    { version: '1.1',
//      width: '16',
//      height: '16',
//      viewBox: '0 0 16 16',
//      class: 'octicon octicon-alert',
//      'aria-hidden': 'true' },
//   toSVG: [Function] }
```

There will be a key for every icon, with [`toSVG`](#octiconsalerttosvg) and other properties.

#### `octicons.alert.symbol`

Returns the string of the symbol name, same as the key for that icon.

```js
octicons.x.symbol
// "x"
```

#### `octicons.person.path`

Returns the string representation of the path of the icon.

```js
octicons.x.path
// <path d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"></path>
```

#### `octicons.issue.options`

This is an object of all the attributes that will be added to the output tag.

```js
octicons.x.options
// { version: '1.1', width: '12', height: '16', viewBox: '0 0 12 16', class: 'octicon octicon-x', 'aria-hidden': 'true' }
```

#### `octicons.alert.width`

Returns the icon's true width, based on the svg view box width. _Note, this doesn't change if you scale it up with size options, it only is the natural width of the icon._

#### `octicons.alert.height`

Returns the icon's true height, based on the svg view box height. _Note, this doesn't change if you scale it up with size options, it only is the natural height of the icon._

#### `keywords`

Returns an array of keywords for the icon. The data comes from the [data file in lib](../data.json). Consider contributing more aliases for the icons.

```js
octicons.x.keywords
// ["remove", "close", "delete"]
```

#### `octicons.alert.toSVG()`

Returns a string of the `<svg>` tag.

```js
octicons.x.toSVG()
// <svg version="1.1" width="12" height="16" viewBox="0 0 12 16" class="octicon octicon-x" aria-hidden="true"><path d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"/></svg>
```

The `.toSVG()` method accepts an optional `options` object. This is used to add CSS classnames, a11y options, and sizing.

##### class

Add more CSS classes to the `<svg>` tag.

```js
octicons.x.toSVG({ "class": "close" })
// <svg version="1.1" width="12" height="16" viewBox="0 0 12 16" class="octicon octicon-x close" aria-hidden="true"><path d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"/></svg>
```

##### aria-label

Add accessibility `aria-label` to the icon.

```js
octicons.x.toSVG({ "aria-label": "Close the window" })
// <svg version="1.1" width="12" height="16" viewBox="0 0 12 16" class="octicon octicon-x" aria-label="Close the window" role="img"><path d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"/></svg>
```

##### width & height

Size the SVG icon larger using `width` & `height` independently or together.

```js
octicons.x.toSVG({ "width": 45 })
// <svg version="1.1" width="45" height="60" viewBox="0 0 12 16" class="octicon octicon-x" aria-hidden="true"><path d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"/></svg>
```

## License

(c) GitHub, Inc.

When using the GitHub logos, be sure to follow the [GitHub logo guidelines](https://github.com/logos).

[MIT](./LICENSE)  

[primer]: https://github.com/primer/primer
[docs]: http://primercss.io/
[npm]: https://www.npmjs.com/
[install-npm]: https://docs.npmjs.com/getting-started/installing-node
[sass]: http://sass-lang.com/
[Open Iconic v1.1.1](http://useiconic.com/open)
===========

### Open Iconic is the open source sibling of [Iconic](http://useiconic.com). It is a hyper-legible collection of 223 icons with a tiny footprint&mdash;ready to use with Bootstrap and Foundation. [View the collection](http://useiconic.com/open#icons)



## What's in Open Iconic?

* 223 icons designed to be legible down to 8 pixels
* Super-light SVG files - 61.8 for the entire set 
* SVG sprite&mdash;the modern replacement for icon fonts
* Webfont (EOT, OTF, SVG, TTF, WOFF), PNG and WebP formats
* Webfont stylesheets (including versions for Bootstrap and Foundation) in CSS, LESS, SCSS and Stylus formats
* PNG and WebP raster images in 8px, 16px, 24px, 32px, 48px and 64px.


## Getting Started

#### For code samples and everything else you need to get started with Open Iconic, check out our [Icons](http://useiconic.com/open#icons) and [Reference](http://useiconic.com/open#reference) sections.

### General Usage

#### Using Open Iconic's SVGs

We like SVGs and we think they're the way to display icons on the web. Since Open Iconic are just basic SVGs, we suggest you display them like you would any other image (don't forget the `alt` attribute).

```
<img src="/open-iconic/svg/icon-name.svg" alt="icon name">
```

#### Using Open Iconic's SVG Sprite

Open Iconic also comes in a SVG sprite which allows you to display all the icons in the set with a single request. It's like an icon font, without being a hack.

Adding an icon from an SVG sprite is a little different than what you're used to, but it's still a piece of cake. *Tip: To make your icons easily style able, we suggest adding a general class to the* `<svg>` *tag and a unique class name for each different icon in the* `<use>` *tag.*  

```
<svg class="icon">
  <use xlink:href="open-iconic.svg#account-login" class="icon-account-login"></use>
</svg>
```

Sizing icons only needs basic CSS. All the icons are in a square format, so just set the `<svg>` tag with equal width and height dimensions.

```
.icon {
  width: 16px;
  height: 16px;
}
```

Coloring icons is even easier. All you need to do is set the `fill` rule on the `<use>` tag.

```
.icon-account-login {
  fill: #f00;
}
```

To learn more about SVG Sprites, read [Chris Coyier's guide](http://css-tricks.com/svg-sprites-use-better-icon-fonts/).

#### Using Open Iconic's Icon Font...


##### …with Bootstrap

You can find our Bootstrap stylesheets in `font/css/open-iconic-bootstrap.{css, less, scss, styl}`


```
<link href="/open-iconic/font/css/open-iconic-bootstrap.css" rel="stylesheet">
```


```
<span class="oi oi-icon-name" title="icon name" aria-hidden="true"></span>
```

##### …with Foundation

You can find our Foundation stylesheets in `font/css/open-iconic-foundation.{css, less, scss, styl}`

```
<link href="/open-iconic/font/css/open-iconic-foundation.css" rel="stylesheet">
```


```
<span class="fi-icon-name" title="icon name" aria-hidden="true"></span>
```

##### …on its own

You can find our default stylesheets in `font/css/open-iconic.{css, less, scss, styl}`

```
<link href="/open-iconic/font/css/open-iconic.css" rel="stylesheet">
```

```
<span class="oi" data-glyph="icon-name" title="icon name" aria-hidden="true"></span>
```


## License

### Icons

All code (including SVG markup) is under the [MIT License](http://opensource.org/licenses/MIT).

### Fonts

All fonts are under the [SIL Licensed](http://scripts.sil.org/cms/scripts/page.php?item_id=OFL_web).
