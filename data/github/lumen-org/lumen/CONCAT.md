---
title: 'Lumen: A software for the interactive visualization of probabilistic models together with data'

tags:  
  - Probabilistic Modelling
  - Model Criticism, Model Validation
  - Model Understanding  
  - Visual-Interactive Exploration
  - Web-Interface

authors:
  - name: Philipp Lucas^[corresponding author]
    orcid: 0000-0002-6687-8209
    affiliation: 1
  - name: Joachim Giesen
    affiliation: 2

affiliations:
 - name: Institute of Data Science, German Aerospace Center
   index: 1
 - name: Friedrich-Schiller-University Jena
   index: 2

date: 19th April 2021
bibliography: paper.bib
---

# Summary

Research in machine learning and applied statistics has led to the development of a plethora of different types of models.
*Lumen* aims to make a particular yet broad class of models, namely, probabilistic models, more easily accessible to humans. 
*Lumen* does so by providing an interactive web application for the visual exploration, comparison, and validation of probabilistic models together with underlying data. 
As the main feature of *Lumen* a user can rapidly and incrementally build flexible and potentially complex interactive visualizations of both the probabilistic model and the data that the model was trained on. 

Many classic machine learning methods learn models that predict the value of some target variable(s) given the value of some input variable(s).
*Probabilistic* models go beyond this point estimation by predicting instead of a particular value a probability distribution over the target variable(s).
This allows, for instance, to estimate the prediction's uncertainty, a highly relevant quantity.
For a demonstrative example consider a model predicts that an image of a suspicious skin area does _not_ show a malignant tumor.
Here it would be extremely valuable to additionally know whether the model is sure to 99.99% or just 51%, that is, to know the uncertainty in the model's prediction.

*Lumen* is build on top of the [*modelbase*](https://github.com/lumen-org/modelbase) back-end, which provides a SQL-like interface for querying models and its data [@Lucas:2021:modelbase].

# Statement of need

A major challenge for both the application and development of machine learning/modelling methods is their accessibility to a human analyst, that is, the amount of hurdles that one must overcome to practically make use and benefit from them.
*Lumen* aims to improve accessibility of probabilistic machine learning models with respect to multiple aspects as follows:

_Model Building:_
Building a statistical/machine learning model is often an iterative, analyst-driven process.
This is particularly true for the field of probabilistic programming, a modelling approach where the analyst explicitly declares the likelihood of the observed data as a probability density function.
The analyst typically starts with an exploration of the data.
Based on insights gained from data exploration and on the analyst's domain knowledge, the analyst creates an initial simple model involving only some data.
Subsequently, this model is iteratively made more complex [@Gelman:2013; @Gabry:2019] until it meets the expert's goals.
In particular, the model must be validated after each iteration.
*Lumen* supports this model building process by (i) enabling visual-interactive data exploration, (ii) supporting model validation by means of a visual comparison of data queries to semantically equivalent model queries, and (iii) enabling a direct comparison of model iterates.

_Debugging:_
Even for a machine learning expert it may be hard to know whether a model has been trained on the data as expected.
Possible reasons for artifacts in a model include an inappropriate application of the machine learning method, implementation bugs in the machine learning method, and issues in the training data.
Direct visual inspection of the probabilistic model provides an approach to model debugging that enables the analyst to literally spot model artifacts that may cause degrading performance.
Classical approaches to validation would rely on aggregating measures like information criterions or preditictive accuracy scores.

_Education:_
By its intuitive visual representations of models, *Lumen* aims to promote understanding of the underlying modelling techniques. 
For instance, the effect of varying a parameter value for a modelling method on the probabilistic model can be observed visually rather than remaining an abstract description in a textbook.
Similarily, the differences between models/model types can be visually illustrated  by plotting them side by side.
Also, probabilistic concepts such as conditioning or marginalization, which are often difficult to grasp, can be tried out interactively, providing immediate feedback.

# Software

*Lumen's* interface is inspired by the academic Polaris project and its commercial successor Tableau [@Stolte:2002]. 
However, while Polaris/Tableau is for _data only_, *Lumen* provides a uniform visual language and interactions for both data and probabilistic models.
\autoref{fig:LumenUI} shows an screenshot of *Lumen* to illustrate the user interface. 
The Schema panel (left) contains the random variables of the probabilistic model that the user has currently selected.
Users can drag'n'drop variables onto the visual channels of the Specification panel (middle-left).
This reconfigures the currently active visualization on the dashboard (middle to right), triggers execution of corresponding data and model queries, and finally updates and re-renders the visualization.
To foster comparison of multiple models (for instance from different classes of models or from iterates of an incremental model building process) Lumen allows users to create as many visualizations of as many models as desired.
All visualization support basic interactions like panning, zoom, or selections and are resizable as well as freely movable on the dashboard.

![The Web-based interface of *Lumen* displaying a variety of visualizatons as created in the process of incrementally building a probabilistic model on the socio-economic ALLBUS data set [@Allbus:2016]: 
(1) Marginal data density. 
(2) Marginal model density (pink) versus observed data density (grey). 
(3) Both plots show the same queries but from (a) to (b) the underlying model was improved to better capture the correlation of the `income` variable and the `sex` variable.
Again, data are shown as histograms and model densities as line plots. 
(4) Connected dots show the model's point predictions of `income` given `age` and `sex`.
Marks in the background as well as the marginal plots at the side represent observed data.
(5) Similar to (4) but visualizing the model's predictions of `income` as well as of `happiness` given `age` and place of origin (`eastwest`). 
Again, the background marks show observed data.\label{fig:LumenUI}](joss/example.png){ width=95% }

While *Lumen* handles all user facing aspects (such as visualizations and interactions) most computational aspects (such as execution of model or data queries that are triggered by a user interaction) are delegated to a dedicated back-end. 
The back-end is implemented in the *modelbase* project [@Lucas:2021:modelbase].
This separation follows a classic client-server architecture where *Lumen* is the web-client and *modelbase* the web-service.
For the standard usage scenario both client and server would be installed locally on the same machine. 
However, they can, of course, also be separated on different machines across a network.

*Lumen* is model-agnostic in the sense that it can be used with models of any class of probabilistic models as long as this model class implements the common, abstract API in the *modelbase* back end. 
The API essentially requires that a model class 

 * contains only quantitative and categorical random variables, i.e. Lumen has no native support for images, time series, or vector-valued random variables, 
 * supports marginalization of random variables, i.e. the operation to remove/integrate out any subset of random variables of the model, 
 * supports conditioning of random variables on values of its domain, i.e. the operation to fix the value of random variables to particular values, and
 * supports density queries, i.e. the operation to ask for the value of the model's probability density function at any point of its domain.

In fact *Lumen* does not depend on any specific properties of a particular model class and we regard this genericity as one of *Lumens* major features. 
Among the model classes that we have used *Lumen* with are Sum-Product-Networks [@Poon:2011; @Molina2019:SPFlow], Condional-Gaussian Distributions [@Olkin:1961:CG; @Nussbaum:2020:paper], Probabilistic Progams based on PyMC3 [@Salvatier:2016:PyMC3], and Kernel-Density-Estimators [@Parzen:1962:KDE; @SciPy:2020].

# Acknowledgements

We thank Andreas Goral, Jonas Aaron Gütter, Laines Schmalwasser, Julien Klaus and Christian Lengert for their steady and patient interest in trying out Lumen, for their valuable feedback and our discussions, as well as for the features they contributed to Lumen.
Philipp Lucas was partially supported by Stiftung der Deutschen Wirtschaft (sdw). 

# References


# What is `lumen`?

`lumen` is an interactive web-application for the visualization and exploration of probabilistic machine learning models. 
Its main feature is the ability to rapidly and incrementally build flexible and potentially complex visualizations of both probabilistic machine learning models and the data these models were trained on.

# Using `lumen`

`lumen` aims to make a particular class of machine learning/statistical models, namely *probabilistic* models,  more easily accessible to humans. 
A probabilistic model models a set of target variables by means of a probability distribution.
That is, different to many classic ML methods which predict a particular value of the target variable(s), probabilistic models instead capture the distribution of the target variables. 
`lumen` lets you 'see' your model, understand how it performs, where it 'fails', and compare this to previous versions of the model or alternative models. 

## Manual
A manual-style description of the UI, the visual encodings, `lumens` usage, its features, and available interactions is available [here](doc/manual.md).

## Walk-Through
A walk-through-style introduction to `lumen` is available [here](doc/walkthrough.md). It demonstrates some of the feature for exploration of probabilistic models in `lumen`.

![`lumens` user interface displaying a variety of visualizations of a probabilistic model on a socio-economic data set](doc/img/example_raw_processed.png)

In particular `lumen`  lets you:

 * plot any marginals of your models. 
Here, marginal means not just 1d but also higher dimensional marginals. 
Studying multiple multi-variate 'slices' of a model may help to understand the interactions between the variables.
We believe this also helps you to fix model degrading artifacts. 
Such artifacts may indicate a problem in your model specification, model parameterization or possibly a bug in the machine learning algorithm of your model.
 * plot the model marginals together with data marginals. 
This lets you directly check the models fit to data.
 * plot predictions of your model along side corresponding data aggregations. 
This lets you understand its predictive behaviour, and also compare it observed quantities.
 * combine any of the above 'layers' into a single visualization.
 * change visualizations by flexibly assigning variables/data attributes to visual channels.
 * create as many of these visualizations side by side on an virtually infinite canvas. 
This lets you compare various stages of a model, compare different modelling approaches, and get a better overall understanding by combining many different visualizations of the same model.

## Augmenting Probabilistic Programming

Probabilistic programming language (PPLs), such as PyMC3, BLOG, or Stan, provide a framework to define probabilistic models by explicitly declaring the likelihood of the observed data as a probability density function.
The analyst typically starts with an exploration of the data.
Based on insights gained from data exploration and on the analyst's domain knowledge, the analyst creates an initial simple model involving only some data.
Subsequently, this model is iteratively made more complex until it meets the expert's goals.
In particular, the model must be validated after each iteration.
`lumen` supports this model building process by 
(i) enabling visual-interactive data exploration, 
(ii) supporting model validation by means of a visual comparison of data queries to semantically equivalent model queries, and 
(iii) enabling a direct comparison of model iterates.

## Model Debugging

Even for a machine learning expert it may be hard to know whether a model has been trained on the data as expected. 
Possible reasons for artifacts in a model include an inappropriate application of the machine learning method, implementation bugs in the machine learning method, and issues in the training data. 
Direct visual inspection of the probabilistic model provides an approach to model debugging that enables the analyst to literally spot model artifacts that may cause degrading performance.
Classical approaches to validation would rely on aggregating measures like information criterions or predictive accuracy scores.

## Education / Teaching

By its intuitive visual representations of models, Lumen aims to promote understanding of the underlying modelling techniques.
For instance, the effect of varying a parameter value for a modelling method on the probabilistic model can be observed visually rather than remaining an abstract description in a textbook.
Similarly, the differences between models/model types can be visually illustrated by plotting them side by side.
Also, probabilistic concepts such as conditioning or marginalization, which are often difficult to grasp, can be tried out interactively, providing immediate feedback.

## Data-only exploration

You don't do any Machine Learning but simply would like to conveniently browse, explore, and compare tabular data? 
`lumen` is the right place for you too!
This is not what `lumen` was built for originally, but regard it as your 'free lunch'.

---

# Installing `lumen`

This explains how to install and configure `lumen` and its dependencies.

Note that `lumen` is build on top of the `modelbase` [back-end](https://github.com/lumen-org/modelbase), which provides a SQL-like interface for querying models and its data.

## Requirements

* `lumen` is a web application that requires access to a web-service instance of the Python3-based `modelbase` backend.
`lumen` allows a user to interactively compile data/model queries and visualize the queries results. `modelbase` does the computation and actually answers the queries. 
You can get `modelbase` [here](https://github.com/lumen-org/modelbase) where you also find information on how to set it up and run it as a web-service.

* `lumen` and `modelbase` need to be configured correctly with 'matching' settings. By default (both run locally on the same physical machine) this is the case and you do not need to change these settings:
  * hostname set in the configuration of `lumen` must match the actual hostname of `modelbase`.
  * port must match
  * protocol must match (http or https)

* `lumen` allows you to explore the models and data that are hosted by the `modelbase` backend. 
You can use the `modelbase` Python package to (1) train/create models from data, and then (2) host them by an instance of the `modelbase` web-service.
See the [documentation and introductory jupyter notebooks](https://github.com/lumen-org/modelbase) in the `doc` folder for more information. 
Also, a number of example models are created during the setup process of `modelbase` for your convenience.

## Setup

Clone/download this repository into a folder `<path>` of your choice.
## Updating it

Just pull/download the lasted branch/version you'd like.

## Running it

1. make sure the `modelbase` backend is running and hosting the models that you'd like to explore. 
2. it's dead simple: Open `<path>/index.html` in your browser. If everything is fine you should now see a model dialog that lists the available models. Select one and start exploring it!

Notes:
 * Using *chrome/chromium* as a browser is recommended, since it provides the best performance from our experience. 

---

# Trouble Shooting

If you have any trouble using `lumen`, need some additional explanation, or even just want to provide some feedback, please don't hesitate to contact us at [philipp.lucas@dlr.de](philipp.lucas@dlr.de).
If you encounter any bugs you can also [submit an issue](https://github.com/lumen-org/lumen/issues/new/choose).

## Typical Issues

#### When open `lumen` in my browser I get the error message: "Could not load remote model from server!"
 
 1. Confirm that the backend server actually running
 2. Check the developer console log of the browser where you are loading the front-end. If it shows something like:
 
     ```Failed to load http://127.0.0.1:5000/web-service: Response to preflight request doesn't pass access control check: The 'Access-Control-Allow-Origin' header has a value 'null' that is not equal to the supplied origin. Origin 'null' is therefore not allowed access.```
 
 Then your probably run into some CORS issue because you serve the file directly from the file system, instead from a webserver running locally. See here for the issues:
   * [problem description: answer 1, point 2 ](https://stackoverflow.com/questions/3595515/xmlhttprequest-error-origin-null-is-not-allowed-by-access-control-allow-origin)
 
 Solutions:
   * serve it from a local web-service (preferred)
   * [disable CORS control in chrome](https://stackoverflow.com/questions/3102819/disable-same-origin-policy-in-chrome) (kind of hacky)
 

#### I get the error message: "Could not load remote model 'XXXX' from server 'XXXX' !"
  1. Confirm that the backend server is actually running
  2. Did the backend server load the particular model that you are trying to retrieve? Loaded models are listed in the terminal output of the backend server on its start up.

---

# Contributing

You wanna contribute? Awesome! Let's get in touch: [philipp.lucas@dlr.de](philipp.lucas@dlr.de) !

### Development Setup

This is only for you, if you want to contribute to the project.

1. Do the steps as described in the Setup section above.
2. Install [node-js](https://nodejs.org/en/download/). For questions refer to the [getting started guide](https://docs.npmjs.com/getting-started/what-is-npm).
3. Update npm (part of node-js): `sudo npm install -g npm`
4. Install all npm-dependencies as provided by the projects `package.json`:
    * run from `<path>`: `npm install`
---

# Contact

For any questions, feedback, bug reports, feature requests, spam, rants, etc please contact: [philipp.lucas@dlr.de](philipp.lucas@dlr.de)

# Copyright and Licence

© 2016-2021 Philipp Lucas (philipp.lucas@dlr.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
# Walk-through: a step-wise introduction to `lumen`

# Overview 

This document shows a step-by-step example how to use `lumen` to analyze and explore a given probabilistic model.
Specifically we will have a look a Conditional Gaussian (CG) distribution fitted to the popular *![mpg cars data set](https://github.com/hadley/data-fuel-economy)*.
The data set contains yearly fuel economy data provided by the EPA (http://www.fueleconomy.gov/feg/download.shtml) for the years from 1978 to 2008.
Each data entry corresponds to a car model, and for each car model a number of attributes have been recorded. 

For a more comprehensive explanation of the features and usage of `lumen` please refer to the ![manual](manual.md).

`lumen` requires an instance of `modelbase` as its backend that hosts all models and executes any queries to these models.
If you follow the ![default instructions to install `modelbase`](https://github.com/lumen-org/modelbase) a Conditional Gaussian (CG) model on the mpg cars data set is trained as part of the setup process and should be available right away for exploration when you run `lumen`.
So, if you like you can follow this guide right away and try it out yourself!

In this walk-through we will touch the following aspects

 * learn: understand concepts of Conditional Gaussian (CG) distributions
 * explore: discover and explain interesting structure of the model
 * validate: spot, identify, and explain data or model artifacts

All figures shown here are screenshots of Lumen.
However, to reduce the required screen space we cropped the full screenshot to only show the lastly modified visualization.
Note that the dashboard in Lumen may hold many visualizations at once.

# 01: Start Up Applications

 1. Start the`modelbase` web-service by executing `webservice.py` of the `modelbase` backend.
 2. Run `lumen` by opening the `index.html` in your browser (best is chrome-based). 
 You should see a list of available model (the models provided by the backend) as follows:

 ![walkthrough 01 - startup](img/walkthrough-01-startup.png)

 3. Click on `mpg_cond_gauss` to load that model. 
 An empty specification and empty visualization is now shown and we can start to create our visualizations.
 However, lets first have a look at the Schema panel on the left:

![walkthrough 02 - schema](img/walkthrough-02-schema.png)

 As we can see the model contains a total of 7 random variables, and since we use a Conditional Gaussian (CG) model all of them are backed by observed data and all are 'random' variables modeled by the probabilistic model.

 
 Here is what the variables mean:

  * `year`: the year the vehicle was sold
  * `mpg_city`: measures the average miles per gallon in the city for the vehicle
  * `mpg_highway`: measures the average miles per gallon on the highway for the vehicle
  * `displacement`: measures the engine displacement in cubic centimeters  
  * `transmission`: describes what kind of transmission that car has as categories of `auto`, `lock-up`, and `manual`.
  * `cylinder`: describes how many cylinders the engine has categorized to `few`, `medium`, and `many`.
  * `car_size`: describes the physical size of the car categorized to `small`, `midsize`, and `large`.

# 02: Exploring the Conditional Gaussian Model

Let's have a look at the uni-variate marginal distributions of `mpg_highway`.
To do so, we drag the variable onto the X-Axis shelf. 
Since by default only data facets are activated, this results in a visualization of the data only. To also show model marginals, we check the respective facet `model - marginals`.

![walkthrough 03 - marginals](img/walkthrough-03-marginals.png)

The visualization shows that the data (grey) and the model (pink) fit quite neatly (for this marginal).
In particular, the Conditional Gaussian (CG) model succeeded in capturing the multi-model structure of the data (here two modes, that is local maxima, are visible).

A Conditional Gaussian (CG) model effectively models the joint distribution over all variables (both quantitative and categorical) by fitting a multivariate Gaussian distribution for each combinus to explain the structure of the data and potentially identify semantic clustersation of categorical values. 
Hence, the CG distribution could help . 
Here, we raise the question which of the categorical attributes are linked to the multi-modality?

We can simply try it out. Let's assign the categorical variable `transmission` to the X-Axis-shelf too.
Make sure to drop it as the first element in the shelf.
This creates a hierarchical axis and we get a marginal distribution for each of the values of `transmission`:

![walkthrough 04 - modality](img/walkthrough-04-explaining_modality.png)

We can see that there is a strong correlation between `transmission` and displacement, as the three marginal plots, especially the one for manual transmission, are quite different. 
However, it does not yet explain the modality of the marginals.

Hence, we also add the variable `cylinder` to our specification.
This time, we drag it onto the Y-Axis shelf. 
Also, we enable the facet `model - density` and `data - density` but disable the `data points` facet.
This results in a table-like visualization. 
Let's quickly explain what we see: 
There is a 'central' 3x3 array of plots that show the (marginal) distribution of `displacement` conditioned on all combinations of value for  `cylinder` and `transmission`. 
Note the middle row carries almost no probability and hence resembles a extremely flat curve.
On the very top there are three more marginal distribution plots, one for each column of the 3x3 array of plots.
Finally, the horizontal bar plots visualize another marginal, namely, the (categorical) distribution of `cylinder` given a particular value of `transmission`.

![walkthrough 05 - modality](img/walkthrough-05-explaining_modality.png)

There is a number of interesting observations here:

First, the newly included variable `cylinder` appears to explain the multi-modality: In each column of the table-like visualization we see how the conditional marginal distribution (top) decomposes into three uni-modal parts (rows of the table).

Second, there are almost no cars with many cylinders (middle row).
This may either indicate a poor choice for data preprocessing (here we have categorizes the originally numerical value for `cylinders`), or simply an interesting observation, namely, that there are few car models with many cylinders.

Often the exploration process leads to some insight that, once found, can be visualized explicitly very nicely. 
Here, `cylinder` appears to be an important variable, yet there is some doubts about its validity. 
Let's create a new plot (select the model in the toolbar -> Load Model -> Go!). 
This time we encode `cylinder` as color but plot two quantitative variables (`mpg_city` and `displacement`) on the positional channels (X-Axis and Y-Axis-shelf).
Also, we enable the `data-aggregation` facet. 
It will indicate the average `displacement` and `miles per gallon`.

![walkthrough 06 - exploration](img/walkthrough-06-exploring_further.png)

Note that this is a visualization of the *data* only. 
Yet, we see a very distinct clustering according to the value of `cylinder`.

Let us activate the model facets for aggregation, marginals and density. 
The resulting visualization shows a very good fit of the models distribution (pink contour plot in the background) of `displacement` and `mpg_city` to the observed data (marks with white strokes).
Additionally, the marks of the aggregation facets (marks with black stroke) are very close, confirming the good fit also on a higher level of aggregation.

![walkthrough 07 - exploration](img/walkthrough-07-exploring_further.png)

# Manual

This document provides an manual-style introduction to `lumen` and its user interface.

In a few sentences `lumen` could be summarized as follows:

`lumen` allows you to visually explore probabilistic models and their data.
You, the user, assign attributes of the data / random variables of the model to visual variables by drag-and-drop interactions to specify what part of the model/data you would like to see and how this model/data is visually encoded. 
A number of combinable 'semantic layers' allow you to visualize different aspects of the model and data.

![`lumens` user interface](img/example.png)

# Overview

These are the five main components and their most important function, see also the screenshot.

 1. Toolbar: Load models and create new visualizations
 2. Schema panel: Shows the data attributes/random variables of the model that is shown in the active visualization.
 3. Specification panel: Lets you modify the assignment of data attributes/random variables of the model to visual variables in order to change the active visualization.
 Also, it lets you choose which semantic facets to enable.
 Facets are layers in a visualization that show a particular aspect of model and/or its data.
 4. Dashboard: A pannable container that holds all the visualizations that you have created on an virtually infinite canvas.
 5. Visualization (contained in the dashboard): The visualization as configured in the specification. 
 No limits on how many visualizations you can have at once on the dashboard.

To keep things simple each visualization is associated to exactly one model and its data.
Put differently: you cannot mix multiple models in one visualization.
However, of course, you can create multiple visualizations, each of a different model.

At several points above it said 'active': 
At any point there is exactly one 'active' visualization.
You can recognize the active visualization by its darker frame around it, compared to the rest of the visualizations.
The specification panel always shows the visual configuration of the active visualization.
Hence, when you modify the specification, you always modify the currently active visualization.
Also, the schema panel shows the variables / attributes of the single model that is associated with a visualization.
Hence, the schema shows the variables / attributes for the _active_ model.
As a user you can change the active visualization and corresponding active model by simply clicking on the desired visualization in the dashboard.

The manual is structured as follows:

 * [Probabilistic Models and Random Variables](#probabilistic-models-and-random-variables): explains what a probabilistic model in `lumen` is and how random variables and data attributes are related.
 * [User Interface and Components](#user-interface-and-components): explains the components of the user interface and how to use them.
* [Interaction](#interaction): Briefly lists all available interactions.

---

# Probabilistic Models and Random Variables
## Probabilistic Models in Lumen

The general probabilistic model explorable in `lumen` can be written as

p(X, U | Y, T)

where

 * X is the set of observed random variables,
 * U is the set of latent random variables,
 * Y is the set of observed covariates ('input' variables), and
 * T is the set of latent covariates (typically (hyper-)parameters of the model)

The basic approach of `lumen` to model understanding is that it lets the user visualize any marginal/conditional sub-models of p(X, U | Y, T) together with the 'equivalent' data slices.
Note, this requires that one can effectively compute such marginal and conditional distributions.
Hence, the appropriate model class needs to be implemented in the `modelbase` backend.

For many model classes the above density function simplifies quite a bit, as some of X,U,Y, or T are empty.
For instance, in conditional Gaussian (CG) distributions, all variables are observed random variables, and any conditional Gaussian distribution simplifies in our notation to p(X).

## Random Variables and Data Attributes

`lumen` builds on the idea that random variables and corresponding data attributes can be treated the same way for the purpose of interactive visualization.
Hence, at most places we make no difference between random variables and the corresponding data attributes.
This may seem confusing at first, however, in fact it is quite natural when you consider our scenario where random variables (often) directly correspond to observed data.

This section aims to clarify on random variables / data attributes, how they are represented in `lumen`, and what properties they expose. 
In fact we use the two terms (random variable and data attribute) almost interchangeably.
However, we prefer random variable (or variable for short) to highlight the modeling aspect:

 * __data type__: A variable may either be of the numerical type or of the categorical type. 
 For instance, a variable that measures the ![sepal length of iris flowers](https://en.wikipedia.org/wiki/Iris_flower_data_set) naturally is of the numerical type, while a variable that encode the species of an iris flower, naturally is of the categorical type.
 * __latent / observed__: A variable may or not be observed. 
 If it is observed, then corresponding data is available and you may visualize the data and aggregations of data alongside the variable. 
 If it is not observed (latent), then you can, of course, not do that.
 * __variate / covariate__:  Some variables in a model only serve as an input but are not modeled as part of the model and are not part of the output/prediction of an model.


---

# User Interface and Components

# Toolbar

The toolbar is located on the top edge of the UI.

![Toolbar in Lumen's UI](img/toolbar.png)

## Loading models / creating new visualizations

Most importantly you can create new visualizations using the toolbar.
Go to the drop-down menu on the left, select among the available models and then hit "Go!" to get a brand new, empty visualization of the selected model on the dashboard.
The newly created visualization is automatically activated, that is, the schema now represents the newly loaded model and the specification configures now the  newly created visualization. 

## Clone Button

The clone button duplicates the currently active visualization. 
On click a new visualization with identical content is created, activated, and can be used and modified on its own.

## Clear button

The clear button removes all assignments of random variables to shelves in the specification.

## Query button

Hitting the query button triggers a re-computation and re-rendering of the active visualization. 
It is comparable to forced refreshing a page that skips the cache. 
In case something in `lumen` just went wrong, this may help ;)

## Details

The Details button toggles a panel on the right that shows some more details about the active model and visualization.
For instance, it contains textual description associated to the active model and allows you to download the data associated to the enabled facets as shown in the active visualization.

## Config

The Config button toggles a panel for advanced configurations to change colors, opacity, strokes and much more.
While it is generally safe to use, it is more of a developer tool due to its complexity and missing documentation in the tool.

---

# Schema Panel

The schema lists all variables of the model's variables / data attributes.
It groups the variables by their scale type, namely into 'quantitative' and 'categorical'.
Additionally, it provides information whether it is a observed or latent variable, and whether it is a random variable (a variate variable modeled as being randomly distributed by model) or a independent variable (a covariate, a variate that just serves as a required input to the model).
See also [Random Variables and Data Attributes](##Random-Variables-and-Data-Attributes).

![Schema panel of Lumen's UI](img/schema.png)

In this example, the schema panel lists the variables for a model with name 'Iris_cond_gauss'. 
The model has five variables, four of them quantitative (`sepal_length`, `sepal_width`,`petal_length`,`petal_width`) and one categorical (`species`).
As the the model does not have any latent (non-observed) variables, all variables are labelled as 'observed'.
Also, all of the listed variables are actually modelled by the model, hence they have the label 'random'. 

---

# Specification Panel

The specification consists of a number of so-called shelves (top: *X-Axis* to *Size*), the semantic facet selector (middle: :*Facets*), and some advanced configurations (bottom: *Config*).
See also the example below.

## Shelves

The shelves represent *visual channels*, and allow you to configure two things at once:

 1. what part of the data / what sub-model you want to visualize
 2. how you want to visualize it

In short these two things are done as follows

 * _what_: by assigning variables to the specification at all.
 Variables that are not in the specification will be removed from a model before the model is visualized (-> marginalization and conditioning of models)
 * _how_: by assigning variables to specific shelves.
 Each shelf has its own semantic as a visual variable, as explained in the following.

### X-Axis and Y-Axis shelf

These two shelves represent the positional channels. 
In the above example `sepal_width` is assigned to *X-Axis*, hence, in the visualization the x-axis encodes the values of `sepal_width`. Quite simple, eh?
Not surprisingly, it works the same for *Y-Axis*.

### Color, Shape and Size shelf

This is quite straight forward too: 
Assigning a variable here, will cause the visualization to use color/shape/size to encode the value of this variable.

Have a look at the following 5-dimensional (!) visualization:

![Advanced Example for Specification](img/advanced_example.png)

It encodes all of the information in the iris data (no model is shown at all). 
Each visual mark encodes one data point.
Obviously, some information is easier to recognize than other, but that's a trade-off we always have to deal with.

### Filter shelf

Using the filter shelf you can further specialize the data and model to show.
It allows you to restrict the values that variables may take. 
After assigning a variable to the shelf, you can click on the shelf item to open a modal (pop-up) dialog.
Here, you can restrict the interval (for quantitative) and set (for categorical) of allowed values. 

![Filter Shelf in Detail](img/filter.png)

In the shown example the value of `petal_width` is restricted to values small than 1.85 (see the modal dialog), and the values for `species` are restricted to 'setosa' and 'virginica' only excluding 'versicolor' (not shown explicitly, but note that green *species* field in the filter shelf). 
As you can see in the visualization, under these conditions, 'setosa'  and 'virginica' can be separated almost perfectly just using the single variable `sepal_length`.

### Details shelf

In the beginning of this section it said that shelves serve two purposes.
All the shelves described so far serve both purposes: select variables to include in the visualization (the 'what?') and assign visual channels for their visualization (the 'how?').
The Details shelf, however, only adds a variable to the visualization (what), but *without* assigning a visual channel (how).
Read on in the very next section on roles of variable usages.

## Aggregating and Grouping of Variables

Actually, for each usage of a variable (that is any assignment of a variables to a shelf in the specification) you configure a third aspect, namely the *role* of that variable when aggregating data or model:

For each variable `X` you choose whether:

 * `X` is grouped by (shown as a light blue fill of the variables usage in shelf of its specification), or
 * `X` is aggregated/predicted (shown as a  light yellow fill of the variable usage in shelf of its specification).

Note that the details shelf and the configuration whether to group or aggregate is only relevant for the *aggregation facets*, see below.

The distinction between grouping and aggregation is easy to understand with some examples.
The screenshot below shows an arrangement of 2x3 visualizations.
*All specifications are identical except for the assigned roles of the variables.*
The bottom right visualization is the active one here.

![Roles of Variables in Lumen](img/roles_of_variables.png)
 
 * Top-left: Aggregate (predict) all three variables. The shown dot represents the single point that is most likely for *sepal_length*, *sepal_width*, and *species*.
 * Top-middle: aggregate (predict) *sepal_width* and *species*, but group by *sepal_length*.
 For each value of the grouping an aggregation/prediction is made.
 Here, *sepal_length* is grouped into 5 values equidistantly chosen along its observed extent of values.
 * Top-right: Aggregate (predict) *sepal_length* and *species*, but group by *sepal_width*.
 * Bottom-left: Aggregate *sepal_length* and *sepal_width* and group by *species*. Hence, the most likely value for the two quantitative variables is computed and shown *per category value* of species.
 * Bottom-middle: Aggregate (predict) *species* grouped by both *sepal_width* and *sepal_length*. For each raster element the most likely species is shown.
 * Bottom-right: Given *sepal_length* and *species* predict *sepal_length*.

To change wether a variable is grouped by or aggregated/predicted simply click on the tiny icon that shows up on the right when hovering on a variable in a shelf.

## Visual Defaults

The specification panel in `lumen` gives a lot of flexibility to you, the user, which may be a bit overwhelming at first. 
But it's quick to get used to, don't worry.

Also, to make it easy for you to keep facets apart and understand plots, `lumen` already took many 'good default decisions' for you too to ensure that plots are consistent and easy to understand:
The default visual encodings are as follows:

 * model-related marks are pink, on a pink scale, and/or square
 * data-related marks are grey, on a grey scale, and/or circle
 * aggregation marks have a black stroke and high opacity
 * data points / samples from models have a white stroke and low opacity
 * within one visualization all encodings are comparable, which is, they are on the same scale

Note that many of the defaults can be overridden by using the specification.

## Facets

Facets act as 'semantic layers', that is, they represent different aspects.
`lumen` provides eight facets organizes in two columns and four rows:

### Columns
 
 * model column: facets here represent aspects of the _model_
 * data column: facets here represent aspects of the _data_ that the model was trained on

### Rows

 * aggregation: Adds marks for aggregated (summarized, predicted) values to the visualization.
 * data points: Adds marks for data points (data) and samples (drawn from the model) to the visualization.
 * marginals: Adds marginal distribution plots to the visualization. 
 * density: Adds density distribution plots to the visualization

For an illustration, see this 2-by-4 arrangement of eight individual visualizations. 
All have the identical specification, however, each has exactly one facet only  activated.
Rows are data and model, and columns are aggregation, data points, marginals and density facets, respectively.

![Semantic Facets](img/facets.png)

 
## Tabular Visualizations: X-Axis and Y-Axis shelves revisited
 
The above description of the X-Axis and Y-Axis shelves only told half the story,
because you can actually have multiple variables assigned to each of these shelves at the same time.
This allows you to create 'tabular arrangements' of plots within a single visualization.

Let's have a look at an example:

![Positional Shelves Revisited 01](img/positional_shelves_revisited_01.png)

For the left visualization left we created a scatter plot of model samples drawn from the probabilistic model by dropping `sepal_length` on the x-axis shelf and `petal_length` on the y-axis shelf.
For the right visualization (the active one), we cloned the visualization, dragged the `species` variable from the schema, and dropped it on the x-axis shelf. 
The visualization now contains three scatter plots, namely one for each value of `species`, instead of only one.
Note how all individual plots share both the x and y-axis.
Here, `species` is used to group/split the single plot into individual ones, creating an additional hierarchical x-axis for `species` on the bottom.

Instead of a creating a hierarchy you can also just 'add' (concatenate) another variable to the horizontal or vertical layout.
For the following visualization we dropped `petal_width` to the y-axis shelf (and resized the plot).
Notice how there is _no_ hierarchical axis and instead `petal_width` is just added next to `petal_length`.

![Positional Shelves Revisited 02](img/positional_shelves_revisited_02.png)

How do you specify whether to create an hierarchical axis or just concatenate?
Here, we reuse the assignment of variables usage to 'aggregating' or 'grouping', see Section [Aggregating and Grouping of Variables](##Aggregating-and-Grouping-of-Variables).
In short, 'blue' shelf items create hierarchies and 'yellow' shelf items concatenate axis. 
You can swap between 'blue' and 'yellow' by hovering on a shelf item and clicking the yellow/blue button.

There is no explicit limit on how many variables you may add to the positional shelves. 
Here is two more examples that illustrate useful applications.
The following visualization contains all data and model marginals of the `iris_cond_gauss` model in one visualization:

![Positional Shelves Revisited 03](img/positional_shelves_revisited_03.png)

As another example the following visualization shows several facets for `age` over `fare` for all combinations of the four variables `sex`, `embarked`, `passenger class(Pclass)`, and `Survived`. 
This kind of a plot is often referred to as 'small multiples'.

![Positional Shelves Revisited 04](img/positional_shelves_revisited_04.png)

---

# Interaction

In Lumen the mouse is the primary interaction device. 
It's made with the idea to be (hopefully) intuitive, and most of the available interaction have already been used or mentioned on in the descriptions above.

Anyway, here is a concise list of what you can do:

## Working with Variables
* assign variables to visual channels or other shelves by dragging and dropping them between the various shelves in the Schema panel and the Specification panel
* remove the assignment of a variable to a shelf by dragging that variable and dropping it 'anywhere'
* reorder variables within a shelf by dragging and dropping them
* open modal dialog to further configure a variable usage in the specification by clicking on it
## Facets
* enable and disable a facet by checking or un-checking it in the specification panel
## Handling Visualizations in the Dashboard
* resize visualizations by dragging the edges of it
* move visualization on the Dashboard by dragging the background of a visualization
* pan the dashboard by dragging the background of the Dashboard
## Zooming and Data in a Visualization
* zoom into a visualization by selecting the desired subregion with the mouse. This work on the marginals as well.
* zoom out to default by double clicking on a visualization
* hovering over marks to show their variable values
