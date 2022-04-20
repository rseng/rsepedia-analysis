## Installation and Usage Guidelines ##

Follow these guidelines to set up ADaPT-ML on your machine and to see how you can add new classification tasks to the system. Each header links to the appropriate file created for the [Example Use Case](./README.md#example-usage) so you can see an example of these instructions implemented.

### Step 1: Review System Requirements ###

#### Linux and MacOS Setup ####
1. Download and Install Docker and Docker Compose:
   - [Docker Engine v19.03.0+](https://docs.docker.com/)
   - [Docker Compose v1 1.29.2](https://docs.docker.com/compose/) (This software has not been tested with a newer version of Docker Compose)
3. Ensure CrateDB will pass the bootstrap checks by following [these instructions](https://crate.io/docs/crate/howtos/en/latest/admin/bootstrap-checks.html#linux), as the host system must be configured correctly to use CrateDB with Docker. 

#### Windows Setup ####
1. Download and Install Docker Desktop (Docker Compose is included):
   - [Docker Desktop](https://docs.docker.com/desktop/windows/install/)
   - If prompted to do so, download and install the [Linux kernel update package](https://docs.microsoft.com/en-us/windows/wsl/install-manual#step-4---download-the-linux-kernel-update-package). Complete steps 4-6 in the linked article.
3. Ensure CrateDB will pass the bootstrap checks by following [these instructions](https://stackoverflow.com/questions/69214301/using-docker-desktop-for-windows-how-can-sysctl-parameters-be-configured-to-per), copied below:
   1. In your Windows `%userprofile%` directory (usually `C:\Users\<username>`), create or edit the file `.wslconfig` with the following:
   ```
   [wsl2]
   kernelCommandLine = "sysctl.vm.max_map_count=262144"
   ```
   2. Exit any WSL instance through Command Prompt, `wsl --shutdown`, and restart your computer.

#### **It is recommended that at this point, you test ADaPT-ML by following these [instructions](./test/README.md). Additionally, if you came straight here without following the [Example Use Case](./README.md#example-usage), please consider doing so now so that you can see the following steps implemented.** ####

### Step 2: [Set up the environment variables for Docker Compose](./.env) ###
Make a copy of the `.env` file in the repository's root directory and call it `.env.dev`. Review the `.env.dev` file, and edit the variables according to their descriptions.

### Step 3: Changes to [label-studio](./label-studio) ###

Most of the setup for Label Studio is done through the UI that is launched at http://localhost:8080 by default, but there are a few things within this project directory to take note of, especially if you plan on using Label Studio's API.

#### (a) [Format your Labeling Config file](./label-studio/config/example_config.xml) ####

This configures how each component of a datapoint will be displayed to the annotators. This file can be copied and pasted into the Label Studio Labeling Configuration UI, or set for a certain project [using the API](https://labelstud.io/api#operation/api_projects_create).

#### (b) [Define your classification task name and classes](./label-studio/ls/__init__.py) ####

Until there is one configuration file for defining the classification task name and classes across all steps in the ADaPT-ML pipeline (see [Contributing](#community-guidelines)), you will need to update the `CLASSIFICATION_TASKS` variable with your new task name and corresponding classes.

#### Please note: ####

The Example Use Case demonstrates how to add a new classification task with only a text component for each datapoint. Therefore, it may be necessary to make changes to the [task sampling](./label-studio/ls/sample_tasks.py), [annotation processing](./label-studio/ls/process_annotations.py), and/or [annotator agreement](./label-studio/ls/annotator_agreement.py) modules if Label Studio's JSON import and export format is different according to the datapoint's number of components (e.g. both text and image), number of annotators, etc. See [Contributing](#community-guidelines).

### Step 4: Changes to [data-programming](./data-programming) ###

Setting up the data-programming project within ADaPT-ML to work with a new classification task requires adding new Python modules and editing some existing files.

#### (a) [Define your class names](./data-programming/label/lfs/__init__.py) ####

Until there is one configuration file for defining the classification task name and classes across all steps in the ADaPT-ML pipeline (see [Contributing](#community-guidelines)), this is where you need to define the Class that will hold both the name of each class and the number representing that class, which the Labeling Functions will use to vote, and which will ultimately make up the Label Matrix. **NOTE**: if your task is specifically a **binary** task, then you need to use the suffix `_pos` for the positive class (and optionally `_neg` for the negative class) in order to have the correct binary classification metrics downstream.

#### (b) [Write your Labeling Functions](./data-programming/label/lfs/example.py) ####

Create a module within [./data-programming/label/lfs](data-programming/label/lfs). This module that you can name after your new classification task is where you will write your Labeling Functions, and create a function called `get_lfs` that will produce an iterable containing all of the Labeling Functions you have defined.

#### (c) [Create your main function as an MLflow endpoint](./data-programming/label/example.py) ####

Create a module within [./data-programming/label](./data-programming/label). This is the main module for your new task. You will need to import the Class you defined in [Step 4(a)](#a-define-your-class-namesdata-programminglabellfs__init__py) and the `get_lfs` function defined in [Step 4(b)](#b-write-your-labeling-functionsdata-programminglabellfsexamplepy). You will also need to create a name for the Label Model that will be specific to your new task, and a dictionary with the names of the columns holding the features extracted for use with the Labeling Functions you defined as keys and any functions necessary to properly transform or unpack the featurized data point as values. You will also need to specify the path within the Label Studio annotations directory to the DataFrame that holds the annotated development data. Here you can add additional arguments to the argument parser if your Labeling Functions need them, like thresholds.

#### (d) [Add your MLflow endpoint to the MLproject file](./data-programming/MLproject) ####

This file is where you will specify the default hyperparameters for training the Label Model, additional parameters for your Labeling Functions, the type of classification your new task falls under (multiclass or multilabel), and the path to the main module you created in [Step 4(c)](#c-create-your-main-function-as-an-mlflow-endpointdata-programminglabelexamplepy). If you perform hyperparameter tuning and find a configuration that works well for your task, then change the defaults here!

### Step 5: Changes to [modelling](./modelling) (including model deployment) ###

There is not much that you have to edit in this project directory unless you need a machine learning algorithm other than a multi-layer perceptron (MLP), but if you do add a new algorithm, please see [Contributing](#community-guidelines)! For now, all of the edits are to the FastAPI app.

#### (a) [Add your class response format and endpoint to the deployment app](./modelling/app/main.py) ####

Until there is one configuration file for defining the classification task name and classes across all steps in the ADaPT-ML pipeline (see [Contributing](#community-guidelines)), you will need to add a response model that validates the output from your prediction endpoint. You will also need to create and set environment variables in [Step 2](#step-2-set-up-the-environment-variables-for-docker-composeenv) for your new End Model and add functions to load them. You can add an element to the `loaded_models_dict` for your model, so you will know if it loaded successfully by visiting the root page. Finally, you will need to add an endpoint to get predictions for new datapoints from your model. This endpoint can return a JSON response in the format of your specified response model, or directly update the data in CrateDB with the predictions. 

### Step 6: Start ADaPT-ML ###

Once you have your new classification task ready to go by completing Steps 1-5, all you need to do is:
```shell
cd ADaPT-ML/
docker-compose --env-file .env.dev --profile dev up -d
docker-compose ps
```
Once you see Docker Compose report this:
```
         Name                        Command                  State                                  Ports                            
--------------------------------------------------------------------------------------------------------------------------------------
crate-db                  /docker-entrypoint.sh crat ...   Up             0.0.0.0:4200->4200/tcp,:::4200->4200/tcp, 4300/tcp, 5432/tcp
dp-mlflow                 /bin/bash                        Up                                                                         
dp-mlflow-db              /entrypoint.sh mysqld            Up (healthy)   3306/tcp, 33060/tcp, 33061/tcp                              
dp-mlflow-server          mlflow server --backend-st ...   Up             0.0.0.0:5000->5000/tcp,:::5000->5000/tcp                    
label-studio-dev          /bin/bash                        Up                                                                         
label-studio-web          ./deploy/docker-entrypoint ...   Up             0.0.0.0:8080->8080/tcp,:::8080->8080/tcp                    
modelling-mlflow          /bin/bash                        Up                                                                         
modelling-mlflow-db       /entrypoint.sh mysqld            Up (healthy)   3306/tcp, 33060/tcp, 33061/tcp                              
modelling-mlflow-deploy   /start.sh                        Up             0.0.0.0:80->80/tcp,:::80->80/tcp                            
modelling-mlflow-server   mlflow server --backend-st ...   Up             0.0.0.0:5001->5000/tcp,:::5001->5000/tcp 
```
Then it's ready! Import your data into a table in CrateDB and refer to the [Example Usage](./README.md) and this [script](./example_data/example_data_import.py) for an example of how to manipulate the data so that it's ready for ADaPT-ML. How you load the data, featurize it, and sample from it to create your unlabeled training data is up to you -- ADaPT-ML does not perform these tasks. However, there may be an opportunity for certain sampling methods to become a part of the system; see [Contributing](#community-guidelines).

### Optional: Create a dev/test dataset using Label Studio ###

If you have two or more domain experts available to label some datapoints in order to create a gold dev/test dataset, then you can follow these steps to use Label Studio to accomplish this.

#### (a) Sample some data from CrateDB ####

Use this module to sample N random datapoints from a table in CrateDB, making sure to include the columns that contain the data that the domain exports will use during annotation.
```shell
docker exec label-studio-dev python ./ls/sample_tasks.py --help
```
```
usage: sample_tasks.py [-h] [--filename FILENAME] table columns [columns ...] n {example}

Sample a number of data points from a table to annotate.

positional arguments:
  table                Table name that stores the data points.
  columns              column name(s) of the data point fields to use for annotation.
  n                    Number of data points to sample.
  {example}            What classification task is this sample for?

optional arguments:
  -h, --help           show this help message and exit
  --filename FILENAME  What would you like the task file to be called?
```

#### (b) Use Label Studio's UI or API to label the sampled datapoints ####

Please refer to these guides to create an account ([UI](https://labelstud.io/guide/signup.html#Create-an-account) or [API](https://labelstud.io/api#operation/api_users_create)), create a project ([UI](https://labelstud.io/guide/setup_project.html#Create-a-project) or [API](https://labelstud.io/api#operation/api_projects_create)), load the tasks file created in [Optional (a)](#optional-create-a-devtest-dataset-using-label-studio) ([UI](https://labelstud.io/guide/tasks.html#Import-data-from-the-Label-Studio-UI) or [API](https://labelstud.io/guide/tasks.html#Import-data-using-the-API)), label the tasks ([UI](https://labelstud.io/guide/labeling.html) or [API](https://labelstud.io/api#tag/Tasks)), and export the resulting annotations ([UI](https://labelstud.io/guide/export.html#Export-using-the-UI-in-Label-Studio) or [API](https://labelstud.io/guide/export.html#Export-using-the-API)).

#### (c) Process the exported annotations ####

Once you have exported the annotations and moved the file to `${LS_ANNOTATIONS_PATH}`, the annotations need to be processed using this module:
```shell
docker exec label-studio-dev python ./ls/process_annotations.py --help
```
```
usage: process_annotations.py [-h] filename {example} gold_choice

Format exported annotations into DataFrames ready for downstream functions.

positional arguments:
  filename     Name of the exported annotations file.
  {example}    Which task is the annotations file for?
  gold_choice  How to settle disagreements between workers. id: Provide the id of the worker whose labels will be chosen every time. random: The least strict. Choose the label that the majority of workers agree on. If they are evenly split, choose a worker label randomly. majority: More strict. Choose the
               label that the majority of workers agree on. If they are evenly split, drop that datapoint. drop: The most strict. If workers disagree at all, drop that datapoint.

optional arguments:
  -h, --help   show this help message and exit
```
This will save three DataFrames in `${LS_ANNOTATIONS_PATH}/[CLASSIFICATON_TASK]`, where `CLASSIFICATION_TASK` is the name of your new classification task defined in [Step 3(b)](#b-define-your-classification-task-name-and-classeslabel-studiols__init__py):
- `ann_df.pkl` contains all of the datapoints that were initially exported from Label Studio with a column for each annotator's label set.
- `task_df.pkl` contains only the datapoints that were labeled by all annotators working on the project (e.g., if worker 1 labeled 50 datapoints and worker 2 labeled 45, then this DataFrame will contain 45 datapoints.)
- `gold_df.pkl` contains the final gold label set that was compiled according to the method selected using the `gold_choice` argument.

#### (d) Calculate interannotator agreement ####

Before moving on with the gold labels in `gold_df.pkl`, this module should be used to determine the level of agreement between all of the annotators:
```shell
docker exec label-studio-dev ./ls/annotator_agreement.py --help
```
```
usage: annotator_agreement.py [-h] {example}

Compute the inter-annotator agreement for completed annotations.

positional arguments:
  {example}   Task to calculate agreement for.

optional arguments:
  -h, --help  show this help message and exit
```
This will log and print a report using Krippendorff's alpha.

### Step 7: Run a data programming experiment ###

Once you have determined how you will sample some of your data for training an End Model, you need to save it as a pickled Pandas DataFrame with columns `id` and `table_name`, and optionally other columns if you need them. `table_name` needs to have the name of the table in CrateDB where the datapoint is stored. Once this DataFrame is in the directory `$DP_DATA_PATH/unlabeled_data`, you can run this command to label your data:
```shell
docker exec dp-mlflow sh -c ". ~/.bashrc && wait-for-it dp-mlflow-db:3306 -s -- mlflow run --no-conda -e [ENTRYPOINT] --experiment-name [EXP_NAME] -P train_data=/unlabeled_data/[DATA] -P dev_data=[0,1] -P task=[TASK] ."
```
where `ENTRYPOINT` is the name of the entrypoint you specified in [Step 4(d)](#d-add-your-mlflow-endpoint-to-the-mlproject-filedata-programmingmlproject), `EXP_NAME` is a name for the experiment of your choosing, `DATA` is the name of the pickled Pandas DataFrame holding your unlabeled data, `[0, 1]` is the flag to set indicating that you have done [Optional (c)](#c-process-the-exported-annotations) to create the `${LS_ANNOTATIONS_PATH}/[CLASSIFICATION_TASK]/gold_df.pkl`for your classification task (`1`), or that you do not have a gold dataset available from Label Studio (`0`), and `TASK` is the type of classification that is appropriate for your new task (multiclass or multilabel). You can then check http://localhost:5000 to access the MLflow UI and see the experiment log, Labeling Function evaluation, artifacts, metrics, and more. Your labeled data will be stored in the directory `${DP_DATA_PATH}/mlruns/EXP_ID/RUN_ID/artifacts/training_data.pkl` where `EXP_ID` is the id corresponding to `EXP_NAME`, and `RUN_ID` is a unique id created by MLflow for the run.

### Step 8: Run a modelling experiment ###

Once you have run some experiments and are happy with the resulting labeled data, take note of the `EXP_ID` and `RUN_ID` from [Step 7](#step-7-run-a-data-programming-experiment) within the filepath to the `training_data.pkl` and `development_data.pkl` (or, if you don't have a gold dataset from Label Studio, then instead of `development_data.pkl` you will use the DataFrame you split off of `training_data.pkl` and saved in the `artifacts` folder). Then you can run this command to train and evaluate an MLP model:
```shell
docker exec modelling-mlflow sh -c ". ~/.bashrc && wait-for-it modelling-mlflow-db:3306 -s -- mlflow run --no-conda -e mlp --experiment-name [EXP_NAME] -P train_data=/dp_mlruns/[EXP_ID]/[RUN_ID]/artifacts/training_data.pkl -P test_data=/dp_mlruns/[EXP_ID]/[RUN_ID]/artifacts/[TEST_DATA] -P features=FEATURE ."
```
where `EXP_NAME` is a name for the experiment of your choosing, `EXP_ID` and `RUN_ID` are from your evaluation from [Step 7](#step-7-run-a-data-programming-experiment), `TEST_DATA` is either `development_data.pkl` or the name of the Pandas DataFrame holding your testing data split off of `training_data.pkl`, and `FEATURE` is a list of column names holding the feature vectors in CrateDB. You can then check http://localhost:5001 to access the MLflow UI and see the experiment log, artifacts, metrics, and more.

### Step 9: Deploying your model ###

After you are satisfied with the performance of an End Model created in [Step 8](#step-8-run-a-modelling-experiment), take note of the `EXP_ID` and `RUN_ID` for the End Model, and update your End Model's [environment variable](#step-2-set-up-the-environment-variables-for-docker-composeenv) to `/mlruns/[EXP_ID]/[RUN_ID]/artifacts/mlp/python_model.pkl`. Then, edit the `environment` section of the `m_deploy` service in [docker-compose.yml](./docker-compose.yml) so that it has your End Model's environment variable.

Now you can reload the deployment API by running these commands:
```shell
docker-compose stop
docker-compose --profile dev up -d
```
and visit http://localhost:80/docs to see the deployment API. You can use this API to get predictions on unseen datapoints, and take note of the `curl` command to get predictions. It should look something like this:
```shell
curl -X 'POST' \
  'http://localhost/[ENDPOINT]' \
  -H 'accept: application/json' \
  -H 'Content-Type: application/json' \
  -d '{
  "table_name": [
    [TABLE_NAME]
  ],
  "id": [
    [ID]
  ]
}'
```
where `ENDPOINT` is the one you created in [Step 5(a)](#step-5-changes-to-modellingmodelling-including-model-deployment), `TABLE_NAME` is a list of names of the table(s) containing the datapoints you need predictions for, and `ID` is the list of ids for the datapoints.


You now have predicted labels for your data and can perform any downstream analyses you need!

## Additional Installation Notes: ##

Check this out if you are hosting CrateDB or another SQLAlchemy-based database on a remote server:
https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_sql.html

If you want to train the Label Model using CUDA tensors, then please refer to these resources:
https://developer.nvidia.com/cuda-toolkit
https://pytorch.org/docs/stable/cuda.html

## Community Guidelines ##

### Contribution ###
Follow these guidelines to see where you can contribute to expand the system's functionality and adaptability. The following items are on ADaPT-ML's "wish list":
- a configuration file that can be used by the label-studio, data-programming, and modelling projects to automatically create the classification task directory for label studio, a coding schema for annotators, the Enum object that stores values that the Labeling Functions use, the ModelResponse schema for deployment, and anything else where it is important to have consistency and maintainability in the classification task name and classes.
- a main UI with links to all of the different UIs, buttons that can run commands to sample data and run end-to-end experiments by returning the `EXP_ID` and `RUN_ID` within mlruns for a successful and performant Label Model and End Model, forms for submitting new classification tasks, an interface that makes writing labeling functions easier, etc.
- implement some algorithms that can take a representative sample of a table in CrateDB for training data creation.
- implement classification algorithms in addition to the MLP.
- determine the best method for updating the CrateDB tables with worker labels, gold labels, Label Model labels and probabilities, and End Model predictions and probabilities.
- a separate project for creating a flexible feature store.

Please [open an issue](https://github.com/U-Alberta/ADaPT-ML/issues) if you would like to propose an approach to adding these features.

### Report Issues and Seek Support ###
If you find a problem with the software or if you need help with any of the steps in this document or the testing document, please [open an issue](https://github.com/U-Alberta/ADaPT-ML/issues) and I will try to address your concerns.
[![status](https://joss.theoj.org/papers/e846e08311cee3886d33101209166f4c/status.svg)](https://joss.theoj.org/papers/e846e08311cee3886d33101209166f4c)

# ADaPT-ML #
## A Data Programming Template for Machine Learning ##

Often when studying natural phenomena by creating data-driven models, processing the data becomes the largest challenge. Without a framework to build upon and implement one's ideas, researchers are forced to hastily build inflexible programs from the ground up. When hypotheses need to be reworked or modelling a new aspect of the phenomena becomes necessary, even more time is spent on the program before finally being able to test out new ideas. This inherently causes problems, with additional problems arising such as including internal and external validation steps as an afterthought rather than a checkstop in the pipeline.

ADaPT-ML aims to be the flexible framework upon which researchers can implement their understanding of the phenomena under study. This software was created especially for any researcher with:

* Some programming experience or interest in learning how to write code based off of examples. 

* Access to large amounts of unlabeled data that is constantly changing, such as social media data. 

* Domain expertise or an intuition about how they would follow rules, heuristics, or use knowledge bases to annotate the unlabeled data. 

ADaPT-ML takes as much of the development work as possible out of creating novel models of phenomenon for which we have well-developed theories that have yet to be applied to big data.

## Introduction ##

ADaPT-ML is composed of a number of open-source tools and libraries, as shown in this system diagram. To familiarize yourself with these components, please review the tools and libraries linked to below the diagram.

![System Diagram](./graphics/system.png)

- [Docker](https://docs.docker.com/) and [Docker Compose](https://docs.docker.com/compose/)
  - working with environment variables, volumes, Dockerfiles, and the main docker-compose commands
- [Pandas](https://pandas.pydata.org/)
  - DataFrames
- [Label Studio v1.0](https://labelstud.io/)
  - Labeling Config files
- [Snorkel v0.9.7](https://www.snorkel.org/get-started/)
  - writing Labeling Functions, Label Matrix, Label Model
- [MLflow v1.19.0](https://mlflow.org/)
  - MLflow Projects, MLflow Tracking, MLflow Models, Model Registry
- [FastAPI v0.68.1](https://fastapi.tiangolo.com/)
  - endpoints, Pydantic, requests, JSON

Now that you are familiar with the concepts, terminology, and tools that make up ADaPT-ML, let's look at the example use case included in this repository. Once you have an understanding of how ADaPT-ML works and want to get started with your own use case, please refer to these instructions for [testing](test/README.md) ADaPT-ML on your machine, and the [usage guidelines](./usage.md), including how to contribute to this project.

## Example Usage ##

Our Example Use Case is to develop a model that can predict whether a data point is about a cat, dog, bird, horse, or snake. Although intuitively this is purely a multilabel task where it is reasonable to assume that one or more animals could be mentioned in one datapoint, this task has been divided into a multiclass setting, where there is only one possible class that the data point can belong to, and a multilabel setting, where one data point can belong to one or many classes, to demonstrate how to handle both tasks (it is not necessary for you to also divide your new classification task into multiclass and multilabel settings). 

All of the directories and files mentioned in the following steps exist in the locations specified in the `.env` file of this repository. To follow along using the various UIs, complete [Step 1](usage.md#step-1-review-system-requirements) and these [tests](test/README.md) to get ADaPT-ML running on your host machine, and go to the following addresses in your web browser of choice:
1. `localhost:4200` for CrateDB
2. `localhost:8080` for Label Studio
3. `localhost:5000` for data programming MLflow
4. `localhost:5001` for modelling MLflow
5. `localhost:81/docs` for FastAPI

### Step 1: obtain and featurize some data ###

We do not have an existing annotated dataset for this classification task, so the first step will be to create one. When you first get started, you will need to gather the appropriate data for your task, and featurize it in two ways:
1. Decide on which features you would pull out to assist in annotation if you were going to manually assign classes to the datapoints.
2. Decide on how you would represent the datapoints as feature vectors for the End Model.
Again, to keep it simple for this use case, our first feature set is simply lemmatized tokens. Our second feature set is the output from the Universal Sentence Encoder, given the raw text as input.

In this use case, data points were manually created with only a text component to keep it simple, but consider the tweets **1a**-**1e** in the diagram below. 

![Step One](./graphics/step_1.png)

Many of them have both text and images that can provide information for more accurate classification. Let's run through each datapoint:
- **1a** has a textual component with the keyword "Birds", and an image component that is a painting where birds can be identified.
- **1b** has a textual component with the keyword "horse", and an image component that does not on its own provide information that it is related to horses.
- **1c** only has a textual component with the keyword "dogs".
- **1d** has a reference to cats in its textual component with a :pouting_cat: emoji and hashtags. Its image component is an obvious picture of a cat.
- **1e** has a textual component that has the keyword "snake", but it is actually not about the animal. The image component does not on its own provide information that is related to snakes.

This diagram demonstrates the process of setting up the example use case data in a table in CrateDB so that it is ready for ADaPT-ML. You can refer to [this script](example_data/example_data_import.py) to see how this was accomplished in detail. As long as each table has these essential columns, you can combine multiple tables to create your training and testing data:
- column **1f** is the _essential_ `id` column: your data table must have this column to work with ADaPT-ML. It can be any combination of numbers and letters.
- column **1g** has the unprocessed, raw text component of each datapoint. Note the column name `txt`: this is the column used in Label Studio for annotation, so this name appears in the [Labeling Config file](./label-studio/config/example_config.xml): in `<Text name="txt" value="$txt"/>` and `<Choices name="topic" toName="txt" choice="multiple" showInLine="false">`
- column **1h** has the first featurization step -- features that will be used by the Labeling Functions. Note the column name `txt_clean_lemma`: this column name is specified in our [data programming MLflow endpoint](./data-programming/label/example.py) for this example task, with no loader function: `LF_FEATURES = {'txt_clean_lemma': None}`. For this use case, as mentioned before, our Labeling Functions will have access to the text in _1g_ that has been lemmatized. Consider, though, a more rich set of features for the text component such as word embeddings, emoji normalization, hashtag splitting, and so on. If our datapoint had an image component, then we could pull out features such as the prediction from an image classifier that has been trained on one or more categories included in our task, the output from an image classifier that can detect language embedded within the image, and so on. The output from all of these different featurization methods would be held in different columns in this table.
- column **1i** has the second featurization step -- features that will be used to train and get predictions from the End Model. For this use case, we have put the text in _1g_ through the Universal Sentence Encoder to obtain a semantic representation of the text. Unlike the first featurization step, these feature columns can _only be real-number vectors_. You can create multiple feature vectors for the datapoint's different components, and ADaPT-ML will concatenate these vectors in the order given, not necessarily the order it appears in the table. 

### Step 2: create a gold dataset using Label Studio ###

We are now ready to annotate a sample of data in Label Studio! Because we only have a total of 15 datapoints for the multiclass setting and 15 for the multilabel setting, they were all annotated manually, but in a real-world application of this classification task, it is likely we would have hundreds of thousands of datapoints. In this case, we would instruct two or more annotators to manually label a few hundred datapoints for a few purposes:
1. Gather feedback from the annotators to inform how we can update or create new Labeling Functions
2. Estimate the class balance and make it available to the Label Model during training
3. Perform an empirical evaluation of the Labeling Functions and Label Model
4. Validate the End Model

The first step to annotate data using Label Studio is to set up the project using the Label Studio UI. For this example use case, we enter `localhost:8080` (if you changed the port in `docker-compose.yml`, replace `8080` with what you entered) in a web browser. Create an account, and set up the project (we simply called it "example").

The second step is to sample some data from CrateDB. The sampling method implemented currently in ADaPT-ML is a random N, so this commannd was used to sample all 30 datapoints for the multiclass and multilabel settings:
```shell
docker exec label-studio-dev python ./ls/sample_tasks.py example_data txt 30 example --filename example_tasks.json
```
This module will format the data in the column names provided so that it can be read by Label Studio, and save a file in the `$LS_TASKS_PATH` directory. The diagram below shows the process of using Label Studio to import the sampled data, annotate it, and export it.

![Step Two](./graphics/step_2.png)

- **2a** shows all the created projects, allowing for multiple classification problems. 
- **2b** has the popup for importing our `example_tasks.json` file.
- **2c** shows the imported tasks, which can be clicked on to arrive at the labeling page.
- **2d** is the main labeling interface determined by the [example config file](label-studio/config/example_config.xml), with different tabs for multiple annotators.
- **2e** once annotation is complete, they are exported in JSON format to the host machine's Downloads folder by default. For this demonstration, they were manually moved and renamed to `$LS_ANNOTATIONS_PATH/annotations/example_annotations.json`

Now that we have labeled all of our sample data and exported the results, we need to process the JSON file back into the Pandas DataFrames that ADaPT-ML can use. Because we had multiple annotators label each datapoint, we need to decide how we want to compile these labels into one gold label set. These two tasks are accomplished through this command:
```shell
docker exec label-studio-dev python ./ls/process_annotations.py example_annotations.json example 1
```
The following DataFrames are saved in `$LS_ANNOTATIONS_PATH/example`:
- `ann_df.pkl` contains all of the datapoints that were initially exported from Label Studio with a column for each annotator's label set.
- `task_df.pkl` contains only the datapoints that were labeled by all annotators working on the project (e.g., if worker 1 labeled 50 datapoints and worker 2 labeled 45, then this DataFrame will contain 45 datapoints.)
- `gold_df.pkl` contains the final gold label set that was compiled according to the method selected using the `gold_choice` argument. In this case, worker_1's labels were selected.
Finally, before we use the gold labels, we need to check the level of agreement between the annotators. To do this, we use this command:

```shell
docker exec label-studio-dev python ./ls/annotator_agreement.py example
```
This module uses `task_df.pkl` to calculate Krippendorff's alpha. For demonstration purposes, worker_2 intentionally disagreed with worker_1 on several datapoints, where worker_1 made all correct choices. Between worker_1 and worker_2, the agreement report looks like this:
```
TASK: example
NOMINAL ALPHA: 0.43847133757961776
RESULT: 0.43847133757961776 < 0.667. Discard these annotations and start again. 
```
This would normally prompt an iteration on the labelling process, but we are choosing only worker_1's labels for the gold dataset.

### Step 3: use data programming to create a labeled dataset ###

Now that we have our gold labels, we are ready to perform data programming to label more data. We have followed [these instructions](./usage.md) to modify ADaPT-ML for our example classification task. We also sampled some data from CrateDB that we want to use as training data; for this example use case, we have one DataFrame with the multiclass datapoints and one DataFrame with the multilabel datapoints, and both DataFrames only have the columns `id` and `table_name`. The DataFrames are called `multiclass_df.pkl` and `multilabel_df.pkl`, and both are stored in `$DP_DATA_PATH/unlabeled_data`.

Once we run the commands in the following code block...
```shell
docker exec dp-mlflow sh -c ". ~/.bashrc && wait-for-it dp-mlflow-db:3306 -s -- mlflow run --no-conda -e example --experiment-name eg -P train_data=/unlabeled_data/multiclass_df.pkl -P dev_data=1 -P task=multiclass -P seed=8 ."

docker exec dp-mlflow sh -c ". ~/.bashrc && wait-for-it dp-mlflow-db:3306 -s -- mlflow run --no-conda -e example --experiment-name eg -P train_data=/unlabeled_data/multilabel_df.pkl -P dev_data=1 -P task=multilabel -P seed=8 ."
```
...we can check out the results using the MLflow UI, as seen in the diagram below. 

![Step 3](./graphics/step_3.png)

- **3a** shows the landing page for MLflow, where the Experiment ID (`EXP_ID`) is shown, and runs can be filtered and selected.
- **3b** has some metadata for the run, including the run command showing that default parameters were used. The Run ID (`RUN_ID`) is shown at the top of the page.
- **3c** shows the metrics available under the multilabel setting.
- **3d** the result of Labeling Function evaluation against the gold labels in the `gold_df` are seen in `lf_summary_dev.html`.
- **3e** `development_data.pkl` has the aforementioned columns in addition to the `gold_label` column from the `gold_df`.
- **3f** `training_data.pkl` has the Label Model's predicted labels in the column `label` and the probability distribution over all classes in the column `label_probs`.

Once we have experimented with the Label Model parameters, Labeling Functions, and datasets to our satisfaction, we can make note of the experiment ID (`EXP_ID`) and run ID (`RUN_ID`) to access the `training_data.pkl` and `development_data.pkl` that we want to use in End Model training and evaluation. For ease of demonstration, these artifacts have been placed in `${DP_DATA_PATH}/mlruns`, but normally these artifacts would be found in `${DP_DATA_PATH}/mlruns/EXP_ID/RUN_ID/artifacts`.

### Step 4: create an End Model ###

Now that we have our training data labeled by the Label Model and testing data with gold labels, we can create an End Model that, given a DataFrame containing only the `id` and `table_name` columns, will look up the appropriate features for each datapoint in CrateDB and produce a DataFrame with a binary encoding for the predicted class(es), and the probability distribution over all classes. Currently, ADaPT-ML only has one machine learning algorithm, [scikit-learn's Multi-layer Perceptron (MLP)](https://scikit-learn.org/stable/modules/generated/sklearn.neural_network.MLPClassifier.html), a classifier that optimizes the log-loss function using LBFGS or stochastic gradient descent.

After running the commands in this code block...
```shell
docker exec modelling-mlflow sh -c ". ~/.bashrc && wait-for-it modelling-mlflow-db:3306 -s -- mlflow run --no-conda -e mlp --experiment-name eg -P train_data=/dp_mlruns/multiclass_training_data.pkl -P test_data=/dp_mlruns/multiclass_development_data.pkl -P features=txt_use -P solver=lbfgs -P random_state=8 ."

docker exec modelling-mlflow sh -c ". ~/.bashrc && wait-for-it modelling-mlflow-db:3306 -s -- mlflow run --no-conda -e mlp --experiment-name eg -P train_data=/dp_mlruns/multilabel_training_data.pkl -P test_data=/dp_mlruns/multilabel_development_data.pkl -P features=txt_use -P solver=lbfgs -P random_state=8 ."
```
...we can check out the results in MLflow, as shown in the diagram below.

![Step 4](./graphics/step_4.png)

- **4a** shows all the training parameters for the MLP model, in addition to the specific features used for the given training data. Note that the `txt_use` parameter is the column 1i in the [diagram for Step 1](#step-1-obtain-and-featurize-some-data).
- **4b** if the test data has gold labels in addition to the Label Model's labels, then the model will be evaluated against both, and they will show up here with "GOLD" and "LM" prefixes. If there are no gold labels, then the metrics will only have "LM".
- **4c** is the confusion matrix for the gold labels and predicted labels in the multiclass setting.
- **4d** is the confusion matrix for the gold labels and predicted labels in the multilabel setting.
- **4e** is the DataFrame with the model's predicted label(s) and the probability distribution over all classes.

Once we have experimented with the MLP parameters, and possibly iterated more on the data programming step if necessary, we can prepare our models for deployment by simply updating the model environment variables in `.env` and the `environment` section of the `m_deploy` service in [docker-compose.yml](./docker-compose.yml) to point to `python_model.pkl`. For this example use case, multiclass and multilabel models were copied and renamed to `${MODELLING_DATA_PATH}/mlruns/multiclass_model.pkl` and `${MODELLING_DATA_PATH}/mlruns/multilabel_model.pkl`.

### Step 5: deploy the End Model ###

This diagram shows the FastAPI UI for the deployed models.

![Step 5](./graphics/step_5.png)

- **5a** on the landing page, we can see the endpoints for predicting multiclass and multilabel.
- **5b** below the endpoints is a schema of the expected format for posting data, and the data format that the model will return in response.
- **5c** this shows an example of how to use the UI to get the `curl` command for the endpoint, and what the response looks like.

Now we can get multiclass predictions...
```shell
curl -X 'POST' \
  'http://localhost/predict_multiclass_example' \
  -H 'accept: application/json' \
  -H 'Content-Type: application/json' \
  -d '{
  "table_name": [
    "example_data"
  ],
  "id": [
    "20"
  ]
}'
```
```
{
  "table_name": [
    "example_data"
  ],
  "id": [
    "20"
  ],
  "cat": [
    0
  ],
  "dog": [
    0
  ],
  "bird": [
    0
  ],
  "horse": [
    0
  ],
  "snake": [
    1
  ],
  "prob_cat": [
    5.6850715594352195e-8
  ],
  "prob_dog": [
    0.0001963686969921083
  ],
  "prob_bird": [
    8.922841061481865e-8
  ],
  "prob_horse": [
    8.82467128837139e-9
  ],
  "prob_snake": [
    0.9998034763992105
  ]
}
```
...and multilabel predictions...
```shell
curl -X 'POST' \
  'http://localhost/predict_multilabel_example' \
  -H 'accept: application/json' \
  -H 'Content-Type: application/json' \
  -d '{
  "table_name": [
    "example_data"
  ],
  "id": [
    "03"
  ]
}'
```
```
{
  "table_name": [
    "example_data"
  ],
  "id": [
    "03"
  ],
  "cat": [
    1
  ],
  "dog": [
    1
  ],
  "bird": [
    0
  ],
  "horse": [
    0
  ],
  "snake": [
    0
  ],
  "prob_cat": [
    0.999976879069893
  ],
  "prob_dog": [
    0.9999725147168369
  ],
  "prob_bird": [
    2.061596293323691e-8
  ],
  "prob_horse": [
    1.7205732529738035e-7
  ],
  "prob_snake": [
    2.0265644234853424e-8
  ]
}
```
...for any datapoint that has the `txt_use` feature set in CrateDB. We have successfully created a model for our new example use case!
---
title: 'ADaPT-ML: A Data Programming Template for Machine Learning'
tags:
  - Python
  - Docker
  - MLOps
  - data programming
  - machine learning
  - model lifecycle
authors:
  - name: Andrea M. Whittaker^[first author] # note this makes a footnote saying 'first author'
    orcid: 0000-0001-8512-1810
    affiliation: 1
affiliations:
  - name: University of Alberta
    index: 1
date: 21 October 2021
bibliography: paper.bib
---

# Summary

Classification is a task that involves making a prediction about which class(es) a data point 
belongs to; this data point can be text, an image, audio, or can even be multimodal. This task can become intractable for many reasons, including:

* Insufficient training data to create a data-driven model; available training data may not be appropriate for the domain being studied, it may not be of the right type (e.g. only text but you want text and images), it may not have all of the categories you need, etc.

* Lack of available annotators with domain expertise, and/or resources such as time and money to label large amounts of data.

* Studying a phenomenon that changes rapidly, so what constitutes a class may change over time, making the available training data obsolete.

ADaPT-ML (\autoref{fig:system}) is a multimodal-ready MLOps system that covers the data processing, data labelling, model design, model training and optimization, and endpoint deployment, with the particular ability to adapt to classification tasks that have the aforementioned challenges. ADaPT-ML is designed to accomplish this by:

* Using `Snorkel` [@Snorkel] as the data programming framework to create large, annotated, multimodal datasets that can easily adapt to changing classification needs for training data-driven models.

* Integrating `Label Studio` [@LabelStudio] for annotating multimodal data.

* Orchestrating the Labelling Function / Label Model / End Model development, testing, and monitoring using `MLflow` [@MLflow].

* Deploying all End Models using `FastAPI` [@FastAPI]

![Overview of how ADaPT-ML orchestrates data labelling and the end model lifecycle.\label{fig:system}](graphics/system.png){ width=70% }

# Statement of Need

Often when studying natural phenomena by creating data-driven models, processing the data becomes the largest challenge. Without a framework to build upon and implement one's ideas, researchers are forced to hastily build inflexible programs from the ground up. When hypotheses need to be reworked or modelling a new aspect of the phenomena becomes necessary, even more time is spent on the program before finally being able to test out new ideas. This inherently causes problems, with additional problems arising such as including internal and external validation steps as an afterthought rather than a checkstop in the pipeline.

ADaPT-ML aims to be the flexible framework upon which researchers can implement their understanding of the phenomena under study. This software was created especially for any researcher with:

* Some programming experience or interest in learning how to write code based off of examples. 

* Access to large amounts of unlabelled data that is constantly changing, such as social media data. 

* Domain expertise or an intuition about how they would follow rules, heuristics, or use knowledge bases to annotate the unlabelled data. 

ADaPT-ML takes as much of the development work as possible out of creating novel models of phenomenon for which we have well-developed theories that have yet to be applied to big data.

# Related Work
An early version of this software supported the modelling of universal personal values and was complementary to the software architecture described in @Gutierrez2021. 

During the development of this software, `Snorkel` progressed into `Snorkel Flow` [@SnorkelFlow], a proprietary MLOps system that incorporates data cleaning, model training and deployment, and model evaluation and monitoring into its existing data programming framework. Unlike ADaPT-ML's focus on adding new classification tasks by updating existing and creating new Python modules, `Snorkel Flow` allows users to perform complex operations with a push-button UI, like creating Labelling Functions from ready-made builders.

`Neu.ro` [@Neuro] and `MLRun` [@MLRun] are MLOps platforms that, like ADaPT-ML, orchestrate open-source tools like `Label Studio` and `MLflow`. However, supporting dynamic and iterative training data creation through data programming is a core feature of ADaPT-ML, but frameworks like `Snorkel` are not integrated out-of-the-box with `Neu.ro` and `MLRun`. Additionally, collaboration with the `Neu.ro` team is necessary to create sandboxes with the needed tools, whereas ADaPT-ML is completely open-source and hands-on.

# Acknowledgements
We acknowledge contributions from Mitacs and The Canadian Energy and Climate Nexus / Le Lien Canadien de L’Energie et du Climat for funding during the early stages of this project's development.

# References
# Testing #

The purpose of this document is to provide instructions for testing ADaPT-ML on your machine.
Please **make sure you have not modified the [`.env`](../.env) file**. These tests will bring all the containers up using docker compose, test core functionality using the data for the example use case included in [`example_data`](../example_data), and will **leave all the containers running** so that you can view the test experiment runs through the MLflow UIs and/or investigate any issues using the containers' respective log files.

## Unix: Run all tests ##

If you are using MacOS / Linux on your host machine, then run these commands:

```shell
cd ADaPT-ML/
bash ./test/all-unix.sh
```

## Windows: Run all tests ##

1. Open the Command Prompt in the ADaPT-ML root directory, and bring all the containers up:
```shell
docker-compose --env-file .env --profile dev up -d
```
2. Make sure that none of the containers have exited or are restarting and wait until you see that the MLflow databases are "healthy", like this:
```shell
docker-compose ps
```
```
NAME                      COMMAND                  SERVICE             STATUS              PORTS
crate-db                  "/docker-entrypoint.…"   cratedb             running             0.0.0.0:4200->4200/tcp
dp-mlflow                 "/bin/bash"              dp                  running
dp-mlflow-db              "/entrypoint.sh mysq…"   dp_db               running (healthy)   33060-33061/tcp
dp-mlflow-server          "mlflow server --bac…"   dp_web              running             0.0.0.0:5000->5000/tcp
label-studio-dev          "/bin/bash"              ls                  running
label-studio-web          "./deploy/docker-ent…"   ls_web              running             0.0.0.0:8080->8080/tcp
modelling-mlflow          "/bin/bash"              m                   running
modelling-mlflow-db       "/entrypoint.sh mysq…"   m_db                running (healthy)   33060-33061/tcp
modelling-mlflow-deploy   "/start.sh"              m_deploy            running             0.0.0.0:8088->80/tcp
modelling-mlflow-server   "mlflow server --bac…"   m_web               running             0.0.0.0:5001->5000/tcp
```
3. Run the test batch file:
```shell
test\all-windows.bat
```

## Optional Clean-up ##

If you wish to bring the containers down after testing, then:

```shell
docker-compose down
```

## Troubleshooting ##

- If the `crate-db` container's status shows `Restarting`, please check its logs, `docker logs crate-db`. If you see something like this in the logs:
    ```
  [1] bootstrap checks failed
  [1]: max virtual memory areas vm.max_map_count [65530] is too low, increase to at least [262144] by adding `vm.max_map_count = 262144` to `/etc/sysctl.conf` or invoking `sysctl -w vm.max_map_count=262144`
    ```
  then check that you completed the appropriate [bootstrap check steps](../usage.md#step-1-review-system-requirements).
- If after the `docker-compose ... up -d` command, an error indicating that there is no space left occurs while pulling the images, then try these steps:
  1. Check the amount of reclaimable space for images: `docker system df`
  2. If reclaimable is greater than zero, clean up unused containers, networks, images, and the build cache: `docker system prune -a`
  3. Try the testing procedure again from the beginning
  4. If there still is not enough space, consider relocating the root directory ([Linux](https://linuxconfig.org/how-to-move-docker-s-default-var-lib-docker-to-another-directory-on-ubuntu-debian-linux), [Windows](https://dev.to/kimcuonthenet/move-docker-desktop-data-distro-out-of-system-drive-4cg2))

For all other concerns with ADaPT-ML's tests, especially if any fail, please create an [Issue](https://github.com/U-Alberta/ADaPT-ML/issues).
