<br/> <br/>

<p align="center">
  <a href="https://github.com/SPARC-FAIR-Codeathon/SPARClink">
    <img src="https://github.com/SPARC-FAIR-Codeathon/SPARClink/blob/main/docs/images/logo.svg" alt="SPARC link logo" height="150">
  </a>
  <br/>
  <h3 align="center">
    Visualizing the Impact of SPARC
  </h3>
</p>

<p align="center">
  <a href="https://github.com/SPARC-FAIR-Codeathon/SPARClink/blob/main/LICENSE" alt="GitHub license">
    <img src="https://img.shields.io/github/license/SPARC-FAIR-Codeathon/SPARClink" />
  </a>
  <a href="https://github.com/SPARC-FAIR-Codeathon/SPARClink/stargazers" alt="GitHub stars">
    <img src="https://img.shields.io/github/stars/SPARC-FAIR-Codeathon/SPARClink" />
  </a>
  <a href="https://github.com/SPARC-FAIR-Codeathon/SPARClink/network" alt="GitHub forks">
    <img src="https://img.shields.io/github/forks/SPARC-FAIR-Codeathon/SPARClink" />
  </a>
  <a href="https://github.com/SPARC-FAIR-Codeathon/SPARClink/issues" alt="GitHub issues">
    <img src="https://img.shields.io/github/issues/SPARC-FAIR-Codeathon/SPARClink" />
  </a>
  <a href="https://github.com/SPARC-FAIR-Codeathon/SPARClink/graphs/contributors">
    <img src="https://img.shields.io/github/contributors/SPARC-FAIR-Codeathon/SPARClink" alt="GitHub contributors">
  </a>
  <a href="#">
    <img src="https://img.shields.io/github/last-commit/SPARC-FAIR-Codeathon/SPARClink" alt="GitHub last commit">
  </a>
<!--   <a href="#">
    <img src="https://img.shields.io/tokei/lines/github/SPARC-FAIR-Codeathon/SPARClink" alt="Lines of code">
  </a> -->
  <a href="#">
    <img src="https://badgen.net/badge/Open%20Source%20%3F/Yes%21/blue?icon=github" alt="Open Source? Yes!">
  </a>
  <br/> 
</p>
 <br/> <br/>
<p align="center">
  <img src="https://github.com/SPARC-FAIR-Codeathon/SPARClink/blob/main/docs/images/2021-07-25%2013-47-30.gif"/>
</p>
<br/> <br/>
 
 
## Table of content
- [What is SPARClink?](#what-is-sparclink)
  - [NIH SPARC](#nih-sparc)
  - [FAIR Data](#fair-data)
  - [Defining Impact](#defining-impact)
  - [Origin Story](#origin-story)
  - [Goal](#goal)
- [How it works](#how-it-works)
- [Run the project](#run-the-project)
  - [Testing](#testing)
  - [Firebase Backend Implementation](#firebase-backend-implementation)
  - [Visualization Web App](#visualization-web-app)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)
- [Further Reading](#further-reading)

## What is SPARClink?
### NIH SPARC
The NIH Common Fund’s Stimulating Peripheral Activity to Relieve Conditions (SPARC) program aims to transform our understanding of nerve-organ interactions with the intent of advancing bioelectronic medicine towards treatments that change lives. [Learn more about SPARC](https://sparc.science/about)

### FAIR Data
By employing a [FAIR](https://www.nature.com/articles/sdata201618) (Findable, Accessible, Interoperable and Reusable) first approach SPARC datasets, protocols and publications generated via the SPARC program is intended to be able to be used by researchers globally with reproducible results. However, at the current moment, there is no real tangible way to show or visualize the usage of SPARC data in outside projects and publications. 

### Origin Story
The SPARClink project was first born as an idea at the 2021 NIH SPARC Codeathon ([more details here](https://sparc.science/help/2021-sparc-fair-codeathon)). The idea behind the topic was created as a method of visualizing citation data on datasets, protocols and publications to determine the degree of use of SPARC material outside of the official channels.

### Defining Impact
The word 'Impact' can have many different meanings depending on the context that it is viewed in. Within the SPARClink project, we consider impact to be the frequency of citations of SPARC funded resources. The SPARC program intends to advance medical understanding by providing datasets, maps and computational studies that follow FAIR principles and is used by researchers all around the world. The usage of SPARC resouces by platforms and programs ouside SPARC is what we view as the meaning of the term 'Impact'.

### Goal
The goal of SPARClink is to provide a system that will query all external publications using open source tools and platforms and create an interactable visualization that is helpful to any person (researcher or otherwise) to showcase the impact that SPARC has on the overall scientific research community. These impact measurements are meant to be used as a showcase of the concept of FAIR data and how good data generation practices and methods are useful in advancing the field of bioelectronic medicine.

However, datasets and protocols are not referenced similar to prior research in manuscripts. Dataset and protocol identifiers or urls are only mentioned in text or under supplementary materials, making this a difficult task to accomplish.

## How it works?
Metadata information on datasets and protocols are extracted from [Pennsieve](https://app.pennsieve.io/), SPARC Airtable database, and [Protocols.io](https://www.protocols.io/workspaces/sparc). This information is queried against the [NIH RePORTER](https://api.reporter.nih.gov/), [NCBI](https://www.ncbi.nlm.nih.gov/), and [Google Scholar](https://serpapi.com/google-scholar-api) to extract citations and create a well connected graph using [d3.js](https://d3js.org/). 

<p align="center">
  <!--<img src="https://user-images.githubusercontent.com/21206996/125478715-d5f83b6f-8a6d-4ef8-a845-952baa27d8da.png" />-->
  <img src="https://github.com/SPARC-FAIR-Codeathon/SPARClink/blob/main/docs/images/sparclink_block_diagram-01.png" width="500"/>
  <br/>
  <span> SPARClink workflow </span>
</p>

## Run the project
Clone or download the repository.
``` bash
git clone https://github.com/SPARC-FAIR-Codeathon/SPARClink.git
```

The development environment uses [Anaconda](https://www.anaconda.com/products/individual) to keep track of the python dependencies. Download Anaconda here: [Anaconda Individual Edition](https://www.anaconda.com/products/individual).

The following would create a new `conda` environment with the dependencies required to run the project.
``` bash
cd SPARClink
conda env create -f environment.yml --prefix ./env 
conda activate ./env
```
The application uses [python-dotenv](https://github.com/theskumar/python-dotenv) to load configuration information from a `.env` file. Create a `.env` file with the following information.
``` bash
PROTOCOLS_IO_KEY="<protocols.io api key>"
SERPAPI_KEY="<serpapi api key>"
```
A public API key for protocols.io can be obtained by signing up as [shown here](https://www.protocols.io/developers). SERP api key is not required at the moment. To integrate google scholar results, an API key can be obtained as [shown here](https://serpapi.com/).

### Testing
Unit tests to verify external APIs are written in Python unittest framework. The tests can be run as shown below:
``` bash
python -m unittest -v tests/test_NIH_NCBI.py
```

### Firebase Backend Implementation
Currently, the central database is implemented as a [Firebase](https://firebase.google.com/) real-time database. The database can be updated by running `FirebaseImplementation.py`. However, this requires a username and a password.

To use your own Firebase instance, setup a Firebase web app as [shown here](https://firebase.google.com/docs/web/setup), and update `firebaseConfig` in `FirebaseImplementation.py` with the new API keys. [Setup a new user](https://firebase.google.com/docs/auth/web/password-auth), and configure the [real-time database](https://firebase.google.com/docs/database/web/start). It is recommended to limit the database write permission to authenticated users. Run `FireabaseImplementation.py` and enter user's email/password when prompted.

<p align="center">
  <!--<img src="https://user-images.githubusercontent.com/21206996/125478715-d5f83b6f-8a6d-4ef8-a845-952baa27d8da.png" />-->
  <img src="https://github.com/SPARC-FAIR-Codeathon/SPARClink/blob/main/docs/images/backend_flow_chart-01.png" width="500"/>
  <br/>
  <span>Backend Flow Chart: Shows the methods implemented in the backend to gather citations of datasets, protocols, and SPARC publications.</span>
</p>

### ML Data Indexing Engine
We have setup a [Flask](https://flask.palletsprojects.com/en/2.0.x/) server on [pythonanywhere](https://www.pythonanywhere.com/) to handle all our machine learning operations. If you would like to setup a backend for your own fork, please setup a flask server on any hosting service of your choice and modify the approriate endpoints in the `flask_app.py` file. To learn more about the techniques we used, refer to the [Further Reading](https://github.com/SPARC-FAIR-Codeathon/SPARClink#further-reading) section. 

### Visualization Web App
The vizualizations created from the realtime database can be viewed directly from our [demo page](https://sparclink-f151d.web.app/sparclink) or by running the local version of our frontend. We use [Vue.js](https://vuejs.org/) and [Tailwind CSS](https://tailwindcss.com/) to render the demo webpage. The interactive force directed graph is created via [d3.js](https://d3js.org/) using data requested from our Firebase real-time database. Within the SPARClink demo page we use the HTML canvas element to render the visualization. In order to get your forked repo frontend to run locally, use the following commands:
```bash
cd frontend
npm install
npm run serve
```
You can now open your browser and visit the url [http://localhost:8080/sparclink](http://localhost:8080/sparclink) to view the webpage. 


`Note:` To use the smart word filter, please refer to the frontend available in the [`smart_filter`](https://github.com/SPARC-FAIR-Codeathon/SPARClink/tree/smart_filter/frontend) branch. This feature will lead to slower render times on the graph visualization so we have not included it in the main branch.

<p align="center">
  <img src="https://github.com/SPARC-FAIR-Codeathon/SPARClink/blob/main/docs/images/2021-07-25 14-44-36.gif" />
  <br/>
  <span> SPARClink smart filter </span>
</p>

<!--Keep track of the project [here](https://github.com/SPARC-FAIR-Codeathon/SPARClink/projects/1)-->

## Maintainers
* [Sanjay Soundarajan](https://github.com/megasanjay)
* [Monalisa Achalla](https://github.com/a-monalisa)
* [Jongchan Kim](https://github.com/Kim-Jongchan)
* [Ashutosh Singh](https://github.com/Ashutosh1712)
* [Sachira Kuruppu](https://github.com/rsachira-abi)

## Contributing
If you would like to suggest an idea to this project, please let us know in the [issues](https://github.com/SPARC-FAIR-Codeathon/SPARClink/issues) page and we will take a look at your suggestion. Please use the `enhacement` tag to label your suggestion. 

If you would like to add your own feature, feel free to fork the project and send a pull request our way. This is an open source project so we will welcome your contributiobs with open arms. 
Refer to our [Contributing Guildeines](./docs/CONTRIBUTING.md) and [Code of Conduct](./docs/CODE_OF_CONDUCT.md) for more information. Add a [GitHub Star](https://github.com/SPARC-FAIR-Codeathon/SPARClink) to support active development!

## License
SPARClink is an open source project and distributed under the  MIT License. See [LICENSE](./LICENSE) for more details.

## Further Reading
- [External APIs](./ExternalAPIs/README.md)
- [SPARC APIs](./SPARC/README.md)
- [Visualization and frontend](./frontend/README.md)
- [ML Data Indexing Engine](./MLDataIndexingEngine/README.md)
# Contributing

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before making a change. 

Please note we have a code of conduct, please follow it in all your interactions with the project.

## Pull Request Process

1. Ensure any install or build dependencies are removed before the end of the layer when doing a 
   build.
2. Update the README.md with details of changes to the interface, this includes new environment 
   variables, exposed ports, useful file locations and container parameters.
3. Increase the version numbers in any examples files and the README.md to the new version that this
   Pull Request would represent. The versioning scheme we use is [SemVer](http://semver.org/).
4. You may merge the Pull Request in once you have the sign-off of two other developers, or if you 
   do not have permission to do that, you may request the second reviewer to merge it for you.

## Code of Conduct

### Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

### Our Standards

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

### Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

### Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at [INSERT EMAIL ADDRESS]. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
{"mode":"full","isActive":false}
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
ssoundarajan@calmi2.org.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# sparc-portal-demo

This demo page is built using [VueJS](https://vuejs.org/) and [Tailwind CSS](https://tailwindcss.com/). It has been developed to mimic the SPARC portal found at [sparc.science](https://sparc.science). In order to get this project running please refer to the [Project setup](https://github.com/SPARC-FAIR-Codeathon/SPARClink/blob/main/frontend/README.md#project-setup) section. 

## Project setup
``` bash
npm install
```

### Compiles and hot-reloads for development
``` bash 
npm run serve
```

### Compiles and minifies for production
``` bash
npm run build
```

### Lints and fixes files
``` bash
npm run lint
```

### Customize configuration
See [Configuration Reference](https://cli.vuejs.org/config/).


## Requesting data from the real-time database
SPARClink uses a Firebase real-time database as its intermediary data storage server. This allows for the backend citation extraction system to reside in a separate server and allow the frontend to pull data asynchronously. If you prefer another data storage mechanism please add and modify the appropriate end point in the SPARClink component found in the `frontend/src/components/SparcLink/SparcLink.vue` file. The url endpoint for the GET request can be found in the `organizeData` function in the `methods` section. Our backend systems only allows authenticated users to write to the database so you will need to refer to the appropriate User ID when referencing the object.

## Adjusting the physics of the visualization
[d3.js](https://d3js.org/) provides extensive documentation on the physics that we implement on our visualizations. If you would like to modify how the graph behaves, the best place to start is at the force simulations in the `drawCanvas` function. d3.js also provides the option of rendering graphs in SVG elements but the limit of nodes that a simulation can accept, before there is a considerable performance drop, is low.
```javascript
this.simulation = d3
  .forceSimulation(nodes)
  .force(
    "link",
    d3.forceLink(links).id((d) => d.id)
  )
  .force("x", d3.forceX())
  .force("y", d3.forceY())
  .force("charge", d3.forceManyBody().strength(this.strength))
  .force("center", d3.forceCenter(WIDTH / 2, HEIGHT / 2));
```
You may define your own forces to act on the nodes or use d3's default forces. Please be aware that the canvas itself has a limit on the amount of nodes that can be drawn in its context before performance takes a hit. Within our own testing 8,000 to 10,000 nodes is the limit of acceptable performance. If your database is larger than this amount, please see if you can filter out any nodes that are present in the data.

## Adjusting the wordcloud
The [d3-cloud](https://github.com/jasondavies/d3-cloud) library is used to generate the visualizations. To modify the rendered image, edit the lines of code below. The data for the function is returned from the backend api found [here](https://github.com/SPARC-FAIR-Codeathon/SPARClink/blob/main/MLDataIndexingEngine/README.md#word-cloud).
```javascript
var layout = cloud()
  .size([300,300])
  .words(
    keywords.map(function (d, i) {
      return { text: d, size: 10 + values[i] * 90 };
    })
  )
  .padding(5)
  .rotate(0)
  .font("Asap, Verdana, Arial, Helvetica, sans-serif")
  .fontSize(function (d) {
    return d.size;
  })
  .on("end", draw);
```
# ML Data Indexing Engine

This module is used to run additional analytics on the search terms for a more smarter SPARClink experience. There are three submodules in this section that references the three key features that are used by SPARClink.

The machine learning engine provides smart search features to the end user of our demo portal. This module also contains algorithms that learn vector embeddings of the descriptors of the elements present in the SPARClink database. Based on these vector embeddings the algorithms compute the similarity between the vectors representation of each word in the vocabulary with the vector representing the whole dataset and finds keywords that would describe the resource. We show these keywords based on their relevance in our word cloud. We have multiple options for searching for keywords in our codebase each providing a different level of performance in terms of time complexity and explainability. We use the [Symspell](https://github.com/wolfgarbe/symspell) algorithm trained on the dataset present in the scikit learn package as well as the vocabulary built using SPARClink database. We use delete-only edit candidate generation for generating different combinations of spelling errors and use both character level embedding and word embedding for recommending the most probable correct spelling. The output of the spell correction algorithm is used to generate sentence level embedding and is then compared with the embeddings of different discriptors of the items in the dataset. We obtain a ranking of all the items in the dataset based on their similarity with the searched string. The top 10 are chosen to be shown on the frontend. 

## Keyword Extraction
This module is used to generate keyword using keyBERT pretrained model. It generates top 50 key words associated with the whole document. It also make use of Maximal marginal relevence algorithm to pick key words that have higher distance among them. This is to ensure diversity among the chosen key words.

## SPARC Search
This module is primarily used for the Top 10 recommended related items in the bottom left panel. Filter keywords are passed to the function via a https request and the server will return a list of most relevant publications, protocols and datasets. The second query parameter is used for the smart filter feature that can be found in the [`smart_filter`](https://github.com/SPARC-FAIR-Codeathon/SPARClink/tree/smart_filter) branch. A list of what we think is the right word is returned in the response.

```javascript
var axios = require('axios');
var data = JSON.stringify({
  "inputString": "Identification of peripheral neural cercuit",
  "fullModel": false,
  "recommendation": true
});

var config = {
  method: 'post',
  url: 'https://megasanjay.pythonanywhere.com/sparcsearch',
  headers: { 
    'Content-Type': 'application/json'
  },
  data : data
};

axios(config)
.then(function (response) {
  console.log(JSON.stringify(response.data));
})
.catch(function (error) {
  console.log(error);
});

```

## Word Cloud
This module is used to generate a list of keywords for the word map visualization that is shown in the left control panel. The api receives a list of sentences and  analyzes the content before suggesting the list of best keywords ranked in descending order. On the front end we use the top 20 terms returned to generate the interactive wordcloud. 

```javascript
var axios = require('axios');
var data = JSON.stringify({
  "sentences": [
    "Annotated heart scaffold for pig available for registration of segmented neural anatomical-functional mapping of cardiac neural circuits.",
    "Data from the innervation of intact mice hearts and the mapping of parasympathetic and sympathetic neural circuits which control heart rate. This data set identifies the cholinergic and noradrenergic neurons which project to the sinoatrial node."
  ]
});

var config = {
  method: 'post',
  url: 'https://megasanjay.pythonanywhere.com/wordcloud',
  headers: { 
    'Content-Type': 'application/json'
  },
  data : data
};

axios(config)
.then(function (response) {
  console.log(JSON.stringify(response.data));
})
.catch(function (error) {
  console.log(error);
});
```
# Keyword Extractor

Description - 
Uses keyBERT to perform key word extraction. This code uses maximal marginal relevance to change the diversity among the keywords.
`main(number of keywords)` takes `int` argument signifying numbe of keywords and return a dictionary in json format where `key` are keywords and `value` are importance factor.

## Requirements
```
!pip install keybert[all]
```
## Demo:
Top 60 key words generated by keeping ```mmr = .3``` within the keyBERT function.
```python 
[('cardiometabolic', 1.5051),
 ('heart', 1.5151000000000001),
 ('piezo1', 1.516),
 ('electrophysiology', 1.52),
 ('macrophage', 1.5461),
 ('myocardial', 1.6327999999999998),
 ('hypertension', 1.6348),
 ('neuroepigenome', 1.6387),
 ('neurobiological', 1.6447999999999998),
 ('vivo', 1.686),
 ('receptor', 1.7062),
 ('neuroscientific', 1.7069),
 ('cnt', 1.7097000000000002),
 ('rhythmogenic', 1.7109),
 ('sympathetic', 1.7194),
 ('synaptic', 1.7195999999999998),
 ('neurophysiology', 1.7406000000000001),
 ('neuroinflammatory', 1.7571),
 ('neurogenic', 1.7622),
 ('noradrenergic', 1.7868),
 ('autonomic', 1.8639000000000001),
 ('synuclein', 1.8853000000000002),
 ('calcitonin', 1.8932),
 ('etiologic', 1.8958),
 ('gcamp6s', 1.9038),
 ('neurological', 2.019),
 ('nerve', 2.0361000000000002),
 ('ischemia', 2.0863),
 ('channelrhodopsin', 2.1014),
 ('neurochemical', 2.125),
 ('tachyarrhythmias', 2.1407999999999996),
 ('ventricle', 2.1545),
 ('glutamatergic', 2.1639999999999997),
 ('coughing', 2.1729),
 ('asthma', 2.2398000000000002),
 ('lidocaine', 2.2516000000000003),
 ('pathophysiology', 2.344),
 ('cgrpergic', 2.3728),
 ('bungarotoxin', 2.4015),
 ('gabaergic', 2.4482999999999997),
 ('pharmacodynamic', 2.4634),
 ('biomarker', 2.4732),
 ('angiotensin', 2.565),
 ('cardioprotection', 2.6090000000000004),
 ('atrioventricular', 2.6842),
 ('nociceptors', 2.6895000000000002),
 ('neuroinflammation', 2.6972000000000005),
 ('microglia', 2.7223),
 ('biomarkers', 2.7470999999999997),
 ('fibrillation', 2.752),
 ('neuro', 3.1549),
 ('neurotransmission', 3.375),
 ('neuronal', 3.3883),
 ('cardiac', 3.573),
 ('acetylcholine', 3.6888),
 ('cholinergic', 3.7815),
 ('antiarrhythmic', 3.7904000000000004),
 ('neurocardiac', 3.8373),
 ('vasoactive', 3.9165000000000005),
 ('neurotransmitter', 5.672)]
# Word Cloud
This module uses frequency of occurrence for all words present in the vocabulary and ranks them based on their frequency. It also provides a method that takes in these words and their frequencies to generate word cloud in the shape of the image `sparc_lung.jpg`. Higher the frequency larger the size of the word in the word cloud.
 
 ## Demo
 Using human lungs shape silhouette image to generate mask
 
 <p align="center">
  <!--<img src="https://github.com/SPARC-FAIR-Codeathon/SPARClink/blob/main/MLDataIndexingEngine/WordCloud/sparc_lungs1.png" />-->
  <img src="https://github.com/SPARC-FAIR-Codeathon/SPARClink/blob/main/MLDataIndexingEngine/WordCloud/sparc_lungs1.png" width="500"/>
  <br/>
  <span> Lung shaped word cloud </span>
</p>
# SPARC Search and Recommender

# Contents:
* [Description](#Description)
* [Requirement](#Requirement)
* [Steps for Usage](#Steps-for-Usage)
* [Demo](#Demo)
* [References](#References)

## Description:
This is a recommendation system that was built with the aim to make searching for documents using keywords simpler and user friendly. The two mechanisms used in this module are explained below. 

First is a word recommendation or spelling recommendation which takes a given string looks for the most probable word in the vocabulary. 

The second part uses the spelling recommendation to look for publications and datasets present in the Firebase real-time daatabase to get most relatable articles with respect to the search terms. This is based on a word embedding generated using a word2vec model. GloVe is an unsupervised algorithm that attempts to obtain high-dimensional vector representations of words using global word-word co-occurrence. This has been shown to encode interesting semantic information and to encode semantically related words nearby one another. (the code also provides option of using GloVe embedding). 

Word2vec model generates a vocabulary and vector representation for each article based on their title and description. When a string is serched its vector representation is generated and the cosine similarity between that vector and the vector representation of each article title/description is measured.

We pick the top 10 of these to show on the frontend.

For further reading please refer to references.

## Requirement:
To run this function standalone follow the instructions below. These commands are already handled if you are using the SPARClink conda environment.
```python
!python -m nltk.downloader stopwords
!python -m pip install -U symspellpy
!pip install --upgrade gensim
```

## Steps for usage:
The main function takes three inputs: 
1. `string`: This parameter takes in a string. The output from the search bar goes here.
2. `full_model`: If 'True' then uses Glove model. (Read below for further instructions)
3. `recommendation`: If `True` then use the word2vec model to return the article recomendation system.


```python
SparcSearch(string, full_model = False, recomendation = True)
```


Returns the publication and dataset ids in order of relevance to search field (empty list in case of `recommendation = False`) and spelling recommendation.
1. To switch off the recommendation system for papers set `recommendation` to `False`. This will then only run the spelling recommendation. 

2. If `recommendation = True` the code will perform spelling recommendation followed by paper recommendation.

3. To train the model using the Glove dataset :
   1. Download https://nlp.stanford.edu/data/glove.6B.zip 
   2. `!unzip glove.6B.zip`
   3. `get_glove2wv()`
   4. Set `full_model` to True

3. To run the code only on the words and tags associated with the dataset keep full model = Flase


## Demo:
```
python Sparcsearch.py "neoral cercuit" False True
```
```python
output :
['103389fnins201900897',
  '101016jtins202009011',
  '101016jconb201911011',
  '1014797FVDN2224',
  '101016jneuron202007010',
  '101093europaceeuy134',
  '101038s41593-021-00828-2',
  '101016jconb201804006',
  '101007s10827-019-00709-5',
  '101038s41587-019-0198-8',
  '101038s41586-020-2474-7',
  '101038s41586-021-03413-6',
  '103389fncel201800469',
  '101210endocrbqab087',
  '101016jexpneurol2020113256',
  '101002mds27321',
  '103389fphys2020563372',
  'https:dxdoiorg1026275yo5c-etlo',
  '101016jcell201705034',
  '101152ajpheart006352019',
  'https:dxdoiorg1026275yum2-z4uf',
  'https:dxdoiorg1026275eyik-qjhm',
  '101128IAI00928-19',
  '101523JNEUROSCI2158-192020',
  'https:dxdoiorg1026275dv4h-izxs',
  '101186s42234-019-0030-2',
  '101016jneuron202009031',
  '101038s41575-020-0271-2',
  '101016jjacc201910046',
  '101016jconb201712011',
  'https:dxdoiorg1026275dn1d-owj9',
  'https:dxdoiorg1026275iprt-7m5c',
  'https:dxdoiorg1026275dqpf-gqdt',
  'https:dxdoiorg1026275higx-q8hs',
  'https:dxdoiorg1026275uxjv-kbrz',
  '101152jn004422020',
  'https:dxdoiorg1026275zpju-kpjd',
  '101016jcell202005029',
  '103389fnins2020619275',
  '101523JNEUROSCI0743-192019',
  '101016jomtm202102012',
  '101146annurev-cancerbio-030419-033413',
  ...]
```
Top 15 paper titles :

```
A Student’s Guide to Neural Circuit Tracing
Neural Circuits of Interoception
Parenting — a paradigm for investigating the neural circuit basis of behavior
Neural Mechanisms and Therapeutic Opportunities for Atrial Fibrillation
Viral vectors for neural circuit mapping and recent advances in trans-synaptic anterograde tracers
Neural ablation to treat ventricular arrhythmias
An amygdala-to-hypothalamus circuit for social reward
The neural circuits of thermal perception
Emerging techniques in statistical analysis of neural data.
Next-generation interfaces for studying neural function
Microbes modulate sympathetic neurons via a gut-brain circuit
An amygdala circuit that suppresses social engagement.
Methods for Three-Dimensional All-Optical Manipulation of Neural Circuits
The Effects of Estrogens on Neural Circuits That Control Temperature
Targeted activation of spinal respiratory neural circuits
...
```
## References:
1. For Symspell look below:
    1. https://github.com/wolfgarbe/symspell
    2. https://wolfgarbe.medium.com/1000x-faster-spelling-correction-algorithm-2012-8701fcd87a5f
2. [GloVe](https://nlp.stanford.edu/projects/glove/)
# metadata_extraction.py

This module contains API implementations to communicate with SPARC Pennsieve and https://www.protocols.io/workspaces/sparc.

## SPARC Pennsieve

### get_list_of_datasets_with_metadata([])

Return all the dataset from SPARC Pennsieve. 

Return object:
``` python
[..., {
    'datasetId': 64, 
    'version': 4, 
    'name': 'Quantified Morphology of the Pig Vagus Nerve'
    'model': 'award', 
    "description": "This dataset contains recordings of compound nerve action potentials from the pelvic nerve as well as the pudendal nerves and its branches in response to electrical stimulation on the epidural surface of the sacral spinal cord in anesthetized cats.",
    'properties': {
        'description': 'Experiments to map physiological functions of autonomic nerves and the continued advance of bioelectronictherapies are limited by inadequate activation or block of targeted nerve fibers and unwanted co-activation orblock of non-targeted nerve fibers. More fundamentally, the relationship between applied stimuli and the nervefibers that are activated or blocked, how this relationship varies across individuals and species, and how theserelationships can be controlled remain largely unknown. We will develop, implement and validate an efficientcomputational pipeline for simulation of electrical activation and block of different nerve fiber types withinautonomic nerves. The pipeline will include segmentation of microanatomy from fixed nerve samples, three-dimensional finite-element models of electrodes positioned on nerves, and non-linear cable models of differentnerve fiber types, enabling calculation of quantitative input-output maps of activation and block of specific nervefibers. As key benchmarks of pipeline development and for the proposed analysis and design efforts, we willimplement models of the cervical (VNc) and abdominal (VNa) vagus nerves in rat, in a SPARC-identified animalmodel, and in human. The VNc is an excellent test bed as it contains a broad spectrum of nerve fiber types,there are experimental data to facilitate model validation, and there are multiple applications of VNc stimulationwhere a lack of fiber selectivity limits the therapeutic window. The VNa is an excellent complement to the cervicalVNc, as a prototypical autonomic nerve of a size comparable to many of the small autonomic nerves targetedby SPARC projects. We will use the models that emerge from the pipeline to achieve analysis and design goalsto address critical gaps identified as SPARC priorities. Specifically, we will quantify of the effects of intra-speciesdifferences in nerve morphology on activation and block by building individual sample-specific models for eachnerve and specie. These models will also be used to quantify inter-species differences in nerve fiber activationand block and to identify electrode designs and stimulation parameters that produce equivalent degrees ofactivation and block across species. We will combine the resulting models with engineering optimization todesign approaches to increase the selectivity and efficiency of activation and block of different nerve fiber types.The outcomes will be a pipeline for modeling autonomic nerves, electrode geometries, and stimulationparameters, as well as tools that address the limitations of nerve stimulation selectivity and efficiency that hinderthe continued advance of physiological mapping studies and the development of bioelectronic therapies.', 
        'id': 'db090035-8a57-4474-8c3d-2536dee9499f', 
        'principal_investigator': 'GRILL, WARREN M.', 
        'title': 'MODELING ACTIVATION AND BLOCK OF AUTONOMIC NERVES FOR ANALYSIS AND DESIGN', 
        'award_id': 'OT2OD025340'}, 
    'versionPublishedAt': '2020-10-01T15:41:57.651749Z', 
    'datasetDOI': 'https://dx.doi.org/10.26275/maq2-eii4', 
    'tags': [
        'vagus nerve stimulation', 
        'neural anatomy', 
        'vagus nerve morphology', 
        'autonomic nervous system'], 
    'contributors': [
        {'firstName': 'Nicole', 'middleInitial': 'A', 'lastName': 'Pelot', 'degree': 'Ph.D.', 'orcid': '0000-0003-2844-0190'}, 
        {'firstName': 'Gabriel', 'middleInitial': 'B', 'lastName': 'Goldhagen', 'degree': None, 'orcid': None}, 
        {'firstName': 'Jake', 'middleInitial': 'E', 'lastName': 'Cariello', 'degree': None, 'orcid': None}, 
        {'firstName': 'Warren', 'middleInitial': 'M', 'lastName': 'Grill', 'degree': 'Ph.D.', 'orcid': '0000-0001-5240-6588'}], 
    'originatingArticleDOI': [], 
    'protocolsDOI': ['dx.doi.org/10.17504/protocols.io.6bvhan6']}, ...]
```

## Protocols.io API

### parsing_protocols(authorization_key)
- `authorization_key` : API key obtained from protocols.io

Returns all the protocols in the SPARC workspace. 

Return object:
``` Python
[...{
  "authors" : [ {
    "affiliation" : "UCLA",
    "blocked_by_you" : false,
    "blocked_you" : false,
    "hide_following" : false,
    "image" : {
      "placeholder" : "/img/avatars/011.png",
      "source" : "/img/avatars/011.png",
      "webp_source" : ""
    },
    "is_verified_user" : false,
    "name" : "Pradeep  Rajendran",
    "username" : "pradeep-rajendran",
    "verified" : 0
  }, {
    "affiliation" : "UCLA",
    "blocked_by_you" : false,
    "blocked_you" : false,
    "hide_following" : false,
    "image" : {
      "placeholder" : "",
      "source" : "",
      "webp_source" : ""
    },
    "is_verified_user" : false,
    "name" : "John Tompkins",
    "username" : "",
    "verified" : 0
  }, {
    "affiliation" : "UCLA",
    "blocked_by_you" : false,
    "blocked_you" : false,
    "hide_following" : false,
    "image" : {
      "placeholder" : "",
      "source" : "",
      "webp_source" : ""
    },
    "is_verified_user" : false,
    "name" : "Kalyanam Shivkumar",
    "username" : "",
    "verified" : 0
  } ],
  "doi" : "10.17504/protocols.io.2n9gdh6",
  "title" : "Pig- Heart Neuonal and Fiber Immunocytochemistry",
  "url" : "https://www.protocols.io/view/pig-heart-neuonal-and-fiber-immunocytochemistry-2n9gdh6"
}...]
```


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
 - OS: [e.g. Windows, macOS, Linux]
 - Browser [e.g. chrome, safari]

**Additional context**
Add any other context about the problem here.
# External APIs

This folder contains the implementations for APIs to communicate with external resources.

## NIH_NCBI.py
API implementations to communicate with [NIH RePORTER](https://api.reporter.nih.gov/) and [NCBI](https://www.ncbi.nlm.nih.gov/home/develop/api/).

### NCBI API uses

#### getCitedBy (id_type, id)
- `id_type`: Type of the given `id`. 'pm_id' for PubMed articles or 'pmc_id' for PubMed Central articles. 
- `id`     : Identifier for the article.

Returns a dictionary of the articles that cite the given article, with doi as the key.
```
title       : Title of the paper
journal     : Name of the journal
year        : Publication year
author_list : Names of authors
doi         : DOI of the paper
pm_id       : PubMed id (if available)
pmc_id      : PubMed Central id of the paper
```

#### getPublicationsWithSearchTerm (search_term)
- `search_term`: Search term to look for in PubMed Central.

Returns a dictionary of the publications that matches the given search term, with the doi as the key. 
```
title       : Title of the paper
journal     : Name of the journal
year        : Publication year
author_list : Names of authors
doi         : DOI of the paper
pm_id       : PubMed id (if available)
pmc_id      : PubMed Central id of the paper
```
e.g. To find all the articles that mention the doi `10.26275/DUZ8-MQ3N`
``` python
getPublicationsWithSearchTerm('"10.26275/DUZ8-MQ3N"')
```
To find the article with the doi `10.26275/DUZ8-MQ3N`
``` python
getPublicationsWithSearchTerm('10.26275/DUZ8-MQ3N[doi]')
```

### NIH RePORTER API uses

#### generateRecord(getProjectFundingDetails (project_no)):
- `project_no`: Project number (also known as the award number) of NIH funding associated with Pennsieve datasets.

Returns a dictionary of the funding details of the grants, with the project number as the key. A project can have multiple grant applications (e.g. extensions, ammendments)
```
appl_id   : Application identifier
institute : Name of the organization that received funding
country   : Country of the organization
amount    : Amount received
year      : Fiscal year
keywords  : Keywords of the project topic.
```

#### getPublications (appl_id)
- `appl_id` : Application identifier of a grant.

Returns a dictionary of all the publications associated with the grant, with doi as the key.
```
title       : Paper title
journal     : Name of the journal
year        : Published year
author_list : Names of authors.
url         : URL to the paper
pm_id       : PubMed id of the paper
doi         : Paper DOI
```






