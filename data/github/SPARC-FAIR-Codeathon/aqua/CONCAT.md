<p align="center">
  <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/logo_aqua-1.jpg" alt="interface" width="420" height="250"> 
  <br/> 
  </img>
</p>

## Table of Contents

* [About AQUA](#about-aqua)
* [The Problem](#the-problem)
* [AQUA solution](#aqua-solution)
* [How it works](#how-it-works)
* [How to install](#how-to-install)
* [Documentation](#documentation)
* [API](#api)
* [Developers](#developers)

## About AQUA

AQUA (Advanced QUery Architecture for the SPARC Portal) is an application that aims at improving the search capabilities of the [SPARC Portal](https://sparc.science/). In particular, we are looking to make the search engine smarter at reading and understanding user input as search keywords. We also enhance the result display feature of the SPARC Portal by making it more user-friendly and providing users with more sophisticated result filtering and sorting options. Our end goal is to improve exponentially the visibility of the SPARC datasets. This in turn will benefit the SPARC community as a whole since their datasets will be more discoverable for reuse and subsequent collaboration. This project was created during the 2021 SPARC FAIR Codeathon.

## The Problem

Currently, the search feature of the SPARC Portal is very limited: 

1) It does not recognize nearby words (typos and close-matches) or synonyms.

2) The result display is limited. E.g.: Limited result filtering and sorting (only by Published Date or Alphabetical Ordered Titles).

## AQUA solution

__1) Apply Artificial Intelligence tools__ (Natural Language Processing) to the processing of users’ search keywords and to the implementation of predictive typing (suggestion-based typing). 

- In details, in addition to lemmatization, other NIH tools (e.g: NIF Ontology) will be used to derive origins of words and make autocomplete suggestions for users as they type. This will help AQUA standardize various user inputs and return the most datasets possible that match the search keywords.
- AQUA also fixes typos and close matches and suggests corrected search keywords.

__2) Enhance the current result display by:__

- Bolding/highlighting matched texts in results for easy lookup

- Add a more sophisticated Dataset results sorting and filtering functionality (based on Relevance, Date of Publication, and other customized filtering) to the current portal.

- Add a “Notify me when related datasets are published”. This will allow users to enter their email to be stored by the SPARC Portal for future alerts. 

## How it works

<p align="left">
  <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/workflow_new.jpg" alt="interface" width="900" height="500"> 
  <br/> 
  </img>
</p>
 
## How to install

**Step 1:** Git clone the AQUA project: `git clone https://github.com/SPARC-FAIR-Codeathon/aqua.git`

**Step 2:** Go into the `aqua` directory and run the following commands:

```
# install dependencies
$ yarn install

# serve with hot reload at localhost:3000
$ yarn dev

# build for production and launch server
$ yarn build
$ yarn start

```

## Documentation

For a detailed user documentation of our application, please visit [:arrow_forward: Documentation](https://github.com/SPARC-FAIR-Codeathon/aqua/blob/main/Documentation/Documentation.md).

## API

To read the AQUA API refer to: [AQUA API](https://github.com/SPARC-FAIR-Codeathon/aqua/tree/main/aqua_docker/docs/api).

## Developers

- [Tram Ngo](https://github.com/tramngo1603) (Lead)
- [Laila Rasmy](https://github.com/lrasmy) (Sysadmin)
- [Niloofar Shahidi](https://github.com/Niloofar-Sh) (Technical writer)
- [Yuda Munarko](https://github.com/napakalas) (Sysadmin)
- [Xuanzhi](https://github.com/marcusLXZ) (Front-end)
# NotifyMe
### An email notification functionality

## Main Purpose

The NotifyMe option is to send emails that summarize search results against exact keywords. NotifyMe sends the email just once at least one exact hitting match exists. Therefore, NotifyMe can be used for:

     1.	Notify users when a dataset gets published against keywords that don’t retrieve any results yet.
     2.	Emailing the current search results in a tabular format, which can be found helpful for users.

Additionally, NotifyMe stores all requests in an SQLite database, which can be further analyzed by the SPARC team to understand the search pattern and get more insights on the demand. For example, the SPARC team can find out the most common keywords searched with no existing matches and decide on actions to fulfill such needs. 

## How it works

 We can summarize NotifyMe actions as follow:
 
       1.	Add email requests with keywords 
 
       2.	Scan for existing search hits and send email
 
       3.	Moving pending requests to a waiting list, that’s scanned daily
 
       4.	Moving fulfilled requests to an archived list
 
       5.	Any failed requests (that already have matching hits) will remain in the waiting list for one month,
          during which NotifyMe will try to send the email on daily basis. Afterward, if the email still failing,
          it will move to the archived list with a final failed status for efficiency.
 


<p align="left">
  <img src="./NotifyMe.jpeg" alt="interface" width="900" height="550"> 
  <br/> 
  </img>
</p>

## How to run

1. First, update the [properties.ini](./properties.ini) with the required information like sending email password and the scicrunch api-key

2. Run [notifyme_api.py](./notifyme_api.py), in order fetch the user email and search keywords.
 
 an example call:  http://localhost:5432/aqua/notifyme?email="<email>"&keywords="<keywords>"

3. In order to schedule the keywords search and sending emails, you need to run [notifyme_sched.py](./notifyme_sched.py)

   The current setting is scheduling emails to be sent daily at 2 am

4. The request are saved in a SQLite database. the description of the database tables is available in [NotifyMe_Database.pdf](./NotifyMe_Database.pdf)
 
5. A sample analytics visualization can run through [NotifyMe_analytics_visual.ipynb](https://nbviewer.jupyter.org/github/lrasmy/aqua/blob/main/NotifyMe/NotifyMe_analytics_visual.ipynb)
 

## Required Packages
- configparser
- flask_restplus
- numpy
- pandas
- schedule
- smtplib
- sqlite3

and plotly express for visualization examples
# Application name:

<p align="center">
  <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/logo_aqua-1.jpg" alt="interface" width="500" height="300"> 
  <br/> 
  </img>
</p>

# Contents:
* [Analysis of the existing SPARC search portal](#chart_with_upwards_trend-analysis-of-the-existing-sparc-search-portal)
* [AQUA objectives for the new SPARC search portal](#bulb-aqua-objectives-for-the-new-sparc-search-portal)
* [AQUA Components](#iphone-aqua-components)
  * [AQUA Backend](#1-aqua-backend)
  * [AQUA UI](#2-aqua-ui)
* [How to use AQUA?](#information_desk_person-how-to-use-aqua)
* [Installation](#hammer_and_wrench-installation)
* [Testing](#mag_right-testing)
* [Examples](#round_pushpin-examples)
* [Ideas?](#speech_balloon-ideas)


# :chart_with_upwards_trend: Analysis of the existing SPARC search portal

**1. Limited search feature of the SPARC Portal:** It does not recognize nearby words (in case of **typos**). As an example, if we type "rattis" (typo) instead of "rattus", it does not recognize it or give any suggestion (Fig 1).
<br/>
<p align="center">
   <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/rattis_current_result.jpg" width="800" height="450"></br>
  <i>Fig 1. Search results for a typo (rattis) on the SPARC portal. </i>
</p>
  
   <br/>
   
**2. Vague result display:** You need to enter the exact correct keywords in the search bar (Fig 2) and yet, it does not **bold/highlight the search keywords** among the search results.

<p align="center">
   <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/rattus_current_result.jpg" width="800" height="450"></br>
  <i>Fig 2. Search results for a correct keyword (rattis) on the SPARC portal. Keywords are not bolded/highlighted in the results.</i>
</p>
   
**3. Limited result filterings:** The website currently refines the results by either being "Public" or "Embargoed" (Fig 3).

**4. Limited result sorting:** The website currently sorts the results by either "Title" (listed alphabetically) or "Published date" (Fig 3).

<p align="center">
   <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/rattus_filters_1.JPG" width="800" height="450"></br>
  <i>Fig 3. Limited filterings for the results display on the SPARC portal.</i>
</p>
   <br/>
   
**5. No option for Email alerts :** In the event of No results being found for the queries, there is no way for the user to get informed on whether a **new dataset/resource** (related to their query) has been published on the SPARC portal or not.
<br/>

# :bulb: AQUA objectives for the new SPARC search portal
Specific features of AQUA are listed below:
<br/>
* :point_right: __Query refinement__
   * Auto-completion:<br/>
      Based on the term, our tool automatically completes the queries if it partially/completely matches any keywords. It then sends the selected keyword to AQUA backend.
   * Suggestion:<br/>
      If no exact matches are found, it finds close-matches and suggests them to the user with popping up the phrase: *"Did you mean ...?"*. Otherwise, it will send the raw and uncorrected query to AQUA backend.

* :point_right: __Results filtering__
   * Sort by:<br/>
      When the results for the query are displayed, user will have the option of sorting them based on the *Relevance*, *Date published*, and *Alphabetical order*.
    * Filter by:<br/>
      The results can also be filtered based on *Keyword*, *Author*, *Category*, and *Publication date*.
    * Matched text bolded/highlighted
    
    
* :point_right: __Alert me about future datasets related to search query__

     At the end, if no results are returned by the AQUA backend, our tool asks the user if they want to get notified when a related resource is published or not. For a given email address, the tool checks for its validity and then stores it using SQLite. Thereafter, it will check for any updated/uploaded related resource on the SPARC portal everyday at 2AM EDT. In case of the requested resource availability, it sends a notification email to the registered user. 
<br/>
<p align="center">   
  <b>The AQUA pipeline:</b>
</p>
<br/>
<p align="center">
   <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/aqua_pipeline.jpg" width="550" height="1000"></br>
  <i>Fig 4. AQUA pipeline including three major sections: Query refinement, Results filtering, and Notify me.</i>
</p>
   <br/>
   


# :iphone: AQUA Components
The AQUA application integrates Python libraries, data mining tools, a SQL database engine, and Document Object Model (DOM) API to mimic an environment similar to the SPARC search portal with an improved functionality in multiple ways. In general, the AQUA platform consists of a presentation layer as the User Interface (UI) (referred to as frontend) and a server-side data-access layer (referred to as backend). The **AQUA UI** and the **AQUA Backend** bridge between the user and the Knowledge Management Core (K-Core) database. K-Core is the SPARC knowledge graph database.


## 1. AQUA Backend

The AQUA backend includes querying the K-Core database for information, delivering data to the frontend, and processing any logic that the AQUA UI requires. The main tools utilised for the AQUA Backend are Python (Jupyter Lab), Docker, SQLite, and SciGraph.
<br/>
The AQUA backend focuses on two main features:
<br/>

:sparkles: __Behind the scenes of AQUA's Query refinement__ 

<p align="center">
   <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/blob/main/src/assets/images/Suggestion_AutoComplete.jpg" alt="interface" width="780" height="700"></br>
  <i>Fig 5. Query refinement by Auto-completion/Suggestions.</i>
</p>
<br/>

AQUA utilises SciGraph for auto-completion and suggestion. However, we found that SciGraph’s suggestions do not deal with query problems such as error spelling and continuous script (*scriptio continua*). Therefore, we have added a new auto-correction feature to segment queries with missing spaces and fix error spelling by creating a pipeline to [SymSpellPy](https://pypi.org/project/symspellpy/). The auto-correction result is combined with the suggestion results and then executed as the final query search terms.

(To read more, please visit: ["Query refinement" Readme](https://github.com/SPARC-FAIR-Codeathon/aqua/tree/main/Documentation/QueryRefinement.md))  

<br/>
<br/>

:sparkles: __Behind the scenes of AQUA's Email notification__ 
<br/>
<p align="center">
   <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/blob/main/src/assets/images/NotifyMe.jpeg" alt="interface" width="780" height="500"></br>
  <i>Fig 6. "Notify me" module.</i>
</p>
 <br/>
The "Notify Me" option is to send emails that summarize search results against exact keywords. "Notify Me" module sends the email just once at least one exact hitting match exists. Therefore, "Notify Me" can be used for:<br/>
<br/>

1.	Notify users when a dataset gets published against keywords that do not retrieve any results yet.
2.	Emailing the current search results in a tabular format, which can be found helpful for users.
Additionally, the "Notify Me" module stores all requests in an SQLite database, which the SPARC team can further analyze to understand the search pattern and get more insights on demand. For example, the SPARC team can find the most common keywords searched with no existing matches and decide to fulfil such needs. 

(To read more, please visit: ["Notify me" Readme](https://github.com/SPARC-FAIR-Codeathon/aqua/tree/main/Documentation/NotifyMe.md))
<br/>

## 2. AQUA UI
<br/>

AQUA UI receives the user's queries, formulates them, and transfers to the AQUA Backend module. When the response from the AQUA backend is received, the AQUA UI interprets it and displays the content on the screen. Like the SPARC portal web application, the AQUA UI is implemented by the HTML-CSS-JS trio using: [VueJS](https://vuejs.org/) and [NuxtJS](https://nuxtjs.org/). Nuxt.js is an upper-level framework that is built over Vue.js to design and create highly advanced web applications [[1](https://cubettech.com/resources/blog/nuxt-js-and-vue-js-reasons-why-they-differ-and-when-do-they-combine/)].



<p align="center">
   <img src="https://github.com/Niloofar-Sh/aqua/blob/main/src/assets/images/vuejsTOnuxtjs.jpg" alt="interface" width="500" height="300"></br>
  <i> Vue.js and Nuxt </i>
</p>
   <br/>
   
# :information_desk_person: How to use AQUA? 

How to use the 7 features added to the existing SPARC Portal Search engine by AQUA:

#### 1. Predictive search typing
AQUA provides autocompletion for user's query as they type. This feature is powered by training data from the NIF Ontologies and Scigraph. To avoid too many results being returned that can slow down the application, we only show autocompletion after users type 3 letters and more. (Figure 7)

<p align="center">
   <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/UI-autocompletion.png" alt="interface" width="850" height="580"></br>
  <i>Fig 7. Predictive typing interface.</i>
</p>

#### 2. Advanced search options

There are currently 2 options for user's search query: "Exact match" or "Any of the words match". The default is "Any of the words match". If users want to return datasets for their exact search phrase, they can do that by clicking on "Advanced search" under the search box (Figure 8).

<p align="center">
   <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/advanced-search.png" alt="interface" width="850" height="500"></br>
  <i>Fig 8. Advanced search interface.</i>
</p>

#### 3. Advanced Sorting
The existing SPARC Portal allows sorting based on dataset titles (alphabetically) and by published date. AQUA adds a "Relevance" sorting criterion that returns results based on how relevant the results are to their search query. This is set as the default sorting option (Figure 9).

<p align="center">
   <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/UI-sorting.png" alt="interface" width="700" height="500"></br>
  <i>Fig 9. Advanced sorting interface.</i>
</p>

#### 4. Advanced Filtering
The existing SPARC Portal only allows for filtering based on "Dataset status", which is either Published or Embargoed. Aqua adds more sophisticated filtering options: by "Published Date", "Keyword", "Author", and "Category". Users can filter datasets by one or several keywords, authors, and categories. Hit "Enter" after each "keyword", "author", or "category" in their respective box to register it. After the entries are registered, click "Apply" to filter dataset results (Filter 10).

<p align="center">
   <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/UI-filtering-feature.png" alt="interface" width="740" height="500"></br>
  <i>Fig 10. Advanced filtering interface.</i>
</p>

#### 5. Email notifications for new matched datasets
Users can opt in to receive emails about new datasets that match their search query. AQUA believes this is a much needed option for users to stay updated about their search and SPARC datasets. Simply click on "Create alerts" under the search box and enter an email. AQUA will trigger an email send when newly added dataset(s) that match the search query are published by SPARC. This is a one-time only email subscription. Options to be alerted more than once can be added in the future. (Figure 11)

<p align="center">
   <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/UI-email.png" alt="interface" width="780" height="500"></br>
  <i>Fig 11. Email Notification interface.</i>
</p>

#### 6. Bold matched texts in result display
When a dataset is returned, any matched text in the dataset title and description will be bolded for easy and convenient lookup (Figure 12). 

<p align="center">
   <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/highlighted.png" alt="interface" width="850" height="500"></br>
  <i>Fig 12. Bolded matched text in result display.</i>
</p>

#### 7. View type
AQUA adds view type to the existing SPARC Portal to enhance user experience with the website. The default option is List view, which is the SPARC Portal's existing view type. AQUA proposes to add a gallery view option in the future. (Figure 13)

<p align="center">
   <img src="https://github.com/SPARC-FAIR-Codeathon/aqua/raw/main/src/assets/images/UI-viewtype.png" alt="interface" width="700" height="500"></br>
  <i>Fig 13. View type interface.</i>
</p>

## :hammer_and_wrench: Installation

**Step 1**: Git clone the AQUA project by running the following command:
`git clone https://github.com/SPARC-FAIR-Codeathon/aqua.git`

**Step 2**: Make sure you have Yarn installed on your computer. Yarn is a package manager (an alternative to [npm](https://www.npmjs.com/)). Installation instructions for Yarn is detailed [here](https://classic.yarnpkg.com/lang/en/docs/install/#windows-stable). Specifically, you can choose which operating system you are installing Yarn for to get a custom installation instructions. 

**Step 3**: Go into the `aqua` directory and run the following commands:
```
# install dependencies
$ yarn install

# serve with hot reload at localhost:3000
$ yarn dev

# build for production and launch server
$ yarn build
$ yarn start
```

# :mag_right: Testing

:point_right: __Auto-completion experiment:__

Comparing **Scigraph** and **fast auto-complete**

We have 465142 of 1-gram and 2-grams extracted from NIFS Ontology and datasets.

```DataTest1``` (test_autocomplete_pure.json) consists of 200 words selected randomly from the n-grams. The selection criteria is:

* word with length between 3 and 15
* word does not contain number

Thereafter, we created a second dataset by changing one character in ```DataTest1``` randomly at position 2 until 15 with * , named ```DataTestWithTypo``` (test_autocomplete_typo.json).

The experiment is set to return 10 completion in maximum for each query.

:point_right: __Execution Time Analysis:__

<p align="center">
   <img src="https://github.com/Niloofar-Sh/aqua/blob/main/src/assets/images/execution_time.png" alt="interface" width="750" height="500"></br>
  <i>Fig 14. Execution time analysis.</i>
</p>

Example results plus the execution rates:

1) SciGraph + ```DataTest1``` : 0.8019 second
`
2) SciGraph + ```DataTestWithTypo```: 0.81809 second

3) fast-autocomplete + ```DataTest1``` : 0.03280 second

4) fast-autocomplete + ```DataTestWithTypo```: 0.07471 second


In general longer words will need longer execution times.

:point_right: __The number of completions:__

<p align="center">
   <img src="https://github.com/Niloofar-Sh/aqua/blob/main/src/assets/images/return_number.png" alt="interface" width="750" height="500"></br>
  <i>Fig 15. Number of completions returned.</i>
</p>

1) SciGraph + ```DataTest1``` : 
 ```
 'sinonasal': [],
 'acid-oxo': [],
 'chondrichthyes': ['chondrichthyes'],
 'turbellarian': 
  ['turbellarian platyhelminths']
  ```
2) SciGraph + ```DataTestWithTypo```:

```
'sinon*sal': [],
 'acid-o*o': [],
 'chondricht*yes': [],
 'turbel*arian': []
```
3) fast-autocomplete + ```DataTest1``` : 

```
'sinonasal': ['sinonasal', 'sinonasal papilloma', 'sinonasal undifferentiated'],
'acid-oxo':['acid oxo',
   'acid oxoanion',
   'acid oxo group',
   'acid oxoacid',
   'acid oxoanions',
   'acid oxoacids',
   'acid oxoglutarate',
   'acid oxoglutarate dehydrogenase',
   'acid oxoacyl-',
   'acid oxoacyl- and'],
'chondrichthyes': ['chondrichthyes'],
'turbellarian': ['turbellarian', 'turbellarian platyhelminths']
```

4) fast-autocomplete + ```DataTestWithTypo```:

```
'sinon*sal': ['sinonasal', 'sinonasal papilloma', 'sinonasal undifferentiated'],
 'acid-o*o': ['acid',
   'acid oxidase',
   'acid oxidation',
   'acid o-linked',
   'acid omega-hydroxylase',
   'acid or',
   'acid o-methyltransferase',
   'acid oxygenase',
   'acid optimum',
   'acid oxidoreductase'],
 'chondricht*yes':  ['chondrichthyes'],
 'turbel*arian': ['turbellarian', 'turbellarian platyhelminths']
```

__SciGraph__ returns a smaller number of completions. When there is a typo, SciGraph returns *almost zero completion*.
Longer words cause a reduce in the completion number. Typo tends to increase the number of completion for __fast-autocomplete__.

:triangular_flag_on_post: __In all cases, fast-autocomplete can return results.__
 
The code for the experiments can be found at [AQUA experiments](https://github.com/Niloofar-Sh/aqua/tree/main/experiment)
# :speech_balloon: Ideas?
To share your ideas, feedback, and comments, please either report an issue via Github Issue or contact any of our team members. Thank you!

# "Notify Me" module
__(Email notification functionality)__

## Main Purpose

The primary purpose of this functionality is to notify users whenever a new dataset is published against their search terms if nothing is available at the moment.

However, users can still use the same function to receive a summary table including basic information and links to all datasets currently matching their keywords.

Additionally, as requests are saved in a database, this information can be further accessed and analysed to get further content improvement.


## How to run

1. First, update the [properties.ini](https://github.com/Niloofar-Sh/aqua/blob/main/NotifyMe/properties.ini) with the required information like sending email password and the scicrunch api-key
 
 
2. Run [notifyme_api.py](https://github.com/Niloofar-Sh/aqua/blob/main/NotifyMe/notifyme_api.py), in order fetch the user email and search keywords.
   An example call:  
   ```
   http://localhost:5432/aqua/notifyme?email="<email>"&keywords="<keywords>"
   ```
3. In order to schedule the keywords search and sending emails, you need to run [notifyme_sched.py](https://github.com/Niloofar-Sh/aqua/blob/main/NotifyMe/notifyme_sched.py) <br/>

The current setting is scheduling emails to be sent daily at 2 AM EDT


 4. The request are saved in a SQLite database. the description of the database tables has come in [How it works](##how-it-works).
 
 5. A sample analytical visualization can run through [NotifyMe_analytics_visual.ipynb](https://nbviewer.jupyter.org/github/lrasmy/aqua/blob/main/NotifyMe/NotifyMe_analytics_visual.ipynb)
 

## How it works

<p align="left">
  <img src="https://github.com/Niloofar-Sh/aqua/blob/main/src/assets/images/NotifyMe.jpeg" alt="interface" width="900" height="550"> 
  <br/> 
  </img>
</p>

<br/>

We can summarize the "Notify Me" actions as follow:
1.	Add email requests with keywords.
2.	Scan for existing search hits and send email
3.	Moving pending requests to a waiting list that’s scanned daily
4.	Moving fulfilled requests to an archived list
5.	Any failed requests (that already have matching hits) will remain on the waiting list for one month, during which the "Notify Me" module will try to send the email daily. Afterwards, if the email still failing, it will move to the archived list with a final failed status for efficiency.

NotifyMe.db is the "Notify Me" SQLite database file. The database will be automatically created during the first call of the "Notify Me" module.  If one needs to drop and recreate the database for any reason, they can call the create_notifyme_db function from Notifyme_utils and set the clean option to True.

The database consists of 4 main tables:

1. __NEW_REGISTER__

| Column Name        | Description           | 
| ------------- |:-------------:| 
| entry_id     | A unique identifier for email requests (autoincrement integer) | 
| email      |   The user entered email (validated at the front end)    |   
| register_date |   The system date and time corresponding to the record creation, which is the time of the request initialization     |    
| keywords     | The user-entered search keywords | 
| status      |  All records in this table should show a ‘New’ status  |   
| last_modified |  In case the record get modified for any reason, this record is representing the last modification date and time for the corresponding record |  


2. __WAITING_LIST__

| Column Name        | Description           | 
| ------------- |:-------------:| 
| entry_id     | A unique identifier for email requests (automatically created in the new_register table) | 
| email      |   The user entered email (validated at the front end)    |   
| register_date |   The date and time of the request     |    
| keywords     | The user-entered search keywords | 
| status      | Request current status. Can be either ‘New’ if no hits still matching, or ‘Failed’ if the last attempt to send an email failed, the detailed error will be stored in Failed_emails |   
| last_modified | The date and time for the last modification of this record |  
| hits| The current number of matching hits exists against the search keywords. This should be ‘0’ for records remaining in ‘New’ status, and an actual number for ‘Failed’ records|
|failed_reqid|This is the reference for the latest corresponding Failed_emails record, showing the exact error that explains why this request failed|

3. __ARCHIEVED_LIST__

| Column Name        | Description           | 
| ------------- |:-------------:| 
| entry_id     | A unique identifier for email requests (automatically created in the new_register table) | 
| email      |   The user entered email (validated at the front end)    |   
| register_date |   The date and time of the request     |    
| keywords     | The user-entered search keywords | 
| status      | Request current status. Can be either ‘Sent’ for successfully sent emails,	‘Duplicate’ for the case in which the request identified earlier to be a duplicate request/entry,	‘Failed’ in case the email request raises an error consistently for more than one month. |   
| last_modified | The date and time for the last modification of this record. Should be corresponding to the time the email is sent if the status is ‘Sent’. |  
| hits| The number of matching hits sent against the search keywords. In case of failed requests, it will be the number of hits that exists at the time the record moved from the waiting_list table to here|
|failed_reqid|This is the reference for the latest corresponding Failed_emails record, showing the exact error that explains why this request failed|

4. __FAILED_EMAILS__

| Column Name        | Description           | 
| ------------- |:-------------:| 
|failed_reqid | A unique identifier for an error recorded against an email request (autoincrement integer)| 
|corresponding_entry_id | The corresponding entry_id of the email request will normally find either in the waiting_list table or the archieved_list table if we fail to send an email for more than a month|   
| register_date | The date and time of the request |    
|error_message | The detailed error message leading for email posting failure | 
|error_date | The system date and time when the error was triggered  |   

<br/>
<br/>


## A sample analytical visualization

```
import NM_analytics_utils as ut
ut.plot_most_frequent_wait()
```
<p align="left">
  <img src="https://github.com/Niloofar-Sh/aqua/blob/main/src/assets/images/not_matched.png" alt="interface" width="700" height="400"> 
  <br/> 
  </img>
</p>

<br/>

```
ut.plot_most_freq_search_term(10)
```
<p align="left">
  <img src="https://github.com/Niloofar-Sh/aqua/blob/main/src/assets/images/topTen.png" alt="interface" width="900" height="550"> 
  <br/> 
  </img>
</p>

<br/>

```
ut.plot_key_hits_pie()
```

<p align="left">
  <img src="https://github.com/Niloofar-Sh/aqua/blob/main/src/assets/images/numOfHits.png" alt="interface" width="700" height="400"> 
  <br/> 
  </img>
</p>

<br/>

```
ut.plot_key_timetomatch()
```
<p align="left">
  <img src="https://github.com/Niloofar-Sh/aqua/blob/main/src/assets/images/maxTime.png" alt="interface" width="750" height="450"> 
  <br/> 
  </img>
</p>

<br/>

## Required Packages
- configparser
- flask_restplus
- numpy
- pandas
- schedule
- smtplib
- sqlite3

and plotly express for visualization examples


<p align="center">
   <img src="https://github.com/Niloofar-Sh/aqua/raw/main/src/assets/images/Suggestion%26AutoComplete.jpg" alt="interface" width="780" height="700"></br>
  <i> Query refinement by Auto-completion/Suggestions. Purple box: suggestions, yellow box: auto-completion.</i>
</p>
<br/>

# What is the SPARC dataset metadata?
Metadata is the “Data about data”, i.e., additional information provided about datasets. The **SPARC dataset metadata** includes information such as title, description, techniques, as well as the number of the files, formats, licenses, etc. ([SPARC dataset metadata](https://staging.sparc.science/help/3vcLloyvrvmnK3Nopddrka#metadata)).
<br/>
<br/>

# What is the NIFS ontology?
NIF Standard ontology (**NIFS ontology**) is a neuroscience ontology that maintains an extensive set of terms and concepts important for the domains of neuroscience and biology ([NIFS Ontology](https://github.com/SciCrunch/NIF-Ontology)). This a set of community ontologies used by SPARC to annotate data and models.
<br/>
<br/>

# Sections of the "AQUA Query refinement" module



## Suggestions path (purple box):
If there is a typo or removed space between the words of a query, the ElasticSeach might return either no results or irrelevant results. In this case, we need a suggestion and auto-correction feature to improve the quality of the query. 
Merely using SciGraph is not sufficient because SciGraph returns alternative queries/suggestions without correcting the initial query. To improve this process, we have implemented an auto-correction pipeline along with SciGraph to correct the queries before giving suggestions. This includes an *Auto-correction n-gram model* and a Python library *symspellpy*. 
<br/>

### SciGraph

Represents ontologies and ontology-encoded knowledge in a Neo4j graph. SciGraph reads ontologies with owlapi and ingests ontology formats available to owlapi (OWL, RDF, OBO, TTL, etc) ([SciGraph](https://github.com/SciGraph/SciGraph)).

### Auto-correction n-gram model

In spelling correction task, an n-gram is a contiguous sequence of n letters from a given sample of text. An n-gram model is utilised to compare strings and compute the similarity between two words, by counting the number of similar n-grams they share. This technique is language independent. The more similar n-grams between two words exist the more similar they are ([Ahmed et al.](http://www.scielo.org.mx/pdf/poli/n40/n40a7.pdf)). 

### symspellpy

symspellpy is a Python port of SymSpell for spelling correction, fuzzy search and approximate string matching ([symspellpy](https://pypi.org/project/symspellpy/),[SymSpell](https://github.com/wolfgarbe/SymSpell)).

* __Word segmentstion__

``` 
word_segmentation(phrase, max_edit_distance=None, max_segmentation_word_length=None, ignore_token=None)
``` 

word_segmentation divides a string into words by inserting missing spaces at the appropriate positions misspelled words are corrected and do not affect segmentation existing spaces are allowed and considered for optimum segmentation.

* __Spelling correction__

```
lookup_compound(phrase, max_edit_distance, ignore_non_words=False, transfer_casing=False, split_phrase_by_space=False, ignore_term_with_digits=False)
```

Supports compound aware automatic spelling correction of multi-word input strings with three cases:

1. mistakenly inserted space into a correct word led to two incorrect terms
2. mistakenly omitted space between two correct words led to one incorrect combined term
3. multiple independent input terms with/without spelling errors <br/>

Find suggested spellings for a multi-word input string (supports word splitting/merging).

## Auto-completion path (yellow box):
It is an added feature to auto-complete the queries while the user is typing. The idea of auto-completion is to prevent typos occuring and to give a better user experience in the SPARC Portal. We have created an n-gram model for auto-completion and utilised a Python library *fast-autocomplete* ([fast-autocomplete](https://pypi.org/project/fast-autocomplete/)).

### Auto-completion model
The format of the n-gram model needs to be in the following format:

``` 
{
    phrase: [
        context,
        display value,
        count
    ]
}
``` 

* "phrase" can be 1-2 words. 
* The "context" is related to the context of words, for example Anatomy, chemical reactions, proteins, etc. 
* "display value" defines the standard display of the phrase based on the context.
* "count" is the appearance of the phrase in the SPARC dataset and the NIFS ontology.


### Fast auto-complete

The Elasticsearch's Autocomplete suggestor is not fast enough and does not do everything that we need. Consequently, we have utilised fast-autocomplete library which provides us with a much faster process (reducing the average latency from 900 ms to 30 ms).

Elasticsearch's Autocomplete suggestor does not handle any sort of combination of the words you have put in. For example Fast Autocomplete can handle ```brainstem neuron in rat``` when the words ```brainstem```, ```neuron```, ```in```, ```rat``` are seperately fed into it. While Elasticsearch's autocomplete needs that whole sentence to be fed to it to show it in Autocomplete results.

# Packages:

1. symspellpy
2. scigraph
3. fast-autocomplete
# AQUA - Experiment
## About
Experiments of autocomplete and spelling error correction for AQUA SPARC

## How to reproduce experiment results?
[Experiment](Experiment.ipynb)

## How to generate new test collections?
[Test Collections](Data_Test_Generator.ipynb)

## Autocompletion

### Test Collections
We generated two type of autocomplete test collections, [no typo collection](test_autocomplete_pure.json) and [with typo collection](test_autocomplete_typo.json), consists of 200 phrases which are selected from SPARC datasets and NIFS Ontology. The selection process basically is random but still we implement filters to allow more than three character and non numeric phrases only.

### Experiment setup
We compared the performance of AQUA and SPARC's Scigraph to provide autocompletion in the terms of the ability to give suggestion and the execution time.

### Results

<p align="center">
   <img src="https://github.com/napakalas/aqua/blob/aqua_docker/experiment/return_number.png" alt="interface" width="750" height="500"></br>
  <i>Fig 1. Number of completions returned.</i>
</p>

<p align="center">
   <img src="https://github.com/napakalas/aqua/blob/aqua_docker/experiment/execution_time.png" alt="interface" width="750" height="500"></br>
  <i>Fig 2. Execution time analysis.</i>
</p>



## Spelling Error Correction

### Test Collections
We generated test collections by randomly select keywords and authors related to the SPARC data sets, including [biological keywords](query_datatest.json) and [authors](author_datatest.json). Hence, a test collection is a pair of query and a list of the corresponding data sets. Then, we differentiated the test collections based on the number of terms in the query and mimicking typos by performing insertion, deletion, replacement, and spaces removal. In total,there are 31 test collections consisting 50 pairs of query and a list of datasets.

Here are the type of the test collection number of terms in the query we used:
- biological keyword
    - 1 term query
        - no typo
        - 1 deletion
        - 1 insertion
        - 1 replacement
    - 2 terms query
        - no typo
        - 1 deletion
        - 1 insertion
        - 1 replacement
        - no space
        - no space with 1 typo
        - no space with 2 typos
        - no space with 3 typos
        - 3 typos
    - 3 terms query
        - no typo
        - 1 deletion
        - 1 insertion
        - 1 replacement
        - no space
        - no space with 1 typo
        - no space with 2 typos
        - no space with 3 typos
        - 3 typos
- author
    - 1 term query
        - no typo
        - 1 deletion
        - 1 insertion
        - 1 replacement
    - 2 terms query
        - no typo
        - 1 deletion
        - 1 insertion
        - 1 replacement
        - no space

### Experiment setup
We found that the current SPARC search tool uses only the basic features of elastic search, so it can only handle single or double range misspellings such as 'rattusnorvegicus' and 'neromodulation'. Elasticsearch (ES) has a fuzzy search feature that analyzes the query and automatically executes the corrected query. However, its performance seems to suffer when the query consists of more than one typo type. Therefore, AQUA adds an extra layer by creating a pipeline using SymSpellPy to suggest updated queries and then sends them to the ES fuzzy search.

Henceforth, we do not compare the performance of AQUA against the current SPARC search tool because it will be unfair. The tool cannot deal with a complex misspelling such as 'rattunorvegicus', which has one deletion and one replacement, or 'alan garny' and 'marcello bosta', which are the full name of the author. Instead, we compared AQUA to ES and observed whether additional layers could improve retrieval performance. For each test collection, we calculate the Mean Average Precision (MAP) to AQUA and ES. MAP is a standard measurement to represent the performance of an information retrieval system with a range value between 0 and 1, where 1 indicates a perfect system.

### Results


Table 1, shows the Mean Average Precision (MAP) of AQUA and Elasticsearch (ES) over 22 test collections consist of biological keywords as queries. AQUA improves retrieval performance when misspelling occurs in multiple terms queries or misspelling occurs several times. ES is excellent at fixing queries that lose space, even if there is one other typo. However, if there is more than one typo, the performance will be deficient, with a MAP value of 0.057 and 0.18 if there are two or three additional typos in a row. This performance is quite far behind AQUA, which has a MAP of 0.48 and 0.45 for the same addition. Interestingly, AQUA's performance for queries without typos is generally better than ES's.

| Typo            |   ('1 term', 'AQUA') |   ('1 term', 'ES') |   ('2 terms', 'AQUA') |   ('2 terms', 'ES') |   ('3 terms', 'AQUA') |   ('3 terms', 'ES') |
|:----------------|---------------------:|-------------------:|----------------------:|--------------------:|----------------------:|--------------------:|
| 0 typo          |             **0.714785** |           0.711452 |              **0.569673** |           **0.569673**  |              **0.680431** |          0.677097   |
| 1 del           |             0.635935 |           **0.677184** |              **0.555371** |           0.505849  |              **0.668609** |          0.653644   |
| 1 insert        |             0.704785 |           **0.742356** |              0.56559  |           **0.572663**  |              **0.680431** |          0.661312   |
| 1 replace       |             0.644126 |           **0.772202** |              0.548968 |           **0.568364**  |              **0.680431** |          0.646185   |
| no space        |           nan        |         nan        |              0.568006 |           **0.987667**  |              0.667097 |          **0.816667**  |
| no space 1 typo |           nan        |         nan        |              0.559696 |           **0.995918**  |              **0.670508** |          0.0561224  |
| no space 2 typo |           nan        |         nan        |              **0.484005** |           0.0566667 |              **0.644305** |          0.0102041  |
| no space 3 typo |           nan        |         nan        |              **0.446296** |           0.184211  |              **0.589903** |          0.00347222 |
| 3 typo          |           nan        |         nan        |              **0.540761** |           0.481212  |              **0.646919** |          0.621238   |


Table 2 shows the Mean Average Precision (MAP) of AQUA and Elasticsearch (ES) over 9 test collections consisting of authors as queries. AQUA is better to correct a misspelling that appears in a two-term author query. A striking performance difference is AQUA's ability to fix author as a query that loses space where AQUA's MAP is 0.92 while ES is only 0.12.

| Typo      |   ('1 term', 'AQUA') |   ('1 term', 'ES') |   ('2 terms', 'AQUA') |   ('2 terms', 'ES') |
|:----------|---------------------:|-------------------:|----------------------:|--------------------:|
| 0 typo    |             0.863212 |           **0.897673** |              0.926911 |            **0.952778** |
| 1 del     |             0.613025 |           **0.675974** |              **0.818579** |            0.797889 |
| 1 insert  |             0.843871 |           **0.914193** |              0.926944 |            **0.96**     |
| 1 replace |             0.822374 |           **0.867786** |              **0.913039** |            **0.913265** |
| no space  |           nan        |         nan        |              **0.926911** |            0.1245   |
# AQUA - Docker
# About
Docker module for AQUA SPARC

# Docker deployment:
## Initial setup:
  - git clone the project
  - go to the containing directory (the one that having docker-compose.yml):
    - `cd aqua_docker`
  - create `.env` file, set the value:

        SANIC_LOGO="OVERRIDE LOGO USING CONFIG"
        ES_API_KEY = "api-key to acces SciCrunch"
        NM_EMAIL_USR = "email address for notifyme"
        NM_EMAIL_PWD = "email password for notifyme"
        DB_NAME = "sqlite3 file path"


## Deployment using docker-machine
  - Preparation
    - Clone the project and go to the containing directory (the one that having docker-compose.yml):
      - `cd aqua_docker`
    - Make sure docker-machine is installed:
      - `docker-machine version`
    - If not, refer to [here](https://docs.docker.com/machine/install-machine/) to install docker-machine.
    - Install a virtual machine application to create a docker machine. We use [VirtualBox](https://www.virtualbox.org/wiki/Downloads) to deploy locally.
  - Create a new docker machine:
      - `docker-machine create --driver virtualbox --virtualbox-memory 2048 aqua;`
        - creating a docker machine on VirtualBox
        - the minimum VirtualBox allocation is 2048 MB
  - Point docker client to `aqua`:
    - Mac OSX & Linux:
      - `eval "$(docker-machine env aqua)"`
    - Windows:
      - Show the environment of `aqua`
        - `docker-machine env --shell cmd aqua`
      - Set the environment by running a command under
        - `Run this command to configure your shell: `
          - If you use CMD, the command can be:
            - `@FOR /f "tokens=*" %i IN ('docker-machine env --shell cmd aqua') DO @%i`
          - If you use Powershell:
            - `docker-machine env --shell powershell aqua | Invoke-Expression`
  - Build and start the containers:
    - `docker-compose up -d --build`
  - Get the machine IP
    - `docker-machine ip aqua`
  - Now you can access AQUA via web browser with given docker api
    - http://DOCKER-MACHINE-IP/
  - Rerun after computer restart:
    - `docker-machine start aqua`

## Cloud deployment using docker context
  - Create a virtual machine instance in the cloud with minimum RAM 2048 MB
    (for our purpose we create a Ubuntu 18.04 LTS (Bionic) amd64 VM in [Nectar](https://ardc.edu.au/services/nectar-research-cloud/) cloud services)
    - before creating the VM, make sure you have created RSA key pair,
      - setup the public key in the cloud and the VM
      - here are links how to create RSA key pair in different OS:
        - Windows 10: https://phoenixnap.com/kb/generate-ssh-key-windows-10
        - Mac OSX: https://www.siteground.com/kb/how_to_generate_an_ssh_key_pair_in_mac_os/
        - Ubuntu: https://help.ubuntu.com/community/SSH/OpenSSH/Keys
  - Install Docker to the VM using ssh:
    - Connect to the VM:
      - `ssh ubuntu@VM-PUBLIC-IP` (if you are not use ubuntu, change the username)
    - Install Docker:
      - `sudo apt-get remove docker docker-engine docker.io containerd runc`
      - `sudo apt-get update`
      - `sudo apt-get install docker-ce docker-ce-cli containerd.io`
      - `sudo systemctl start docker`
    - Make sure you can run docker
      - `docker run hello-world`
    - If cannot run docker
      - `sudo groupadd docker`
      - `sudo usermod -aG docker $USER`
      - `newgrp docker`
    - Logout from the VM:
      - `exit`

  - Create context
    - `docker context create --docker host=ssh://ubuntu@VM-PUBLIC-IP remote_aqua`
  - Set remote_aqua as default context
    - `docker context use remote_aqua`
  - Build and run the containers in the VM:
    - `docker-compose --context remote_aqua up -d --build`
  - Now you can access AQUA via web browser with your VM public IP
    - http://VM-PUBLIC-IP/


## Cloud rebuild
  - Set remote_aqua as default context
    - `docker context use remote_aqua`
  - Make sure the available container is down
    - `docker-compose --context remote_aqua down`
  - Rebuild, create, and start container
    - `docker-compose --context remote_aqua up --build`

MIT - **Free Software, Enjoy!**

[//]: #URLs
   [sanic]: <https://github.com/channelcat/sanic>
   [nginx]: <https://www.nginx.com/resources/wiki/>
# AQUA - API

# About

API module for AQUA SPARC

## Searching for datasets
GET datasets based on the provided query
#### Format:

    GET /search?query={query}&force={yes|no}&match={yes|no}

-   query (optional): is a query terms.
-   force (optional): to force the server to execute a refined query or not.
    -   yes : execute the refined query.
    -   no (default): execute the original query.
-   match : to signal the server to identify datasets with full query phrase or not.
    -   yes : return datasets having the query phrase.
    -   no (default): return dataset having at least one term of the query phrase.

#### Examples:

-   A GET request without additional arguments will return all datasets
    -   `http://130.216.216.55/search`
-   Example of query: `rat`
    -   `http://130.216.216.55/search?query=rat`
        -   since `rat` is available, query `rat` is executed
-   Example of query: `rattis`
    -   `http://130.216.216.55/search?query=ratis`
        -   since `rattis` is not available, it is refined to `rattus`
        -   `rattus` is executed
    -   `http://130.216.216.55/search?query=ratis&force=yes`
        -   `rattis` is executed since `force` is set to `yes`
        -   no datasets returned since `rattis` is not available in all datasets
-   Example of query: `rattisnorvegicu`
    -   `http://130.216.216.55/search?query=rattis%20norvegicu&match=yes`
        -   `rattisnovergicu` is segmented to `rattis norvegicu`
        -   `rattis norvegicu` is refined to `rattus norvegicus`
        -   `rattus norvegicus` is executed
        -   since `match` is `yes`, only datasets having phrase `rattus norvegicus` are returned.

#### Return example:

    {
      "query": "rattus",
      "executed": "rattus",
      "force": "no",
      "match": "no",
      "suggestions": ["rattle", "rattus sp.", "rousettus"],
      "total": 27,
      "filters": {
        "keywords": {
          "vagus nerve stimulation": ["60", "16", "9"],
          },
        "authors": {
          "Terry Powley": ["10", "12", "9", "123", "107", "90", "11", "24", "121"],
          ....
          }
        },
        "status": {
          "public": ["10", "12", "9", ...],
          "embargoed": ["137"]
        }
      "sorts": {
        "ranking": ["29", "60", "20", "16", "21", "10", "12", "9", "139", "123", "77", "51", "107", "88", "90", "130", "37", "31", "11", "106", "48", "24", "121", "52", "151", "62", "105"],
        "date": ["151", "139", "105", "130", "123", "121", "107", "106", "90", "60", "88", "77", "37", "16", "51", "29", "12", "11", "10", "62", "20", "21", "52", "24", "31", "9", "48"],
        "name": ["51", "31", "48", "77", "52", "130", "123", "105", "62", "106", "9", "37", "107", "121", "29", "90", "151", "20", "21", "139", "16", "60", "88", "10", "11", "12", "24"]
      },
      "hits": {
        "29": {
          "url": "https:\/\/sparc.science\/datasets\/29",
          "banner": "https:\/\/assets.discover.pennsieve.io\/dataset-assets\/29\/6\/revisions\/1\/banner.jpg",
          "_id": "DOI:10.26275\/xmyx-rnm9",
          "_score": 1.1826954,
          "firstPublishedAt": "2019-07-19T23:20:23.12942Z",
          "updatedAt": "2021-06-23T07:25:53.105281Z",
          "name": "Molecular Phenotype Distribution of Single Rat ICN Neurons",
          "description": "We developed an approach to appreciating ... .",
          "readme": {
            "description": "**Study purpose:** The purpose of this study ... ."
          }
          "samples": "151",
          "subjects": "1",
          "anatomy": ["heart"],
          "organisms": ["Rattus norvegicus"],
          "publication": "true",
          "techniques": ["euthanasia technique", "dissection technique", "tissue fixation technique", "staining technique", "computational technique"],
          "highlight": {
            "name":[...],
            "description":[...],
            "readme":[...]
            }
        },
        "60": { .... },
        ....
      }
    }

## Get query autocomplete

### Autocomplete using fast_autocomplete
GET autocomplete from fast-autocomplete
#### Format:

    GET /autocomplete?query=query&limit={limit}

-   query (mandatory): is a query terms
-   limit (optional): the number of returns (default = 10)

#### Example:

-   Example of query: `ratt`
    -   `http://130.216.216.55/autocomplete?query=ratt&limit=10`
        -   get the autocomplete of ratt as a list
        -   return 10 phrases
    -   `http://130.216.216.55/autocomplete?query=ratt`
        -   get the autocomplete of ratt as a list

#### Return example:
    ["rattus","rattlesnake","rattle","rattle virus","ratti","rattus rattus","rattus norvegicus","rattay","rattails","rattus sp"]

### Autocomplete using SciGraph
GET autocomplete from SciGraph
#### Format:

    GET /autocomplete_sc?query=query&limit={limit}

-   query (mandatory): is a query terms
-   limit (optional): the number of returns (default = 10)

#### Example:

-   Example of query: `ratt`
    -   `http://130.216.216.55/autocomplete_sc?query=ratt&limit=10`
        -   get the autocomplete of ratt as a list
        -   return 10 phrases
    -   `http://130.216.216.55/autocomplete_sc?query=ratt`
        -   get the autocomplete of ratt as a list

#### Return example:
    ["rattails","rattle","rattus","rattus leucopus","rattus muelleri","rattus norvegicus","rattus norvegicus genome","rattus norvegicus8","rattus norwegicus"]

## Get query suggestions
GET autocorrection from symspellpy pipeline and suggestions from SciGraph
#### Format:

    GET /suggestions?query=query&limit={limit}

-   query (mandatory): is a query terms
-   limit (optional): the number of returns (default = 10)

#### Example:

-   Example of query: `rattus`
    -   `http://130.216.216.55/suggestions?query=rattus&limit=10`
        -   get the suggestion of rattus as a list
        -   return 10 phrases
    -   `http://130.216.216.55/suggestions?query=rattus`
        -   get the suggestion of rattus as a list

#### Return example:
    ["rattus","rattle","rattus sp.","rousettus"]

## NotifyMe
POST email and keywords to get notification for new datasets
#### Format:

    curl -d “email=email_address&keywords=keywords”

- email: the registrating email
- keywords: the topic keywords to match to new datasets

#### Return:
      - {'success':'true|false'}
