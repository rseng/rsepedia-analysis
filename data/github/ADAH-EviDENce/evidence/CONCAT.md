# Troubleshooting

## Permisson issue on Windows

When using Windows 10 Pro, the Docker containers may not be able to read or write to the folder with the corpus and model. In this case, you get output like the following:

```powershell
PS C:\Users\JohnDoe\evidence-master\evidence-master> docker-compose up --build
Creating network "evidence-master_default" with the default driver
Building server
Step 1/28 : FROM golang:1.12-stretch as buildserver
 ---> 563601c9e3b2
...
Creating evidence-master_elasticsearch_1 ... done                                                                               Creating evidence-master_indexer_1       ... error
ERROR: for evidence-master_indexer_1  Cannot create container for service indexer: status code not OK but 500: {"Message":"Unhandled exception: Filesharing has been cancelled",...

ERROR: for evidence-master_server_1  Cannot create container for service server: status code not OK but 500: {"Message":"Unhandled exception: Filesharing has been cancelled",...

ERROR: for indexer  Cannot create container for service indexer: status code not OK but 500: {"Message":"Unhandled exception: Filesharing has been cancelled",...

ERROR: for server  Cannot create container for service server: status code not OK but 500: {"Message":"Unhandled exception: Filesharing has been cancelled",...
ERROR: Encountered errors while bringing up the project.
```

This can be solved by explicitly giving Docker access to the folder. You can do this by:

1. opening the Docker dashboard by right-clicking the Docker icon in the Windows taskbar
2. go to Settings/Resources/FILE SHARING and add the folder where you extracted the files before
3. try running the command again:

    ```powershell
    $Env:EXPERIMENT="getuigenverhalen"
    docker-compose up --build
    ```
# evidence

`evidence` --  a doc2Vec-based assisted close reading tool with support for abstract concept-based search and context-based search.

| Five recommendations for fair software from [fair-software.nl](https://fair-software.nl) | Badges |
| --- | --- |
| 1. Code repository | [![GitHub badge](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/ADAH-EviDENce/evidence/) |
| 2. License | [![License badge](https://img.shields.io/github/license/ADAH-EviDENce/evidence)](https://github.com/ADAH-EviDENce/evidence/) |
| 3. Community registry | [![Research Software Directory](https://img.shields.io/badge/rsd-evidence-00a3e3.svg)](https://www.research-software.nl/software/evidence) |
| 4. Enable citation | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3954885.svg)](https://doi.org/10.5281/zenodo.3954885) |
| 5. Checklist | N/A |
| **Other best practices** | |
| Test model generation | [![Test model generation](https://github.com/ADAH-EviDENce/evidence/workflows/Test%20model%20generation/badge.svg)](https://github.com/ADAH-EviDENce/evidence/actions?query=workflow%3A%22Test+model+generation%22) |
| Frontend | [![Frontend](https://github.com/ADAH-EviDENce/evidence/workflows/Frontend/badge.svg)](https://github.com/ADAH-EviDENce/evidence/actions?query=workflow%3A%22Frontend%22) |
| docker-compose | [![docker-compose](https://github.com/ADAH-EviDENce/evidence/workflows/docker-compose/badge.svg)](https://github.com/ADAH-EviDENce/evidence/actions?query=workflow%3Adocker-compose) |
| GitHub Super Linter| [![Lint Code Base](https://github.com/ADAH-EviDENce/evidence/workflows/Lint%20Code%20Base/badge.svg)](https://github.com/ADAH-EviDENce/evidence/actions?query=workflow%3A%22Lint+Code+Base%22) |
| Markdown Link Checker| [![Check Markdown links](https://github.com/ADAH-EviDENce/evidence/workflows/Check%20Markdown%20links/badge.svg)](https://github.com/ADAH-EviDENce/evidence/actions?query=workflow%3A%22Check+Markdown+links%22) |

## Machine-supported research in humanities
While research in the humanities has been able to leverage the digitization of text corpora and the development of computer based text analysis tools to its benefit, the interface current systems provide the user with is incompatible with the proven method of scholarly close reading of texts which is key in many research scenarios pursuing complex research questions.

What this boils down to, is the fact that it is often restrictive and difficult, if not impossible, to formulate adequate selection criteria, in particular for more complex or abstract concepts, in the framework of a keyword based search which is the standard entry point to digitized text collections.

## Querying by example - close reading with tailored suggestions
`evidence` provides an alternative, intuitive entry point into collections by leveraging the doc2vec framework. Using doc2vec `evidence` learns abstract representations of the theme and content of the elements of the user's corpus. Then, instead of trying to translate the scientific query into keywords, after compiling a set of relevant elements as starting points, i.e. examples of the concept the user is interested in, the user can query the corpus based on these examples of their concept of interest. Specifically, `evidence` retrieves elements with similar abstract representations and presents them to the user, using the users feedback to refine its retrieval.
Furthermore, this concept-based query mode is complemented by the ability to perform additional retrieval using `more-like-this` context based retrieval function provided by `elasticsearch`.
Together, this enables a user to combine the power of a close-reading approach with that of a large digitized corpus, selecting elements from the entire corpus which are likely to be of interest, but leaving the decision up to the user as to what evidence they deem useful.


## Documentation for users

## Running the demo

The repository contains a demonstration including a corpus and a model. The demonstration allows you the explore the features of this software without supplying your own corpus.

Prerequisites:

- [docker](https://docs.docker.com/engine/install/)
- [docker-compose](https://docs.docker.com/compose/install/)

### Step 1

First test that the docker installation is working. Depending on your system, you need to use either a PowerShell (on Windows) or a terminal (on Linux or on MacOS).

- For Windows:

    Open a Powershell prompt (press Windows+S and type Powershell) and run:

    ```powershell
    docker run hello-world
    ```

- For Linux/MacOs:

    Open a terminal and run:

    ```shell
    docker run hello-world
    ```

This should show a message that your Docker installation is working correctly. If so, we can proceed to the installation of ``evidence``, otherwise we suggest to check the [Docker troubleshooting page](https://docs.docker.com/docker-for-windows/troubleshoot/).

### Step 2

[Download](https://github.com/ADAH-EviDENce/evidence/archive/master.zip) a copy of evidence archive and extract its contents on your machine.

Alternatively, if you have [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) installed, you can also clone the repository.

```shell
git clone https://github.com/ADAH-EviDENce/evidence.git
```

### Step 3

- For Windows:
  - Open a Powershell prompt
  - Change your current working directory to where you extracted the files. For instance:

      ```powershell
      cd C:\Users\JohnDoe\Downloads\evidence-master\evidence-master
      ```

- Linux/MacOS:
  - Open a terminal
  - Change your current working directory to where you extracted the files. For instance:

      ```shell
      cd /home/JohnDoe/Downloads/evidence
      ```

The demo can be started with the commands below. Keep this PowerShell/Terminal window open and running during the demo.

- Set the experiment name

    For Windows:

    ```powershell
    $Env:EXPERIMENT="getuigenverhalen"
    ```

    For Linux/MacOS:

    ```shell
    export EXPERIMENT="getuigenverhalen"
    ```

- Start the demo

    ```shell
    docker-compose up --build
    ```

The command above downloads necessary Docker images, builds all the Docker images and starts the demo.

The command prints many log messages. If all goes well, the last lines of the output should be:

```shell
...
indexer_1        | Indexing done.
evidence-master_indexer_1 exited with code 0
```

Check [troubleshooting](troubleshooting.md) if you have any issues about this step.

### Step 4

Go to the following URL in your web browser: [http://localhost:8080/](http://localhost:8080/ui/search/).

### Step 5

Once you are done with exploring the demo, you can stop it by selecting the PowerShell/Terminal that is still running the demo and press Ctrl+C.

## Generating a model

### Prerequisites

Verify that your ``docker-compose`` version is at least 1.25.4. (Earlier versions may work).

```shell
docker-compose --version
```

Verify that your ``docker`` version is at least 19.03.12. (Earlier versions may work).

```shell
docker --version
```

If you want to use your own corpus, refer to [./experiments/README.md](./experiments/README.md) for notes on the required format and directory layout.

### Define which corpus to use

Define the name of the dataset/experiment. Here we choose 'getuigenverhalen'. The corpus files should reside under ``/experiments/<EXPERIMENT>/corpus``, see sample corpora.

```shell
export EXPERIMENT=getuigenverhalen
```

### Building the model generation image

Be aware that building can take a couple of minutes.

```shell
# (starting from the repo root directory)
docker-compose --file generate-model.yml build generate-model
```

## Generating the doc2vec model

```shell
# (starting from the repo root directory)
docker-compose --file generate-model.yml run --user $(id -u):$(id -g) generate-model
```

### Build the user interface web application and start it

```shell
# (starting from the repo root directory)
export EXPERIMENT=getuigenverhalen
docker-compose build
docker-compose up
```

Frontend should now be usable at [``http://localhost:8080``](http://localhost:8080).

> We strongly suggest not making the frontend available publicly as there is no authentication. Anyone with the url will have access to the frontend.
Running it on a local network, for example a university network, should be protected from most evil-doers.

Besides interaction with a web browser you can also interact with the frontend from the command line see [here](ui#elastic-search-example-queries) and [here](ui#doc2vec-example-queries) for examples.

## Optional: manage frontend users

The first page of the frontend forces you to select a user or 'gebruiker' in Dutch.
A user called `demo` exists and can be selected.

### Change initial user

The initial user named ``demo`` can be renamed by setting the `FRONTEND_USER` environment variable before running `docker-compose up`.

For example to have `myinitialusername` as a user name, do the following:

```shell
# (starting from the repo root directory)
export EXPERIMENT=getuigenverhalen
export FRONTEND_USER=myinitialusername
docker-compose up
```

### Add additional users

If the existing user is not enough, you can add a new user to the frontend with the following command
(you can choose your own username by replacing `mynewusername` value in the command below):

```shell
export EXPERIMENT=getuigenverhalen
export FRONTEND_USER=mynewusername
docker-compose run usercreator
```

To add more users, repeat the command with different values for `FRONTEND_USER`.

---

## Documentation for developers

### Checking the MarkDown links

When updating the documentation, you can check if the links are all working by running:

```shell
npm install
npm run mlc
```

### Related repositories

[https://github.com/ADAH-EviDENce/EviDENce_doc2vec_docker_framework](https://github.com/ADAH-EviDENce/EviDENce_doc2vec_docker_framework)

[https://github.com/ADAH-EviDENce/evidence-gui](https://github.com/ADAH-EviDENce/evidence-gui)
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
reported by contacting the project team at generalization@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org), version 1.4,
available at [https://www.contributor-covenant.org/version/1/4/code-of-conduct.html](https://www.contributor-covenant.org/version/1/4/code-of-conduct.html).

For answers to common questions about this code of conduct, see
[https://www.contributor-covenant.org/faq](https://www.contributor-covenant.org/faq)
# Preparing a corpus for use with evidence

The aim of the `evidence` tool is to enable and support researchers interested in
applying a close-reading methodology for their scientific question to efficiently
make use of large digitized text corpora. Of course, different researchers will
be interested in different corpora based on their scientific focus. Rather than
provide support for (a selection of) specific corpora, `evidence` is therefore designed
to ingest and make accessible via its user interface (UI) a corpus of the user's choice.

The ability to do so, however, depends on the user presenting `evidence` the corpus in
question in a manner the tool understands. Accordingly, this document outlines what is
understood as a corpus in the context of `evidence`, and in what format the user should
supply their corpus in order to make use of `evidence`. Given the plethora of formats in
which corpora may exist in digital archives we do not provide general instructions on
how to prepare any specific corpus to meet the required format, instead focusing only on
the structure `evidence` expects.

## A corpus

In the context of `evidence`, we refer to a corpus as a collection of individual documents
or texts.

### Documents
There are no strict upper or lower limits to the number of documents a collection
can/should contain, however, as `evidence` trains a machine learning model on the corpus,
too small corpora may give rise to poor performance. Typically, a corpus containing several
million words, respectively several ten thousand paragraphs (see below) should suffice.

### Paragraphs
The expectation is that each document will consist of a (varying) number of paragraphs.
There is no upper bound on the number of paragraphs an individual document may have, nor is there
a strict lower bound. However, as a rule of thumb, at least 2 paragraphs (with of order 150 words) is
a good lower bound.

### Languages
Preferably, all documents within a corpus should be written in one language. However, this is
not a strict requirement. If multiple languages are present, the fractions of the total corpus in
each language should ideally be balanced. If this is not the case, `evidence` will still work,
however results for the minority language(s) may be untrustworthy, especially if the corpus itself
is small ( < 10.000 fragments).

## Structure within `evidence`

While the constituent documents provide a natural subdivision of a corpus, they are not sufficiently
granular for accurately identifying sections of text with a common topic or theme. Therefore, rather
than processing the documents of a corpus as a whole, `evidence` makes use of sub-sections of
documents of a corpus.

### Fragments
Referred to as fragments, these sub-sections represent the atomic unit of a corpus within `evidence` and
the set of all fragments makes up the corpus. Each fragment has a unique id, which simultaneously also links it to its parent document. Together the fragments form a flat hierarchical layer on which search queries are executed.

We strongly suggest a fragment length of ca. 150 words, while not splitting inside of a paragraph. This choice
is meant to restrain a fragment to one main topic, while providing enough length/context to account for more
complex topics. It is, however, up to the user to decide the fragment length for their corpus.


## Corpus preparation

In order to ingest a corpus into the evidence tool, the user must provide it as a nested file/folder structure.

Within a top-level folder named `corpus` the user places a folder per corpus document, with a unique name
identifying the document.

Within each of these folders all fragments constituting a document are placed. Each
fragment should have a name identifying its place in the sequence of fragments constituting the parent document
and adhering to the format `AAA_paragraph_XX-YY_BBB.txt`, where
`AAA` can be any string (not containing `paragraph_`), `XX` and `YY` should be two numbers
with `YY` > `XX` denoting the paragraphs contained in the fragment, and `BBB` can again be any string
(not containing `paragraph_`). Valid examples are:

`paragraph_9-12.txt`

`mycol_paragraph_99-112_clean.txt`

This filename must be unique WITHIN a document, but need not be unique within the corpus.

Each of these fragment files should be formatted as a single line of flat ascii text
containing the full content of the contained paragraphs.

## Corpus ingestion
The aforementioned `corpus` folder with its substructures and contents can then be copied to the folder the
user has created in the `experiment` folder and ingested into the evidence tool as detailed in the instructions
(see [README](../README.md)).
# Evidence user interface

User interface consists of

* a web service written in Go language in current folder.
* a single page web application written in React in `ui/` folder.

## Run

See [/README.md#build-the-user-interface-web-application-and-start-it](/README.md#build-the-user-interface-web-application-and-start-it).

## Elasticsearch example queries

List of snippets:

    http://localhost:8080/es/snippets/_search

Single snippet:

    http://localhost:8080/es/snippets/snippet/$id

More like this:

    curl -XGET -H 'Content-Type: application/json' \
        http://localhost:8080/es/snippets/_search -d '{
        "query": {
            "more_like_this": {
                "fields": ["text", "lemma"],
                "boost_terms": 1,
                "max_query_terms": 150,
                "min_doc_freq": 1,
                "min_term_freq": 1,
                "like": [{
                    "_index": "snippets",
                    "_type": "snippet",
                    "_id": "'${id}'"
                }]
            }
        }
    }'

## Doc2Vec example queries

More like this:

    curl http://localhost:8080/doc2vec/$id

Paging in `doc2vec` is done using request parameters `from` (default `0`) and `size` (default `10`)

    curl http://localhost:8080/doc2vec/$id?from=3&size=8

## Legal matters

Copyright 2016-2019 Koninklijke Nederlandse Academie van Wetenschappen

Distributed under the terms of the GNU General Public License, version 3.
See the file LICENSE for details.
# UI

User interface of Evidence.

## Development

- Start containers as described in ../README.md
- Run in ./ui: `npm start`
- Open: [http://localhost:3000](http://localhost:3000)

### First time

- Run in ./ui: `npm install`

## Build & run

See [/README.md#build-the-user-interface-web-application-and-start-it](/README.md#build-the-user-interface-web-application-and-start-it).

To run the Go web service (one level up) and the React app from different hosts, the React app must be built with environment variable `REACT_APP_HOST` set to the host URL of the Go web service.

## React Starter App

This project was bootstrapped with [Create React App](https://github.com/facebook/create-react-app).
# Instructions
See [/README.md#generating-a-model-from-the-corpus](/README.md#generating-a-model-from-the-corpus).
