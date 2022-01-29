# cookiecutter-cdp-deployment

[![Cookiecutter Check Status](https://github.com/CouncilDataProject/cookiecutter-cdp-deployment/workflows/Build%20Example%20Repo/badge.svg)](https://github.com/CouncilDataProject/cookiecutter-cdp-deployment/tree/example-build)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03904/status.svg)](https://doi.org/10.21105/joss.03904)

Cookiecutter template for creating new Council Data Project deployments.

---

## Council Data Project

Council Data Project is an open-source project dedicated to providing journalists, activists, researchers, and all members of each community we serve with the tools they need to stay informed and hold their Council Members accountable.

For more information about Council Data Project, please visit [our website](https://councildataproject.org/).

## About

This repository is "cookiecutter template" for an entirely new Council Data Project (CDP) Instance. By following the steps defined in
the [Usage](#usage) section, our tools will create and manage all the database, file storage, and processing infrastructure needed to serve the CDP web application.

While our tools will setup and manage all processing and storage infrastructure, you (or your team) must provide and maintain the custom Python code to gather event information and handle billing for the costs of the deployment.

For more information about costs and billing, see [Cost](#cost).

### CDP Instance Features

-   Plain text search of past events and meeting items<br>
    _(search for "missing middle housing" or "bike lanes")_
-   Filter and sort event and meeting item search results<br>
    _(filter by date range, committee, etc.)_
-   Automatic timestamped-transcript generation<br>
    _(jump right to a specific public comment or debate)_
-   Meeting item and amendment tracking<br>
    _(check for amendment passage, upcoming meetings, etc.)_
-   Share event at timepoint<br>
    _(jump right to the point in the meeting you want to share)_
-   Full event minutes details<br>
    _(view all documents and presentations related to each event)_

See the current [Seattle CDP Instance](https://councildataproject.org/seattle/#/) for a live example.

_Note: Some features are dependent on how much data is provided during event gather. More information see our [ingestion models documentation](https://councildataproject.org/cdp-backend/ingestion_models.html)._

## Usage

Regardless of your deployment strategy, you may find reading the
[Things to Know](#things-to-know) section helpful prior to deployment.

### Deploying Under the councildataproject.org Domain

If you want your deployment under the councildataproject.org domain (i.e. https://councildataproject.org/seattle),
you will need to fill out the ["New Instance Deployment" Issue Form](https://github.com/CouncilDataProject/cookiecutter-cdp-deployment/issues/new/choose).

The Council Data Project team will help you along in the process on the issue from there.

### Deploying Under Your Own Domain

If you want to host your deployment under a different domain (i.e. Your-Org-Name.github.io/your-municipality),
you will need to install `cookiecutter` and use this template.

[**Follow along with the video walkthrough**](https://youtu.be/xdRhh-ocSfc)

Before you begin, please note that you will need to install or have available the following:

-   [gcloud](https://cloud.google.com/sdk/docs/install)
-   [pulumi](https://www.pulumi.com/docs/get-started/install/)
-   [gsutil](https://cloud.google.com/storage/docs/gsutil_install)
-   [Python 3.6+](https://www.python.org/downloads/) (Any Python version greater than or equal to 3.6)

Once all tools are installed, the rest of the infrastructure setup process
should take about 15 minutes.

In a terminal with Python 3.6+ installed:

```bash
pip install cookiecutter
cookiecutter gh:CouncilDataProject/cookiecutter-cdp-deployment
```

Follow the prompts in your terminal and fill in the details for the instance deployment. At the end of the process a new directory will have been created with all required files and further instructions to set up your new deployment.

For more details and examples on each parameter of this cookiecutter template, see [Cookiecutter Parameters](#cookiecutter-parameters)

Follow the next steps in the generated repository's "Initial Repository Setup" section of the generated `README.md` file with the `SETUP` directory.

For more details on what is created from using this cookiecutter template, see [Cookiecutter Repo Generation](#cookiecutter-repo-generation)

The short summary of setup tasks remaining are:

-   The creation of a new GitHub repository for the instance.
-   Logging in or creating accounts for Google Cloud and Pulumi.
-   Initialize the basic infrastructure.
-   Assign a billing account to the created Google Cloud project.
-   Generate credentials for the Google Project for use in automated scripts.
-   Attach credentials as secrets to the GitHub repository.
-   Push the cookiecutter generated files to the GitHub repository.
-   Setup web hosting through GitHub Pages.
-   Enable open access for data stored by Google Cloud and Firebase.
-   Write an event scraper for your municipality (it may be useful to
    build off of [cdp-scrapers](https://github.com/CouncilDataProject/cdp-scrapers))

You can also see an example generated repository and the full steps listed
[here](https://github.com/CouncilDataProject/cookiecutter-cdp-deployment/tree/example-build/SETUP).

### Cookiecutter Parameters

| Parameter                      | Description                                                                                                                 | Example 1                                      | Example 2                                         |
| ------------------------------ | --------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------- | ------------------------------------------------- |
| municipality                   | The name of the municipality (town, city, county, etc.) that this CDP Instance will store data for.                         | Seattle                                        | King County                                       |
| governing_body_type            | What type of governing body this instance is for.                                                                           | city council                                   | county council                                    |
| municipality_slug              | The name of the municipality cleaned for use in the web application and parts of repository naming.                         | seattle                                        | king-county                                       |
| python_municipality_slug       | The name of the municipality cleaned for use in specifically Python parts of the application.                               | seattle                                        | king_county                                       |
| infrastructure_slug            | The name of the municipality cleaned for use in specifically application infrastructure. Must be globally unique to GCP.    | cdp-seattle-abasjkqy                           | cdp-king-county-uiqmsbaw                          |
| maintainer_or_org_full_name    | The full name of the primary maintainer or organization that will be managing this instance deployment.                     | Jackson Maxfield Brown                         | Council Data Project                              |
| hosting_github_username_or_org | The GitHub username or organization that will host this instance's repository. (Used in the web application's domain name)  | JacksonMaxfield                                | CouncilDataProject                                |
| hosting_github_repo_name       | A specific name to give to the repository. (Used in the web application's full address)                                     | cdp-seattle                                    | king-county                                       |
| hosting_github_url             | From the provided information, the expected URL of the GitHub repository.                                                   | https://github.com/JacksonMaxfield/cdp-seattle | https://github.com/CouncilDataProject/king-county |
| hosting_web_app_address        | From the provided information, the expected URL of the web application.                                                     | https://jacksonmaxfield.github.io/cdp-seattle  | https://councildataproject.org/king-county        |
| firestore_region               | The desired region to host the firestore instance. ([Firestore docs](https://firebase.google.com/docs/firestore/locations)) | us-west1                                       | europe-central2                                   |

### Things to Know

Much of Council Data Project processing and resource management can be handled for free and purely on GitHub. However we do rely on a select few resources outside of GitHub to manage all services and applications.

The only service that will require a billing account to manage payment for resources used, is [Google Cloud](#google-cloud). Google Cloud will manage all databases, file storage, and (if needed) [speech-to-text](#speech-to-text) for transcription. You can see more about the average monthly cost of running a CDP Instance in [Cost](#cost).

[Pulumi](#pulumi) is a service to manage and track infrastructure deployment state. For those familiar with [Terraform](https://www.terraform.io/), the two are quite similar. Pulumi's purpose is to ensure that we can move from infrastructure upgrade to infrastructure upgrade without breaking anything (and skipping things that don't need to be done).

For more details see [Cookiecutter Repo Generation](#cookiecutter-repo-generation). _After creating the repo, the following steps will have instructions and links specific to your deployment in the generated repository's README._

### Cookiecutter Repo Generation

`Cookiecutter` is a Python package to generate templated projects. This repository is a template for `cookiecutter` to generate a CDP deployment repository which contains following:

-   A directory structure for your project
-   A directory for your web application to build and deploy from
-   A directory for infrastructure management
-   A directory for your Python event gather function and it's requirements
-   Continuous integration
    -   Preconfigured for your web application to fully deploy
    -   Preconfigured to deploy all required CDP infrastructure
    -   Preconfigured to run CDP pipelines using GitHub Actions

To generate a new repository from this template, in a terminal with Python 3.5+ installed, run:

```bash
pip install cookiecutter
cookiecutter gh:CouncilDataProject/cookiecutter-cdp-deployment
```

_Note: This will only create the basic repository. You will still need to setup Google Cloud and Pulumi accounts._

### Google Cloud

All of your deployments data and some data processing will be done using Google Cloud Platform (GCP).

-   Your deployment's provided and generated data (meeting dates, committee names, councilmember details, etc) will live in [Firestore](https://cloud.google.com/firestore).
-   Your deployment's generated files (audio clips, transcripts, etc.) will live in [Filestore](https://cloud.google.com/filestore).
-   When provided a video without closed captions, the audio from the provided video will be processed using [Speech-to-Text](https://cloud.google.com/speech-to-text).

All of these resources will be set up for you using [Pulumi](#pulumi) but you will need to create both Google Cloud and Pulumi accounts. More information on these services and the steps for account creation can be found in the generated repository's README.

### Pulumi

Pulumi allows CDP developers and Instance maintainers to create, deploy, and manage infrastructure on any cloud using familiar programming languages and tools. It additionally, stores and tracks the _state_ of the CDP infrastructure, i.e. how many and which file storage, database, and processing resources are available.

For CDP Instance maintainers, this simply means, the infrastructure management is packaged up as a part of `cdp-backend`, _and_ the infrastructure will never be incompatible with the pipelines as they are versioned together.

Pulumi is free, and generally, you as an instance maintainer should never have to interact with Pulumi other than during the CDP Instance creation and setup process.

## Cost

CDP was created and maintained by a group of people working on it in their free time. We didn't want to pay extreme amounts of money so why should you?

To that end, we try to make CDP as low cost as possible. Many of the current features are entirely free as long as the repo is open source:

Free Resources and Infrastructure:

-   Event Processing (GitHub Actions)
-   Event and Legislation Indexing (GitHub Actions)
-   Web Hosting (GitHub Pages)
-   Infrastructure State Management (Pulumi)

The backend resources and processing are the only real costs and depend on usage. The more users that use your web application, the more the database and file storage cost. The CDP-Seattle monthly averages below are for the most utilized months of its existance so take these as close to upper-bounds.

Billed Resources and Infrastructure:

-   [Cloud Firestore Pricing](https://firebase.google.com/pricing/)
    _CDP-Seattle monthly average: ~$8.00_
-   [Google Storage Pricing](https://cloud.google.com/storage/pricing#price-tables)
    _CDP-Seattle monthly average: ~$3.00_
-   [Google Speech-to-Text Pricing](https://cloud.google.com/speech-to-text/pricing)
    _CDP-Seattle monthly average: ~$22.00_

**Total Average Monthly Cost**: $33.00

### Speech-to-Text

You may not need to use speech-to-text! In the case your municipalicity provides closed caption files in a format we support parsing and cleaning, we can use those files instead of using speech-to-text. When using closed caption files for transcription generation, CDP-Seattle speech-to-text costs dropped to ~$2.00 / month.

To use closed captions files instead of generating a transcript from audio, your event gather function can attach a [`closed_caption_uri`](https://councildataproject.org/cdp-backend/cdp_backend.pipeline.html#cdp_backend.pipeline.ingestion_models.Session) to the `Session` object.

With speech-to-text cost removed, total average monthly cost for CDP-Seattle is ~$13.00.

### Future Processing Features

As we add more features to CDP that require additional processing or resources we will continue to try to minimize their costs wherever possible. Further, if a feature is optional, we will create a flag that maintainers can set to include or exclude the additional processing or resource usage. See [Upgrades and New Features](#upgrades-and-new-features) for more information.

## Upgrades and New Features

In general, all updates and upgrades are handled easily for you through automated systems that run using GitHub Actions on your repository.

Upgrades are delivered either just before processing time (during package installation) or as a part of weekly system checks and auto-deployments.

#### Backend Pipeline Upgrades

Every time CDP pipelines are initiated, event processing will also install the latest non-breaking-changes version of the [cdp-backend](https://github.com/CouncilDataProject/cdp-backend) package.

In this way, every time your pipeline runs to gather events, index events, or index pieces of legislation, you can be sure you are up-to-date.

#### Backend Infrastructure Upgrades

Every week, your repository will automatically check the status of your infrastructure and deploy updates as detected.

#### Frontend Web Application Upgrades

Every week, your repository will automatically build and deploy the CDP web application, during the build process any non-breaking changes and new features will automatically be pulled in.

### Breaking Changes

In the case we release a new version of [cdp-backend](https://github.com/CouncilDataProject/cdp-backend) or [cdp-frontend](https://github.com/CouncilDataProject/cdp-frontend) that includes breaking changes we will also release an upgrade guide that you will be able to find in this repository.

In the case we can fully automate the upgrade, we will include a script to do so that CDP Instance maintainers will simply have to run for their project.

### Notification of Updates and Release Notes

General CDP updates will be posted to our main website: [https://councildataproject.org](https://councildataproject.org).

However, to receive notifications for individual CDP front and back-end application updates and to receive release notes you can "watch" the two primary repositories.

-   [cdp-frontend](https://github.com/CouncilDataProject/cdp-frontend)
-   [cdp-backend](https://github.com/CouncilDataProject/cdp-backend)

See [GitHub's Documentation on Watching individual repositories](https://docs.github.com/en/github/managing-subscriptions-and-notifications-on-github/configuring-notifications#configuring-your-watch-settings-for-an-individual-repository).

## Citation

If you have found CDP software, data, or ideas useful in your own work, please consider citing us:

Brown et al., (2021). Council Data Project: Software for Municipal Data Collection, Analysis, and Publication. Journal of Open Source Software, 6(68), 3904, https://doi.org/10.21105/joss.03904

```bibtex
@article{Brown2021,
  doi = {10.21105/joss.03904},
  url = {https://doi.org/10.21105/joss.03904},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {68},
  pages = {3904},
  author = {Jackson Maxfield Brown and To Huynh and Isaac Na and Brian Ledbetter and Hawk Ticehurst and Sarah Liu and Emily Gilles and Katlyn M. f. Greene and Sung Cho and Shak Ragoler and Nicholas Weber},
  title = {{Council Data Project: Software for Municipal Data Collection, Analysis, and Publication}},
  journal = {Journal of Open Source Software}
}
```

## License

[MIT](./LICENSE)
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
reported by contacting any of the maintainers of this project and
we will attempt to resolve the issues with respect and dignity.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

## Cookiecutter Contribution vs Application Contribution

Please note that this repository is only the cookiecutter and not the entire
CDP tooling and infrastructure ecosystem. This repository ties all of our
tooling together into a single repository that is easy to deploy and maintain.
Contributions to this repository should largely be documentation,
devops, bugfixes, or similar.

If you experience a bug or incorrect documentation while using the cookiecutter
please do send us a Pull Request! If you want to add or fix an auto-deployment
bot, or add more GitHub Actions to the produced repository, all such contributions
welcome and appreciated.

Examples of these types of contributions include:

-   adding more instance admin documentation to the generated repository
-   updating the `cdp-backend` and `cdp-frontend` version pins
-   upgrading or fixing and auto-deployment GitHub Action
-   adding new GitHub Actions to the generated repository

For contributions to the major pipelines and infrastructure that are used by all
CDP deployments, please see:
[cdp-backend](https://github.com/councildataproject/cdp-backend)

For contributions to the web application which is used by all CDP deployments, please
see: [cdp-frontend](https://github.com/councildataproject/cdp-frontend)

For contributions to the existing event scrapers used by some CDP deployments, please
see: [cdp-scrapers](https://github.com/councildataproject/cdp-scrapers)

## Get Started!

Ready to contribute? Here's how to set up `cookiecutter-cdp-deployment` for local development.

1. Fork the `cookiecutter-cdp-deployment` repo on GitHub.

2. Clone your fork locally:

    ```bash
    git clone git@github.com:{your_name_here}/cookiecutter-cdp-deployment.git
    ```

3. Install `cookiecutter`. (It is also recommended to work in a virtualenv or anaconda environment):

    ```bash
    cd cookiecutter-cdp-deployment/
    pip install cookiecutter
    ```

4. Create a branch for local development:

    ```bash
    git checkout -b {your_development_type}/short-description
    ```

    Ex: feature/read-tiff-files or bugfix/handle-file-not-found<br>
    Now you can make your changes locally.

5. When you're done making changes, check that the cookiecutter still generates
   properly:

    ```bash
    cookiecutter . --no-input
    ```

6. Commit your changes and push your branch to GitHub:

    ```bash
    git add .
    git commit -m "Resolves gh-###. Your detailed description of your changes."
    git push origin {your_development_type}/short-description
    ```

7. Submit a pull request through the GitHub website.
---
title: "Council Data Project: Software for Municipal Data Collection, Analysis, and Publication"
tags:
    - Python
    - JavaScript
    - open government
    - open data
    - open infrastructure
    - municipal governance
    - data archival
    - civic technology
    - natural language processing
authors:
    - name: Jackson Maxfield Brown
      orcid: 0000-0003-2564-0373
      affiliation: 1
    - name: To Huynh
      orcid: 0000-0002-9664-3662
      affiliation: 2
    - name: Isaac Na
      orcid: 0000-0002-0182-1615
      affiliation: 3
    - name: Brian Ledbetter
      affiliation: 1
    - name: Hawk Ticehurst
      affiliation: 1
    - name: Sarah Liu
      affiliation: 4
    - name: Emily Gilles
      affiliation: 4
    - name: Katlyn M. F. Greene
      affiliation: 4
    - name: Sung Cho
      affiliation: 4
    - name: Shak Ragoler
      affiliation: 4
    - name: Nicholas Weber
      orcid: 0000-0002-6008-3763
      affiliation: 1

affiliations:
    - name: University of Washington Information School, University of Washington, Seattle
      index: 1
    - name: University of Washington, Seattle
      index: 2
    - name: Washington University, St. Louis
      index: 3
    - name: Independent Contributor
      index: 4

date: 29 October 2021
bibliography: paper.bib
---

# Summary

Cities, counties, and states throughout the USA are bound by law to archive recordings of public meetings. Most local governments comply with these laws by posting documents, audio, or video recordings online. As there is no set standard for municipal data archives however, parsing and processing such data is typically time consuming and highly dependent on each municipality. Council Data Project (CDP) is a set of open-source tools that improve the accessibility of local government data by systematically collecting, transforming, and re-publishing this data to the web. The data re-published by CDP is packaged and presented within a searchable web application that vastly simplifies the process of finding specific information within the archived data. We envision this project being used by a variety of groups including civic technologists hoping to promote government transparency, researchers focused on public policy, natural language processing, machine learning, or information retrieval and discovery, and many others.

# Statement of Need

Comparative research into municipal governance in the USA is often prohibitively difficult due to a broad federal system where states, counties, and cities divide legislative powers differently. This has contributed to the lack of large-scale quantitative studies of municipal government, and impeded necessary research into effective procedural elements of administrative and legislative processes [@Trounstine2009]. Council Data Project enables large-scale quantitative studies by generating standardized municipal governance corpora - including legislative voting records, timestamped transcripts, and full legislative matter attachments (related reports, presentations, amendments, etc.).

## Related Work

Work in extracting and repackaging government data into machine-readable and experiment ready datasets has historically happened in fields with highly structured data, such as meteorology [@Sparks2017] and legal review and monitoring [@courtlistener]. Notably, there has been prior work in extracting and repackaging municipal government data with [Councilmatic](https://github.com/codeforamerica/councilmatic) [@councilmatic]. However, this work largely aims to make municipal data more accessible to a general public, and does not add any specific data processing to expand the research capabilities of the produced dataset. Recent advances in natural language processing have made it possible to conduct large-scale transcript-based studies on the effects of gender, ideology, and seniority in Supreme Court oral argument [@jacobi2017] and the effects that information communication technology has on civic participation [@einstein2021].

## CDP Architecture

Council Data Project consists of three primary tools:

1. [cookiecutter-cdp-deployment](https://github.com/CouncilDataProject/cookiecutter-cdp-deployment): A Python [cookiecutter](https://cookiecutter.readthedocs.io/) [@cookiecutter] template to assist users in fully deploying a new CDP instance. A "CDP Instance" is a unique deployment of CDP software and tools. For example, there is an "instance" of CDP for the ["Seattle City Council"](https://councildataproject.org/seattle/#/) and an instance of CDP for the "King County Council". Each instance is comprised of its own repository, database, file storage bucket, processing pipelines, and web application.

2. [cdp-backend](https://github.com/CouncilDataProject/cdp-backend): A Python package containing CDP's database schema definition, a file format for transcripts generated by speech-to-text algorithms, an infrastructure specification, and processing pipelines. This package currently contains an event gather and processing workflow that will parse event details, generate a transcript for the event using either the provided closed caption file, or using Google Speech-to-Text from the provided event video, and finally, generate and store event metadata (voting records, thumbnails, minutes items, etc.) This package additionally provides a workflow for generating a TF-IDF based event index for weighted term search. The processing workflows and all utilities and schemas are separate from any one CDP instance so that all CDP instances can be easily upgraded whenever there is a new version of `cdp-backend` released.

3. [cdp-frontend](https://github.com/CouncilDataProject/cdp-frontend): A TypeScript and React-based component library and web application. The web application allows for simple data exploration and sharing, and as such, acts as a method to interactively explore the data produced by the backend pipelines. The web application and the component library are separate from any single CDP instance so that all CDP instances can be easily upgraded whenever there is a new version of `cdp-frontend` released.

## Cookiecutter and the Produced Repository

`cookiecutter-cdp-deployment` will generate all necessary files for an entirely new CDP instance as well as additional setup documentation for the user to follow to fully complete the instance deployment process.

Utilizing [GitHub Actions](https://github.com/features/actions) and [GitHub Pages](https://pages.github.com/), data processing and web hosting are entirely free as long as the user sets their instance's GitHub repository visibility to public.

Deploying a CDP instance incurs some small primary costs by using:

1. [Google Speech-to-Text](https://cloud.google.com/speech-to-text/) for transcript generation.
2. [Firebase Cloud Firestore](https://firebase.google.com/docs/firestore/) for event metadata storage and access.
3. [Firebase Storage](https://firebase.google.com/docs/storage) for file storage and access.

![CDP Core Infrastructure and Pipelines. A CDP instance's event gather and processing pipeline can be triggered by providing the GitHub Action a custom datetime range to process, providing a pre-constructed JSON object (useful for special events like debates and such), or the pipeline will automatically run every 6 hours. Once initiated, the event gather pipeline will get the events for the provided or default date range, create a transcript and extra metadata objects (thumbnails, and more), then will finally archive the event to Firebase. A CDP instance's event indexing pipeline can be triggered by running the pipeline manually or will automatically run every two days.Once initiated, the event index pipeline will in parallel, generate and store unigrams, bigrams, and trigrams from all event transcripts and store them to Firebase for query. Finally, the web application is built from a manual trigger or once a week and simply runs a standard NPM build process then publishes the generated web application to GitHub Pages. Once built and published, the deployment website fetches data from Firebase for search and web access.\label{fig:core-infra}](./assets/cdp_core_infrastructure.png)

CDP tools allow for decentralized control over the management and deployment of each CDP instance while producing a standardized open-access dataset for both research and for municipal transparency and accessibility.

## Data Access

### Web

Once data is processed by a CDP instance, it is available through that instance's interactive web application.

![CDP Web Application. Screenshot of a single event's page. Navigation tabs for basic event details such as the minutes items, the entire transcript, and voting information. Additionally both the transcript search and the full transcript have links to jump to a specific sentence in the meeting. This example event page can be found on our Seattle City Council "staging" instance: http://councildataproject.org/seattle-staging/#/events/0ec08c565d45 \label{fig:event-page}](./assets/event-page-screenshot.png)

### Python

For users who want programmatic access, each instance's repository README includes a programmatic quickstart guide and our database schema is automatically generated and stored in our `cdp-backend` [documentation](https://councildataproject.org/cdp-backend/database_schema.html).

```python
from cdp_backend.database import models as db_models
from cdp_backend.pipeline.transcript_model import Transcript
import fireo
from gcsfs import GCSFileSystem
from google.auth.credentials import AnonymousCredentials
from google.cloud.firestore import Client

# Connect to the database
fireo.connection(client=Client(
    project="cdp-test-deployment-435b5309",
    credentials=AnonymousCredentials()
))

# Read from the database
five_people = list(db_models.Person.collection.fetch(5))

# Connect to the file store
fs = GCSFileSystem(project="cdp-test-deployment-435b5309", token="anon")

# Read a transcript's details from the database
transcript_model = list(db_models.Transcript.collection.fetch(1))[0]

# Read the transcript directly from the file store
with fs.open(transcript_model.file_ref.get().uri, "r") as open_resource:
    transcript = Transcript.from_json(open_resource.read())

# OR download and store the transcript locally with `get`
fs.get(transcript_model.file_ref.get().uri, "local-transcript.json")
# Then read the transcript from your local machine
with open("local-transcript.json", "r") as open_resource:
    transcript = Transcript.from_json(open_resource.read())
```

# Acknowledgements

We wish to thank the many volunteers that have contributed code, design, conversation, and ideas to the project. We wish to thank DemocracyLab and Open Seattle for helping build a civic technology community. From DemocracyLab, we would specifically like to thank Mark Frischmuth for the continued support and helpful discussions. We wish to thank the University of Washington Information School for support. We wish to thank Code for Science and Society and the Digital Infrastructure Incubator for providing guidance on developing a sustainable open source project.

# References
# CDP - {{ cookiecutter.municipality }}

[![Infrastructure Deployment Status]({{ cookiecutter.hosting_github_url }}/workflows/Infrastructure/badge.svg)]({{ cookiecutter.hosting_github_url }}/actions?query=workflow%3A%22Infrastructure%22)
[![Event Processing Pipeline]({{ cookiecutter.hosting_github_url }}/workflows/Event%20Gather/badge.svg)]({{ cookiecutter.hosting_github_url }}/actions?query=workflow%3A%22Event+Gather%22)
[![Event Index Pipeline]({{ cookiecutter.hosting_github_url }}/workflows/Event%20Index/badge.svg)]({{ cookiecutter.hosting_github_url }}/actions?query=workflow%3A%22Event+Index%22)
[![Web Deployment Status]({{ cookiecutter.hosting_github_url }}/workflows/Web%20App/badge.svg)]({{ cookiecutter.hosting_web_app_address }})
[![Repo Build Status]({{ cookiecutter.hosting_github_url }}/workflows/Build%20Main/badge.svg)]({{ cookiecutter.hosting_github_url }}/actions?query=workflow%3A%22Build+Main%22)

---

## Council Data Project

Council Data Project is an open-source project dedicated to providing journalists, activists, researchers, and all members of each community we serve with the tools they need to stay informed and hold their Council Members accountable.

For more information about Council Data Project, please visit [our website](https://councildataproject.org/).

## Instance Information

This repo serves the municipality: **{{ cookiecutter.municipality }}**

### Python Access

Install:

`pip install cdp-backend`

Quickstart:

```python
from cdp_backend.database import models as db_models
from cdp_backend.pipeline.transcript_model import Transcript
import fireo
from gcsfs import GCSFileSystem
from google.auth.credentials import AnonymousCredentials
from google.cloud.firestore import Client

# Connect to the database
fireo.connection(client=Client(
    project="{{ cookiecutter.infrastructure_slug }}",
    credentials=AnonymousCredentials()
))

# Read from the database
five_people = list(db_models.Person.collection.fetch(5))

# Connect to the file store
fs = GCSFileSystem(project="{{ cookiecutter.infrastructure_slug }}", token="anon")

# Read a transcript's details from the database
transcript_model = list(db_models.Transcript.collection.fetch(1))[0]

# Read the transcript directly from the file store
with fs.open(transcript_model.file_ref.get().uri, "r") as open_resource:
    transcript = Transcript.from_json(open_resource.read())

# OR download and store the transcript locally with `get`
fs.get(transcript_model.file_ref.get().uri, "local-transcript.json")
# Then read the transcript from your local machine
with open("local-transcript.json", "r") as open_resource:
    transcript = Transcript.from_json(open_resource.read())
```

-   See the [CDP Database Schema](https://councildataproject.org/cdp-backend/database_schema.html)
    for a Council Data Project database schema diagram.
-   See the [FireO documentation](https://octabyte.io/FireO/)
    to learn how to construct queries using CDP database models.
-   See the [GCSFS documentation](https://gcsfs.readthedocs.io/en/latest/index.html)
    to learn how to retrieve files from the file store.

## Contributing

If you wish to contribute to CDP please note that the best method to do so is to contribute to the upstream libraries that compose the CDP Instances themselves. These are detailed below.

-   [cdp-backend](https://github.com/CouncilDataProject/cdp-backend): Contains all the database models, data processing pipelines, and infrastructure-as-code for CDP deployments. Contributions here will be available to all CDP Instances. Entirely written in Python.
-   [cdp-frontend](https://github.com/CouncilDataProject/cdp-frontend): Contains all of the components used by the web apps to be hosted on GitHub Pages. Contributions here will be available to all CDP Instances. Entirely written in TypeScript and React.
-   [cookiecutter-cdp-deployment](https://github.com/CouncilDataProject/cookiecutter-cdp-deployment): The repo used to generate new CDP Instance deployments. Like this repo!
-   [councildataproject.org](https://github.com/CouncilDataProject/councildataproject.github.io): Our landing page! Contributions here should largely be text changes and admin updates.

## Instance Admin Documentation

You can find documentation on how to customize, update, and maintain this CDP instance
in the
[admin-docs directory]({{ cookiecutter.hosting_github_url }}/tree/main/admin-docs).

## License

CDP software is licensed under a [MIT License](./LICENSE).

Content produced by this instance is available under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).
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
reported by contacting any of the maintainers of this project and
we will attempt to resolve the issues with respect and dignity.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# CDP Instance Setup

This document outlines the steps necessary to finish initializing this CDP Instance.

## Before You Begin

Install the command line tools that will help shorten the setup process

1. Install [gcloud](https://cloud.google.com/sdk/docs/install)
2. Install [pulumi](https://www.pulumi.com/docs/get-started/install/)
3. Install [gsutil](https://cloud.google.com/storage/docs/gsutil_install)

## Initial Repository Setup

There are additional tasks required after generating this repository.

1.  Create the GitHub repository for this deployment to live in.

    [Create a new Repository](https://github.com/new) with the following parameters:

    -   Set the repo name to: **{{ cookiecutter.hosting_github_repo_name }}**
    -   Set the repo owner to: **{{ cookiecutter.hosting_github_username_or_org }}**
    -   Set the repo visibility to: "Public"
    -   Do not initialize with any of the extra options
    -   Click "Create repository".

1.  Login to both Google Cloud and Pulumi.

    During this process Pulumi will provide a token to use for authentication.
    Keep this token available for use in a later step.

    This step should be run while within the `SETUP` directory (`cd SETUP`).

    Run:

    ```bash
    make login
    ```

1.  Initialize the basic project infrastructure.

    This step should be run while within the `SETUP` directory (`cd SETUP`)

    Run:

    ```bash
    make init
    ```

1.  Create (or re-use) a
    [Google Cloud billing account](https://console.cloud.google.com/billing/linkedaccount?project={{ cookiecutter.infrastructure_slug }})
    and attach it to the newly created project ({{ cookiecutter.infrastructure_slug }}).

    For more details on the cost of maintaining a CDP Instance, see our [estimated cost breakdown](https://github.com/CouncilDataProject/cookiecutter-cdp-deployment#cost).

1.  Generate a Google Service Account JSON Key for your Google Cloud Project.

    This will create a directory called `.keys` within this `SETUP` directory and
    add a file called `{{ cookiecutter.infrastructure_slug }}.json` to it
    (i.e. `.keys/{{ cookiecutter.infrastructure_slug }})`. This file will be used later on.

    Run:

    ```bash
    make gen-key
    ```

1.  Attach the Pulumi Access Token and the
    Google Service Account JSON as GitHub Repository Secrets.

    1. Pulumi Access Token -- Create a [new secret]({{ cookiecutter.hosting_github_url }}/settings/secrets/actions/new)

    -   Set the name to: **PULUMI_ACCESS_TOKEN**
    -   Set the value to: The token you kept from step #2
    -   Click "Add secret"

    2. Google Service Account JSON -- Create a [new secret]({{ cookiecutter.hosting_github_url }}/settings/secrets/actions/new)

    -   Set the name to: **GOOGLE_CREDENTIALS**
    -   Set the value to: the contents of the file `.keys/{{ cookiecutter.infrastructure_slug }}.json`
    -   Click "Add secret"

1.  Initialize and push the local repository to GitHub.

    This step should be run while within the base directory of the repository (`cd ..`).

    To initialize the repo locally, run:

    ```bash
    git init
    git add -A
    git commit -m "Initial commit"
    git branch -M main
    ```

    To setup a connection to our GitHub repo, run either:

    ```bash
    git remote add origin {{ cookiecutter.hosting_github_url }}.git
    ```

    Or (with SSH):

    ```bash
    git remote add origin git@github.com:{{ cookiecutter.hosting_github_username_or_org }}/{{ cookiecutter.hosting_github_repo_name }}.git
    ```

    Finally, to push this repo to GitHub, run:

    ```bash
     git push -u origin main
    ```

    Now refresh your repository's dashboard to ensure that all files were pushed.

1.  Once the
    ["Web App" GitHub Action Successfully Complete]({{ cookiecutter.hosting_github_url }}/actions?query=workflow%3A%22Web+App%22)
    configure GitHub Pages.

    Go to your repository's [GitHub Pages Configuration]({{ cookiecutter.hosting_github_url }}/settings/pages)

    -   Set the source to: "gh-pages"
    -   Set the folder to: `/ (root)`
    -   Click "Save"

1.  Once the
    ["Infrastructure" GitHub Action Successfully Completes]({{ cookiecutter.hosting_github_url }}/actions?query=workflow%3A%22Infrastructure%22)
    set the CORS policy for your Storage Bucket.

    This step should be run while within the `SETUP` directory (`cd SETUP`)

    Run:

    ```bash
    make set-cors
    ```

1.  Once the
    ["Infrastructure" GitHub Action Successfully Completes]({{ cookiecutter.hosting_github_url }}/actions?query=workflow%3A%22Infrastructure%22)
    enable data-logging for the Google Speech-to-Text service.

    [Direct Link to Enable](https://console.cloud.google.com/apis/api/speech.googleapis.com/data_logging?project={{ cookiecutter.infrastructure_slug }})

    If the above direct link doesn't work, follow the instructions from
    [Google Documentation](https://cloud.google.com/speech-to-text/docs/enable-data-logging).

1.  Once the
    ["Infrastructure" GitHub Action Successfully Completes]({{ cookiecutter.hosting_github_url }}/actions?query=workflow%3A%22Infrastructure%22)
    configure Firebase Security Rules.

    -   Navigate to [Firebase Console](https://console.firebase.google.com),
        login to the Google Account you used during step #2, select the `{{ cookiecutter.infrastructure_slug }}` Firebase project
        -   Navigate to "Firestore Database", select the "Rules" tab, paste the following in:
            ```
            rules_version = '2';
            service cloud.firestore {
                match /databases/{database}/documents {
                    match /{document=**} {
                        allow read;
                    }
                }
            }
            ```
        -   Click "Publish"
        -   Navigate to "Storage", select the "Rules" tab, paste the following in:
            ```
            rules_version = '2';
            service firebase.storage {
                match /b/{bucket}/o {
                    match /{allPaths=**} {
                        allow read;
                    }
                }
            }
            ```
        -   Click "Publish"

**If all steps complete successful your web application will be viewable at: {{ cookiecutter.hosting_web_app_address }}**

## Data Gathering Setup

Once your repository, infrastructure, and web application have been set up, you will need to write an event data gathering function.

Navigate and follow the instructions in the the file: `python/cdp_{{ cookiecutter.python_municipality_slug }}_backend/scraper.py`.

As soon as you push your updates to your event gather function (`get_events`) to your GitHub repository, everything will be tested and configured for the next pipeline run.

There are some optional configurations for the data gathering pipeline which can be added to `python/event-gather-config.json`. No action is needed for a barebones pipeline run, but the optional parameters can be checked in the [CDP pipeline config documentation](https://councildataproject.org/cdp-backend/cdp_backend.pipeline.html#module-cdp_backend.pipeline.pipeline_config). Note that `google_credentials_file` and `get_events_function_path` should not be modified and will populate automatically if you have followed the steps above.

Be sure to review the [CDP Ingestion Model documentation](https://councildataproject.github.io/cdp-backend/ingestion_models.html) for the object definition to return from your `get_events` function.

Once your function is complete and pushed to the `main` branch, feel free to delete this setup directory.

## Other Documentation

For more documentation on adding data to your new CDP instance and maintainer or customizing your instance
please see the "admin-docs" directory.
# Manually Triggering Event Gather

### Background

While event gather and processing is automated and runs on a schedule,
you may want to manually trigger the event gather and processing pipeline:

1. To backfill data into the infrastructure
    - I.e. to gather the weeks or months from before the instance
      was created.
2. To reprocess data
    - I.e. a prior automated run failed for some reason and you wish to rerun
      that time period.
    - I.e. new and more data in available from the scraper and as such updated and
      new data can enter into the CDP infrastructure.
3. To add a custom or special event
    - I.e. a press conference, debate, forum, etc.

## Backfilling and Reprocessing

To backfill or rerun the pipeline for a specific datetime range go to the
[Event Gather GitHub Action Page]({{ cookiecutter.hosting_github_url }}/actions/workflows/event-gather-pipeline.yml).

Once there, you can add the begin and end datetimes as parameters to the workflow run.

![screenshot of "Run workflow" for event gather pipeline](./resources/backfill-event-gather.png)

See the Python
[`datetime.fromisoformat` documentation](https://docs.python.org/3/library/datetime.html#datetime.datetime.fromisoformat)
for examples of the allowed string patterns for these two parameters.

**Note:** If you chose a datetime range that results in many of events, the pipeline may
error as the pipeline can only run for six hours for one job. So if you wish to
backfill many months of data, consider breaking up the whole time period into
small enough datetime ranges that no single pipeline run lasts longer than six hours.

## Adding Custom Events

To add custom or special events to the CDP infrastructure go to the
[Process Special Event GitHub Action Page]({{ cookiecutter.hosting_github_url }}/actions/workflows/process-special-event.yml).

Once there, you can add the full event details as a parameter to the workflow run.

![screenshot of "Run workflow" for special event pipeline](./resources/special-event-gather.png)

The safest method to construct the event JSON string is using the
`cdp-backend` Python API:

1. `pip install cdp-backend`
2. Create a Python file with your Event definition, at the end,
   print the JSON string of the event
   (including the wrapping `'` characters -- `repr` does this for you).

```python
from datetime import datetime

from cdp_backend.pipeline import ingestion_models

# Define your event
event = ingestion_models.EventIngestionModel(
    body=ingestion_models.Body(name="2021 Mayoral Debates"),
    sessions=[
        ingestion_models.Session(
            video_uri="https://video.seattle.gov/media/council/brief_091321_2012171V.mp4",
            session_datetime=datetime(2021, 9, 13),
            session_index=0,
        ),
    ],
)

# Print out the JSON string in it's full form
print(repr(event.to_json()))
```

3. Run your file

```bash
python your-file.py
```

4. Run the workflow

    Copy and past the output from your terminal into the "Run workflow" parameter and
    then click the "Run workflow" button itself. If the data provided is valid and
    publically available (i.e. videos are publically downloadable), the pipeline will
    run and your custom event will make it's way into CDP infrastructure.
# Customizing Automated Pipelines

Council Data Pipelines all run for you automatically on set schedules
or on push to the repository and generally you should never need to change the
default configuration.

However, you may want to for example, change the interval, or the time that pipelines
run.

To do so, navigate to the `.github` directory of your repository
([GitHub link]({{ cookiecutter.hosting_github_url }}/tree/main/.github)),
and then to the `workflows` sub-directory
([GitHub link]({{ cookiecutter.hosting_github_url }}/tree/main/.github/workflows)).

![image of workflows sub-directory](./resources/workflows.png)

All automated and manual workflows recide in this `.github/workflows` sub-directory.

The following sections will detail what is "safe" to edit.
If there isn't a section for the workflow you are wish to customize,
it is likely because we haven't thought of a reason for needing to customize the
pipeline. In general however, please refer to
[GitHub Actions Documentation](https://docs.github.com/en/actions) for more
information on how all of these workflows are constructed.

## Event Gather and Processing Pipeline

This pipeline ([including the manual trigger](./manual-event-gather.md))
can be found under `event-gather-pipeline.yml`
([GitHub link]({{ cookiecutter.hosting_github_url }}/tree/main/.github/workflows/event-gather-pipeline.yml)).

The common reasons for customizing this pipeline are:

1. to change the automated schedule
2. to add required extra OS level dependencies (such as language and tool installations)

### Customizing Schedule

The pipeline schedule is handled by the
[GitHub Action CRON string](https://docs.github.com/en/actions/reference/events-that-trigger-workflows#scheduled-events).

You may want to change the CRON string to something that more closely matches the
time intervals when you know your municipality tends to post events.

Be careful, don't add too many intervals as if you try to scrape or request
from a resource too often, your pipeline may start to fail because you are
requesting _too_ often.

### Adding Extra OS Dependencies

The pipeline runs on an Ubuntu server and as such to install many OS level
dependencies you can use `apt` to install more packages. In general,
it is safe to add as many of these extra dependencies as you need and there is already
a section where we add dependencies like this to the pipeline.

See the "Install Packages" task of this workflow file and add any more
packages you may need there.

For more programming language support, look into
[GitHub's existing "setup x" actions](https://github.com/actions) and add them
to the pipeline just like the current "actions/setup-python" task.

## Event Indexing Pipeline

This pipeline can be found under `event-index-pipeline.md`
([GitHub link]({{ cookiecutter.hosting_github_url }}/tree/main/.github/workflows/event-index-pipeline.yml)).

The common reasons for customizing this pipeline are:

1. to increase or decrease the frequency a new index is created

### Changing Indexing Frequency

This is controlled by the
[GitHub Action CRON string](https://docs.github.com/en/actions/reference/events-that-trigger-workflows#scheduled-events).

You may want to increase the frequency to generate a fresh index more often or decrease
the frequency because events only happen once a week (or less) and in doing so you can
reduce the cost of running the instance.

The more frequent the index pipeline runs and the more events that it is indexing,
the database cost increases simply due to how many times you are writing many thousands
of documents to the database per pipeline run.
# Adding Analytics

Council Data Project is setup with [Plausible Analytics](https://plausible.io/about)
by default.

## CDP Org Hosted Instances

If your CDP instance is hosted under the "councildataproject.org" domain,
you shouldn't have to change anything from the default cookiecutter settings
as the `data-domain` value for the Plausible Analytics script in `web/public/index.html`
should be set to `councildataproject.github.io` for you.

If it isn't however, please update the `data-domain` value to `councildataproject.github.io`.

Once done, the analytics for your CDP instance should be publicly available on our
[Plausible Dashboard](https://plausible.io/councildataproject.github.io?page=%2F{{ cookiecutter.municipality_slug }}%2F**).

## Self Hosted Instances

If you want to use a different analytics platform or service,
simply replace the Plausible Analytics script in `web/public/index.html` with the service
of your choosing and setup the rest as you normally would.
# Updating to New Cookiecutter Releases

As we have abstracted away the `cdp-backend` for pipeline, infrastructure,
and, database functionality, and `cdp-frontend` for web-app functionality,
and because every time the pipelines which utilize those packages run, the
pipeline pulls in the latest (non-breaking) versions of each, you should
rarely need to update anything yourself.

However, in the instance that a new version of
_this cookiecutter generated repository_ becomes available
(as seen on our [releases page](https://github.com/CouncilDataProject/cookiecutter-cdp-deployment/releases)),
you may want to pull in the changes as they might change default pipeline configuration,
enable new features, add more documentation, and more.

To do so, you should feel comfortable with the command line and git.

### Steps to Upgrade

1. Clone (or fetch and pull) your repository:

    - `git clone {{ cookiecutter.hosting_github_url }}.git`

    OR if you have previously cloned your repository and other changes have occurred:

    - `git checkout main`
    - `git fetch`
    - `git pull main`

2. Install `cookiecutter`:

    You will need Python 3.6+ installed for this.

    - `pip install cookiecutter`

3. _Partially_ run `cookiecutter`:

    We simply want to download the latest version of this cookiecutter template,
    to do so:

    - `cookiecutter gh:CouncilDataProject/cookiecutter-cdp-deployment`

    Re-download the latest version and then when prompted for parameter input,
    escape or exit the process (`control+c`).

    ![screenshot of terminal after cookiecutter interrupt](./resources/cookiecutter-interrupt.png)

4. Update your repo using the latest cookiecutter version:

    This will pull in the current cookiecutter updates as changes to the current
    repository.

    Ensure you are in the repository directory:

    - `cd {{ cookiecutter.municipality_slug }}`

    Then run:

    - `make update-from-cookiecutter`

5. Select the desired changes to commit:

    In an editor or in the terminal you can now review the changes to commit.

    With git in a terminal:

    - `git status`
    - `git add {some-file}`
    - OR `git restore {some-file}`

    ![screenshot of make update and resulting git status](./resources/update-and-git-status.png)

    In VS Code:

    ![screenshot of source control pane in vs code](./resources/vs-code-status.png)

6. Commit and push:

    With the desired changes selected, commit and push to update your repository.
    It is recommended to include the cookiecutter version in your commit message.

    - `git commit -m "Update to cookiecutter v0.1.0"`
    - `git push origin main`

After all steps are complete, you should see updates on your repository and
any pipelines that run on `push` they will automatically trigger.
<!--
  Thank you for submitting a pull request!

  ⚠️⚠️ Please do the following before submitting: ⚠️⚠️

  - Ensure that the code is up-to-date with the `main` branch.
  - Provide or update documentation for any feature added by your pull request.
  - Provide relevant tests for your feature or bug fix.

  ❗️ Also: ❗️

  Please name your pull request {development-type}/{short-description}.
  For example: feature/read-tiff-files
-->

### Link to Relevant Issue

This pull request resolves #

### Description of Changes

_Include a description of the proposed changes._
---
name: Bug Report
about: Create a report to help us improve cookiecutter-cdp-deployment
labels: bug
---

<!--
  ⚠️⚠️ Please do the following before submitting: ⚠️⚠️

  📖 Please read our Code of Conduct.
  🔎 Please search existing issues to avoid creating duplicates.
-->

### Describe the Bug

_A clear and concise description of the bug._

### Expected Behavior

_What did you expect to happen instead?_

### Reproduction

_Steps to reproduce the behavior and/or a minimal example that exhibits the behavior._

### Environment

_Any additional information about your environment._

-   OS Version: _[e.g. macOS 11.3.1]_
-   Cookiecutter Version: _[e.g. 0.5.0]_
---
name: Feature Request
about: Suggest a feature for cookiecutter-cdp-deployment
labels: enhancement
---

<!--
  ⚠️⚠️ Please do the following before submitting: ⚠️⚠️

  📖 Please read our Code of Conduct.
  🔎 Please search existing issues to avoid creating duplicates.
-->

### Feature Description

_A clear and concise description of the feature you're requesting._

### Use Case

_Please provide a use case to help us understand your request in context._

### Solution

_Please describe your ideal solution._

### Alternatives

_Please describe any alternatives you've considered, even if you've dismissed them._
