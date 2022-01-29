---
title: 'Beiwe: A data collection platform for high-throughput digital phenotyping'
tags:
  - high-throughput
  - digital phenotyping
  - multi-language
  - smartphone
  - AWS
authors:
  - name: Jukka-Pekka Onnela^[Principal Investigator and corresponding author.]
    orcid: 0000-0001-6613-8668
    affiliation: 1
  - name: Caleb Dixon
    affiliation: 2
  - name: Keary Griffin
    affiliation: 3
  - name: Tucker Jaenicke
    affiliation: 2
  - name: Leila Minowada
    affiliation: 2
  - name: Sean Esterkin
    affiliation: 2
  - name: Alvin Siu
    affiliation: 2
  - name: Josh Zagorsky
    affiliation: 2
  - name: Eli Jones
    affiliation: 2
affiliations:
 - name: Department of Biostatistics, Harvard T.H. Chan School of Public Health, Harvard University
   index: 1
 - name: Zagaran, Inc.
   index: 2
 - name: Rocket Farm Studios
   index: 3
date: 21 April 2021
bibliography: paper.bib

---

# Summary

Beiwe is a high-throughput data collection platform for smartphone-based digital phenotyping. It has been in development and use since 2013. Beiwe consists of two native front-end applications: one for Android (written in Java and Kotlin) and one for iOS (written in Swift and Objective-C). The Beiwe back-end, which is based on Amazon Web Services (AWS), has been implemented primarily in Python 3.6, but it also makes use of Django for ORM and Flask for API and web servers. It uses several AWS services, such as S3 for flat file storage, EC2 virtual servers for data processing, Elastic Beanstalk for orchestration, and RDS for PostgreSQL database engine. Most smartphone applications use software development kits (SDKs) that generate unvalidated behavioral summary measures using closed proprietary algorithms. These applications do not meet the high standards of reproducible science, and often require researchers to modify their scientific questions based on what data happens to be available. In contrast, Beiwe collects raw sensor and phone use data, and its data collection parameters can be customized to address specific scientific questions of interest. Collection of raw data also improves reproducibility of studies and enables re-analyses of data and pooling of data across studies. Every aspect of Beiwe data collection is fully customizable, including which sensors to sample, how frequently to sample them, whether to add Gaussian noise to GPS location, whether to use Wi-Fi or cellular data for uploads, how frequently to upload data, specification of surveys and their response options, and skip logic. All study settings are captured in a JSON-formatted configuration file, which can be exported from and imported to Beiwe to enhance transparency and reproducibility of studies.

# Statement of need

**Background.** The phenotype of an organism is a collection of traits, such as enzyme activity, hormone levels, and behavior. Many researchers are increasingly advocating a more substantial role for large scale phenotyping as a route to advances in the biomedical sciences [@houle2010phenomics; @delude2015deep; @bilder2009phenomics; @robinson2012deep]. Of the many different phenotypes, social, behavioral, and cognitive phenotypes are difficult to study due to their temporal and contextual dependence. The standard approach to learning about something as complex as human behavior is to use surveys, which are cross-sectional, subjective, and often burdensome. The ubiquity and capability of smartphones—when coupled with appropriate data analytic techniques—can be part of the solution. We coined the term Digital Phenotyping to refer to the “moment-by-moment quantification of the individual-level human phenotype *in situ* using data from personal digital devices, in particular smartphones” [@onnela2016harnessing; @torous2016new]. We developed Beiwe specifically for use in smartphone-based digital phenotyping research. In addition to enabling more objective measurement of existing phenotypes, the approach can also give rise to entirely new phenotypes.

**State of the field.** Social and behavioral phenotypes have traditionally been studied using either self-administered or investigator-administered surveys in research settings and self-administered or clinician-administered surveys in clinical settings. For example, the Amyotrophic Lateral Sclerosis Functional Rating Scale - Revised (ALSFRS-R) includes 12 questions (items) each scored on 0 (no function) to 4 (full function) scale, and it has been used for both diagnosing patients and for measuring disease progression [@mora2017edaravone]. In observational studies and clinical trials, it may be administered every six weeks, with smaller within-subject standard deviation in the score if self-administered by the patient on the smartphone rather than administered by a clinician at the clinic on paper [@berry2019design]. To eliminate recall bias, some of these items may be actually measured objectively in free-living settings. For example, two of the items in ALSFRS-R are related to physical activity: walking (Item 8) and climbing stairs (Item 9), both of which can be estimated using smartphone accelerometer and gyroscope data [@straczkiewicz2021systematic].

The development of Beiwe was driven by several key considerations: (1) Beiwe is open source software under a permissive 3-clause BSD license; (2) it has been developed and continues to be maintained by professional software engineers to ensure high quality of code base; (3) it supports both Android and iOS devices to allow access to nearly all  smartphones; (4) it collects raw sensor data to allow reproducible research and pooling and re-analyses of data; (5) it is compliant with HIPAA for maintaining the privacy of individually identifiable information; (6) it enables full configurability to accommodate the data collection needs of different studies; (7) it emphasizes study replicability and reproducibility by capturing all data collection settings in a JSON configuration file that can be imported into / exported from Beiwe; (8) it uses scalable cloud-based architecture.

Digital phenotyping is related to sensing in computer science. There is a Wikipedia page that compares different mobile phone sensing platforms [@wiki_sensing], but this page appears to be out of date. The comparison includes 26 pieces of software, but these platforms are mostly not intended for digital phenotyping, which among others requires continuous (or essentially continuous) collection and banking of raw passive data. In addition, for representative cohorts and equitable participation in research, it is important to support both Android and iOS users given the socioeconomic differences across their user base; 14 of these platforms support both. To support reproducibility and replicability, being open source is important, which leaves 8 software systems. Some of these systems have not been updated in years, which leaves us with 4 pieces of software: AWARE, Beiwe, EARS, and mindLamp. As far as we know, AWARE has not been used in biomedical research; EARS is not actually open source as it has significant portions of code redacted [@github_ears], although the Wikipedia site has it listed as open source. Finally, mindLamp appears to have been released fairly recently; although it can collect some passive data, it appears to be geared mostly towards clinical use, e.g., it allows patients to track medication, keep a journal, and perform guided meditations.

An important development for the field has been the introduction of software development kits (SDKs) for smartphones, such as Apple’s ResearchKit and Google’s Research-Stack, which has facilitated writing of software for these devices. Use of prepackaged software however limits what types of data can be collected, which then limits the types of analyses that may be performed [@onnela2021opportunities]. For example, Apple’s ResearchKit does not support background sensor data collection [@researchkit]; Apple's HealthKit supports background sensor data collection for selected sensors only [@sensorkit]; and the Core Motion framework makes it possible to collect raw accelerometer data in the background but only for up to 12 hours at a time [@cmsensor]. The algorithms that underlie HealthKit metrics, such as step count [@sensorkit_stepcount], appear all to be proprietary. The use of closed algorithms, which may change at any time without notice, makes it hard or impossible to compare data collected at different times (using possibly different versions of these algorithms) or data collected using different SDKs.

**Workflow.** The Beiwe platform is used primarily in health research settings to collect active and passive data. When using open-source Beiwe, the first step is to deploy the AWS-based system back-end. While historically it's been possible to use either a single server deployment or a scalable server cluster deployment, currently only the server cluster deployment is documented and supported. It's worth pointing out that the back-end deployment is non-trivial, and among other things it requires creating an AWS account, setting up a user with sufficient permissions, generating credentials for programmatic access, obtaining a domain name, and creating a Sentry account for monitoring errors. Once deployed, researchers can then login to the web-based Beiwe portal where they create a study, which includes specifying surveys and configuring passive data collection settings. Note that a single back-end can support tens or even hundreds of studies each with their own surveys and passive data collection settings. Investigators then create a collection of Beiwe user IDs, which are randomly generated 8-character strings (e.g., abcd1234), and distribute them to study subjects together with their temporary passwords. Subjects then download the Beiwe2 smartphone application from Apple App Store (iOS) or Google Play Store (Android) and enter their Beiwe ID, password, and the name of the study server (e.g., my.beiwe.org). The main goal of Beiwe is to collect high-throughput passive data in the background, and as such the idea is for the application to intervene as little as possible in the daily life of the subject. Therefore, the only time the user would actively use the Beiwe smartphone application is when entering surveys or contributing audio diary entries; at all other times the application is simply running in the background. Note that Beiwe never asks the subject for her name, phone number, or any other identifier; the only piece of information that links a user to the Beiwe account is the Beiwe user ID.

**Data elements.** The Beiwe platform can collect both active data (subject input required) and passive data (subject input not required). Currently supported active data types for both Android and iOS are text surveys and audio diary entries and their associated metadata. Passive data can be further divided into two groups: phone sensor data and phone logs. Beiwe collects raw sensor data and raw phone logs, which is absolutely crucial in scientific settings, yet this point remains underappreciated. Relying on generic SDKs for data collection is convenient but for many reasons ineffective: the data generated by different SDKs are not comparable; the algorithms used to generate data are proprietary and hence unsuitable for reproducible research; new data summaries cannot be implemented post data collection; the composition of the metrics collected by SDKs changes in time, making it difficult or impossible to make comparisons across subjects when they are enrolled at different points in time; and data cannot be pooled across studies for later re-analyses or meta-analyses. Currently supported passive data types are the following: accelerometer, gyroscope, magnetometer, GPS, call and text message logs on Android devices (metadata only, no content), proximity,  device motion, reachability, Wi-Fi, Bluetooth, and power state. For each sensor, such as GPS, data collection alternates between an on-cycle (data collected) and an off-cycle (data not collected); logs are collected without sampling if their collection is specified. The investigators specify what data is collected based on the scientific question: for example, a study on mobility might choose to collect accelerometer and gyroscope data used in human activity recognition. The text that appears within the application is also customizable for each study. Data streams that contain identifiers, such as phone numbers in communication (call and text message) logs, are anonymized on the device; the "fuzzy" GPS feature if enabled adds randomly generated noise to the GPS coordinates on the device. Finally, study meta settings are also customizable, and include items such as frequency of uploading data files to the back-end (typically 1 hour) and duration before auto logout from the application (typically 10 minutes).

**Privacy and security.** All Beiwe data are encrypted while stored on the phone awaiting upload and while in transit, and they are re-encrypted for storage on the study server while at rest. More specifically, during study registration the platform provides the smartphone app with the public half of a 2048-bit RSA encryption key. With this key the device can encrypt data, but only the server, which has the private key, can decrypt it. Thus, the Beiwe application cannot read its own data that it stores temporarily, and therefore there is no way for a user (or anyone else) to export the data. The RSA key is used to encrypt symmetric Advanced Encryption Standard (AES) keys for bulk encryption. These keys are generated as needed by the app and must be decrypted by the study server before data recovery. Data received by the cloud server is re-encrypted with the study master key provided and then stored on the cloud. Some of the collected data contain identifiers: communication logs on Android devices contain phone numbers, and Wi-Fi and Bluetooth scans contain media access control (MAC) address. If the study is configured to collect these data, the identifiers in them are anonymized on the phone, and only anonymized versions of the data are uploaded to the back-end server. Briefly, the Beiwe front-end application generates a unique cryptographic code, called a salt, during the Beiwe registration process, and then uses the salt to encrypt phone numbers and other similar identifiers. The salt never gets uploaded to the server and is known only to the phone for this purpose. Using the industry standard SHA-256 (Secure Hash Algorithm) and PBKDF2 (Password-Based Key Derivation Function 2) algorithms, an identifier is transformed into an 88-character anonymized string that can then be used in data analysis.

**Use cases.** At the time of writing, Beiwe is or has been used in tens of scientific studies on three continents across various fields, and there are likely several additional studies we are not aware of. Smartphone-based digital phenotyping is potentially very promising in behavioral and mental health [@onnela2016harnessing], and new research tools like Beiwe are especially needed in psychiatry [@torous2016new], where in the context of schizophrenia it has been used to predict patient relapse [@barnett2018relapse], compare passive and active estimates of sleep [@staples2017comparison], and characterize the clinical relevance of digital phenotyping data quality [@torous2018characterizing]. The platform has also been used to assess depressive symptoms in a transdiagnostic cohort [@pelligrini2021estimating] and to capture suicidal thinking during the COVID-19 pandemic [@fortgang2020increase]. There is an increasing amount of research on the use of Beiwe in neurological disorders, such as in the quantification of ALS progression [@berry2019design] and behavioral changes in people with ALS during the COVID-19 pandemic [@beukenhorst2021smartphone]. The platform has been used in the context of cancer to assess postoperative physical activity among patients undergoing cancer surgery [@panda2021smartphone], to capture novel recovery metrics after cancer surgery [@panda2020using], to enhance recovery assessment after breast cancer surgery [@panda2020smartphone], and to enhance cancer care [@wright2018hope]. Digital phenotyping and Beiwe have also been applied to quantifying mobility and quality of life of spine patients [@cote2019digital] and to study psychosocial well-being of individuals after spinal cord injury [@mercier2020digital].


# Acknowledgements

The Principal Investigator, Jukka-Pekka Onnela, is extremely grateful for his NIH Director’s New Innovator Award in 2013 (DP2MH103909) for enabling the crystallization of the concept of digital phenotyping and the construction of the Beiwe platform. He is also grateful to the members of the Onnela Lab. The authors also wish to acknowledge the contributions of the following individuals at Zagaran: Aaron Klein, Kevin Fan, Christopher McCarthy, Dor Samet, Alicia Fan, Rebecca Magazine Malamud, Jacob Klingensmith, and Benjamin Zagorsky.

# References
![Beiwe logo](https://github.com/onnela-lab/beiwe/blob/master/images/beiwe%20logo.png)
# Beiwe Research Platform (www.beiwe.org)
The Beiwe Research Platform is an open-source digital phenotyping platform designed for the collection and analysis of research-grade raw data from smartphone sensors and usage logs.  It was developed for biomedical research  and funded largely by a 2013 NIH Director’s New Innovator Award to Dr. Jukka-Pekka Onnela at the Harvard T.H. Chan School of Public Health.  More than just the Android and iOS apps, the Beiwe Research Platform also consists of three cloud-based components for collecting data, managing studies and performing data analysis which we call the backend.  Beiwe was released as open source at the end of 2017 and this document is to help you understand how Beiwe works, what's available and what's planned.  For background information, please see https://www.hsph.harvard.edu/onnela-lab/research/ and www.beiwe.org.

Of the Beiwe open source repositories, users need only implement the backend since apps, called Beiwe2,  are ready to go on the iOS App Store and the Google Play Store.   The Beiwe2 apps will prompt your study participants for your server name & will connect them to the correct study on your platform.  

Beiwe consists of five main github repositories as follows. 

# Repositories
## [beiwe](https://github.com/onnela-lab/beiwe)
This repository has no code but includes the wiki with the Beiwe documentation. Please visit the Beiwe wiki here: https://github.com/onnela-lab/beiwe/wiki

## [beiwe-backend](https://github.com/onnela-lab/beiwe-backend)
The Beiwe backend is the AWS-based source code for collecting data, managing studies and performing data analysis.  The backend is a modular, scalable system where: 

 - studies are created, configured and managed
 - patient IDs are created for enrolling study participants in individual studies
 - surveys, if used, are set up for each study
 - data is collected from all study participants' phones
 - data is analysed if the data analysis pipeline is configured 
 - data is downloaded for analysis

The data analysis pipeline infrastructure has not been fully tested at this time but is included in this repo.  It is designed to create containers to run specific code from the data analysis pipeline repo on specified studies and is completely configurable to your particular research.

For documentation on the Beiwe back-end and building your own backend server, [start here](https://github.com/onnela-lab/beiwe-backend/wiki/Deployment-instructions).  

For user documentation on Beiwe,  [start here](https://github.com/onnela-lab/beiwe/wiki).  
## [beiwe-android](https://github.com/onnela-lab/beiwe-android)
This repository contains the Beiwe Android app code. The Beiwe2 app is also available on the Google Play store to use with open source builds of the Beiwe backend.  The Beiwe2 app prompts for your server name and patient IDs from studies on your server. 
## [beiwe-ios](https://github.com/onnela-lab/beiwe-ios)
This repository contains the Beiwe iOS app code. The Beiwe2 app is also available on the Apple App store to use with open source builds of the Beiwe backend.  The Beiwe2 app prompts for your server name and patient IDs from studies on your server. 
## [Beiwe-Analysis](https://github.com/onnela-lab/Beiwe-Analysis)
This repository contains an evolving code base for running various types of analyses on data collected with the Beiwe research platform for smartphone-based digital phenotyping. We will merge methods from several distinct projects into this repository by the end of 2018. Note that this repository is distinct from the Beiwe data anlaysis pipeline. 

## beiwedata (https://github.com/onnela-lab/beiwedata)
This repository is provided by a doctoral student and contains a set of Python scripts designed to help analyze, and manipulate data generated by the Beiwe application.

> Written with [StackEdit](https://stackedit.io/).
<p align="left">
  <img width="264" height="99" src="https://github.com/onnela-lab/beiwe-backend/blob/main/beiwe-logo-color.png">
</p>

# Beiwe
The Onnela Lab at the Harvard T.H. Chan School of Public Health has developed the Beiwe (bee-we) research platform to collect smartphone-based high-throughput digital phenotyping data. The fully configurable open-source platform supports collection of a range of social, behavioral, and cognitive data, including spatial trajectories (via GPS), physical activity patterns (via accelerometer and gyroscope), social networks and communication dynamics (via call and text logs), and voice samples (via microphone). The platform consists of a front-end smartphone application for iOS (by Apple) and Android (by Google) devices and a back-end system, which supports a web-based study management portal for data processing and storage, based on Amazon Web Services (AWS) cloud computing infrastructure. Data analysis is increasingly identified as the main bottleneck; our data analysis platform, [Forest](https://github.com/onnela-lab/forest), makes sense of the data collected by Beiwe.

Beiwe can collect active data (participant input required) and passive data (participant input not required). Currently supported active data types for both Android and iOS include text surveys and audio recordings and their associated metadata. The questions, answers, and skip logic can be configured from the web-based dashboard. Passive data can be further divided into two groups: phone sensor data (e.g., GPS) and phone logs (e.g., communication logs). Beiwe collects raw sensor data and phone logs, which is crucial in scientific research settings. Beiwe has two front-end applications, one for Android (written in Java and Kotlin) and another for iOS (written in Swift & Objective-C). The Beiwe back-end, which is based on Amazon Web Services (AWS) cloud computing infrastructure, and runs on Python 3.8 using the Django webserver and ORM framework. It also uses several AWS services: primary S3 (for flat file storage), EC2 (servers), Elastic Beanstalk (scaling servers), and RDS (PostgreSQL).

Every aspect of data collection is fully customizable, including which sensors to sample, sampling frequency, addition of Gaussian noise to GPS location, use of Wi-Fi or cellular data for uploads, data upload frequency, and specification of surveys and their response options. Study participants simply download the Beiwe application from the app store and enter three pieces of information: a system-generated 8-character user ID, a system-generated temporary password, and an IP address of the back-end server. If no active data is being collected in the study (i.e., no surveys), this is the only time the participant will interact with the application. However, most studies make use of occasional self-reports or EMA, and some use the audio diary feature to collect rich data on lived experiences.

All Beiwe data is encrypted while stored on the phone awaiting upload and while in transit, and are re-encrypted for storage on the study server. During study registration, Beiwe provides the smartphone app with the public half of a 2048-bit RSA encryption key. With this key, the device can encrypt data, but only the server, which has the private key, can decrypt it. Thus, the Beiwe application cannot read its own temporarily stored data, and the study participant (or somebody else) cannot export the data. The RSA key is used to encrypt a symmetric Advanced Encryption Standard (AES) key for bulk encryption. These keys are generated as needed by the app and must be decrypted by the study server before data recovery. Data received by the cloud server is re-encrypted with the study master key and then stored.

Some of the data collected by Beiwe contain identifiers, such as phone numbers. The Beiwe app generates a unique cryptographic code, called a salt, during the Beiwe registration process, and then uses the salt to encrypt phone numbers and other similar identifiers. The salt never gets uploaded to the server and is known only to the phone for this purpose. Using the industry-standard SHA-256 (Secure Hash Algorithm) and PBKDF2 (Password-Based Key Derivation Function 2) algorithms, an identifier is transformed into an 88-character anonymized string that can then be used in data analysis.

A recent study found that 65% of medical studies were inconsistent when retested, and only 6% were completely reproducible. Reproducibility of studies using mobile devices may be even lower given the variability of devices, heterogeneity in their use, and lack of standardized methods for data analysis. All Beiwe study data collection settings, from sensors to surveys, are captured in a human readable JSON file. These files can be imported into Beiwe and exported out of Beiwe. To replicate a study, the investigator can simply upload an existing configuration file.

Cite the code: [![DOI](https://zenodo.org/badge/53344506.svg)](https://zenodo.org/badge/latestdoi/53344506)

# Setup instructions

## Configuring SSL
Because Beiwe often deals with sensitive data covered under HIPAA, it's important to add an SSL certificate so that web traffic is encrypted with HTTPS.

The setup script [uses AWS Certificate Manager to generate an SSL certificate](http://docs.aws.amazon.com/acm/latest/userguide/gs-acm-request.html).  AWS Certificate Manager [will check that you control the domain by sending verification emails](http://docs.aws.amazon.com/acm/latest/userguide/gs-acm-validate.html) to the email addresses in the domain's WHOIS listing.

## Configuring Firebase
To initialize the Firebase SDK, [generate a private key file](https://firebase.google.com/docs/admin/setup#initialize-sdk).
Rename the file firebase_cloud_messaging_credentials.json and place it in the project root.

***

# Configuration settings

### Mandatory Settings

If any of these environment options are not provided, Beiwe will not run. Empty strings and None  are considered invalid.

```
    FLASK_SECRET_KEY - a unique, cryptographically secure string
    AWS_ACCESS_KEY_ID - AWS access key for S3
    AWS_SECRET_ACCESS_KEY - AWS secret key for S3
    S3_BUCKET - the bucket for storing app-generated data
    SYSADMIN_EMAILS - a comma separated list of email addresses for recipients of error reports. (whitespace before and after addresses will be ignored)
    RDS_DB_NAME - postgress database name (the name of the database inside of postgres)
    RDS_USERNAME - database username
    RDS_PASSWORD - database password
    RDS_HOSTNAME - database IP address or url
    S3_ACCESS_CREDENTIALS_USER - the user id for s3 access for your deployment
    S3_ACCESS_CREDENTIALS_KEY - the secret key for s3 access for your deployment
```

### Optional Settings
There are additional settings that you will find documented in the [config/settings.py](https://github.com/onnela-lab/beiwe-backend/blob/main/config/settings.py) file.

We _strongly_ recommend adding Sentry DSNs to all your Beiwe servers.  Without these there is very little data to work with when something goes wrong, and we won't be able to assist.

***

# Development setup
How to set up beiwe-backend running on a development machine (NOT a production instance!  For a production instance,
see https://github.com/onnela-lab/beiwe-backend/wiki/Deployment-Instructions---Scalable-Deployment)

#### Before starting:
While it is possible to run your development environment inside of the system Python environment, this practice is _strongly discouraged_.  We recommend familiarizing yourself with one of the following: Python's [venv](https://docs.python.org/3/tutorial/venv.html) library (basic virtual environments), [Pyenv](https://github.com/pyenv/pyenv) (allows for compiling particular target versions of Python, plus some quality-of-life command-line shell integrations), or [Conda](https://docs.conda.io/en/latest/) (another option, includes integrations with non-Python libraries).  Note also that the codebase expects at least Python version 3.8.

1. `sudo apt-get update; sudo apt-get install postgresql libpq-dev`
2. `pip install --upgrade pip setuptools wheel`
3. `pip install -r requirements.txt`
4. Create a file for your environment variables that contains at least these:
    ```
    export DOMAIN_NAME="localhost://8080"
    export FLASK_SECRET_KEY="asdf"
    export S3_BUCKET="a"
    export SYSADMIN_EMAILS="sysadmin@localhost"
    ```
    I usually store it at `private/environment.sh`.  Load up these environment variables by running `source private/environment.sh` at the Bash prompt.

For additional tips on running a local development enironment please see [this wiki page](https://github.com/onnela-lab/beiwe-backend/wiki/Tips-For-Local-Beiwe-Development).  If you are having difficulty getting started, or believe you could assist with any issues of documentation, please [post an issue with a documentation tag](https://github.com/onnela-lab/beiwe-backend/labels/documentation).

### Local Celery setup
**Update**: it is no longer necessary to use Celery for local testing, though you still need it to be installed in your Python environment in order to avoid import errors.  A full test of Celery requires the full setup below, including installing `rabbitmq`, but as long as the file for the rabbitmq host server IP and password (`manager_ip` in the root of the repository) is missing you will instead be presented with output similar to the example shell session below, indicating a that you are running in a _much_ more convenient single-threaded local testing mode:

```
In [1]: from services.celery_data_processing import *
task declared, args: (), kwargs:{'queue': 'data_processing'}
Instantiating a FalseCeleryApp for celery_process_file_chunks.
```

For those souls brave enough to run the entire broker queue and Celery task dispatch machinery locally, here are our best instructions.  Caveat: this configuration is based on a one that is known to work on Ubuntu 18.04, and is potentially incompatible with the version of RabbitMQ provided in Ubuntu 20.04. Also, due to the use of the system `service` command it is incompatible with the varient of Ubuntu for use on the Windows Subsystem for Linux.  Have at:

1. Install RabbitMQ (https://docs.celeryproject.org/en/latest/getting-started/backends-and-brokers/rabbitmq.html#broker-rabbitmq)
    1. Edit `/etc/rabbitmq/rabbitmq-env.conf` and add the line `NODE_PORT=50000`
    2. Restart RabbitMQ like this in the Bash shell: `time sudo service rabbitmq-server restart` (`time` isn't necessary, but it tells you that the command has finished, and how long the command took to execute... which can be... random and excessive?)
2. `pip install -r requirements_data_processing.txt` (this will install Celery using pip)
3. Create a file called `manager_ip` in the top level of your `beiwe-backend` repo, and enter these two lines in it.  Do not provide a trailing new-line character.
    ```
    127.0.0.1:50000
    [YOUR DESIRED PASSWORD]
    ```
    Where the password is the one you set when setting up RabbitMQ
4. run this command to create a user for rabbitmq: `rabbitmqctl add_user beiwe [INSERT THAT PASSWORD HERE]`
5. run this command to allow the beiwe user access to the appropriate queues: `sudo rabbitmqctl set_permissions -p / beiwe ".*" ".*" ".*"`
6. If you intend to test Firbase push notifications you will need to upload functional firebase credentials on the local website interface.
7. To execute push notification tasks run this command _while inside the root of the repo_: `celery -A services.celery_push_notifications worker -Q push_notifications --loglevel=info -Ofair --hostname=%%h_notifications --concurrency=20 --pool=threads`
8. To run data processing tasks run this command _while inside the root of the repo_: `celery -A services.celery_data_processing worker -Q data_processing --loglevel=info -Ofair --hostname=%%h_processing`
9. To run forest tasks run this comand _while inside the root of the repo_: `celery -A services.celery_forest worker -Q forest_queue --loglevel=info -Ofair --hostname=%%h_forest` (Forest is still in beta.)
10. Run this command to dispatch new tasks, which will then be consumed by the Celery processes, _while inside the root of the repo_. `python services/cron.py five_minutes`



### Forest

Warning: The Forest integration is still in beta; running Forest may cause significant data processing costs.
