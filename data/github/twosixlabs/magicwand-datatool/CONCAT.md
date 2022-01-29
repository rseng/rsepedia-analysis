# Magicwand Tests

*Welcome to the Magicwand tests!*

## Running Magicwand Tests

To run the tests locally, first install the tests-specific requirements with `pip` using the `requirements-test.txt` file:

```
$ pip install -r magicwand-data-generator/requirements-tests.txt
```


Tests can then be run as follows from the project `root`:

```bash
$ make test
```

The Makefile uses the `pytest` runner and testing suite as well as the coverage library.

These are automated tests that are run before each new release
# CHANGELOG

## v0.4.0

- Removed view runs
- Added CIC convert command
- CIC conversion done post run automatically
- Documentation updated
- General code clean up


## v0.3.0
- Code refactor release, no new features 

## v0.2.0

- Updated documentation
- Quality of Life Improvements
- Variable for Locust stagger start time

## v0.1.0

Lots of changes, this is a release that would be safe for external users to test.

### Added new attacks
* httpflood
* synflood

### Code refactored
* less yaml files
* smaller json config files
* removed dead code and files

### Locust configuration 
* added random keep alive toggle
* added random start time for locust clients

### SUT configuration
* can configure snap length for pcap capture

### Documentation 
* Updated mkdocs

### Data Verification 
* added data tests for apachekill
* added ip check tests

### Data viewer updated
 * Better readability
 * always up to date



## v0.0.2

- Added mkdocs-based documentation	
- Added tests	
  - Test that `init` executed properly	
  - Test that an apachekill `run` created the correct files	
- Added new attacks 	
  - sockstress	
  - slowhttptest (slowloris, rudy, slowread)	
  - goloris	
- Added version tagging for data generated	
- Added local image building
# Magicwand Data Tool

**PCAP generation platform to facilitate machine learning model development for identifying benign traffic vs. malicious "low volume" denial of service traffic.**

## What is Magicwand

Magicwand is a platform to provide high-quality, reliable, and reproducible data sets for low-and-slow DDoS attacks. With the use of Docker images and customizable JSON files, users can generate a multitude of network traffic PCAPS.

For complete documentation on Magicwand, a gallery of available attacks, the configuration guide, tutorials, and teaching resources, frequently asked questions, and more, please read our [documentation](https://magicwand-datatool.readthedocs.io/en/latest/).

### What is a "low and slow" DDoS

Low and slow attacks stealthily degrade server performance through cleverly crafted transmissions of data. A normal DDoS uses volumetric transmissions to overwhelm a server.

## Why develop Magicwand

In the network security space, questions like these are common amongst researchers...

> ![Magicwand Apachekill Run](docs/images/why_magicwand.png)


Network security datasets are hard to come by due to data privacy concerns. This struggle to capture live traffic and use it for research purposes is limited to very static and unreproducible datasets provided on an ad-hoc basis. To fix the stale dataset problem, we have developed Magicwand, as we to provide researchers with high quality data they can use for network security analysis.

## Installing Magicwand

### Dependencies

The following are dependencies need to properly run Magicwand

* docker (https://docs.docker.com/get-docker/)
* docker-compose (https://docs.docker.com/compose/install/)
* python3.6+ (https://www.python.org/downloads/)
* tshark (https://www.wireshark.org/docs/man-pages/tshark.html)

### Hardware Requirements

* >= 8GB of RAM
* >= 2 CPUS

### Installing Magicwand

Magicwand is compatible with Python 3.6 or later. The simplest way to install Magicwand and its dependencies is from PyPI with pip, Python's preferred package installer.

**Note** Depending on your permissions (Docker/Python), you may need to run some commands as sudo (e.g. `sudo bash bash scripts/pull_images.sh`). It is recommended to run without root first, and only run as root if necessary.

```bash
python3 -m virtualenv ./venv
source venv/bin/activate
```

It is also recommended to utilize a vitural environment when installing python packages to avoid compatibility issues.

```bash
pip install magicwand
```

Note that Magicwand is an active project and routinely publishes new releases. In order to upgrade Magicwand to the latest version, use pip as follows

```bash
pip install -U magicwand
```

Magicwand can also be installed from source.

```bash
make -C magicwand-data-generator/ install
```

In addition to the python package, Magicwand leverages prebuilt docker images to run experiments

You can pull from docker hub

```bash
bash scripts/pull_images.sh
```

### Quick Start

Here's how you can quickly use magicwand

### 1. Install Magicwand CLI Tool

```bash
pip install magicwand
```

### 2. Create Test Folder

```bash
magicwand init --project test
cd test
```

### 3. Run Calibration Command
``` bash
magicwand calibrate --attack apachekill
```

### 4. Start Runs
```bash
magicwand run --config configs/mw_locust-apachekill.json --count 1 --data_version test_runs
```

For To get started using the Magicwand Data Generator, please visit our [documentation](https://magicwand-datatool.readthedocs.io/en/latest/).

## Contributing to Magicwand

Magicwand is an open source project that is supported by a community who will gratefully and humbly accept any contributions you might make to the project. Large or small, any contribution makes a big difference; and if you've never contributed to an open source project before, we hope you will start with Magicwand!

If you are interested in contributing, check out our contributor's guide. Here are some of the many ways to contribute:

* Submit a bug report or feature request on GitHub Issues.
* Assist us with user testing.
* Add a new attack to our repository
* Add to the documentation or help with our website,
* Write unit or integration tests for our project.
* Answer questions on our issues, mailing list, Stack Overflow, and elsewhere.
* Translate our documentation into another language.
* Write a blog post, tweet, or share our project with others.
* Teach someone how to use Magicwand.

As you can see, there are lots of ways to get involved and we would be very happy for you to join us! The only thing we ask is that you abide by the principles of openness, respect, and consideration of others as described in the Python Software Foundation Code of Conduct.

For more information, checkout the [CONTRIBUTING.md](CONTRIBUTING.md) file in the root of the repository.

## Magicwand Datasets

Example datasets have been published to [Zenodo](https://zenodo.org/record/4774258#.YKV-1pNKjmE) following the steps documented in the `Quick Start`

## Citing Magicwand

We would be glad if you used Magicwand in your scientific publications! If you do, please cite us using the citation guidelines.

## Affiliations


[<img src="https://www.twosixlabs.com/wp-content/uploads/TwoSix-Logo-Media-Kit-Reverse-Red.jpg">](https://www.twosixlabs.com/)
# Contributing to Magicwand

**NOTE: This document is a "getting started" summary for contributing to the Magicwand project.** Please make sure to read this page carefully to ensure the review process is as smooth as possible and to ensure the greatest likelihood of having your contribution be merged.


## How to Contribute

Magicwand is an open source project that is supported by a community who will gratefully and humbly accept any contributions you might make to the project. Large or small, any contribution makes a big difference; and if you've never contributed to an open source project before, we hope you will start with Magicwand!

Principally, Magicwand development is about the addition and creation of new traffic generation images. We'll discuss in detail how to add new images later.

Beyond creating images, there are many ways to contribute:

* Submit a bug report or feature request on GitHub Issues.
* Assist us with user testing.
* Add a new attack to our repository
* Add to the documentation or help with our website, 
* Write unit or integration tests for our project.
* Answer questions on our issues, mailing list, Stack Overflow, and elsewhere.
* Translate our documentation into another language.
* Write a blog post, tweet, or share our project with others.
* Teach someone how to use Magicwand.

As you can see, there are lots of ways to get involved and we would be very happy for you to join us! The only thing we ask is that you abide by the principles of openness, respect, and consideration of others as described in the [Python Software Foundation Code of Conduct](https://www.python.org/psf/codeofconduct/).

## Getting Started on GitHub

Magicwand is hosted on GitHub at https://github.com/twosixlabs/magicwand-datatool

The typical workflow for a contributor to the codebase is as follows:

1. **Discover** a bug or a feature by using Magicwand.
2. **Discuss** with the core contributes by [adding an issue](https://github.com/twosixlabs/magicwand-datatool/issues).
3. **Fork** the repository into your own GitHub account.
4. Create a **Pull Request** first thing to [connect with us](https://github.com/twosixlabs/magicwand-datatool/pulls) about your task.
5. **Code** the feature, write the documentation, add your contribution.
6. **Review** the code with core contributors who will guide you to a high-quality submission.
7. **Merge** your contribution into the Magicwand codebase.

We believe that *contribution is collaboration* and therefore emphasize *communication* throughout the open source process. We rely heavily on GitHub's social coding tools to allow us to do this. For instance, we use GitHub's [milestone](https://help.github.com/en/articles/about-milestones) feature to focus our development efforts for each Magicwand semester, so be sure to check out the issues associated with our [current milestone](https://github.com/twosixlabs/magicwand-datatool/milestones)!

Once you have a good sense of how you are going to implement the new feature (or fix the bug!), you can reach out for feedback from the maintainers by creating a [pull request](https://github.com/twosixlabs/magicwand-datatool/pulls). Please note that if we feel your solution has not been thought out in earnest, or if the PR is not aligned with our [current milestone](https://github.com/twosixlabs/magicwand-datatool/milestones) goals, we may reach out to ask that you close the PR so that we can prioritize reviewing the most critical feature requests and bug fixes.

Ideally, any pull request should be capable of resolution within 6 weeks of being opened. This timeline helps to keep our pull request queue small and allows Magicwand to maintain a robust release schedule to give our users the best experience possible. However, the most important thing is to keep the dialogue going! And if you're unsure whether you can complete your idea within 6 weeks, you should still go ahead and open a PR and we will be happy to help you scope it down as needed.

If we have comments or questions when we evaluate your pull request and receive no response, we will also close the PR after this period of time. Please know that this does not mean we don't value your contribution, just that things go stale. If in the future you want to pick it back up, feel free to address our original feedback and to reference the original PR in a new pull request.

### Forking the Repository

The first step is to fork the repository into your own account. This will create a copy of the codebase that you can edit and write to. Do so by clicking the **"fork"** button in the upper right corner of the Magicwand GitHub page.

Once forked, use the following steps to get your development environment set up on your computer:

1. Clone the repository.

    After clicking the fork button, you should be redirected to the GitHub page of the repository in your user account. You can then clone a copy of the code to your local machine.

    ```
    $ git clone https://github.com/twosixlabs/magicwand-datatool.git
    $ cd magicwand
    ```

    Optionally, you can also [add the upstream remote](https://help.github.com/articles/configuring-a-remote-for-a-fork/) to synchronize with changes made by other contributors:

    ```
    $ git remote add upstream https://github.com/twosixlabs/magicwand-datatool
    ```

    See "Branching Conventions" below for more on this topic.


2. Switch to the develop branch.

    The Magicwand repository has a `develop` branch that is the primary working branch for contributions. It is probably already the branch you're on, but you can make sure and switch to it as follows:

    ```
    $ git fetch
    $ git checkout develop
    ```

At this point you're ready to get started writing code!

### Branching Conventions

The Magicwand repository is set up in a typical production/release/development cycle as described in "[A Successful Git Branching Model](http://nvie.com/posts/a-successful-git-branching-model/)." The primary working branch is the `develop` branch. This should be the branch that you are working on and from, since this has all the latest code. The `master` branch contains the latest stable version and release, _which is pushed to PyPI_. No one but maintainers will push to master or develop.

**NOTE:** All pull requests should be into the `magicwand/develop` branch from your forked repository.

You should work directly in your fork and create a pull request from your fork's develop branch into ours. We also recommend setting up an `upstream` remote so that you can easily pull the latest development changes from the main Magicwand repository (see [configuring a remote for a fork](https://help.github.com/articles/configuring-a-remote-for-a-fork/)). You can do that as follows:

```
$ git remote add upstream https://github.com/https://github.com/twosixlabs/magicwand-datatool`
$ git remote -v
origin    https://github.com/YOUR_USERNAME/YOUR_FORK.git (fetch)
origin    https://github.com/YOUR_USERNAME/YOUR_FORK.git (push)
upstream  https://github.com/https://github.com/twosixlabs/magicwand-datatool (fetch)
upstream  https://github.com/https://github.com/twosixlabs/magicwand-datatool (push)
```

When you're ready, request a code review for your pull request. Then, when reviewed and approved, you can merge your fork into our main branch. Make sure to use the "Squash and Merge" option in order to create a Git history that is understandable.

**NOTE to maintainers**: When merging a pull request, use the "squash and merge" option and make sure to edit the both the subject and the body of the commit message so that when we're putting the changelog together, we know what happened in the PR. I recommend reading [Chris Beams' _How to Write a Git Commit Message_](https://chris.beams.io/posts/git-commit/) so we're all on the same page!

Core contributors and those who are planning on contributing multiple PRs might want to consider using feature branches to reduce the number of merges (and merge conflicts). Create a feature branch as follows:

```
$ git checkout -b feature-myfeature develop
$ git push --set-upstream origin feature-myfeature
```

Once you are done working (and everything is tested) you can submit a PR from your feature branch. Synchronize with `upstream` once the PR has been merged and delete the feature branch:

```
$ git checkout develop
$ git pull upstream develop
$ git push origin develop
$ git branch -d feature-myfeature
$ git push origin --delete feature-myfeature
```

Head back to Github and checkout another issue!

## Developing Magicwand Components 

In this section, we'll discuss the basics of developing magicwand components. This of course is a big topic, but hopefully these simple tips and tricks will help make sense.

### Magicwand Components

There are four basic types of components

- **Attacks** are components that attempt to shut down the SUT with malicious traffic
- **Benign** are components that send benign traffic to the SUT
- **SUT** The System under test that retrieves the attack and benign traffic
- **Sensors** are components that do auxiliary tasks such as RTT monitoring 

Creating a new component requires the following 

- Docker Image,
- Docker-compose script
- JSON configuration file,
- MwComponent python class
- MwRunner Updates

A barebones implementation is as follows...


#### Docker Image

This image will run a python script on start to send traffic to the SUT. You will need a Dockerfile, a begin-test.sh script, a python script, and a docker-compose.yml.  
The image can be placed in the `images` folder

```Dockerfile
FROM ubuntu:18.04
RUN apt-get update && apt-get install net-tools python3 python3-pip iproute2 vim -y
RUN pip3 install --upgrade pip; pip install requests


ADD begin-test.sh /usr/local/bin/
ADD example.py /usr/local/bin/
RUN chmod a+x /usr/local/bin/begin-test.sh
WORKDIR /usr/local/bin/
CMD ["./begin-test.sh"]
```


**begin-test.sh**
```bash
#!/bin/sh

_term() {
    echo [begin-test.sh] Example Caught SIGTERM signal!
    pkill -f python3.6\ example
}

trap _term TERM

#CLIENT OR ATTACK 
LOCAL_IP=`ip route show | grep src | cut -f 9  -d ' '`
echo "ip,type,subtype" >> /home/ip_map_rtt.csv
echo "$LOCAL_IP,client,rtt" >> /home/ip_map_client.csv

echo "timestamp,rtt" >> /home/rtt_stats.csv

echo "[begin-test.sh] python rtt_tracker.py"
python3.6 rtt_tracker.py

tail -f /dev/null
```

**example.py**
```python
import requests
import os
import time
import datetime


TARGET_URL = os.environ['TEST_TARGET']
DUR = os.environ['TEST_DURATION']
MY_PARM = os.environ['MY_PARM']
def main():

    counter = 0
    start = time.time()

    while counter < int(DUR):

        #curr_time = datetime.datetime.now().strftime("%H:%M:%S")

        try:
            r = requests.get("http://"+TARGET_URL,timeout=2)
            roundtrip = r.elapsed.total_seconds()            
        except:
            roundtrip = 2

        time.sleep(int(MY_PARM))
        counter += 1

if __name__ == "__main__":
    main()

```


#### MwComponent Folder

The following files will be placed inside the appropriate MwComponent folder in `magicwand-data-generator/magicwand/magicwand_components`   
This example will be for benign so our folder iwll be `magicwand-data-generator/magicwand/magicwand_components/benign/example_client/`


##### Docker-compose script

This is an example docker-compose script, that will start the component during runtime

**exmaple.yml**
```
  mw-client-example:
    image: MY_IMAGE:latest
    privileged: true
    depends_on:
      - "mw-sut-apachewp"
    environment:
      - TEST_DURATION=${CURR_TEST_DURATION}
      - MY_PARAM=${CURR_MY_PARAM}
      - TEST_TARGET=mw-sut-apachewp:80
    volumes:
      - ./${CURR_RUN}:/home/
```
Add this file to the `magicwand-data-generator/magicwand/magicwand_components/benign/example_client/example_client.yml`

##### MwComponent

This module contains the execution logic for our example client. It inherits from the MwComponent abstract class. Here is a full example
```python
"""
    Purpose:
        This file contains the class for example client
"""
import logging
import os
import time

from typing import Any, Dict
from magicwand.magicwand_components.mw_component import MwComponent
from magicwand.magicwand_utils.magicwand_utils import get_logger
from magicwand.magicwand_utils.magicwand_utils import load_json


class ExampleClient(MwComponent):

    name = "example_client"

    def __init__(self, log_level=logging.INFO) -> None:
        """
        Purpose:
            Init Component
        Args:
            log_level: Verbosity of the logger
        Returns:
            self: The MwComponent
        """
        # get logger
        # self._logger = get_logger("mw-log", log_level)
        super().__init__(log_level=log_level)

        # set config to expected spot
        self._config = load_json("magicwand_components/benign/example_client.json")

    @property
    def config(self) -> Dict[str, Any]:
        """
        Purpose:
            Get config
        Args:
            N/A
        Returns:
            config: The json config for the component
        """
        return self._config

    @config.setter
    def config(self, val: Dict[str, Any]) -> None:
        """
        Purpose:
            Set config
        Args:
            val: value of the config
        Returns:
            N/A
        """
        self.config = val

    def set_env_variables(self) -> int:
        """
        Purpose:
            Set environment variables 
        Args:
            N/A
        Returns:
            stats: 0 if worked, -1 if failed
        """
        try:
            os.environ["MY_PARM"] = str(self.config["MY_PARM"])
        except Exception as error:
            self.logger.error(error)
            return -1

        return 0

    def verify(self, config: Dict[str, Any], post_run_data: Dict[str, Any]) -> bool:
        """
        Purpose:
            Verify if the component worked during the run
        Args:
            config: Run config options
            post_run_data: Data for verifications
        Returns:
            passed: True if passed, False if failed
        """
        return True

```

Add this file to the `magicwand-data-generator/magicwand/magicwand_components/benign/example_client/example_client.py`

##### JSON configuration file

Here is an example JSON file

```json
{"benign":"example_client", "compose-file": "magicwand_components/benign/example_client.yml", "MY_PARM": 2}
```

Add this file to the `magicwand-data-generator/magicwand/magicwand_components/benign/example_client/example_client.json`


##### Final steps

We need to add some code to extra spots as well.

Create a `__init__.py` and add

```python
from .example_client import *
```

Also add this to benign's `__init__.py`

```python
from .example_client import *
```

In mw_runner.py add to the component dictionaries

```python
mw_components: Mapping[str, Type[MwComponent]] = {
    "apachekill": Apachekill,
    "sockstress": Sockstress,
    "goloris": Goloris,
    "sht_rudeadyet": Sht_rudeadyet,
    "sht_slowread": Sht_slowread,
    "sht_slowloris": Sht_slowloris,
    "httpflood": Httpflood,
    "synflood": Synflood,
    "mw_locust": MwLocust,
    "mw_rtt_sensor": MwRttSensor,
    "mw_apache_wp": MwApacheWp,
    "example_client": ExampleClient    <------ Update
}


valid_values = {
    "attack": [
        "apachekill",
        "sockstress",
        "goloris",
        "sht_rudeadyet",
        "sht_slowread",
        "sht_slowloris",
        "httpflood",
        "synflood",
    ],
    "benign": ["mw_locust","example_client"],          <------ Update
    "sut": ["mw_apache_wp"],
    "rtt": ["mw_rtt_sensor"],
}
```


#### Test your component

Once all of these are in place and your image is built. You can rebuild Magicwand

```bash
pip install --editable .
```

Create a new project

```bash
magicwand init --project new_example
```

Start your example

```bash
#make new config file

{
  "benign": "example_client",
  "sut": "mw_apache_wp",
  "rtt": "mw_rtt_sensor",
  "run_type": "example_client-only",
}


magicwand run --config configs/example_client.json
```

I know this is daunting and may require some give and take. Feel free to leave an issue on the GitHub page if you need help creating a new image or attack. We will also improve the process to add new images as well.

---
title: 'Magicwand: A platform to provide high-quality, reliable, and reproducible data sets for low-and-slow DDoS attacks.'
tags:
  - Python
  - Cybersecurity
  - Network Traffic
  - Distributed Denial of Service
  - DDoS
  - LSDDoS
  - Machine Learning
authors:
  - name: Banjo Obayomi
    affiliation: 1
  - name: Christopher H. Todd
    affiliation: 1
  - name: Lucas Cadalzo
    affiliation: 1
  - name: David Killian
    affiliation: 1
  - name: Anthony C. Wong
    affiliation: 2
affiliations:
 - name: Two Six Labs
   index: 1
 - name: Unaffiliated
   index: 2
date: xx May 2021
bibliography: paper.bib
---

# Summary

`Magicwand` is a Python platform to provide high-quality, reliable, and replicable data sets for low-and-slow distributed denial-of-service (LSDDoS) attacks. Using a Command line interface (CLI), researchers can generate packet capture (.pcap) data with labels for attackers and benign clients. With the use of Docker, and JavaScript Object Notation (JSON) files, `Magicwand` provides a customizable and extensible data generation framework to generate a multitude of network traffic PCAPs.

# Statement of Need

Distributed denial-of-service (DDoS) attacks remain a pervasive cybersecurity threat [@krebs:2016], [@dota:2015], and [@russia:2008]. While various commercial services such as Akamai and CloudFlare are effective in mitigating traditional, high-volume DDoS attacks, LSDDoS attacks present a different kind of threat [@makrushin:2013]. These attacks stealthily degrade server performance through cleverly crafted transmissions of data, rather than brute-force floods of traffic. The adversary sends low-bandwidth requests that slowly transmit data and consume a disproportionate amount of the server's resources. The relatively small amount of traffic needed to cause a denial-of-service has two important consequences: these attacks are both easier to launch and more difficult to detect, in comparison to high-volume DDoS.

In order to create a robust LSDDoS defense, researchers in the security community must have access to sufficient data capturing how these attacks behave. The most widely-adopted public source of data for LSDDoS was generated by the Canadian Institute for Cybersecurity (CIC) at the University of New Brunswick [@sharafaldin:2018]. However, the vast majority of instances in this dataset are from high-volume attacks, and the dataset does not contain any benign traffic directed towards the victim host. Additionally, the range of LSDDoS attacks employed is limited. In an effort to facilitate the generation of more thorough LSDDoS datasets, we are introducing our traffic generation tool `Magicwand`.

Through its development, these key goals informed the design and implementation:

* The user must not be required to be an expert in low-and-slow DDoS attacks.
* Data generated must be representative of low-and-slow DDoS attack behavior and their impact on victim servers.
* Data generated should be validated through a series of automated tests and checks.
* The software must be able to generate data on any hardware that meets a minimum set of requirements.
* The software should be extensible to enable new components in the future.
* The software should allow for user flexibility in configuration of LSDDoS scenarios.

# Extensibility

`Magicwand` was designed as a framework that can be extended based on project needs; to this end, the tool is comprised of extensible *components* that can be considered a foundation for expanding the tool's capabilities. The four core components of the system are the attack clients, benign clients, systems under threat (SUTs), and sensors. Each of these components was designed to be operationalized into any `Magicwand` experiment.

## Attack Clients

The set of implemented attack images consists of Slowloris [@goloris:2014], R U Dead Yet [@sht:2011], Slow Read [@sht:2011], Apache Killer [@killapache:2011], Sockstress [@sockstress:2012], and a second implementation of Slowloris [@goloris:2014].  The JSON files provide configurations that can scale the intensity of attack and adjust other attack-specific-specific parameters (e.g. threads per connection for Apache Killer).

## Benign Clients

For benign traffic generation, `Magicwand` uses customized agents, based on Locust [@locust:2019], that continuously browse the SUT, traversing randomly from link to link within the SUT. The JSON files allow for the configuration of parameters such as clients per container, connection runtime, rate of activity, and more.

## System Under Threat (SUT)

The primary SUT included with `Magicwand` is an Apache WordPress server. We choose this stack for two reasons: first, it is widely used on the web, ensuring our system has wide applicability. Second, there is a substantial set of known LSDDoS attack implementations for this stack.

## Sensors

`Magicwand` includes sensors that use tcpdump to capture real-time network traffic and provide a PCAP for each experiment. An additional sensor captures Round-Trip-Time (RTT) data of the SUT to understand the impact of the attack and what a client would experience during the experiment.

## Developing A New Component

Developing a new `Magicwand` component requires extending the Python `MwComponent` abstract base class. The critical methods/attributes are `validate`, `run`, and `config` and must be implemented for a new component. The developer must also create a `configuration.json` file that defines the parameters for the component and default values.

## Component Containerization

`Magicwand` experiments rely on Docker containers to run any component. For each component, a user will need to create a `Dockerfile` that will run the component during the experiment, a `docker-compose.yml` file that will handle parameters and runtime-config for the component during the experiment.

# Experiments

Experiments are defined by the components that will run, the timeline for each component to execute, and customizable configuration for each component. Each experiment will run locally and generate data that includes PCAPs of the network traffic generated, data from each sensor component, and metadata from the experiment.

# System Calibration

To ensure the data is representative of LSDDoS attacks and realistic benign clients, `Magicwand` includes calibration functionality.

For benign components, `Magicwand` ensures that the amount of benign traffic flowing towards the SUT is not so large that it causes congestion, absent any malicious traffic.

For attack components, `Magicwand` finds the configuration required to ensure an attack is strong enough to induce a denial-of-service condition on the SUT. For example, Apache Killer calibration process will determine the parameters needed to yield a mean RTT greater than 2x the mean RTT of experiments with only benign traffic.

Attack calibration works by launching two experiments: one with only benign traffic and the other with only attack traffic. We then compare various metrics between the two experiments that detail the state of the SUT, such as RTT and memory consumption. The SUT parameters are held constant between each experiment to compare how the attack and benign traffic changes based on the inputs.  If the ratio between the metrics of the benign and attack experiments indicates that the attack succeeded in inducing a denial-of-service condition, then the JSON configuration file is saved.

# Data Validation

We validate each experiment to ensure representative and expected data was produced. After each experiment, `Magicwand` conducts a series of checks based on the input parameters against the generated data and run logs. Each component has its own validation functionality.

Attack components verify that each attack yielded expected signatures in the data; an example would be Apache Killer generating packets that send overlapping bytes in the resulting PCAP. Each attack has specific validations that are tailored to how the attack works.

Benign client validation confirms that the clients performed the tasks they were configured to do; an example is verifying that an experiment with 10 benign clients attempted to spin up 10 connections during the experiment.

Sensor validation confirms that the sensors ran and yielded the data each sensor is responsible for and that the data is not corrupted. SUT validation confirms that SUT receives all expected trafic and that it ran as expected through the experiment.

While we have been systematic and methodical in crafting these checks to make them as comprehensive as possible, these are not exhaustive. This is where researchers can define validation methods to ensure their `Magicwand` experiments are producing data that they would expect.

# Acknowledgements

This research was developed with funding from the Defense Advanced Research Projects Agency (DARPA) under Contract # HR0011-16-C-0060. This document was cleared for release under Distribution Statement” A” (Approved for Public Release, Distribution Unlimited). The views, opinions, and/or findings expressed are those of the authors and should not be interpreted as representing the official views or policies of the Department of Defense of the U.S. Government.

# References# Magicwand Docker Images

This directory contains the code for the magicwand CLI.  
To learn more about the CLI please visit our documentation at https://magicwand-datatool.readthedocs.io/en/latest/
# CSL Magicwand User Guide

**This guide will explain how to generate datasets for CSL**

## What is Magicwand 
Magicwand is a platform to provide high-quality, reliable, and reproducible data sets for low-and-slow DDoS attacks. With the use of Docker images and customizable JSON files, users can generate a multitude of network traffic PCAPS, that are converted to CSVs labeled with attack or benign network flows.

## Installing Magicwand

For full details checkout the official [README.md](../README.md#installing-magicwand) file in the root of the repository.

### Quick install
```
git clone https://github.com/twosixlabs/magicwand-datatool
cd images
./build_all.sh
cd ..
cd magicwand-data-generator
pip install --editable .
```

### Create CSL datasets

Here's how you can begin creating CSL datasets


**1. Create Dataset Folder**  

This step initializes a magicwand 'environment' to begin running your experiments 

```bash
localhost$ magicwand init --project csl_datasets
localhost$ cd csl_datasets
```

**2. Start an experiment**  

Once your environment is initialized, you can begin creating data using the run command


localhost$ magicwand run --config configs/mw-locust-apachekill.json --count 1 –data_version test_runs


The run config file specifies what components to use for your experiment  
Here is an example of a run config file

```
{
    "attack": "apachekill",
    "benign": "mw_locust",
    "sut": "mw_apache_wp",
    "rtt": "mw_rtt_sensor",
    "run_type": "mw-locust-apachekill"
}
```

The count parameter specifies how many times to run the experiment  
The data_version parameter specifies where to save the data



**3. Evaluate data**  

Each run will output metrics of the final CSV benign vs attack split for the cic_flow_labeled.csv

```
2020-10-21 13:44:29,005 - INFO -  mw-log - Benign Stats:
2020-10-21 13:44:29,005 - INFO -  mw-log - 744 flows out of 1750 flows: 0.43
2020-10-21 13:44:29,005 - INFO -  mw-log - Attack Stats:
2020-10-21 13:44:29,005 - INFO -  mw-log - 1005 flows out of 1750 flows: 0.57
```

If the resulting CSV does not have the desired stats, then we need to edit the component config 

**4. Edit Config files**  

### Editing Attack Configs
Each component config has its own set of parameters to configure in order to affect the number of flows in each experiment.

All components configs are located in the magicwand_components folder.
In this example we will edit the Apachekill Config located here `magicwand_components/attacks/apachekill.json`

``` bash
cat magicwand_components/attacks/apachekill.json | jq
{
  "attack_options": {
    "ak_num_threads": 50,              <- Play around with
    "ak_num_ips": 20,                  <- Play around with
    "attack_delay": 15,                <- Play around with
    "ak_duration": 30                  <- Play around with.  (My experience so far time is biggest factor)
  },
  "attack": "apachekill",
  "compose-file": "magicwand_components/attacks/apachekill.yml"
}

```

ak_num_threads -> How many apachekill threads to start, effect how strong the attack is  
ak_num_ips -> How many unique IPS to use  
attack_delay -> How long to wait before attack starts  
ak_duration -> How long the attack runs for  

The goal of tuning this config is to find out which of these parameters can result in a more even split of data. My experience so far is how long the attack runs determine how many attack flows get generated.

**Note** Since apachekill is a memory consumption attack a high value of ak_num_threads (>100) can fully utilize the RAM on weaker systems.


### Editing Benign Configs

In this example we will edit the MW_Locust Config located here `magicwand_components/benign/mw_locust.json `


```bash
{
  "client_options": {
    "stagger": "ON",                     <- Play around with
    "num_ips": 20,                       <- Play around with (My experience so far num_ips is biggest factor)
    "client_duration": 300,              <- Play around with
    "locust_duration": 300,              <- Play around with
    "wait_max": 60                       <- Play around with
    "keepalive": "ON"                    <- Play around with
    "traffic_behavior": "default"        <- See Traffic behavior section below

  },
  "benign": "mw_locust",
  "compose-file": "magicwand_components/benign/mw_locust.yml"
}

```

stagger -> Determines if clients should randomly start up or start all at once. Can either be "ON" or "OFF"  
num_ips -> Determines how many unique IPs to use  
client_duration -> How long each client is sending connections  
locust_duration -> How long the container is up for  
wait_max -> The max wait time of the random start up interval
keepalive -> Requests will send http keepalive flag can be ON,OFF, or RANDOM
traffic_behavior: The traffic pattern the clients will use

Benign clients don’t generate as much flows as attack traffic, so increasing the number of clients has the biggest effect on how many flows are generated.

### Traffic behavior

Our Locust container comes with a "smart profile" that can be configured to give clients a more "realistic" traffic pattern.

The pattern is configurable to act on a simulated "per hour basis" using our time compression algorithm which determines how much delay to have between requests. The algorithm works like this...    

1. Determine where in the behavior pattern our current time is
2. Select the behavior pattern value just before, and just after our current time (bounding pattern)
3. Interpolate the two bounding patterns values into approximate requests per hour
4. Interpolate requests per hour into seconds between requests
5. Linearly interpolate to determine our current seconds between requests
6. Scale seconds between requests by the ratio of seconds_per_hour and 3600

There are 8 possible behavior states for the clients 
 
```bash
"requests_per_hour": {
        "A": 1,
        "B": 2,
        "C": 4,
        "D": 8,
        "E": 16,
        "F": 32,
        "G": 64,
        "H": 128
}
```

#### Example traffic behavior patterns
Here are some example behavior patterns to try. Each client in the experiment will use the same pattern.

* Default Config
The default config is 
`AAAAAAAAAAAAAAAAAAAAAAAA`   

This means that at every "hour" the clients are sending a simulated 1 request per hour.

* Super Max Config
The "strongest" config would be all H for each of the 24 hours  
`HHHHHHHHHHHHHHHHHHHHHHHH`  


This means that at every "hour" the clients are sending a simulated 128 requests per hour.

* Mix and Match

Here is a config that is "sleeping" during off hours, and has lots of activity during "working hours"

`AAAAAAABHDEFGBAAAAAAAA`  

### Editing Experiment Duration

To make experiments longer or shorter, you need to change the run duration of the SUT(System Under Test)   
The config file is located here `magicwand_components/suts/mw_apache_wp.json`


```bash
{
  "pcap_snap_length": 400,
  "max_clients": 250,                              <-  Max number of concurrent clients
  "run_duration": 300,                             <-  Run duration 
  "sut": "mw_apache_wp",
  "compose-file": "magicwand_components/suts/mw_apache_wp.yml"
}
```

Once the SUT stops, all other containers will also shutdown

### Editing Component Resource Usage

If we wanted to limit how much resources a component uses, you can edit the component docker.yml file  

For example, we can edit the SUT be editing this file `magicwand_components/suts/mw_apache_wp.yml`

```
version: '2.2'
services:

  mw-sut-apachewp:
    image: twosixlabsmagicwand/mw-sut-apachewp:latest
    privileged: true
    ports:
      - "80:80"
    environment:
      - TEST_DURATION=${CURR_TEST_DURATION}
      - MAX_CLIENTS=${CURR_MAX_CLIENTS}
      - PCAP_PACKET_SNAPLEN=${CURR_PCAP_PACKET_SNAPLEN}
    volumes:
      - ./${CURR_RUN}:/home/
    #mem_limit: 16G     <- Uncomment and set how much memory to use
    #cpus: 8  <- Uncomment and set how many CPUs to use
```

### Editing Advice

While there are lots of editable features, for the purposes of a CSL dataset, Attack duration, Attack Strength, number of clients IPS, and client duration seem to have the biggest effect on how many flows are generated. Each Attack will have a sweet spot to play around with.

**5. Repeat Runs**  
Once everything is configured you can do another run, to see if you get the results you desired. if Not back to step 4# Calibration

Calibration attempts to find the parameters needed for the runs to exhibit the desired effects of the attack based on current hardware resources.

It does this by following this sequence.


!!! note 
    As of v0.4.0, only `apachekill` can be calibrated

1. Start a run with only client traffic
2. Start a run with only attack traffic
3. Compare various metrics such as RTT and SUT resources
4. If the ratio between the client and attack metrics is good, then save a tuned json file


Example Use
```bash
magicwand calibrate --attack apachekill
```
# MAGICWAND Data Generator

<!-- ![Magicwand Apachekill Run](images/mw-ak-calibrate.png) -->

## What is Magicwand

Magicwand is a platform to provide high-quality, reliable, and reproducible data sets for low-and-slow DDoS attacks. With the use of Docker images and customizable JSON files, users can generate a multitude of network traffic PCAPS.

Data is generated by orchestrating configurable test with both benign and attack traffic targeting a system under test (SUT) using Docker on the user's machine.

Resulting data includes PCAP captures, experiment metadata, the SUT resource/consumption metrics, and more; data is stored locally with machine learning model development in mind.

## Getting Started

This guide will walk through the basics of starting your first Magicwand experiment and producing your first PCAP.

Before you get started, you’re going to need a few things:


* [docker](https://docs.docker.com/get-docker/)
* [docker-compose](https://docs.docker.com/compose/install/)
* [python3.6+](https://www.python.org/downloads/)
* [tshark](https://www.wireshark.org/docs/man-pages/tshark.html)


### Hardware Requirements
Magicwand experiments consume a considerable amount of system resources, so we recommend that your machine has at least 8GB of RAM and 2 CPUS


## Installing Magicwand
Magicwand is compatible with Python 3.6 or later. The simplest way to install Magicwand and its dependencies is from PyPI with pip, Python's preferred package installer.

```bash
$ pip install magicwand
```
Note that Magicwand is an active project and routinely publishes new releases. In order to upgrade Magicwand to the latest version, use pip as follows

```bash
$ pip install -U magicwand
```
In addition to the python package, Magicwand leverages prebuilt docker images to run experiments


 You can build the images by running

```bash
$ bash images/build_all.sh
```


## Using Magicwand

The Magicwand CLI provides a streamlined methodology to generate reproducible network traffic data. This section will walk through how we setup an experiment and produce data.

### Initialize your Experiment

To begin we must first initialize your experiment directory. 

```bash
$ magicwand init --folder test_magicwand
$ cd test_magicwand
```

This will create a new folder called test_magicwand with all the required configuration/orchestration files needed to run experiments on your local machine.

To learn more visit the [CLI page](cli.md##initialize-project-init) for further details.


### Calibrate Your Machine

After our experiment directory is created, we can now calibrate the attacks to create the desired effect based on your hardware resources. This is an optional step, but it does allow us to verify that your machine can run magicwand experiments. To learn more about the calibration process and methodology view the [Calibration page](calibration.md) for further details.

To start simply run

```bash
$ magicwand calibrate --attack apachekill
```

This command will roughly take 5 minutes to complete. It will run two experiments one with benign traffic only and one with attack traffic only. The goal is to see if the apache server is being affected in the desired way, with only one type of traffic. Once complete the calibrate JSON will be saved in the tuned_jsons folder. 

This is an example output from calibration.

```bash
#Example Output
2020-10-28 17:29:27,197 | INFO : Running calibration for: apachekill
2020-10-28 17:29:27,201 - INFO -  mw-log - Starting runs for data version: mw_calibrate_runs
2020-10-28 17:29:27,204 - INFO -  mw-log - Running for 150 seconds
Creating network "suts_default" with the default driver
Creating suts_mw-sut-apachewp_1 ... done
Creating suts_mw-client-rtt-tracker_1 ... done
Creating suts_mw-client-smartswarm_1  ... done
Stopping suts_mw-client-smartswarm_1  ... done
Stopping suts_mw-client-rtt-tracker_1 ... done
Stopping suts_mw-sut-apachewp_1       ... done
Removing suts_mw-client-smartswarm_1  ... done
Removing suts_mw-client-rtt-tracker_1 ... done
Removing suts_mw-sut-apachewp_1       ... done
Removing network suts_default
Running cfm
cic.cs.unb.ca.ifm.Cmd You select: /home/tcpdump.pcap
cic.cs.unb.ca.ifm.Cmd Out folder: /home
cic.cs.unb.ca.ifm.Cmd CICFlowMeter received 1 pcap file
Working on... tcpdump.pcap
tcpdump.pcap is done. total 360 flows 
Packet stats: Total=8570,Valid=8545,Discarded=25
-------------------------------------------------------------------------------
Done
2020-10-28 17:32:30,338 - INFO -  mw-log - Benign Stats:
2020-10-28 17:32:30,339 - INFO -  mw-log - 358 flows out of 359 flows: 1.0
2020-10-28 17:32:30,339 - INFO -  mw-log - Attack Stats:
2020-10-28 17:32:30,339 - INFO -  mw-log - 0 flows out of 359 flows: 0.0
2020-10-28 17:32:30,361 - INFO -  mw-log - Finished Run 1 Successfully
2020-10-28 17:32:30,368 - INFO -  mw-log - Starting runs for data version: mw_calibrate_runs
2020-10-28 17:32:30,370 - INFO -  mw-log - Running for 150 seconds
Creating network "suts_default" with the default driver
Creating suts_mw-sut-apachewp_1 ... done
Creating suts_mw-attack-apachekill_1  ... done
Creating suts_mw-client-rtt-tracker_1 ... done
Stopping suts_mw-attack-apachekill_1  ... done
Stopping suts_mw-client-rtt-tracker_1 ... done
Stopping suts_mw-sut-apachewp_1       ... done
Removing suts_mw-attack-apachekill_1  ... done
Removing suts_mw-client-rtt-tracker_1 ... done
Removing suts_mw-sut-apachewp_1       ... done
Removing network suts_default
Running cfm
cic.cs.unb.ca.ifm.Cmd You select: /home/tcpdump.pcap
cic.cs.unb.ca.ifm.Cmd Out folder: /home
cic.cs.unb.ca.ifm.Cmd CICFlowMeter received 1 pcap file
Working on... tcpdump.pcap
tcpdump.pcap is done. total 4461 flows 
Packet stats: Total=28584,Valid=28559,Discarded=25
-------------------------------------------------------------------------------
Done
2020-10-28 17:35:43,252 - INFO -  mw-log - Benign Stats:
2020-10-28 17:35:43,253 - INFO -  mw-log - 234 flows out of 4460 flows: 0.05
2020-10-28 17:35:43,253 - INFO -  mw-log - Attack Stats:
2020-10-28 17:35:43,253 - INFO -  mw-log - 4225 flows out of 4460 flows: 0.95
2020-10-28 17:35:43,517 - INFO -  mw-log - Finished Run 1 Successfully
2020-10-28 17:35:43,555 - INFO -  mw-log - Checking results...
2020-10-28 17:35:43,555 - INFO -  mw-log - Calibration Completed Successfully
2020-10-28 17:35:43,562 - INFO -  mw-log - Calibration Process Finished
```

### Run Your Experiment

Now we are ready run out first experiment. To start simply invoke the ```run``` command as follows

```bash
$ magicwand run --config configs/mw_locust-apachekill.json --count 1 --data_version apachekill_runs 
```

This will start one run using the apachekill attack configuration and save the results to ```data_runs/apachekill_runs```
This command will roughly take 5 minutes to complete, here is some example output

```bash
2020-09-01 13:03:02,463 - INFO -  mw-log - Starting runs for data version: test_runs
2020-09-01 13:03:02,463 - INFO -  mw-log - Running for 300 seconds
Creating network "run_configs_default" with the default driver
Creating run_configs_mw-sut-apachewp_1 ... done
Creating run_configs_mw-attack-apachekill_1  ... done
Creating run_configs_mw-client-rtt-tracker_1 ... done
Creating run_configs_mw-client-smartswarm_1  ... done
Stopping run_configs_mw-client-rtt-tracker_1 ... done
Stopping run_configs_mw-client-smartswarm_1  ... done
Stopping run_configs_mw-attack-apachekill_1  ... done
Stopping run_configs_mw-sut-apachewp_1       ... done
Removing run_configs_mw-client-rtt-tracker_1 ... done
Removing run_configs_mw-client-smartswarm_1  ... done
Removing run_configs_mw-attack-apachekill_1  ... done
Removing run_configs_mw-sut-apachewp_1       ... done
Removing network run_configs_default
Running cfm
cic.cs.unb.ca.ifm.Cmd You select: /home/tcpdump.pcap
cic.cs.unb.ca.ifm.Cmd Out folder: /home
cic.cs.unb.ca.ifm.Cmd CICFlowMeter received 1 pcap file
Working on... tcpdump.pcap
tcpdump.pcap is done. total 4461 flows 
Packet stats: Total=28584,Valid=28559,Discarded=25
-------------------------------------------------------------------------------
Done
2020-10-28 17:35:43,252 - INFO -  mw-log - Benign Stats:
2020-10-28 17:35:43,253 - INFO -  mw-log - 234 flows out of 4460 flows: 0.05
2020-10-28 17:35:43,253 - INFO -  mw-log - Attack Stats:
2020-10-28 17:35:43,253 - INFO -  mw-log - 4225 flows out of 4460 flows: 0.95
2020-09-01 13:04:04,998 - INFO -  mw-log - Finished Run 1 Successfully
2020-09-01 13:05:26,552 - INFO -  mw-log - ALL IPs present
2020-09-01 13:06:36,705 - INFO -  mw-log - All verifications passed
2020-09-01 13:07:41,509 - INFO -  mw-log - Run Process Finished
```

Once finished we can view the output files in the folder ```data_runs/apachekill_run/*/```

```
$ ls data_runs/apachekill_run/*/
apache_stats.csv        ip_map_client.csv       mem_stats.csv           tcpdump.pcap
cic_flow_labeled.csv    ip_map_rtt.csv          rtt_stats.csv           tcpdump.pcap_Flow.csv
ip_attr_map.csv         ip_map_sut.csv          run_config.json         tcpdump_verify.csv
ip_map_attack.csv       locust_4ce7ab68e5df.csv run_parms.json          verify_run.json
```

To learn more about each file visit the [Data page](data.md#) for further details.
To learn more about the run command visit the [CLI page](cli.md#execute-run-run) for further details.


## Advanced Usage

The follow section will walk through running advanced experiments. 

### Changing Configurations

The configurations files are very flexible allowing you to change up the parameters to fit your use case. To edit the overall run, you will edit the JSON file in the configs folder. To edit each component, you edit the files in the magicwand_components directory.


This section will show you some example configurations to produce a variety of PCAPs. To see the full configuration options visit the [Configuration page](config.md) for further details.

#### No Attack Experiment

To generate a PCAP without attack traffic simply remove the attack parameter from the JSON file. 
The JSON file provided is an example to do a no attack experiment.

`configs/mw_locust-only.json`
```
{
    "benign": "mw_locust",
    "sut": "mw_apache_wp",
    "rtt": "mw_rtt_sensor",
    "run_type": "mw-locust-only"
}
```

#### 60 second Sockstress attack with staggered clients

This configuration will have Sockstress run for the first 60 seconds, and then stop, while the clients continue to ping server.

`configs/mw_locust-sockstress.json`
```
{
    "attack": "sockstress",
    "benign": "mw_locust",
    "sut": "mw_apache_wp",
    "rtt": "mw_rtt_sensor",
    "run_type": "mw-locust-sockstress"
}
```

`magicwand_components/benign/mw_locust.json`
```
{"client_options": {"stagger":"ON","num_ips": 20, "client_duration": 300, "locust_duration": 60, "wait_max": 60, "seed": "None"}, "benign": "mw_locust", "compose-file": "magicwand_components/benign/mw_locust.yml"}
```


`magicwand_components/attacks/sockstress.json`

```
{"attack_options": {"packet_delay": 5000, "attack_duration": 60, "attack_delay": 15}, "attack": "sockstress", "compose-file": "magicwand_components/attacks/sockstress.yml"}

```

#### Client Keepalive Random Experiment Slowread 

This configuration will have a Slowread attack, and locust clients that will randomly send HTTP keepalive on or off messages to the server. 


```configs/mw_locust-sht_slowread.json```
```
{
    "attack": "sht_slowread",
    "benign": "mw_locust",
    "sut": "mw_apache_wp",
    "rtt": "mw_rtt_sensor",
    "run_type": "mw-locust-sht_slowread"
}
```

`magicwand_components/benign/mw_locust.json`
```
{"client_options": {"stagger":"ON","num_ips": 20, "client_duration": 300, "locust_duration": 60, "wait_max": 60 "keepalive":"RANDOM"}, "benign": "mw_locust", "compose-file": "magicwand_components/benign/mw_locust.yml"}
```


`magicwand_components/attacks/sht_slowread.json`
```
{"attack_options": {"connections_per_sec": 256, "attack_duration": 300, "num_connections": 16384,"attack_delay": 15}, "attack": "sht_slowread", "compose-file": "magicwand_components/attacks/sht_slowread.yml"}
```
# Attacks

Magicwand uses "low volume" DDOS attacks in the data generation runs. The following are all the attacks included in the tool. 

## Apachekill

The ApacheKill attack abuses the HTTP protocol by requesting that the target web server return the requested URL content in a huge number of individual chunks, or byte ranges. It uses up a significant amount of memory on the target web server.

More info: <https://github.com/tkisason/KillApachePy>


## Sockstress

Sockstress is a Denial of Service attack on TCP services, it works by using RAW sockets to establish many TCP connections to a listening service. Because the connections are established using RAW sockets, connections are established without having to save any per-connection state on the attacker's machine.

Like SYN flooding, Sockstress is an asymmetric resource consumption attack: It requires very little resources (time, memory, and bandwidth) to run a Sockstress attack, but uses a lot of resources on the victim's machine. Because of this asymmetry, a weak attacker (e.g. one bot behind a cable modem) can bring down a rather large web server.

More info: <https://github.com/defuse/sockstress>

 
## Goloris

Goloris works by occupying and keeping many TCP connections open to the victim. It uses low network bandwidth to prolong the attack. The goal is to eventually consume all available connections to the victim, so no other client could connect to it.

More info: https://github.com/valyala/goloris

## SlowHTTPTest

SlowHTTPTest is a highly configurable tool that simulates some Application Layer Denial of Service attacks by prolonging HTTP connections in different ways.

It implements most common low-bandwidth Application Layer DoS attacks, such as Slowloris, Slow HTTP POST, Slow Read attack (based on TCP persist timer exploit) by draining concurrent connections pool, as well as Apache Range Header attack by causing very significant memory and CPU usage on the server.

Slowloris and Slow HTTP POST DoS attacks rely on the fact that the HTTP protocol, by design, requires requests to be completely received by the server before they are processed. If an HTTP request is not complete, or if the transfer rate is very low, the server keeps its resources busy waiting for the rest of the data. If the server keeps too many resources busy, this creates a denial of service. This tool is sending partial HTTP requests, trying to get denial of service from target HTTP server.

More info: <https://tools.kali.org/stress-testing/slowhttptest>
# Configuration

The following explains the fields mean for each Magicwand component JSONs, and how you can customize it for your experiments.


## Experiment 

Each Magicwand experiment is started from a base config in the configs folder. Experiments can have an attack, benign or both components set.

* attack - Attack component to use (optional)
* benign - Benign component to use (optional)
* sut - SUT component to use (required)
* rtt - rtt component to use (required)
* run_type - Name prefix for run folder (required)


### Example Config

```
{
    "attack": "apachekill",
    "benign": "mw_locust",
    "sut": "mw_apache_wp",
    "rtt": "mw_rtt_sensor",
    "run_type": "mw-locust-apachekill"
}
```

## Benign

### `mw_locust.json`

This is our specialized benign locust client

* benign - Name of benign client
    - client_options - Options for client image
        * num_ips - Number of unique IPs to use
        * client_duration - How long image should send traffic
        * locust_duration - How long individual locusts send traffic
        * wait_max - Maximum time to wait before starting a client
        * stagger - Options to stagger clients
            * ON - Will randomly stagger between 0 and WAIT_MAX seconds each client
            * OFF - No stagger, will use default locust settings
        * keepalive - Options to set http keepalive for clients
            * ON - http keepalive is on, this is default behavior 
            * OFF - http keepalive is off
            * RANDOM - http keepalive will randomly be on or off for each client request
        * seed - Unsigned 32 bit int used as a seed for all python random module generations (None if don’t want to set)
        * traffic_behavior -  24 character string of either ABCDEFGH representing traffic intensity per hour on a scale from 1 - 128  (default if don’t want to change)

   * compose-file - The docker compose file for this component

#### Example Config

``` 
{
  "client_options": {
    "stagger": "ON",
    "num_ips": 20,    
    "client_duration": 300,
    "locust_duration": 300,
    "wait_max": 60              
    "keepalive": "ON"          
    "traffic_behavior": "default"
    "seed": "None"
  },
  "benign": "mw_locust",
  "compose-file": "magicwand_components/benign/mw_locust.yml"
}

```

## Attacks

### `apachekill.json`

This is the apachekill attack client

* attack - Name of attack
    - attack_options - Options for attack image
        * ak_num_threads - Number of threads to spawn
        * ak_num_ips - How many IPs to spawn
        * attack_delay - How long to wait before attack starts
        * ak_duration: How long the attack will go on for
   * compose-file - The docker compose file for this component

#### Example Config

``` 
{"attack_options": {"ak_num_threads": 50, "ak_num_ips": 20, "attack_delay": 15, "ak_duration": 300}, "attack": "apachekill", "compose-file": "magicwand_components/attacks/apachekill.yml"}
```

### `sockstress.json`

This is the Sockstress attack client

* attack - Name of attack
    - attack_options - Options for attack image
        * packet_delay: How long to wait before requests (lower is stronger)
        * attack_duration: How long the attack will go on for 
        * attack_delay - How long to wait before attack starts
   * compose-file - The docker compose file for this component

#### Example Config

``` 
{"attack_options": {"packet_delay": 5000, "attack_duration": 300, "attack_delay": 15}, "attack": "sockstress", "compose-file": "magicwand_components/attacks/sockstress.yml"}
```




### `goloris.json`

This is the Goloris attack client

* attack - Name of attack
    - attack_options - Options for attack image
        * worker_count: How many processes to spawn,
        * ramp_up_interval: How quickly processes spawn
        * attack_duration: How long the attack will go on for
        * attack_delay - How long to wait before attack starts
   * compose-file - The docker compose file for this component

#### Example Config

``` 
{"attack_options": {"worker_count": 16, "ramp_up_interval": 1, "attack_delay": 15, "attack_duration": 300}, "attack": "goloris", "compose-file": "magicwand_components/attacks/goloris.yml"}
```


### `sht_rudeadyet.json`

This is the sht_rudeadyet attack client

* attack - Name of attack
    - attack_options - Options for attack image
        * connections_per_sec: How many connections to spawn per second
        * attack_duration: How long the attack will go on for
        * num_connections": Number of connections to make
        * attack_delay - How long to wait before attack starts
   * compose-file - The docker compose file for this component

#### Example Config

``` 
{"attack_options": {"connections_per_sec": 256, "attack_duration": 300, "num_connections": 16384,"attack_delay": 15}, "attack": "sht_rudeadyet", "compose-file": "magicwand_components/attacks/sht_rudeadyet.yml"}
```


### `sht_slowloris.json`

This is the sht_slowloris attack client

* attack - Name of attack
    - attack_options - Options for attack image
        * connections_per_sec: How many connections to spawn per second
        * attack_duration: How long the attack will go on for
        * num_connections": Number of connections to make
        * attack_delay - How long to wait before attack starts
   * compose-file - The docker compose file for this component

#### Example Config

``` 
{"attack_options": {"connections_per_sec": 256, "attack_duration": 300, "num_connections": 16384,"attack_delay": 15}, "attack": "sht_rudeadyet", "compose-file": "magicwand_components/attacks/sht_slowloris.yml"}
```

### `sht_slowread.json`

This is the sht_slowread attack client

* attack - Name of attack
    - attack_options - Options for attack image
        * connections_per_sec: How many connections to spawn per second
        * attack_duration: How long the attack will go on for
        * num_connections": Number of connections to make
        * attack_delay - How long to wait before attack starts
   * compose-file - The docker compose file for this component

#### Example Config

``` 
{"attack_options": {"connections_per_sec": 256, "attack_duration": 300, "num_connections": 16384,"attack_delay": 15}, "attack": "sht_rudeadyet", "compose-file": "magicwand_components/attacks/sht_slowread.yml"}
```



## SUT

### `mw_apache_wp.json`

This is the mw_apache_wp SUT

* sut - Name of SUT
* max_clients - Maximum number of httpd processes that can spawn
* run_duration - How long Magicwand run is
* pcap_snap_length - How many bytes to capture per packet from tcpdump
* compose-file - The docker compose file for this component

#### Example Config

```
{"pcap_snap_length": 400 , "max_clients": 250,  "run_duration": 300, "sut": "mw_apache_wp", "compose-file": "magicwand_components/suts/mw_apache_wp.yml"}
```


## Sensors


### `mw_rtt_sensor.json`

This is the rtt sensor

* rtt - Name of rtt sensor
* compose-file - The docker compose file for this component
* timeout - How many seconds for request to wait before timeout


#### Example Config

```
{"rtt":"mw_rtt_sensor", "compose-file": "magicwand_components/sensors/mw_rtt_sensor.yml", "timeout": 2}
```
## Data

Each Magicwand experiment will generate data that can be used for modeling purposes to identify attack and client traffic.
This document will go over how each of these files is generated. How you use them is up to you

```
apache_stats.csv        ip_map_client.csv       mem_stats.csv           tcpdump.pcap
cic_flow_labeled.csv    ip_map_rtt.csv          rtt_stats.csv           tcpdump.pcap_Flow.csv
ip_attr_map.csv         ip_map_sut.csv          run_config.json         tcpdump_verify.csv
ip_map_attack.csv       locust_4ce7ab68e5df.csv run_parms.json          verify_run.json   
```

### tcpdump.pcap

This is the PCAP generated for each experiment. PCAPs hold the raw network traffic data and are the basis for most network security datasets.

This is the command used to collect the PCAP for each experiment.
```
tcpdump -vvv -i any -s $PCAP_SNAP_LEN -w /home/tcpdump.pcap
```

### tcpdump_verify.csv


To verify the data for each experiment, we used tshark to convert the raw PCAP into a CSV. With this CSV we can do some data validation steps to ensure that the data is correct. This is the tshark command run automatically post-experiment to convert the PCAP to a CSV.

```
tshark -r tcpdump.pcap -i wlan1 -T fields -E header=y -E separator=, -E quote=d -e _ws.col.No. -e _ws.col.Time -e _ws.col.Source -e _ws.col.Destination -e _ws.col.Protocol -e _ws.col.Length -e _ws.col.Info  -e http.user_agent -e http.connection -e http.request -e http.request.line > tcpdump_verify.csv
```
To learn more about how you can use tshark to convert the PCAP for your needs [click here](https://tshark.dev/)  



This is a glimpse of what the tcpdump_verify.csv file looks like.
```
_ws.col.No.,_ws.col.Time,_ws.col.Source,_ws.col.Destination,_ws.col.Protocol,_ws.col.Length,_ws.col.Info,http.user_agent,http.connection,http.request,http.request.line
"1","0.000000","192.168.32.4","192.168.32.2","TCP","74","39684 → 80 [SYN] Seq=0 Win=1152 Len=0 MSS=1460 SACK_PERM=1 TSval=1823557209 TSecr=0 WS=1",,,,
"2","0.001309","192.168.32.4","192.168.32.2","TCP","74","39246 → 80 [SYN] Seq=0 Win=1152 Len=0 MSS=1460 SACK_PERM=1 TSval=1823557211 TSecr=0 WS=1",,,,
"3","0.001316","192.168.32.4","192.168.32.2","TCP","74","39244 → 80 [SYN] Seq=0 Win=1152 Len=0 MSS=1460 SACK_PERM=1 TSval=1823557211 TSecr=0 WS=1",,,,
```


### verify_run.json

This JSON file depicts the output of the experiment verifier.    
Each component has its own verify function and will return true or false if the experiment checks passed.

```
{
  "MW_Apache_WP": true,
  "Apachekill": true,
  "MW_Locust": true,
  "MW_RTT_Sensor": true,
  "MW_Global": true
}
```


### rtt_stats.csv

This is a CSV file produced by the Round-Trip Time (RTT) tracker that collects how long each request took. This file can be used to measure how responsive the SUT was during each experiment.

* timestamp - Seconds since experiment started
* rtt - Round Trip time in seconds for a request to be sent

```
timestamp,rtt
0,0.007287
1,0.002221
2,0.003041
3,0.001833
4,0.002691
```

### ip_attr_map.csv

This is an IP attribution csv that keeps track of what each IP is i.e. is it an attack or benign client. This file is an aggregate of the following files

```
ip_map_sut.csv ip_map_attack.csv ip_map_client.csv ip_map_rtt.csv
```

The following explains what each field means.  

* index - Number in csv
* ip - IP address in the experiment
* type - Type of traffic attack, client or from the SUT
* subtype - More granular view of IP such as attack name 

```
index,ip,type,subtype
0,192.168.1.1,attack,apachekill
1,192.168.1.2,attack,apachekill
2,192.168.1.3,attack,apachekill
3,192.168.1.4,attack,apachekill
4,192.168.1.5,attack,apachekill
5,192.168.2.1,attack,apachekill
6,192.168.2.2,attack,apachekill
7,192.168.2.3,attack,apachekill
8,192.168.2.4,attack,apachekill
```


### locust_HOSTNAME.csv

This is csv is produced by locust. It catalogs the requests and fails at each timestamp
The following explains what each field means.  


* timestamp -local time for image
* requests - Number of successful connections
* fails - Number of failed connections

```
timestamp,requests,fails
14:47:33,0,0
14:47:34,0.0,0
14:47:35,0.0,0
14:47:36,0.0,0
14:47:37,2.0,0
14:47:38,2.6666666666666665,0
```

### apache_stats.csv 

This csv is produced by apache. it catalogs various system metrics during the run.
The following explains what each field means.  


* timestamp - Seconds since experiment started
* total_access - Total number of requests 
* total_traffic - Size of traffic transferred to clients 
* cpu_load - CPU usage of Apache image
* requests per second - Number of requests per second
* httpd mem used(kb) - Memory used from httpd processes
* num_httpd - Number of httpd processes spawned
* system memory free - Total system memory free

```
timestamp,total_access,total_traffic,cpu_load,requests per second,httpd mem used(kb),num_httpd,system memory free
1.0,94,369 kB,.2%,18.8,22728,9,3499
2.0,99,375 kB,.167%,16.5,28528,11,3498
3.0,105,382 kB,.286%,15,40088,15,3497
4.0,112,389 kB,.375%,14,63168,23,3495
5.0,120,398 kB,.333%,13.3,63112,23,3494
6.0,129,407 kB,.3%,12.9,63140,23,3494
7.0,139,418 kB,.273%,12.6,63112,23,3493
9.0,150,430 kB,.25%,12.5,66116,25,3492
10.0,162,442 kB,.231%,12.5,72000,26,3492

```

### mem_stats.csv

This csv shows system memory being used by the SUT as reported from docker stats.
The following explains what each field means.  

* timestamp - Seconds since experiment started
* memory_percent - Memory usage reported by docker stats

```
timestamp,memory_percent
3,2.21
10,2.36
17,7.12
25,15.24
32,17.63
40,25.41
48,27.44
55,29.6
63,34.9
```

### run_config.json

This is the run config used to start the experiment.  
For details on the component fields view the [configuration page](config.md#) for details on what each field means.  
```
{
    "attack": "apachekill",
    "benign": "mw_locust",
    "sut": "mw_apache_wp",
    "rtt": "mw_rtt_sensor",
    "run_type": "mw-locust-apachekill"
}
```


### run_params.json

This is the JSON file that was used to configure the experiment.   

These are the run metadata fields.  
* "compose_file_string": docker compose file string
* "run_loc": Location of run data
* "version": Version of magicwand used to generate data


For details on the component fields view the [configuration page](config.md#) for details on what each field means.


```
{
  "sut": {
    "pcap_snap_length": 400,
    "max_clients": 250,
    "run_duration": 300,
    "sut": "mw_apache_wp",
    "compose-file": "magicwand_components/suts/mw_apache_wp.yml"
  },
  "attack": {
    "attack_options": {
      "ak_num_threads": 50,
      "ak_num_ips": 20,
      "attack_delay": 15,
      "ak_duration": 300
    },
    "attack": "apachekill",
    "compose-file": "magicwand_components/attacks/apachekill.yml"
  },
  "benign": {
    "client_options": {
      "stagger": "ON",
      "num_ips": 20,
      "client_duration": 300,
      "locust_duration": 300,
      "wait_max": 10,
      "keepalive": "ON",
      "traffic_behavior": "HHHHHHHHHHHHHHHHHHHHHHHH"
    },
    "benign": "mw_locust",
    "compose-file": "magicwand_components/benign/mw_locust.yml"
  },
  "rtt": {
    "rtt": "mw_rtt_sensor",
    "compose-file": "magicwand_components/sensors/mw_rtt_sensor.yml",
    "timeout": 2
  },
  "compose_file_string": "-f magicwand_components/suts/mw_apache_wp.yml -f magicwand_components/attacks/apachekill.yml -f magicwand_components/benign/mw_locust.yml -f magicwand_components/sensors/mw_rtt_sensor.yml ",
  "run_loc": "magicwand_components/suts/runs/mw-locust-apachekill_10_29_2020T10_37_57Z/",
  "version": "beta-1.0.1"
}
```


### cic_flow_labeled.csv

This is a labeled csv has 80 statistical network traffic features such as Duration, Number of packets, Number of bytes, and Length of packets for each bidirectional flow in the PCAP. 

The labels are `attack` for attack clients and `client` for benign clients.

!!! note 
    tcpdump.pcap_Flow.csv is the unlabeled version of the csv

To learn more about the features visit [here](https://github.com/CanadianInstituteForCybersecurity/CICFlowMeter/blob/master/ReadMe.txt#L80)

Here is an example output
```
,Flow ID,Src IP,Src Port,Dst IP,Dst Port,Protocol,Timestamp,Flow Duration,Total Fwd Packet,Total Bwd packets,Total Length of Fwd Packet,Total Length of Bwd Packet,Fwd Packet Length Max,Fwd Packet Length Min,Fwd Packet Length Mean,Fwd Packet Length Std,Bwd Packet Length Max,Bwd Packet Length Min,Bwd Packet Length Mean,Bwd Packet Length Std,Flow Bytes/s,Flow Packets/s,Flow IAT Mean,Flow IAT Std,Flow IAT Max,Flow IAT Min,Fwd IAT Total,Fwd IAT Mean,Fwd IAT Std,Fwd IAT Max,Fwd IAT Min,Bwd IAT Total,Bwd IAT Mean,Bwd IAT Std,Bwd IAT Max,Bwd IAT Min,Fwd PSH Flags,Bwd PSH Flags,Fwd URG Flags,Bwd URG Flags,Fwd Header Length,Bwd Header Length,Fwd Packets/s,Bwd Packets/s,Packet Length Min,Packet Length Max,Packet Length Mean,Packet Length Std,Packet Length Variance,FIN Flag Count,SYN Flag Count,RST Flag Count,PSH Flag Count,ACK Flag Count,URG Flag Count,CWR Flag Count,ECE Flag Count,Down/Up Ratio,Average Packet Size,Fwd Segment Size Avg,Bwd Segment Size Avg,Fwd Bytes/Bulk Avg,Fwd Packet/Bulk Avg,Fwd Bulk Rate Avg,Bwd Bytes/Bulk Avg,Bwd Packet/Bulk Avg,Bwd Bulk Rate Avg,Subflow Fwd Packets,Subflow Fwd Bytes,Subflow Bwd Packets,Subflow Bwd Bytes,FWD Init Win Bytes,Bwd Init Win Bytes,Fwd Act Data Pkts,Fwd Seg Size Min,Active Mean,Active Std,Active Max,Active Min,Idle Mean,Idle Std,Idle Max,Idle Min,Label
0,172.18.0.4-172.18.0.2-37402-80-6,172.18.0.4,37402,172.18.0.2,80,6,18/09/2020 02:58:46 PM,1151,5,3,146.0,334.0,146.0,0.0,29.2,65.29318494299386,334.0,0.0,111.33333333333331,192.83498990933498,417028.6707211121,6950.4778453518675,164.42857142857144,269.1758143808969,748.0,12.0,1151.0,287.75,327.5956603294169,748.0,41.0,337.0,168.5,99.7020561473032,239.0,98.0,0,0,0,0,168,104,4344.0486533449175,2606.429192006951,0.0,334.0,53.33333333333334,115.7972365818805,13409.0,1,2,0,2,7,0,0,0,0.0,60.0,29.2,111.33333333333331,0,0,0,0,0,0,0,18,0,41,29200,235,1,32,0.0,0.0,0.0,0.0,1600441126433346.0,0.0,1600441126433346.0,1600441126433346.0,client
```
# Magicwand Images

These are the docker images needed for magicwand
You can build them with ```build_all.sh``` or pull from docker-hub with ```pull_images.sh```

## SUTS

System Under Test (SUT) is our target server that will be receiving traffic from the clients and attacks.

### Apache SUT

This is a clean Ubuntu 18.04 image with Apache/2.2.11
This version of Apache is vulnerable to attacks.

### Clean WP

This image is built from the Apache SUT and installs WordPress 4.7.3

## Clients

These are images that send non attack data to the SUT

### RTT-Tracker

This is our Round Trip Time (RTT) tracker that will send a request every second to the SUT. It has a timeout of 2 seconds.

### Smart-Swarm

This image uses [Locust](https://locust.io/) to simulate normal traffic sent to the SUT. We have a custom profile that simulates a "9-5 workweek" per 5 minutes, and configuration options to change the behavior of http keepalive.

## Attacks

These images are built to implement the DDOS attacks
Please see the [Attacks page](attacks.md#license) for further details.

## CIC

These images are for the using the CIC (Canadian Institute for Cybersecurity) PCAP to csv pipeline

### Converter 

This image converts PCAPs into a csv with 80 statistical network traffic features such as Duration, Number of packets, Number of bytes, Length of packets, etc.
To learn more visit [CICFlowMeter](https://www.unb.ca/cic/research/applications.html#CICFlowMeter)
# MAGICWAND Data Generator CLI

This document will go over the functions in the Magicwand CLI 

```buildoutcfg
Usage: magicwand [OPTIONS] COMMAND [ARGS]...

  Magicwand data generation tool

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  calibrate  Calibrate the magicwand tool
  convert    Convert any pcap file to a CIC csv file
  init       Create and initialize magicwand folder
  run        Start a run to generate data
```

## Initialize Project (`init`)

The init command creates a magicwand project

```buildoutcfg
Usage: magicwand init [OPTIONS]

  Create magicwand project

Options:
  --project TEXT Folder to create  [required]
  --help         Show this message and exit.
```

Example use
```buildoutcfg
magicwand init --project test
```

## Execute Runs (`run`)

The run command starts the experiment using the JSON file for N number of runs saved to a folder inside data_runs.

```buildoutcfg
Usage: magicwand run [OPTIONS]

  Start a run to generate data

Options:
  --config TEXT        JSON file with run parameters  [required]
  --data_version TEXT  Folder to save runs  [required]
  --count INTEGER      Number of runs
  -v, --verbose        Print more output.
  --help               Show this message and exit.
```

Example use
```bash
magicwand run --config configs/mw_locust-apachekill.json --count 5 --data_version test_runs 
```

## Calibrate Runs (`calibrate`)

The calibrate command tunes attacks to create the desired effects based on hardware resources.

```buildoutcfg
Usage: magicwand calibrate [OPTIONS]

  Calibrate Magicwand

Options:
  --ratio TEXT   Ratio config file
  --attack TEXT  Attack to tune, valid options: apachekill
                 [required]
  --help         Show this message and exit.
```

Example use
```bash
magicwand calibrate --attack apachekill
```


## Convert PCAPs (`convert`)

The convert command turns PCAP files to a CIC CSV file


```buildoutcfg
Usage: magicwand convert [OPTIONS]

  Convert any pcap file to a CIC csv file

Options:
  --pcap TEXT        PCAP to convert  [required]
  -o, --output TEXT  Output .csv filename
  -f, --force        Force output file to be overwritten if it already exists.
  --help             Show this message and exit.
```

Example use
```bash
magicwand convert --pcap my_pcap.pcap
```
# Data Validation

The key purpose of this tool is to provide high-quality, reliable, and reproducible data sets for low-and-slow DDoS attacks. To validate each data generation run, our tool conducts a series of checks. While we have been systematic and methodical in crafting these checks to make them as comprehensive as possible, we imagine that we weren't totally exhaustive.

## Motivating Questions

Through the data validation process, we aim to address the following questions:

1. Is the attack functioning as intended?
2. Are all the components of the experiment run functioning?
3. Are all the collection sensors functioning?
4. Did the run complete? Did we capture all of the expected data?

### Is the attack functioning as intended?

Our primary usage for this dataset was in the development of a machine-learning-based low-and-slow DDoS defense. Each of the attacks we have included could be converted into a volumetric (connection- or packet-based) attack, simply by ramping up the parameters of the attack (e.g. number of threads to use). If so, the generated data would not properly represent a low-and-slow attack.

<!-- TODO: add in discussion of how we approach this --> 

!!! note 
    There are certainly applications of these data sets where the volumetric attacks traits are exhibited. To achieve this behavior, you will need to use a custom [run configuration](data.md#run_tunedjson).

Here are some examples of attack signatures (Different for each attack) we check for during each experiment:

- Did attack work as expected?
    - apachekill
        -  is the attack client sending overlapped bytes?
    - sockstress: is attack client sending RSTs?
        - is the attack client sending RSTs?

- Did we see the expected data sent/received

### Are all of the components of the experiment run functioning?

For each experiment we check for evidence that each of these components is functioning.

!!! note 
    These checks have **NOT** be implemented in 1.0.1


Client Signatures:

- Did the benign client work as expected?
- Did the client do what we configured it to do
- Were the settings, right?

SUT Signatures:

- Did the SUT work as expected?
- Were there any unexpected incoming/outgoing connections?
- Did the SUT appear to send back valid HTTP data?

### Are all the sensors functioning correctly?

For each experiment we check for evidence that our sensor was up and running.

Evidence of sensors functioning

- Is a PCAP generated?

Resources:

- Did we RTT Stats?
- Did we capture Apache server stats (CPU,Memory,etc..)?

### Did the run complete? Did we capture all the expected data?

For each experiment we check for evidence that the run completed and produced data

Test Run: 

- Did we see all expected IPs in the PCAP
- are there IPs in the PCAP that aren't in ip_attr_map.csv?

Higher Level Metrics:

- RTTs as expected

Run Data:

- Were all files generated?
- Was everything captured as expected and not corrupted?
- was there anything unexpected in the data?

<!-- Tune using volumetric attacks?

- Connection-volumetric
- Packet-volumetric
- Data-volumetric (not considered for now) -->
