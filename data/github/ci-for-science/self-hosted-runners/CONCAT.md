# [0.1.0] 10/06/2020

## New
* Add recipe to add a Runner using Ubuntu.
# Self-hosted runners

| Five recommendations for fair software from [fair-software.nl](https://fair-software.nl) | Badges |
| --- | --- |
| 1. Code repository | [![GitHub badge](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/ci-for-science/self-hosted-runners/) |
| 2. License | [![License badge](https://img.shields.io/github/license/ci-for-science/self-hosted-runners)](https://github.com/ci-for-science/self-hosted-runners/) |
| 3. Community registry | [![Ansible Galaxy badge](https://img.shields.io/badge/galaxy-fixme.fixme-660198.svg)](https://galaxy.ansible.com/fixme/fixme) [![Research Software Directory](https://img.shields.io/badge/rsd-self--hosted--runners-00a3e3.svg)](https://www.research-software.nl/software/self-hosted-runners) |
| 4. Enable citation | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3904265.svg)](https://doi.org/10.5281/zenodo.3904265) |
| 5. Checklist | N/A |
| **Other best practices** | |
| Markdown Link Checker| [![Check Markdown links](https://github.com/ci-for-research/self-hosted-runners/workflows/Check%20Markdown%20links/badge.svg)](https://github.com/ci-for-research/self-hosted-runners/actions?query=workflow%3A%22Check+Markdown+links%22) |

## How to set up GitHub Action runners on self-hosted infrastructure

This repository explains how to set up a server for running continuous integration tests on other hardware than what
GitHub provides. This can be useful when the code you want to test has special requirements, for example if

- it needs a GPU to run ([CUDA installation with ansible instructions](/install-cuda/README.md))
- it needs multiple nodes
- testing requires data that needs to stay on-premises for privacy reasons or legal reasons
- testing requires data that is too big to move
- testing requires specific software

This guide distinguishes between the _client_ and the _server_; the client is your own machine; the server is whichever
machine runs the tests. For either side, we'll explain what configuration needs to be done. For people who just want to
try out the instructions but don't have access to remote hardware, we included a few alternatives for running the server
locally as well, through the use of virtualization (with VirtualBox) and containerization (with Docker).

For the client, we included instructions for Linux Ubuntu, Mac, and Windows; the server-side instructions all assume
Linux Ubuntu.

| Status | Client OS | Server hardware | Runner |
| --- | --- | --- | --- |
| :heavy_check_mark: Completed | Linux Ubuntu | local machine via Docker           | [link](/ubuntu-docker/README.md)          |
| :heavy_check_mark: Completed | Linux Ubuntu | local machine via Singularity      | [link](/ubuntu-singularity/README.md)     |
| :heavy_check_mark: Completed | Linux Ubuntu | local machine via Vagrant          | [link](/ubuntu-vagrant/README.md)         |
| :heavy_check_mark: Completed | Linux Ubuntu | local machine via VirtualBox       | [link](/ubuntu-virtualbox/README.md)      |
| :heavy_check_mark: Completed | Linux Ubuntu | remote machine at [SURF HPC Cloud] | [link](/ubuntu-surf-hpc-cloud/README.md)  |
| :hourglass_flowing_sand: WIP | Mac          | local machine via Docker           | -                                         |
| :hourglass_flowing_sand: WIP | Mac          | local machine via Vagrant          | -                                         |
| :hourglass_flowing_sand: WIP | Mac          | local machine via VirtualBox       | -                                         |
| :hourglass_flowing_sand: WIP | Mac          | remote machine at [SURF HPC Cloud] | -                                         |
| :hourglass_flowing_sand: WIP | Windows      | local machine via Docker           | -                                         |
| :heavy_check_mark: Completed | Windows      | local machine via Vagrant          | [link](windows-vagrant/README.md)         |
| :hourglass_flowing_sand: WIP | Windows      | local machine via VirtualBox       | -                                         |
| :heavy_check_mark: Completed | Windows      | remote machine at [SURF HPC Cloud] | [link](/windows-surf-hpc-cloud/README.md) |

## Security

**A warning from GitHub for self-hosted runners in combination with public repositories is shown
[here](https://help.github.com/en/actions/hosting-your-own-runners/about-self-hosted-runners#self-hosted-runner-security-with-public-repositories).
Please take this seriously. It basically means that the combination of a self-hosted runner and a public GitHub
repository is unsafe. However, there was a [recent discussion](https://github.com/actions/runner/issues/494) indicating
that GitHub may add features to make this combination safe in the near future.**

[SURF HPC Cloud]: https://userinfo.surfsara.nl/systems/hpc-cloud

## Documentation for developers


If you want to check if the links in your markdown work, install markdown-link-check

```shell
npm install
```

then run

```shell
find . -name '*.md' -not -path './node_modules/*' -exec markdown-link-check '{}' --config .mlc-config.json ';'
```

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
available at [https://www.contributor-covenant.org/version/1/4/code-of-conduct.html](https://www.contributor-covenant.org/version/1/4/code-of-conduct.html)

For answers to common questions about this code of conduct, see
[https://www.contributor-covenant.org/faq](https://www.contributor-covenant.org/faq)
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/). Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/ci-for-research/self-hosted-runners/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/ci-for-research/self-hosted-runners/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community _before you start working_. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. make sure the existing tests still work by running ``python setup.py test``;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. [push](http://rogerdudler.github.io/git-guide/) your feature branch to (your fork of) the Fortran_Davidson repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
# Setting up a CI server for a GitHub Action runner with Docker from Linux Ubuntu

After following this guide, you'll have a simple GitHub action workflow on a GitHub repository of your choice. When new commits are made to your repository, the workflow delegates work to a server which runs in a [Docker](https://www.docker.com/) container.

This guide distinguishes between the _client_ and the _server_; the client is your own machine; the server is whichever
machine will run the tests. This document describes the case where the server is a Docker container running on your own machine.

For guides on how to configure other features in addition to just the runner, go [here](/README.md).

## Prerequisites

1. Install Docker: [https://docs.docker.com/engine/install/](https://docs.docker.com/engine/install/)
2. Follow post-installation steps [https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user) to manage docker as a non-root user

### Testing your Docker setup

Refence: [https://docs.docker.com/docker-for-windows/#test-your-installation](https://docs.docker.com/docker-for-windows/#test-your-installation)

1. Open a terminal window

2. Run ``docker --version`` to ensure that you have a supported version of Docker:

```shell
> docker --version

Docker version 19.03.11-ce, build 42e35e61f3
```

3. Pull the hello-world image from Docker Hub and run a container:

```shell
> docker run hello-world

Unable to find image 'hello-world:latest' locally
latest: Pulling from library/hello-world
0e03bdcc26d7: Pull complete
Digest: sha256:d58e752213a51785838f9eed2b7a498ffa1cb3aa7f946dda11af39286c3db9a9
Status: Downloaded newer image for hello-world:latest

Hello from Docker!
This message shows that your installation appears to be working correctly.
...
```

4. List the hello-world image that was downloaded from Docker Hub:

```shell
> docker image ls

REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
hello-world         latest              bf756fb1ae65        6 months ago        13.3kB
```

5. List the hello-world container (that exited after displaying “Hello from Docker!”):

```shell
> docker container ls --all

CONTAINER ID        IMAGE               COMMAND             CREATED              STATUS                          PORTS               NAMES
1d624a063f22        hello-world         "/hello"            About a minute ago   Exited (0) About a minute ago                       flamboyant_ramanujan
```

## Server side configuration

### Build image

Now we are ready to build our Docker image. The following command will use [Dockerfile](Dockerfile) to build the image. It will create a system user, install necessary system packages and dependencies for the runner.

```shell
docker build \
    --tag github-actions-runner \
    .
```

The Docker image will create a user called ``tester`` with a password ``password``. If you want to change the default username and the password, you will need to adjust `<username>` and `<user password>` for the user which will be added to the Docker image.

```shell
docker build \
    --tag github-actions-runner \
    --build-arg DOCKER_USER="<username>" \
    --build-arg DOCKER_PASS="<user password>" \
    .
```

## Client side configuration

### Generate an OAuth token

We're almost ready to use our Docker image to set up a GitHub Runner, but first we need to
generate an OAuth token, as follows:

1. Go to [https://github.com/settings/tokens](https://github.com/settings/tokens) and click the ``Generate new token`` button.
2. Provide your GitHub password when prompted
3. Fill in a description for the token, for example _GitHub runner for github.com/&lt;your organization&gt;/&lt;your repository&gt;_
4. Enable the ``repo`` scope and all of its checkboxes, like so:

    ![Token permissions](/images/token_permissions.png)

5. Click ``Generate`` at the bottom. Make sure to copy its value because we'll need it in the next step

### Run the server

#### Daemon mode

The command below will start the Docker container in daemon mode. The Docker container will run in the background and the terminal will be available. To stop the running container you need to run the ``stop`` command, which is explained in the [Cleanup](#cleanup) section.

```shell
docker run -d --restart always --name github-actions-runner \
    --env PERSONAL_ACCESS_TOKEN=<Github OAuth token> \
    --env RUNNER_NAME=<runner name to appear on Github> \
    --env RUNNER_WORKDIR=/tmp/actions-runner-repo \
    --env GITHUB_ORG=<organization or username> \
    --env GITHUB_REPO=<name of the repository> \
    github-actions-runner:latest
```

If you stop the running container, the Github actions runner will also stop.

To stop the running Docker container:

```shell
docker container stop github-actions-runner
```

To start the Docker container again:

```shell
docker container start github-actions-runner
```

#### Temporary mode

The command below will run the docker image and set up the runner. When user presses `CTRL+C`, it automatically removes the runner from GitHub and removes the Docker container as well.

```shell
docker run --rm --name github-actions-runner \
    --env PERSONAL_ACCESS_TOKEN=<personal access token> \
    --env RUNNER_NAME=<runner name to appear on Github> \
    --env RUNNER_WORKDIR=/tmp/actions-runner-repo \
    --env GITHUB_ORG=<organization or username> \
    --env GITHUB_REPO=<name of the repository> \
    github-actions-runner:latest
```

**Warning:** If you use ``--rm`` argument, Docker will remove the container and you will loose your changes. However, this can be useful for testing purposes.

#### Adding a CI workflow on Github

If you now go to GitHub [https://github.com/&lt;your organization&gt;/&lt;your repository&gt;/settings/actions](https://github.com/%3Cyour%20organization%3E/%3Cyour%20repository%3E/settings/actions),
you should see a self-hosted runner with status "Idle":

![Self hosted runner status is Idle](/images/github-self-hosted-runners-status-idle.png)

Add the following simple workflow as ``.github/workflows/self_hosted_ci.yml`` in your repository:

```yaml
name: Self-hosted CI example
on: [push, pull_request]
jobs:
  test:
    name: test
    runs-on: self-hosted
    steps:
      - name: Show directory listing
        shell: bash -l {0}
        run: |
          ls -la
```

Now try making a change to one of the files in your repository to see if you can trigger running the simple workflow
on your self-hosted server. If successful, the status will change to "Active" while the workflow is running. You can
get an overview of previous GitHub actions by navigating to [https://github.com/&lt;your organization&gt;/&lt;your repository&gt;/actions](https://github.com/%3Cyour%20organization%3E/%3Cyour%20repository%3E/actions).


### Cleanup

To stop the running Docker container:

```shell
docker container stop github-actions-runner
```

To start the running Docker container:

```shell
docker container start github-actions-runner
```

To remove a Docker container

```shell
docker container rm github-actions-runner
```

To remove a Docker image

```shell
docker image rm github-actions-runner
```

### Extras

#### Get Docker container details

Use `docker inspect` to display details of the container
```shell
> docker inspect  github-actions-runner

[
    {
        "Id": "4ff7f37a894e351dc1def53c76959ac30c09083347ada25fcb55dc6687f8d295",
        "Created": "2020-07-02T18:29:53.106234584Z",
        "Path": "/bin/sh",
...
```

You can use the command below to only find out the IP address of the container
```shell
> docker inspect --format '{{ .NetworkSettings.IPAddress }}' github-actions-runner

172.17.0.2
```

#### Accessing Docker container

If you need access to a shell on the running Docker container:

```shell
docker exec -ti github-actions-runner /bin/bash
```

### What's next

Find instructions for provisioning additional functionality [here](../README.md).
# Setup a server for a GitHub Action runner with Vagrant

Vagrant is a tool to build a VirtualBox virtual machine (VM).

We will use a [Vagrant](https://www.vagrantup.com) to create a VM and an Ansible playbook to install a [GitHub Action runner](https://help.github.com/en/actions/hosting-your-own-runners) on it. When done a GitHub action workflow configured with `runs-on: self-hosted` will run on that runner in the VM.

## Prerequisites

* [VirtualBox](https://www.virtualbox.org/wiki/Downloads)
* [Vagrant](https://www.vagrantup.com/downloads)
* [Ansible](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html),
    I suggest using a Python virtual environment and `pip install ansible`.

## Start VM

Virtual machine can be started from this (/ubuntu-vagrant) directory with

```shell
vagrant up
```

This will have started a Ubuntu 18.04 virtual machine in VirtualBox.

## Test ssh connection to VM

Login with vagrant ssh and get the hostname with

```shell
vagrant ssh -c hostname
```

This should output `vagrant`, which is the hostname of the VM.

(If you get `Host key verification failed` error then clear previous key with `ssh-keygen -R "[127.0.0.1]:2222"` and try again)

To use other ssh client get the ssh config with

```shell
vagrant ssh-config
```

This will output something like

```shell
Host default
  HostName 127.0.0.1
  User vagrant
  Port 2222
  UserKnownHostsFile /dev/null
  StrictHostKeyChecking no
  PasswordAuthentication no
  IdentityFile /home/verhoes/git/NLESC-JCER/linux_actions_runner/ubuntu-vagrant/.vagrant/machines/default/virtualbox/private_key
  IdentitiesOnly yes
  LogLevel FATAL
```

So to login with ssh client and to get hostname use

```shell
ssh -i .vagrant/machines/default/virtualbox/private_key -p 2222 vagrant@127.0.0.1 hostname
```

It should output `vagrant`, which is the hostname of the VM.

## Configure

To use Ansible you need an [inventory file](https://docs.ansible.com/ansible/latest/user_guide/intro_inventory.html). An example inventory file called `hosts.example` should be copied to `hosts` and updated to reflect your situation.

```shell
cp hosts.example hosts
```

Ansible must be configured for which GitHub account/organization and repository it should setup a runner for.
The repository must be configured in `github_account` and `github_repo` fields in the `hosts` file.
As a repository, you can use a duplicate of [https://github.com/ci-for-science/example-python-1](https://github.com/ci-for-science/example-python-1) repository which has a workflow that runs on a self-hosted runner or any repository which has a GitHub Action workflow that has `runs-on: self-hosted`.

The Ansible playbook uses Personal Access Token for GitHub account to register the runner.
The token needs to have full admin rights for the repo. At the moment the checkbox that needs to be checked is called `repo          Full control of private repositories`. The token can be created [here](https://github.com/settings/tokens).

The token should be set as the `PAT` environment variable.

```shell
export PAT=xxxxxxxxxxxxxxx
```

## Install GitHub Action runner

To install GitHub Action runner we use an Ansible playbook to provision the VM.

Test that Ansible can ping server with

```shell
ansible vagrants -m ping
```

Should output something like

```shell
vagrant | SUCCESS => {
    "changed": false,
    "ping": "pong"
}
```

(If ping fails please check the connection configuration in `hosts` file matches output of `vagrant ssh-config`)

The playbook uses roles from [Ansible galaxy](https://galaxy.ansible.com/), they must be downloaded with

```shell
ansible-galaxy install -r requirements.yml
```

To provision VM use

```shell
ansible-playbook playbook.yml
```

The log of the runner can be viewed with

```shell
vagrant ssh -- journalctl -u actions.runner.*
```

## Destroy VM

First unregister runner with

```shell
ansible-playbook playbook.yml --tags uninstall
```

To get rid of VM use

```shell
vagrant destroy
```
# Setup GitHub Action runner on a VM on SURF HPC Cloud from Windows

Most of the steps in [/ubuntu-surf-hpc-cloud/README.md](/ubuntu-surf-hpc-cloud/README.md) can be reused except for installing Ansible.

Ansible can not be initiated from Windows powershell or command prompt.
You will need to install Ansible in [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/) using any of the [WSL OS choices](https://docs.microsoft.com/en-us/windows/wsl/install-win10#install-your-linux-distribution-of-choice) that [Ansible supports](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html) like Ubuntu 20.04.


## Example run-through

[I](https://github.com/sverhoeven) used the [https://github.com/ci-for-science/example-gpu-houston](https://github.com/ci-for-science/example-gpu-houston) repo and am already started VM in the SURF HPC cloud to run a self hosted GitHub action runner. Through out this run-through I will use my account `sverhoeven` and `example-gpu-houston` as repo name, please replace with your account/repo for your own run-through. Below are screenshots of the run-through.

I was using Ubuntu 20.04 on WSL1 with Python 3.8.2 and Ansible v2.9.10.

![Versions](ci-hpc-versions.png)

> At the time the `sleep` command did not work, so I used [this](https://github.com/microsoft/WSL/issues/4898#issuecomment-642703700) workaround

I [duplicated](https://help.github.com/en/github/creating-cloning-and-archiving-repositories/duplicating-a-repository) the [https://github.com/ci-for-science/example-gpu-houston](https://github.com/ci-for-science/example-gpu-houston) repo to my own account and made it private.

I cloned [https://github.com/ci-for-science/self-hosted-runners/](https://github.com/ci-for-science/self-hosted-runners/) repo locally.

I changed dir to `ubuntu-surf-hpc-cloud/runner/`, because the playbook and other files are there.

I copied `hosts.example` to `hosts` and updated it to connect to the VM.

I set the `PAT` environment variable to my GitHub personal access token.

Provision the runner by running the playbook with

```shell
ansible-playbook --ask-become-pass playbook.yml
```

Fill in the account and repo name.

![Fill in the account and repo name](ci-hpc-prompt.png)

Playbook ran successfully.

![Playbook ran OK](ci-hpc-playbook-end.png)

Now I made a change (commit+push) to the repo ([https://github.com/sverhoeven/example-gpu-houston](https://github.com/sverhoeven/example-gpu-houston)).

Check in [https://github.com/sverhoeven/example-gpu-houston/settings/actions](https://github.com/sverhoeven/example-gpu-houston/settings/actions) (replace with your account/repo) for runner being active and it is.

![Runner status](ci-runner-active.png)

In the actions tab we can see the job ran successfully

![Job ran OK](ci-action.png)
# Setting up a CI server for a GitHub Action runner Singularity from Linux Ubuntu

After following this guide, you'll have a simple GitHub action workflow on a GitHub repository of your choice. When new commits are made to your repository, the workflow delegates work to a server which runs in a [Singlularity](https://sylabs.io/singularity/) container. You can use this Singularity container to on a HPC cluster as a regular user and you will not need root permissions.

This guide distinguishes between the _client_ and the _server_; the client is your own machine; the server is whichever
machine will run the tests. This document describes the case where the server is a Singularity container running on your own machine.

For guides on how to configure other features in addition to just the runner, go [here](/README.md).

## Prerequisites

1. Install Singularity: https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps

2. Follow a [short tutorial](https://sylabs.io/guides/3.5/user-guide/quick_start.html#overview-of-the-singularity-interface) (Optional)

### Testing your Singularity setup

Check Singularity version
```shell
> singularity version

3.5.3
```

Pull an example image

```shell
singularity pull library://sylabsed/examples/lolcow
```

Start a shell in the Singlarity container

```shell
singularity shell lolcow_latest.sif
```

Run some test commands

```shell
Singularity> id

uid=1000(fdiblen) gid=985(users) groups=985(users),98(power),108(vboxusers),972(docker),988(storage),998(wheel)

Singularity> hostname

archlinux
```

## Server side configuration

### Build image

Now we are ready to build our Singularity image. The following command will use [Definition file](github-actions-runner-singularity.def) to build the image. It will create a system install necessary system packages and dependencies for the runner. In order to create a Singularity image, you will need root permission (or sudo) on your system.

```shell
sudo singularity build github-actions-runner-singularity.sif github-actions-runner-singularity.def
```

This command will generate ``github-actions-runner-singularity.sif`` (SIF stands for Singularity Image Format) image which we will use to set up the runner.

## Client side configuration

### Generate an OAuth token

We're almost ready to use our Docker image to set up a GitHub Runner, but first we need to
generate an OAuth token, as follows:

1. Go to [https://github.com/settings/tokens](https://github.com/settings/tokens) and click the ``Generate new token`` button.
2. Provide your GitHub password when prompted
3. Fill in a description for the token, for example _GitHub runner for github.com/&lt;your organization&gt;/&lt;your repository&gt;_
4. Enable the ``repo`` scope and all of its checkboxes, like so:

    ![Token permissions](/images/token_permissions.png)

5. Click ``Generate`` at the bottom. Make sure to copy its value because we'll need it in the next step

### Run the server

#### Preperation
Before using the Singularity image we need to set some environment variables. The Singularity container will use these environment variables to set up the runner.

```shell
export SINGULARITYENV_PERSONAL_ACCESS_TOKEN="<Github OAuth token>"
export SINGULARITYENV_RUNNER_NAME="<runner name to appear on Github>"
export SINGULARITYENV_RUNNER_WORKDIR="/tmp/actions-runner-repo"
export SINGULARITYENV_GITHUB_ORG="<organization or username>"
export SINGULARITYENV_GITHUB_REPO="<name of the repository>"
```

Create an envionment file which will be user by the runner to save some variables.

```shell
cp env.template env
```

#### Instance mode

Alternatively, you can start it as an instance (service).

```shell
singularity instance start github-actions-runner-singularity.sif github-actions-runner --writable-tmpfs --bind ./env:/opt/actions-runner/.env
```

For more information about Singularity services see [this link](https://sylabs.io/guides/3.5/user-guide/running_services.html).


To list the running instances:

```shell
singularity instance list
```

To stop the running Singularity instance:

```shell
singularity instance stop github-actions-runner
```

To start the Singularity instance again:

```shell
singularity instance start github-actions-runner
```

#### Temporary mode

Now we can run Singularity container with the following command.

```shell
singularity run \
    --writable-tmpfs \
    --bind ./env:/opt/actions-runner/.env \
    github-actions-runner-singularity.sif
```

Singularity containers by-default starts in ``read-only`` mode so you cannot make changes. While setting up the runner, some scripts needs to create a few files so we need a write access. This is achieved by adding ``--writable-tmpfs`` argument.

If you stop the running container or interrupt it by pressing to ``CTRL+C``, the Github actions runner will stop and it will be unregistered from your Github repository.

#### Accessing the logs

The singularity instances save the logs in
`~/.singularity/instances/logs/INSTANCE_NAME` folder.

## Using on a HPC Cluster

In most of the cases, the HPC user does not have root persmissions s it wont be possible to build the Singularity image. The Singularity image can be built locally and copied to
the cluster. We will use `scp` command to copy the singularity image.

```shell
scp github-actions-runner-singularity.sif USERNAME@CLUSTER_IP_ADDRESS:$REMOTE_FOLDER
```

The `$REMOTE_FOLDER` is typically your home folder which you can figure out using the command below.

```shell
echo $HOME
```

In oder to use the Singularity image on a HPC cluster, we first need to create a jobscript. The job script will be handled by the scheduler and eventually it will set up the runner when it is executed.


Example:
```
#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 480

cd $HOME/work

export SINGULARITYENV_PERSONAL_ACCESS_TOKEN="<Github OAuth token>"
export SINGULARITYENV_RUNNER_NAME="<runner name to appear on Github>"
export SINGULARITYENV_RUNNER_WORKDIR="/tmp/actions-runner-repo"
export SINGULARITYENV_GITHUB_ORG="<organization or username>"
export SINGULARITYENV_GITHUB_REPO="<name of the repository>"

srun singularity run \
    --writable-tmpfs \
    github-actions-runner-singularity.sif
```

To submit the job script using ``sbatch``:

```shell
sbatch jobscript
```

## GPU support

See [Singularity documentation](https://sylabs.io/guides/3.5/user-guide/gpu.html)

## Limitations

GitHub owned servers have some pre-installed software (see [here](https://docs.github.com/en/actions/reference/software-installed-on-github-hosted-runners)). Self hosted runners are able to download and use these software. However, this is only possible if self-hosted runner is being run on one of the supported platforms such as Ubuntu. Otherwise, users have to install the required software by themselves.
Also actions based on Docker images will not work as running Docker containers inside a Singularity container does not work.

### Extras

#### Get Singularity image details

Use `singularity inspect` to display details of the image

```shell
> singularity inspect github-actions-runner-singularity.sif*
WARNING: No SIF metadata partition, searching in container...
org.label-schema.build-date: Monday_6_July_2020_16:55:58_CEST
org.label-schema.schema-version: 1.0
org.label-schema.usage.singularity.deffile.bootstrap: library
org.label-schema.usage.singularity.deffile.from: ubuntu:19.10
org.label-schema.usage.singularity.deffile.mirrorurl: http://us.archive.ubuntu.com/ubuntu/
org.label-schema.usage.singularity.deffile.osversion: eoan
org.label-schema.usage.singularity.deffile.stage: build
org.label-schema.usage.singularity.version: 3.5.3
```

#### Accessing Singularity container

If you need access to a shell on the running Singularity container:

```shell
singularity shell \
    --writable-tmpfs \
    github-actions-runner-singularity.sif
```

### What's next

Find instructions for provisioning additional functionality [here](../README.md).
# Setting up a CI server for a GitHub Action runner with Vagrant from Windows

This guide distinguishes between the _client_ and the _server_; the client is your own machine; the server is whichever
machine will run the tests. This document describes the case where the server is VirtualBox Vm, setup through Vagrant, running on localhost.

For guides on how to configure other features in addition to just the runner, go [here](/README.md).

Vagrant is a tool to build a VirtualBox virtual machine (VM).
We will use a [Vagrant](https://www.vagrantup.com) to create a VM and an Ansible playbook to install a [GitHub Action runner](https://help.github.com/en/actions/hosting-your-own-runners) on it. When done a GitHub action workflow configured with `runs-on: self-hosted` will run on that runner in the VM.

## Prerequisites

* [VirtualBox](https://www.virtualbox.org/wiki/Downloads)
* [Vagrant](https://www.vagrantup.com/downloads)
* [Windows subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
* [Ansible installed within WSL](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html),
    I suggest using a Python virtual environment within WSL and `pip install ansible`.

## Starting a VM

Virtual machine can be started from this (/window-vagrant) directory with

```shell
vagrant up
```

This will have started a Ubuntu 18.04 virtual machine in VirtualBox.

## Login to VM

To login to VM with ssh use

```shell
vagrant ssh
```

To use other ssh client get the ssh config with

```shell
vagrant ssh-config
```

This will output something like

```shell
Host default
  HostName 127.0.0.1
  User vagrant
  Port 2222
  UserKnownHostsFile /dev/null
  StrictHostKeyChecking no
  PasswordAuthentication no
  IdentityFile /<path-to-repo/windows-vagrant/.vagrant/machines/default/virtualbox/private_key
  IdentitiesOnly yes
  LogLevel FATAL
```

**We will use Windows Subsystem for Linux (WSL) from from here on.**
Clone the repo again, this time within the disk space managed by WSL. Because of specific file permissions required, it is easier to work from the disk space with posix permissions, than to use the existing clone on the mounted Windows drive.

To have access to the machine from WSL, the private key needs to be copied from Windows to the WSL repo. To do this, copy the `/window-vagrant/.vagrant` into WSL. Correct the permissions with chmod

```shell
chmod go-rwx .vagrant/machines/default/virtualbox/private_key
```

When this is done, login with ssh using

```shell
ssh -i .vagrant/machines/default/virtualbox/private_key -p 2222 vagrant@127.0.0.1
```

## Client side configuration


## Ansible

Ansible is a tool with which you can do so-called _provisioning_, i.e. automated system administration of remote
machines. We'll use it to set up the GitHub Actions runner.

Ansible does not support Windows. Installing Ansible, and the major part of the rest of this guide, should be done from a Windows Subsystem for Linux (WSL) shell.

Install Ansible.
We recommend installing ansible through pip/PyPi. For more installation options see the list below.

- PyPI: ``pip install ansible``
- default repositories for the OS if they're available (apt, ndf, apk, homebrew)
- PPA (for Ubuntu)
- other options see [docs](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html#selecting-an-ansible-version-to-install)

Make sure your Ansible version is 2.9.9 or later with:
```shell
ansible --version
```

(Find more information [here](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html#installing-ansible-on-ubuntu)).

An example inventory file called `hosts.example` should be copied to `hosts` and updated to reflect your situation.

```shell
cp hosts.example hosts
```

Ansible must be configured for which GitHub account/organization and repository it should setup a runner.
Edit the `hosts` file and set `github_account` key and repository `github_repo` key.

The Ansible playbook uses personal Access Token for GitHub account to register the runner.
The token has to have admin rights for the repo.
Token can be created [here](https://github.com/settings/tokens).

The token should be set as the `PAT` environment variable.

```shell
export PAT=xxxxxxxxxxxxxxx
```

## Install GitHub Action runner using the playbook

To install GitHub Action runner we use an Ansible playbook to provision the VM.

### Test connection with server using ``ansible``

We're about ready to test if we can connect to the server using Ansible. For this we will use the ``ping`` module, and
we'll instruct Ansible to run the module on all hosts as defined in the inventory, as follows:

```shell
ansible all -m ping
```

Which should return:


```shell
vagrant | SUCCESS => {
    "changed": false,
    "ping": "pong"
}
```

For more complicated tasks than ``ping``, it's often inconvenient having to put everything on the command line. Instead,
a better option is to create a so-called _playbook_ containing all the steps that you want to include in your
provisioning. The playbook is a YAML file that defines a series of ``tasks``. When creating new tasks, one can start
from scratch, or make use of tasks that have been published by others (see [https://galaxy.ansible.com/](https://galaxy.ansible.com/)).



The playbook uses roles from [Ansible galaxy](https://galaxy.ansible.com/), they must be downloaded with

```shell
ansible-galaxy install -r requirements.yml
```

To provision VM use

```shell
ansible-playbook playbook.yml
```

The log of the runner can be viewed with

```shell
vagrant ssh -- journalctl -u actions.runner.*
```

## Destroy VM

First unregister runner with

```shell
ansible-playbook playbook.yml --tags uninstall
```

To get rid of VM use

```shell
vagrant destroy
```

# Using SURF HPC Cloud

This guide assumes that the user already has a VM running on [SURF HPC Cloud](https://www.surf.nl/en/hpc-cloud-your-flexible-compute-infrastructure).

We will use an Ansible playbook to install a [GitHub Action runner](https://help.github.com/en/actions/hosting-your-own-runners) on it. When done a GitHub action workflow configured with `runs-on: self-hosted` will run on that runner in the VM.

## Prerequisites

* VM running on SURF HPC Cloud
* [Ansible](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html) version 2.9.10 or later
  * Option 1 - Using a Python virtual environment and
    `pip install ansible`.
  * Option 2 - Using package manager (tested on Ubuntu 20.04)
    ```shell
        sudo apt update
        sudo apt install ansible
    ```

## Create a VM

To deploy the runner on SURF HPC Cloud, we need to create a Virtual machine (VM).
You can find the instructions [here](https://doc.hpccloud.surfsara.nl/).

## Set SSH keys for the VM

In order to access to VM, you will need to create ssh keys. Please see [this link](https://doc.hpccloud.surfsara.nl/SSHkey)

## Login to VM (optional)

If you want to test the connection with the VM, you can use the command below

```shell
ssh -i $YOUR_KEY_PATH -p 22 $USER@$HOSTNAME
```

In this command,

- `$HOSTNAME` is the IP address of the VM you created,
- `$YOUR_KEY_PATH` is the ssh key you generated in previous step.
- `$USER` is the posix user on the vm associated with the ssh key.

## Configuration

### Step-1 Creating the Ansible inventory file

To use Ansible you need an [inventory file](https://docs.ansible.com/ansible/latest/user_guide/intro_inventory.html). An example inventory file called `hosts.example` should be copied to `hosts` and updated to reflect your situation.

```shell
cp hosts.example hosts
```

Ansible playbook will ask for which GitHub account/organization and repository it should setup a runner for.

When Ansible command is executed, the Ansible playbook will ask for

- the user or rganization name
- the repository name

As a repository, you can use a clone of [https://github.com/ci-for-science/example-python-1](https://github.com/ci-for-science/example-python-1) or any repository which has a GitHub Action workflow that has [`runs-on: self-hosted`](https://github.com/ci-for-science/example-python-1/blob/4dea9c4f32a9bfcfcf166eb631c7aed3b2097d6c/.github/workflows/ci.yml#L15).

### Step-2 Generating a Github Personal Access Token

The Ansible playbook uses Personal Access Token for GitHub account to register the runner.
The token needs to have full admin rights for the repo. The only scope needed is `repo          Full control of private repositories`.

![Token permissions](/images/token_permissions.png)

The token can be created [here](https://github.com/settings/tokens).

The generated token should be set as the `PAT` environment variable.

```shell
export PAT=xxxxxxxxxxxxxxx
```

## Install GitHub Action runner

### Step 1- Testing the connection with the server
To install GitHub Action runner we use an Ansible playbook to provision the VM.

To test the connection with the server, Ansible can run ping command.

```shell
ansible all -m ping
```

If it successfully connects to the server, the output should be something like

```shell
hpc | SUCCESS => {
    "changed": false,
    "ping": "pong"
}
```

### Step 2- Installing required Ansible dependencies

The playbook uses roles from [Ansible galaxy](https://galaxy.ansible.com/), they must be downloaded with

```shell
ansible-galaxy install -r requirements.yml
```

### Step 3- Provisioning (installation on the server)

To provision VM use

```shell
ansible-playbook --ask-become-pass playbook.yml
```

To view the log of the runner, you can connect to the server via ssh and run

```shell
journalctl -u actions.runner.*
```

## Uninstalling the runner

First unregister runner with

```shell
ansible-playbook --ask-become-pass playbook.yml --tags uninstall
```

## Examples

- [Simple Python example](https://github.com/ci-for-science/example-python-1)
- [Simple GPU example](https://github.com/ci-for-science/example-gpu-houston)
# Title

Title format:

- Setting up a CI server for a GitHub Action runner with [Docker|Virtualbox|Vagrant] from [Linux Ubuntu|MacOS|Windows]
- Setting up a CI server for a GitHub Action runner on [HPC Cloud|other hardware] from [Linux Ubuntu|MacOS|Windows]

After following this guide, you'll have a simple GitHub action workflow on a GitHub repository of your choice. When new
commits are made to your repository, the workflow delegates work to a server which runs **<in a Virtual Machine on your own
computer.>**

This guide distinguishes between the _client_ and the _server_; the client is your own machine; the server is whichever
machine will run the tests. This document describes the case where the server is **<something something, e.g. a HPC
cloud machine, a VirtualBox Vm running on localhost, etc.>.**

For guides on how to configure other features in addition to just the runner, go [here](/README.md).

## Prerequisites

_Describe things that users need to install on their system to follow the guide. Things like VirtualBox, Vagrant,
Docker, Windows Subsystem for Linux, where to get iso images, how to get an account for remote hardware, etc. Out of
scope for this section: ssh, putty, Ansible_

## Server side configuration

E.g. how to configure VirtualBox, how to run docker container, how to configure HPC cloud machine

## Client side configuration

### Install Ansible

Ansible is a tool with which you can do so-called _provisioning_, i.e. automated system administration of remote
machines. We'll use it to set up the GitHub Actions runner.

Install Ansible:

- default repositories for the OS if they're available (apt, ndf, apk, homebrew)
- PPA (for Ubuntu)
- PyPI: ``pip install ansible``
- other options see [docs](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html#selecting-an-ansible-version-to-install)

Make sure your Ansible version is 2.9.9 or later with:
```shell
ansible --version
```

(Find more information [here](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html#installing-ansible-on-ubuntu)).

### Install SSH Client

e.g.

- sudo apt install openssh-client
- install putty
- homebrew install ssh

### Generate SSH key pair

Generate a key pair (files ``id_rsa`` and ``id_rsa.pub``) in directory
[``something/something``](something/something) using RSA encryption:

**Note: ``id_rsa`` is the private half of the SSH key pair; don't share it with anybody else.**

**e.g.**

```shell
cd something/something/
ssh-keygen -t rsa -f ./id_rsa -N ''
```
Make sure that the permissions are set correctly:

```shell
chmod 600 id_rsa
chmod 644 id_rsa.pub
```

Note you can use ``stat``'s ``%a`` option to see a file's permissions as an octal number, e.g.

```shell
stat -c "%a %n" <filename>
stat -c "%a %n" `ls -1`
```


### Copy the key pair to server

Copy the public half of the key pair (i.e. ``id_rsa.pub``) to the server.

**e.g.**

```shell
ssh-copy-id -i ./id_rsa.pub -p 2222 tester@127.0.0.1
```


### Test connection with server using ``ssh``

Test if you can SSH into the server using the other half of the key pair (i.e. ``id_rsa``)

**e.g.**

```shell
ssh -i ./id_rsa -p 2222 tester@127.0.0.1
```

If you get a ``Host key verification failed`` error, clear the existing key with

```shell
ssh-keygen -R "[127.0.0.1]:2222"
```

and try again.


Log out of the server with

```shell
exit
```

### Troubleshooting SSH

Getting SSH connections to work can be tricky. Check out [this document](/docs/troubleshooting-ssh.md) if you're
experiencing difficulties.

### The inventory file

Ansible uses so-called _inventory_ files to define how to connect to remote machines. The inventory file is typically
called  ``hosts``. The following inventory file is equivalent to the ``ssh`` command line we just used:

```yaml
all:
  hosts:
    ci-server:
      ansible_connection: ssh
      ansible_host: 127.0.0.1
      ansible_port: 2222
      ansible_ssh_private_key_file: ./id_rsa
      ansible_user: tester
```

This inventory file defines a group ``all`` with just one machine in it, which we labeled ``ci-server``. ``ci-server``
has a bunch of variables that define how to connect to it. For more information on inventory files, read
[Ansible's documentation](https://docs.ansible.com/ansible/latest/user_guide/intro_inventory.html).

### The Ansible configuration file

In addition to the inventory file, it's often convenient to use a configuration file. The default filename for this file
is ``ansible.cfg``, and it can be used to specify Ansible's behavior. The configuration option documentation can be
found [here](https://docs.ansible.com/ansible/latest/reference_appendices/config.html#the-configuration-file).

### Test connection with server using ``ansible``

We're about ready to test if we can connect to the server using Ansible. For this we will use the ``ping`` module, and
we'll instruct Ansible to run the module on all hosts as defined in the inventory, as follows:

```shell
ansible all -m ping
```

Which should return:

```text
ci-server | SUCCESS => {
    "changed": false,
    "ping": "pong"
}
```

### Install the runner using the playbook

For more complicated tasks than ``ping``, it's often inconvenient having to put everything on the command line. Instead,
a better option is to create a so-called _playbook_ containing all the steps that you want to include in your
provisioning. The playbook is a YAML file that defines a series of ``tasks``. When creating new tasks, one can start
from scratch, or make use of tasks that have been published by others (see [https://galaxy.ansible.com/](https://galaxy.ansible.com/)).

We're almost ready to use ``ansible-playbook`` to set up a GitHub Runner on your own server, but first we need to
generate an OAuth token, as follows:

1. Go to [https://github.com/settings/tokens](https://github.com/settings/tokens) and click the ``Generate new token`` button.
1. Provide your GitHub password when prompted
1. Fill in a description for the token, for example _GitHub runner for github.com/&lt;your organization&gt;/&lt;your repository&gt;_
1. Enable the ``repo`` scope and all of its checkboxes, like so:

    ![Token permissions](/images/token_permissions.png)

1. Click ``Generate`` at the bottom. Make sure to copy its value because we'll need it in the next step

Configuring your server such that it can run continuous integration requires 4 pieces of information, for which you will be prompted:

1. Because our playbook requires elevated permissions, the command uses the ``--ask-become-pass`` option to prompt for
the root password. Fill in the password ``password`` to become ``root`` in the server.
1. Fill in the GitHub organization (which might be simply your GitHub user name) and ...
1. ...the repository name for which you want to run workflows on a self-hosted server
1. Finally, you need to supply the Personal Access Token

Now run this command to provision the GitHub Action runner on your server:

```shell
ansible-playbook playbook.yml --ask-become-pass
```

If you now go to GitHub [https://github.com/&lt;your organization&gt;/&lt;your repository&gt;/settings/actions](https://github.com/%3Cyour%20organization%3E/%3Cyour%20repository%3E/settings/actions),
you should see a self-hosted runner with status "Idle":

![Self hosted runner status is Idle](/images/github-self-hosted-runners-status-idle.png)

Add the following simple workflow as ``.github/workflows/self_hosted_ci.yml`` in your repository:

```yaml
name: Self-hosted CI example

on: [push, pull_request]

jobs:
  test:
    name: test
    runs-on: self-hosted
    steps:
      - name: Show directory listing
        shell: bash -l {0}
        run: |
          ls -la
```

Now try making a change to one of the files in your repository to see if you can trigger running the simple workflow
on your self-hosted server. If successful, the status will change to "Active" while the workflow is running. You can
get an overview of previous GitHub actions by navigating to [https://github.com/&lt;your organization&gt;/&lt;your repository&gt;/actions](https://github.com/%3Cyour%20organization%3E/%3Cyour%20repository%3E/actions).


### Monitoring the runner service's logs

The log of the runner can be viewed with

```shell
ssh -i <keyfile> -p <port> <username>@<hostname>
```

Then

```shell
journalctl -u actions.runner.*
```

### Start the runner each time the machine boots

```shell
ansible-playbook playbook.yml --tags enable
```

### Start the runner

```shell
ansible-playbook playbook.yml --tags start
```

### Managing the runner service through the playbook

```shell
ansible-playbook playbook.yml --tags start
ansible-playbook playbook.yml --tags stop
ansible-playbook playbook.yml --tags restart
ansible-playbook playbook.yml --tags status
ansible-playbook playbook.yml --tags enable
ansible-playbook playbook.yml --tags disable
```

Uninstalling the runner

```shell
ansible-playbook playbook.yml --tags uninstall
```

### Verify that your newly configured runner is triggered

Add the following simple workflow as ``.github/workflows/self_hosted_ci.yml`` in your repository
[https://github.com/&lt;your organization&gt;/&lt;your repository&gt;](https://github.com/%3Cyour%20organization%3E/%3Cyour%20repository%3E):

```yaml
name: Self-hosted CI example

on: [push, pull_request]

jobs:
  test:
    name: test
    runs-on: self-hosted
    steps:
      - name: Show directory listing
        shell: bash -l {0}
        run: |
          ls -la
```

With this workflow in place, new pushes and new pull requests should trigger your self-hosted server.
Try making a change to one of the files in your repository to see if you can trigger running the simple workflow
on your self-hosted server. If successful, the status will change to "Active" while the workflow is running.
You can see a record of past and current GitHub Actions by pointing your browser to
[https://github.com/&lt;your organization&gt;/&lt;your repository&gt;/actions?query=workflow:"Self-hosted+CI+example"](https://github.com/%3Cyour%20organization%3E/%3Cyour%20repository%3E/actions?query=workflow%3A%22Self-hosted+CI+example%22).


### What's next

Find instructions for provisioning additional functionality [here](../README.md).
# Troubleshooting ``SSH``

Note: these troubleshooting instructions assume Linux as the operating system on both the client side and the server
side.

## File permissions

The first thing to check is whether your system has the correct permissions on the following files (you can check the
octal representation of the file permission with: ``stat -c %a <filename>``):

```shell
# client-side
chmod go-w $HOME
chmod 700 $HOME/.ssh
chmod 600 $HOME/.ssh/config       # (documentation varies 644, 600, 400)
chmod 600 $HOME/.ssh/id_rsa       # (private keys, rsa and other types)
chmod 644 $HOME/.ssh/id_rsa.pub   # (public keys, rsa and other types)
chmod 600 $HOME/.ssh/known_hosts  # (not documented)

# server side:
chmod 600 <other system''s home dir>/.ssh/authorized_keys
```

## Ownership

All files and directories under `~/.ssh`, as well as `~/.ssh` itself, should be owned by the user.

```shell
chown -R $(id -u):$(id -g) ~/.ssh
```

## Verbosity

Increase ``ssh``'s verbosity using the ``-vvvv`` option (more ``v``'s means higher verbosity), e.g.

```shell
ssh -vvv username@hostname
```

Another useful option is to ask ``ssh`` for a list of its configuration options and their values with the ``-G`` option,
e.g.

```shell
ssh -G anyhost
ssh -G user@some.system.com
```

Sometimes, a connection cannot be set up because of a configuration problem on the server side. If you have access to
the server through another way, running

```shell
sshd -T
```

might help track the problem down. Note that the results may be user-dependent, for example the result may be different
for ``root`` or for a user.

## Configuration settings

1. client-side, global user configuration: ``/etc/ssh/ssh_config``
1. client-side, local user configuration ``$HOME/.ssh/config``
1. server-side, global system configuration for ssh server daemon ``/etc/ssh/sshd_config``

## Related to `known_hosts` file


1. host name hashed or not ``hashKnownHosts``
1. strict host key checking ``strictHostKeyChecking``
1. removing a given host's key goes like this

    ```shell
    ssh-keygen -R [localhost]:10022
    ```

## Encrypted ``/home``

Using encryption on your home drive may negatively affect things: <https://help.ubuntu.com/community/SSH/OpenSSH/Keys>
# Intro

This short guide shows how to install Nvidia drivers and CUDA for GRID K2 hardware. For more information about GPUs on SURF HPC Cloud please visit [SURF HPC Documentation](https://doc.hpccloud.surfsara.nl/gpu-attach).

We have 2 methods to install Nvidia drivers and CUDA.
- Using ansible playbook
- Following the steps and installing manually

## Method 1: Installation using ansible playbook

The command below runs the ansible-playbook and installs all necassary software. This ansible-playbook is specifically written for SURF HPC Cloud platform and GRID K2 hardware. However, it can easily be adapted for different platforms and graphics cards.

Copy example hosts file and make necassary changes.

```shell
cp hosts.example hosts
```

```shell
docker run --rm -ti -v $PWD:/data --workdir=/data ansible/ansible-runner ansible-playbook playbook-install-cuda-gridk2.yml
```

``id_rsa`` is your rivate ssh key and ``hosts`` is the file which has the connection details of the server.

## Method 2: Manual installation

### Requirements

CUDA currently officially supports only two versions of Ubuntu: 18.04 and 16.04. This instructions were tested on Ubuntu 18.04.

### System info

Distribution info

```shell
lsb_release -a
No LSB modules are available.
Distributor ID: Ubuntu
Description:    Ubuntu 18.04.4 LTS
Release:        18.04
Codename:       bionic
```

Kernel version

```shell
uname -a
Linux packer-Ubuntu-18 4.15.0-101-generic #102-Ubuntu SMP Mon May 11 10:07:26 UTC 2020 x86_64 x86_64 x86_64 GNU/Linux
```

GPU hardware information

```shell
lspci | grep -i nvidia
01:01.0 VGA compatible controller: NVIDIA Corporation GK104GL [GRID K2] (rev a1)
```

### Pre install

For Grid K2 card we will need CUDA 8.0. CUDA 8.0 only works with only gcc 5.0 so it should be installed before.

To decide what version of CUDA and Nvidia drivers you need, please check the links below.

See what drivers you need:
[https://www.nvidia.com/Download/index.aspx?lang=en-us](https://www.nvidia.com/Download/index.aspx?lang=en-us)

Check compatibility first:
[https://docs.nvidia.com/deploy/cuda-compatibility/index.html](https://docs.nvidia.com/deploy/cuda-compatibility/index.html)

```shell
apt install gcc-5 g++-5
```

### Install Nvidia drivers

Download Nvidia driver (version 367) and install it.

```shell
wget http://us.download.nvidia.com/XFree86/Linux-x86_64/367.134/NVIDIA-Linux-x86_64-367.134.run
sh ./NVIDIA-Linux-x86_64-367.134.run --accept-license  -s
```

### Install CUDA

Download CUDA 8.0 installer
```shell
# wget https://developer.nvidia.com/compute/cuda/8.0/Prod2/local_installers/cuda_8.0.61_375.26_linux-run
```

Download Patch release:

```shell
wget https://developer.nvidia.com/compute/cuda/8.0/Prod2/patches/2/cuda_8.0.61.2_linux-run
```

While installing CUDA, we had some issues related to Perl scripts.
See: [https://forums.developer.nvidia.com/t/cant-locate-installutils-pm-in-inc/46952/10](https://forums.developer.nvidia.com/t/cant-locate-installutils-pm-in-inc/46952/10)


These commands solves the Perl issues.

```shell
sh ./cuda_8.0.61_375.26_linux-run  --tar mxvf
cp InstallUtils.pm /usr/lib/x86_64-linux-gnu/perl-base/
export $PERL5LIB
rm -rf InstallUtils.pm cuda-installer.pl run_files uninstall_cuda.pl
```

After fixing the Perl issue, we can install CUDA.

```shell
sh ./cuda_8.0.61_375.26_linux-run --silent --samples --toolkit --override --verbose
```

### Environment variables

In order to be able to use CUDA, we need to change our environment variables.

Add the lines below to .profile file.

```shell
export PATH=$PATH:/usr/local/cuda-8.0/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-8.0/lib64
```

## Test the installation

### CUDA compiler

Check the CUDA compiler version.

```shell
nvcc --version
nvcc: NVIDIA (R) Cuda compiler driver
Copyright (c) 2005-2017 NVIDIA Corporation
Built on Fri_Nov__3_21:07:56_CDT_2017
Cuda compilation tools, release 9.1, V9.1.85
```

### Example code

See [gpu-houston](https://github.com/ci-for-research/example-gpu-houston) for a simple example code.
# Setting up a CI server with a GitHub Action runner using VirtualBox, from Linux Ubuntu

After following this guide, you'll have a simple GitHub action workflow on a GitHub repository of your choice. When new
commits are made to your repository, the workflow delegates work to a server which runs in a Virtual Machine on your own
computer.

This guide distinguishes between the _client_ and the _server_; the client is your own machine; the server is whichever
machine runs the tests. This document describes the case where the server is a virtual machine, running on your own
physical machine. For guides on how to configure other features in addition to just the runner, go [here](/README.md).

## TL;DR

1. create a virtual machine with an SSH server
1. enable access to the server via SSH keys
1. ``ansible-playbook --ask-become-pass playbook.yml``

## Prerequisites

1. Install VirtualBox on the client: [https://www.virtualbox.org/wiki/Linux_Downloads](https://www.virtualbox.org/wiki/Linux_Downloads)
1. Download an Ubuntu iso image from [https://ubuntu.com/#download](https://ubuntu.com/#download). Both the desktop and the server variant are
suitable --choose whichever you're comfortable with.

## Server side configuration

1. Create a new virtual machine in VirtualBox. It's recommended to give it at least 4 GB memory, 2 CPUs, and 20 GB disk space (dynamically allocated).
1. For the new virtual machine, go to _Settings_ > _Storage_, then under _IDE controller_ select the item marked _Empty_. Then click the icon to load something into the virtual optical disk, then select the Ubuntu iso file.
1. For the new virtual machine, go to _Settings_ > _Network_
    1. On the _Adapter 1_ tab,
        - make sure that the _Enable Network Adapter_ checkbox is checked
        - set the _Attached to_ dropdown menu to _NAT_
        - Click _Advanced_, then _Port Forwarding_
        - Add a new rule, with _Protocol_ TCP, _HostIP_ 127.0.0.1, _Host Port_ 2222, leave _Guest IP_ empty, and _Guest Port_ 22
1. Start the Virtual Machine for the first time.
1. In Ubuntu's install wizard, call the user ``tester``
1. In Ubuntu's install wizard, set the user's password to ``password``
1. Update packages

    ```shell
    sudo apt update
    sudo apt upgrade
    ```

1. Configure an SSH server (OpenSSH) for remote connection; check permissions on relevant files and directories:

    ```shell
    sudo apt install openssh-server
    chmod go-w /home/tester
    mkdir /home/tester/.ssh
    chmod 700 /home/tester/.ssh
    touch /home/tester/.ssh/known_hosts && chmod 644 /home/tester/.ssh/known_hosts
    touch /home/tester/.ssh/config      && chmod 600 /home/tester/.ssh/config
    chown -R tester:tester /home/tester/.ssh
    ```

    Note you can use ``stat``'s ``%a`` option to see a file's permissions as an octal number, e.g.

    ```shell
    stat -c "%a %n" <filename>
    ```

## Client side configuration

### Install Ansible

Ansible is a tool with which you can do so-called _provisioning_, i.e. automated system administration of remote
machines. We'll use it to set up the GitHub Actions runner.

Install Ansible from Ubuntu's repositories:

```shell
sudo apt install ansible
```

Make sure your Ansible version is 2.9.9 or later with:
```shell
ansible --version
```

(Find more information [here](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html#installing-ansible-on-ubuntu)).

### Install SSH client

To be able to connect to remote machines via SSH, we'll need an SSH client. We'll use OpenSSH. Install it from Ubuntu's
repositories with:

```shell
sudo apt install openssh-client
```

### Generate SSH key pair

Generate a key pair (files ``id_rsa`` and ``id_rsa.pub``) in directory
[``/ubuntu-virtualbox/``](/ubuntu-virtualbox/) using RSA encryption:

**Note: ``id_rsa`` is the private half of the SSH key pair; don't share it with anybody else.**

```shell
cd ubuntu-virtualbox/
ssh-keygen -t rsa -f ./id_rsa -N ''
```

Make sure that the permissions are set correctly:

```shell
chmod 600 id_rsa
chmod 644 id_rsa.pub
```

Note you can use ``stat``'s ``%a`` option to see a file's permissions as an octal number, e.g.

```shell
stat -c "%a %n" <filename>
stat -c "%a %n" `ls -1`
```

### Copy the key pair to the server

Copy the public half of the key pair (i.e. ``id_rsa.pub``) to the server.

```shell
ssh-copy-id -i ./id_rsa.pub -p 2222 tester@127.0.0.1
```

### Test connection with server using ``ssh``

Test if you can SSH into the server using the other half of the key pair (i.e. ``id_rsa``)

```shell
ssh -i ./id_rsa -p 2222 tester@127.0.0.1
```

If you get a ``Host key verification failed`` error, clear the existing key with

```shell
ssh-keygen -R "[127.0.0.1]:2222"
```

and try again.

Log out of the server with

```shell
exit
```

### Troubleshooting SSH

Getting SSH connections to work can be tricky. Check out [this document](/docs/troubleshooting-ssh.md) if you're
experiencing difficulties.

### The inventory file

Ansible uses so-called _inventory_ files to define how to connect to remote machines. The inventory file is typically
called  ``hosts``. The following inventory file is equivalent to the ``ssh`` command line we just used:

```yaml
all:
  hosts:
    ci-server:
      ansible_connection: ssh
      ansible_host: 127.0.0.1
      ansible_port: 2222
      ansible_python_interpreter: /usr/bin/python3
      ansible_ssh_private_key_file: ./id_rsa
      ansible_user: tester
```

This inventory file defines a group ``all`` with just one machine in it, which we labeled ``ci-server``. ``ci-server``
has a bunch of variables that define how to connect to it. For more information on inventory files, read
[Ansible's documentation](https://docs.ansible.com/ansible/latest/user_guide/intro_inventory.html).

### The Ansible configuration file

In addition to the inventory file, it's often convenient to use a configuration file. The default filename for this file
is ``ansible.cfg``, and it can be used to specify Ansible's behavior. The configuration option documentation can be
found [here](https://docs.ansible.com/ansible/latest/reference_appendices/config.html#the-configuration-file).

### Test connection with server using ``ansible``

We're about ready to test if we can connect to the server using Ansible. For this we will use the ``ping`` module, and
we'll instruct Ansible to run the module on all hosts as defined in the inventory, as follows:

```shell
ansible all -m ping
```

Which should return:

```text
ci-server | SUCCESS => {
    "changed": false,
    "ping": "pong"
}
```

### Install the runner using the playbook

For more complicated tasks than ``ping``, it's often inconvenient having to put everything on the command line. Instead,
a better option is to create a so-called _playbook_ containing all the steps that you want to include in your
provisioning. The playbook is a YAML file that defines a series of ``tasks``. When creating new tasks, one can start
from scratch, or make use of tasks that have been published by others (see [https://galaxy.ansible.com/](https://galaxy.ansible.com/)).

We're almost ready to use ``ansible-playbook`` to set up a GitHub Runner on your own server, but first we need to
generate an OAuth token, as follows:

1. Make a copy of the template file. We will store your token in the copied file momentarily.

    ```shell
    cp secret.yml.template secret.yml
    ```
1. Go to [https://github.com/settings/tokens](https://github.com/settings/tokens) and click the ``Generate new token`` button.
1. Provide your GitHub password when prompted
1. Fill in a description for the token, for example _Token for self-hosted GitHub runners_
1. Enable the ``repo`` scope and all of its checkboxes, like so:

    ![Token permissions](/images/token_permissions.png)

1. Click ``Generate`` at the bottom, and update the value of ``PERSONAL_ACCESS_TOKEN`` in ``secret.yml``. **Don't share the contents of ``secret.yml``.**

Configuring your server such that it can run continuous integration requires 4 pieces of information, for which you will be prompted:

1. Because our playbook requires elevated permissions, the command uses the ``--ask-become-pass`` option to prompt for
the root password. Fill in the password ``password`` to become ``root`` in the server.
1. Fill in the GitHub organization (which might be simply your GitHub user name) and ...
1. ...the repository name for which you want to run workflows on a self-hosted server
1. Specify how you want the runner to show up in the GitHub interface


Now run this command to provision the GitHub Action runner on your server:

```shell
ansible-playbook --ask-become-pass playbook.yml
```

If you now go to GitHub [https://github.com/&lt;your organization&gt;/&lt;your repository&gt;/settings/actions](https://github.com/%3Cyour%20organization%3E/%3Cyour%20repository%3E/settings/actions),
you should see a self-hosted runner with status "Idle":

![Self hosted runner status is Idle](/images/github-self-hosted-runners-status-idle.png)

### Monitoring the runner service's logs

The log of the runner can be viewed with

```shell
ssh -i <keyfile> -p <port> <username>@<hostname>
```

Then

```shell
journalctl -u actions.runner.*
```

### Start the runner each time the machine boots

```shell
ansible-playbook --ask-become-pass playbook.yml --tags enable
```

### Start the runner

```shell
ansible-playbook --ask-become-pass playbook.yml --tags start
```

### Managing the runner service through the playbook

```shell
ansible-playbook --ask-become-pass playbook.yml --tags start
ansible-playbook --ask-become-pass playbook.yml --tags stop
ansible-playbook --ask-become-pass playbook.yml --tags restart
ansible-playbook --ask-become-pass playbook.yml --tags status
ansible-playbook --ask-become-pass playbook.yml --tags enable
ansible-playbook --ask-become-pass playbook.yml --tags disable
```

Uninstalling the runner

```shell
ansible-playbook --ask-become-pass playbook.yml --tags uninstall
```

### Verify that your newly configured runner is triggered

Add the following simple workflow as ``.github/workflows/self_hosted_ci.yml`` in your repository
[https://github.com/&lt;your organization&gt;/&lt;your repository&gt;](https://github.com/%3Cyour%20organization%3E/%3Cyour%20repository%3E):

```yaml
name: Self-hosted CI example

on: [push, pull_request]

jobs:
  test:
    name: test
    runs-on: self-hosted
    steps:
      - name: Show directory listing
        shell: bash -l {0}
        run: |
          ls -la
```

With this workflow in place, new pushes and new pull requests should trigger your self-hosted server.
Try making a change to one of the files in your repository to see if you can trigger running the simple workflow
on your self-hosted server. If successful, the status will change to "Active" while the workflow is running.
You can see a record of past and current GitHub Actions by pointing your browser to
[https://github.com/&lt;your organization&gt;/&lt;your repository&gt;/actions?query=workflow:"Self-hosted+CI+example"](https://github.com/%3Cyour%20organization%3E/%3Cyour%20repository%3E/actions?query=workflow%3A%22Self-hosted+CI+example%22).

### What's next

Find instructions for provisioning additional functionality [here](/README.md).
