# Badges

| fair-software.nl recommendations | Badge |
|:-|:-:|
| [1. Code Repository](https://fair-software.nl/recommendations/repository) | [![GitHub](https://img.shields.io/github/last-commit/nlesc/lokum)](https://img.shields.io/github/last-commit/nlesc/lokum) |
| [2. License](https://fair-software.nl/recommendations/license) | [![License](https://img.shields.io/github/license/nlesc/lokum)]((https://img.shields.io/github/license/nlesc/lokum)) |
| [3. Community Registry](https://fair-software.nl/recommendations/registry) | [![Research Software Directory](https://img.shields.io/badge/rsd-lokum-00a3e3.svg)](https://www.research-software.nl/software/lokum) |
| [4. Enable Citation](https://fair-software.nl/recommendations/citation) | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3482939.svg)](https://doi.org/10.5281/zenodo.3482939) |
| [5. Code Quality Checklist](https://fair-software.nl/recommendations/checklist) | [![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/3754/badge)](https://bestpractices.coreinfrastructure.org/projects/3754)  |


# Lokum

This repository includes scripts to deploy a cluster with GlusterFS, Docker, Spark and JupyterHub services on bare-metal. Currently, it only supports Opennebula platform.

Lokum uses [emma](https://github.com/nlesc-sherlock/emma) ansible playbooks to deploy services.

## Technologies & Tools

- [Terraform Client](https://www.terraform.io)
- [Runtastic Terraform Opennebula provider](https://github.com/runtastic/terraform-provider-opennebula)
- [Ansible](https://www.ansible.com/)
- [emma](https://github.com/nlesc-sherlock/emma)

## Usage

### 1-Pull the Docker image from Docker Hub

```bash
docker pull nlesc/lokum:latest
```

### 2-Settings

#### 2.1 VM configuration (template)

Edit **config/opennebula_k8s.tpl** to adjust the following VM settings:

    CPU = "2.0"
    VCPU = "2"
    IMAGE_ID = "YOUR_IMAGE_ID"
    MEMORY = "4096"
    NIC = [
      NETWORK = "INTERNAL_NETWORK_NAME",
      NETWORK_UNAME = "NETWORK_USERNAME" ]

There are two **SIZE** variables. The first one is for the cluster itselft and the second one is for the persistent storage. The default values are about 15G and 30G.

#### 2.2 Credentials

Edit **config/variables.tf** and set user credentials.

### 3-Deploy the cluster

```bash
docker run --rm --net=host -it \
  -v $(pwd)/config:/lokum/config \
  -v $(pwd)/deployment:/lokum/deployment \
  nlesc/lokum:latest
```

Confirm the planned changes by typing **yes**

Configuration and the ssh-keys of each deployed cluster will be stored under **deployment/clusterX** folder.

## Connecting to the nodes

### ssh to nodes

You can connect to the nodes using generated ssh keys. For example:

```bash
ssh -i ./deployment/cluster0/id_rsa_lokum_root.key root@SERVER_IP
or
ssh -i ./deployment/cluster0/id_rsa_lokum_ubuntu.key ubuntu@SERVER_IP
```

## Starting the services manually

Run the following command in main lokum directory.

```bash
docker run --rm --net=host -it \
  -v $(pwd)/config:/lokum/config \
  -v $(pwd)/deployment:/lokum/deployment \
  nlesc/lokum:latest
```

```bash
DEPLOYMENT_DIR=/lokum/deployment/cluster0; ANSIBLE_HOST_KEY_CHECKING=False; export CLUSTER_NAME=lokum; cd /lokum/emma/vars; sh ./create_vars_files.sh; cd /lokum/emma; ansible-playbook -i ${DEPLOYMENT_DIR}/hosts.yaml --extra-vars 'CLUSTER_NAME=lokum' start_platform.yml --skip-tags 'jupyterhub,cassandra' --private-key=${DEPLOYMENT_DIR}/id_rsa_lokum_ubuntu.key -v
```

To check Apache Spark open the link below in a browser:
http://NODE_1_IP:8080/

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
reported by contacting the project team at . All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
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
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
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
