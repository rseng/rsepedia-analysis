# Badges

| fair-software.nl recommendations | Badge |
|:-|:-:|
| [1. Code Repository](https://fair-software.nl/recommendations/repository) | [![GitHub](https://img.shields.io/github/last-commit/NLeSC-GO-common-infrastructure/marzipan)](https://img.shields.io/github/last-commit/NLeSC-GO-common-infrastructure/marzipan) |
| [2. License](https://fair-software.nl/recommendations/license) | [![License](https://img.shields.io/github/license/NLeSC-GO-common-infrastructure/marzipan)]((https://img.shields.io/github/license/NLeSC-GO-common-infrastructure/marzipan)) |
| [3. Community Registry](https://fair-software.nl/recommendations/registry) | [![Research Software Directory](https://img.shields.io/badge/rsd-marzipan-00a3e3.svg)](https://www.research-software.nl/software/marzipan) |
| [4. Enable Citation](https://fair-software.nl/recommendations/citation) | [![DOI](https://zenodo.org/badge/303491969.svg)](https://zenodo.org/badge/latestdoi/303491969) |
| [5. Code Quality Checklist](https://fair-software.nl/recommendations/checklist) | [![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/3754/badge)](https://bestpractices.coreinfrastructure.org/projects/3754)  |


# marzipan
Automated instantiation and deployment of (clusters of) virtual machine(s) on bare metal using the OpenNebula platform, as well as subsequent provisioning and deployment of services incl., e.g. Dask.

`marzipan` consists of the core `marzipan.py` python [module](https://github.com/NLeSC-GO-common-infrastructure/marzipan/tree/improve-documentation/Marzipan) providing a high level interface to the OpenNebula cloud,
as well as an accompanying [Docker framework](https://github.com/NLeSC-GO-common-infrastructure/marzipan/tree/improve-documentation/Docker) and configurable [deployment scripts](https://github.com/NLeSC-GO-common-infrastructure/marzipan/tree/improve-documentation/marzipan_scripts) providing a fully automated instantiation and provisioning environment.

For provisioning marzipan makes use of the [`emma_marzipan` fork](https://github.com/NLeSC-GO-common-infrastructure/emma/tree/emma_marzipan) ansible playbooks.

`marzipan` is based off and strongly draws from [`Lokum`](https://github.com/NLeSC/lokum), but is updated to make use of current versions of Ansible as well as python 3, and circumvents recurrent synchronicity and timeout issues arsing from the interplay of terraform, the runtastic OpenNebula provider for terraform, and various (legacy) OpenNebula versions.

`marzipan` has been tested on the SURFsara HPC cloud, but should work for any OpenNebula platform.

## Technologies and tools

- [OpenNebula](https://opennebula.io)
- [Docker](https://www.docker.com)
- [Ansible](https://www.ansible.com)
- [emma_marzipan fork](https://github.com/NLeSC-GO-common-infrastructure/emma/tree/emma_marzipan) of [emma](https://github.com/nlesc-sherlock/emma)



## Usage

### 1 Clone repository
To make use of `marzipan` the user should clone this repository to their local system. Further instructions on the use of `marzipan` assume a full replica of the repository on the users local system.

### 2 Adjust configuration and template
The user should modify the `ClusterConf.ini` file located in the [`config`](./config) subdirectory, as well as the `opennebula_goera.tpl` file in the [`templates`](./templates) subdirectory to match their requirements.

#### 2.1 configuration
The `ClusterConf.ini` file enables the user to set desired configuration values such as the number of nodes, the name of the VMs, the OpenNebula endpoint and their credentials.
The [`config`](./config) subdirectory of the repository includes a file `ClusterConf.ini.example` which can be appropriately modified and subsequently renamed.

#### 2.2 template
The user must supply a template file specifiying the desired configuration for the VM(s) to be created.
An example, `opennebula_goera.tpl`, is provided in the [`templates`](./templates) subfolder of the repository.
In particular, the following fields will require modification:
```
CONTEXT = [
    GROUP = "your_group"]
DISK = [
    DATASTORE = "nameOfYourBaseImageDataStore"
    DATASTORE_ID = "IDOfYourBaseImageDataStore"
    IMAGE_ID ="IDOfYourBaseImage"
]
```

Please bear in mind, that the base image for the OS disk must be made available for the user (with the credentials being used) before executing `marzipan`. This is up to the user and can be accomplished using the OpenNebula user interface.


### 3 Build Docker image
Change directories to the [`Docker`](./Docker) subdirectory.
Build the `nlesc/marzipan` docker image by running
```bash
./build_marzipan.sh
```
This creates the image with tag set to `latest`.

### 4 Run Docker framework to instantiate and provision a (cluster of) VM(s)
Change back to the root directory of the repository.
The docker framework can then be used to instantiate and provision a cluster of VMs by running
```bash
./deployCluster.sh
```
. The `root` and `ubuntu` user ssh keys generated for the cluster, as well as the `hosts.yaml` file enabling provisioning with the `emma` platform leveraging `ansible` are written to the `deployments` subdirectory (is created on execution) in a subfolder with the clusters name. They can be used to subsequently interact with the cluster.

The user can adapt the provisioning by modifying the `marzipan_deploy.py` script in the `marzipan_scripts` subdirectory in the section below
```
"""
emma based provisioninig 
"""
```
 The user is referred to the [`emma_marzipan` fork](https://github.com/NLeSC-GO-common-infrastructure/emma/tree/emma_marzipan) for supported options.

 __NOTE__: changes to the `marzipan_deploy.py` script require the [docker image](#3-build-docker-image) to be rebuilt before taking effect.


## Access to the cluster
The cluster VMs can be accesed via ssh as `ubuntu` or `root` user, using the generated keys. For example:
```bash
ssh -i ./deployments/<clustername>/id_rsa_marzipan_root.key root@SERVER_IP
or
ssh -i ./deployments/<clustername>/id_rsa_marzipan_ubuntu.key ubuntu@SERVER_IP
```


## The marzipan OpenNebula interface

The core `marzipan.py` module (located in the [`Marzipan`](./Marzipan) subfolder) provides the `one_interface` class. The class' methods provide a high level interface to set up a cluster of VMs on the (SURFsara) OpenNebula cloud.

`marzipan.py` can be imported, providing access to the class methods, or run as a script to fully automatedly set up a cluster of VMs. `marzipan.py` also provides a high level class method `deploy_cluster()` which corresponds to the execution as a script.

`marzipan` is complemented by a `ClusterConf.ini` file (in [`config`](./config)), where the user can set desired configuration values such as the number of nodes, the name of the VMs, the OpenNebula endpoint and their credentials. The repository includes a file `ClusterConf.ini.example` which can be appropriately modified and subsequently renamed.

Furthermore, the user should supply a template file specifying the desired VM configuration. An example, `opennebula_goera.tpl`, is provided in the [`templates`](./templates) subfolder of the repository. Please bear in mind, that the base image for the OS disk must be made available for the user (with the credentials being used) before executing `marzipan`.

Finally, the user can provide a public ssh key file for the `root` user to be included with the template. Alternatively, the ssh key can be supplied in the `ClusterConf.ini` file.
NOTE: if using the Docker framework, the use should refrain from, or take great care in, changing the settings relating to the ssh keys, as these are autogenerated during deployment.

When run as a script or by invoking the full deployment method `mazipan` will construct a VM template in OpenNebula, and deploy the requested number of VMs based on this template. `marzipan` will then monitor the deployment, only reporting successful execution when all VMs are in the `RUNNING` LCM_STATE. If this has failed to complete after 120 seconds marzipan will exit, notifying the user of failure.


## Reference/Documentation

[OpenNebula 5.2 Documentation for the XML-RPC API](http://docs.opennebula.io/5.2/integration/system_interfaces/api.html#actions-for-templates-management)

[PYONE bindings documentation](http://docs.opennebula.io/5.12/integration/system_interfaces/python.html)
# ClusterConf.ini configuration file

The `ClusterConf.ini` file allows the user to set configuration parameters for instantiating a (cluster of) VM(s) on bare metal using the OpenNebula platform.

The parameters to be supplied are divided into three setions:

## one (OpenNebula)
This section covers access to the OpenNebula platform. The user should specify the endpoint, as well as their login credentials. __Note__: Please set access to the ClusterConf.ini file at an appropriate level (e.g. `chmod 600`). DO NOT push this file to a remote repository. The `.gitignore` file currently precludes this happening.

## ssh
This section specifies the path to the root public ssh key to be added to each VM created. The current setting:
```
root_public_key = /marzipan/deployments/tmp/id_rsa_marzipan_root.key.pub
```
ensures comaptibility with the automatic deployment options. The user may modify this HOWEVER, they are then responsible for ensuring compatability with the automatic deployment opions if desired.

## cluster
This section specifies which template file should be used and how many VM(s) will be created. Please ensure that he specifed templaye exists.
VM(s) will be named using the `basename` supplied with a numerial suffix starting at `1`.
# Marzipan scripts

this subdirectory contains two helper scripts for the use of `marzipan` in automatically deploying a (cluster of) VM(s).

## `marzipan_deploy.py`
This is the `ENTRYPOINT` script for the Docker framework. It deals with the fully automated instantiation deployment and provisioning of the cluster, leveraging `marzipan` for the instantiation, ansible for futher set-up and `emma_marzipan` for the (software) provisioning. It also auto-generates ssh key pairs for the `root` and `ubuntu` users via the [`generate_keys.sh` script](#generate_keyssh).


## `generate_keys.sh`
This script generates pairs of SSH keys for the `ubuntu` and `root` users. and makes them available to be included in the VM(s) being spun up.
# Marzipan ansible playbooks

The `ansible_playbooks` folder in this subdirectory contains several ansible playbooks specifying further setup actions performed across all created VMs/nodes. Ansible provides a tool to script such provisioning actions for execution across distributed systems.

Of particular importance are the `firewall.yml`, `set_ssh_keys.yml`, and `update_hosts_file.yml`. However, no modifications are required.
The first defines firewall settings, the second copys ssh keys to all VMs enabling access, and the last creates/upates a file with all cluster hosts.
The hosts file is/can be used by subsequent ansible playbooks (e.g. those in `emma_marzipan`) to provision all VMs/nodes with desired software, as well as to start and stop cluster sevices such as, e.g. Dask.# VM template(s)

This subdirectory contains VM template configurations to be deployed on the OpenNebula platform.
An example is provided by `opennebula_goera.tpl`.

This template can and must be adapted according to the users needs. If renamed the path to the appropriate template must be set in the `ClusterConf.ini` config file.

Relevant sections that should be modified include
```
CONTEXT = [
  DNS_HOSTNAME = "YES",
  NETWORK = "YES",
  SSH_PUBLIC_KEY = "replace_root_key",
  USERNAME = "root",
  GROUP = "eratosthenes-uu"]
```
The group should be changed to the user's. __Note__ the `SSH_PUBLIC_KEY` is filled with a placeholder that filled by the automatically generated root public ssh key when using fully automatic deployment. Should the user wish deploy more manually this field can/should be modified.

```
DISK = [
  DATASTORE = "local_images_ssd",
  DATASTORE_ID = "104",
  IMAGE_ID = "25621",
  SIZE = "15360",
  TYPE = "fs"]
```
This specifies the OS disk. The user MUST ensure that a base image (e.g. and Ubuntu server) exists in the specified datastore of the OpenNebula Platform with the specified `IMAGE_ID`. The size corresponds to the OS partition.

```
DISK = [
  DATASTORE =  "ceph",
  DATASTORE_ID = "106",
  FORMAT = "raw",
  SIZE = "71680",
  TARGET = "vdb",
  TYPE = "fs" ]
```
This section specifies the local HDD of the VM.

```
DISK = [
  DATASTORE = "local_system_ssd",
  DATASTORE_ID = "103",
  DISK_TYPE = "FILE",
  SIZE = "49152",
  TARGET = "vdc",
  TYPE = "swap"]
```
This section defines a swap partition.


```
MEMORY = "32768"
MEMORY_UNIT_COST = "MB"
```
This specifies the VM's RAM allocation, anmd finally the following defines it's network connection. `NETWORK` and `NETWORK_UNAME` can be adpated, but this is not required.

```
NIC = [
  NETWORK = "internet",
  NETWORK_UNAME = "oneadmin" ]
```  

# Docker framework for marzipan

Docker file to create image and run container from which to perform instantiation and provisioning
# marzipan core module

The core `marzipan.py` module leverages the [PYONE](http://docs.opennebula.io/5.12/integration/system_interfaces/python.html) python bindings to the OpenNebula XML-RPC API to provide higer level abstrations for the steps required in instantiating a (cluster of) VM(s) on bare metal using the OpenNebula platform

`marzipan.py` provides the `one_interface` class with methods enabling the creation, deletion, instantiation, and termination  of templates and virtual machines.
In additon, the module provides a function `deploy_cluster()` which can be called and will automatically handle spinning up a cluster as defined in the `ClusterConf.ini` and `opennebula_goera.tpl` configuration and template files.

Finally `marzipan.py` can also be run as a script invoking the `deploy_cluster()` function.

Further documentation is provided inline

## References/Documentation

[OpenNebula 5.2 Documentation for the XML-RPC API](http://docs.opennebula.io/5.2/integration/system_interfaces/api.html#actions-for-templates-management)

[PYONE bindings documentation](http://docs.opennebula.io/5.12/integration/system_interfaces/python.html)
# Inventory builder

Python module to create and inventorize hosts (VMs/nodes) for different roles as created and provisioned by ansible playbooks (from `emma_marzipan`).
