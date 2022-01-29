# Minio

[Minio](https://www.minio.io/) is a distributed object storage server built for cloud applications and devops.
To use minio in distributed mode and have redundancy there are some pre-requisites. To understand them you should read the [distributed minio quickstart guide](https://docs.minio.io/docs/distributed-minio-quickstart-guide).
Before you set a minio cluster, make sure you set minio global variables using the template under *vars/*.
Once initialized a web GUI will be available at *http://${CLUSTER_NAME}0.${HOST_DOMAIN}:9091*, or any other host part of the *minio* group.

For unit tests please read the README under *roles/minio/tests/*.

# Ansible
Ansible is used to setup environment in a set of machines.

## Setup environment:
Every time the user opens a bash console on Windows or Linux the environment is set through the following commands:
```
#Windows
cd <path_to_emma>/emma
. env_linux.sh
./env_windows.sh

#Linux
cd <path_to_emma>/emma
. env_linux.sh
```
The user should check if each role has extra steps to setup environment. For example, the Hadoop role requires the generation of an extra ssh-key for the Hadoop user.

## Install Ansible
The recommended version is Ansible 2.3. However, we try to have the playbooks aligned with the latest version.

Ansible should then be install using Ubuntu package manager **apt-get**. To install the latest Ansible the user should follow the [installation instructions from the Ansible web-site](http://docs.ansible.com/ansible/latest/intro_installation.html#latest-releases-via-apt-ubuntu).

## Provision

Create the `hosts` file see `hosts.template` for template. To create a default four nodes cluster run the following command:
```
./create_hosts.sh
```

To change ansible configurations the user should edit ansible.cfg at the root directory of this repository. A diff between it and the file under /etc/ansible/ansible.cfg shows the additions to the default version. 

Now use ansible to verify login.
```
ansible all -u ubuntu -i hosts -m ping
```

For cloud based setup, skip this when deploying to vagrant. The disk (in example /dev/vdb) for /data/local can be partitioned/formatted/mounted (also sets ups ssh keys for root) with:
```
ansible-playbook -e datadisk=/dev/vdb -e host_name=$CLUSTER_NAME prepcloud-playbook.yml
```

If a apt is auto updating the playbook will fail. Use following commands to clean on the host:
```
kill <apt process id>
dpkg --configure -a
```

The roles defined for Ansible will create a platform with the following features:

* [GlusterFS](gluster.md)
* [Minio](minio.md)
* [Spark](spark.md) Standalone cluster
* [Docker Swarm](dockerswarm.md)
* [JupyterHub](jupyterhub.md)
* [Dask](roles/dask/README.md)

Each role will deploy the respective service or system on the node's group specified in **playbook.yml**. Each node's group is defined in the inventory file **hosts**. The **playbook.yml** defines on which order the roles are executed. 

The global variables for each role are defined under *vars/* in a file with the same name as the role. Their values should be set before running any
playbook. All templates contain the default values for each role's variable. Note: each variable defined as global will over-write the ones defined
under each role in the **defaults/** dir. Hence, to change or extend the a variable definition use the global variable definition. To create **vars'
files** with the default values, with exception of **minio_vars.yml** which needs to be edited to set the secret and access key, just run:
```
cd vars/
./create_vars_files.sh
```

Once all variables are defined the platform is installed with the following command:
```
ansible-playbook install_platform.yml
```

The platform only needs to be installed once. Once it is installed the services, e.g., Hadoop and Spark, are started using the following command:
```
ansible-playbook start_platform.yml
```

To shutdown the platform just run the following command:
```
ansible-playbook shutdown_platform.yml
```

Through [ansible-tags](http://docs.ansible.com/ansible/playbooks_tags.html) it is possible to have fine grained control over the task execution. To know all existent tags run:
```
ansible-playbook install_platform_light.yml --list-tags
```

Currently we have the following tags (if some tag is missing please fill in an issue):
* **common**: All tasks for role *common*.
* **extra_python_packages**: All tasks to install python packages using pip.
* **extra_system_packages**: All tasks to install system packages using apt-get.
* **firewall**: All tasks to update firewall

* **minio**: All tasks to install/start/stop services related with minio role.

* **hadoop**: All tasks to install/start/stop services related with hadoop role.

* **spark**: All tasks to install/start/stop services related with spark role.

* **jupyterhub**: All tasks to install/start/stop services related with jupyterhub role.
* **jupyter_modules**: All tasks to install extra modules for jupyterhub.
* **dask**: All tasks to install/start/stop services related with dask role.

If you wanted to just to update firewall instead of run the entire installation, you could do this:
```
ansible-playbook install_platform.yml --tags "firewall"
```
On the other hand, if you want to start the platform with exception of Minio service, you could do this:
```
ansible-playbook start_platform.yml --skip-tags "minio"
```

### Update an existing platform
In case a cluster is already installed and the user wants update the **Hadoop** or **Spark** cluster to an older or newer version the user should do the following:
```
#edit vars/hadoop_vars.yml and set
hadoop_prev_version: "<current_version>"
hadoop_version: "<new_version>"

#edit vars/spark_vars.yml and set
spark_prev_version: "<current_version>"
spark_version: "<new_version>"

#Run installation script just for Hadoop and Spark
ansible-playbook install_platform.yml --tags "hadoop,spark"

#In case the user does not want to format HDFS
ansible-playbook install_platform.yml --tags "hadoop,spark" --skip-tags "hdfs_format"
```

### Extend an existing platform
In case a Spark cluster is already defined and the user wants to add new nodes, it is not necessary to re-create the cluster. The user only has to install Spark and Hadoop (for HDFS) in the new nodes and add them to HDFS as *data nodes* and to Spark as *workers*. The first step is to create a new inventory file, simply make a copy of the [current one](https://github.com/nlesc-sherlock/emma/blob/master/ansible.md#provision) and add the new nodes address under **[hadoop-datanode]** and **[spark-worker]**, keep the ones under **[hadoop-namenode]** and **[spark-master]** and remove all other ones.
```
cd emma
cp hosts hosts_new
#edit hosts_new
vim hosts_new
```

The second step is to install all systems listed in the playbook `install_platform_light.yml`. Using `--tags` the user decides which systems should be installed.
```
#Shutdown cluster
ansible-playbook shutdown_platform.yml

#Create data partitions, only run this command if the nodes are not Vagrant boxes.
ansible-playbook -i hosts_new -e datadisk=/dev/vdb -e host_name=$CLUSTER_NAME prepcloud-playbook.yml

#Install Spark and Hadoop.
ansible-playbook -i hosts_new install_platform_light.yml --skip-tags "hdfs_format"
```

Once *Spark* and *Hadoop* is installed in the new nodes add their address, in the original inventory file (hosts), under **[hadoop-datanode]** and **[spark-worker]**. Once it is done, the user should start the cluster and check if they are listed as *Data nodes* (visit <**hadoop-namenode_url**>:50070/dfshealth.html#tab-datanode) and *Spark workers* ( visit <**spark-master_url**>:8080).
```
ansible-playbook start_platform.yml
```

A final step is to re-balance data in HDFS using **hdfs balancer**. It can be done at one of the cluster's nodes or remotely from the user's laptop. For the latter the user should check the instructions to [install Hadoop binaries and configure it](https://github.com/phenology/infrastructure/blob/master/platform/README.md#hadoop-binaries). The command to re-balance data is:
```
cd hadoop-2.8.1/bin
sudo ./hdfs balancer -threshold 5
```

### Add new modules
For a new application an user might need to install a new module, such as a python module, using either **pip* or **apt-get**. The library to be installed needs to be listed in **emma/vars/common_vars.yaml** either under the *python_packages* or *system_packages* variables. To install them the user only needs to run an Emma's Ansible playbook:
```
#Pip
ansible-playbook install_platform_light.yml --tags "extra_python_packages"

#Apt-get
ansible-playbook install_platform_light.yml --tags "extra_system_packages"
```

It is also possible to copy *user defined modules (UDM)* such as a **python library**, a **R library** and a **scala jar** to the cluster. To achieve that the user needs to copy the module to the respective folder **emma/files/<python | scala | r>** and call:
```
ansible-playbook install_platform_light.yml --tags "user_defined_modules"
```

The *UDM* is then available at the path **{{ jupyterhub_modules_dir }}/< python | scala | r>**, the default path is */data/local/jupyterhub/modules/*.

## Remote command execution
In case the user needs to run a specific command at one of the remote hosts, such as restarting the network service, the user should use **ansible**. Such option is interesting when the user simply wants to restart a service. However, if it is to install a system package the user should follow the steps in [**Add new module**](https://github.com/nlesc-sherlock/emma/blob/master/ansible.md#add-new-modules) and in case of a new system the user should create a new [*Ansible role*](http://docs.ansible.com/ansible/latest/playbooks_reuse_roles.html).

For fine-grain access the [**inventory file**](https://github.com/nlesc-sherlock/emma/blob/master/ansible.md#provision) (*often called hosts*) should have an entry per node. To extend the inventory file with these entries the user should run the following command, but only once:
```
. env_linux.sh
for f in `seq 0 $(calc $NUM_HOSTS-1)`; do echo $f; echo [$CLUSTER_NAME$f] >> hosts; echo $CLUSTER_NAME.$HOST_DOMAIN >> hosts; done
```

To execute a remote command the user needs to call *ansible <host_group> -a "<command>" -u ubuntu*. For example, if the user wants to restart the network service at host number zero (s)he should do the following:
```
ansible ${CLUSTER_NAME}0 -a "sudo systemctl restart networking.service" -u ubuntu
```


## Demo deployment

A demo deployment which uses the platform set by the above playbooks is done using the demos for the Sherlock project.
```
ansible-playbook demo.yml
```
Once deployed, a website is available on http://\<docker-swarm-manager\> (\<docker-swarm-manager\> is defined in the hosts file as described in [ansible.md](ansible.md).
# JupyterHub

JupyterHub connects with Spark using [toree](http://toree.apache.org).
# Spark

Spark is installed in `/data/shared/spark` directory as Spark Standalone mode.
* For master use `spark://<spark-master>:7077`
* The UI on http://<spark-master>:8080
* The JupyterHub on http://<spark-master>:8000

To get shell of Spark cluster run:
```
spark-shell --master spark://<spark-master>:7077
```

Before staring interaction with Spark, we recommend the read of [RDD, DataFrame and Data Set](https://indatalabs.com/blog/data-engineering/convert-spark-rdd-to-dataframe-dataset) and its [programming guide](http://spark.apache.org/docs/latest/programming-guide.html).

## SparkML
The installation followed the instructions for [SparkML using Spark 2.1.1](http://spark.apache.org/docs/latest/ml-guide.html).
All the dependencies should be installed. However, since Spark is a **pre-built tarball** the user should be aware of the [netlib-java dependency](http://spark.apache.org/docs/latest/ml-guide.html#dependencies) which requires an include.
```
com.github.fommil.netlib:all:1.1.2
```

## GeoTrellis
The installation followed the information available at [GeoTrellis GitHub](https://github.com/locationtech/geotrellis). The same repository has a **doc** dir where there is a [guide on how to use the library in Spark](https://github.com/locationtech/geotrellis/blob/master/docs/guide/spark.rst). When using sbt, maven or gradle the Geotrellis jars should be downloaded from [org.locationtech.geotrellis](https://mvnrepository.com/artifact/org.locationtech.geotrellis).

Before starting using GeoTrellis, we recommend the read of [Core concepts for Geo-trellis](https://geotrellis.readthedocs.io/en/1.0/guide/core-concepts/).

## Examples

A python notebook to do [Unsupervised classification of imagery using scikit-learn](http://nbviewer.jupyter.org/gist/om-henners/c6c8d40389dab75cf535). This example shows how to classify imagery (for example from LANDSAT) using scikit-learn. There are many classification methods available, but for this example they use K-Means as it's simple and fast. It uses [**numpy**](http://www.numpy.org/) and [**rasterio**](https://github.com/mapbox/rasterio).

## Debug mode
All information to debug set Spark for remote debugging and performance tuning.

### Remote debugging
To set Spark for remote debugging the user should set *spark_debug_mode* to **true** in **emma/vars/spark_vars/.yml** and reconfigure Spark to have only an executor per worker. In **emma/vars/spark_vars/.yml** the user should set *spark_executor_cores* equal to *spark_worker_cores* and *spark_executor_memory* equal to *spark_worker_memory*. Such setting allows us to open a single debugger per executor, i.e., one per node.

By default **driver's debugging port** is *5005* while the **worker's debugging port** is *5006*. They are defined in **emma/vars/spark_vars/.yml** by the variables *worker_debug_port* and *driver_debug_port*. To have either a worker or driver *waiting on startup* the variables *worker_waiting_on_startup* and *driver_waiting_on_startup* should be set to **y** (yes), by default they are set to **n** (no). 

For the configuration take place the user should restart Spark and Jupyterhub.
```
ansible-playbook shutdown_platform.yml --tags "spark,jupyterhub"
ansible-playbook start_platform.yml --tags "spark,jupyterhub"
```

### Tuning garbage collection
Spark is written in Scala, therefore, each process is a Java virtual machine (JVM). JVM memory management is done by a Garbage Collector (GC). The type of GC to be used is defined by variable *gc_type*. If the user wants to obtain GC debug information (s)he should set in **emma/vars/spark_vars/.yml** *gc_debug* to:
```
-XX:+PrintFlagsFinal -XX:+PrintReferenceGC -verbose:gc -XX:+PrintGCDetails -XX:+PrintGCTimeStamps -XX:+PrintAdaptiveSizePolicy
```

[Tuning Java Garbage Collection for Apache Spark Applications](https://databricks.com/blog/2015/05/28/tuning-java-garbage-collection-for-spark-applications.html) is a blog post which explains how memory management in a JVM is done and how we can tune GC for Spark applications.
# Hadoop
The platform uses Hadoop HDFS as a file system for Spark. For web-ui access to it use [the hadoop-name node address defined in the inventory file](https://github.com/nlesc-sherlock/emma/blob/master/ansible.md#provision) which is listening on port **50070**.

For setting up Hadoop on a cluster of machines, the master should be able to do a password-less ssh to start the daemons on all the slaves. Hence, the user should generate a ssh-key with the name **hadoop_id_rsa** and the key should be stored under the directory **files** located at the root directory.
```
cd files/
ssh-keygen -t rsa -f hadoop_id_rsa
```
# Vagrant

Vagrant is a tool to emulate a cluster of machines. Version **2.0.1** is recommended.

## Requirements

The user needs to have [VirtualBox](https://www.virtualbox.org/wiki/Downloads) installed. Version **5.2.4** is recommended.

## Installation

For Linux systems a simple package installation is enough.
```
#Ubuntu
wget https://releases.hashicorp.com/vagrant/2.0.1/vagrant_2.0.1_x86_64.deb
sudo dpkg -i vagrant_2.0.1_x86_64.deb
```

For Windows, despite the [Ubuntu environment](#windows) was set to run Ansible, vagrant needs to be installed as if it was to be executed using the CMD console. To install it download *msi* file from: https://www.vagrantup.com/downloads.html. Sometimes there are directories ownership issues with vagrant installation. To solve it is required to click in properties and claim ownership of the directory so the installation can proceed. Despite it is installed to be used on the CMD console vagrant.exe can be called from using [Ubuntu environment](#windows). Before doing that some environment variables need to be set. Create *env_linux.sh* and run *env_windows.sh* on [Ubuntu environment](#windows) before using *vagrant.exe*. It is important to make sure home directory for vagrant has write permissions for Windows users, not only for Windows root/administrator.
```
#create and edit env_linux.sh.template
cp env_linux.sh.template env_linux.cmd

#edit it
vim env_linux.cmd

#On Ubuntu bash for Windows run, it is required to restart all consoles to have the environment variables set.
./env_windows.sh
```

**NOTE:** despite the [Ubuntu environment](#windows) is used to run Vagrant the user should call the vagrant's command as this:
```
vagrant.exe
```

### Plugins
Vagrant needs two plugins and they will be installed in *VAGRANT\_HOME*.
```
#Plugin for persistent storage
vagrant(.exe) plugin install vagrant-persistent-storage

#Plugin to manage hosts
vagrant(.exe) plugin install vagrant-hostmanager

#It is recommended to install vbguest plugin
vagrant(.exe) plugin install vagrant-vbguest
```

## VMs management

On Windows to run Vagrant's commands simply use Ubuntu bash console.
To create the virtual machines or start them use:
```
vagrant(.exe) up
```

To update guest machines */etc/hosts* the user after a *vagrant up* should always run:
```
vagrant(.exe) hostmanager
```

On Linux the host machine */etc/hosts* will automatically be updated. On Windows because the [Ubuntu environment](#windows) has its own */etc/hosts* the guest nodes IPs need to be retrieved by hand. After *vagrant hostmanager* run:
```
sh getHosts.sh
```

To halt all VMs
```
vagrant(.exe) halt
```

To destroy all VMs
```
vagrant(.exe) destroy
```

In case vagrant needs to be set using a private network due to issues in getting IPs in the public network the option **--network-type=private_network** should be used.
```
vagrant(.exe) --network-type=private_network up
```

If not used, vagrant will set a public network by default. To switch between a public and private network and vice versa it is required a **vagrant(.exe) halt** and then **vagrant(.exe) up**, it is not recommended to use **vagrant(.exe) reload**.

## Check

Verify login for *N* hosts.
```
ssh -i ${CLUSTER_NAME}.key root@${CLUSTER_NAME}0.$HOST_DOMAIN uptime
...
ssh -i ${CLUSTER_NAME}.key root@${CLUSTER_NAME}N.$HOST_DOMAIN uptime
```
# Emma

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.996308.svg)](https://doi.org/10.5281/zenodo.996308)

Emma is a project to create a platform for development of application for Spark and DockerSwarm clusters. The platform runs on an infra-structure composed by virtual machines and Ansible playbooks are used to create a storage layer, processing layer and [JupyterHub](https://jupyter-notebook.readthedocs.io/en/latest/index.html) services. The storage layer offers two flavors of storage, file-base by [GlusterFS](https://www.gluster.org/) and [Hadoop Distributed File System (HDFS)](http://hadoop.apache.org/), and object-based using [Minio](https://www.minio.io). The processing layer has a [Apache Spark cluster](http://spark.apache.org/) and a [Docker Swarm](https://docs.docker.com/engine/swarm/) sharing the storage instances.

## Deployment
At the moment the deployment of clusters with **emma** was tested with two OSs environment, Linux and Windows.

### Linux as host OS
For Linux it was tested using Ubuntu 14.04 and Ubuntu 16.04, however, it is recommended to use the latter. For Linux no special environment setup is required.
Before deploying a cluster user needs to clone **emma** repository:
```
git clone https://github.com/nlesc-sherlock/emma
```

### Windows as host OS
For Windows, **emma** was only tested on Windows 10 using the embedded Ubuntu 16.04 environment. It setup is straight forward and it simple follows the steps listed in the [installation guide](https://msdn.microsoft.com/en-us/commandline/wsl/install_guide).

Once installed it recommended to verify if the version 16.04 is installed. For that the user simply needs to type **bash** in the search windows at the left bottom of Windows 10 Desktop-environment and press enter.
Once the bash console is open the user should then type:
```
lsb_release -a
```

If not 16.04, then do the following:
```
#Powershell as administrator, and enter the command
lxrun /uninstall /full .

lxrun /install /y

#Verify again the version:
lsb_release -a
```

After the installation the Ubuntu environment is accessible through the bash command of Windows.
To add Windows executables to your Ubuntu *$PATH*, add the following to your **~/.bashrc** in the bash console:
```
vim ~/bashrc
export PATH=$PATH:/mnt/c/Windows/System32/
```

Note the *C* drive will be mounted with the files owned by *root* and file permissions set to *777*.
This means ssh keys will to open for Ansible. Hence, before you run ansible you need to call getHosts.sh.

For windows users need to change the line-ending setting for **git**, they should configure **git** to *Checkout as-is, commit as-is*:
```
git config --global core.autocrlf false
```

Before deploying a cluster user needs to clone **emma** repository:
```
git clone https://github.com/nlesc-sherlock/emma
```

The Windows Subsystem Linux will mount the C: drive as `/mnt/c`. 
By default the file permissions will be read/write/executable for all, this is bad because to use ssh-keys stricter permissions are required.
To solve this create the `/etc/wsl.conf` file with the following content and restart the shell.
```ini
[automount]
enable=true
root=/mnt
options="metadata"
```

### Setup environment

Independently of which OS the user is using the following steps need to be done after the repository is cloned.
```
cd emma

#create and edit env_linux.sh
cp env_linux.sh.template env_linux.sh
vim env_linux.sh

#Linux environments (also in the embedded Ubuntu environment in Windows).
#On each bash
. env_linux.sh

# Key used by root, do not set passphrase when asked
ssh-keygen -f ${CLUSTER_NAME}.key
# Key used by ${CLUSTER_NAME} user
ssh-keygen -f files/${CLUSTER_NAME}.key
```

Every time the user opens a bash console on Windows or Linux the environment is set through the following commands:
```
#Windows
cd <path_to_emma>/emma
. env_linux.sh
./env_windows.sh

#Linux
cd <path_to_emma>/emma
. env_linux.sh
```

### Infra-structure

With the environment set, the next step is to setup the infra-structure. The infra-structure, physical place where the platform runs, is composed by a set of virtual machines with the following characteristics:
1. Ubuntu 16.04 OS
2. Public network interface
3. OS disk, 200Mb for software + enough room in /tmp
4. Passwordless login as root with `${CLUSTER_NAME}.key` private key.
5. XFS Partition mounted at /data/local (used for swapfile, GlusterFS brick, Docker root)
6. Python2 to run Ansible tasks


The infrastructure is a collection of machines which must be reachable by ssh. The machines must be prepared/constructed by either [preparing cloud virtual machine](cloud.md) or [constructing using Vagrant boxes](vagrant.md).
Once the machines are prepared the servers are provisioned using [Ansible](https://www.ansible.com/), an automation tool for IT infra-structure. The roles defined for Ansible will create a platform with the following features:

* [GlusterFS](glusterfs.md)
* [Minio](minio.md)
* [Hadoop](hadoop.md)
* [Spark](spark.md) Standalone cluster
* [Docker Swarm](dockerswarm.md)
* [JupyterHub](jupyterhub.md)

Preceding the platform's installation, the user should click on each feature to understand the setup requirements for each of them. Once all the requirements have been fulfilled, the user should follow the platform's installation steps listed in **[ansible.md](ansible.md)**.
# Cloud

Once the virtual machines are instantiated in the cloud give some time for their kernel get updated. For cloud based setup, skip this when deploying to vagrant. The disk (in example /dev/vdb) for /data/local can be partitioned/formatted/mounted (also sets ups ssh keys for root) with (make sure first you read [ansible.md](ansible.md) to have Ansbile up and running):
```
#First define the environment, use env_linux.sh.template. Then run
. env_linux.sh

ansible-playbook -e datadisk=/dev/vdb -e host_name=$CLUSTER_NAME prepcloud-playbook.yml
```

If the first run fails because the **apt-get update** fails the solution is to reboot the machines and then log in into each of them using the web-ui and do the following:
```
sudo dpkg --configure -a

#A graphical window will prompt and the first option should be selected, i.e., install the version from the package manager.
```

On [SURFsara HPC cloud each guest node gets an DNS entry](https://doc.hpccloud.surfsara.nl/access-your-VM) as ${vmname}.${projectname}.surf-hosted.nl. They can be used to access web-uis of installed systems such as JupyterHub.
# Docker swarm

All nodes have a Docker daemon running.

The Docker swarm endpoint is at `<docker_manager_ip>` IP address (Set in `hosts` file).
Howto see https://docs.docker.com/engine/swarm/swarm-tutorial/deploy-service/

To use Swarm login on `docker-swarm-manager` host as configured in `hosts` file.

# GlusterFS

See http://gluster.readthedocs.io
All nodes have a xfs partition which is available as `gv0` volume and mounted as /data/shared on all nodes.
The volume is configured (replicas/stripes/transport/etc) in `roles/glusterfs/tasks/create-volume.yml` file.

# Global file location

Location to store files specific to each user running different playbooks. One example are the ssk keys.
# Scala modules location

Location to store Scala modules specific to each user. To install them run:
```
ansible-playbook install_platform_light.yml --tags "user_defined_modules"
```
# Python modules location

Location to store Python modules specific to each user. To install them run:
```
ansible-playbook install_platform_light.yml --tags "user_defined_modules"
```
# R modules location

Location to store R modules specific to each user. To install them run:
```
ansible-playbook install_platform_light.yml --tags "user_defined_modules"
```
#Minio

To mount a local dir to a bucket in minio you need to mirror in two directions with option -w to watch for changes. The mirror should run in background.
```
#All these commands should be run in all hosts.

#register a alias in ~/.mc/config for a S3 storage service
mc config host add s3 <host_ip>:9091 <access_key> <secret_key> S3v4
example: mc config host add s3 http://145.100.116.245:9091 A24H1RIGV4RKFGXJTEMS 5jd7ARCOi/XVjLzXqT5wA1NSgjmUo9mYJBgyGyIh S3v4

#Lets get what is in the S3 my bucket into /data/shared/
mc mirror -w /data/shared s3/mybucket &

#Lets get what is in the S3 my bucket into /data/shared/
mc mirror -w s3/mybucket /data/shared/ &

To mirror sub-directories you need to mirror them with a different bucket. 
For example, to mirror /data/shared/scratch

mc mirror -w /data/shared/scratch s3/superbucket &
mc mirror -w s3/superbucket /data/shared/scratch/ &
```

Configuration file for s3cmd .s3cfg
```
host_base = 145.100.116.153:9091
host_bucket = 145.100.116.153:9091
access_key = A24H1RIGV4RKFGXJTEMS
secret_key = 5jd7ARCOi/XVjLzXqT5wA1NSgjmUo9mYJBgyGyIh
use_https = False
list_md5 = False
use_mime_magic = False
#Make sure the region is the same as the one used by minio
bucket_location = us-east-1

```

Example of commands:
```
s3cmd  ls s3://files

s3cmd get s3://files/sonnets.txt romulo.txt
```

To upload data to a sub-directory follow this example:
```
cd <root_dir>/<sub_dir> ; for f in `ls *`; do s3cmd put $f s3://<root_dir>/<sub_dir>/$f; done
```
For unit tests we will use [testinfra](https://github.com/philpep/testinfra).

To install dependencies and testinfra do the following: 
```
sudo apt-get install python-pip python-dev libffi-dev libssl-dev libxml2-dev libxslt1-dev libjpeg8-dev zlib1g-dev
sudo pip install testinfra
```

It is necessary to define ssh-config.
```
#On windows you need to replace the windows path by the Ubuntu environment path.
vagrant ssh-config > ./roles/minio/tests/ssh-config
```

To run a test to verify if minio is installed and up and running run the following:
```
testinfra --ansible-inventory=hosts --ssh-config=./roles/minio/tests/ssh-config --connection=ssh --sudo ./roles/minio/tests/test_minio.py
```
#Spark-shell
To connect spark-shell you need to do:
```
spark-shell
```

Download data and load it into hdfs:
```
#Create dir in HDFS, if you login as Ubuntu
hadoop dfs -mkdir /user/ubuntu/files

#Download sample_kmeans_data.txt file
wget https://raw.githubusercontent.com/apache/spark/master/data/mllib/sample_kmeans_data.txt

#Upload the sample_kmeans_data.txt to HDFS
hadoop dfs -put ./sample_kmeans_data.txt /user/ubuntu/files/sample_kmeans_data.txt
```
#Spark-shell
To connect spark-shell with drivers to read from S3 you need to do:
```
spark-shell
```

Download data and load it into S3:
```
#Create a bucket in S3
s3cmd mb s3://files

#Download sonnets.txt files
wget http://www.gutenberg.org/cache/epub/1041/pg1041.txt -o sonnets.txt

#Upload the sonnets.txt to S3 storage
s3cmd put sonnets.txt s3://files/
```
Dask
====

Shared [Dask](http://dask.pydata.org) cluster.

The Dask scheduler is run on the `dask-scheduler` host on port 9091.

To use in Python
```python
from dask.distributed import Client
client = Client('<dask-scheduler>:9091')
```

Requirements
------------

* Cluster is running as `ubuntu` posix user, so user should exist.
* The firewall should allow incoming connections to the Dask scheduler port on the scheduler host.

Role Variables
--------------

This role has the following variables:
```yaml
# Port on which dask scheduler is running
dask_scheduler_port: 9091
# Work directory of Dask scheduler
dask_scheduler_dir: "/home/ubuntu/"
# Work directory of Dask worker
dask_worker_dir: "/home/ubuntu/"
```

This role has the following host groups:
* `dask-scheduler`, Host on which Dask scheduler is run
* `dask-worker`, Hosts on which Dask workers are run

Dependencies
------------

No dependencies on other Ansible roles.

Example Playbook
----------------

Including an example of how to use your role (for instance, with variables passed in as parameters) is always nice for users too:

    - hosts: servers
      roles:
         - { role: username.rolename, x: 42 }

License
-------

Apache v2.0

Author Information
------------------

An optional section for the role authors to include contact information, or a website (HTML is not allowed).
#Example in on to manipulate GeoTiffs
The examples were obtained from [GeoTrellis documentation](http://geotrellis.readthedocs.io/en/latest/guide/spark.html#example-use-cases).

##Stitching Tiles into a single GeoTiff
This example show how to [stich tiles into a single GeoTiff](https://github.com/locationtech/geotrellis/blob/master/docs/guide/spark.rst#stiching-tiles-into-a-single-geotiff).
