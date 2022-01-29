# [DARE Platform](https://project-dare.gitlab.io/dare-platform/)

## Welcome 
The DARE platform is the integrated technical outcome of the [H2020 DARE project (#777413)
](http://project-dare.eu). The platform is documented in greater detail in its [dedicated GitLab microsite](https://project-dare.gitlab.io/dare-platform/).

The DARE platform is designed to address the requirements of research engineers and scientists, and empowers them to overcome technical and collaboration hurdles. It achieves these goals by exploiting the elasticity of the Cloud in resource acquisition and allocation. Crucially, it combines data provenance with other metadata and registries regarding datasets, experiments and resources to best enact workflows on present Cloud resources but also on remote and diverse platforms.

Using the DARE platform, users can design, implement and make available scientific workflows. They can subsequently execute them via the DARE API in order to power user-facing products, such as Virtual Research Environments, Web portals, etc. These workflows and their constituent processing elements are registered and available on the DARE repositories. These repositories can be used to find, communicate, replicate or debug existing methods. During execution DARE offers tunable provenance, which can help research engineers debug methods as well as site administrators locate bottlenecks and optimise local installations.

Whilst the DARE platform is being used in a number of use-cases, it is still a research product, and the DARE gitlab group an active collaboration tool between consortium members. As such, not all DARE repositories are being used as part of the deployable DARE platform. The ones used and their addresses are clearly listed below. 

DARE is a platform that should be deployed on a [Kubernetes](https://kubernetes.io)-enabled resource. In an operational setting it would require an administrator for deployment and day-to-day maintenance. Interactive use with the platform is typically provided via a [Jupyter Hub](https://jupyter.org). E.g. it could:

* support a Web application on a commercial cloud, or 
* it could form the basis of an Open Research solution at a research or academic site, or
* a combination of the above.

### More info
More information on the DARE platform, its components and how they relate, as well as its overall architecture and direction can be found in the following papers:

* I. Klampanos et al., "DARE: A Reflective Platform Designed to Enable Agile Data-Driven Research on the Cloud," 2019 15th International Conference on eScience (eScience), San Diego, CA, USA, 2019, pp. 578-585, doi: 10.1109/eScience.2019.00079. [[Read on IEEE Xplore]](https://ieeexplore.ieee.org/document/9041753)
 
* M. Atkinson et al., "Comprehensible Control for Researchers and Developers Facing Data Challenges," 2019 15th International Conference on eScience (eScience), San Diego, CA, USA, 2019, pp. 311-320, doi: 10.1109/eScience.2019.00042. [[Read on IEEE Xplore]](https://ieeexplore.ieee.org/document/9041709)


### Contributing to the DARE platform

The DARE platform and its components are published on an [Apache 2.0 License](https://www.apache.org/licenses/LICENSE-2.0). Everyone is welcome to download, deploy, and modify the source code, as well as to propose bug fixes and changes, either by creating issues or by contributing source code.
The most straightforward way to contribute code to the DARE platform and to its component repositories is by working on a fork and creating a pull request. 


## Deployment Instructions

The DARE platform consists of additional components, present in additional repositories in the DARE GitLab group. 

It is assumed that these components live in the same Kubernetes-managed cluster. Some components expose REST APIs to the outside world, while for some this is achieved via the DARE Execution API (see below). With some exceptions, the components are decoupled, with the DARE-platform and the Execution API repositories providing the glue. For instance, one could download and use dispel4py locally, but the Execution API makes it available as a service and manages a self-spawning MPI cluster on the Cloud.

The DARE components are documented in our [GitLab microsite](https://project-dare.gitlab.io/dare-platform/api-documentation/) with
instructions on how to use the DARE API. Additionally, there are available [installation instructions](https://project-dare.gitlab.io/dare-platform/configuration/).

In order to perform and integration test to the platform, check our toy example [here](https://project-dare.gitlab.io/dare-platform/demo/)

## Components

The DARE platform consists of a number of components developed by the partners of the DARE project, as well as of 3rd-party components. The DARE platform integrates the following internally-provided components.

### dispel4py
dispel4py is a Python library for describing abstract stream-based workflows for distributed data-intensive applications. It enables users to focus on their scientific methods, avoiding distracting details and retaining flexibility over the computing infrastructure they use.

It delivers mappings to diverse computing infrastructures, including cloud technologies, HPC architectures and  specialised data-intensive machines, to move seamlessly into production with large-scale data loads. More information can be found at the [dispel4py repository](https://gitlab.com/project-dare/dispel4py).

### s-ProvFlow
s-ProvFlow implements the P4 aspects of the DARE platform. It is a provenance framework for storage and access of data-intensive streaming lineage. It offers a a web API and a range of dedicated visualisation tools based on the underlying provenance model, S-PROV, which utilises and extends PROV and ProvONE models.

S-PROV addresses aspects of mapping between logical representation and concrete implementation of a workflow until its enactment onto a target computational resource.  The model captures aspects associated with the distribution of the computation, runtime changes and support for flexible metadata management and discovery for the data products generated by the execution of a data-intensive workflow.

Complete Documentation for the component can be found at the [s-ProvFlow repository](https://gitlab.com/project-dare/s-ProvFlow).

### dispel4py Registry
The dispel4py Registry is a RESTful Web service providing functionality for registering workflow entities, such as processing elements (PEs), functions and literals, while encouraging sharing and collaboration via groups and workspaces. More information is provided in the [dispel4py Registry repository](https://gitlab.com/project-dare/d4p-registry).

### CWL Workflow Registry
The CWL Workflow Registry provides a similar functionality as the Dispel4py Registry, with the difference that it is associated with CWL workflows. More information is provided at the [CWL workflow registry repository](https://gitlab.com/project-dare/workflow-registry).

### Data catalogue - Semantic Data Discovery
The Data catalogue is part of the DARE Knowledge Base and manages information related to the data elements processed via the platform. It exposes a RESTful API for registering new data sources and retrieving information on data previously registered or data results produced by processes executed over DARE. For
more information visit the [Data catalogue repository](https://gitlab.com/project-dare/semantic-data-discovery).

### DARE Execution API
The DARE Execution API enables the distributed and scalable execution of numerical simulations (now using SPECFEM3D code) and dispel4py workflows (e.g. used to describe the steps of RA except for simulations), which can be extended to other execution contexts. The Execution API also offers services such as uploading/downloading and referencing of data and process monitoring. For more information check out the [Execution API repository](https://gitlab.com/project-dare/exec-api).

### DARE playground 
The purpose of the playground is to provide an environment for testing and debugging purposes, especially dispel4py workflows. This helps users debug their methods before making them available for exection on the platform. More information is provided in the [DARE playground repository](https://gitlab.com/project-dare/playground).

---
title: 'DARE Platform\: a Developer-Friendly and Self-Optimising Workflows-as-a-Service Framework for e-Science on the Cloud'
tags:
  - Python
  - CWL
  - Worfklows-as-a-Service
  - Kubernetes
  - Cloud
  - e-infrastructures
authors:
  - name: Iraklis A. Klampanos^[Corresponding author.]
    orcid: 0000-0003-0478-4300
    affiliation: "1"
  - name: Chrysoula Themeli
    orcid: 0000-0002-6759-4136
    affiliation: "1"
  - name: "Alessandro Spinuso"
    affiliation: "2"
  - name: "Rosa Filgueira"
    affiliation: "3"
  - name: "Malcolm Atkinson"
    affiliation: "3"
  - name: "André Gemünd"
    affiliation: "4"
  - name: Vangelis Karkaletsis
    affiliation: "1"
affiliations:
 - name: National Centre for Scientific Research “Demokritos”, Greece
   index: 1
 - name: Koninklijk Nederlands Meteorologisch Instituut, the Netherlands
   index: 2
 - name: The University of Edinburgh, UK
   index: 3
 - name: Fraunhofer-Institut für Algorithmen und Wissenschaftliches Rechnen (SCAI), Germany
   index: 4
date: 24 July 2020
bibliography: bibliography.bib
---

# Statement of need

In recent years, science has relied more than ever on large-scale data as well as on distributed computing and human resources. Scientists and research engineers in fields such as climate science and computational seismology, constantly strive to make good use of remote and largely heterogeneous computing resources (HPC, Cloud, institutional or local resources, etc.), process, archive and analyse results stored in different locations and collaborate effectively with other scientists.

The DARE platform enables the seamless development and reusability of scientific workflows and applications, and the reproducibility of the experiments. Further, it provides Workflow-as-a-Service (WaaS) functionality and dynamic loading of execution contexts in order to hide technical complexity from its end users. This paper introduces the software implementing the DARE platform. More information on the H2020 DARE project is provided in @klampanos2019dare, @atkinson2019comprehensible, and  @atkinson_malcolm_2020_3697898.


# The DARE platform

The DARE platform is designed to live in-between user applications and the underlying computing resources. It is built on top of containerisation as well as parallelisation technologies, e.g., Kubernetes and MPI. Interfacing with client systems and end-users is achieved via RESTful APIs. The execution of scientific workflows is achieved via a Workflows-as-a-service layer, which can handle workflows described in either the dispel4py Python library [@Filgueira2017], or in the Common Workflow Language (CWL) [@amstutz2016common].

## The software

The DARE platform consists of a number of largely independent software components developed by the partners of the DARE project. All core software components are provided via the [DARE GitLab group](https://gitlab.com/project-dare). The [DARE Platform repository](https://gitlab.com/project-dare/dare-platform) provides pointers to all relevant repositories, documentation and more. 
Installation instructions and API documentation are provided in a [GitLab page](https://project-dare.gitlab.io/dare-platform/). A demo is available in the [DARE Execution API GitLab Repository](https://gitlab.com/project-dare/exec-api/-/tree/master/examples/mySplitMerge), which can also be used as an integration test.

The DARE platform and its components are published with the Apache 2.0 License. Everyone is welcome to download, deploy, and modify the source code, as well as to propose bug fixes and changes, either by creating issues or by contributing source code.
The most straightforward way to contribute code to the DARE platform and to its component repositories is by working on a fork and creating a pull request.

The core DARE platform components are the following:

### dispel4py
dispel4py is a Python library for describing abstract stream-based workflows for distributed data-intensive applications. It can translate higher-level workflows to diverse computing contexts, such as Apache Storm, MPI and plain shared-memory multi-core, to enable moving seamlessly into production with large-scale data loads. More information can be found at the [dispel4py repository](https://gitlab.com/project-dare/dispel4py).

### s-ProvFlow
s-ProvFlow implements a provenance framework for storage and access of data-intensive streaming lineage. It offers a web API and a range of dedicated visualisation tools based on the underlying provenance model, S-PROV, which utilises and extends PROV and ProvONE models.
Complete documentation for this component can be found at the [s-ProvFlow repository](https://gitlab.com/project-dare/s-ProvFlow).

### dispel4py Registry
The dispel4py Registry is a RESTful Web service providing functionality for registering workflow entities, such as processing elements (PEs), functions and literals, while encouraging sharing and collaboration via groups and workspaces. More information is provided in the [dispel4py Registry repository](https://gitlab.com/project-dare/d4p-registry).

### CWL Workflow Registry
The CWL Workflow Registry provides a similar functionality as the Dispel4py Registry, with the difference that it is associated with CWL workflows. More information is provided at the [CWL workflow registry repository](https://gitlab.com/project-dare/workflow-registry).


### DARE Execution API
The DARE Execution API enables the distributed and scalable execution of dispel4py and CWL workflows, and is extensible to other contexts. The Execution API also offers services such as uploading/downloading and referencing of data and process monitoring. More information is provided in the [Execution API repository](https://gitlab.com/project-dare/exec-api).

### DARE playground 
The purpose of the playground is to provide an environment for testing and debugging purposes, especially dispel4py workflows. This helps users debug their methods before making them available for execution on the platform. More information is provided in the [DARE playground repository](https://gitlab.com/project-dare/playground).


## Characteristics of the DARE platform

1. It interfaces with users and external systems via a comprehensive RESTful API.
2. It facilitates the development of modular, reusable and shareable data-intensive solutions.
3. It combines two different workflow approaches, dispel4py and CWL, within the same platform and development environment.
4. Via its execution API, it orchestrates the dynamic spawning and closing of MPI clusters on the cloud for MPI-enabled components.
5. It provides a flexible environment, which local administrators can parametrise, by supporting custom docker-based environments and user interfaces.
6. It supports the collection, mining and visualisation of provenance information.




# DARE platform use cases

The DARE platform is currently used in the following domain applications:

1. Seismology: [Rapid Assessment (RA) of ground motion parameters during large earthquakes](https://gitlab.com/project-dare/WP6_EPOS).
2. Seismology: [Moment Tensor 3D (MT3D) for ensemble-type of seismic modelling](https://gitlab.com/project-dare/WP6_EPOS).
3. Volcanology: [Ash fall hazard modelling](https://gitlab.com/project-dare/wp6_volcanology).
4. Climate-change: [Extending Climate4Impact with efficient and transparent access to diverse computing resources](https://gitlab.com/project-dare/WP7_IS-ENES_Climate4Impact).
5. Atmospheric sciences: [Cyclone tracking and visualisation application](https://gitlab.com/project-dare/wp7_cyclone-tracking).
 

# State of the field

The DARE platform implements research coming from multiple areas. This section is therefore not meant to be exhaustive but rather to provide basic state-of-the-field information for further study.
The need for unifying underlying e-infrastructures and platforms via higher-level interfaces, programmatic or interactive, is especially pronounced in Europe due to the widespread policy and technological diversity.  Generic technological solutions, such as the ones produced by the [COLA](https://project-cola.eu) project [@cola-kiss], move towards providing unifying low-level views of underlying infrastructures. However, to raise the level of abstraction for researchers also requires automation powered by tighter integration of heterogeneous components. Much of this functionality is powered by shared catalogues within and outside proposed technological solutions. 

Using shared catalogues as a basis for integration is central to projects, such as [VRE4EIC project](https://vre4eic.ercim.eu), which has developed research environments for collaborating research communities [@Martin2019]. Similar to DARE, the [SWITCH project](https://cordis.europa.eu/project/id/643963) has demonstrated using knowledge-bases for supporting enactment-target selection, optimisation, mapping and coping with heterogeneity [@Stefanic2019], focusing on time-critical applications. 

In terms of leveraging the Cloud paradigm to raise the abstraction level, the project [DEEP-Hybrid-DataCloud](https://deep-hybrid-datacloud.eu/)  makes use of underlying data representation and transformation functionality to provide machine learning as a service to a variety of target user groups [@deep-lopez-joss]. DEEP focuses on the exposure of computational resources, e.g. GPU clusters over federated Clouds. The [PROCESS project](https://www.process-project.eu/) has built a set of services and tools to enable extreme scale data processing in  scientific and advanced industry settings. Similar to the DARE platform, PROCESS offers a set of composable services covering from data processing to workflow specification and enactment. However DARE places more weight on supporting reflection via catalogues and registries to aid automation and optimisation.


# Acknowledgements

This work has been supported by the EU H2020 research and innovation programme under grant agreement No 777413.

# References
# joss-dare-platform---
title: "Architecture"
---

![Platform Architecture](../images/dare_architecture.png)

# DARE Components

### dispel4py
Dispel4py is a free and open-source Python library for describing abstract stream-based workflows for distributed data-intensive applications. It enables users to focus on their scientific methods, avoiding distracting details and retaining flexibility over the computing infrastructure they use. It delivers mappings to diverse computing infrastructures, including cloud technologies, HPC architectures and specialised data-intensive machines, to move seamlessly into production with large-scale data loads. The dispel4py system maps workflows dynamically onto multiple enactment systems, such as MPI, STORM and Multiprocessing, without users having to modify their workflows.

More information on dispel4p:

* [dispel4py documentation on pythonhosted.org](https://pythonhosted.org/dispel4py/)
* [Filgueira Rosa and others. dispel4py: A Python framework for data-intensive scientific computing, The International Journal of High Performance Computing Applications 2017, Vol. 31(4) 316–334](https://journals.sagepub.com/doi/pdf/10.1177/1094342016649766?casa_token=O03gkqs1Um8AAAAA:Mdg0KZoyTH4qIUQfw9jHzWFiu-T5ifgwxVGr1Nj_4X4_iEKI8GyyJGShr81pTS4T4wUD8TbJU9A)

### s-ProvFlow
s-ProvFlow implements the P4 aspects of the DARE platform. It is a provenance framework for storage and access of data-intensive streaming lineage. It offers a a web API and a range of dedicated visualisation tools based on the underlying provenance model, S-PROV, which utilises and extends PROV and ProvONE models.

S-PROV addresses aspects of mapping between logical representation and concrete implementation of a workflow until its enactment onto a target computational resource.  The model captures aspects associated with the distribution of the computation, runtime changes and support for flexible metadata management and discovery for the data products generated by the execution of a data-intensive workflow.

Complete Documentation for the component can be found at the [relevant repository](https://gitlab.com/project-dare/s-ProvFlow/tree/master).

### Dispel4py Workflow Registry 

The dispel4py Registry is a RESTful Web service providing functionality for registering workflow entities, such as processing elements (PEs), functions and literals, while encouraging sharing and collaboration via groups and workspaces.

The DARE users should register their workflows in the dispel4py registry before accessing the DARE Execution API in order to run them. Once the workflows are registered, they can be retrieved by name. Each workflow is uniquely identified by its workspace, package and PE name. For more information on the API, check the [Documentation](api).

### CWL Workflow Registry

Similarly to the dispel4py workflows, CWLs should also be registered in the respective registry in order to be accessible for execution. This component allows the registration of docker execution environments, which are associated with CWL worklfows. 

The platform admins should create and build docker environments and then register them in the CWL Workflow Registry. Afterwards, the research developers can list the existing dockers in the Registry, download their files etc in order to find a suitable environment for their application.

Once they have found an execution environment, they can register their application and associate it with a docker. After the registration, the workflows can be retrieved with their name and version. The DARE execution API will request only this information (i.e. name and version) in order to retrieve the workflow and execute it. If the application used inside the CWL workflow supports MPI, the users can request their docker to run in multiple containers.

### DARE Execution API

The DARE Execution API enables the distributed and scalable execution of numerical simulations (now using SPECFEM3D code), dispel4py workflows (e.g. used to describe the steps of RA except for simulations), which can be extended to other execution contexts, and CWL workflow (e.g. used in the Cyclone use case). Execution API also offers services such as uploading/downloading and referencing of data and process monitoring.
More information can be found at the relevant [repository page](https://gitlab.com/project-dare/exec-api).

### Data catalogue - Semantic Data Discovery
The Data catalogue is part of the DARE Knowledge Base and manages information related to the data elements processed via the platform. It exposes a RESTful API for registering new data sources and retrieving information on data previously registered or data results produced by processes executed over DARE. For
more information visit the [respective repository](https://gitlab.com/project-dare/semantic-data-discovery)

### Testing environment - Playground

DARE platform provides a "playground" - testing environment to research developers. During the testing phase of the workflow development, users can simulate a dispel4py workflow execution as well as simulate a local dispel4py run using the playground module of the DARE platform. Additional information and API description is available in the [corresponding repository](https://gitlab.com/project-dare/playground)
---
title: "Demo"
---
# Toy Example

Let's discuss on how to interact with the DARE Platform!

Before trying this toy example, we recommend to read the features and API documentation sections in order to have an idea of the components
that constitute the DARE Platform.

DARE platform exposes a RESTful API for workflow execution to research developers, in other words provides workflow execution as a service.
Let's discuss the necessary steps that our users need to implement. 

First of all, users should prepare their python script with their dispel4py workflow, that makes also use of provenance! 
Before the official execution, the DARE platform provides a testing environment for workflow development and testing! 
Research developers can test their code directly inside the platform and make the necessary adjustments using our playground module.

Once the workflow is ready for execution, the research developers can execute their workflow using the official execution API. 
In order to interact with the DARE platform, we provide a script named [helper_manager.py](https://gitlab.com/project-dare/exec-api/-/blob/master/client/helper_manager.py),
which wraps the necessary API calls to the DARE platform. In the below example code, we make use of this script in order to authenticate a user,
create workspace & register workflows, execute a workflow and monitor the execution etc.

The necessary steps to execute a workflow are listed below.

1) Authentication

```python
BASE_URL = "https://platform.dare.scai.fraunhofer.de/"

import json
from os import getcwd
from os.path import join, exists

import requests

# Download the DARE platform client - helper function library
hf_scripts = requests.get("https://gitlab.com/project-dare/exec-api/-/raw/master/client/helper_manager.py")
if hf_scripts.status_code == 200:
    with open("helper_manager.py", "w") as f:
        f.write(hf_scripts.text)
from helper_manager import DareManager

# Running notebook locally

credentials_file = "credentials.yaml"
if not exists(credentials_file):
    credentials_yaml = requests.get("https://gitlab.com/project-dare/exec-api/-/raw/master/client/example_credentials.yaml")

    if credentials_yaml.status_code == 200:
        with open("credentials.yaml", "w") as f:
            f.write(credentials_yaml.text)

# if you are working locally and not in the JupyterHub, pass the parameter:
# config_file="credentials.yaml
dm = DareManager(dare_platform_url=BASE_URL, config_file=credentials_file)
# only if you are working locally
print(dm.token)

# Running notebook in DARE platform's JupyterHub
dm = DareManager(dare_platform_url=BASE_URL)
print(dm.token)

```
2) D4p Information Registry: create a workspace and register your workflow

```python
code = requests.get('https://gitlab.com/project-dare/exec-api/-/raw/master/examples/mySplitMerge/scripts/mySplitMerge_prov.py')
code = str(code.text)

# TODO provide a name for your workflow
name = "mySplitMerge"
workspace_id, impl_id = dm.register_d4p_workflow(name=name, code=code)
print("Your workspace ID is: {}".format(workspace_id))
print("Your PE ID is: {}".format(impl_id))
```

3) Use execution API

a) Execute a workflow

```python
dm.exec_d4p(nodes=6, no_processes=6, iterations=1, 
            reqs='https://gitlab.com/project-dare/exec-api/-/raw/master/examples/mySplitMerge/scripts/reqs.txt')
```

b) Monitor the execution

```python
dm.monitor_job()
```

c) Upload files in the platform

```python
remote_path = "d4p-input"
filename = "input.json"
dm.upload_file(remote_path=remote_path, filename=filename)
```
d) List your folders & files

```python
dm.list_workspace()
dm.list_exec_folder()
```

e) Download a file from the platform

```python
dm.download_file(filename="logs.txt")
```

f) Share files in B2DROP

```python
# Case 1 file
kind = "file"
dare_path_kind = "run" # you can also use upload if you want to use some file from the uploads directory
# use the dare_directory parameter if you want a file from a different run directory than the one stored in the session
filename = ""
remote_dir_name = None
dm.b2drop_share(kind=kind, filename=filename, dare_path_kind=dare_path_kind, remote_dir=remote_dir_name)

# Case 2 directory
kind = "directory"
dare_path_kind = "run"
remote_dir_name = None
# Again use the dare_directory if you want to upload another run_dir and not the one in the session
dm.b2drop_share(kind=kind, dare_path_kind=dare_path_kind, remote_dir=remote_dir_name)
```

For hands-on practice with the DARE platform, for both Dispel4py and CWL workflow, we provide a tutorial jupyter notebook
in order to get familiar with the platform. We include the above examples as well as additional material, which can
be found in our [GitLab repository](https://gitlab.com/project-dare/exec-api/-/tree/master/examples/tutorial). You can download
the tutorial folder and use the Jupyter Notebook to interact with the platform. 

We provide client-side helper functions in order to make easier the interaction with the platform. You can find the relevant documentation
[here](https://project-dare.gitlab.io/exec-api/client.html).

Contact our team to request an account to our JupyterLab to access our demos and tutorials directly into the platform!
---
title: DARE platform January 2020 Release
date: 2020-01-14
image: "images/dare_transparent.png"
---

DARE platform January 2020 release is now available at https://testbed.project-dare.eu/

A new testing environment is now included in the platform. Use our toy example to get familiar with the new functionalities.

Additionally, the DARE platform's shared file system is re-organized and more user-oriented. 
Users can view/download their generated files by interacting with the DARE's execution API. 
Read more information check the Features section. To download the v2.1 release of the DARE platform
visit our [GitLab Repository](https://gitlab.com/project-dare/dare-platform/tree/v2.1)
 
---
title: DARE platform February 2020 Release
date: 2020-02-07
image: "images/dare_transparent.png"
---

DARE platform February 2020 release is now available at https://testbed.project-dare.eu/

In this version a queue service added to unburden the workflowexecutions/insert sprov-api end point.
Additionally, CWL provenance to s-prov mapping has been improved and sprovflow-viewer search options 
have been expanded. You can now add more search terms in the filter dialog.

To download the v2.2 release of the DARE platform visit our 
[GitLab Repository](https://gitlab.com/project-dare/dare-platform/-/tree/v2.2)
 
---
title: DARE platform v3.6 release
date: 2020-10-23
image: "images/dare_transparent.png"
---

DARE platform v3.6 release is now available at https://platform.dare.scai.fraunhofer.de/

The v3.6 release contains various fixes and updates in the Execution API component. It's written almost from scratch,
has improved features and easier API. The two Shared File Systems are now integrated into a single one and all the DARE
use cases are updated and follow the same folder structure.

To download the v3.6 Release of the DARE platform visit our 
[GitLab Repository](https://gitlab.com/project-dare/dare-platform/-/tree/v3.6)
 
---
title: DARE platform v3.3 release
date: 2020-06-24
image: "images/dare_transparent.png"
---

DARE platform v3.3 release is now available at https://platform.dare.scai.fraunhofer.de/

Thanks to new use cases, we fixed multiple issues, especially regarding the CWL support. Visit our toy example [here](https://project-dare.gitlab.io/dare-platform/demo/)

To download the v3.3 Release of the DARE platform visit our 
[GitLab Repository](https://gitlab.com/project-dare/dare-platform/-/tree/v3.3)
 
---
title: DARE platform v3.1 release
date: 2020-05-04
image: "images/dare_transparent.png"
---

DARE platform v3.1 release is now available at https://testbed.project-dare.eu/

From this version, apart from dispel4py, we also support CWL workflows.
We have a new component, deployed in the platform, i.e. the workflow-registry. Domain developers
can register their docker containers and CWL workflows in the CWL Workflow Registry and refer to them
by name and version. For more information, check the technical documentation of the Workflow registry
[here](https://project-dare.gitlab.io/workflow-registry/). 

We provide client-side functions 
[here](https://gitlab.com/project-dare/workflow-registry/-/blob/master/workflow_client/helper_functions.py),
which wrap all the endpoints provided by the workflow-registry component. You can use these functions
to store, update, retrieve and delete dockers and workflows. Download and modify our demo jupyter notebook
from [here](https://gitlab.com/project-dare/workflow-registry/-/tree/master/workflow_client) in order to 
register your own dockers and workflows.

For the CWL execution, research developers can use the /run-cwl endpoint of the Execution API,
specifying the name and version of the workflow. Thus, the DARE platform will deploy your
execution environment and run your workflow. As usual, the logs and output files will be stored in the
platform's Shared File System.

To download the v3.1 Release of the DARE platform visit our 
[GitLab Repository](https://gitlab.com/project-dare/dare-platform/-/tree/v3.1)
 
---
title: DARE platform 3rd major Release
date: 2020-04-03
image: "images/dare_transparent.png"
---

DARE platform v3.0 release is now available at https://testbed.project-dare.eu/

In this version, there are changes in almost all major platform components. First of all,
the new testing environment (playground module) is updated to run in an nginx server
(instead of gunicorn) and is completely integrated with the AAI. In a similar way,
the Execution API is now based on nginx image and it is also integrated with the AAI.
The Provenance component have also some important changes, i.e. a more resilient queue,
unit tests on ProvenanceStore class, new search expression on API and in GUI. Finally,
the dispel4py library includes the following changes:

* The -f and -d arguments can be used to define the workflow inputs as shown in the workflow provenance. 
The -f argument accepts a path to a file containing the input dataset in JSON format. 
The -d argument accepts the input dataset in JSON format.
* The -f argument has priority over the -d argument.
* Fixed the issue inline provenance definition in the workflow script is used to create the workflow provenance when the argument --provenance-config is present. 
The --provenance-config has priority over the inline provenance definition now.
* Attention: "s-prov:WFExecutionsInputs" in --provenance-config is deprecated.

Finally, we have a completely new component integrated in the platform, the Execution Registry.
This component serves as backend to the Execution API handling the Shared File System, and more specifically,
the experiments & runs of the users as well as their uploads.

To download the v3.0 Release of the DARE platform visit our 
[GitLab Repository](https://gitlab.com/project-dare/dare-platform/-/tree/v3.0)
 
---
title: DARE platform v3.2 release
date: 2020-05-25
image: "images/dare_transparent.png"
---

DARE platform v3.2 release is now available at https://testbed.project-dare.eu/

In the v3.2 release, we have completely integrated the Keycloak AAI to the Dispel4py Information Registry and performed
various updates in our Execution API. We have also deployed a new component, as a central login point, which handles
all the sign in actions and validates the provided access tokens. The dare-login component uses Keycloak as backend
to issue tokens and validate the users. The rest DARE components use the dare-login for all AAI actions.

The new DARE Login component is available [here](https://gitlab.com/project-dare/dare-login). As usual, we provide
client-side helper functions to interact with the component as well as a GitLab page for technical documentation.

To download the v3.2 Release of the DARE platform visit our 
[GitLab Repository](https://gitlab.com/project-dare/dare-platform/-/tree/v3.2)
 
---
title: DARE platform October 2019 Release
date: 2019-10-21
image: "images/dare_transparent.png"
---

DARE platform October 2019 release is now available at https://testbed.project-dare.eu/. 

The DARE Execution API is now integrated in the DARE platform. All the necessary files for
Docker and Kubernetes are added in the new DARE platform 
release in our [GitLab Repository](https://gitlab.com/project-dare/dare-platform/tree/v2.0).

Visit our [toy example](https://project-dare.gitlab.io/dare-platform/demo/) to interact with the DARE platform!
---
title: "Latest releases"
subtitle: ""
# meta description
description: ""
draft: false
---
---
title: DARE platform v3.4 release
date: 2020-07-27
image: "images/dare_transparent.png"
---

DARE platform v3.4 release is now available at https://platform.dare.scai.fraunhofer.de/

The v3.4 release contains various fixes and updates in workflow execution and monitoring.
The dare-login and exec-api components are improved to handle much more requests at the 
same time during the monitoring of the workflow execution.

To download the v3.4 Release of the DARE platform visit our 
[GitLab Repository](https://gitlab.com/project-dare/dare-platform/-/tree/v3.4)
 
---
title: DARE platform v3.5 release
date: 2020-10-14
image: "images/dare_transparent.png"
---

DARE platform v3.5 release is now available at https://platform.dare.scai.fraunhofer.de/

The v3.5 release contains various fixes and updates in the provenance component for the CWL workflows.
It also contains fixes in the Execution API and its helper functions as well as integration tests for
the components' APIs.

To download the v3.5 Release of the DARE platform visit our 
[GitLab Repository](https://gitlab.com/project-dare/dare-platform/-/tree/v3.5)
 
---
title: "Installation"
---
# DARE Deployment instructions

* Deployment on ubuntu

## Kubernetes setup

```bash
    # Install docker 
	sudo apt install docker.io
	
    # Enable docker 
	sudo systemctl enable docker

    # Add Kubernetes signing key 
	curl -s https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key add
	
    # Add Xenial Kubernetes Repository 
	sudo apt-add-repository "deb http://apt.kubernetes.io/ kubernetes-xenial main"
	
    # Install Kubeadm 
	sudo apt install kubeadm=1.15.3-00 kubectl=1.15.3-00 kubelet=1.15.3-00
	
    # Initialize Kubernetes on the master node
	sudo kubeadm init

    # start using your cluster 
	mkdir -p $HOME/.kube & sudo cp -i /etc/kubernetes/admin.conf $HOME/.kube/config & sudo chown $(id -u):$(id -g) $HOME/.kube/config

    # Deploy a Pod Network through the master node 
	kubectl apply -f "https://cloud.weave.works/k8s/net?k8s-version=$(kubectl version | base64 | tr -d '\n')"

    # Run 
	kubectl taint nodes --all node-role.kubernetes.io/master-
```

## MPI-operator (v 0.1.0)

* wget https://github.com/kubeflow/mpi-operator/archive/0.1.0.tar.gz
* tar -xvzf 0.1.0.tar.gz 
* cd mpi-operator-0.1.0/deploy/
* In 3-mpi-operator.yaml, change the two images (mpi-operator & kubectl-delivery) from latest to 0.1.0
```yaml
	- name: mpi-operator
        image: mpioperator/mpi-operator:0.1.0
        args: [
          "-alsologtostderr",
          "--gpus-per-node", "8",
          "--kubectl-delivery-image",
          "mpioperator/kubectl-delivery:0.1.0"
```
* Deploy mpi-operator:
```bash
   kubectl create -f 0-crd.yaml
   kubectl create -f 1-namespace.yaml
   kubectl create -f 2-rbac.yaml
   kubectl create -f 3-mpi-operator.yaml
```

## Rook Shared file system (release-0.8)

* git clone https://github.com/rook/rook.git 
* cd rook/
* git checkout release-0.8
* cd cluster/examples/kubernetes/ceph
```bash
   kubectl create -f operator.yaml
   kubectl create -f cluster.yaml
   kubectl create -f filesystem.yaml
   kubectl create -f storageclass.yaml
```
* NOTE: rook-ceph-block storageclass is not default. To set as default: 
```bash
   kubectl patch storageclass rook-ceph-block -p '{"metadata": {"annotations":{"storageclass.kubernetes.io/is-default-class":"true"}}}'
```

## Ingress

* kubectl apply -f  https://raw.githubusercontent.com/kubernetes/ingress-nginx/130af33510882ae62c89277f2ad0baca50e0fafe/deploy/static/mandatory.yaml
* Check ingress controller: kubectl get pods -n ingress-nginx
* mkdir ingress-deployment && cd ingress-deployment
* vi nginx-ingress.yaml
* Copy-Paste inside the file the below:
```yaml
kind: Service
apiVersion: v1
metadata:
  name: ingress-nginx
  namespace: ingress-nginx
  labels:
    app.kubernetes.io/name: ingress-nginx
    app.kubernetes.io/part-of: ingress-nginx
spec:
  externalTrafficPolicy: Local
  type: LoadBalancer
  selector:
    app.kubernetes.io/name: ingress-nginx
    app.kubernetes.io/part-of: ingress-nginx
  ports:
    - name: http
      port: 80
      targetPort: http
    - name: https
      port: 443
      targetPort: https
```

* kubectl apply -f nginx-ingress.yaml
* Check the created service 
```
kubectl get svc -n ingress-nginx
```

## Install Helm & Tiller
DARE uses the Helm package manager (https://helm.sh/) for Kubernetes to install and manage some of the external packages it uses. This facilitates the installation and upgrade of external components and prevents duplicated work on Kubernetes descriptors for well-known applications. To use Helm, we need to install the helm command and corresponding service (Tiller) first.
Summary:
* Download 3.1.1 release from https://github.com/helm/helm/releases
* Install helm binary to /usr/local/bin
* Initialize helm and Tiller
* Update package sources
```bash
## Download & unpack release package
$ wget https://get.helm.sh/helm-v3.1.1-linux-amd64.tar.gz
$ tar xf helm-v3.1.1-linux-amd64.tar.gz
## move to path
$ sudo mv linux-amd64/helm /usr/local/bin/
## create service account
$ kubectl create serviceaccount -n kube-system tiller
## create role for RBAC
$ kubectl create clusterrolebinding tiller-binding --clusterrole=cluster-admin --serviceaccount kube-system:tiller
## Update package sources
$ helm repo update
```

## Install cert-manager
To maintain certificates for the externally reachable services and pages, DARE uses the Let's Encrypt Certification Authority (https://letsencrypt.org/) through its ACME protocl (https://tools.ietf.org/html/rfc8555). The Cert-Manager addon (https://cert-manager.io) automates this process even further so that certificates are automatically issued, configured in the ingress and updated based on annotations in our Kubernetes descriptors. For installation, we use the official Helm package from Helm hub (https://hub.helm.sh/charts/jetstack/cert-manager). 
Summary:
* Install the custom resource definitions required for cert-manager
* Activate the Jetstack Helm repository
* Install the Cert-Manager package
* Add a clusterissue for letsencrypt

```bash
## install custom resource definitions
$ kubectl apply --validate=false -f https://raw.githubusercontent.com/jetstack/cert-manager/release-0.14/deploy/manifests/00-crds.yaml
## Add the Jetstack Helm repository
$ helm repo add jetstack https://charts.jetstack.io
## Install the cert-manager helm chart
$ helm install cert-manager --version v0.14.0 jetstack/cert-manager
## Install the letsencrypt ClusterIssuer (careful: needs customization!)
$ kubectl apply -f https://raw.githubusercontent.com/kubernetes/k8s.io/master/cert-manager/letsencrypt-prod.yaml
```

## Install Keycloak (WIP)

```bash
## Add codecentric package source and update package sources.
$ helm repo add codecentric https://codecentric.github.io/helm-charts
$ helm repo update
## Install Keycloak Helm Chart
$ helm install keycloak -f keycloak-values.yaml --version 8.0.0 codecentric/keycloak 

## To redeploy keycloak:
$ helm upgrade keycloak -f keycloak-values.yaml codecentric/keycloak
```
* Note: in order to find keycloak-values.yaml, go to k8s directories

Once Keycloak is deployed:
* Login to the Admin panel in keycloak UI and create a new realm named dare
* In the new realm, register the dare-login component as client. Use confidential strategy and after saving the client copy the secret key
* Paste the secret key in the dare-login-dp.yaml
* Create user accounts
* Create a user with username admin and some password. Copy the password in the dare-login-dp.yaml in readiness and liveness prob

## DARE platform
* git clone https://gitlab.com/project-dare/dare-platform.git
* cd dare-platform/k8s/
* execute the script deploy.sh
* Expose deployments:
```bash
    kubectl expose deployment d4p-registry --type=NodePort --name=d4p-registry-public
    kubectl expose deployment dare-login --type=NodePort --name=dare-login-public
    kubectl expose deployment exec-api --type=NodePort --name=exec-api-public
    kubectl expose deployment exec-registry --type=NodePort --name=exec-registry-public
    kubectl expose deployment playground --type=NodePort --name=playground-public
    kubectl expose deployment semantic-data --type=NodePort --name=semantic-data-public
    kubectl expose deployment workflow-registry --type=NodePort --name=workflow-registry-public
```
* If volumes are not mounted --> systemctl restart kubelet

## JupyterHub

This is an optional step. If you have multiple users and you want them to create an environment with existing notebooks,
you can create a JupyterHub. We have a demo yaml in k8s (jupyterhub-config.yaml) but since it contained sensitive
information we have comments to the fields that need to be completed.

* Create secret token for jupyterhub-config.yaml

```bash
    openssl rand -hex 32
```
* Copy the token and paste it in the aforementioned yaml file in the secretToken field
* Run the following

```bash
    helm repo add jupyterhub https://jupyterhub.github.io/helm-chart/
    helm repo update
    kubectl create namespace jupyterhub
    helm upgrade --install jupyterhub jupyterhub/jupyterhub --namespace jupyterhub --version=0.9.0 --values jupyterhub-config.yaml
```
---
title: "API Documentation"
---

In this page, information on DARE API is provided. DARE Platform follows the Microservices software architecture, therefore the DARE API is
constitute by multiple individual APIs. The following sections list the API endpoints of each DARE service. For high level description of 
the components visit the features section of our site.

The main component APIs described here are:

1. dare-login service which interacts with Keycloak in the backend
2. d4p-registry, API to interact with the dispel4py workflow registry
3. workflow-registry, API to interact with the CWL workflow registry
4. exec-api, the DARE execution API
5. s-prov, the provenance API
6. playground, which enables a development environment for scientists to write their workflows
7. semantic-data-discovery, API to retrieve data from the Data Catalogue

Finally, we provide documentation on the dispel4py library.

## DARE Login API

DARE login service interacts between DARE components and Keycloak. It exposes functionality for sign in to the platform, refreshing and validating a token etc. The main API calls are:

<table class="table">
	
<thead>
	<tr>
		<th scope="col"><b>HTTP method</b></th>
        <th scope="col"><b>Endpoint</b></th>
        <th scope="col"><b>Description</b></th>
        <th scope="col"><b>Content Type</b></th>
        <th scope="col"><b>Parameters</b></th>
	</tr>
</thead>
	<tbody>
		<tr>
			<td>POST</td>
			<td>/auth/</td>
			<td>Authenticates a user performing HTTP call <br>to the Keycloak service. <br>After having successfully authenticated <br>the user, Dispel4py Registry, <br>CWL registry and Execution API are <br>notified to check if the <br>user already exists in their local DBs</td>
			<td>application/json</td>
			<td>data (body), example: <br>{ <br>"username": "string", <br>"password": "string", <br>"requested_issuer": "string"}</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/validate-token/</td>
			<td>Validates a token using the Keycloak Service</td>
			<td>application/json</td>
			<td>data (body),example: <br>{<br>"access_token": "string"<br>}</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/delegation-token/</td>
			<td>Issues a token for internal application use</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>"access_token": "string"<br>}</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/refresh-token/</td>
			<td>Uses the refresh token to issue a new token for a user</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>"refresh_token": "string",<br>"issuer": "string"<br>}</td>
		</tr>
	</tbody>
</table>

The technical documentation of the dare-login component can be found [here](https://project-dare.gitlab.io/dare-login/).

## D4p Information Registry

The Dispel4py Workflow Registry enables the research developers to register their dispel4py workflows, re-use and share them. The fowlloing table shows the available API endpoints of this component.

<table class="table">
<thead>
	<tr>
            <th scope="col"><b>HTTP method</b></th>
            <th scope="col"><b>Endpoint</b></th>
            <th scope="col"><b>Description</b></th>
            <th scope="col"><b>Content Type</b></th>
            <th scope="col"><b>Parameters</b></th>
        </tr>
</thead>
<tbody>
	<tr>
		<td>GET</td>
		<td>/connections/</td>
		<td>Retrieves all the available PE Connection resources. <br>A PE Connection resource allows the <br>addition and manipulation <br>of PE connections. <br>Connections are associated with PEs and <br>are not themselves workspace items</td>
		<td>application/json</td>
		<td>No parameters</td>
	</tr>
	<tr>
		<td>POST</td>
		<td>/connections/</td>
		<td>Creates a new PE Connection resource, <br>which allows the addition <br>and manipulation of PE connections. <br>Connections are associated <br>with PEs and are not themselves <br>workspace items</td>
		<td>application/json</td>
		<td>data (body) example: {<br>  "comment" : "string" ,  <br>"kind" : "string" ,  <br>"modifiers" : "string" ,  <br>"name" : "string" ,  <br>"is_array" : true ,  <br>"s_type" : "string" ,  <br>"d_type" : "string" ,  <br>"pesig" : "string" }</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/connections/{id}/</td>
		<td>Retrieves a specific PE Connection resource. <br>A PE Connection resource allows the addition <br>and manipulation of PE connections. <br>Connections are associated with PEs and <br>are not themselves workspace items.</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>PUT</td>
		<td>/connections/{id}/</td>
		<td>Updates an existing PE Connetion resource. <Br>A PE Connection resource allows the addition <br>and manipulation of PE connections. <br>Connections are associated with PEs and <br>are not themselves workspace items.</td>
		<td>application/json</td>
		<td>-id (integer) <br>-data (body) example: <br>{<br>  "comment" : "string" ,  <br>"kind" : "string" ,  <br>"modifiers" : "string" ,  <br>"name" : "string" ,  <br>"is_array" : true ,  <br>"s_type" : "string" ,  <br>"d_type" : "string" ,  <br>"pesig" : "string" <br>}</td>
	</tr>
	<tr>
		<td>DELETE</td>
		<td>/connections/{id}/</td>
		<td>Deletes an existing PE Connection resource <br>from the DB</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/fnimpls/</td>
		<td>Retrieve all the available function <br>implementation resources</td>
		<td>application/json</td>
		<td>No parameters</td>
	</tr>
	<tr>
		<td>POST</td>
		<td>/fnimpls/</td>
		<td>Creates a new Function Implementation</td>
		<td>application/json</td>
		<td>data (body), example: <br>{<br>"code" : "string", <br>"parent_sig": "string", <br>"description" : "string", <br>"pckg" : "string", <br>"workspace": "string" , <br>"clone_of": "string" , <br>"name": "string" <br>}</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/fnimpls/{id}/</td>
		<td>Retrieves a specific Function implementation resource</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>PUT</td>
		<td>/fnimpls/{id}/</td>
		<td>Updates an existing function implementation</td>
		<td>application/json</td>
		<td>-id (integer) <br>-data (body), example: <br>{ <br>"code": "string", <br>"parent_sig": "string", <br>"description": "string", <br>"pckg": "string", <br>"workspace": "string", <br>"clone_of": "string", <br>"name": "string" <br>}</td>
	</tr>
	<tr>
		<td>DELETE</td>
		<td>/fnimpls/{id}/</td>
		<td>Deletes an existing Function Implementation</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/fnparams/</td>
		<td>Retrieves all the available Function Parameters</td>
		<td>application/json</td>
		<td>No parameters</td>
	</tr>
	<tr>
		<td>POST</td>
		<td>/fnparams/</td>
		<td>Creates a new Function Parameters</td>
		<td>application/json</td>
		<td>data (body), example: <br>{<br>"parent_function": "string", <br>"param_name": "string", <br>"param_type" : "string" <br>}</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/fnparams/{id}/</td>
		<td>Retrieves a specific Function Parameters</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>PUT</td>
		<td>/fnparams/{id}/</td>
		<td>Updates an existing Function <br>Parameters entry</td>
		<td>application/json</td>
		<td>-id (integer) <br>-data (body) example: <br>{<br>"parent_function": "string", <br>"param_name": "string", <br>"param_type": "string" <br>}</td>
	</tr>
	<tr>
		<td>DELETE</td>
		<td>/fnparams/{id}/</td>
		<td>Deletes an existing Function Parameters</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/functions/</td>
		<td>Retrieves all the Function resources <br>from the DB</td>
		<td>application/json</td>
		<td>No parameters</td>
	</tr>
	<tr>
		<td>POST</td>
		<td>/functions/</td>
		<td>Creates a new Function Resource</td>
		<td>application/json</td>
		<td>data (body), example: <br>{ <br>"description" : "string", <br>"parameters" : ["string"], <br>"fnimpls": ["string"], <br>"pckg": "string", <br>"workspace": "string", <br>"return_type": "string", <br>"clone_of": "string", <br>"name": "string" <br>}</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/functions/{id}/</td>
		<td>Retrieves an existing Function resource</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>PUT</td>
		<td>/functions/{id}/</td>
		<td>Updates an existing function resouce</td>
		<td>application/json</td>
		<td>-id (integer) <br>-data(body), example: <br>{ <br>"description": "string", <br>"parameters": ["string"], <br>"fnimpls": ["string"], <br>"pckg": "string", <br>"workspace": "string", <br>"return_type": "string", <br>"clone_of": "string", <br>"name": "string" <br>}</td>
	</tr>
	<tr>
		<td>DELETE</td>
		<td>/functions/{id}/</td>
		<td>Deletes an existing function resource</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/groups/</td>
		<td>Retrieves all the available groups</td>
		<td>application/json</td>
		<td>No parameters</td>
	</tr>
	<tr>
		<td>POST</td>
		<td>/groups/</td>
		<td>Creates a new user group</td>
		<td>application/json</td>
		<td>data (body), example: { "name": "string" }</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/groups/{id}/</td>
		<td>Retrieves a user group</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>PUT</td>
		<td>/groups/{id}/</td>
		<td>Updates an existing user group</td>
		<td>application/json</td>
		<td>-id (integer) <br>-data (body), <br>example: { "name": "string" }</td>
	</tr>
	<tr>
		<td>DELETE</td>
		<td>/groups/{id}/</td>
		<td>Removes a user group</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/literals/</td>
		<td>Retrieves all the literal entities</td>
		<td>application/json</td>
		<td>No parameters</td>
	</tr>
	<tr>
		<td>POST</td>
		<td>/literals/</td>
		<td>Creates a new Literal Entities</td>
		<td>application/json</td>
		<td>data (body), example: <br>{ <br>"description": "string", <br>"value": "string", <br>"name": "string", <br>"pckg": "string", <br>"workspace": "string", <br>"clone_of" : "string" <br>}</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/literals/{id}/</td>
		<td>Retrieves a Literal Entities</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>PUT</td>
		<td>/literals/{id}/</td>
		<td>Updates an existing Literal Entities</td>
		<td>application/json</td>
		<td>-id (integer) <br>-data (body), example: <br>{ <br>"description": "string", <br>"value": "string", <br>"name": "string", <br>"pckg": "string", <br>"workspace": "string", <br>"clone_of": "string" <br>}</td>
	</tr>
	<tr>
		<td>DELETE</td>
		<td>/literals/{id}/</td>
		<td>Deletes a Literal Entities</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/peimpls/</td>
		<td>Retrieves all the available <br>PE Implementation</td>
		<td>application/json</td>
		<td>No parameters</td>
	</tr>
	<tr>
		<td>POST</td>
		<td>/peimpls/</td>
		<td>Creates a new PE Implementation</td>
		<td>application/json</td>
		<td>data (body), example: <br>{ <br>"code": "string", <br>"parent_sig": "string", <br>"description": "string", <br>"pckg": "string", <br>"workspace": "string", <br>"clone_of": "string", <br>"name": "string" <br>}</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/peimpls/{id}/</td>
		<td>Retrieves a specific <br>PE Implementation</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>PUT</td>
		<td>/peimpls/{id}/</td>
		<td>Updates an existing PE Implementation</td>
		<td>application/json</td>
		<td>-id (integer) <br>-data(body), example: <br>{ <br>"code": "string", <br>"parent_sig": "string", <br>"description": "string", <br>"pckg": "string", "workspace": "string", <br>"clone_of" : "string", <br>"name": "string" <br>}</td>
	</tr>
	<tr>
		<td>DELETE</td>
		<td>/peimpls/{id}/</td>
		<td>Deletes an existing PE Implementation</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/pes/</td>
		<td>Retrieves all the available PE resources</td>
		<td>application/json</td>
		<td>No parameters</td>
	</tr>
	<tr>
		<td>POST</td>
		<td>/pes/</td>
		<td>Creates a new PE</td>
		<td>application/json</td>
		<td>data (body), example: <br>{<br>"description": "string", <br>"name": "string", <br>"connections": ["string"], <br>"pckg": "string", <br>"workspace": "string", <br>"clone_of": "string", <br>"peimpls": ["string"] <br>}</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/pes/{id}/</td>
		<td>Retrieves a specific PE</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>PUT</td>
		<td>/pes/{id}/</td>
		<td>Updates an existing PE</td>
		<td>application/json</td>
		<td>-id(integer) <br>-data(body), example: <br>{ <br>"description": "string", <br>"name": "string", <br>"connections": ["string"], <br>"pckg": "string", <br>"workspace": "string", <br>"clone_of": "string", <br>"peimpls": ["string"] <br>}</td>
	</tr>
	<tr>
		<td>DELETE</td>
		<td>/pes/{id}/</td>
		<td>Deletes an existing PE</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/registryusergroups/</td>
		<td>Retrieves all the available <br>registry user groups</td>
		<td>application/json</td>
		<td>No parameters</td>
	</tr>
	<tr>
		<td>POST</td>
		<td>/registryusergroups/</td>
		<td>Creates a new Registry user group</td>
		<td>application/json</td>
		<td>data (body), example: <br>{ <br>"description": "string", <br>"group_name": "string" <br>}</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/registryusergroups/{id}/</td>
		<td>Retrieves a specific Registry user group</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>PUT</td>
		<td>/registryusergroups/{id}/</td>
		<td>Updates an existing Registry user group</td>
		<td>application/json</td>
		<td>-id(integer) <br>-data(body), example: <br>{ <br>"owner": "string", <br>"description": "string", <br>"group_name": "string" <br>}</td>
	</tr>
	<tr>
		<td>DELETE</td>
		<td>/registryusergroups/{id}/</td>
		<td>Deletes an existing Registry user group</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/users/</td>
		<td>Retrieves all the existing users</td>
		<td>application/json</td>
		<td>No parameters</td>
	</tr>
	<tr>
		<td>POST</td>
		<td>/users/</td>
		<td>Creates a new user</td>
		<td>application/json</td>
		<td>data (body), example: <br>{ <br>"username": "string", <br>"password": "string", <br>"first_name": "string", <br>"last_name": "string", <br>"email": "string" <br>}</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/users/{id}/</td>
		<td>Retrieves a specific user</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>PUT</td>
		<td>/users/{id}/</td>
		<td>Updates a specific user</td>
		<td>application/json</td>
		<td>-id(integer) <br>-data(body), example: <br>{ <br>"username": "string", <br>"password": "string", <br>"first_name": "string", <br>"last_name": "string", <br>"email": "string" <br>}</td>
	</tr>
	<tr>
		<td>DELETE</td>
		<td>/users/{id}/</td>
		<td>Deletes a specific user</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>/workspaces/</td>
		<td>Retrieves all the available workspaces</td>
		<td>application/json</td>
		<td>parameters:<br> <br>-name: name <br>-description: The name of the workspace we want to display <br>-paramType: query <br><br>-name: username <br>-description: The username the workspace is associated with  <br>-paramType: query <br><br>-name: search <br>- description: perform a simple full-text on <br>descriptions and names of workspaces. <br>-paramType: query</td>
	</tr>
	<tr>
		<td>POST</td>
		<td>/workspaces/</td>
		<td>Create or clone a new workspace</td>
		<td>application/json</td>
		<td>parameters: <br> <br>- name: name <br>- description: the name of the workspace. <br><br>- name: description <br>- description: a textual description of the workspace. <br><br>- name: clone_of <br>- description: indicates that a cloning operation is requested. <br>- paramType: query <br>- type: long</td>
	</tr>
	<tr>
		<td>GET</td>
		<td>Retrieves a specific workspace</td>
		<td></td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
	<tr>
		<td>PUT</td>
		<td>/workspaces/{id}/</td>
		<td>Updates an existing workspace</td>
		<td>application/json</td>
		<td>-id (integer) <br>-data (body), example: <br>{ <br>"clone_of": "string", <br>"name": "string", <br>"description": "string" <br>}</td>
	</tr>
	<tr>
		<td>DELETE</td>
		<td>/workspaces/{id}/</td>
		<td>Deletes an existing workspace</td>
		<td>application/json</td>
		<td>id (integer)</td>
	</tr>
</tbody>
</table>

## CWL Workflow Registry

This component is a Django Web Service exposing an API for CWLs and dockers registration. The technical documentation for the CWL workflow registry is available in the project's micro site [here](https://project-dare.gitlab.io/workflow-registry/). The following table shows all the available API calls in the CWL workflow registry.

<table>
	<thead>
		<tr>
	        <th scope="col"><b>HTTP method</b></th>
	        <th scope="col"><b>Endpoint</b></th>
	        <th scope="col"><b>Description</b></th>
	        <th scope="col"><b>Content Type</b></th>
	        <th scope="col"><b>Parameters</b></th>
	    </tr>
	</thead>
	<tbody>
		<tr>
			<td>POST</td>
			<td>/docker/</td>
			<td>Creates a new Docker Environment. <br>The environment consist of a Dockerfile and <br>can be associated with one or <br>multiple DockerScript entries <br>(which represent bash or python scripts)</td>
			<td>application/json</td>
			<td>data (body), example:<br>{<br>
        	"docker_name": "name",<br>
       		"docker_tag": "tag",<br>
      		"script_names": ["script1.sh", "script2.sh"]<br>
        	"files": <br>{<br>
            	"dockerfile": "string",<br>
            	"script1.sh": "string.",<br>
            	"script2.sh": "string"<br>
        	},<br>
        	"access_token": "token"<br>}
    		</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/docker/update_docker/</td>
			<td>Updates an existing Docker Environment</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>
				"docker_name": "name",<br>
				"docker_tag" "tag",<br>
				"update": {"tag": "v2.0"},<br>
				"files": {"dockerfile": "string"},<br>
				"access_token": "token"<br>}
			</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/docker/provide_url/</td>
			<td>Updates an existing Docker environment’s url field. <br> Once the docker image is built <br>and pushed in a public repository, <br>the relevant Docker entry should <br>be updated with the URL.</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>
				"docker_name": "name",<br>
				"docker_tag": "tag",<br>
				"docker_url": "url",<br>
				"access_token": "token"<br>}
			</td>
		</tr>
		<tr>
			<td>DELETE</td>
			<td>/docker/delete_docker/</td>
			<td>Deletes an existing docker environment.</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>
				"docker_name": "name",<br>
				"docker_tag": "tag",<br>
				"access_token": "token"<br>}
			</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/docker/bynametag/</td>
			<td>Retrieves a Docker Environment using <br>its name and tag.</td>
			<td>application/json</td>
			<td>-docker_name (string)<br>-docker_tag(string)</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/docker/byuser/</td>
			<td>Retrieves all the registered Docker <br>environments by user</td>
			<td>application/json</td>
			<td>-requested_user(string) if exists, <br>otherwise it uses the user that <br>performed the request</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/docker/download/</td>
			<td>Downloads in a zip file the Dockerfile <br>and the relevant scripts of a <br>Docker Environment.</td>
			<td>application/json</td>
			<td>-docker_name (string)<br>-docker_tag(string)</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/scripts/add/</td>
			<td>Adds a new script in an existing <br>Docker Environment</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>
				"docker_name": "name",<br>
				"docker_tag": "tag",<br>
				"script_name": "entrypoint.sh",<br>
				"files": {"entrypoint.sh": "string"},<br>
				"access_token": "token"<br>}
			</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/scripts/edit/</td>
			<td>Edits an existing script of a <br>Docker Environment</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>
				"docker_name": "name",<br>
				"docker_tag": "tag",<br>
				"script_name": "entrypoint.sh",<br>
				"files": {"entrypoint.sh": “string”},<br>
				"access_token": "token"<br>}
			</td>
		</tr>
		<tr>
			<td>DELETE</td>
			<td>/scripts/delete/</td>
			<td>Deletes an existing script from a docker environment</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>
				"docker_name": "name",<br>
				"docker_tag": "tag",<br>
				"script_name": "entrypoint.sh",<br>
				"access_token": "token"<br>}
			</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/scripts/download</td>
			<td>Downloads a specific script from <br>a Docker Environment</td>
			<td>application/json</td>
			<td>-docker_name(string)<br>-docker_tag(string)<br>-script_name(string)</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/scripts/byname</td>
			<td>Retrieves a specific script based on the name <br>& tag of the Docker Environment <br>and on the name of the script.</td>
			<td>application/json</td>
			<td>-docker_name(string)<br>-docker_tag(string)<br>-script_name(string)</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/workflows/</td>
			<td>Creates a new CWL workflow of <br>class Workflow</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>
				"workflow_name": "demo_workflow.cwl",<br>
				"workflow_version": "v1.0",<br>
				"spec_file_name": "spec.yaml",<br>
				"docker_name": "name",<br>
				"docker_tag": "tag",<br>
				"workflow_part_data":<br>[{<br>"name":arguments.cwl”,<br>
				"version":"v1.0", <br>"spec_name": "arguments.yaml"<br>},<br>{<br>"name": "tar_param.cwl", <br>"version":"v1.0", <br>"spec_name": "tar_param.yaml"<br>}],<br>
				"files": <br>{<br>"demo_workflow.cwl":"string",<br>
				"spec.yaml": "string",<br>
				"arguments.cwl": "string",<br>
				"arguments.yaml": "string",<br>
				"tar_param.cwl": "string",<br>
				"tar_param.yaml": "string"<br>},<br>
				"access_token": "token"<br>
				}
			</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/workflows/update_workflow/</td>
			<td>Updates an existing CWL workflow of <br>class Workflow</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>
				"workflow_name":"demo_workflow.cwl",<br>
				"workflow_version": "v1.0",<br>
				"files": {"workflow_file": "string",<br>
				"spec_file": "string",},<br>
				"update": {"version":"v1.1"},<br>
				"access_token": "token"<br>}
			</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/workflows/update_docker/</td>
			<td>Associate a CWL workflow of class <br>Workflow with a different Docker <br>Environment.</td>
			<td>application/json</td>
			<td>data (body), example:<br>
				{<br>
				"workflow_name":"demo_workflow.cwl",<br>
				"workflow_version": "v1.0",<br>
				"docker_name": "test",<br>
				"docker_tag": "v1.0",<br>
				"access_token": "token"<br>
				}
			</td>
		</tr>
		<tr>
			<td>DELETE</td>
			<td>/workflows/delete_workflow/</td>
			<td>Deletes an existing CWL workflow <br>(class Workflow) and all the <br>associated Workflow parts <br>(class CommandLineTool).</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>
				"workflow_name":"demo_workflow.cwl",<br>
				"workflow_version": "v1.0",<br>
				"access_token": "token"<br>
				}
			</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/workflows/bynameversion/</td>
			<td>Retrieve a CWL workflow of class Workflow <br>and its associated workflow parts <br>as well as the related docker <br>environment, based on the workflow <br>name and version.</td>
			<td>application/json</td>
			<td>-workflow_name(string)<br>-workflow_version(string)</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/workflows/download</td>
			<td>Downloads in a zip file all the CWL files <br>(Workflow and CommandLineTool) as well as <br>the relevant Dockerfile and scripts <br>(if the parameter dockerized is provided)</td>
			<td>application/json</td>
			<td>-workflow_name(string)<br>-workflow_version(string)<br>-dockerized(boolean)</td>
		</tr>
		<tr>
			<td></td>
			<td></td>
			<td></td>
			<td>application/json</td>
			<td></td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/workflow_parts/add/</td>
			<td>Adds a new CommandLineTool CWL in an existing CWL workflow</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>
				"workflow_name":"demo_workflow.cwl",<br>
				"workflow_version": "v1.0",<br>
				"workflow_part_name":"arguments.cwl",<br>
				"workflow_part_version": "v1.0",<br>
				"spec_name": "arguments.yaml",<br>
				"files": {"arguments.cwl": "string",<br>
				"arguments.yaml": "string"},<br>
				"access_token": "token"<br>
				}
			</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/workflow_parts/edit/</td>
			<td>Edits an existing CommandLineTool CWL workflow</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>
				"workflow_name":"demo_workflow.cwl",<br>
				"workflow_version": "v1.0",<br>
				"workflow_part_name":"arguments.cwl",<br>
				"workflow_part_version": "v1.0",<br>
				"spec_name": "arguments.yaml",<br>
				"files": <br>{<br>"arguments.cwl":"string",<br>
				"arguments.yaml": "string"<br>},<br>
				"update": <br>{<br>"version":"v1.1”<br>},<br>
				"access_token": "token"<br>
				}
			</td>
		</tr>
		<tr>
			<td>DELETE</td>
			<td>/workflow_parts/delete/</td>
			<td>Deletes an existing CommandLIneTool CWL workflow</td>
			<td>application/json</td>
			<td>data (body), example:
				<br>{<br>
				"workflow_name":"demo_workflow.cwl",<br>
				"workflow_version": "v1.0",<br>
				"workflow_part_name": "arguments.cwl",<br>
				"workflow_part_version": "v1.0",<br>
				"access_token": "token"<br>
				}
			</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/workflow_parts/bynameversion</td>
			<td>Retrieves a specific CommandLineTool CWL based <br>on its parent (name & version) and <br>its own name and version.</td>
			<td>application/json</td>
			<td>-workflow_name(string)<br>-workflow_version(string)<br>-workflow_part_name(string)<br>-workflow_part_version(string)</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/workflow_parts/download/</td>
			<td>Downloads a specific CWL of class CommandLineTool</td>
			<td>application/json</td>
			<td>-workflow_name(string)<br>-workflow_version(string)<br>-workflow_part_name(string)<br>-workflow_part_version(string)</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/accounts/login/</td>
			<td>Authenticates a user (login) used by the <br>dare-login component described above <br>when a user calls the /auth/ endpoint of the dare-login. <br>If the user does not exist in the CWL workflow <br>registry’s local DB, <br>it creates a new user.</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>
				"username": "string",<br>
				"password":"string",<br>
				"access_token":"string",<br>
				"email":"string",<br>
				"given_name":"string",<br>
				"family_name":"string"<br>
				}
			</td>
		</tr>
	</tbody>
</table>

## Execution API

### General

Execution API provides endpoints for multiple execution contexts:

* Dispel4py: dynamically creates containers to execute Dispel4py workflows
* CWL: execution environments spawned on-demand to execute CWL workflows.
* Specfem: creates containers in a dynamic way to execute Specfem executable. This endpoint is deprecated, Specfem is now executed via the CWL endpoint.

### API calls

<table>
	<thead>
		<tr>
	        <th scope="col"><b>HTTP method</b></th>
	        <th scope="col"><b>Endpoint</b></th>
	        <th scope="col"><b>Description</b></th>
	        <th scope="col"><b>Content Type</b></th>
	        <th scope="col"><b>Parameters</b></th>
	    </tr>
	</thead>
	<tbody>
		<tr>
			<td>POST</td>
			<td>/create-folders/</td>
			<td>Endpoint used by the /auth/ endpoint of dare-login. <br>Checks if the user’s workspace in the DARE platform <br>is available, otherwise it creates the <br>necessary folder structure</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>	"username": "string"<br>}</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/d4p-mpi-spec/</td>
			<td>Used internally by the dispel4py execution environment <br>in order to retrieve the respective <br>PE Implementation and spec.yaml</td>
			<td>application/json</td>
			<td>data (body),example: <br>{<br>
				"pe_imple": "name",<br>
				"nodes": 3,<br>
				"input_data": {}<br>
				}
			</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/run-d4p/</td>
			<td>Creates a new dispel4py execution environment, <br>using the Kubernetes API. <br>Generates a new run directory, <br>stored under the user’s “runs” folder<br>(i.e. /<home>/<username>/runs/). <br>All the execution results are <br>stored in the generated run directory.</td>
			<td>application/json</td>
			<td>data (body), example: <br>{<br>
				"access_token": "string",<br>
				"workspace": "string",<br>
				"pckg": "string",<br>
				"pe_name":"string",<br>
				"target":"string",<br>
				"nodes":1<br>
				}
			</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/run-specfem/</td>
			<td>Deprecated endpoint. Use /run-cwl instead.</td>
			<td>application/json</td>
			<td>-</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/run-cwl/</td>
			<td>Endpoint to instantiate an execution environment <br>for CWL workflow execution. <br>The environment to be instantiated <br>is retrieved from the CWL using the <br>CWL Workflow Registry. Generates a new run <br>directory, stored under the user’s “runs” folder<br>(i.e. /<home>/<username>/runs/). All the execution <br>results are stored in the <br>generated run directory.</td>
			<td>application/json</td>
			<td>data (body), example:
				<br>{<br>
				"access_token": "string",<br>
				"nodes":12,<br>
				"workflow_name":"string",<br>
				"workflow_version": "string",<br>
				"input_data": <br>{<br>
				"example1":"string"<br>
				}<br>
				}
			</td>
		</tr>
		<tr>
			<td>POST</td>
			<td>/upload/</td>
			<td>Endpoint used to upload files in the DARE platform. <br>The files are stored under the user’s home directory. <br>The home directory is named after his/hers username and <br>inside there are 3 folders, i.e. uploads, debug and runs.<br>All the uploaded files are stored under the user’s “uploads” directory</td>
			<td>application/json</td>
			<td>data (body), example:
				<br>{<br>
				"dataset_name": "string",<br>
				"path": "string",<br>
				"access_token": "string",<br>
				"files": [<file list>]<br>
				}
			</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/my-files/</td>
			<td>Lists all the users’ directories under the “uploads”, <br>“runs” and “debug” folders. If the parameter <br>num_run_dirs is present, the response <br>is limited to the most recent directories based on <br>the number provided in the aforementioned parameter</td>
			<td>application/json</td>
			<td>-access_token(string)<br>-num_run_dirs(integer)</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/list/</td>
			<td>Lists all the files inside a specific directory. <br>This directory could be retrieved from the previous endpoint</td>
			<td>application/json</td>
			<td>-access_token(string)<br>-path(string)</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/download/</td>
			<td>Downloads a specific file from the DARE platform. <br>To find the file’s full path use the two previous endpoints</td>
			<td>application/json</td>
			<td>-access_token(string)<br>-path(string)</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/send2drop/</td>
			<td>Uploads files from the dare platform to B2DROP</td>
			<td>application/json</td>
			<td>-access_token(string)<br>-path(string)</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/cleanup/</td>
			<td>Clears the user’s folders (uploads, runs, debug)</td>
			<td>application/json</td>
			<td>-access_token(string)<br>-runs(boolean)<br>-uploads(boolean)<br>-debug(boolean)</td>
		</tr>
		<tr>
			<td>GET</td>
			<td>/my-pods</td>
			<td>List the running jobs of a user</td>
			<td>application/json</td>
			<td>data example: {"access_token": "string"}</td>
		</tr>
	</tbody>
</table>

Technical documentation of the component is also available [here](https://project-dare.gitlab.io/exec-api/)

## Provenance

**Version:** v1

### /data
---
##### ***GET***
**Description:** The data is selected by specifying a query string. Query parameters allow to search by attribution to a component or to an implementation, generation by a workflow execution and by combining more metadata and parameters terms with their min and max valuesranges. Mode of the search can also be indicated (mode ::= (OR | AND). It will apply to the search upon metadata and parameters values-ranges

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| usernames | query | csv list of users the Workflows Executons are associated with | No | string |
| terms | query | csv list of metadata or parameter terms. These relate positionally to the maxvalues and the minvalues | No | string |
| functionNames | query | csv list of functions the Data was generated with | No | string |
| wasAttributedTo | query | csv list of Component or Component Instances involved in the generation of the Data | No | string |
| minvalues | query | csv list of metadata or parameters minvalues. These relate positionally to the terms and the minvalues | No | string |
| rformat | query | unimplemented: format of the response payload (json,json-ld) | No | string |
| start | query | index of the starting item | Yes | integer |
| limit | query | max number of items expected | Yes | integer |
| maxvalues | query | csv list of metadata or parameters maxvalues. These relate positionally to the terms and the minvalues | No | string |
| formats | query | csv list of data formats (eg. mime-types) | No | string |
| types | query | csv list of data types | No | string |
| wasGeneratedBy | query | the id of the Invocation that generated the Data | No | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /data/filterOnAncestor
---
##### ***POST***
**Description:** Filter a list of data ids based on the existence of at least one ancestor in their data dependency graph, according to a list of metadata terms and their min and max values-ranges. Maximum depth level and mode of the search can also be indicated (mode ::= (OR | AND)

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| body | body |  | No | object |

**Responses**

| Code | Description |
| ---- | ----------- |

### /data/{data_id}
---
##### ***GET***
**Description:** Extract Data and their DataGranules by the Data id

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| data_id | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /data/{data_id}/derivedData
---
##### ***GET***
**Description:** Starting from a specific data entity of the data dependency is possible to navigate through the derived data or backwards across the element's data dependencies. The number of traversal steps is provided as a parameter (level).

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| level | query | level of depth in the data derivation graph, starting from the current Data | Yes | string |
| data_id | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /data/{data_id}/export
---
##### ***GET***
**Description:** Export of provenance information PROV-XML or RDF format. The S-PROV information returned covers the whole workflow execution or is restricted to a single data element. In the latter case, the graph is returned by following the derivations within and across runs. A level parameter allows to indicate the depth of the resulting trace

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| format | query | export format of the PROV document returned | No | string |
| rdfout | query | export rdf format of the PROV document returned | No | string |
| creator | query | the name of the user requesting the export | No | string |
| level | query | level of depth in the data derivation graph, starting from the current Data | Yes | string |
| data_id | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /data/{data_id}/wasDerivedFrom
---
##### ***GET***
**Description:** Starting from a specific data entity of the data dependency is possible to navigate through the derived data or backwards across the element's data dependencies. The number of traversal steps is provided as a parameter (level).

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| level | query | level of depth in the data derivation graph, starting from the current Data | Yes | string |
| data_id | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /instances/{instid}
---
##### ***GET***
**Description:** Extract details about a single instance or component by specifying its id. The returning document will indicate the changes that occurred, reporting the first invocation affected. It support the specification of a list of runIds the instance was wasAssociateFor, considering that the same instance could be used across multiple runs

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| start | query | index of the starting item | Yes | integer |
| limit | query | max number of items expected | Yes | integer |
| wasAssociateFor | query | cvs list of runIds the instance was wasAssociateFor (when more instances are reused in multiple workflow executions) | No | string |
| instid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /invocations/{invocid}
---
##### ***GET***
**Description:** Extract details about a single invocation by specifying its id

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| invocid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /summaries/collaborative
---
##### ***GET***
**Description:** Extract information about the reuse and exchange of data between workflow executions based on terms' valuesranges and a group of users. The API method allows for inclusive or exclusive (mode ::= (OR j AND) queries on the terms' values. As above, additional details, such as running infrastructure, type and name of the workflow can be selectively extracted by assigning these properties to a groupBy parameter. This will support the generation of grouped views

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| wasAssociatedWith | query | csv lis of Components involved in the Workflow Executions | No | string |
| usernames | query | csv list of users the Workflows Executons are associated with | No | string |
| terms | query | csv list of metadata or parameter terms. These relate positionally to the maxvalues and the minvalues | No | string |
| functionNames | query | csv list of functions that are executed by at least one workflow's components | No | string |
| level | query | level of depth in the data derivation graph, starting from the current Data | Yes | string |
| minvalues | query | csv list of metadata or parameters minvalues. These relate positionally to the terms and the minvalues | No | string |
| rformat | query | unimplemented: format of the response payload (json,json-ld) | No | string |
| maxvalues | query | csv list of metadata or parameters maxvalues. These relate positionally to the terms and the minvalues | No | string |
| formats | query | csv list of data formats (eg. mime-types) | No | string |
| clusters | query | csv list of clusters that describe and group one or more workflow's component | No | string |
| groupby | query | express the grouping of the returned data | No | string |
| types | query |  | No | string |
| mode | query | execution mode of the workflow in case it support different kind of concrete mappings (eg. mpi, simple, multiprocess, etc.. | No | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /summaries/workflowexecution
---
##### ***GET***
**Description:** Produce a detailed overview of the distribution of the computation, reporting the size of data movements between the workflow components, their instances or invocations across worker nodes, depending on the specified granularity level. Additional information, such as process pid, worker, instance or component of the workflow (depending on the level of granularity) can be selectively extracted by assigning these properties to a groupBy parameter. This will support the generation of grouped views

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| level | query | level of depth in the data derivation graph, starting from the current Data | Yes | string |
| mintime | query | minimum start time of the Invocation | No | string |
| groupby | query | express the grouping of the returned data | No | string |
| runId | query | the id of the run to be analysed | No | string |
| maxtime | query | maximum start time of the Invocation | No | string |
| maxidx | query | maximum iteration index of an Invocation | No | integer |
| minidx | query | minimum iteration index of an Invocation | No | integer |

**Responses**

| Code | Description |
| ---- | ----------- |

### /terms
---
##### ***GET***
**Description:** Return a list of discoverable metadata terms based on their appearance for a list of runIds, usernames, or for the whole provenance archive. Terms are returned indicating their type (when consistently used), min and max values and their number occurrences within the scope of the search

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| aggregationLevel | query | set whether the terms need to be aggreagated by runId, username or across the whole collection (all) | No | string |
| usernames | query | csv list of usernames | No | string |
| runIds | query | csv list of run ids | No | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions
---
##### ***GET***
**Description:** Extract documents from the bundle collection according to a query string which may include usernames, type of the workflow, the components the run wasAssociatedWith and their implementations. Data results' metadata and parameters can also be queried by specifying the terms and their min and max values-ranges and data formats. Mode of the search can also be indicated (mode ::= (OR j AND). It will apply to the search upon metadata and parameters values of each run

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| wasAssociatedWith | query | csv lis of Components involved in the Workflow Executions | No | string |
| usernames | query | csv list of users the Workflows Executons are associated with | No | string |
| terms | query | csv list of metadata or parameter terms. These relate positionally to the maxvalues and the minvalues | No | string |
| functionNames | query | csv list of functions that are executed by at least one workflow's components | No | string |
| minvalues | query | csv list of metadata or parameters minvalues. These relate positionally to the terms and the minvalues | No | string |
| rformat | query | unimplemented: format of the response payload (json,json-ld) | No | string |
| start | query | index of the starting item | Yes | integer |
| limit | query | max number of items expected | Yes | integer |
| maxvalues | query | csv list of metadata or parameters maxvalues. These relate positionally to the terms and the minvalues | No | string |
| formats | query | csv list of data formats (eg. mime-types) | No | string |
| clusters | query | csv list of clusters that describe and group one or more workflow's component | No | string |
| types | query |  | No | string |
| mode | query | execution mode of the workflow in case it support different kind of concrete mappings (eg. mpi, simple, multiprocess, etc.. | No | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions/insert
---
##### ***POST***
**Description:** Bulk insert of bundle or lineage documents in JSON format. These must be provided as encoded stirng in a POST request

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| body | body |  | No | object |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions/import
---
##### ***POST***
**Description:** Import of provenance output which is not yet mapped to the s-ProvFlowMongoDB format. 
The files provided in the archive will be mapped to s-ProvFlowMongoDB if they are in one of the supported formats.

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| archive | form | Zip archive of provenance output, which will be mapped to s-ProvFlowMongoDB and stored. Currently only files in the CWLProv format are supported | Yes | file |
| format | form | Format of the provenance output to be imported. | Yes | String |


**Responses**

| Code | Description |
| ---- | ----------- |



### /workflowexecutions/{run_id}/export
---
##### ***GET***
**Description:** Export of provenance information PROV-XML or RDF format. The S-PROV information returned covers the whole workflow execution or is restricted to a single data element. In the latter case, the graph is returned by following the derivations within and across runs. A level parameter allows to indicate the depth of the resulting trace

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| rdfout | query | export rdf format of the PROV document returned | No | string |
| creator | query | the name of the user requesting the export | No | string |
| format | query | export format of the PROV document returned | No | string |
| run_id | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions/{runid}
---
##### ***DELETE***
**Description:** Extract documents from the bundle collection by the runid of a WFExecution. The method will return input data and infomation about the components and the libraries used for the specific run

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| runid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

##### ***GET***
**Description:** Extract documents from the bundle collection by the runid of a WFExecution. The method will return input data and infomation about the components and the libraries used for the specific run

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| runid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions/{runid}/delete
---
##### ***POST***
**Description:** Delete a workflow execution trace, including its bundle and all its lineage documents

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| runid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions/{runid}/edit
---
##### ***POST***
**Description:** Update of the description of a workflow execution. Users can improve this information in free-tex

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| body | body |  | No | object |
| runid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions/{runid}/showactivity
---
##### ***GET***
**Description:** Extract detailed information related to the activity related to a WFExecution (id). The result-set can be grouped by invocations, instances or components (parameter level) and shows progress, anomalies (such as exceptions or systems' and users messages), occurrence of changes and the rapid availability of accessible data bearing intermediate results. This method can also be used for runtime monitoring

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| start | query | index of the starting item | Yes | integer |
| limit | query | max number of items expected | Yes | integer |
| level | query | level of aggregation of the monitoring information (component, instance, invocation, cluster) | No | string |
| runid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |


## Testing environment

The purpose of this component is to provide a DARE environment for test and debugging purposes. The component exposes two endpoints:

* The /playground endpoint: this simulates the dispel4py execution in DARE and prints the logs and output content directly to the user
* The /run-command endpoint: accepts any bash command, which is executed and returns the result to the user

### Use in notebook

* For the first endpoint, you need to execute the first steps as always: login, create workspace, register the workflow
* For the second endpoint, you need to provide the endpoint, the token, the command and the output file name if exists

### Update helper_functions

Add the below two methods in helper_functions:

* For the first endpoint
```python 
def debug_d4p(hostname, impl_id, pckg, workspace_id, pe_name, token, reqs=None, output_filename="output.txt",
              **kw):
    # Prepare data for posting
    data = {
        "impl_id": impl_id,
        "pckg": pckg,
        "wrkspce_id": workspace_id,
        "n_nodes": 1,
        "name": pe_name,
        "access_token": token,
        "output_filename": output_filename,
        "reqs": reqs if not (reqs is None) else "None"
    }
    d4p_args = {}
    for k in kw:
        d4p_args[k] = kw.get(k)
    data['d4p_args'] = d4p_args
    r = requests.post(hostname + '/playground', data=json.dumps(data))
    if r.status_code == 200:
        response = json.loads(r.text)
        if response["logs"]:
            print("Logs:\n========================")
            for log in response["logs"]:
                print(log)
        if response["output"]:
            print("Output content:\n==============================")
            for output in response["output"]:
                print(output)
    else:
        print('Playground returns status_code: \
                ' + str(r.status_code))
        print(r.text)
```

* For the second endpoint:
```python
import requests
import json

def exec_command(hostname, token, command, run_dir="new", output_filename="output.txt"):
    data = {
        "access_token": token,
        "command": command,
        "run_dir": run_dir,
        "output_filename": output_filename
    }

    r = requests.post(hostname + '/run-command', data=json.dumps(data))
    if r.status_code == 200:
        response = json.loads(r.text)
        if response["logs"]:
            print("Logs:\n========================")
            for log in response["logs"]:
                print(log)
        if response["output"]:
            print("Output content:\n==============================")
            for output in response["output"]:
                print(output)
        if response["run_dir"]:
            print("Run directory is: ")
            print(response["run_dir"])
    else:
        print('Playground returns status_code: \
                ' + str(r.status_code))
        print(r.text)
```

### Update the jupyter notebook

* For the /playground endpoint:

```python
F.debug_d4p(impl_id=impl_id, pckg="mysplitmerge_pckg", workspace_id=workspace_id, pe_name="mySplitMerge", 
            token=F.auth(), creds=creds, no_processes=6, iterations=1,
            reqs='https://gitlab.com/project-dare/dare-api/raw/master/examples/jupyter/requirements.txt')
```

* For the /run-command endpoint:

```python
    F.exec_command(PLAYGROUND_API_HOSTNAME, F.auth(), "pip install --user numpy")
```

Technical documentation of the component is also available [here](https://project-dare.gitlab.io/playground/)

## Semantic Data Discovery

The API documentation of the Semantic Data Discovery component is available in our [testbed environment](https://testbed.project-dare.eu/semantic-data)

# Dispel4py Documentation

dispel4py is a free and open-source Python library for describing abstract stream-based workflows for distributed data-intensive applications. It enables users to focus on their scientific methods, avoiding distracting details and retaining flexibility over the computing infrastructure they use.  It delivers mappings to diverse computing infrastructures, including cloud technologies, HPC architectures and  specialised data-intensive machines, to move seamlessly into production with large-scale data loads. The dispel4py system maps workflows dynamically onto multiple enactment systems, and supports parallel processing on distributed memory systems with MPI and shared memory systems with multiprocessing, without users having to modify their workflows.

## Dependencies

dispel4py has been tested with Python *2.7.6*, *2.7.5*, *2.7.2*, *2.6.6* and Python *3.4.3*, *3.6*, *3.7*.

The following Python packages are required to run dispel4py:

- networkx (https://networkx.github.io/)

If using the MPI mapping:

- mpi4py (http://mpi4py.scipy.org/)

## Installation

Clone this repository to your desktop. You can then install from the local copy to your python environment by calling:

```
python setup.py install
```

from the dispel4py root directory.

## Docker

The Dockerfile in the dispel4py root directory builds a Debian Linux distribution and installs dispel4py and OpenMPI.

```
docker build . -t dare-dispel4py
```

Start a Docker container with the dispel4py image in interactive mode with a bash shell:

```
docker run -it dare-dispel4py /bin/bash
```

For the EPOS use cases obspy is included in a separate Dockerfile `Dockerfile.seismo`:

```
docker build . -f Dockerfile.seismo -t dare-dispel4py-seismo
```

# Provenance in Dispel4py

lean_empty
-----------

.. code-block:: python

   clean_empty(d)

Utility function that given a dictionary in input, removes all the properties that are set to None.
It workes recursively through lists and nested documents

total_size
----------

.. code-block:: python

   total_size(o, handlers={}, verbose=False)

Returns the approximate memory footprint an object and all of its contents.

Automatically finds the contents of the following builtin containers and
their subclasses:  tuple, list, deque, dict, set and frozenset.
To search other containers, add handlers to iterate over their contents:

handlers = {SomeContainerClass: iter,
            OtherContainerClass: OtherContainerClass.get_elements}

write
-----

.. code-block:: python

   write(self, name, data)

Redefines the native write function of the dispel4py SimpleFunctionPE to take into account
provenance payload when transfering data.

getDestination_prov
-------------------

.. code-block:: python

   getDestination_prov(self, data)

When provenance is activated it redefines the native dispel4py.new.process getDestination function to take into account provenance information
when redirecting grouped operations.

commandChain
------------

.. code-block:: python

   commandChain(commands, envhpc, queue=None)

Utility function to execute a chain of system commands on the hosting oeprating system.
The current environment variable can be passed as parameter env.
The queue parameter is used to store the stdoutdata, stderrdata of each process in message

ProvenanceType
--------------

.. code-block:: python

   ProvenanceType(self)

A workflow is a program that combines atomic and independent processing elements
via a specification language and a library of components. More advanced systems
adopt abstractions to facilitate re-use of workflows across users'' contexts and application
domains. While methods can be multi-disciplinary, provenance
should be meaningful to the domain adopting them. Therefore, a portable specification
of a workflow requires mechanisms allowing the contextualisation of the provenance
produced. For instance, users may want to extract domain-metadata from a component
or groups of components adopting vocabularies that match their domain and current
research, tuning the level of granularity. To allow this level of flexibility, we explore
an approach that considers a workflow component described by a class, according to
the Object-Oriented paradigm. The class defines the behaviour of its instances as their
type, which specifies what an instance will do in terms of a set of methods. We introduce
the concept of *ProvenanceType*\ , that augments the basic behaviour by extending
the class native type, so that a subset of those methods perform the additional actions
needed to deliver provenance data. Some of these are being used by some of the preexisting
methods, and characterise the behaviour of the specific provenance type, some
others can be used by the developer to easily control precision and granularity. This approach,
tries to balance between automation, transparency and explicit intervention of the developer of a data-intensive tool, who
can tune provenance-awareness through easy-to-use extensions.

The type-based approach to provenance collection provides a generic *ProvenanceType* class
that defines the properties of a provenance-aware workflow component. It provides
a wrapper that meets the provenance requirements, while leaving the computational
behaviour of the component unchanged. Types may be developed as **Pattern Type** and **Contextual Type** to represent respectively complex
computational patterns and to capture specific metadata contextualisations associated to the produce output data.

The *ProvenanceType* presents the following class constants to indicate where the lineage information will be stored. Options include a remote
repository, a local file system or a *ProvenanceSensor* (experimental).


* \_SAVE_MODE\ *SERVICE='service'*
* \_SAVE_MODE\ *FILE='file'*
* \_SAVE_MODE\ *SENSOR='sensor'*

The following variables will be used to configure some general provenance capturing properties


* \_PROV\ *PATH*\ : When _SAVE_MODE\ *SERVICE* is chosen, this variable should be populated with a string indicating a file system path where the lineage will be stored.
* \_REPOS\ *URL*\ : When _SAVE_MODE\ *SERVICE* is chosen, this variable should be populated with a string indicating the repository endpoint (S-ProvFlow) where the provenance will be sent.
* \_PROV_EXPORT_URL: The service endpoint from where the provenance of a workflow execution, after being stored, can be extracted in PROV format.
* \_BULK\ *SIZE*\ : Number of lineage documents to be stored in a single file or in a single request to the remote service. Helps tuning the overhead brough by the latency of accessing storage resources.


getProvStateObjectId
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   ProvenanceType.getProvStateObjectId(self, name)

Return the id of a named object stored in the provenance state

apply_derivation_rule
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   ProvenanceType.apply_derivation_rule(self, event, voidInvocation, oport=None, iport=None, data=None, metadata=None)

In support of the implementation of a *ProvenanceType* realising a lineage *Pattern type*. This method is invoked by the *ProvenanceType* each iteration when a decision has to be made whether to ignore or discard the dependencies on the ingested stream
and stateful entities, applying a specific provenance pattern, thereby creating input/output derivations. The framework invokes this method every time the data is written on an output port (\ *event*\ : *write*\ ) and every
time an invocation (\ *s-prov:Invocation*\ ) ends (\ *event*\ : \_end_invocation\ *event*\ ). The latter can be further described by  the boolean parameter *voidInvocation*\ , indicating whether the invocation terminated with any data produced.
The default implementation provides a *stateless* behaviour, where the output depends only from the input data recieved during the invocation.

getInputAt
^^^^^^^^^^

.. code-block:: python

   ProvenanceType.getInputAt(self, port='input', index=None)

Return input data currently available at a specific *port*. When reading input of a grouped operator, the *gindex* parameter allows to access exclusively the data related to the group index.

addNamespacePrefix
^^^^^^^^^^^^^^^^^^

.. code-block:: python

   ProvenanceType.addNamespacePrefix(self, prefix, url)

In support of the implementation of a *ProvenanceType* realising a lineage *Contextualisation type*.
A Namespace *prefix* can be declared with its vocabulary *url* to map the metadata terms to external controlled vocabularies.
They can be used to qualify the metadata terms extracted from the *extractItemMetadata* function,
as well as for those terms injected selectively at runtime by the *write* method. The namespaces will be used
consistently when exporting the lineage traces to semantic-web formats, such as RDF.

extractItemMetadata
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   ProvenanceType.extractItemMetadata(self, data, port)

In support of the implementation of a *ProvenanceType* realising a lineage *Contextualisation type*.
Extracts metadata from the domain specific content of the data (s-prov:DataGranules) written on a components output *port*\ , according to a particular vocabulary.

ignorePastFlow
^^^^^^^^^^^^^^

.. code-block:: python

   ProvenanceType.ignorePastFlow(self)

In support of the implementation of a *ProvenanceType* realising a lineage **Pattern type**.

It instructs the type to ignore the all the inputs when the method \_apply_derivation\ *rule* is invoked for a certain event."

ignoreState
^^^^^^^^^^^

.. code-block:: python

   ProvenanceType.ignoreState(self)

In support of the implementation of a *ProvenanceType* realising a lineage **Pattern type**.

It instructs the type to ignore the content of the provenance state when the method \_apply_derivation\ *rule* is invoked for a certain event."

discardState
^^^^^^^^^^^^

.. code-block:: python

   ProvenanceType.discardState(self)

In support of the implementation of a *ProvenanceType* realising a lineage **Pattern type**.

It instructs the type to reset the data dependencies in the provenance state when the method \_apply_derivation\ *rule* is invoked for a certain event.
These will not be availabe in the following invocations."

discardInFlow
^^^^^^^^^^^^^

.. code-block:: python

   ProvenanceType.discardInFlow(self, wlength=None, discardState=False)

In support of the implementation of a *ProvenanceType* realising a lineage **Pattern type**.

It instructs the type to reset the data dependencies related to the component''s inputs when the method \_apply_derivation\ *rule* is invoked for a certain event.
These will not be availabe in the following invocations."

update_prov_state
^^^^^^^^^^^^^^^^^

.. code-block:: python

   ProvenanceType.update_prov_state(self, lookupterm, data, location='', format='', metadata={}, ignore_inputs=False, ignore_state=True, stateless=False, **kwargs)

In support of the implementation of a *ProvenanceType* realising a lineage *Pattern type* or inn those circumstances where developers require to explicitly manage the provenance information within the component''s logic,.

Updates the provenance state (\ *s-prov:StateCollection*\ ) with a reference, identified by a *lookupterm*\ , to a new *data* entity or to the current input. The *lookupterm* will allow developers to refer to the entity when this is used to derive new data.
Developers can specify additional *medatata* by passing a metadata dictionary. This will enrich the one generated by the *extractItemMetadata* method.
Optionally the can also specify *format* and *location* of the output when this is a concrete resource (file, db entry, online url), as well as instructing the provenance generation to 'ignore_input' and 'ignore_state' dependencies.

The *kwargs* parameter allows to pass an argument *dep* where developers can specify a list of data *id* to explicitly declare dependencies with any data in the provenance state (\ *s-prov:StateCollection*\ ).

write
^^^^^

.. code-block:: python

   ProvenanceType.write(self, name, data, **kwargs)

This is the native write operation of dispel4py triggering the transfer of data between adjacent
components of a workflow. It is extended by the *ProvenanceType* with explicit provenance
controls through the *kwargs* parameter. We assume these to be ignored
when provenance is deactivated. Also this method can use the lookup tags to
establish dependencies of output data on entities in the provenance state.

The *kwargs* parameter allows to pass the following arguments:


* *dep* : developers can specify a list of data *id* to explicitly declare dependencies with any data in the provenance state (\ *s-prov:StateCollection*\ ).
* *metadata*\ : developers can specify additional medatata by passing a metadata dictionary.
* _ignore\ *inputs*\ : instructs the provenance generation to ignore the dependencies on the current inputs.
* *format*\ : the format of the output.
* *location*\ : location of the output when this is a concrete resource (file, db entry, online url).

checkSelectiveRule
^^^^^^^^^^^^^^^^^^

.. code-block:: python

   ProvenanceType.checkSelectiveRule(self, streammeta)

In alignement with what was previously specified in the configure_prov_run for the Processing Element,
check the data granule metadata whether its properies values fall in a selective provenance generation rule.

checkTransferRule
^^^^^^^^^^^^^^^^^

.. code-block:: python

   ProvenanceType.checkTransferRule(self, streammeta)

In alignement with what was previously specified in the configure_prov_run for the Processing Element,
check the data granule metadata whether its properies values fall in a selective data transfer rule.

extractDataSourceId
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   ProvenanceType.extractDataSourceId(self, data, port)

In support of the implementation of a *ProvenanceType* realising a lineage *Pattern type*. Extract the id from the incoming data, if applicable,
to reuse it to identify the correspondent provenance entity. This functionality is handy especially when a workflow component ingests data represented by
self-contained and structured file formats. For instance, the NetCDF attributes Convention includes in its internal metadata an id that can be reused to ensure
the linkage and therefore the consistent continuation of provenance tracesbetween workflow executions that generate and use the same data.

AccumulateFlow
--------------

.. code-block:: python

   AccumulateFlow(self)

A *Pattern type* for a Processing Element (\ *s-prov:Component*\ ) whose output depends on a sequence of input data; e.g. computation of periodic average.

Nby1Flow
--------

.. code-block:: python

   Nby1Flow(self)

A *Pattern type* for a Processing Element (\ *s-prov:Component*\ ) whose output depends
on the data received on all its input ports in lock-step; e.g. combined analysis of multiple
variables.

SlideFlow
---------

.. code-block:: python

   SlideFlow(self)

A *Pattern type* for a Processing Element (\ *s-prov:Component*\ ) whose output depends
on computations over sliding windows; e.g. computation of rolling sums.

ASTGrouped
----------

.. code-block:: python

   ASTGrouped(self)

A *Pattern type* for a Processing Element (\ *s-prov:Component*\ ) that manages a stateful operator
with grouping rules; e.g. a component that produces a correlation matrix with the incoming
coefficients associated with the same sampling-iteration index

SingleInvocationFlow
--------------------

.. code-block:: python

   SingleInvocationFlow(self)

A *Pattern type* for a Processing Element (\ *s-prov:Component*\ ) that
presents stateless input output dependencies; e.g. the Processing Element of a simple I/O
pipeline.

AccumulateStateTrace
--------------------

.. code-block:: python

   AccumulateStateTrace(self)

A *Pattern type* for a Processing Element (\ *s-prov:Component*\ ) that
keeps track of the updates on intermediate results written to the output after a sequence
of inputs; e.g. traceable approximation of frequency counts or of periodic averages.

IntermediateStatefulOut
-----------------------

.. code-block:: python

   IntermediateStatefulOut(self)

A *Pattern type* for a Processing Element (\ *s-prov:Component*\ ) stateful component which produces distinct but interdependent
output; e.g. detection of events over periodic observations or any component that reuses the data just written to generate a new product

ForceStateless
--------------

.. code-block:: python

   ForceStateless(self)

A *Pattern type* for a Processing Element (\ *s-prov:Component*\ ). It considers the outputs of the component dependent
only on the current input data, regardless from any explicit state update; e.g. the user wants to reduce the
amount of lineage produced by a component that presents inline calls to the \_update_prov\ *state*\ , accepting less accuracy.

get_source
----------

.. code-block:: python

   get_source(object, spacing=10, collapse=1)

Print methods and doc strings.
Takes module, class, list, dictionary, or string.

injectProv
----------

.. code-block:: python

   injectProv(object, provType, active=True, componentsType=None, workflow={}, **kwargs)

This function dinamically extend the type of each the nodes of the graph
or subgraph with ProvenanceType type or its specialisation.

configure_prov_run
------------------

.. code-block:: python

   configure_prov_run(graph, provRecorderClass=None, provImpClass=<class 'dispel4py.provenance.ProvenanceType'>, input=None, username=None, workflowId=None, description=None, system_id=None, workflowName=None, workflowType=None, w3c_prov=False, runId=None, componentsType=None, clustersRecorders={}, feedbackPEs=[], save_mode='file', sel_rules={}, transfer_rules={}, update=False, sprovConfig=None, sessionId=None, mapping='simple')

In order to enable the user of a data-intensive application to configure the lineage metadata extracted from the execution of their
worklfows we adopt a provenance configuration profile. The configuration is used at the time of the initialisation of the workflow to prepare its provenance-aware
execution. We consider that a chosen configuration may be influenced by personal and community preferences, as well as by rules introduced by institutional policies.
For instance, a certain RI would require to choose among a set of contextualisation types, in order to adhere to
the infrastructure's metadata portfolio. Thus, a provenance configuration profile play
in favour of more generality, encouraging the implementation and the re-use of fundamental
methods across disciplines.

With this method, the users of the workflow provide general provenance information on the attribution of the run, such as *username*\ , *runId* (execution id),
*description*\ , *workflowName*\ , and its semantic characterisation *workflowType*. It allows users to indicate which provenance types to apply to each component
and the belonging conceptual provenance cluster. Moreover, users can also choose where to store the lineage (\_save\ *mode*\ ), locally in the file system or in a remote service or database.
Lineage storage operations can be performed in bulk, with different impacts on the overall overhead and on the experienced rapidity of access to the lineage information.


* **Configuration JSON**\ : We show here an example of the JSON document used to prepare a worklfow for a provenance aware execution. Some properties are described inline. These are defined by terms in the provone and s-prov namespaces.

.. code-block:: python

       {
               'provone:User': "aspinuso",
               's-prov:description' : "provdemo demokritos",
               's-prov:workflowName': "demo_epos",
               # Assign a generic characterisation or aim of the workflow
               's-prov:workflowType': "seis:preprocess",
               # Specify the unique id of the workflow
               's-prov:workflowId'  : "workflow process",
               # Specify whether the lineage is saved locally to the file system or remotely to an existing serivce (for location setup check the class prperties or the command line instructions section.)
               's-prov:save-mode'   : 'service'         ,
               # Assign the Provenance Types and Provenance Clusters to the processing elements of the workflows. These are indicated by the name attributed to their class or function, eg. PE_taper. The 's-prov:type' property accepts a list of class names, corrisponding to the types' implementation. The 's-prov:cluster' is used to group more processing elements to a common functional section of the workflow.
               's-prov:componentsType' :
                                  {'PE_taper': {'s-prov:type':["SeismoPE"]),
                                                's-prov:prov-cluster':'seis:Processor'},
                                   'PE_plot_stream':    {'s-prov:prov-cluster':'seis:Visualisation',
                                                      's-prov:type':["SeismoPE"]},
                                   'StoreStream':    {'s-prov:prov-cluster':'seis:DataHandler',
                                                      's-prov:type':["SeismoPE,AccumulateFlow"]}
                                   }}


* **Selectivity rules**\ : By declaratively indicating a set of Selectivity rules for every component ('s-prov:sel_rules'), users can respectively activate the collection
  of the provenance for particular Data elements or trigger transfer operations of the data to external locations. The approach takes advantage of the contextualisation
  possibilities offered by the provenance *Contextualisation types*. The rules consist of comparison expressions formulated in JSON that indicate the boundary
  values for a specific metadata term. Such representation is inspired by the query language and selectors adopted by a popular document store, MongoDB.
  These can be defined also within the configuration JSON introduced above.

Example, a Processing Element *CorrCoef* that produces lineage information only when the *rho* value is greater than 0:

.. code-block:: python

       { "CorrCoef": {
           "rules": {
               "rho": {
                   "$gt": 0
       }}}}


* ** Command Line Activation**\ : To enable proveance activation through command line dispel4py should be executed with specific command line instructions. The following command will execute a local test for the provenance-aware execution of the MySplitAndMerge workflow.**

.. code-block:: python

   dispel4py --provenance-config=dispel4py/examples/prov_testing/prov-config-mysplitmerge.json --provenance-repository-url=<url> multi dispel4py/examples/prov_testing/mySplitMerge_prov.py -n 10

* The following command instead stores the provenance files to the local filesystem in a given directory. To activate this mode, the property *s-prov:save_mode* of the configuration file needs to be set to 'file'.

.. code-block:: python

    dispel4py --provenance-config=dispel4py/examples/prov_testing/prov-config-mysplitmerge.json --provenance-path=/path/to/prov multi dispel4py/examples/prov_testing/mySplitMerge_prov.py -n 10
    
    
ProvenanceSimpleFunctionPE
--------------------------

.. code-block:: python

   ProvenanceSimpleFunctionPE(self, *args, **kwargs)

A *Pattern type* for the native  *SimpleFunctionPE* of dispel4py

ProvenanceIterativePE
---------------------

.. code-block:: python

   ProvenanceIterativePE(self, *args, **kwargs)

A *Pattern type* for the native  *IterativePE* Element of dispel4py
---
title: "FAQ"
---

# Easy Setup (Hugo + Netlify + Forestry)
Build your website with bigspring theme by following this easy steps (No Coding Required)

<a href="http://bit.ly/meghna-hugo-installation" target="_blank" title="meghna hugo installation" rel="nofollow"><img width="100%" src="https://user-images.githubusercontent.com/37659754/70844354-4028be00-1e6a-11ea-8d84-02e9a25e7db8.png"></a>

In this tutorial we will show you to make your website live without buying any hosting and touching a single line of code. We made this tutorial based on [meghna hugo](https://github.com/themefisher/meghna-hugo) but you can setup everithing like this.

### What you need !!

1. Git acccount (Ex: Github, Gitlab etc ) . In our case we use github.
2. [Netlify](https://bit.ly/netlify-account) account to host files and add custom domain .
3. [Forestry](https://bit.ly/forestry-account) account to maintain whole project without code.


### Step 1 : Fork or Clone repository

First we will fork this [bigspring](https://github.com/themefisher/bigspring-hugo-startup-theme) template.

### Step 2 : Add your repository in Forestry

Go to your [forestry](https://bit.ly/forestry-account)  account and click on `import your site now`. declare your config.toml file [`exampleSite`] and fill up basic settings .

**Or just click this button for one click installation** [![import to forestry](https://assets.forestry.io/import-to-forestryK.svg)](https://app.forestry.io/quick-start?repo=themefisher/bigspring-hugo-startup-theme&engine=hugo&version=0.60.1&config=exampleSite)

Now mark everything as done, then go to configuration to change the base url . You can put any url but this have to similar as netlify . So for now put a name which you are going to put in netlify as netlify subdomain.

### Step 3 : Setup and host website with Netlify

Here comes the last step . Go to your [netlify](https://bit.ly/netlify-account) account and click add new site . Choose your git repository to import your website in netlify .  And now you can see the forked `bigspring` theme. select it and follow the steps. Then go to `site settings` for change the site name and put your subdoamin name here what you puted on forestry as base url. save it and go to `deploy` from top menu, Wait a while and click on `site preview` or just simply go to the subdomain you puted as base url. **BOOM! Your site is live.** Now you can go to forestry and add, remove or customize every setting and content.

> If you face any issue regarding the installation feel free to onen [open a new issue](https://github.com/themefisher/bigspring-hugo-startup-theme/issues)


## Table of Contents

- [Demo](#demo)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Reporting Issues](#reporting-issues)
- [Technical Support or Questions](#technical-support-or-questions-(paid))
- [Licensing](#licensing)
- [More Hugo Themes](https://themefisher.com/hugo-themes/)

## Demo

[Live Preview](http://demo.themefisher.com/bigspring-hugo/).

## Quick Start
Quick start options:

- Clone the repo: `git clone https://github.com/themefisher/bigspring-hugo-startup-theme.git`.
- [Download from Github](https://github.com/themefisher/bigspring-hugo-startup-theme/archive/master.zip).
- [Download from gethugothemes website](https://gethugothemes.com/products/bigspring/).

## Installation
At the top we have shown an easy hugo installation. but still if you think you want to go with the traditional way then use the following commands:

```
$ git clone git@github.com:themefisher/bigspring-hugo-startup-theme.git
$ cd bigspring/exampleSite/
$ hugo server --themesDir ../..
```


## Reporting Issues

We use GitHub Issues as the official bug tracker for the **bigspringe**. Please Search [existing issues](https://github.com/themefisher/bigspring-hugo-startup-theme/issues). It’s possible someone has already reported the same problem.
If your problem or idea is not addressed yet, [open a new issue](https://github.com/themefisher/bigspring-hugo-startup-theme/issues/new)

## Technical Support or Questions (Paid)

If you have questions or need help integrating the product please [contact us](mailto:mehedi@themefisher.com) instead of opening an issue.

## Licensing

- Copyright 2020 Designed by [Themefisher](https://themefisher.com/) & Developed by [Gethugothemes](https://gethugothemes.com/)
- Licensed under MIT (https://github.com/themefisher/bigspring-hugo-startup-theme/blob/master/LICENSE)

### Donate Us (Author) 
This project needs you! If you would like to support this project's further development, the creator of this project or the continuous maintenance of this project, feel free to donate. Your donation is highly appreciated . Thank you!

* **[Donate $10](https://www.paypal.me/themefisher/10USD)**: Thank's for creating this project, here's a tea (or some juice) for you!
* **[Donate $20](https://www.paypal.me/themefisher/20USD)**: Wow, I am stunned. Let me take you to the movies!
* **[Donate $30](https://www.paypal.me/themefisher/30USD)**: I really appreciate your work, let's grab some lunch!
* **[Donate $40](https://www.paypal.me/themefisher/40USD)**: That's some awesome stuff you did right there, dinner is on me!
* **[Donate $50](https://www.paypal.me/themefisher/50USD)**: I really really want to support this project, great job!
* **[Donate $100](https://www.paypal.me/themefisher/100USD)**: You are the man! This project saved me hours (if not days) of struggle and hard work, simply awesome!
* **[Donate $1500](https://www.paypal.me/themefisher/1500USD)**: Go buddy, buy Macbook Pro for yourself!

Of course, you can also choose what you want to donate, all donations are awesome !


## Premium Themes

| [![Mega-Bundle-HUGO](https://gethugothemes.com/wp-content/uploads/edd/2019/09/Mega-Bundle-HUGO.png)](https://themefisher.com/products/hugo-mega-bundle/) | [![Phantop](https://gethugothemes.com/wp-content/uploads/edd/2019/06/Phantom.jpg)](https://gethugothemes.com/products/phantom-hugo-theme/) | [![redlab](https://gethugothemes.com/wp-content/uploads/edd/2019/09/redlab-hugo-thumbnail.jpg)](https://gethugothemes.com/products/redlab-hugo/) |
|:---:|:---:|:---:|
| **Hugo Mega Bundle**  | **Phantom**  | **Red Lab**  |
| [![northendlab](https://gethugothemes.com/wp-content/uploads/2019/11/Blogplate-Blog-Template.png)](https://gethugothemes.com/products/northendlab/) | [![Influencer](https://gethugothemes.com/wp-content/uploads/2019/11/Influencer.png)](https://gethugothemes.com/products/influencer-hugo/) | [![Vex](https://gethugothemes.com/wp-content/uploads/edd/2019/07/Vex.jpg)](https://gethugothemes.com/products/vex-hugo-theme/) |
| **Northendlab** | **Influencer** | **Vex** |
| [![Timer](https://gethugothemes.com/wp-content/uploads/edd/2019/07/Timer.jpg)](https://gethugothemes.com/products/timer-hugo-theme/) | [![Parsa](https://gethugothemes.com/wp-content/uploads/edd/2019/07/parsa-768x576.jpg)](https://gethugothemes.com/products/parsa-hugo-theme/) | [![all](https://gethugothemes.com/wp-content/uploads/2019/12/get-more-hugo-themes.png)](https://gethugothemes.com/shop/) |
| **Timer** | **Parsa** | **More Hugo Themes** |
---
title: "{{ replace .Name "-" " " | title }}"
date: {{ .Date }}
draft: true
---

# DARE platform testing

## General

The basic components of the DARE platform include: the dare login service, which uses Keycloak in the backend to 
authenticate the users, the two workflow registries, where users can register their workspaces and workflows and
then share it with other users, the Execution API, which depends on the two registries so as to retrieve a workflow and 
then execute it and finally the sprov API which is used so as to create/retrieve provenance logs. Below, 
we include some instructions on how to run the provided tests on those workflows.

## DARE login service

Although this service is tested through all the other tests, as no interaction with the platform can be performed
without a session token, we provide a separate test of the dare-login component. Before start testing, update
the credentials.yaml file only if you have a local installation of dare. In that case, use a username/password
of a user already registered in Keycloak and update the ip with the DARE platform's machine IP. Then, you need to find 
the port of the dare-login public service. To help you with this, we provide the dareports.sh script. Use it and
copy the port that is printed. Then paste it in the login_port field of the credentials.yaml

You can use the login test as a quick test, to check if the services are up and running. The login endpoint of the
DARE login service calls the d4p-registry, the cwl workflow registry, the execution API to initialize the user's
 environment as well as the Keycloak Service for authentication. The only main component that is not used in the login
 endpoint is the sprov API.

## CWL Workflow Registry component

Enter the workflow_registry folder to find the necessary files so as to test the Workflow Registry API.
If you have a local installation of DARE platform, you need to modify the credentials.yaml. Update the username/password fields
using a registered user in your local Keycloak. Then, you need to update the ip, login_port and workflow_port fields.
The ip field should contain the machine's IP while the login_port and workflow_port fields should contain the ports
of the public Kubernetes services for dare-login and workflow-registry. To find the ports use the dareports.sh script.
In case you do not have a local installation, the test will use the URL to the DARE's operational environment.

Once the credentials.yaml script is updated, execute the test_workflow_registry.py to test the endpoints of the
Workflow Registry API. When the script finishes, it will print the number of passed and failed tests. To test a CWL
execution go to the cwl_exec folder under the test_dare_components directory and execute the test_cwl.py. Of course,
modify again the yaml configuration file.

### CWL execution

The tests for the CWL execution are under the ```test_dare_components/cwl_exec``` directory.
The test_cwl.py script combines some functionality of the Workflow Registry API and the execution and monitoring 
functionality of the Execution API (test_cwl.py). 
The test_exec.py provides some additional tests on the Execution API once the CWL execution is finished. 

Keep in mind that in case you want to test a local DARE installation, you need to need to update the 
credentials.yaml file. You need a valid username/password from a registered user in you local Keycloak service
as well as the IP and ports of the services. Regarding the service ports, you can use the dareports.sh script
which prints the ports that you need for the cwl_exec test.

## Dispel4py Workflow Registry & Execution API

Enter the d4p_registry folder and prepare the credentials.yaml in a similar way as before. We provide a new dareports.sh
to find the d4p-registry port in your local installation. If you do not have a local installation, leave blank the ip,
login_port and workflow_port fields in the yaml file, so the test to use the operational DARE platform.

Use the test_d4p_registry to test the Dispel4py Registry API. Using this script, you can also test the execution of 
the dispel4py workflow. The dispel4py library installed is retrieved from [this repository](https://gitlab.com/project-dare/dispel4py) 
using the master branch.

For the testing of the execution result, use the test_exec.py. Modify the run_dir variable of the script with the run
directory that the execution step returned.

## Sprov API

To run the sprov tests, you need to enter the sprov container. To do so:

1. Use ```kubectl get pods``` command to list all the containers.
2. Find the sprov container and copy the name. Copy the simple sprov, not the sprov-db or sprov-cronjob or sprov-view.
3. Execute  ```kubectl exec -it <sprov name> -- apk add bash```
4. Enter the container: ```kubectl exec -it <sprov name> -- bash``` and:
    * cd ../test to enter the folder containing the tests.
    * Follow the instructions [here](https://gitlab.com/project-dare/s-ProvFlow/-/tree/master/provenance-api/src/test) to run manually the tests
 
The provenance API is used via the dispel4py test execution (d4p_registry folder).
To view the provenance traces use the sprov-viewer, i.e. http://IP:port/html/view.jsp or if you do not have a local
installation and you used the operational DARE platform you can enter [here](https://platform.dare.scai.fraunhofer.de/sprovflow-viewer/html/view.jsp).
You need of course to execute the above test in d4p-registry and then check the provenance traces.

Once on the S-ProvFlow viewer click on Open Run and in the popup, choose which of your executions you want to explore.
The viewer will show the workflow processes on the left panel with information such last active time,
data produced, execution node and messages such as logs and errors. Clicking on a process allows to see in detail what data has been produced.
Data are listed in the Data Product Panel. Each item includes the metadata, as specified in the workflow, its internal content, in case of larger
collections, and other general details, such as time of execution, process and associated messages.

Starting from a data element, it is possible to visualise the data-dependency graph. Here, you can explore the whole history (lineage)
of the data and processes that were involved in the generation of the data element. The navigation of the graph proceeds backwards,
from the data element of interest until the beginning of the analysis.
## Synchronize 3d-party repositories
```
git submodule update --init --recursive
```

## Spawn a Vagrant VM to install ansible
```
vagrant up
vagrant ssh
```


## Create SSH Keys to access nodes (optional).
```
ssh-keygen -t rsa
ssh-copy-id 
```

```
ansible-playbook -i inventory/hosts.ini --become --become-user=root --become-method=sudo --ask-sudo-pass -u user --extra-vars bootstrap_os=ubuntu cluster.yml
```
# Local Docker builds

Docker builds that are not pushed in a public registry
are built and pushed to `registry.gitlab.com/project-dare/dare-platform/`

# Work-in-progress Docker image for the Vulcanology use case using Fall3D.

To test manually, follow these steps:

* build image, e.g. `docker build -t exec-context-fall3dpy .`
* download static input data and unzip to current working directory: 
```
 wget https://owncloud.scai.fraunhofer.de/index.php/s/3RJqXi5k7HEA9TB/download
 unzip WP6-Fall3d.zip
```
* run image interactively using a bind mount for the input data
```
 docker run -it -v `pwd`/WP6-Fall3d:/input:ro exec-context-fall3dpy /bin/bash
```
* the rest of the commands are to be run inside the container.
* link the input folders to the working directory:
```
cd /home/mpiuser/docker/wp6_volcanology
ln -s /input/Cartopy_Data/ .
ln -s /input/TsupyDEM/ .
```
* export PYTHONPATH to fix relative imports 
```
export PYTHONPATH=.:$PYTHONPATH
```
* Start a project, e.g. with Stromboli case:
```
python f3_setup.py Eruption_test_stromboli Stromboli
```
* run download step:
```
dispel4py multi f3_dispel4py_download.py -n 10 -f project_parameters.json
```
* run Simulation step: 
```
python FALL3D_Executor.py
```
* create input_parameters.json with plotting parameters (see notebook in wp6)
```
python set-plotting-parameters.py
```
* run postprocessing step: 
```
dispel4py multi f3_dispel4py_plotting -n 10 -f input_parameters.json
```
* Output will be in the Simulations/.../Figures folder of your project (depending on
the name you chose when calling f3_setup.py). 
# Playground

The purpose of this component is to provide a DARE environment for test and debugging purposes. The component exposes two endpoints:

* The /playground endpoint: this simulates the dispel4py execution in DARE and prints the logs and output content directly to the user
* The /run-command endpoint: accepts any bash command, which is executed and returns the result to the user

## Use in notebook

* For the first endpoint, you need to execute the first steps as always: login, create workspace, register the workflow
* For the second endpoint, you need to provide the endpoint, the token, the command and the output file name if exists

## Update helper_functions

Add the below two methods in helper_functions:

* For the first endpoint

```python 
def debug_d4p(impl_id, pckg, workspace_id, pe_name, token, creds, reqs=None, output_filename="output.txt", **kw):
    # Prepare data for posting
    data = {
        "impl_id": impl_id,
        "pckg": pckg,
        "wrkspce_id": workspace_id,
        "n_nodes": 1,
        "name": pe_name,
        "access_token": token,
        "output_filename": output_filename,
        "reqs": reqs if not (reqs is None) else "None"
    }
    d4p_args = {}
    for k in kw:
        d4p_args[k] = kw.get(k)
    data['d4p_args'] = d4p_args
    r = requests.post(creds['PLAYGROUND_API_HOSTNAME'] + '/playground', data=json.dumps(data))
    if r.status_code == 200:
        response = json.loads(r.text)
        if response["logs"]:
            print("Logs:\n========================")
            for log in response["logs"]:
                print(log)
        if response["output"]:
            print("Output content:\n==============================")
            for output in response["output"]:
                print(output)
    else:
        print('Playground returns status_code: \
                ' + str(r.status_code))
        print(r.text)
```
* For the second endpoint:
```python
    def exec_command(hostname, username, command, run_dir="new", output_filename="output.txt"):

    data = {
        "username": username,
        "command": command,
        "run_dir": run_dir,
        "output_filename": output_filename
    }

    r = requests.post(hostname + '/run-command', data=json.dumps(data))
    if r.status_code == 200:
        response = json.loads(r.text)
        if response["logs"]:
            print("Logs:\n========================")
            for log in response["logs"]:
                print(log)
        if response["output"]:
            print("Output content:\n==============================")
            for output in response["output"]:
                print(output)
        if response["run_dir"]:
            print("Run directory is: ")
            print(response["run_dir"])
    else:
        print('Playground returns status_code: \
                ' + str(r.status_code))
        print(r.text)

```

## Update the jupyter notebook

* For the /playground endpoint:

```python
F.debug_d4p(impl_id=impl_id, pckg="mysplitmerge_pckg", workspace_id=workspace_id, pe_name="mySplitMerge", 
            token=F.auth(), creds=creds, no_processes=6, iterations=1,
            reqs='https://gitlab.com/project-dare/dare-api/raw/master/examples/jupyter/requirements.txt')
```

* For the /run-command endpoint:

```python
F.exec_command(PLAYGROUND_API_HOSTNAME, F.auth(), "pip install --user numpy")
```s-prov
======
**Version:** v1

### /data
---
##### ***GET***
**Description:** The data is selected by specifying a query string. Query parameters allow to search by attribution to a component or to an implementation, generation by a workflow execution and by combining more metadata and parameters terms with their min and max valuesranges. Mode of the search can also be indicated (mode ::= (OR | AND). It will apply to the search upon metadata and parameters values-ranges

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| usernames | query | csv list of users the Workflows Executons are associated with | No | string |
| terms | query | csv list of metadata or parameter terms. These relate positionally to the maxvalues and the minvalues | No | string |
| functionNames | query | csv list of functions the Data was generated with | No | string |
| wasAttributedTo | query | csv list of Component or Component Instances involved in the generation of the Data | No | string |
| minvalues | query | csv list of metadata or parameters minvalues. These relate positionally to the terms and the minvalues | No | string |
| rformat | query | unimplemented: format of the response payload (json,json-ld) | No | string |
| start | query | index of the starting item | Yes | integer |
| limit | query | max number of items expected | Yes | integer |
| maxvalues | query | csv list of metadata or parameters maxvalues. These relate positionally to the terms and the minvalues | No | string |
| formats | query | csv list of data formats (eg. mime-types) | No | string |
| types | query | csv list of data types | No | string |
| wasGeneratedBy | query | the id of the Invocation that generated the Data | No | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /data/filterOnAncestor
---
##### ***POST***
**Description:** Filter a list of data ids based on the existence of at least one ancestor in their data dependency graph, according to a list of metadata terms and their min and max values-ranges. Maximum depth level and mode of the search can also be indicated (mode ::= (OR | AND)

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| body | body |  | No | object |

**Responses**

| Code | Description |
| ---- | ----------- |

### /data/{data_id}
---
##### ***GET***
**Description:** Extract Data and their DataGranules by the Data id

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| data_id | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /data/{data_id}/derivedData
---
##### ***GET***
**Description:** Starting from a specific data entity of the data dependency is possible to navigate through the derived data or backwards across the element's data dependencies. The number of traversal steps is provided as a parameter (level).

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| level | query | level of depth in the data derivation graph, starting from the current Data | Yes | string |
| data_id | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /data/{data_id}/export
---
##### ***GET***
**Description:** Export of provenance information PROV-XML or RDF format. The S-PROV information returned covers the whole workflow execution or is restricted to a single data element. In the latter case, the graph is returned by following the derivations within and across runs. A level parameter allows to indicate the depth of the resulting trace

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| format | query | export format of the PROV document returned | No | string |
| rdfout | query | export rdf format of the PROV document returned | No | string |
| creator | query | the name of the user requesting the export | No | string |
| level | query | level of depth in the data derivation graph, starting from the current Data | Yes | string |
| data_id | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /data/{data_id}/wasDerivedFrom
---
##### ***GET***
**Description:** Starting from a specific data entity of the data dependency is possible to navigate through the derived data or backwards across the element's data dependencies. The number of traversal steps is provided as a parameter (level).

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| level | query | level of depth in the data derivation graph, starting from the current Data | Yes | string |
| data_id | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /instances/{instid}
---
##### ***GET***
**Description:** Extract details about a single instance or component by specifying its id. The returning document will indicate the changes that occurred, reporting the first invocation affected. It support the specification of a list of runIds the instance was wasAssociateFor, considering that the same instance could be used across multiple runs

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| start | query | index of the starting item | Yes | integer |
| limit | query | max number of items expected | Yes | integer |
| wasAssociateFor | query | cvs list of runIds the instance was wasAssociateFor (when more instances are reused in multiple workflow executions) | No | string |
| instid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /invocations/{invocid}
---
##### ***GET***
**Description:** Extract details about a single invocation by specifying its id

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| invocid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /summaries/collaborative
---
##### ***GET***
**Description:** Extract information about the reuse and exchange of data between workflow executions based on terms' valuesranges and a group of users. The API method allows for inclusive or exclusive (mode ::= (OR j AND) queries on the terms' values. As above, additional details, such as running infrastructure, type and name of the workflow can be selectively extracted by assigning these properties to a groupBy parameter. This will support the generation of grouped views

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| wasAssociatedWith | query | csv lis of Components involved in the Workflow Executions | No | string |
| usernames | query | csv list of users the Workflows Executons are associated with | No | string |
| terms | query | csv list of metadata or parameter terms. These relate positionally to the maxvalues and the minvalues | No | string |
| functionNames | query | csv list of functions that are executed by at least one workflow's components | No | string |
| level | query | level of depth in the data derivation graph, starting from the current Data | Yes | string |
| minvalues | query | csv list of metadata or parameters minvalues. These relate positionally to the terms and the minvalues | No | string |
| rformat | query | unimplemented: format of the response payload (json,json-ld) | No | string |
| maxvalues | query | csv list of metadata or parameters maxvalues. These relate positionally to the terms and the minvalues | No | string |
| formats | query | csv list of data formats (eg. mime-types) | No | string |
| clusters | query | csv list of clusters that describe and group one or more workflow's component | No | string |
| groupby | query | express the grouping of the returned data | No | string |
| types | query |  | No | string |
| mode | query | execution mode of the workflow in case it support different kind of concrete mappings (eg. mpi, simple, multiprocess, etc.. | No | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /summaries/workflowexecution
---
##### ***GET***
**Description:** Produce a detailed overview of the distribution of the computation, reporting the size of data movements between the workflow components, their instances or invocations across worker nodes, depending on the specified granularity level. Additional information, such as process pid, worker, instance or component of the workflow (depending on the level of granularity) can be selectively extracted by assigning these properties to a groupBy parameter. This will support the generation of grouped views

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| level | query | level of depth in the data derivation graph, starting from the current Data | Yes | string |
| mintime | query | minimum start time of the Invocation | No | string |
| groupby | query | express the grouping of the returned data | No | string |
| runId | query | the id of the run to be analysed | No | string |
| maxtime | query | maximum start time of the Invocation | No | string |
| maxidx | query | maximum iteration index of an Invocation | No | integer |
| minidx | query | minimum iteration index of an Invocation | No | integer |

**Responses**

| Code | Description |
| ---- | ----------- |

### /terms
---
##### ***GET***
**Description:** Return a list of discoverable metadata terms based on their appearance for a list of runIds, usernames, or for the whole provenance archive. Terms are returned indicating their type (when consistently used), min and max values and their number occurrences within the scope of the search

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| aggregationLevel | query | set whether the terms need to be aggreagated by runId, username or across the whole collection (all) | No | string |
| usernames | query | csv list of usernames | No | string |
| runIds | query | csv list of run ids | No | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions
---
##### ***GET***
**Description:** Extract documents from the bundle collection according to a query string which may include usernames, type of the workflow, the components the run wasAssociatedWith and their implementations. Data results' metadata and parameters can also be queried by specifying the terms and their min and max values-ranges and data formats. Mode of the search can also be indicated (mode ::= (OR j AND). It will apply to the search upon metadata and parameters values of each run

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| wasAssociatedWith | query | csv lis of Components involved in the Workflow Executions | No | string |
| usernames | query | csv list of users the Workflows Executons are associated with | No | string |
| terms | query | csv list of metadata or parameter terms. These relate positionally to the maxvalues and the minvalues | No | string |
| functionNames | query | csv list of functions that are executed by at least one workflow's components | No | string |
| minvalues | query | csv list of metadata or parameters minvalues. These relate positionally to the terms and the minvalues | No | string |
| rformat | query | unimplemented: format of the response payload (json,json-ld) | No | string |
| start | query | index of the starting item | Yes | integer |
| limit | query | max number of items expected | Yes | integer |
| maxvalues | query | csv list of metadata or parameters maxvalues. These relate positionally to the terms and the minvalues | No | string |
| formats | query | csv list of data formats (eg. mime-types) | No | string |
| clusters | query | csv list of clusters that describe and group one or more workflow's component | No | string |
| types | query |  | No | string |
| mode | query | execution mode of the workflow in case it support different kind of concrete mappings (eg. mpi, simple, multiprocess, etc.. | No | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions/insert
---
##### ***POST***
**Description:** Bulk insert of bundle or lineage documents in JSON format. These must be provided as encoded stirng in a POST request

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| body | body |  | No | object |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions/import
---
##### ***POST***
**Description:** Import of provenance output which is not yet mapped to the s-ProvFlowMongoDB format. 
The files provided in the archive will be mapped to s-ProvFlowMongoDB if they are in one of the supported formats.

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| archive | form | Zip archive of provenance output, which will be mapped to s-ProvFlowMongoDB and stored. Currently only files in the CWLProv format are supported | Yes | file |
| format | form | Format of the provenance output to be imported. | Yes | String |


**Responses**

| Code | Description |
| ---- | ----------- |



### /workflowexecutions/{run_id}/export
---
##### ***GET***
**Description:** Export of provenance information PROV-XML or RDF format. The S-PROV information returned covers the whole workflow execution or is restricted to a single data element. In the latter case, the graph is returned by following the derivations within and across runs. A level parameter allows to indicate the depth of the resulting trace

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| rdfout | query | export rdf format of the PROV document returned | No | string |
| creator | query | the name of the user requesting the export | No | string |
| format | query | export format of the PROV document returned | No | string |
| run_id | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions/{runid}
---
##### ***DELETE***
**Description:** Extract documents from the bundle collection by the runid of a WFExecution. The method will return input data and infomation about the components and the libraries used for the specific run

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| runid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

##### ***GET***
**Description:** Extract documents from the bundle collection by the runid of a WFExecution. The method will return input data and infomation about the components and the libraries used for the specific run

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| runid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions/{runid}/delete
---
##### ***POST***
**Description:** Delete a workflow execution trace, including its bundle and all its lineage documents

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| runid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions/{runid}/edit
---
##### ***POST***
**Description:** Update of the description of a workflow execution. Users can improve this information in free-tex

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| body | body |  | No | object |
| runid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

### /workflowexecutions/{runid}/showactivity
---
##### ***GET***
**Description:** Extract detailed information related to the activity related to a WFExecution (id). The result-set can be grouped by invocations, instances or components (parameter level) and shows progress, anomalies (such as exceptions or systems' and users messages), occurrence of changes and the rapid availability of accessible data bearing intermediate results. This method can also be used for runtime monitoring

**Parameters**

| Name | Located in | Description | Required | Schema |
| ---- | ---------- | ----------- | -------- | ---- |
| start | query | index of the starting item | Yes | integer |
| limit | query | max number of items expected | Yes | integer |
| level | query | level of aggregation of the monitoring information (component, instance, invocation, cluster) | No | string |
| runid | path |  | Yes | string |

**Responses**

| Code | Description |
| ---- | ----------- |

# dispel4py

dispel4py is a free and open-source Python library for describing abstract stream-based workflows for distributed data-intensive applications. It enables users to focus on their scientific methods, avoiding distracting details and retaining flexibility over the computing infrastructure they use.  It delivers mappings to diverse computing infrastructures, including cloud technologies, HPC architectures and  specialised data-intensive machines, to move seamlessly into production with large-scale data loads. The dispel4py system maps workflows dynamically onto multiple enactment systems, and supports parallel processing on distributed memory systems with MPI and shared memory systems with multiprocessing, without users having to modify their workflows.

## Dependencies

dispel4py has been tested with Python *2.7.6*, *2.7.5*, *2.7.2*, *2.6.6* and Python *3.4.3*, *3.6*, *3.7*.

The following Python packages are required to run dispel4py:

- networkx (https://networkx.github.io/)

If using the MPI mapping:

- mpi4py (http://mpi4py.scipy.org/)

## Installation

Clone this repository to your desktop. You can then install from the local copy to your python environment by calling:

```
python setup.py install
```

from the dispel4py root directory.

## Docker

The Dockerfile in the dispel4py root directory builds a Debian Linux distribution and installs dispel4py and OpenMPI.

```
docker build . -t dare-dispel4py
```

Start a Docker container with the dispel4py image in interactive mode with a bash shell:

```
docker run -it dare-dispel4py /bin/bash
```

For the EPOS use cases obspy is included in a separate Dockerfile `Dockerfile.seismo`:

```
docker build . -f Dockerfile.seismo -t dare-dispel4py-seismo
```# Execution API Documentation

## General

Execution API provides endpoints for multiple execution contexts:

* Dispel4py: dynamically creates containers to execute Dispel4py workflows
* Specfem: creates containers in a dynamic way to execute Specfem executable
* Cyclone: to be created

## API calls

* **create-folders:** endpoint to be used by a new user to generate his/her workspace (folders).
The folder structure follows the below logic:
    - All folders are stored in the shared file system: .i.e. /home/mpiuser/sfs/
    - Inside the above folder, a new directory is generated based on the given username, 
    i.e. /home/mpiuser/sfs/username/
    - The user's folder contains: one folder for uploads (named uploads), one folder for the
    execution (named runs) and one folder for the test executions (named debug).
* **run-d4p:** endpoint to execute dispel4py workflows. It is based on a docker image (can be found in 
DARE platform repository, named exec-context-d4p), which creates new containers that execute the given
workflow. Data to be provided:
    - PE implementation id
    - Workspace id
    - package name
    - PE name
    - number of containers to be created
    - username
    - requirements: txt file containing the requirements of the PE
    - token for authentication
* **run-specfem:** endpoint to execute specfem software. The containers are generated to execute it and
are based in a docker image, which can be found in the platform (exec-context-specfem). Data to be
provided:
    - PE implementation id
    - Workspace id
    - package name
    - PE name
    - username
    - token for authentication
    - output filename
    - requirements: txt file containing the requirements of the PE
 * **run-cyclone:**
 * **upload:** API endpoint to upload files in the DARE platform. All uploaded files are stored in the
 uploads directory of each user (e.g. "../username/uploads")
 * **my-files:** API endpoint to list all the folders inside the three main directories of a user. The
 user directory is named after its username and inside it the three directories are: uploads, runs and
 debug (which is used for testing purposes)
 * **list:** API endpoint which lists all the files inside a specific directory. Use the previous endpoint
 to get the specific directory you need and then use this endpoint to view the files inside the directory.
 * **download:** API endpoint, used to download a specific file. Use the previous endpoint to get the 
 desired file and then use this endpoint to download it locally.
 * **send2drop:** API endpoint which uploads a file in B2Drop
 * **my-pods:** API endpoint to list all the execution containers related to a specific MPI job
