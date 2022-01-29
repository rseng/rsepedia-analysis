# Not yet released
- Xenonflow now cleans up the remote directory when a job completes
- Xenonflow now cleans up the output directory correctly when deleting a job
- Added an option to the sourceDirectory config that controls deleting the input
  to a job when it is completed.

# Version 1.0.1
- Adding cwl compliance check and fixed a bunch of related bugs
- Reduced logging output
- Fixed various bugs
- Fixed a number of crashes related to the tests ran by the compliance test.
- Fixing output bug when parameter has # in its id
- Fixing file and directory array location wrangling
- Fixed directory tests and staging
- Fixing directory path staging for cwl check: nr85 directory_input_docker

# Version 1.0
- Various bugfixes
- Fixed uri's supplied by xenonflow when runing behind a proxy server
- Updated admin interface to connect to the backend using its current location
- You can now add xenonflow_jobid or xenonflow_jobname as input to a workflow to be run
  xenonflow will then supply these values automatically.
- Added XENONFLOW_FILES environment variable for use in the config.yml

# Version 1.0-rc1
- Split SourceFileSystem setting into sourceFileSystem for inputs and cwlFileSystem for workflow storage
- Added check on job submission whether the referenced workflow exists on the cwlFileSystem
- Added /workflows api that supplies a list of available workflows in the cwlFileSystem
- Upgrade admin interface to Angular 11

# Version 0.4-process
- Upgrade to Xenon 3

# Version 0.3-alpha
- Fixed pending jobs throwing an error

# Version 0.2-alpha
Output file serving and cwl parsing updates
- Output files are now served if they are stored on the local filesystem
  * This requires a new config block with a parameter "hosted" set to true:
  targetFileSystem:
   adaptor: file
   location: /home/bweel/Documents/projects/xenon-flow/output/
   hosted: true

- Output file paths are now set to this location if used

- Parsing of cwl has been updated to support Maps in addition to arrays for inputs and outputs

# Version 0.1-alpha
- First pre-release
- Updated to Xenon 2.1.0
[![Build and Test Xenonflow](https://github.com/xenon-middleware/xenon-flow/actions/workflows/build.yml/badge.svg?branch=master)](https://github.com/xenon-middleware/xenon-flow/actions/workflows/build.yml)
[![DOI](https://zenodo.org/badge/63334137.svg)](https://zenodo.org/badge/latestdoi/63334137)

##### CWL Compliance v1.0

![Required](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/required.json?icon=commonwl)
![Workflow](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/workflow.json?icon=commonwl)
![CommanLine](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/command_line_tool.json?icon=commonwl)
![Docker](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/docker.json?icon=commonwl)
![Badge](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/env_var.json?icon=commonwl)
![Badge](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/expression_tool.json?icon=commonwl)
![Badge](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/initial_work_dir.json?icon=commonwl)
![Badge](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/inline_javascript.json?icon=commonwl)
![Badge](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/multiple_input.json?icon=commonwl)
![Badge](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/resource.json?icon=commonwl)
![Badge](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/scatter.json?icon=commonwl)
![Badge](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/schema_def.json?icon=commonwl)
![Badge](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/shell_command.json?icon=commonwl)
![Badge](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/step_input.json?icon=commonwl)
![Badge](https://badgen.net/https/raw.githubusercontent.com/xenon-middleware/xenon-flow/gh-pages/badges/subworkflow.json?icon=commonwl)


# Xenonflow
Run CWL workflows using Xenon through a REST api.

# Usage:
The following diagram shows a rough overview of the interaction when using xenonflow.
The overview shows 3 file systems, all of which can be configured for Xenonflow (See the quick-start guide below).
1. The input filesystem, this should contain all the input files needed for running the cwl workflow
2. The cwl filesystem, this filesystem should contain the cwl workflows you want to run with Xenonflow.
3. The output filesystem, this is where xenonflow will put the output of the workflows.

On the right you can see a compute resource: Xenonflow can be configured to run on a number of computing backends, including the local machine, to actually execute the cwl workflow.

Before making a call to the Xenonflow REST API make sure the data is available on the input filesystem and the workflow is available on the cwl filesystems.
The rest call will return a JSON object which contains some information on the job you just submitted. Such as its current state, what input was provided, a uri to the job (for instance to poll the server for new states) and a uri to the log of the job.

After the workflow is completed the results will be available in the target filesystem

![Xenonflow Usage Pattern](https://user-images.githubusercontent.com/16776108/116083884-718bdd80-a69d-11eb-982a-7351b1e586f3.png)


# Quick-start
## 1. Install the dependencies:
 - Java 11
 - a cwl runner

For the cwl runner you can use the reference implementation called cwltool.
It can be installed by
```bash
pip install cwltool
```
You may need to use pip3 on some systems.
For a full list of available cwl runners check https://www.commonwl.org/#Implementations

After installing the cwl runner it is a good idea to double check that your workflow can be run
using the runner on the command line

## 2. Download Xenon-flow
```
wget https://github.com/xenon-middleware/xenon-flow/releases/V1.0
unzip xenonflow-1.0.zip
```

## 3. Configure the server
Configuration of the server is done by editing the `XENONFLOW_HOME/config/config.yml` file.
As well as the `XENONFLOW_HOME/config/application.properties`.

By default it is set the use the local file system as the source and the local
computer to run workflows.

For information on which filesystems and schedulers can be used refer to the xenon documentation: https://xenon-middleware.github.io/xenon/versions/3.1.0/javadoc/.

### config.yml
Xenon-flow configuration consists of
1. `sourceFileSystem`: Any filesystem supported by Xenon can be used here
2. `targetFileSystem`: Any filesystem supported by Xenon can be used here
3. `cwlFileSystem`: Any filesystem supported by Xenon can be used here
4. `ComputeResources`: A map of compute resource descriptions for Xenon
5. Each resource has the following settings:
    1. `cwlCommand`: A script to run the cwl runner, allowing for python environments to be started first.
    	* Default:
    		```![xenonflow (3)](https://user-images.githubusercontent.com/16776108/116083865-6d5fc000-a69d-11eb-81ca-5e82fa2c6727.png)

    		#!/usr/bin/env bash

    		cwltool $@
    		```
    2. `scheduler`: A Xenon scheduler
    3. `filesystem` A Xenon filesystem
    4. Both the scheduler and filesystem have to following format:
        1. `adaptor`: The name of the Xenon adaptor (for instance slurm for scheduler or sftp for filesystem)
        2. `location`: The URI for the resource
        3. `credential`: Optional credentials (if not supplied the current user and ~/.ssh/id_rsa is used)
        	1. `user`: Username
        	2. `password`: Password in base64 encoded
        4. `properties`: Optional properties (usually not needed)

#### Environment Variables
There are two environment variables that can be set in your environement which can then be
used in the config.yml file: `XENONFLOW_FILES` and `XENONFLOW_HOME`.

 

### application.properties
The application.properties needs configuration for the following things:
1. api-key
	1. `xenonflow.http.auth-token-header-name` controls the name of the header that holds the api key
	2. `xenonflow.http.auth-token` the value of the api key. IMPORTANT you should really change this one
2. The Database Configuration.
	* These settings should be changed!
    	1. `spring.datasource.username` The database username
    	2. `spring.datasource.password` The database password3.
	* The following settings can be left as is.
    	1. `server.port` The port for the server to run on.
    	2. `local.server.address=localhost` The servername.
    	3. `server.http.interface` Set up the server to be publicaly available by setting this to 0.0.0.0


## 4. Start the server
The following command will run the server.
```
./bin/xenonflow
```

## 5. Run a workflow
Put the workflow and any input files and directories to into the location as defined by the `sourceFileSystem` in the config. For instance when using a webdav server, upload the files there.

Send a POST http request with a job description to the server.

### Example:

Given the echo command-line-tool (in yaml notation):

```yaml
cwlVersion: v1.0
class: CommandLineTool
inputs:
  - id: inp
    type: string
    inputBinding: {}

outputs:
  - id: out
    type: string
    outputBinding:
      glob: out.txt
      loadContents: true
      outputEval: $(self[0].contents)

baseCommand: echo
stdout: out.txt
```

The job description looks something like the following.

Note that the input map contains a key `inp` which refers to the corresponding input of the echo command-line-tool.

```json
{
    "name": "My First Workflow",
    "workflow": "cwl/echo.cwl",
    "input": {
        "inp": "Hello CWL Server!"
    }
}
```

```bash
curl -X POST -H "Content-Type: application/json" -H "api-key: <insert api key here>" -d '{"name": "My First Workflow","workflow": "cwl/echo.cwl","input": {"inp": "Hello CWL Server!"}}' http://localhost:8080/jobs
```

### Using the jobid or jobname in a workflow
If you need access to the jobid generated by xenonflow, or the jobname that was used to submit the workflow
then you can add them as inputs to your cwl file as parameters with the ids `xenonflow_jobid` and `xenonflow_jobname` respectively.
Xenonflow will then automatically inject the values into the job-order.json as input to the cwl file.

For example the following cwl file will echo the xenonflow_jobid and xenonflow_jobname:
```yaml
cwlVersion: v1.0
class: CommandLineTool
inputs:
  - id: xenonflow_jobid
    type: string
    inputBinding:
      position: 1
  - id: xenonflow_jobname
    type: string
    inputBinding:
      position: 2

outputs:
  - id: out
    type: string
    outputBinding:
      glob: out.txt
      loadContents: true
      outputEval: $(self[0].contents)

baseCommand: echo
stdout: out.txt
```



### Running Xenonflow behind a proxy server
We recommend running xenonflow behind a proxy server. Both nginx and apache httpd are good candidates for this. In addition both nginx and apache httpd can act as webdav servers which xenonflow can use as a sourceFileSystem.

Doing this requires no changes to the configuration of xenonflow as long as the correct X-forwarded-* headers are set in the proxy server.

To ensure that xenonflow returns the correct uri's for the jobs you should set the following headers:
* X-Forwarded-Host
* X-Forwarded-Server
* X-Forwarded-Proto
* X-Forwarded-Port
* X-Forwarded-Prefix

Below is an example location from a nginx config that correctly proxies a xenonflow instance running at localhost:8080
```nginx
...
location /api/ {

    include cors;
    proxy_pass http://localhost:8080/;
    proxy_redirect off;
    proxy_set_header Host $host;
    proxy_set_header X-Forwarded-Host $host;
    proxy_set_header X-Forwarded-Server $host;
    proxy_set_header X-Forwarded-Proto http;
    proxy_set_header X-Forwarded-Port $server_port;
    proxy_set_header X-Forwarded-Prefix /api/;
}
...
```

### Running Xenonflow in SSL mode
To run xenonflow in ssl (https) mode you can follow the following steps:
1. Please read https://dzone.com/articles/spring-boot-secured-by-lets-encrypt for setup using Letsencrypt.
2. You should now have a certificate with a private key store
3. You should now set the following properties in the application.properties file:
   1. `server.ssl.enabled=true` Enable ssl encryption in the server
   2. `server.ssl.key-store-type` The keystore type (spring boot only supports PKCS12).
   3. `server.ssl.key-store` The store for the certificate files.
   4. `server.ssl.key-store-password` The password to the key store.
   5. `server.ssl.key-alias` The alias as given to the keystore.

### Cleaning Input Data
*Warning:* This will delete the input data on the source directory. It is recommended to set
the input filesystem to a different location than the cwl and output filesystems so files are
not lost by accident.

You can have xenonflow clean up the input files after a task has run by setting the `clearOnJobDone` parameter
to true in the sourceFileSystem.

i.e.

```yaml
sourceFileSystem:
   adaptor: file
   location: ${XENONFLOW_FILES}/input
   clearOnJobDone: true
```
# XenonflowAdmin

This project was generated with [Angular CLI](https://github.com/angular/angular-cli) version 6.1.2.

## Development server

Run `ng serve` for a dev server. Navigate to `http://localhost:4200/`. The app will automatically reload if you change any of the source files.

## Code scaffolding

Run `ng generate component component-name` to generate a new component. You can also use `ng generate directive|pipe|service|class|guard|interface|enum|module`.

## Build

Run `ng build` to build the project. The build artifacts will be stored in the `dist/` directory. Use the `--prod` flag for a production build.

## Running unit tests

Run `ng test` to execute the unit tests via [Karma](https://karma-runner.github.io).

## Running end-to-end tests

Run `ng e2e` to execute the end-to-end tests via [Protractor](http://www.protractortest.org/).

## Further help

To get more help on the Angular CLI use `ng help` or go check out the [Angular CLI README](https://github.com/angular/angular-cli/blob/master/README.md).
