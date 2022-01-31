# CHANGELOG

This is a manually generated log to track changes to the repository for each release. 
Each section should include general headers such as **Implemented enhancements** 
and **Merged pull requests**. Critical items to know are:

 - renamed commands
 - deprecated / removed commands
 - changed defaults
 - backward incompatible changes
 - migration guidance
 - changed behaviour

The versions coincide with releases on pypi.

## [0.0.x](https://github.com/singularityhub/singularity-compose/tree/master) (0.0.x)
 - fix a bug triggered when using startoptions in conjunction with network=false (0.0.15)
 - bind volumes can handle tilde expansion (0.0.14)
 - fix module import used by check command (0.0.13)
 - adding jsonschema validation and check command (0.0.12)
   - implement configuration override feature
   - implement `--preview` argument for the `check` command 
 - add network->enable option on composer file (0.1.11)
 - add network->allocate_ip option on composer file (0.1.10)
 - version 2.0 of the spec with added fakeroot network, start, exec, and run options (0.1.0)
   - stop option added (equivalent functionality to down)   
   - spython version 0.1.0 with Singularity 3.x or greater required
 - removed check for sudo when adding network flags (0.0.21)
 - singularity-compose down supporting timeout (0.0.20)
 - command, ability to associate arguments to the instance's startscript (0.0.19)
 - depends\_on, check circular dependencies at startup and shutdown in reverse order (0.0.18)
 - resolv.conf, etc.hosts generated if needed, network disabled non-sudo users (0.0.17)
 - resolv.conf needs to bind by default (0.0.16)
 - adding run command (0.0.15)
 - ensuring that builds are streamed (0.0.14)
 - adding more build options to build as build-flags (0.0.13)
 - when not using sudo, need to set --network=none, and catching exec error (0.0.12)
 - pyaml version should be for newer kind, and still check for attribute (0.0.11)
 - alpha release with simple (single container) working example (0.0.1)
 - dummy release (0.0.0)
# Singularity Compose

[![status](http://joss.theoj.org/papers/1fc2593b11b5e18df12efb59ed8757a0/status.svg)](http://joss.theoj.org/papers/1fc2593b11b5e18df12efb59ed8757a0)
[![DOI](https://zenodo.org/badge/188852712.svg)](https://zenodo.org/badge/latestdoi/188852712)

This is a simple orchestration library for Singularity containers, akin to
Docker Compose. See the [specification](https://singularityhub.github.io/singularity-compose/#/spec/spec-1.0) 
and the [documentation](https://singularityhub.github.io/singularity-compose) for
details, or more examples below.

## Examples

See our [Singularity Compose Examples](https://www.github.com/singularityhub/singularity-compose-examples) repository for
finding or contributing examples for using scompose.
---
title: 'Singularity Compose: Orchestration for Singularity Instances'
tags:
  - containers
  - singularity
  - linux
  - orchestration
authors:
 - name: Vanessa Sochat
   orcid: 0000-0002-4387-3819
   affiliation: 1
affiliations:
 - name: Stanford University Research Computing, Stanford University, Stanford, CA 94305
   index: 1
date: 24 June 2019
bibliography: paper.bib
---

# Summary

Singularity Compose is an orchestration tool for Singularity container instances.

![Singularity Compose](singularity-compose.png)

The Singularity container technology started to become popular in 2016,
as it offered a more secure option to run encapsulated environments [@Kurtzer2017-xj].
Traditionally, this meant that Singularity users could run a script built into the container
(called a runscript), execute a custom command, or shell into a container. 
Unlike Docker [@Merkel2014-da], these basic interactions simply interacted with processes in the 
foreground (e.g., running a script and exiting) and were not appropriate to run 
background services. This was a task for container instances [@SingularityInstances].

A container instance [@SingularityInstances] equates to running a container in a detached or
daemon mode. Instances allow for running persistent services in the background,
and then interaction with these services from the host and other containers.
Examples of services include databases, web servers, and associated applications
that interact with them. While a container technology can provide command line
and other programmatic interfaces for interaction with instances, what is also needed
is a configuration file for orchestration and customization of several instances.
For sibling container technology Docker, Docker Compose [@DockerCompose] was developed 
for this purpose. For local and production usage, the user could create a `docker-compose.yml` 
file to define services, volumes, ports exposed, and other customizations to networking and environment
[@DockerCompose]. Notably, there was strong incentive for the development of such a tool,
because Docker Compose existed before Kubernetes was available in the middle of 2015 [@Wikipedia_contributors2019-bw].

No equivalent orchestration tool was created for Singularity container
instances. While Singularity has empowered enterprise users to run 
services via platforms such as Kubernetes [@Meyer2019-sd], these platforms come
with privilege. It is often the case that a production Kubernetes cluster is not 
readily available to a user via his or her institution, or that the user 
cannot pay a cloud provider to deploy one. However, this does not imply that
a non enterprise user (e.g., an open source developer
or academic) would not benefit from such an orchestration tool. Unfortunately,
since the current trend and strongest potential for making profits is centered
around encouraging usage of enterprise tools like Kubernetes [@Wikipedia_contributors2019-bw],
there is not any urgent incentive on part of the provider companies to 
invest in a non-enterprise orchestration tool. It is logical, rational, and
understandable that companies exist to make profit, and must make profit
to exist. As the need is unfulfilled, it is the responsibility of the open source community to step up.


## Singularity Compose

The solution for orchestration of container instances from the open source
community is Singularity Compose [@SingularityCompose]. Singularity Compose 
is software for non enterprise users to easily create a configuration file to 
control creation and interaction of Singularity container instances.
It allows for the creation of a `singularity-compose.yml` file, in which
the user can define one or more container services, optionally with exposed ports
and volumes on the host. The user can easily define a container binary
to build or pull from a remote resource, along with custom scripts to
run after creation of the instances. Singularity Compose handles designation
of addresses on a local bridge network for each container, and creation of
resource files to bind to the containers to "see" one another. 
Importantly, by way of adding a Singularity Compose to a repository,
a user is ensuring not just reproducibility of a container recipe, but also
reproducibility of it's build and creation of services. For example, a simplified
version of a sequence of steps to build two containers and bring them up
as instances might look like this:

```bash
$ sudo singularity build app/app.sif app/Singularity
$ sudo singularity build nginx/nginx.sif nginx/Singularity.nginx

$ singularity instance start \
   --bind nginx.conf:/etc/nginx/conf.d/default.conf \
   --bind nginx/uwsgi_params.par:/etc/nginx/uwsgi_params.par \
   --bind nginx/cache:/var/cache/nginx \
   --bind nginx/run:/var/run \
   --bind app:/code \
   --bind static:/var/www/static \
   --bind images:/var/www/images \
   --bind etc.hosts:/etc/hosts \
   --net --network-args "portmap=80:80/tcp" --network-args "IP=10.22.0.2" \
   --hostname nginx --writable-tmpfs nginx/nginx.sif nginx

$ singularity instance start \
   --bind app:/code \
   --bind static:/var/www/static \
   --bind images:/var/www/images \
   --bind etc.hosts:/etc/hosts \
   --net --network-args "portmap=8000:8000/tcp" --network-args "IP=10.22.0.3" \
   --hostname app --writable-tmpfs app/app.sif app

$ singularity instance list
```

This is a complicated set of commands. In the above, we
first build the two containers. There are no checks here if the recipes
exist, or if the containers themselves already exist.
We then start instances for them. If we save these commands in a file,
we need to consistently hard code the paths to the container binaries,
along with the ip addresses, hostnames, and volumes. There are no checks
done before attempting the creation if the volumes meant to be bound
actually exist. We also take for granted that we've already generated an 
`etc.hosts` file to bind to the container at `/etc/hosts`, which will
define the container instances to have the same names supplied with `--hostname`. 
For the networking, we have to be mindful of the default bridge provided by Singularity, 
along with how to specify networking arguments under different conditions. 
This entire practice is clearly tedious. For a user to constantly need to generate and then
re-issue these commands, it's not a comfortable workflow. However, 
with Singularity Compose, the user writes a `singularity-compose.yml`
file once:

```yaml
version: "1.0"
instances:

  nginx:
    build:
      context: ./nginx
      recipe: Singularity.nginx
    volumes:
      - ./nginx.conf:/etc/nginx/conf.d/default.conf
      - ./uwsgi_params.par:/etc/nginx/uwsgi_params.par
      - ./nginx/cache:/var/cache/nginx
      - ./nginx/run:/var/run
    ports:
      - 80:80
    depends_on:
      - app
    volumes_from:
      - app

  app:
    build:
      context: ./app
    volumes:
      - ./app:/code
      - ./static:/var/www/static
      - ./images:/var/www/images
    ports:
      - 8000:8000
```

And then can much more readily see and reproduce generation of the same services.
The user can easily build all non-existing containers, and bring up all services
with one command:

```bash
$ singularity-compose up
```

And then easily bring services down, restart, shell into a container, execute
a command to a container, or run a container's internal runscript.

```bash
$ singularity-compose down                   # stop services
$ singularity-compose restart                # stop and start services
$ singularity-compose shell app              # shell into an instance
$ singularity-compose exec app "Hello!"      # execute a command
$ singularity-compose run app                # run internal runscript
```

These interactions greatly improve both reproducibility and running of
any development workflow that is not appropriate for an enterprise cluster but
relies on orchestration of container instances.

For the interested reader, the complete documentation for Singularity Compose [@SingularityCompose]
is provided, along with the code on GitHub [@SingularityComposeGithub]. For 
additional walkthroughs and complete examples, we direct the reader to the examples 
repository, also on GitHub [@SingularityComposeExamples]. Contribution by way
of additional examples, questions, or requests for development of a new example
are appreciated and welcome.


# References
# Commands

The following commands are currently supported. Remember, you must be in the 
present working directory of the compose file to reference the correct instances.

## check

To do a sanity check of your singularity-compose.yml, you can use `singularity-compose check`

```bash
$ singularity-compose check
singularity-compose.yml is valid.

$ singularity-compose -f singularity-compose.yml \ 
          -f singularity-compose.override.yml check 
singularity-compose.yml is valid.
singularity-compose.override.yml is valid.
```

To view the combined compose files you can use `--preview`.

```bash
$ singularity-compose -f singularity-compose.yml \ 
          -f singularity-compose.override.yml check  --preview
          
version: '2.0'
instances:
  cvatdb:
    start:
      options:
      - containall
    network:
      enable: false
    volumes:
    - ./recipes/postgres/env.sh:/.singularity.d/env/env.sh
    - ./volumes/postgres/conf:/opt/bitnami/postgresql/conf
    - ./volumes/postgres/tmp:/opt/bitnami/postgresql/tmp
    - /home/vagrant/postgres_data:/bitnami/postgresql
    build:
      context: .
      recipe: ./recipes/postgres/main.def
      options:
      - fakeroot

```

## build

Build will either build a container recipe, or pull a container to the
instance folder. In both cases, it's named after the instance so we can
easily tell if we've already built or pulled it. This is typically
the first step that you are required to do in order to build or pull your
recipes. It ensures reproducibility because we ensure the container binary
exists first.

```bash
$ singularity-compose build
```

The working directory is the parent folder of the singularity-compose.yml file.
If the build requires sudo (if you've defined sections in the config that warrant
setting up networking with sudo) the build will instead give you an instruction
to run with sudo.

## up

If you want to both build and bring them up, you can use "up." Note that for
builds that require sudo, this will still stop and ask you to build with sudo.

```bash
$ singularity-compose up
```

### resolv.conf

By default, singularity-compose will generate a resolv.conf for you to bind
to the container, instead of using the host. It's a simple template
that uses Google nameservers:

```
# This is resolv.conf generated by singularity-compose. It is provided
# to provide Google nameservers. If you don't want to have it generated
# and bound by default, use the up --no-resolv argument.

nameserver 8.8.8.8
nameserver 8.8.4.4
```

If you want to disable this:

```bash
$ singularity-compose up --no-resolv
Creating app
```

## create

Given that you have built your containers with `singularity-compose build`,
you can create your instances as follows:

```bash
$ singularity-compose create
```

Akin to "up," you can also disable the generation of the resolv.conf.

```bash
$ singularity-compose create --no-resolv
Creating app
```


## restart

Restart is provided as a convenience to run "down" and then "up." You can
specify most of the same arguments as create or up.

```bash
$ singularity-compose restart --no-resolv
Stopping app
Creating app
```


## ps

You can list running instances with "ps":

```bash
$ singularity-compose ps
INSTANCES  NAME PID     IMAGE
1           app	6659	app.sif
2            db	6788	db.sif
3         nginx	6543	nginx.sif
```

## shell

It's sometimes helpful to peek inside a running instance, either to look at permissions,
inspect binds, or manually test running something.
You can easily shell inside of a running instance:

```bash
$ singularity-compose shell app
Singularity app.sif:~/Documents/Dropbox/Code/singularity/singularity-compose-example> 
```

## exec

You can easily execute a command to a running instance:

```bash
$ singularity-compose exec app ls /
bin
boot
code
dev
environment
etc
home
lib
lib64
media
mnt
opt
proc
root
run
sbin
singularity
srv
sys
tmp
usr
var
```

## run

If a container has a `%runscript` section (or a Docker entrypoint/cmd that 
was converted to one), you can run that script easily:

```bash
$ singularity-compose run app
```

If your container didn't have any kind of runscript, the startscript
will be used instead.


## down

You can bring one or more instances down (meaning, stopping them) by doing:

```bash
$ singularity-compose down
Stopping (instance:nginx)
Stopping (instance:db)
Stopping (instance:app)
```

To stop a custom set, just specify them:


```bash
$ singularity-compose down nginx
```

It is also possible to specify a timeout (as for singularity instance stop)
in order to kill instances after the specified number of seconds:

```bash
singularity-compose down -t 100
```

## logs

You can of course view logs for all instances, or just specific named ones:

```bash
$ singularity-compose logs --tail 10
nginx ERR
nginx: [emerg] host not found in upstream "uwsgi" in /etc/nginx/conf.d/default.conf:22
2019/06/18 15:41:35 [emerg] 15#15: host not found in upstream "uwsgi" in /etc/nginx/conf.d/default.conf:22
nginx: [emerg] host not found in upstream "uwsgi" in /etc/nginx/conf.d/default.conf:22
2019/06/18 16:04:42 [emerg] 15#15: host not found in upstream "uwsgi" in /etc/nginx/conf.d/default.conf:22
nginx: [emerg] host not found in upstream "uwsgi" in /etc/nginx/conf.d/default.conf:22
2019/06/18 16:50:03 [emerg] 15#15: host not found in upstream "uwsgi" in /etc/nginx/conf.d/default.conf:22
nginx: [emerg] host not found in upstream "uwsgi" in /etc/nginx/conf.d/default.conf:22
2019/06/18 16:51:32 [emerg] 15#15: host not found in upstream "uwsgi" in /etc/nginx/conf.d/default.conf:22
nginx: [emerg] host not found in upstream "uwsgi" in /etc/nginx/conf.d/default.conf:22
```

## config

You can load and validate the configuration file (singularity-compose.yml) and
print it to the screen as follows:

```bash
$ singularity-compose config
{
    "version": "1.0",
    "instances": {
        "nginx": {
            "build": {
                "context": "./nginx",
                "recipe": "Singularity.nginx"
            },
            "volumes": [
                "./nginx.conf:/etc/nginx/conf.d/default.conf:ro",
                "./uwsgi_params.par:/etc/nginx/uwsgi_params.par:ro",
                ".:/code",
                "./static:/var/www/static",
                "./images:/var/www/images"
            ],
            "volumes_from": [
                "app"
            ],
            "ports": [
                "80"
            ]
        },
        "db": {
            "image": "docker://postgres:9.4",
            "volumes": [
                "db-data:/var/lib/postgresql/data"
            ]
        },
        "app": {
            "build": {
                "context": "./app"
            },
            "volumes": [
                ".:/code",
                "./static:/var/www/static",
                "./images:/var/www/images"
            ],
            "ports": [
                "5000:80"
            ],
            "depends_on": [
                "nginx"
            ]
        }
    }
}
```

# Global arguments

The following arguments are supported for all commands available.

## debug

Set logging verbosity to debug.

```bash
$ singularity-compose --debug version
```

This is equivalent to passing `--log-level=DEBUG` to the CLI.

```bash
$ singularity-compose --log-level='DEBUG' version
```

## log_level

Change logging verbosity. Accepted values are: `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`

```bash
$ singularity-compose --log-level='INFO' version
```

## file

Specify the location of a Compose configuration file

Default value: `singularity-compose.yml`

Aliases `--file`, `-f`. 

You can supply multiple `-f` configuration files. When you supply multiple files, `singularity-compose`
 combines them into a single configuration. It builds the configuration in the order you supply the
files. Subsequent files override and add to their predecessors.

For example consider this command:

```bash
$ singularity-compose -f singularity-compose.yml -f singularity-compose.dev.yml up
```

The `singularity-compose.yml` file might specify a `webapp` instance:

```yaml
instances:
  webapp:
    image: webapp.sif
    start:
      args: "start-daemon"
    port:
      - "80:80"
    volumes:
      - /mnt/shared_drive/folder:/webapp/data
```

if the `singularity-compose.dev.yml` also specifies this same service, any matching fields override 
the previous files.

```yaml
instances:
  webapp:
    start:
      args: "start-daemon -debug"
    volumes:
      - /home/user/folder:/webapp/data
```

The result of the examples above would be translated in runtime to:

```yaml
instances:
  webapp:
    image: webapp.sif
    start:
      args: "start-daemon -debug"
    port:
      - "80:80"
    volumes:
      - /home/user/folder:/webapp/data
```

## project-name

Specify project name.

Default value: `$PWD`

Aliases `--project-name`, `-p`. 

```bash
$ singularity-compose --project-name 'my_cool_project' up
```

## project-directory

Specify project working directory

Default value: compose file location


```bash
$ singularity-compose --project-directory /home/user/myfolder up
```

[home](/README.md#singularity-compose)
![logo](img/scompose.png)

# Singularity Compose <small>docs</small>

> Singularity Compose for simple container orchestration

- Definition of services
- Interaction with instances
- Logging and Networking


<style>
section.cover .cover-main > p:last-child a:last-child {
    background-color: #ffffff;
    color: black !important;
}

section.cover .cover-main>p:last-child a {
    border: 1px solid #ffffff !important;
    color: white !important;
}

.cover {
    background: linear-gradient(to left bottom, hsl(352, 88%, 48%) 0%,hsl(331, 100%, 19%) 100%) !important;;
    color: white;
}

.cover-main span {
    color: whitesmoke !important;
}
</style>

[GitHub](https://github.com/singularityhub/singularity-compose)
[Get Started](#singularity-compose)
# Singularity Compose

This is a simple orchestration library for Singularity containers, akin to
Docker Compose. It is under development, and working for basic examples.

## Who is this intended for?

Singularity compose is intended to run a small number of container instances
on your host. It is *not* a complicated orchestration tool like Kubernetes,
but rather a controlled way to represent and manage a set of container instances,
or services.

## When do I need sudo?

Singularity compose uses Singularity on the backend, so anything that would require sudo (root)
permissions for Singularity is also required for Singularity compose. This includes most
networking commands (e.g., asking to allocate ports) and builds from recipe files. 
However, if you are using Singularity v3.3 or higher, you can take advantage of 
[fakeroot](https://sylabs.io/guides/3.3/user-guide/fakeroot.html) to try and get around this.
The snippet below shows how to add fakeroot as an option under a build section:

```yaml
    build:
      context: ./nginx
      recipe: Singularity.nginx
      options:
       - fakeroot
```

And a complete example is provided [here](https://github.com/singularityhub/singularity-compose-examples/blob/master/rstudio-simple/singularity-compose.yml).


## Getting Started

### Dependencies

Singularity Compose *must* use a version of [Singularity](https://sylabs.io/guides/latest/user-guide/) 
3.2.1 or greater. It's recommended to use the latest (3.3.0 release at the time of this writing) otherwise there was
a bug with some versions of 3.2.1. Singularity 2.x absolutely will not work.
Python 3 is also required, as Python 2 is at end of life.

### singularity-compose.yml

For a singularity-compose project, it's expected to have a `singularity-compose.yml`
in the present working directory. You can look at a simple example here, here is a 
version 1.0 spec (before we added networking and exec options):

```yaml
version: "1.0"
instances:
  app:
    build:
      context: ./app
    volumes:
      - ./nginx.conf:/etc/nginx/conf.d/default.conf
      - ./app:/code
      - ./static:/var/www/static
      - ./images:/var/www/images
    ports:
      - 80:80
```

and here is a version 2.0 spec that shows adding networking and exec options:

```yaml
```

If you are familiar with [docker-compose](https://docs.docker.com/compose/) 
the file should look very familiar. A key difference is that instead of 
"services" we have "instances." And you guessed correctly - each 
section there corresponds to a 
[Singularity instance](https://sylabs.io/guides/3.2/user-guide/running_services.html)
that will be created. In this guide, we will walk through each of the sections
in detail.

### Instance folders

Generally, each section in the yaml file corresponds with a container instance to be run, 
and each container instance is matched to a folder in the present working directory.
For example, if I give instruction to build an `nginx` instance from
a `nginx/Singularity.nginx` file, I should have the
following in my singularity-compose:

```
  nginx:
    build:
      context: ./nginx
      recipe: Singularity.nginx
```

The above says that I want to build a container and corresponding instance
named `nginx`, and use the recipe `Singularity.nginx` in the context
folder `./nginx` in the present working directory. This gives me the following
directory structure:

```bash
singularity-compose-example
├── nginx
...
│   ├── Singularity.nginx
│   └── uwsgi_params.par
└── singularity-compose.yml

```

Notice how I also have other dependency files for the nginx container
in that folder.  While the context for starting containers with Singularity
compose is the directory location of the `singularity-compose.yml`,
the build context for this container is inside the nginx folder.
We will talk about the [build command](commands.md) soon. First,
as another option, you can just define a container to pull,
and it will be pulled to the same folder that is created if it doesn't exist.

```
  nginx:
    image: docker://nginx
```

This will pull a container `nginx.sif` into a `nginx` context folder:

```bash
├── nginx                    (- created if it doesn't exist
│   └── nginx.sif            (- named according to the instance
└── singularity-compose.yml
```

It's less likely that you will be able to pull a container that is ready to
go, as typically you will want to customize the 
[startscript](https://sylabs.io/guides/3.2/user-guide/definition_files.html#startscript) 
for the instance. Now that we understand the basic organization, let's
bring up some instances.

## Quick Start

For this quick start, we are going to use the 
[singularity-compose-simple](https://www.github.com/singularityhub/singularity-compose-simple) 
example. Singularity has a networking issue that currently doesn't allow communication
between multiple containers (due to iptables and firewall issues) so for now the most we
can do is show you one container. First, install singularity-compose from pip:

```bash
$ pip install singularity-compose
```

Then, clone the repository:

```bash
$ git clone https://www.github.com/singularityhub/singularity-compose-simple
```

cd inside, and you'll see a `singularity-compose.yml` like we talked about.

```bash
$ cd singularity-compose-simple
$ ls
app  images  LICENSE  nginx.conf  README.md  singularity-compose.yml  static
```

Let's take a look at the `singularity-compose.yml`

```yaml
version: "1.0"
instances:
  app:
    build:
      context: ./app
    volumes:
      - ./nginx.conf:/etc/nginx/conf.d/default.conf
      - ./app/nginx/uwsgi_params.par:/etc/nginx/uwsgi_params.par
      - ./app/nginx/cache:/var/cache/nginx
      - ./app/nginx/run:/var/run
      - ./app:/code
      - ./static:/var/www/static
      - ./images:/var/www/images
    ports:
      - 80:80
...
```

It defines a single service, `app`, which has both a Django application and
a nginx server with the [nginx-upload module](https://www.nginx.com/resources/wiki/modules/upload/) enabled.
It tells us right away that the folder `app` is the context folder, and inside we can
see dependency files for nginx and django.

```bash
$ ls app/
manage.py  nginx  requirements.txt  run_uwsgi.sh  Singularity  upload...
```

What we don't see is a container. We need to build that from the Singularity recipe.
Let's do that.

```bash
$ singularity-compose build
```

Will generate an `app.sif` in the folder:

```bash
$ ls app/
app.sif manage.py  nginx  requirements.txt  run_uwsgi.sh  Singularity  upload...
```

And now we can bring up our instance!

```bash
$ singularity-compose up
```

Verify it's running:

```bash
$ singularity-compose ps
INSTANCES  NAME PID     IMAGE
1           app	20023	app.sif
```

And then look at logs, shell inside, or execute a command.

```bash
$ singularity-compose logs app
$ singularity-compose logs app --tail 30
$ singularity-compose shell app
$ singularity-compose exec app uname -a
```

When you open your browser to [http://127.0.0.1](http://127.0.0.1)
you should see the upload interface. 

![img/upload.png](img/upload.png)

If you drop a file in the box (or click
to select) we will use the nginx-upload module to send it directly to the
server. Cool!

![img/content.png](img/content.png)

This is just a simple Django application, the database is sqlite3, and it's
now appeared in the app folder:

```bash
$ ls app/
app.sif  db.sqlite3  manage.py  nginx  requirements.txt  run_uwsgi.sh  Singularity  upload  uwsgi.ini
```

The images that you upload are stored in `images` at the root:

```bash
$ ls images/
2018-02-20-172617.jpg  40-acos.png  _upload 
```

And static files are in `static`.

```bash
$ ls static/
admin  css  js
```

Finally, the volumes that we specified in the `singularity-compose.yml`
tell us exactly where nginx and the application need write. The present
working directory (where the database is written) is bound to the
container at `/code`, and nginx dependencies are bound to locations
in `/etc/nginx`. Notice how the local static and images folder are bound
to locations in the container where we normally wouldn't have write.

```bash
    volumes:
      - ./nginx.conf:/etc/nginx/conf.d/default.conf
      - ./app/nginx/uwsgi_params.par:/etc/nginx/uwsgi_params.par
      - ./app/nginx/cache:/var/cache/nginx
      - ./app/nginx/run:/var/run
      - ./app:/code
      - ./static:/var/www/static
      - ./images:/var/www/images
```

This is likely a prime different between Singularity and Docker compose - Docker doesn't need
binds for write, but rather to reduce isolation. When you develop an application,
a lot of your debug will come down to figuring out where the services need to write
log and similar files, which you might not have been aware of when using Docker.

Continue below to read about networking, and see these commands in detail.

## Networking

When you bring the container up, you'll see generation of an `etc.hosts` file,
and if you guessed it, this is indeed bound to `/etc/hosts` in the container.
Let's take a look:

```bash
10.22.0.2	app
127.0.0.1	localhost

# The following lines are desirable for IPv6 capable hosts
::1     ip6-localhost ip6-loopback
fe00::0 ip6-localnet
ff00::0 ip6-mcastprefix
ff02::1 ip6-allnodes
ff02::2 ip6-allrouters
```

This file will give each container that you create (in our case, just one)
a name on its local network. Singularity by default creates a bridge for
instance containers, which you can conceptually think of as a router,
This means that, if I were to reference the hostname "app" in a second container,
it would resolve to `10.22.0.2`. Singularity compose does this by generating
these addresses before creating the instances, and then assigning them to it.
If you would like to see the full commands that are generated, run the up
with `--debug` (binds and full paths have been removed to make this easier to read).

```bash
$ singularity instance start \
    --bind /home/vanessa/Documents/Dropbox/Code/singularity/singularity-compose-simple/etc.hosts:/etc/hosts \
    --net --network-args "portmap=80:80/tcp" --network-args "IP=10.22.0.2" \
    --hostname app \
    --writable-tmpfs app.sif app
```

Control and customization of these instances is probably the coolest (and not widely
used) feature of Singularity. You can create your own network configurations,
and customie the arguments to the command. Read [here](https://sylabs.io/guides/3.2/user-guide/running_services.html) for more detalis.

## Commands

Read more about the commands shown above [here](commands.md#commands). For the
Python API, see [here](/singularity-compose/api/).

## Specification

The [specification](spec/) describes in more detail the sections of the singularity-compose.yml.
For example, in the quick start above, we have a post command for the app instance
that creates a series of folders on the host. 

## Examples

See our [Singularity Compose Examples](https://www.github.com/singularityhub/singularity-compose-examples) repository for
finding or contributing examples for using scompose.
# Singularity Compose Version 2.0

The second version of the singularity compose spec has added options for
starting an instance, starting, exec'ing and running commands after that. Note that spec v2.0 is
supported for Singularity Compose versions 0.1.0 and later.

## Overview

Here is a simple example of a 2.0 spec. The full example is provided at
[singularity-compose-examples](https://github.com/singularityhub/singularity-compose-examples/tree/master/v2.0/ping).

```yaml
version: "2.0"
instances:
  alp1:
    build:
      context: ./
      options:
        - fakeroot
    ports:
      - "1025:1025"
    start:
      options:
        - fakeroot
    exec:
      options: 
        - "env-file=myvars.env"
      command: printenv SUPERHERO
  alp2:
    build:
      context: ./
      options:
        - fakeroot
    ports:
      - "1026:1026"
    start:
      options:
        - fakeroot
    depends_on:
      - alp1
```

You'll notice that this differs from v1.0 because we've added groups for start, exec,
and run options. Start options and arguments are provided to the instance at start,
while exec options and arguments are exec'd to the instance after creation (only
if an exec argument exists). In the example above, we want to generate
two instances, each with an alpine base, and use fakeroot so that sudo is not required.


## Networking

Singularity Compose 0.1.0 later (version 2.0 of the spec here) supports using
fakeroot for start arguments, as shown above. However, you likely will need to
add the user of interest as follows:

```bash
$ sudo singularity config fakeroot --add $USER
```

The above command would add _your_ username. You'll also need to update the fakeroot
network configuration at `/usr/local/etc/singularity/network/40_fakeroot.conflist`
(or replace with the prefix where you installed Singularity):

```json
{
    "cniVersion": "0.4.0",
    "name": "fakeroot",
    "plugins": [
        {
            "type": "bridge",
            "bridge": "sbr0",
            "isGateway": true,
            "ipMasq": true,
            "ipam": {
                "type": "host-local",
                "subnet": "10.22.0.0/16",
                "routes": [
                    { "dst": "0.0.0.0/0" }
                ]
            }
        },
        {
            "type": "firewall"
        },
        {
            "type": "portmap",
            "capabilities": {"portMappings": true},
            "snat": true
        }
    ]
}
``` 

If you don't want to make these changes, then you won't be able to use fakeroot
as a network (start) option (you might still be able to use it as a build option).

## Network Group

### Allocate IP Address

By default `singularity-compose` will allocate an IP address for every instance in 
the listed yaml file. Binding an IP address to a process requires `sudo` so in certain
scenarios in which access to a privileged user isn't an option, you might want to tell
`singularity-compose` not to allocate an IP address, that way you can run instances 
without `sudo`. 

The example below will run a container that exposes the port `5432` to the host. 

```yaml
  instance1:
    ...
    network:      
      allocate_ip: true | false
    ports:
      - 5432:5432
```

**Observation:** In recent versions of the Singularity CLI, there is the need for tweaking the 
`/etc/singularity/singularity.conf` to allow `fakeroot` to bind to ports otherwise 
an error will be thrown at container execution similar to this:

```
INFO:    Converting SIF file to temporary sandbox...
ERROR:   Network fakeroot is not permitted for unprivileged users.
INFO:    Cleaning up image...
```

To allow fakeroot to bind ports without sudo you need to execute this:

```
echo "allow net networks = bridge, fakeroot" >> /etc/singularity/singularity.conf
```

### Enable/Disable Network

By default `singularity-compose` will always append `--net` to command to be executed in which it
will prompt for either having `--network=none` or `--fakeroot` added.

Depending on your environment's configuration and isolation requirements you may want to be able
to instruct `singularity-compose` not to append `--net` or any network-related params to the command
sent to singularity CLI.

The example below will disable the network features:

```yaml
  instance1:
    ...
    network:      
      enable: true | false
```


## Start Group

Startscript options generally include those for networking, and any other flags
that you want to provide. The previous "command" option is deprecated, and moved to be under the "start"
group as "args," since we are technically providing arguments to the start script. 
As an example in the configuration above, we are starting with options
for `--fakeroot`:

```yaml
  alp1:
    ...
    start:
      options:
        - fakeroot
```

You could also add "args" here within the start group to provide arguments to the start script.


## Environment

While Singularity compose doesn't currently have support for an environment 
section, it's easy to add custom environments by way of binding an environment
file to the instance! For example:

```yaml
    build:
      context: ./app
    volumes:
      - ./env_file.sh:/.singularity.d/env/env_file.sh
```

Any file that is found in the `.singularity.d/env` folder will be sourced.
For example, you could define or export variables like so:

```bash
#!/bin/bash
MYNAME=dinosaur
export MYNAME
```

Make sure to export variables.


## Exec Group

The exec group, if defined, will run an exec command after an instance is started.
There must be a command defined for this to run. For example,
if you want to provide a one off environment variable to an exec (run after start)
then you can do the following:

```yaml
    exec:
      options: 
        - "env-file=myvars.env"
      command: printenv MYNAME
```

As with the above, make sure to export variables.
We are only printing there to show that it works.

## Run Group

Let's say that you want to start the container, exec a command, and then run the container.
This is possible with the run group, and your container must define a runscript.
If you just want the run to happen (without options or arguments) you can do:

```yaml
  alp1:
    run: []
```

And if you want args or options, you can again add them:

```yaml
  alp1:
    run:
      args: "arg1 arg2 arg3"
      options: 
        - "env-file=myvars.env"
```

The run and exec sections are separate to allow you to run both, or either without
the other.

## Instance

An instance currently must be instantiated from a container built 
from a Singularity recipe in a named folder (the example above) 
alongside the singularity-compose.yml:

```
  nginx:
    build:
      context: ./nginx
      recipe: Singularity.nginx
    volumes:
      - ./nginx.conf:/etc/nginx/conf.d/default.conf:ro
      - ./uwsgi_params.par:/etc/nginx/uwsgi_params.par:ro
    volumes_from:
      - uwsgi
    ports:
      - "80"
```

or from a container unique resource identifier (uri) that can be pulled
to a local directory with the same name as the section.

```
  nginx:
    image: docker://vanessa/sregistry_nginx
    volumes:
      - ./nginx.conf:/etc/nginx/conf.d/default.conf:ro
      - ./uwsgi_params.par:/etc/nginx/uwsgi_params.par:ro
    volumes_from:
      - uwsgi
    ports:
      - "80"
```

We build or pull a local container for reproducibility. The first approach,
building from a local file, is recommended as you have full control over
the environment and startscript, and there are likely few containers out in the
wild (and ready to be pulled) that have the correct entry and start scripts
for your application. In the case that you *do* want to pull
a container (and use the existing startscript or entrypoint) you can do that
as follows:

```
  nginx:
    image: docker://vanessa/sregistry_nginx
...
```

Customization of an image (e.g., labels, help, post) is out of
scope for singularity-compose, and you must build from a recipe instead.
The fields for instances are discussed below:


### Fields

|Name| Description |
|----|-------------|
|name|The name of the instance will be "nginx" unless the user provides a "name" field (not defined above).|
|build| a section to define how and where to build the base container from.|
|build.context| the folder with the Singularity file (and other relevant files). Must exist.
|build.recipe| the Singularity recipe in the build context folder. It defaults to `Singularity`|
|build.options| a list of one or more options (single strings for boolean, or key value pairs for arguments) to provide to build.  This is where you could provide fakeroot.|
|start| a section to define start (networking) arguments and options |
|start.options| a list of one or more options for starting the instance |
|start.args| arguments to provide to the startscript when starting the instance |
|run| a section to define run arguments and options (container must have runscript) |
|run.options| a list of one or more options for running the instance after start |
|run.args| arguments to provide when running the instance |
|exec| a section to define an exec directly after instance start (requires a command) |
|exec.options| a list of one or more options for exec'ing to the instance |
|exec.command| the command and arguments to provide the instance exec |
|image| is looked for after a build entry. It can be a unique resource identifier, or container binary. |
|volumes| one or more files or files to bind to the instance when it's started.|
|volumes_from| shared volumes that are defined for other instances|
|ports| currently not sure how I'm going to handle this!|
|post| a section of post commands and arguments, run after instance creation |
|post.commands| a list of commands to run (directly or a script) on the host |
# Singularity Compose Version 1.0

## Overview

The main format of the file is yaml. We must define a version and one or more
instances under "instances." Here is a full example for reference.

```yaml
version: "1.0"
instances:

  nginx:
    build:
      context: ./nginx
      recipe: Singularity.nginx
    volumes:
      - ./nginx.conf:/etc/nginx/conf.d/default.conf:ro
      - ./uwsgi_params.par:/etc/nginx/uwsgi_params.par:ro
    volumes_from:
      - app
    ports:
      - "80"

  db:
    image: postgres:9.4
    volumes:
      - db-data:/var/lib/postgresql/data

  app:
    build:
      context: ./app
    volumes:
      - .:/code
      - ./static:/var/www/static
      - ./images:/var/www/images
    ports:
      - "5000:80"
    depends_on:
      - nginx
```

Each of nginx, uwsgi, and db are instances to be built as containers, and run
as instances. 

## Networking

Since Singularity does not (currently) have control over custom networking,
all instance ports are mapped to the host (localhost) and we don't have any
configuration settings to control this (how to handle ports?)

## Startscript arguments

It is possible to use the `command` option to pass arguments to an instance's
startscript.

The following example shows how to pass the arguments `arg0 arg1 arg2` to the
startscript of instance `app`,

```yaml
  app:
    build:
      context: ./app
    command: "arg0 arg1 arg2"
```

## Environment

While Singularity compose doesn't currently have support for an environment 
section, it's easy to add custom environments by way of binding an environment
file to the instance! For example:

```yaml
  app:
    build:
      context: ./app
    volumes:
      - ./env_file.sh:/.singularity.d/env/env_file.sh
```

Any file that is found in the `.singularity.d/env` folder will be sourced.
For example, you could define or export variables like so:

```bash
#!/bin/bash
MYNAME=dinosaur
export MYNAME
```

Make sure to export variables.

## Instance

An instance currently must be instantiated from a container built 
from a Singularity recipe in a named folder (the example above) 
alongside the singularity-compose.yml:

```
  nginx:
    build:
      context: ./nginx
      recipe: Singularity.nginx
    volumes:
      - ./nginx.conf:/etc/nginx/conf.d/default.conf:ro
      - ./uwsgi_params.par:/etc/nginx/uwsgi_params.par:ro
    volumes_from:
      - uwsgi
    ports:
      - "80"
```

or from a container unique resource identifier (uri) that can be pulled
to a local directory with the same name as the section.

```
  nginx:
    image: docker://vanessa/sregistry_nginx
    volumes:
      - ./nginx.conf:/etc/nginx/conf.d/default.conf:ro
      - ./uwsgi_params.par:/etc/nginx/uwsgi_params.par:ro
    volumes_from:
      - uwsgi
    ports:
      - "80"
```

We build or pull a local container for reproducibility. The first approach,
building from a local file, is recommended as you have full control over
the environment and startscript, and there are likely few containers out in the
wild (and ready to be pulled) that have the correct entry and start scripts
for your application. In the case that you *do* want to pull
a container (and use the existing startscript or entrypoint) you can do that
as follows:

```
  nginx:
    image: docker://vanessa/sregistry_nginx
...
```

Customization of an image (e.g., labels, help, post) is out of
scope for singularity-compose, and you must build from a recipe instead.
The fields for instances are discussed below:


### Fields

|Name| Description |
|----|-------------|
|name|The name of the instance will be "nginx" unless the user provides a "name"
field (not defined above).|
|build| a section to define how and where to build the base container from.|
|build.context| the folder with the Singularity file (and other relevant files). Must exist.
|build.recipe| the Singularity recipe in the build context folder. It defaults to `Singularity`|
|build.options| a list of one or more options (single strings for boolean, or key value pairs for arguments) to provide to build.  This is where you could provide fakeroot.|
|image| is looked for after a build entry. It can be a unique resource identifier, or container binary. |
|volumes| one or more files or files to bind to the instance when it's started.|
|volumes_from| shared volumes that are defined for other instances|
|ports| currently not sure how I'm going to handle this!|
|post| a section of post commands and arguments, run after instance creation |
|post.commands| a list of commands to run (directly or a script) on the host |
# Specifications

This folder contains specifications for each version of a singularity-compose.yml file.

  - [2.0](spec/spec-2.0.md) is version 2 of the specification with added network, start, exec, and run options. See [the file here](https://github.com/singularityhub/singularity-compose/tree/master/docs/spec/spec-2.0.md). You should use version 0.1.0 or later for this spec.
 - [1.0](spec/spec-1.0.md) is version 1 (current) of the under development specification. See [the file here](https://github.com/singularityhub/singularity-compose/tree/master/docs/spec/spec-1.0.md). You should use a version less than 0.1.0 for this version.

../CHANGELOG.mdWelcome to Singularity Compose API's documentation!
===================================================

Singularity Compose is an orchestration tool for Singularity instances.

Contents:

.. toctree::
   :maxdepth: 3

   source/scompose.rst
   changelog.md


Indices and tables
------------------

* :ref:`modindex`
scompose.utils package
======================

Module contents
---------------

.. automodule:: scompose.utils
    :members:
    :undoc-members:
    :show-inheritance:
scompose.templates package
==========================

Module contents
---------------

.. automodule:: scompose.templates
    :members:
    :undoc-members:
    :show-inheritance:
scompose.client package
=======================

Submodules
----------

scompose.client.build module
----------------------------

.. automodule:: scompose.client.build
    :members:
    :undoc-members:
    :show-inheritance:

scompose.client.config module
-----------------------------

.. automodule:: scompose.client.config
    :members:
    :undoc-members:
    :show-inheritance:

scompose.client.create module
-----------------------------

.. automodule:: scompose.client.create
    :members:
    :undoc-members:
    :show-inheritance:

scompose.client.down module
---------------------------

.. automodule:: scompose.client.down
    :members:
    :undoc-members:
    :show-inheritance:

scompose.client.exec module
---------------------------

.. automodule:: scompose.client.exec
    :members:
    :undoc-members:
    :show-inheritance:

scompose.client.logs module
---------------------------

.. automodule:: scompose.client.logs
    :members:
    :undoc-members:
    :show-inheritance:

scompose.client.ps module
-------------------------

.. automodule:: scompose.client.ps
    :members:
    :undoc-members:
    :show-inheritance:

scompose.client.restart module
------------------------------

.. automodule:: scompose.client.restart
    :members:
    :undoc-members:
    :show-inheritance:

scompose.client.shell module
----------------------------

.. automodule:: scompose.client.shell
    :members:
    :undoc-members:
    :show-inheritance:

scompose.client.up module
-------------------------

.. automodule:: scompose.client.up
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: scompose.client
    :members:
    :undoc-members:
    :show-inheritance:
scompose.project package
========================

Submodules
----------

scompose.project.instance module
--------------------------------

.. automodule:: scompose.project.instance
    :members:
    :undoc-members:
    :show-inheritance:

scompose.project.project module
-------------------------------

.. automodule:: scompose.project.project
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: scompose.project
    :members:
    :undoc-members:
    :show-inheritance:
scompose.tests package
======================

Subpackages
-----------

.. toctree::

    scompose.tests.testdata

Submodules
----------

scompose.tests.test\_client module
----------------------------------

.. automodule:: scompose.tests.test_client
    :members:
    :undoc-members:
    :show-inheritance:

scompose.tests.test\_utils module
---------------------------------

.. automodule:: scompose.tests.test_utils
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: scompose.tests
    :members:
    :undoc-members:
    :show-inheritance:
scompose.tests.testdata package
===============================

Module contents
---------------

.. automodule:: scompose.tests.testdata
    :members:
    :undoc-members:
    :show-inheritance:
singularity-compose
===================

.. toctree::
   :maxdepth: 4

   scompose
scompose.logger package
=======================

Submodules
----------

scompose.logger.message module
------------------------------

.. automodule:: scompose.logger.message
    :members:
    :undoc-members:
    :show-inheritance:

scompose.logger.progress module
-------------------------------

.. automodule:: scompose.logger.progress
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: scompose.logger
    :members:
    :undoc-members:
    :show-inheritance:
scompose package
================

Subpackages
-----------

.. toctree::

    scompose.client
    scompose.logger
    scompose.project
    scompose.templates
    scompose.tests
    scompose.utils

Submodules
----------

scompose.version module
-----------------------

.. automodule:: scompose.version
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: scompose
    :members:
    :undoc-members:
    :show-inheritance:
