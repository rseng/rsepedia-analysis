# pending

# 0.4.1

Small release to archive for JOSS acceptance.

- Move `cloud_sql_proxy` installation before code copy (https://github.com/google/caliban/pull/87)

# 0.4.0

The biggest feature in this new release is native support for logging to an
MLFlow tracking server using the [UV
Metrics](http://github.com/google/uv-metrics) project.
(https://github.com/google/caliban/pull/35) This feature is in alpha; expect
documentation soon.

### More features

- minor bugfixes for GKE (https://github.com/google/caliban/pull/85)
- additional tests for gke.{types, util} (https://github.com/google/caliban/pull/84)
- re-order custom apt packages before pip requirements (https://github.com/google/caliban/pull/82)
- modify base image to our more general cloudbuild naming scheme (https://github.com/google/caliban/pull/80)
- updated `google-auth` dependency version to `1.19.0` (https://github.com/google/caliban/pull/79)
- add clearer contribution info (https://github.com/google/caliban/pull/76)
- Update uv-metrics tutorial (https://github.com/google/caliban/pull/74, https://github.com/google/caliban/pull/72)
- add support for running an embedded cloudsql_proxy (https://github.com/google/caliban/pull/60)
- bugfix for #65: do not add resource maxima when quota is < 1 (#67)
- Updated accelerator regions (and globally availabe AI Platform regions to
  match the current state here):
  https://cloud.google.com/ai-platform/training/docs/regions

# 0.3.0

- @ramasesh Added a fix that prevented `pip` git dependencies from working in
  `caliban shell` mode (https://github.com/google/caliban/pull/55) This adds a
  small update to the base image, so be sure to run

```
docker pull gcr.io/blueshift-playground/blueshift:cpu
docker pull gcr.io/blueshift-playground/blueshift:gpu
```

to get access to this fix.

- Thanks to @eschnett, `--docker_run-args` can now deal with arbitrary
  whitespace in the list of arguments, instead of single spaces only.
  (https://github.com/google/caliban/pull/46)

- Caliban now authenticates AI Platform job submissions using the authentication
  provided by `gcloud auth login`, rather than requiring a service account key.
  This significantly simplifies the setup required for a first time user.

- `caliban cloud` now checks if the image exists remotely before issuing a
  `docker push` command on the newly built image
  (https://github.com/google/caliban/pull/36)

- Big internal refactor to make it easier to work on code, increase test
  coverage, add new backends (https://github.com/google/caliban/pull/32)

- add `schema` validation for `.calibanconfig.json`. This makes it much easier
  to add configuration knobs: https://github.com/google/caliban/pull/37

- Custom base image support (https://github.com/google/caliban/pull/39), thanks
  to https://github.com/google/caliban/pull/20 from @sagravat.
  `.calibanconfig.json` now supports a `"base_image"` key. For the value, can
  supply:
  - a Docker base image of your own
  - a dict of the form `{"cpu": "base_image", "gpu": "base_image"}` with both
    entries optional, of course.

  Two more cool features.

  First, if you use a format string, like `"my_image-{}:latest"`, the format
  block `{}` will be filled in with either `cpu` or `gpu`, depending on the mode
  Caliban is using.

  Second, we now have native support for [Google's Deep Learning
  VMs](https://cloud.google.com/ai-platform/deep-learning-vm/docs/introduction)
  as base images. The actual VM containers [live
  here](https://console.cloud.google.com/gcr/images/deeplearning-platform-release/GLOBAL).
  If you provide any of the following strings, Caliban will expand them out to
  the actual base image location:

```
dlvm:pytorch-cpu
dlvm:pytorch-cpu-1.0
dlvm:pytorch-cpu-1.1
dlvm:pytorch-cpu-1.2
dlvm:pytorch-cpu-1.3
dlvm:pytorch-cpu-1.4
dlvm:pytorch-gpu
dlvm:pytorch-gpu-1.0
dlvm:pytorch-gpu-1.1
dlvm:pytorch-gpu-1.2
dlvm:pytorch-gpu-1.3
dlvm:pytorch-gpu-1.4
dlvm:tf-cpu
dlvm:tf-cpu-1.0
dlvm:tf-cpu-1.13
dlvm:tf-cpu-1.14
dlvm:tf-cpu-1.15
dlvm:tf-gpu
dlvm:tf-gpu-1.0
dlvm:tf-gpu-1.13
dlvm:tf-gpu-1.14
dlvm:tf-gpu-1.15
dlvm:tf2-cpu
dlvm:tf2-cpu-2.0
dlvm:tf2-cpu-2.1
dlvm:tf2-cpu-2.2
dlvm:tf2-gpu
dlvm:tf2-gpu-2.0
dlvm:tf2-gpu-2.1
dlvm:tf2-gpu-2.2
```

Format strings work here as well! So, `"dlvm:pytorch-{}-1.4"` is a totally valid
base image.

# 0.2.6

- Prepared for a variety of base images by setting up a cloud build matrix:
  https://github.com/google/caliban/pull/25
- Added better documentation for `gcloud auth configure-docker`
  https://github.com/google/caliban/pull/26
- Added `close()` to `TqdmFile`, preventing an error when piping `stdout`:
  https://github.com/google/caliban/pull/30
- `tqdm` progress bars and other interactive outputs now display correctly in
  `caliban run` outputs. `stdout` flushes properly! Before these changes,
  `stderr` would appear before any `stdout`, making it difficult to store the
  logs in a text file. Now, by default, python processes launched by `caliban
  run` won't buffer. https://github.com/google/caliban/pull/31

![2020-06-26 09 48 50](https://user-images.githubusercontent.com/69635/85877300-2a3e7300-b794-11ea-9792-4cf3ae5e4263.gif)

# 0.2.5

- fixes the python binary that caliban notebook points to (now that we use
  conda)
- adds DEBIAN_FRONTEND=noninteractive to the apt-get command, so that packages
  like texlive won't freeze and wait for you to specify a timezone.

This makes it easy to add, for example, npm and latex support to your caliban
notebook invocations.

# 0.2.4

- fixes a bug with `parse_region` not handling a lack of default.
- converts the build to Github Actions.
- Rolls Caliban back to requiring only python 3.6 support.
- Removes some unused imports from a few files.

# 0.2.3

- Added fix for an issue where large user IDs would crash Docker during the
  build phase. https://github.com/google/caliban/pull/8

# 0.2.2

- Fix for bug with requirements.txt files.

# 0.2.1

- Added support for Conda dependencies
  (https://github.com/google/caliban/pull/5). If you include `environment.yml`
  in your project's folder Caliban will attempt to install all dependencies.
- Caliban now uses a slimmer base image for GPU mode.
- Base images for CPU and GPU modes now use Conda to manage the container's
  virtual environment instead of `virtualenv`.

# 0.2.0

- Caliban now caches the service account key and ADC file; you should see faster
  builds, BUT you might run into trouble if you try to run multiple Caliban
  commands in the same directory in different processes, due to a race condition
  with a temp file dropped in the directory. If you see a failure, try again.
- Private Cloud Source Repositories are now supported as requirements in the
  `requirements.txt` and `setup.py` of projects executed using Caliban.
- `caliban notebook` now installs `jupyter` instead of `jupyterlab`. `caliban
  notebook --lab` of course still uses `jupyterlab`.
- Caliban now works with Python 3.5.
- The default release channel for gke clusters in caliban is now 'regular', as
  node autoprovisioning now works with preemptible instances in the regular
  channel.
- If you provide `--cloud_key` Caliban will now properly use the supplied cloud
  key to submit jobs to AI platform. Previously, Caliban would rely on the
  system's auth method, which made it impossible to point to a different
  project if your current service account key wasn't the owner.
- Changed gke job submission to accept min cpu and min mem arguments instead
  of an explicit machine-type. This allows gke to efficiently schedule jobs
  and prevents issues where jobs can be oversubscribed on compute nodes.
- added support for `.calibanconfig.json`. You can now add this file with an
  `"apt_packages"` entry to specify aptitude packages to install inside the
  container. The value under this key can be either a list, or a dictionary with
  `"gpu"` and `"cpu"'` keys. For example, any of the following are valid:

```
# This is a list by itself. Comments are fine, by the way.
{
     "apt_packages": ["libsm6", "libxext6", "libxrender-dev"]
}
```

This works too:

```
# You can also include a dictionary with different deps
# for gpu and cpu modes. It's fine to leave either of these blank,
# or not include it.
{
    "apt_packages": {
        "gpu": ["libsm6", "libxext6", "libxrender-dev"],
        "cpu": ["some_other_package"]
    }
}
```

# 0.1.15

- `caliban notebook` now attempts to search for the first free port instead of
  failing due to an already-occupied port.

- `pip` is now called with `--no-cache-dir` inside the container; this should
  shrink container sizes with no impact on performance.

- All commands have a new `--no-cache` option; when supplied, Docker will skip
  using its build cache. This is helpful to use if you want to, say, force new
  dependencies to get installed without bumping their versions explicitly.

# 0.1.14

- JSON experiment configuration files can now handle arguments which are varied
  together, by supplying a compound key, of the form e.g. `[arg1,arg2]`.

- better error messages print when a docker command fails.

- Caliban can now handle pushing containers to "domain scoped projects":
  https://cloud.google.com/container-registry/docs/overview#domain-scoped_projects
  The colon in the project name separating domain and project ID is handled
  properly.

# 0.1.13

- 'caliban run' and 'caliban shell' now take an --image_id argument; if
  provided, these commands will skip their 'docker build' phase and use the
  image ID directly.

- AI Platform labels now swap periods for underscores (thanks to vinay@!); this
  means that floating point numbers will no longer have pre- and post- decimal
  components concatenated.

- A new `expansion` script will expand experiment configs into all of the
  individual experiments they'd generate. This command can accept `stdin`, just
  like the `--experiment_config` argument. Options include `--pprint` and
  `--print_flags`. The output of this script can be piped directly into `caliban
  cloud --experiment_config stdin`.

- `caliban shell` will now default to bash if you're using a shell that's not
  `bash` or `zsh` (fish shell, for example) instead of erroring out.

- `caliban shell` has a new `--shell` argument that you can use to override the
  container's default shell.

# 0.1.12

- consolidated gke tpu/gpu spec parsing with cloud types
- modified all commands to accept as the module argument paths to arbitrary
  shell scripts. Any argument of the format "trainer.train" will execute using
  "python -m trainer.train", just as before. If instead you pass a python script
  as a file, like "trainer/train.py", caliban will execute this file inside the
  container using "python trainer/train.py". Any other argument, if it exists in
  the local directory, will be executed as a bash script.

  This allows users to run commands like "caliban cloud my_script.sh" and have
  it all work.

- "caliban run" now supports --experiment_config and --dry_run. These work just
  like they do for "caliban cloud"; the experiment config will expand out and
  execute N jobs on your local machine.
- moved some methods from cluster/cluster.py to gke/utils.py
- added unit tests for some gke/utils.py methods
- Support for ADC credentials! if application_default_credentials.json is
  present on the user's machine, they now get copied into the container.
- if ADC credentials are NOT present but a service account key is we write a
  placeholder. this is required to get ctpu working inside containers.

# 0.1.11

- added tpu driver specification for gke jobs
- added query for getting available tpu drivers for cluster/project

# 0.1.10

- set host_ipc=True for cluster jobs

# 0.1.9

- moved cluster constants to separate file
- moved cluster gpu validation to separate file
- added test for gpu limits validation

# 0.1.8

- TPU and GPU spec now accept validate_count arg to disable count validation.

# 0.1.7

- Fixed a bug where the label for the job name wasn't getting properly
  sanitized - this meant that if you provided an upper-cased job name job
  submission would fail.
- Fixed a bug that prevented parsing experiment config values that were floats.
- experiment config parsing now performs the full expansion at CLI-parse-time
  and validates every expanded config.

# 0.1.6

- `--docker_run_args` allows you to pass a string of arguments directly through
  to `docker run`. This command works for `caliban run`, `caliban notebook` and
  `caliban shell`.

- `docker.py` reorganized, now takes explicit `JobMode` instances throughout.

- new `--cloud_key` argument for all commands. If specified this lets you
  override the `GOOGLE_APPLICATION_CREDENTIALS` value that's usually inspected
  for the service account key.

- fixed a bug in `caliban.util.TempCopy` where a `None`-valued path would fail. This affected environments where `GOOGLE_APPLICATION_CREDENTIALS` wasn't set.

# 0.1.5

- `--experiment_config` can now take experiment configs via stdin (pipes, yay!);
  specify `--experiment_config stdin`, or any-cased version of that, and the
  script will wait to accept your input.

  As an example, this command pipes in a config and also passes `--dry_run` to
  show the series of jobs that WILL be submitted when the `--dry_run` flag is
  removed:

  ```
  cat experiment.json | caliban cloud -e gpu --experiment_config stdin --dry_run trainer.train
  ```

  You could pipe the output of a nontrivial python script that generates a JSON
  list of dicts.

- `--image_tag` argument to `caliban cloud`; if you supply this it will bypass
  the Docker build and push steps and use this image directly. This is useful if
  you want to submit a job quickly without going through a no-op build and push,
  OR if you want to broadcast an experiment to some existing container.

- if you supply a `--tpu_spec` and DON'T supply an explicit `--gpu_spec`,
  caliban will default to CPU mode. `--gpu_spec` and `--nogpu` are still
  incompatible. You can use a GPU and TPU spec together without problems.

- added support for `--tpu_spec` in `caliban cloud` mode. This validates in a
  similar style to `--gpu_spec`; any invalid combination of count, region and
  TPU type will fail.

  (Unlike `--gpu_spec` this mode IS compatible with `--nogpu`. In fact, many
  demos seem to use the non-GPU version of tensorflow here.)

- added a `caliban build` mode that simply builds the image and returns the
  image ID. Useful for checking if your image can build at all with the current
  settings.

# 0.1.4

- the CLI will now error if you pass any caliban keyword arguments AFTER the
  python module name, but before `--`. In previous versions, if you did something like

```bash
caliban cloud trainer.train --nogpu
```

  That final `--nogpu` would get passed on directly to your script, vs getting
  absorbed by Caliban.

- If you pass `--nogpu` mode and have a setup.py file, caliban will
  automatically pass `--extras cpu` and attempt to install an extras dependency
  called `cpu`. Same goes for `gpu` if you DON'T pass `--nogpu`. So, if you have
  a setup.py file, the default pip installation will now look like one of these:

```bash
pip install .[gpu]
pip install .[cpu]
```

Instead of

```
pip install .
```

# 0.1.3.1

- Minor bugfix; I was calling "len" on an iterator, not a list.

# 0.1.3

This version:

- Moves from batch submission to submitting 1 request at a time, so Cloud can handle our rate limiting
- adds support for lists of experiment configs
- adds terminal highlighting
- adds a progress bar for Cloud submissions
- cleans up the terminal output for each job to show the interesting bits, not the entire training spec.
- moves the job index to the end of the jobId, so searching is easier

If you like you can set `-v 1` to see the full spec output.

# 0.1.2

- `caliban.cloud.types` has lots of enums and types that make it easier to code
  well against AI Platform.
- `caliban cloud` has more validations, now.
- `caliban cloud` now supports:
  - `--gpu_spec`, which you can use to configure the GPU count and type for your
    job.
  - `--machine_type` allows you to specify the machine type for all jobs that
    run.
  - `--experiment_config` lets you submit batches of jobs at once.
  - `--force` skips all validations and forces a submission with the specified
    config.
  - `--dry_run` will generate logs showing what WOULD happen if you submit a
    batch of jobs.
- The `caliban cloud --stream_logs` argument is now gone; the command prints,
  so this is easy enough to run without special help, and the argument made
  batch job submission difficult.
- all local Docker commands now run `--ipc host`, which gives `docker run`
  access to all of the host's memory.
- Base images now contain `gsutil`. `docker.py` automatically configures
  `gcloud` and `gsutil` with shared credentials when they're present on the
  machine running caliban commands.
# Committers

These are the folks who can +1 a pull request and approve it for merge.

## Active

| Name            | Handle                                               |
|-----------------|------------------------------------------------------|
| Sam Ritchie     | [@isnotinvain](https://github.com/sritchie)          |
| Ambrose Slone   | [@ajslone](https://github.com/ajslone)               |
| Guy Gur-Ari     | [@guygurari](https://github.com/guygurari)           |
| Vinay Ramasesh  | [@ramasesh](https://github.com/ramasesh)             |


## Emeritus
# Caliban

[![Build status](https://github.com/google/caliban/workflows/build/badge.svg?branch=master)](https://github.com/google/caliban/actions?query=workflow%3Abuild+branch%3Amaster)
[![Codecov branch](https://img.shields.io/codecov/c/github/google/caliban/master.svg?maxAge=3600)](https://codecov.io/github/google/caliban)
[![JOSS](https://joss.theoj.org/papers/c33c8b464103b2fb3b641878722bf8f3/status.svg)](https://joss.theoj.org/papers/c33c8b464103b2fb3b641878722bf8f3)
[![readthedocs](https://img.shields.io/readthedocs/caliban?maxAge=3600)](https://caliban.readthedocs.io/en/latest/?badge=latest)
[![caliban version](https://img.shields.io/pypi/v/caliban?maxAge=3600)](https://pypi.org/project/caliban)

Caliban is a tool that helps researchers launch and track their numerical
experiments in an isolated, reproducible computing environment. It was developed
by machine learning researchers and engineers, and makes it easy to go from a
simple prototype running on a workstation to thousands of experimental jobs
running on Cloud.

With Caliban, you can:

- Develop your experimental code locally and test it inside an isolated (Docker)
  environment
- Easily sweep over experimental parameters
- Submit your experiments as Cloud jobs, where they will run in the same
  isolated environment
- Control and keep track of jobs

## Quickstart

[Install Docker](#docker), make sure it's running, then install Caliban (you'll need [Python >= 3.6](#python-36)):

```bash
pip install caliban
```

Train a simple deep learning model on your local machine:

```bash
git clone https://github.com/google/caliban.git && cd caliban/tutorials/basic
caliban run --nogpu mnist.py
```

Sweep over learning rates to find the best one (flags are specified in JSON format):

```bash
echo '{"learning_rate": [0.01, 0.001, 0.0001]}' | caliban run --experiment_config stdin --nogpu mnist.py
```

**Next**:

- See how to submit the experiment to Cloud and use other Caliban features in ["Getting Started with Caliban"](#getting-started-with-caliban)
- See [Installation](#installation-and-prerequisites) for detailed installation instructions
- Read the [Command Overview](#command-overview) for info on Caliban commands.

Full documentation for Caliban lives at [Read The Docs](https://caliban.readthedocs.io/en/latest).

### Dramatic Interlude

<p>
<img style="float: right;" align="right" src="https://upload.wikimedia.org/wikipedia/commons/a/ad/Stephano%2C_Trinculo_and_Caliban_dancing_from_The_Tempest_by_Johann_Heinrich_Ramberg.jpg" width="350">

> “Be not afeard; the isle is full of noises, \
> Sounds, and sweet airs, that give delight and hurt not. \
> Sometimes a thousand twangling instruments \
> Will hum about mine ears; and sometime voices, \
> That, if I then had waked after long sleep, \
> Will make me sleep again: and then, in dreaming, \
> The clouds methought would open, and show riches \
> Ready to drop upon me; that, when I waked, \
> I cried to dream again.”
>
> -- <cite>Shakespeare, The Tempest</cite>
</p>

## Installation and Prerequisites

Caliban's prequisites are [Docker](#docker) and [Python >= 3.6](#python-36).

Make sure your Python is up to date:

```bash
$ python --version
Python 3.6.9 # should be >=3.6.0
```

If not, visit ["Installing Python 3.6"](#python-36) before proceeding.

Next, install Caliban via [pip](https://pypi.org/project/caliban/):

```bash
pip install -U caliban
```

check if your installation worked by navigating to an empty folder and running
`caliban --help`. You should see the usage dialogue:

```bash
$ caliban --help
usage: caliban [-h] [--helpfull] [--version]
               {shell,notebook,build,run,cloud,cluster,status,stop,resubmit}
               ...
```

### Docker

Caliban executes your code inside a "container", managed by
[Docker](https://hub.docker.com/editions/community/docker-ce-desktop-mac). To get Docker:

- On MacOS, follow the installation instructions at [Docker
  Desktop](https://hub.docker.com/editions/community/docker-ce-desktop-mac) and
  start the newly-installed Docker Desktop application.
- On Linux, visit the [Docker installation
  instructions](https://docs.docker.com/engine/install/ubuntu/#installation-methods).
  (It's important that you configure [sudo-less
  Docker](https://caliban.readthedocs.io/en/latest/getting_started/prerequisites.html#docker)
  and start Docker running on your machine.)

Make sure Docker is correctly installed, configured and running by executing the
following command:

```bash
docker run hello-world
```

You should see output that looks like this:

```text
...
Hello from Docker!
This message shows that your installation appears to be working correctly.
...
```

### Python 3.6

Make sure your Python version is up to date:

```bash
$ python --version
Python 3.6.9 # should be >=3.6.0
```

If you need to upgrade:

- On MacOS, install the latest Python version from
  [python.org](https://www.python.org/downloads/mac-osx) ([direct
  link](https://www.python.org/ftp/python/3.8.3/python-3.8.3-macosx10.9.pkg)).
- On Linux, run `sudo apt-get update && sudo apt-get install python3.7`.

### Cloud Submission and GPUs

Caliban's [Read the Docs](https://caliban.readthedocs.io/) documentation has
instructions on:

- [Installing the `nvidia-docker2`
  runtime](https://caliban.readthedocs.io/en/latest/getting_started/prerequisites.html#docker-and-cuda),
  so you can use Caliban to run jobs that use your Linux machine's GPU.
- [Setting up a Google Cloud
  account](https://caliban.readthedocs.io/en/latest/getting_started/cloud.html)
  so you can submit your code to Google's [Cloud AI
  Platform](https://cloud.google.com/ai-platform) with `caliban cloud`.

## Getting Started with Caliban

In this section we will use Caliban to train an image classification network
(implemented in
[TensorFlow](https://www.tensorflow.org/tutorials/quickstart/beginner)). We
will:

- Train a neural network on the local machine
- Increase the model's accuracy by changing the [learning
  rate](https://medium.com/octavian-ai/which-optimizer-and-learning-rate-should-i-use-for-deep-learning-5acb418f9b2)
  with a command-line flag
- Sweep across a range of learning rates with Caliban's [experiment
  broadcasting](https://caliban.readthedocs.io/en/latest/explore/experiment_broadcasting.html)
  feature
- Train the model in the Cloud on Google's [AI
  Platform](https://cloud.google.com/ai-platform)
- Develop code interactively using `caliban shell` in the exact same
  environment.

### Preparing your Project

Create an empty directory and use `curl` to download a [python
script](https://github.com/google/caliban/blob/master/tutorials/basic/mnist.py#L16)
that trains a basic neural network.

```
mkdir demo && cd demo
curl --output mnist.py https://raw.githubusercontent.com/google/caliban/master/tutorials/basic/mnist.py
```

Create a file called `requirements.txt` to declare `tensorflow-cpu` as a dependency:

```bash
echo "tensorflow-cpu" > requirements.txt
```

Caliban will automatically make any entry in `requirements.txt` available when
you run your code. See ["Declaring
Requirements"](https://caliban.readthedocs.io/en/latest/explore/declaring_requirements.html)
for more information.

### Training the Network

Run this command to train your first ML model:

```bash
caliban run --nogpu mnist.py
```

You should see a stream of output ending in this:

```text
Training model with learning rate=0.1 for 3 epochs.
Epoch 1/3
1875/1875 - 3s - loss: 2.0989 - accuracy: 0.2506
Epoch 2/3
1875/1875 - 3s - loss: 1.9222 - accuracy: 0.2273
Epoch 3/3
1875/1875 - 3s - loss: 2.0777 - accuracy: 0.1938
Model performance:
313/313 - 0s - loss: 2.0973 - accuracy: 0.1858
```

Your model was able to recognize digits from the
[MNIST](https://en.wikipedia.org/wiki/MNIST_database) dataset with 18.58%
accuracy. Can we do better?

### Improving the Model

The default learning rate is `0.1`. Run the code again with a smaller learning
rate by passing a command-line flag, separated from your original command by
`--`:

```bash
$ caliban run --nogpu mnist.py -- --learning_rate 0.01

<<elided>>

Training model with learning rate=0.01 for 3 epochs.
Epoch 1/3
1875/1875 - 4s - loss: 0.2676 - accuracy: 0.9221
Epoch 2/3
1875/1875 - 4s - loss: 0.1863 - accuracy: 0.9506
Epoch 3/3
1875/1875 - 4s - loss: 0.1567 - accuracy: 0.9585
Model performance:
313/313 - 0s - loss: 0.1410 - accuracy: 0.9642
```

96% accuracy! Much better! Can we do better still?

### Experiment Broadcasting

Caliban's [experiment
broadcasting](https://caliban.readthedocs.io/en/latest/explore/experiment_broadcasting.html)
feature will allow us to run many jobs with different sets of arguments.

Create a file called `experiment.json` with a
[JSON](https://www.json.org/json-en.html) dictionary of the format
`{"flag_name": ["list", "of", "values"]}`:

```bash
echo '{"learning_rate": [0.01, 0.001, 0.0001]}' > experiment.json
```

Pass the config with `--experiment_config` and run again:

```bash
caliban run --experiment_config experiment.json --nogpu mnist.py
```

You should see accuracies of roughly `0.9493`, `0.9723` and `0.9537`. Looks like
`0.001` is a nice choice.

### Submitting to Cloud AI Platform

Now it's time to submit the job to [Cloud AI
Platform](https://cloud.google.com/ai-platform).

(**NOTE**: This section requires a Google Cloud account. You can create a free
account with $300 of credit to get started. Follow Caliban's ["Getting Started
with Google
Cloud"](https://caliban.readthedocs.io/en/latest/getting_started/cloud.html)
documentation, then come back here to proceed.)

Submit the job to AI Platform by changing the word `run` to `cloud`:

```bash
caliban cloud --nogpu mnist.py -- --learning_rate 0.01
```

You should see output like this:

```bash
I0615 19:57:43.354172 4563361216 core.py:161] Job 1 - jobId: caliban_totoro_1, image: gcr.io/research-3141/974a776e6037:latest
I0615 19:57:43.354712 4563361216 core.py:161] Job 1 - Accelerator: {'count': 0, 'type': 'ACCELERATOR_TYPE_UNSPECIFIED'}, machine: 'n1-highcpu-32', region: 'us-central1'
I0615 19:57:43.355082 4563361216 core.py:161] Job 1 - Experiment arguments: ['--learning_rate', '0.01']
I0615 19:57:43.355440 4563361216 core.py:161] Job 1 - labels: {'gpu_enabled': 'false', 'tpu_enabled': 'false', 'job_name': 'caliban_totoro', 'learning_rate': '0_01'}

I0615 19:57:43.356621 4563361216 core.py:324] Submitting request!
I0615 19:57:45.078382 4563361216 core.py:97] Request for job 'caliban_totoro_20200615_195743_1' succeeded!
I0615 19:57:45.078989 4563361216 core.py:98] Job URL: https://console.cloud.google.com/ai-platform/jobs/caliban_totoro_20200615_195743_1?projectId=totoro-project
I0615 19:57:45.079524 4563361216 core.py:100] Streaming log CLI command: $ gcloud ai-platform jobs stream-logs caliban_totoro_20200615_195743_1
Submitting caliban_totoro_1: 100%|####################################################################################################################################################################################| 1/1 [00:02<00:00,  2.65s/requests]
I0615 19:57:45.405600 4563361216 core.py:673]
I0615 19:57:45.405819 4563361216 core.py:676] Visit https://console.cloud.google.com/ai-platform/jobs/?projectId=research-3141 to see the status of all jobs.
I0615 19:57:45.405959 4563361216 core.py:677]
```

This output means that Caliban has:

- built a Docker container with all of your code
- Pushed that container up to Google Cloud's [Container
  Registry](https://cloud.google.com/container-registry)
- Submitted the job to [AI Platform](https://cloud.google.com/ai-platform).

You can now visit the link in the output that looks like:
https://console.cloud.google.com/ai-platform/jobs/caliban_totoro_20200615_195743_1?projectId=totoro-project
to see all of your job's logs.

#### Why do I need Cloud?

With Google Cloud, you can use on-demand
[GPUs](https://caliban.readthedocs.io/en/latest/cloud/gpu_specs.html) and
[TPUs](https://caliban.readthedocs.io/en/latest/cloud/ai_platform_tpu.html) and
train models on large datasets at very high speeds. You can also customize the
[machine
type](https://caliban.readthedocs.io/en/latest/cloud/gpu_specs.html#custom-machine-types)
that AI Platform uses to run your job. You might need high memory or more CPU,
for example.

See Caliban's ["Customizing Machines and
GPUs"](https://caliban.readthedocs.io/en/latest/cloud/gpu_specs.html#) for more
information.

### Interactive Development with `caliban shell`

[`caliban
shell`](https://caliban.readthedocs.io/en/latest/cli/caliban_shell.html) lets
you develop code interactively inside of the exact same environment that your
code will have available, locally during `caliban run` or in the Cloud with
`caliban cloud`.

Run the following command to activate the shell:

```bash
caliban shell --nogpu
```

You should see Caliban's terminal:

```
I0611 12:33:17.551121 4500135360 docker.py:911] Running command: docker run --ipc host -w /usr/app -u 735994:89939 -v /Users/totoro/code/example:/usr/app -it --entrypoint /bin/bash -v /Users/totoro:/home/totoro ab8a7d7db868
   _________    __    ________  ___    _   __  __  __
  / ____/   |  / /   /  _/ __ )/   |  / | / /  \ \ \ \
 / /   / /| | / /    / // __  / /| | /  |/ /    \ \ \ \
/ /___/ ___ |/ /____/ // /_/ / ___ |/ /|  /     / / / /
\____/_/  |_/_____/___/_____/_/  |_/_/ |_/     /_/ /_/

You are running caliban shell as user with ID 735994 and group 89939,
which should map to the ID and group for your user on the Docker host. Great!

[totoro@6a9b28990757 /usr/app]$
```

You're now living in an isolated [Docker
container](https://www.docker.com/resources/what-container) with your
`tensorflow-cpu` dependency available (and any others [you've
declared](https://caliban.readthedocs.io/en/latest/explore/declaring_requirements.html)).

Run the `python` command and check that `tensorflow` is installed:

```bash
$ python
Python 3.6.9 (default, Nov  7 2019, 10:44:02)
[GCC 8.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import tensorflow as tf
>>> tf.__version__
'2.2.0'
```

Your home directory and the folder where you ran the command are both mounted
into this isolated environment, so any changes you make to either of those
directories will be reflected immediately.

Any code you add to the current folder and edit on your computer will be
available in this special Caliban shell. Run the example from before like this:

```
python mnist.py --learning_rate 0.01
```

If your code runs in `caliban shell`, you can be almost certain that your code
will execute in a Cloud environment, with potentially many GPUs attached and
much larger machines available.

### What next?

Read the [Overview](#overview) for more information on Caliban's subcommands,
then head over to [Caliban's documentation
site](https://caliban.readthedocs.io/en/latest/) and check out the links on the
sidebar.

If you find anything confusing, please feel free to [create an
issue](https://github.com/google/caliban/issues) on our [Github Issues
page](https://github.com/google/caliban/issues), and we'll get you sorted out.

## Command Overview

Caliban provides seven subcommands that you run inside some project directory on
your machine:

* [`caliban
  shell`](https://caliban.readthedocs.io/en/latest/cli/caliban_shell.html)
  generates a Docker image containing any dependencies you've declared in a
  `requirements.txt` and/or `setup.py` in the directory and opens an interactive
  shell in that directory. The `caliban shell` environment is ~identical to the
  environment that will be available to your code when you submit it to AI
  Platform; the difference is that your current directory is live-mounted into
  the container, so you can develop interactively.

* [`caliban
  notebook`](https://caliban.readthedocs.io/en/latest/cli/caliban_notebook.html)
  starts a Jupyter notebook or lab instance inside of a Docker image containing
  your dependencies; the guarantee about an environment identical to AI Platform
  applies here as well.

* [`caliban run`](https://caliban.readthedocs.io/en/latest/cli/caliban_run.html)
  packages your directory's code into the Docker image and executes it locally
  using `docker run`. If you have a GPU, the instance will attach to it by
  default - no need to install the CUDA toolkit. The Docker environment takes
  care of all that. This environment is truly identical to the AI Platform
  environment. The Docker image that runs locally is the same image that will
  run in AI Platform.

* [`caliban
  cloud`](https://caliban.readthedocs.io/en/latest/cli/caliban_cloud.html)
  allows you to [submit jobs to AI
  Platform](https://caliban.readthedocs.io/en/latest/getting_started/cloud.html)
  that will run inside the same Docker image you used with `caliban run`. You
  can submit hundreds of jobs at once. Any machine type, GPU count, and GPU type
  combination you specify will be validated client side, so you'll see an
  immediate error with suggestions, rather than having to debug by submitting
  jobs over and over.

* [`caliban
  build`](https://caliban.readthedocs.io/en/latest/cli/caliban_build.html) builds
  the Docker image used in `caliban cloud` and `caliban run` without actually
  running the container or submitting any code.

* [`caliban
  cluster`](https://caliban.readthedocs.io/en/latest/cli/caliban_cluster.html)
  creates GKE clusters and submits jobs to GKE clusters.

* [`caliban
  status`](https://caliban.readthedocs.io/en/latest/cli/caliban_status.html)
  displays information about all jobs submitted by Caliban, and makes it easy to
  interact with large groups of experiments. Use `caliban status` when you need
  to cancel pending jobs, or re-build a container and resubmit a batch of
  experiments after fixing a bug.

## Disclaimer

This is a research project, not an official Google product. Expect bugs and
sharp edges. Please help by trying out Caliban, [reporting
bugs](https://github.com/google/caliban/issues), and letting us know what you
think!

## Get Involved + Get Support

Pull requests and bug reports are always welcome! Check out our [Contributor's
Guide](CONTRIBUTING.md) for information on how to get started contributing to
Caliban.

The TL;DR; is:

- send us a pull request,
- iterate on the feedback + discussion, and
- get a +1 from a [Committer](COMMITTERS.md)

in order to get your PR accepted.

Issues should be reported on the [GitHub issue
tracker](https://github.com/google/caliban/issues).

If you want to discuss an idea for a new feature or ask us a question,
discussion occurs primarily in the body of [Github
Issues](https://github.com/google/caliban/issues), though the project is growing
large enough that we may start a Gitter channel soon.

The current list of active committers (who can +1 a pull request) can be found
here: [COMMITTERS.md](COMMITTERS.md)

A list of contributors to the project can be found at the project's
[Contributors](https://github.com/google/caliban/graphs/contributors) page.

## Citing Caliban

If Caliban helps you in your research, please consider citing Caliban's
associated academic paper:

```
@article{Ritchie2020,
  doi = {10.21105/joss.02403},
  url = {https://doi.org/10.21105/joss.02403},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {53},
  pages = {2403},
  author = {Sam Ritchie and Ambrose Slone and Vinay Ramasesh},
  title = {Caliban: Docker-based job manager for reproducible workflows},
  journal = {Journal of Open Source Software}
}
```

## License

Copyright 2020 Google LLC.

Licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).
# How to Contribute

So you want to add some code to Caliban. Excellent!

Pull requests and bug reports are always welcome! Check out our [Contributor's
Guide](CONTRIBUTING.md) for information on how to get started contributing to
Caliban.

The TL;DR; is:

- send us a pull request,
- iterate on the feedback + discussion, and
- get a +1 from a [Committer](COMMITTERS.md)

in order to get your PR accepted.

Issues should be reported on the [GitHub issue
tracker](https://github.com/google/caliban/issues).

If you want to discuss an idea for a new feature or ask us a question,
discussion occurs primarily in the body of [Github
Issues](https://github.com/google/caliban/issues), though the project is growing
large enough that we may start a Gitter channel soon.

The current list of active committers (who can +1 a pull request) can be found
here: [COMMITTERS.md](COMMITTERS.md)

A list of contributors to the project can be found at the project's
[Contributors](https://github.com/google/caliban/graphs/contributors) page.

## Contributor License Agreement

Contributions to this project must be accompanied by a Contributor License
Agreement. You (or your employer) retain the copyright to your contribution;
this simply gives us permission to use and redistribute your contributions as
part of the project. Head over to <https://cla.developers.google.com/> to see
your current agreements on file or to sign a new one.

You generally only need to submit a CLA once, so if you've already submitted one
(even if it was for a different project), you probably don't need to do it
again.

## Developing in Caliban

We use [pre-commit](https://pre-commit.com/) to manage a series of git
pre-commit hooks for the project; for example, each time you commit code, the
hooks will make sure that your python is formatted properly. If your code isn't,
the hook will format it, so when you try to commit the second time you'll get
past the hook.

All hooks are defined in `.pre-commit-config.yaml`. To install these hooks,
install `pre-commit` if you don't yet have it. I prefer using
[pipx](https://github.com/pipxproject/pipx) so that `pre-commit` stays globally
available.

```bash
pipx install pre-commit
```

Then install the hooks with this command:

```bash
pre-commit install
```

Now they'll run on every commit. If you want to run them manually, you can run either of these commands:

```bash
pre-commit run --all-files

# or this, if you've previously run `make build`:
make lint
```

## Documentation

We use Sphinx to generate docs. If you want to live-preview your changes to the
documentation as you are editing, you can use
[sphinx-reload](https://pypi.org/project/sphinx-reload/). To get this working:

```bash
pipx install sphinx-reload
```

Then, inside the caliban folder:

```bash
make build
sphinx-reload docs
```

If all goes well, `sphinx-reload` will tell you it is serving the documentation
on a port, which you can listen into from your browser.

## Publishing Caliban

- First, run `make build` to get your virtual environment set up.
- Make sure that you're on the master branch!
- add a new tag, with `git tag 0.2.3` or the equivalent
- run `make release` to push the latest code and tags to all relevant
  repositories.
---
title: 'Caliban: Docker-based job manager for reproducible workflows'
tags:
  - python
  - docker
  - machine learning
  - reproducibility
authors:
  - name: Sam Ritchie
    orcid: 0000-0002-0545-6360
    affiliation: 1
  - name: Ambrose Slone
    affiliation: 1
  - name: Vinay Ramasesh
    orcid: 0000-0003-0625-3327
    affiliation: 1
affiliations:
 - name: Google, United States of America
   index: 1
date: 22 June 2020
bibliography: paper.bib
---

# Summary

Caliban is a command line tool that helps researchers launch and track their
numerical experiments in an isolated, reproducible computing environment. It was
developed by machine learning researchers and engineers, and makes it easy to go
from a simple prototype running on a workstation to thousands of experimental
jobs running in a Cloud environment.

# Motivation

Modern machine learning research typically requires a researcher to execute code
in multiple computing environments. To investigate some property of a machine
learning model, the researcher has to write a script that can accept a path to
some dataset, train a model using some set of configurable parameters, and
generate measurements or a serialized model for later analysis.

Writing and debugging model training code is fastest on a local workstation.
Running the script to generate measurements almost always takes place on some
much more powerful machine, typically in a Cloud environment. The [imagenet
dataset](https://www.tensorflow.org/datasets/catalog/imagenet2012) [@deng2009imagenet]
is 144 GiB, for example, far too large to process on a stock laptop.

Moving between these environments almost always requires a nontrivial amount of
effort on the part of the researcher, unrelated to the original research
question. In Python, a common language for machine learning research, the
researcher installs their script's dependencies in a system-wide package
registry. The Cloud environment can have different dependencies available, or
different versions of the same dependency; different versions of Python itself;
or different software drivers, which elicit different behavior from the script
that seemed to work locally.

This environment mismatch introduces friction into the research process. The
only way to debug a script that succeeds locally and fails in a remote Cloud
environment is to attempt to interpret the often-cryptic error logs.

## Docker

One solution to this problem is Docker [@merkel2014docker]. A researcher can
package their code and dependencies inside of a Docker container and execute
this container on different platforms, each with different hardware options
available but with a consistent software environment.

Many Cloud services (Google's [Cloud AI
Platform](https://cloud.google.com/ai-platform), [Amazon
Sagemaker](https://aws.amazon.com/sagemaker/) etc) allow users to submit and
execute Docker containers that they've built locally.

Packaging research code inside of a Docker container has many benefits for
reproducibility [@cito2016]. But the process of building a Docker container is
difficult, error-prone and orthogonal to the skill set of a machine learning
researcher. The friction of debugging between local and Cloud environments is
solved, but only by accepting a not-insignificant baseline level of toil into
the local development experience.

Projects like [MLFlow](https://mlflow.org/docs/latest/projects.html)
[@zaharia2018accelerating] attempt to streamline the container creation process,
but still force the researcher to absorb much of Docker's mental model.

# Caliban and Reproducible Research

Caliban is a command line tool that solves this problem by providing execution
modes with opinionated, intuitive interfaces for each phase of machine learning
research - interactive development, local execution, cloud execution and data
analysis in a notebook environment.

With Caliban, the researcher executes all code using Caliban's various
subcommands. This process is, for the researcher, identical to the process of
executing code directly. To prepare a research environment, all they need to do
is specify required packages in a `requirements.txt` file and Caliban will
automatically make those dependencies available inside the container.

Behind the scenes, all software execution has moved inside of a Docker
container. This makes it transparent to move that execution from a local
environment to Cloud, enabling a research project to grow from a simple
prototype running on a workstation to thousands of experimental jobs running on
Cloud with little to no cognitive load.

This removal of friction allows a researcher the freedom to be creative in ways
that their psychology would resist, given the typical pain caused by moves
between environments.

In addition, Caliban makes it easy to launch multiple jobs with varying
command-line arguments with a single command, using experiment configuration
files.

# Impact

Before we introduced Caliban in our lab, researchers reported spending multiple
weeks learning how to run their experiments at scale. Caliban allows a new
researcher to reach this level of proficiency in under an hour. Multiple papers
currently in preparation contain research results generated from thousands of
experiments executed using Caliban. Given the limited tenure of a typical
machine learning internship and residency, the efficiency boost offered by
Caliban can materially change the scope of project that a researcher would be
willing to take on.

In addition, any research conducted with Caliban that the researcher open
sources is trivially executable by any interested party. The Docker environment
managed by Caliban guarantees that anyone with access to the shared source code
repository will be able to run experiments in an environment identical to the
environment used in the original research program.

# Caliban's Execution Environments

Caliban provides a suite of execution engines that can execute the containers
built by Caliban.

[`caliban
shell`](https://caliban.readthedocs.io/en/latest/cli/caliban_shell.html)
generates a Docker image containing any dependencies declared in a
`requirements.txt` and/or `setup.py` in a project's directory and opens an
interactive shell. Any update to code in the project's folder will be reflected
immediately inside the container environment. The `caliban shell` environment is
identical to the environment available during Cloud execution, up to access to
different hardware.

[`caliban
notebook`](https://caliban.readthedocs.io/en/latest/cli/caliban_notebook.html)
starts a Jupyter notebook or lab instance inside of a Docker image containing
the project's dependencies. As with `caliban shell`, the environment available
to the notebook is identical to the Cloud environment.

[`caliban run`](https://caliban.readthedocs.io/en/latest/cli/caliban_run.html)
packages a project's code into a Docker container and executes it locally using
`docker run`. If the local machine has access to a GPU, the instance will attach
to it by default, with no need to configure drivers.

[`caliban
cloud`](https://caliban.readthedocs.io/en/latest/cli/caliban_cloud.html) allows
a researcher to [submit a research script to Google's AI
Platform](https://caliban.readthedocs.io/en/latest/getting_started/cloud.html).
The code will run inside the same Docker container available with all other
subcommands. Researchers can submit hundreds of jobs at once. Any machine type,
GPU count, and GPU type combination specified will be validated client side,
instead of requiring a round trip to the server.

[`caliban
build`](https://caliban.readthedocs.io/en/latest/cli/caliban_build.html) builds
the Docker image used in `caliban cloud` and `caliban run` without actually
running the container or submitting any code.

[`caliban
cluster`](https://caliban.readthedocs.io/en/latest/cli/caliban_cluster.html)
allows a researcher to create and submit jobs to a Kubernetes cluster. This
environment is superficially similar to the Cloud environment offered by
Google's AI Platform and `caliban cloud`, but a different range of hardware and
pricing is available to the researcher.

[`caliban
status`](https://caliban.readthedocs.io/en/latest/cli/caliban_status.html)
displays information about all jobs submitted by Caliban, and allows a
researcher to cancel, inspect or resubmit large groups of experiments.

# Related Work

**Reproducible Environments** [Binder](https://mybinder.org/)
[@Forde2018ReproducingML] is a project that allows a researcher to open up
Jupyter notebooks hosted in a git repository in a software environment described
by the repository's `requirements.txt` file. The motivation is similar to
`caliban notebook`, with the added benefit of being able to interact with a
notebook in the browser, without any burden of configuring a local machine.
[repo2docker](https://github.com/jupyter/repo2docker) offers the same ability to
execute a repository of notebooks in its required environment. This environment
can be [customized and
configured](https://repo2docker.readthedocs.io/en/latest/config_files.html#dockerfile-advanced-environment)
in similar ways to containers build by Caliban.

**Scientific Cloud Computing** A related series of approaches to lowering the
friction of scientific computing on the Cloud expose primitives to the user that
can execute both locally, or in parallel on many Cloud machines. A "pure" Python
function is a function depends only on its inputs. If a scientific experiment is
built out of pure functions, it becomes possible to write a framework that can
execute that function in parallel on a large set of distinct inputs in a Cloud
environment. [Pywren](http://pywren.io/) [@DBLP:journals/corr/JonasVSR17] is a
project that makes this possible on [AWS
Lambda](https://aws.amazon.com/lambda/).
[Cloudknot](https://github.com/nrdg/cloudknot)
[@adam_richie-halford-proc-scipy-2018] extends this model to functions requiring
more computational resources. Cloudknot packages code into a Docker image,
providing the same potential level of configurability as Caliban.

# References
# Caliban Tutorials

This directory contains a number of tutorials that show off various aspects of
[Caliban](https://github.com/google/caliban).

The `basic` directory contains the code for the ["Getting Started with
Caliban"](https://github.com/google/caliban#getting-started-with-caliban)
tutorial on the main page of [Caliban's github
repository](https://github.com/google/caliban).

More coming soon!
# UV + MLFlow Tutorial [ALPHA!]

This directory contains a demo of a model training workflow that uses the
[uv-metrics](https://github.com/google/uv-metrics) library to persist metrics to
an [MLFlow](https://mlflow.org/) tracking server.

This is mostly here for testing and reference. Check back for a documentation
update once the API settles down.

## Prerequisites

Right now we are supporting logging metrics to a sql-based backing store only
in this tutorial, but we will update things to allow for local storage in the
future. For now you will need to have a google cloud sql instance configured
for this, and you will need an MLFlow server set up to serve results from
this instance.

To run this tutorial, you will need to edit the `.calibanconfig.json`
file in this directory to reflect your database settings so that the training
script can connect to the database and log metrics. The specific entries to
edit here are in the `mflow_config` entry in `.calibanconfig.json`:

```
{
  "apt_packages" : ["openssh-client", "curl"],
  "mlflow_config" : {"project": <your gcp project where your cloudsql db lives>,
                     "region": <the region where your database lives>,
                     "db": <the name of your mlflow database>,
                     "user": <connect as this database user>,
                     "password": <the database password for the above user>,
                     "artifact_root": <the location to store artifacts, typically a gs bucket>,
                     "debug" : false}
}
```

One note here is that currently artifact storage is not working completely, but
please specify this entry and we will update this tutorial once that is working properly.

Once you have set these parameters properly, you should be able to run the tutorial code.

## Sanity Check (optional)

A quick sanity check to test your database connection is to set the `debug` flag in
the `.calibanconfig.json` file to `true`, and then use Caliban to run the `hello_world.sh`
script. This script simply prints "hello, world", but by enabling the `debug` flag, we
can check the status of the database connection.

To run this test:

```
caliban run --nogpu hello_world.sh
```

If your database settings are configured properly, you should see output like the following:

```
Successfully built 5eb8dcef14ce
I0807 13:02:53.008464 139963939288896 tqdm.py:90] Restoring pure python logging
I0807 13:02:53.010536 139963939288896 run.py:74]
I0807 13:02:53.010816 139963939288896 run.py:75] Job 1 - Experiment args: []
I0807 13:02:53.010974 139963939288896 run.py:198] Running command: docker run --ipc host -e PYTHONUNBUFFERED=1 -e COLUMNS=211 -e LINES=19 5eb8dcef14ce ...
2020/08/07 20:02:53 current FDs rlimit set to 1048576, wanted limit is 8500. Nothing to do here.
2020/08/07 20:02:53 using credential file for authentication; path="/home/<username>/.config/gcloud/application_default_credentials.json"
2020/08/07 20:02:54 Listening on /tmp/cloudsql/<project>:<region>:<db>/.s.PGSQL.5432 for <project>:<region>:<db>
2020/08/07 20:02:54 Ready for new connections
INFO:root:/bin/bash hello_world.sh
hello, world
I0807 13:03:04.015075 139963939288896 run.py:111] Job 1 succeeded!
```

As long as you see `Ready for new connections`, then your configuration should be ok, and you
can disable the `debug` flag and continue with the rest of the tutorial.

## Running a Job

In the Caliban repository:

```
git checkout master && git pull
cd tutorials/uv-metrics
```

Run a single job:

```
caliban run --nogpu trainer.train
```

Name the experiment group and run 3:

```
caliban run --experiment_config experiment.json --xgroup mlflow_tutorial --nogpu trainer.train
```

## Check the MLFlow UI

You may need to refresh, but the UI should now show multiple experiments. You can view the
status and metrics for your jobs from the UI while your jobs are in progress, which is
useful for long-running jobs.
# Basic Tutorial

This directory contains the code for the ["Getting Started with
Caliban"](https://github.com/google/caliban#getting-started-with-caliban)
tutorial on the main page of [Caliban's github
repository](https://github.com/google/caliban).

Visit ["Getting Started with
Caliban"](https://github.com/google/caliban#getting-started-with-caliban) for
the full tutorial, and instructions on how to run the code in this folder.
Caliban
=======

Caliban is a tool for developing research workflow and notebooks in an isolated
Docker environment and submitting those isolated environments to Google Compute
Cloud.

For a short tutorial introduction to Caliban, see the `GitHub page
<https://github.com/google/caliban#caliban>`_.

Overview
--------

Caliban provides five subcommands that you run inside some directory on your
laptop or workstation:

* :doc:`/cli/caliban_shell` generates a Docker image containing any dependencies
  you've declared in a ``requirements.txt`` and/or ``setup.py`` in the directory
  and opens an interactive shell in that directory. The ``caliban shell``
  environment is ~identical to the environment that will be available to your
  code when you submit it to AI Platform; the difference is that your current
  directory is live-mounted into the container, so you can develop
  interactively.

* :doc:`/cli/caliban_notebook` starts a Jupyter notebook or lab instance inside
  of a docker image containing your dependencies; the guarantee about an
  environment identical to AI Platform applies here as well.

* :doc:`/cli/caliban_run` packages your directory's code into the Docker image
  and executes it locally using ``docker run``. If you have a workstation GPU,
  the instance will attach to it by default - no need to install the CUDA
  toolkit. The docker environment takes care of all that. This environment is
  truly identical to the AI Platform environment. The docker image that runs
  locally is the same image that will run in AI Platform.

* :doc:`/cli/caliban_cloud` allows you to submit jobs to AI Platform that will
  run inside the same docker image you used with ``caliban run``. You can submit
  hundreds of jobs at once. Any machine type, GPU count, and GPU type
  combination you specify will be validated client side, so you'll see an
  immediate error with suggestions, rather than having to debug by submitting
  jobs over and over.

* :doc:`/cli/caliban_build` builds the docker image used in ``caliban cloud``
  and ``caliban run`` without actually running the container or submitting any
  code.

* :doc:`/cli/caliban_cluster` creates GKE clusters and submits jobs to GKE
  clusters.

* :doc:`/cli/caliban_status` displays information about all jobs submitted by
  Caliban, and makes it easy to interact with large groups of experiments. Use
  :doc:`/cli/caliban_status` when you need to cancel pending jobs, or re-build a
  container and resubmit a batch of experiments after fixing a bug.

These all work from :doc:`your Macbook Pro <explore/mac>`. (Yes, you can build
and submit GPU jobs to Cloud from your Mac!)

The only requirement for the directory where you run these commands is that it
declare some set of dependencies in either a ``requirements.txt`` or
``setup.py`` file. See the :doc:`requirements docs
<explore/declaring_requirements>` for more detail.

The rest of this document contains detailed information and guides on Caliban's
various modes. If you want to get started in a more interactive way, head over
to `the Caliban tutorials
directory <https://github.com/google/caliban/blob/master/tutorials/README.md>`_.

Caliban's code lives on `Github <https://github.com/google/caliban>`_.

Using Caliban
-------------

If you want to practice using Caliban with a proper getting-started style guide,
head over to `Caliban's tutorials
<https://github.com/google/caliban/blob/master/tutorials/README.md>`_ (Coming
Soon!).

See the sidebar for information on the subcommands exposed by Caliban and a
whole series of tutorials and guides that you might find interesting as you work
with Caliban.

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   getting_started/prerequisites
   getting_started/getting_caliban
   getting_started/cloud

.. toctree::
   :maxdepth: 1
   :caption: Using Caliban

   cli/caliban_shell
   cli/caliban_notebook
   cli/caliban_run
   cli/caliban_cloud
   cli/caliban_build
   cli/caliban_cluster
   cli/caliban_status
   cli/caliban_stop
   cli/caliban_resubmit
   cli/expansion

.. toctree::
   :maxdepth: 1
   :caption: Exploring Further

   explore/why_caliban
   explore/base_image
   explore/custom_docker_run
   explore/declaring_requirements
   explore/experiment_groups
   explore/calibanconfig
   explore/custom_script_args
   explore/experiment_broadcasting
   explore/exp_stdin
   explore/script_vs_module
   explore/gcloud
   explore/mac

.. toctree::
   :maxdepth: 1
   :caption: Common Recipes

   recipes/flagfile
   recipes/single_gpu
   recipes/local_dir
   recipes/dockerignore

.. toctree::
   :maxdepth: 1
   :caption: Cloud-Specific Tutorials

   cloud/labels
   cloud/gpu_specs
   cloud/ai_platform_tpu
   cloud/rate_limit
   cloud/service_account
   cloud/adc
   cloud/bucket

.. toctree::
   :maxdepth: 1
   :caption: Caliban + GKE

   gke/concepts
   gke/prereq
   gke/cluster_management
   gke/job_submission
=======
Caliban reference documentation
===================================

Composable transformations of Python+NumPy programs: differentiate, vectorize,
JIT to GPU/TPU, and more.

For an introduction to Caliban, start at the `Caliban GitHub page
<https://github.com/google/caliban>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. argparse::
   :module: caliban.cli
   :func: caliban_parser
   :prog: caliban

.. argparse::
   :module: caliban.expansion
   :func: expansion_parser
   :prog: expansion


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
caliban notebook
^^^^^^^^^^^^^^^^

This command generates the same isolated environment as the other commands, but
instead of running your code or dropping you into a shell, runs a local instance
of Jupyter based in the folder where you execute the command.

``caliban notebook`` supports the following arguments:

.. code-block:: text

   usage: caliban notebook [-h] [--helpfull] [--nogpu] [--cloud_key CLOUD_KEY]
                           [--extras EXTRAS] [--docker_run_args DOCKER_RUN_ARGS]
                           [-p PORT] [-jv JUPYTER_VERSION] [--lab] [--bare]

   optional arguments:
     -h, --help            show this help message and exit
     --helpfull            show full help message and exit
     --nogpu               Disable GPU mode and force CPU-only.
     --cloud_key CLOUD_KEY
                           Path to GCloud service account key. (Defaults to
                           $GOOGLE_APPLICATION_CREDENTIALS.)
     --extras EXTRAS       setup.py dependency keys.
     --docker_run_args DOCKER_RUN_ARGS
                           String of args to add to Docker.
     -p PORT, --port PORT  Port to use for Jupyter, inside container and locally.
     -jv JUPYTER_VERSION, --jupyter_version JUPYTER_VERSION
                           Jupyterlab version to install via pip.
     --lab                 run 'jupyter lab', vs the default 'jupyter notebook'.
     --bare                Skip mounting the $HOME directory; run an isolated
                           Jupyter lab.

By default ``caliban notebook`` runs ``jupyter notebook`` inside the container. To
run Jupyterlab, pass the ``--lab`` flag:

.. code-block:: bash

   caliban notebook --lab

As with the other commands, the only python dependencies available in the
container will be dependencies that you declare explicitly in either:


* a ``requirements.txt`` file
* a ``setup.py`` file.

Your setup file can declare groups of dependencies using the setuptools
`extras_require
<https://setuptools.readthedocs.io/en/latest/setuptools.html#declaring-extras-optional-features-with-their-own-dependencies>`_
feature. (See the :doc:`../explore/declaring_requirements` docs for more detail
on how to use ``extras_require`` to create separate environments for GPU and
CPU.)

Mounted Home Directory
~~~~~~~~~~~~~~~~~~~~~~

``caliban notebook`` mounts your ``$HOME`` directory into the container, which
allows your Jupyter settings to persist across sessions. If you don't want this
for some reason, run the command with the ``--bare`` flag.

Custom Jupyer Port
~~~~~~~~~~~~~~~~~~

If you'd like to run ``notebook`` using a different port, use the ``--port`` option:

.. code-block:: bash

   caliban notebook --lab --port 8889

On the Mac you'll have to pass ``--nogpu`` to ``notebook``\ , as the NVIDIA runtime
isn't supported on non-Linux machines.
caliban build
^^^^^^^^^^^^^

This command builds the Docker image used in :doc:`caliban_run`,
:doc:`caliban_cloud` and friends, without actually executing the container or
submitting it remotely.

``caliban build`` supports the following arguments:

.. code-block:: text

   usage: caliban build [-h] [--helpfull] [--nogpu] [--cloud_key CLOUD_KEY]
                        [--extras EXTRAS] [-d DIR]
                        module

   positional arguments:
     module                Code to execute, in either 'trainer.train' or
                           'trainer/train.py' format. Accepts python scripts,
                           modules or a path to an arbitrary script.

   optional arguments:
     -h, --help            show this help message and exit
     --helpfull            show full help message and exit
     --nogpu               Disable GPU mode and force CPU-only.
     --cloud_key CLOUD_KEY
                           Path to GCloud service account key. (Defaults to
                           $GOOGLE_APPLICATION_CREDENTIALS.)
     --extras EXTRAS       setup.py dependency keys.
     --no_cache            Disable Docker's caching mechanism and force a
                           rebuild of the container from scratch.
     -d DIR, --dir DIR     Extra directories to include. List these from large to
                           small to take full advantage of Docker's build cache.
caliban stop
^^^^^^^^^^^^^^^^^^^^

This command allows you to stop running jobs submitted using caliban.

For example, suppose you submit a group of experiments to GKE using an
experiment config file like the following:

.. code-block::

   $ caliban cluster job submit --xgroup my-xgroup ... --experiment_config exp.json cpu.py --

After a bit, you realize that you made a coding error, so you'd like to stop
these jobs so that you can fix your error without wasting cloud resources (and
money). The ``caliban stop`` command makes this relatively simple:

.. code-block::

   $ caliban stop --xgroup my-xgroup
   the following jobs would be stopped:
   cpu.py --foo 3 --sleep -1
       job 61       RUNNING        GKE 2020-05-28 11:55:04 container: gcr.io/totoro-project/0f6d8a3ddbee:latest name: job-stop-test-57pr9
   cpu.py --foo 3 --sleep 2
       job 62       RUNNING        GKE 2020-05-28 11:55:04 container: gcr.io/totoro-project/0f6d8a3ddbee:latest name: job-stop-test-s67jt
   cpu.py --foo 3 --sleep 600
       job 63       RUNNING        GKE 2020-05-28 11:55:04 container: gcr.io/totoro-project/0f6d8a3ddbee:latest name: job-stop-test-gg9zm

   do you wish to stop these 3 jobs? [yN]: y

   stopping job: 61       RUNNING        GKE 2020-05-28 11:55:04 container: gcr.io/totoro-project/0f6d8a3ddbee:latest name: job-stop-test-57pr9
   stopping job: 62       RUNNING        GKE 2020-05-28 11:55:04 container: gcr.io/totoro-project/0f6d8a3ddbee:latest name: job-stop-test-s67jt
   stopping job: 63       RUNNING        GKE 2020-05-28 11:55:04 container: gcr.io/totoro-project/0f6d8a3ddbee:latest name: job-stop-test-gg9zm

   requested job cancellation, please be patient as it may take a short while for this status change to be reflected in the gcp dashboard or from the `caliban status` command.

This command will stop all jobs that are in a ``RUNNING`` or ``SUBMITTED`` state,
and checks with you to make sure this is what you *really* intend, as
accidentally stopping a job that has been running for days is a particularly
painful experience if your checkpointing is less than perfect. Similar to other
caliban commands, you can use the ``--dry_run`` flag to just print what jobs would
be stopped.

This command supports the following arguments:

.. code-block::

   $ caliban stop --help
   usage: caliban stop [-h] [--helpfull] [--xgroup XGROUP] [--dry_run]

   optional arguments:
     -h, --help       show this help message and exit
     --helpfull       show full help message and exit
     --xgroup XGROUP  experiment group
     --dry_run        Don't actually submit; log everything that's going to
                      happen.
caliban shell
^^^^^^^^^^^^^

This command is designed for fast, iterative workflows on scripts in an
environment that's guaranteed to match the environment available to your code on
Cloud.

``caliban shell`` supports the following arguments:

.. code-block:: text

   usage: caliban shell [-h] [--helpfull] [--nogpu] [--cloud_key CLOUD_KEY]
                        [--extras EXTRAS] [--image_id IMAGE_ID]
                        [--docker_run_args DOCKER_RUN_ARGS] [--shell {bash,zsh}]
                        [--bare]

   optional arguments:
     -h, --help            show this help message and exit
     --helpfull            show full help message and exit
     --nogpu               Disable GPU mode and force CPU-only.
     --cloud_key CLOUD_KEY
                           Path to GCloud service account key. (Defaults to
                           $GOOGLE_APPLICATION_CREDENTIALS.)
     --extras EXTRAS       setup.py dependency keys.
     --image_id IMAGE_ID   Docker image ID accessible in the local Docker
                           registry. If supplied, Caliban will skip the 'docker
                           build' step and use this image.
     --docker_run_args DOCKER_RUN_ARGS
                           String of args to add to Docker.
     --shell {bash,zsh}    This argument sets the shell used inside the container
                           to one of Caliban's supported shells. Defaults to the
                           shell specified by the $SHELL environment variable, or
                           'bash' if your shell isn't supported.
     --bare                Skip mounting the $HOME directory; load a bare shell.

Running ``caliban shell`` in any directory will generate a Docker image
containing the minimal environment necessary to execute Python ML workflows and
drop you into an interactive shell inside of that image.

Caliban will copy in your Cloud credentials and set the required
``$GOOGLE_APPLICATION_CREDENTIALS`` env variable, so all Cloud interaction from
Python should Just Work. (See the :doc:`guide on gcloud authentication
<../explore/gcloud>` for more detail.)

The base Caliban images also have ``gcloud`` installed; all ``gcloud`` and ``gsutil``
commands will work with the same permissions granted to the key found at
``$GOOGLE_APPLICATION_CREDENTIALS``.

.. NOTE:: If you run ``caliban shell --bare``\ , your gcloud and gsutil will
   have the same permissions that they'll have in the cloud - the permissions
   granted by your JSON key file. If you just run ``caliban shell``\ , which
   mounts your home directory, ``gcloud`` and ``gsutil`` will preferentially
   load the config you have on your local machine.

The only python dependencies available in the container will be dependencies
that you declare explicitly in either:


* a ``requirements.txt`` file
* a ``setup.py`` file.

Your setup file can declare groups of dependencies using the setuptools
`extras_require
<https://setuptools.readthedocs.io/en/latest/setuptools.html#declaring-extras-optional-features-with-their-own-dependencies>`_
feature. (See the :doc:`../explore/declaring_requirements` docs for more detail
on how to use ``extras_require`` to create separate environments for GPU and
CPU.)

By default your home directory will mount into the container, along with the
folder you're in when you run ``caliban shell``. This means that:


* your default ``bash`` (or ``zsh``\ ) environment will be available to you at the
  ``caliban shell``.
* Any changes you make to files in the mounted directory will be immediately
  available to you to run with, say, ``python -m trainer.train`` or some similar
  command.

On the Mac you'll have to pass ``--nogpu`` to ``shell``\ , as the NVIDIA runtime isn't
supported on non-Linux machines. If you forget ``caliban`` will remind you and
prevent you from getting too far.

.. NOTE:: Caliban currently supports ``bash`` and ``zsh`` shells. The command
   will use your ``$SHELL`` environment variable to pick a default; to override
   the default, you can always pass the ``--shell`` argument, like this:
   ``caliban shell --shell bash``.

.. WARNING:: One potential issue resulting from the fact that your home directory will mount
    into the container is that some binaries from your ``$HOME``  directory might
    leak into the container.  For example, we have seen a case in which, in trying
    to run a CUDA binary to communicate with the GPU, ``caliban shell`` called a
    binary from the home directory rather than the one which the container should
    have used. This issue can be mitigated simply by using the ``--bare`` option,
    which will not mount the ``$HOME``  directory inside the container.
caliban cloud
^^^^^^^^^^^^^

This command bundles your code and any other directories you specify into an
isolated Docker container and runs the resulting Python code on Google's
`AI Platform <https://cloud.google.com/ai-platform/>`_.

``caliban cloud`` supports the following arguments:

.. code-block:: text

   usage: caliban cloud [-h] [--helpfull] [--nogpu] [--cloud_key CLOUD_KEY]
                        [--extras EXTRAS] [-d DIR]
                        [--experiment_config EXPERIMENT_CONFIG] [--dry_run]
                        [--image_tag IMAGE_TAG] [--project_id PROJECT_ID]
                        [--region REGION] [--machine_type MACHINE_TYPE]
                        [--gpu_spec NUMxGPU_TYPE] [--tpu_spec NUMxTPU_TYPE]
                        [--force] [--name NAME] [-l KEY=VALUE]
                        module ...

   positional arguments:
     module                Code to execute, in either 'trainer.train' or
                           'trainer/train.py' format. Accepts python scripts,
                           modules or a path to an arbitrary script.

   optional arguments:
     -h, --help            show this help message and exit
     --helpfull            show full help message and exit
     --nogpu               Disable GPU mode and force CPU-only.
     --cloud_key CLOUD_KEY
                           Path to GCloud service account key. (Defaults to
                           $GOOGLE_APPLICATION_CREDENTIALS.)
     --extras EXTRAS       setup.py dependency keys.
     -d DIR, --dir DIR     Extra directories to include. List these from large to
                           small to take full advantage of Docker's build cache.
     --experiment_config EXPERIMENT_CONFIG
                           Path to an experiment config, or 'stdin' to read from
                           stdin.
     --dry_run             Don't actually submit; log everything that's going to
                           happen.
     --image_tag IMAGE_TAG
                           Docker image tag accessible via Container Registry. If
                           supplied, Caliban will skip the build and push steps
                           and use this image tag.
     --project_id PROJECT_ID
                           ID of the GCloud AI Platform project to use for Cloud
                           job submission and image persistence. (Defaults to
                           $PROJECT_ID; errors if both the argument and
                           $PROJECT_ID are empty.)
     --region REGION       Region to use for Cloud job submission and image
                           persistence. Must be one of ['us-east4', 'us-west1',
                           'us-west2', 'us-central1', 'us-east1', 'europe-west4',
                           'europe-west1', 'europe-north1', 'asia-northeast1',
                           'asia-east1', 'asia-southeast1']. (Defaults to $REGION
                           or 'us-central1'.)
     --machine_type MACHINE_TYPE
                           Cloud machine type to request. Must be one of
                           ['n1-standard-8', 'n1-highmem-4', 'n1-highcpu-96',
                           'n1-highcpu-16', 'n1-highcpu-32', 'n1-highcpu-64',
                           'n1-standard-96', 'n1-highmem-2', 'n1-highmem-16',
                           'n1-highmem-64', 'n1-standard-32', 'n1-standard-16',
                           'n1-standard-4', 'n1-highmem-8', 'n1-highmem-96',
                           'n1-standard-64', 'n1-highmem-32', 'cloud_tpu'].
                           Defaults to 'n1-standard-8' in GPU mode, or
                           'n1-highcpu-32' if --nogpu is passed.
     --gpu_spec NUMxGPU_TYPE
                           Type and number of GPUs to use for each AI Platform
                           submission. Defaults to 1xP100 in GPU mode or None if
                           --nogpu is passed.
     --tpu_spec NUMxTPU_TYPE
                           Type and number of TPUs to request for each AI
                           Platform submission. Defaults to None.
     --force               Force past validations and submit the job as
                           specified.
     --name NAME           Set a job name for AI Platform jobs.
     -l KEY=VALUE, --label KEY=VALUE
                           Extra label k=v pair to submit to Cloud.

   pass-through arguments:
     -- YOUR_ARGS          This is a catch-all for arguments you want to pass
                           through to your script. any arguments after '--' will
                           pass through.

.. NOTE:: To use ``caliban cloud`` you'll need to make sure your machine is
   configured for Cloud access. To verify that you're set up, visit
   :doc:`../getting_started/cloud`.

Specifically, you'll need to make sure the following environment variables are
set:

* ``$PROJECT_ID``\ : The ID of the Cloud project where you'll be submitting jobs.
* ``$GOOGLE_APPLICATION_CREDENTIALS``\ : a local path to your JSON Cloud
  credentials file.

``caliban cloud`` works almost exactly like ``caliban run`` (by design!). Thanks to
Docker, the container environment available to your job on AI Platform will look
exactly like the environment available on your local machine.

This means that if you can get your job running in ``caliban local`` mode you can
be quite sure that it'll complete in Cloud as well. The advantages of Cloud mode
are:


#. The machines are much bigger
#. Multi-GPU machines, clusters and TPUs are available
#. Cloud can execute many jobs in parallel, and will pipeline jobs for you as
   jobs complete and resources become available on your project.

See the ``caliban run`` docs for a detailed walkthrough of most options available
to ``caliban cloud``.

This mode has many features explored in the "Cloud-Specific Tutorials" section
in the left-hand menu. Read on here for a description of each keyword argument
supported by ``caliban cloud``.

Arguments as Labels
~~~~~~~~~~~~~~~~~~~

As with ``caliban run``\ , any arguments you pass to your script after ``--``\ :

.. code-block:: bash

   caliban cloud trainer.train -- --epochs 2

Will be passed directly through to your script.

In cloud mode, all user arguments will be passed to cloud as labels, which means
that you can filter by these labels in the AI platform jobs UI.

Keyword Arguments
~~~~~~~~~~~~~~~~~

The additional options available to ``caliban cloud`` are:


* **image_tag**\ : If you supply the tag of a Docker image accessible from your
  project, caliban will bypass the Docker build and push steps and use this
  image tag directly for AI Platform job submission. This is useful if you want
  to submit a job quickly without going through a no-op build and push, or if
  you want to :doc:`broadcast an experiment
  <../explore/experiment_broadcasting>` using some existing container. Note that
  this flag will cause ``caliban cloud`` to ignore any ``--extras`` or ``--dir``
  arguments, as no ``docker build`` step will be executed.

* **project_id**\ : This is the ID of the Cloud project that Caliban will use to
  push Docker containers and to submit AI platform jobs. By default Caliban will
  examine your environment for a ``$PROJECT_ID`` variable; if neither is set and
  you attempt to run a Cloud command, Caliban will exit.

* **region**\ : The Cloud region you specify with this flag is used for AI
  Platform job submission. Any value listed in `AI Platform's region docs
  <https://cloud.google.com/ml-engine/docs/regions>`_ is valid. If you don't
  specify a region Caliban will examine your environment for a ``$REGION``
  variable and use this if supplied; if that's not set it will default to
  ``"us-central1"``. See ``caliban cloud --help`` for all possible arguments.

* **machine_type**\ : Specifies the type of machine to use for each submitted AI
  platform job. See ``caliban cloud --help`` for all possible values. See
  :doc:`../cloud/gpu_specs` for more detail.

* **gpu_spec**\ : optional argument of the form GPU_COUNTxGPU_TYPE. See
  ``caliban cloud --help`` for all possible GPU types, and for the default.
  Usually 1, 2, 4 or 8 of each are supported, though this depends on the machine
  type you specify. Caliban will throw a validation error and give you a
  suggestion for how to proceed if you supply a combination that's not possible
  on AI Platform. See :doc:`../cloud/gpu_specs` for more details.

* **tpu_spec**\ : optional argument of the form TPU_COUNTxTPU_TYPE. See
  ``caliban cloud --help`` for all supported TPU types. As of December 2019,
  ``8xV2`` and ``8xV3`` are the only available options. TPUs are compatible with
  GPUs specified using ``--gpu_spec``. See :doc:`../cloud/ai_platform_tpu` for
  more details.

*
  **--force**\ : If supplied, this flag will disable all validations on
  combinations of region, machine type, GPU count and GPU type and force
  caliban to submit the job to AI Platform as specified. This is useful in
  case some new GPU was added to a region or machine type and caliban hasn't
  yet been updated.

* **name**\ : If you pass a string via this optional flag, ``caliban cloud``
  will submit your job with a job id of ``"{name}_{timestamp}"`` and add a
  ``job_name:{name}`` label to your job. It's useful to pass the same name for
  MANY jobs and use this field to group various experiment runs. Experiment
  broadcasting (the next flag, keep reading!) will do this for you
  automatically.

* **experiment_config**\ : If you pass the location (relative or absolute) of a
  local JSON file of the proper format, caliban will generate many jobs using
  this experiment config and submit them all in batch to AI platform. The
  formatting rules are - keys must be strings, values can be list, int, boolean
  or string. If the value is a list, caliban will generate N copies of the
  experiment config, 1 for each entry in the list, and submit a job for each.
  The total number of jobs submitted is the cardinality of the cartesian product
  of all lists in the experiment config. Lists of valid dicts are also allowed.
  See :doc:`../explore/experiment_broadcasting` for more details.

* **label**\ : You can use this flag to pass many labels to ``caliban cloud``\ ;
  just pass the flag over and over. Labels must be of the form ``k=v``\ ;
  ``--label epochs=2``\ , for example. If you pass any labels identical to your
  flags these labels will take precedence. See :doc:`../cloud/labels` below for
  more detail.

* **dry_run**\ : this flag will force logging output of all jobs that caliban
  will submit without the ``--dry_run`` flag. Docker will also skip an actual
  build and push. Use this to check that your other arguments are well formatted
  before submitting a potentially very large batch of jobs (depending on your
  experiment config).
expansion
^^^^^^^^^

The ``expansion`` script allows you to expand the :doc:`experiment.json
<../explore/experiment_broadcasting>` files accepted by the
``--experiment_config`` flags in subcommands like :doc:`caliban_run` and
:doc:`caliban_cloud` .


``expansion`` supports the following arguments:

.. code-block:: text

   usage: expansion [-h] [--helpfull] [--version] [--pprint] [--print_flags]
                    experiment_config

   Experiment config expander. For documentation, visit https://github.com/google/caliban

   positional arguments:
     experiment_config  Path to an experiment config, or 'stdin' to read from
                        stdin.

   optional arguments:
     -h, --help         show this help message and exit
     --helpfull         show full help message and exit
     --version          show program's version number and exit
     --pprint           Pretty-print the config to stdout.
     --print_flags      Print the actual flags generated by each experiment in
                        the expansion, one per line.

Given a file called ``experiment.json`` with these contents:

.. code-block:: json

   [    {
           "epochs": [1,2, 3, 4],
           "batch_size": [64, 128],
           "constant_arg": "something",
           "important_toggle": [true, false]
       },
       {
           "epochs": 1000,
           "batch_size": 1
       }
   ]

Running ``expansion experiment.json`` will generate a NON-pretty-printed string
containing every combination of flag key and value described by the structure
above (See the :doc:`../explore/experiment_broadcasting` page for details):

.. code-block:: json

   [{"epochs": 1, "batch_size": 64, "constant_arg": "something", "important_toggle": true}, {"epochs": 1, "batch_size": 64, "constant_arg": "something", "important_toggle": false}, {"epochs": 1, "batch_size": 128, "constant_arg": "something", "important_toggle": true}, {"epochs": 1, "batch_size": 128, "constant_arg": "something", "important_toggle": false}, {"epochs": 2, "batch_size": 64, "constant_arg": "something", "important_toggle": true}, {"epochs": 2, "batch_size": 64, "constant_arg": "something", "important_toggle": false}, {"epochs": 2, "batch_size": 128, "constant_arg": "something", "important_toggle": true}, {"epochs": 2, "batch_size": 128, "constant_arg": "something", "important_toggle": false}, {"epochs": 3, "batch_size": 64, "constant_arg": "something", "important_toggle": true}, {"epochs": 3, "batch_size": 64, "constant_arg": "something", "important_toggle": false}, {"epochs": 3, "batch_size": 128, "constant_arg": "something", "important_toggle": true}, {"epochs": 3, "batch_size": 128, "constant_arg": "something", "important_toggle": false}, {"epochs": 4, "batch_size": 64, "constant_arg": "something", "important_toggle": true}, {"epochs": 4, "batch_size": 64, "constant_arg": "something", "important_toggle": false}, {"epochs": 4, "batch_size": 128, "constant_arg": "something", "important_toggle": true}, {"epochs": 4, "batch_size": 128, "constant_arg": "something", "important_toggle": false}, {"epochs": 1000, "batch_size": 1}]

You can also pipe the ``experiment.json`` file in via stdin by passing ``stdin``
instead of a filename:

.. code-block:: bash

   cat experiment.json | expansion stdin

The :doc:`../explore/experiment_broadcasting` page, specifically its section on
pipes, has more information on how to use this function with other caliban
commands to generate very complex sets of experiments with ease.

expansion --pprint
~~~~~~~~~~~~~~~~~~

The ``--pprint`` argument forces a more sane expansion. Here are the results of
``expansion experiment.json --pprint``\ :

.. code-block:: json

   [
     {
       "epochs": 1,
       "batch_size": 64,
       "constant_arg": "something",
       "important_toggle": true
     },
     {
       "epochs": 1,
       "batch_size": 64,
       "constant_arg": "something",
       "important_toggle": false
     },
     {
       "epochs": 1,
       "batch_size": 128,
       "constant_arg": "something",
       "important_toggle": true
     },
     // etc etc
   ]

expansion --print_flags
~~~~~~~~~~~~~~~~~~~~~~~

The ``--print_flags`` argument goes one step further and prints the actual
argparse flags that correspond to each of the expanded JSON objects. Here are
the results of ``expansion experiment.json --print_flags``\ :

.. code-block:: text

   --epochs 1 --batch_size 64 --constant_arg something --important_toggle
   --epochs 1 --batch_size 64 --constant_arg something
   --epochs 1 --batch_size 128 --constant_arg something --important_toggle
   --epochs 1 --batch_size 128 --constant_arg something
   --epochs 2 --batch_size 64 --constant_arg something --important_toggle
   --epochs 2 --batch_size 64 --constant_arg something
   --epochs 2 --batch_size 128 --constant_arg something --important_toggle
   --epochs 2 --batch_size 128 --constant_arg something
   --epochs 3 --batch_size 64 --constant_arg something --important_toggle
   --epochs 3 --batch_size 64 --constant_arg something
   --epochs 3 --batch_size 128 --constant_arg something --important_toggle
   --epochs 3 --batch_size 128 --constant_arg something
   --epochs 4 --batch_size 64 --constant_arg something --important_toggle
   --epochs 4 --batch_size 64 --constant_arg something
   --epochs 4 --batch_size 128 --constant_arg something --important_toggle
   --epochs 4 --batch_size 128 --constant_arg something
   --epochs 1000 --batch_size 1
caliban resubmit
^^^^^^^^^^^^^^^^^^^^^^^^

Often one needs to re-run an experiment after making code changes, or to run the
same code with a different random seed. Caliban supports this with its
``resubmit`` command.

This command allows you to resubmit jobs in an experiment group without having
to remember or re-enter all of the parameters for your experiments. For example,
suppose you run a set of experiments in an experiment group on CAIP:

.. code-block::

   caliban cloud --xgroup resubmit_test --nogpu --experiment_config experiment.json cpu.py -- --foo 3

You then realize that you made a coding error, causing some of your jobs to
fail:

.. code-block::

   $ caliban status --xgroup resubmit_test
   xgroup resubmit_test:
   docker config 1: job_mode: CPU, build url: ~/sw/cluster/caliban/tmp/cpu, extra dirs: None
     experiment id 37: cpu.py --foo 3 --sleep 2
       job 69       SUCCEEDED     CAIP 2020-05-29 10:53:41 container: gcr.io/totoro-project/cffd1475aaca:latest name: caliban_totoro_20200529_105340_2
     experiment id 38: cpu.py --foo 3 --sleep 1
       job 68       FAILED        CAIP 2020-05-29 10:53:40 container: gcr.io/totoro-project/cffd1475aaca:latest name: caliban_totoro_20200529_105338_1

You then go and modify your code, and now you can use the ``resubmit`` command to
run the jobs that failed:

.. code-block::

   $ caliban resubmit --xgroup resubmit_test
   the following jobs would be resubmitted:
   cpu.py --foo 3 --sleep 1
     job 68       FAILED        CAIP 2020-05-29 10:53:40 container: gcr.io/totoro-project/cffd1475aaca:latest name: caliban_totoro_20200529_105338_1

    do you wish to resubmit these 1 jobs? [yN]: y
   rebuilding containers...
   ...
   Submitting request!
   ...

Checking back in with ``caliban status`` shows that the code change worked, and
now all of the experiments in the group have succeeded, and you can see that the
container hash has changed for the previously failed jobs, reflecting your code
change:

.. code-block::

   $ caliban status --xgroup resubmit_test
   xgroup resubmit_test:
   docker config 1: job_mode: CPU, build url: ~/sw/cluster/caliban/tmp/cpu, extra dirs: None
     experiment id 37: cpu.py --foo 3 --sleep 2
       job 69       SUCCEEDED     CAIP 2020-05-29 10:53:41 container: gcr.io/totoro-project/cffd1475aaca:latest name: caliban_totoro_20200529_105340_2
     experiment id 38: cpu.py --foo 3 --sleep 1
       job 70       SUCCEEDED     CAIP 2020-05-29 11:03:01 container: gcr.io/totoro-project/81b2087b5026:latest name: caliban_totoro_20200529_110259_1

The ``resubmit`` command supports the following arguments:

.. code-block::

   $ caliban resubmit --help
   usage: caliban resubmit [-h] [--helpfull] [--xgroup XGROUP] [--dry_run] [--all_jobs] [--project_id PROJECT_ID] [--cloud_key CLOUD_KEY]

   optional arguments:
     -h, --help            show this help message and exit
     --helpfull            show full help message and exit
     --xgroup XGROUP       experiment group
     --dry_run             Don't actually submit; log everything that's going to happen.
     --all_jobs            resubmit all jobs regardless of current state, otherwise only jobs that are in FAILED or STOPPED state will be resubmitted
     --project_id PROJECT_ID
                           ID of the GCloud AI Platform/GKE project to use for Cloud job submission and image persistence. (Defaults to $PROJECT_ID; errors if both the argument and $PROJECT_ID are empty.)
     --cloud_key CLOUD_KEY
                           Path to GCloud service account key. (Defaults to $GOOGLE_APPLICATION_CREDENTIALS.)
caliban cluster
^^^^^^^^^^^^^^^

This subcommand allows you to create and submit jobs to a GKE cluster using
caliban's packaging and interface features.

``caliban cluster ls``
~~~~~~~~~~~~~~~~~~~~~~~~~~

This command lists the clusters currently available in your project.

.. code-block:: text

   usage: caliban cluster ls [-h] [--helpfull] [--project_id PROJECT_ID]
                             [--cloud_key CLOUD_KEY] [--zone ZONE]

   list clusters

   optional arguments:
     -h, --help            show this help message and exit
     --helpfull            show full help message and exit
     --project_id PROJECT_ID
                           ID of the GCloud AI Platform/GKE project to use for
                           Cloud job submission and image persistence. (Defaults
                           to $PROJECT_ID; errors if both the argument and
                           $PROJECT_ID are empty.) (default: None)
     --cloud_key CLOUD_KEY
                           Path to GCloud service account key. (Defaults to
                           $GOOGLE_APPLICATION_CREDENTIALS.) (default: None)
     --zone ZONE           zone (default: None)

Here you may specify a specific project, credentials file, or cloud zone to
narrow your listing. If you do not specify these, caliban tries to determine
these from the system defaults.

``caliban cluster create``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This command creates a new cluster in your project. Typically if you are going
to use GKE in your project, you will create a single long-running cluster in
your project first, and leave it running across many job submissions. In caliban
we configure the cluster to take advantage of autoscaling wherever possible.

In GKE, there are two types of autoscaling. The first is known as
`'cluster autoscaling' <https://cloud.google.com/kubernetes-engine/docs/concepts/cluster-autoscaler>`_.
This mode automatically increases the number of nodes in your cluster's node
pools as job demand increases. In caliban, we configure this automatically and
we query your cpu and accelerator quota to configure the cluster autoscaling
limits. In this way, your cluster will automatically add nodes when you need
them, and then automatically delete them when they are no longer needed. This
is, of course, quite useful for keeping your costs low.

The second type of autoscaling in GKE is
`'node autoprovisioning' <https://cloud.google.com/kubernetes-engine/docs/how-to/node-auto-provisioning>`_.
This form of autoprovisioning addresses the issue that accelerator-enabled
instances must be allocated from a node pool of instances where the particular
cpu/memory/gpu configuration is fixed. For simple configurations where you
support only a small number of node configurations, you can manually create
autoscaling node pools. If, however, you wish to support several, or in
caliban's case, general, configurations, then this becomes more difficult. Node
autoprovisioning automatically creates autoscaling node pools based on the
requirements of the jobs submitted to the cluster, and also deletes these node
pools once they are no longer needed. In caliban we enable node autoprovisioning
so you can specify your gpu- and machine- types on a per-job basis, and the
kubernetes engine will automatically create the appropriate node pools to
accomodate your jobs.

The syntax for this command is as follows:

.. code-block:: text

   totoro@totoro:$ caliban cluster create --help
   usage: caliban cluster create [-h] [--helpfull] [--project_id PROJECT_ID]
                                 [--cloud_key CLOUD_KEY]
                                 [--cluster_name CLUSTER_NAME] [--zone ZONE]
                                 [--dry_run]
                                 [--release_channel ['UNSPECIFIED', 'RAPID', 'REGULAR', 'STABLE']]
                                 [--single_zone]

   create cluster

   optional arguments:
     -h, --help            show this help message and exit
     --helpfull            show full help message and exit
     --project_id PROJECT_ID
                           ID of the GCloud AI Platform/GKE project to use for
                           Cloud job submission and image persistence. (Defaults
                           to $PROJECT_ID; errors if both the argument and
                           $PROJECT_ID are empty.) (default: None)
     --cloud_key CLOUD_KEY
                           Path to GCloud service account key. (Defaults to
                           $GOOGLE_APPLICATION_CREDENTIALS.) (default: None)
     --cluster_name CLUSTER_NAME
                           cluster name (default: None)
     --zone ZONE           for a single-zone cluster, this specifies the zone for
                           the cluster control plane and all worker nodes, while
                           for a multi-zone cluster this specifies only the zone
                           for the control plane, while worker nodes may be
                           created in any zone within the same region as the
                           control plane. The single_zone argument specifies
                           whether to create a single- or multi- zone cluster.
                           (default: None)
     --dry_run             Don't actually submit; log everything that's going to
                           happen. (default: False)
     --release_channel ['UNSPECIFIED', 'RAPID', 'REGULAR', 'STABLE']
                           cluster release channel, see
                           https://cloud.google.com/kubernetes-
                           engine/docs/concepts/release-channels (default: REGULAR)
     --single_zone         create a single-zone cluster if set, otherwise create
                           a multi-zone cluster: see
                           https://cloud.google.com/kubernetes-
                           engine/docs/concepts/types-of-
                           clusters#cluster_availability_choices (default: False)

You can use the ``--dry_run`` flag to see the specification for the cluster that
would be submitted to GKE without actually creating the cluster.

A typical creation request (with ``--dry_run``\ ):

.. code-block:: text

   totoro@totoro:$ caliban cluster create --zone us-central1-a --cluster_name newcluster --dry_run
   I0303 13:07:34.257717 140660011796288 cli.py:160] request:
   {'cluster': {'autoscaling': {'autoprovisioningNodePoolDefaults': {'oauthScopes': ['https://www.googleapis.com/auth/compute',
                                                                                     'https://www.googleapis.com/auth/cloud-platform']},
                                'enableNodeAutoprovisioning': 'true',
                                'resourceLimits': [{'maximum': '72',
                                                    'resourceType': 'cpu'},
                                                   {'maximum': '4608',
                                                    'resourceType': 'memory'},
                                                   {'maximum': '8',
                                                    'resourceType': 'nvidia-tesla-k80'},
                                                   {'maximum': '1',
                                                    'resourceType': 'nvidia-tesla-p100'},
                                                   {'maximum': '1',
                                                    'resourceType': 'nvidia-tesla-v100'},
                                                   {'maximum': '1',
                                                    'resourceType': 'nvidia-tesla-p4'},
                                                   {'maximum': '4',
                                                    'resourceType': 'nvidia-tesla-t4'}]},
                'enable_tpu': 'true',
                'ipAllocationPolicy': {'useIpAliases': 'true'},
                'locations': ['us-central1-a',
                              'us-central1-b',
                              'us-central1-c',
                              'us-central1-f'],
                'name': 'newcluster',
                'nodePools': [{'config': {'oauthScopes': ['https://www.googleapis.com/auth/devstorage.read_only',
                                                          'https://www.googleapis.com/auth/logging.write',
                                                          'https://www.googleapis.com/auth/monitoring',
                                                          'https://www.googleapis.com/auth/service.management.readonly',
                                                          'https://www.googleapis.com/auth/servicecontrol',
                                                          'https://www.googleapis.com/auth/trace.append']},
                               'initialNodeCount': '3',
                               'name': 'default-pool'}],
                'releaseChannel': {'channel': 'RAPID'},
                'zone': 'us-central1-a'},
    'parent': 'projects/totoro-project/locations/us-central1-a'}

Cluster creation can take a while to complete (often on the order of five
minutes). When you use caliban to create a cluster, caliban will provide a link
to the relevant GCP dashboard page where you can monitor the progress of your
cluster creation request. Caliban will also monitor your creation request, and
when your cluster is created, it will apply a
`daemonset <https://kubernetes.io/docs/concepts/workloads/controllers/daemonset/>`_
to your cluster to automatically apply nvidia drivers to any gpu-enabled nodes
that get created, as described
`here <https://cloud.google.com/kubernetes-engine/docs/how-to/gpus#installing_drivers>`_.

``caliban cluster delete``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This command simply deletes an existing cluster. Typically you will leave your
cluster running, but the cluster does consume some resources even when idle, so
if you are not actively using the cluster you may want to shut it down to save
money.

The syntax of this command:

.. code-block:: text

   totoro@totoro:$ caliban cluster delete --help
   usage: caliban cluster delete [-h] [--helpfull] [--project_id PROJECT_ID]
                                 [--cloud_key CLOUD_KEY]
                                 [--cluster_name CLUSTER_NAME] [--zone ZONE]

   delete cluster

   optional arguments:
     -h, --help            show this help message and exit
     --helpfull            show full help message and exit
     --project_id PROJECT_ID
                           ID of the GCloud AI Platform/GKE project to use for
                           Cloud job submission and image persistence. (Defaults
                           to $PROJECT_ID; errors if both the argument and
                           $PROJECT_ID are empty.) (default: None)
     --cloud_key CLOUD_KEY
                           Path to GCloud service account key. (Defaults to
                           $GOOGLE_APPLICATION_CREDENTIALS.) (default: None)
     --cluster_name CLUSTER_NAME
                           cluster name (default: None)
     --zone ZONE           zone (default: None)

As with most caliban commands, if you do not specify arguments, then caliban
does its best to determine them from defaults. For example, if you have only a
single cluster in your project, you can simply type ``caliban cluster delete``.

``caliban cluster job submit``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most of the cli arguments for ``caliban cluster job submit`` are the same as
those for :doc:`../cli/caliban_cloud`:

.. code-block:: text

   totoro@totoro:$ caliban cluster job submit --help
   usage: caliban cluster job submit [-h] [--helpfull]
                                  [--cluster_name CLUSTER_NAME] [--nogpu]
                                  [--cloud_key CLOUD_KEY] [--extras EXTRAS]
                                  [-d DIR] [--image_tag IMAGE_TAG]
                                  [--project_id PROJECT_ID]
                                  [--min_cpu MIN_CPU] [--min_mem MIN_MEM]
                                  [--gpu_spec NUMxGPU_TYPE]
                                  [--tpu_spec NUMxTPU_TYPE]
                                  [--tpu_driver TPU_DRIVER]
                                  [--nonpreemptible_tpu] [--force]
                                  [--name NAME]
                                  [--experiment_config EXPERIMENT_CONFIG]
                                  [-l KEY=VALUE] [--nonpreemptible]
                                  [--dry_run] [--export EXPORT]
                                  [--xgroup XGROUP]
                                  module ...

   submit cluster job(s)

   positional arguments:
     module                Code to execute, in either trainer.train' or
                           'trainer/train.py' format. Accepts python scripts,
                           modules or a path to an arbitrary script.

   optional arguments:
     -h, --help            show this help message and exit
     --helpfull            show full help message and exit
     --cluster_name CLUSTER_NAME
                           cluster name (default: None)
     --nogpu               Disable GPU mode and force CPU-only. (default: True)
     --cloud_key CLOUD_KEY
                           Path to GCloud service account key. (Defaults to
                           $GOOGLE_APPLICATION_CREDENTIALS.) (default: None)
     --extras EXTRAS       setup.py dependency keys. (default: None)
     -d DIR, --dir DIR     Extra directories to include. List these from large to
                           small to take full advantage of Docker's build cache.
                           (default: None)
     --image_tag IMAGE_TAG
                           Docker image tag accessible via Container Registry. If
                           supplied, Caliban will skip the build and push steps
                           and use this image tag. (default: None)
     --project_id PROJECT_ID
                           ID of the GCloud AI Platform/GKE project to use for
                           Cloud job submission and image persistence. (Defaults
                           to $PROJECT_ID; errors if both the argument and
                           $PROJECT_ID are empty.) (default: None)
     --min_cpu MIN_CPU     Minimum cpu needed by job, in milli-cpus. If not
                           specified, then this value defaults to 1500 for
                           gpu/tpu jobs, and 31000 for cpu jobs. Please note that
                           gke daemon processes utilize a small amount of cpu on
                           each node, so if you want to have your job run on a
                           specific machine type, say a 2-cpu machine, then if
                           you specify a minimum cpu of 2000, then your job will
                           not be schedulable on a 2-cpu machine as the daemon
                           processes will push the total cpu needed to more than
                           two full cpus. (default: None)
     --min_mem MIN_MEM     Minimum memory needed by job, in MB. Please note that
                           gke daemon processes utilize a small amount of memory
                           on each node, so if you want to have your job run on a
                           specific machine type, say a machine with 8GB total
                           memory, then if you specify a minimum memory of
                           8000MB, then your job will not be schedulable on a 8GB
                           machine as the daemon processes will push the total
                           memory needed to more than 8GB. (default: None)
     --gpu_spec NUMxGPU_TYPE
                           Type and number of GPUs to use for each AI
                           Platform/GKE submission. Defaults to 1xP100 in GPU
                           mode or None if --nogpu is passed. (default: None)
     --tpu_spec NUMxTPU_TYPE
                           Type and number of TPUs to request for each AI
                           Platform/GKE submission. Defaults to None. (default:
                           None)
     --tpu_driver TPU_DRIVER
                           tpu driver (default: 1.14)
     --nonpreemptible_tpu  use non-preemptible tpus: note this only applies to
                           v2-8 and v3-8 tpus currently, see:
                           https://cloud.google.com/tpu/docs/preemptible
                           (default: False)
     --force               Force past validations and submit the job as
                           specified. (default: False)
     --name NAME           Set a job name for AI Platform or GKE jobs. (default:
                           None)
     --experiment_config EXPERIMENT_CONFIG
                           Path to an experiment config, or 'stdin' to read from
                           stdin. (default: None)
     -l KEY=VALUE, --label KEY=VALUE
                           Extra label k=v pair to submit to Cloud. (default:
                           None)
     --nonpreemptible      use non-preemptible VM instance: please note that you
                           may need to upgrade your cluster to a recent
                           version/use the rapid release channel for preemptible
                           VMs to be supported with node autoprovisioning:
                           https://cloud.google.com/kubernetes-
                           engine/docs/release-notes-rapid#december_13_2019
                           (default: False)
     --dry_run             Don't actually submit; log everything that's going to
                           happen. (default: False)
     --export EXPORT       Export job spec(s) to file, extension must be one of
                           ('.yaml', '.json') (for example: --export my-job-
                           spec.yaml) For multiple jobs (i.e. in an experiment
                           config scenario), multiple files will be generated
                           with an index inserted (for example: --export my-job-
                           spec.yaml would yield my-job-spec_0.yaml, my-job-
                           spec_1.yaml...) (default: None)
     --xgroup XGROUP       This specifies an experiment group, which ties
                           experiments and job instances together. If you do not
                           specify a group, then a new one will be created. If
                           you specify an existing experiment group here, then
                           new experiments and jobs you create will be added to
                           the group you specify. (default: None)

   pass-through arguments:
     -- YOUR_ARGS          This is a catch-all for arguments you want to pass
                           through to your script. any arguments after '--' will
                           pass through.

Again, this command very closely mirrors :doc:`../cli/caliban_cloud`.

You can export job requests created with caliban as a ``yaml`` or ``json`` file
using the ``--export`` flag. You can then use this file with ``caliban cluster job
submit_file`` or
`\ ``kubectl`` <https://kubernetes.io/docs/concepts/workloads/controllers/jobs-run-to-completion/#running-an-example-job>`_
to submit the same job again.

``caliban cluster job submit_file``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This command submits a kubernetes k8s job file to your cluster. This can be
useful if you have a job that you run regularly, as you can create the job
initially with ``caliban cluster job submit`` and use the ``--export`` option to
save the job spec file. Then you can use this command to submit the job again
without having to specify all of the cli arguments.

The syntax of this command:

.. code-block:: text

   totoro@totoro:$ caliban cluster job submit_file --help
   usage: caliban cluster job submit_file [-h] [--helpfull]
                                          [--cluster_name CLUSTER_NAME]
                                          [--cloud_key CLOUD_KEY]
                                          [--project_id PROJECT_ID] [--dry_run]
                                          job_file

   submit gke job from yaml/json file

   positional arguments:
     job_file              kubernetes k8s job file ('.yaml', '.json')

   optional arguments:
     -h, --help            show this help message and exit
     --helpfull            show full help message and exit
     --cluster_name CLUSTER_NAME
                           cluster name (default: None)
     --cloud_key CLOUD_KEY
                           Path to GCloud service account key. (Defaults to
                           $GOOGLE_APPLICATION_CREDENTIALS.) (default: None)
     --project_id PROJECT_ID
                           ID of the GCloud AI Platform/GKE project to use for
                           Cloud job submission and image persistence. (Defaults
                           to $PROJECT_ID; errors if both the argument and
                           $PROJECT_ID are empty.) (default: None)
     --dry_run             Don't actually submit; log everything that's going to
                           happen. (default: False)

Thus a common invocation would resemble:

.. code-block:: text

   caliban cluster job submit_file my-job.yaml
caliban run
^^^^^^^^^^^

This command bundles your code and any other directories you specify into an
isolated Docker container and runs the resulting Python code on your local
machine, but inside of the Docker environment.

``caliban run`` supports the following arguments:

.. code-block:: text

   usage: caliban run [-h] [--helpfull] [--nogpu] [--cloud_key CLOUD_KEY]
                      [--extras EXTRAS] [-d DIR]
                      [--experiment_config EXPERIMENT_CONFIG] [--dry_run]
                      [--image_id IMAGE_ID] [--docker_run_args DOCKER_RUN_ARGS]
                      module ...

   positional arguments:
     module                Code to execute, in either 'trainer.train' or
                           'trainer/train.py' format. Accepts python scripts,
                           modules or a path to an arbitrary script.

   optional arguments:
     -h, --help            show this help message and exit
     --helpfull            show full help message and exit
     --nogpu               Disable GPU mode and force CPU-only.
     --cloud_key CLOUD_KEY
                           Path to GCloud service account key. (Defaults to
                           $GOOGLE_APPLICATION_CREDENTIALS.)
     --extras EXTRAS       setup.py dependency keys.
     -d DIR, --dir DIR     Extra directories to include. List these from large to
                           small to take full advantage of Docker's build cache.
     --experiment_config EXPERIMENT_CONFIG
                           Path to an experiment config, or 'stdin' to read from
                           stdin.
     --dry_run             Don't actually submit; log everything that's going to
                           happen.
     --image_id IMAGE_ID   Docker image ID accessible in the local Docker
                           registry. If supplied, Caliban will skip the 'docker
                           build' step and use this image.
     --docker_run_args DOCKER_RUN_ARGS
                           String of args to add to Docker.

   pass-through arguments:
     -- YOUR_ARGS          This is a catch-all for arguments you want to pass
                           through to your script. any arguments after '--' will
                           pass through.

Because the container is completely isolated, to get any results from ``caliban
run`` you'll have to depend on either:


* ``stdout``\ , if you're just interested in checking if the job is running at all
  before submission to Cloud, for example, or
* Cloud buckets for persistence.

Your credentials are set up inside the container and available via the required
``$GOOGLE_APPLICATION_CREDENTIALS`` environment variable, so all Cloud access
via Python should Just Work. (See the :doc:`guide on gcloud authentication
<../explore/gcloud>` for more detail.)

The base Caliban images also have ``gcloud`` installed; all ``gcloud`` and ``gsutil``
commands will work with the same permissions granted to the key found at
``$GOOGLE_APPLICATION_CREDENTIALS``.

As with the other commands, the only python dependencies available in the
container will be dependencies that you declare explicitly in either:


* a ``requirements.txt`` file
* a ``setup.py`` file.

Your setup file can declare groups of dependencies using the setuptools
`extras_require
<https://setuptools.readthedocs.io/en/latest/setuptools.html#declaring-extras-optional-features-with-their-own-dependencies>`_
feature. (See the :doc:`../explore/declaring_requirements` docs for more detail
on how to use ``extras_require`` to create separate environments for GPU and
CPU.)


Executing a Python script
~~~~~~~~~~~~~~~~~~~~~~~~~

The most basic way to trigger a run is by passing a file path or a module name.
Any of the following will work:

.. code-block:: bash

   caliban run trainer.train -- --epochs 2
   caliban run trainer/train.py
   caliban run mycode.py
   caliban run mycode -- --learning_rates '[1,2,3]'

Any flags or commands you pass to the command after ``--`` will be passed along,
untouched, to your Python code. By configuring your job with flags you can get a
large range of behavior out of the same module.

If you specify a Python module inside of a folder, ``caliban run`` will copy only
that folder into the Docker environment. For example, if you have a
``trainer/train.py`` file and run either of the following:

.. code-block:: bash

   caliban run trainer.train
   caliban run trainer/train.py

Caliban will copy only the ``trainer`` directory into the container.

If your script lives in the root of the directory, as in the ``mycode.py`` example
above, the entire current working directory will be copied in.

This could be inefficient if your directory has lots of data you don't want, or
a folder of notebooks; if you want a smaller build image you can move your
script into a folder. Make sure to create ``__init__.py`` inside the folder to
make it a proper module.

In addition to the required module name, ``caliban run`` supports many optional
arguments. All of these must be supplied **before** the module name.

Jobs run in GPU mode by default. To toggle GPU mode off, use ``--nogpu``.

Extra Directories
~~~~~~~~~~~~~~~~~

If you want to make extra directories available inside your container, pass them
like this:

.. code-block:: bash

   caliban -d data -d models/stored trainer.train

This invocation will copy the ``data`` and ``models/stored`` directories into the
container, where they can be accessed using a relative path. All directories
must exist relative to the directory where you run ``caliban run``.
caliban status
^^^^^^^^^^^^^^^^^^^^^^

The ``caliban status`` command allows you to check on the status of jobs submitted
via caliban. There are two primary modes for this command. The first returns
your most recent job submissions across all experiment groups:

.. code-block::

   $ caliban status --max_jobs 5
   most recent 5 jobs for user totoro:

   xgroup totoro-xgroup-2020-05-28-11-33-35:
     docker config 1: job_mode: CPU, build url: ~/sw/cluster/caliban/tmp/cpu, extra dirs: None
      experiment id 28: cpu.py --foo 3 --sleep 2
        job 56       STOPPED        GKE 2020-05-28 11:33:35 container: gcr.io/totoro-project/0f6d8a3ddbee:latest name: job-stop-test-rssqq
      experiment id 29: cpu.py --foo 3 --sleep 600
        job 57       STOPPED        GKE 2020-05-28 11:33:36 container: gcr.io/totoro-project/0f6d8a3ddbee:latest name: job-stop-test-c5x6v

   xgroup totoro-xgroup-2020-05-28-11-40-52:
     docker config 1: job_mode: CPU, build url: ~/sw/cluster/caliban/tmp/cpu, extra dirs: None
       experiment id 30: cpu.py --foo 3 --sleep -1
         job 58       STOPPED       CAIP 2020-05-28 11:40:54 container: gcr.io/totoro-project/0f6d8a3ddbee:latest name: caliban_totoro_20200528_114052_1
       experiment id 31: cpu.py --foo 3 --sleep 2
         job 59       STOPPED       CAIP 2020-05-28 11:40:55 container: gcr.io/totoro-project/0f6d8a3ddbee:latest name: caliban_totoro_20200528_114054_2
       experiment id 32: cpu.py --foo 3 --sleep 600
         job 60       RUNNING       CAIP 2020-05-28 11:40:56 container: gcr.io/totoro-project/0f6d8a3ddbee:latest name: caliban_totoro_20200528_114055_3

Here we can see five jobs that we recently submitted, in two experiment groups.
The first experiment group has jobs submitted to GKE, while the second has jobs
submitted to CAIP. You can specify the maximum number of jobs to return using
the ``--max_jobs`` flag.

The second mode for the ``caliban status`` command returns jobs in a given
experiment group, using the ``--xgroup`` flag:

.. code-block::

   $ caliban status --xgroup xg2 --max_jobs 2
   xgroup xg2:
   docker config 1: job_mode: CPU, build url: ~/sw/cluster/caliban/tmp/cpu, extra dirs: None
     experiment id 1: cpu.py --foo 3 --sleep -1
       job 34       FAILED        CAIP 2020-05-08 18:26:56 container: gcr.io/totoro-project/e2a0b8fca1dc:latest name: caliban_totoro_1_20200508_182654
       job 37       FAILED        CAIP 2020-05-08 19:01:08 container: gcr.io/totoro-project/e2a0b8fca1dc:latest name: caliban_totoro_1_20200508_190107
     experiment id 2: cpu.py --foo 3 --sleep 2
       job 30       SUCCEEDED    LOCAL 2020-05-08 09:59:04 container: e2a0b8fca1dc
       job 35       SUCCEEDED     CAIP 2020-05-08 18:26:57 container: gcr.io/totoro-project/e2a0b8fca1dc:latest name: caliban_totoro_2_20200508_182656
     experiment id 5: cpu.py --foo 3 --sleep 600
       job 36       STOPPED       CAIP 2020-05-08 18:26:58 container: gcr.io/totoro-project/e2a0b8fca1dc:latest name: caliban_totoro_3_20200508_182657
       job 38       SUCCEEDED     CAIP 2020-05-08 19:01:09 container: gcr.io/totoro-project/e2a0b8fca1dc:latest name: caliban_totoro_3_20200508_190108

Here we can see the jobs that have been submitted as part of the ``xg2``
experiment group. By specifying ``--max_jobs 2`` in the call, we can see the two
most recent job submissions for each experiment in the group. In this case, we
can see that experiment 2 was submitted both locally and to CAIP at different
times. We can also see that experiment 1 failed (due to an invalid parameter),
and that the first submision to CAIP of experiment 5 was stopped by the user.

Another interesting thing to note here is that the container hash is the same
for each of these job submissions, so we can tell that the underlying code did
not change between submissions.

This command supports the following arguments:

.. code-block::

   $ caliban status --help
   usage: caliban status [-h] [--helpfull] [--xgroup XGROUP]
                         [--max_jobs MAX_JOBS]

   optional arguments:
     -h, --help           show this help message and exit
     --helpfull           show full help message and exit
     --xgroup XGROUP      experiment group
     --max_jobs MAX_JOBS  Maximum number of jobs to view. If you specify an
                          experiment group, then this specifies the maximum
                          number of jobs per experiment to view. If you do not
                          specify an experiment group, then this specifies the
                          total number of jobs to return, ordered by creation
                          date, or all jobs if max_jobs==0.
TPUs on AI Platform
^^^^^^^^^^^^^^^^^^^

.. NOTE:: This documentation is currently quite sparse; expect a tutorial soon.

.. IMPORTANT:: Unlike on Cloud, TPUs on AI Platform only support (as of
   Dec 2019) Tensorflow versions 1.13 and 1.14. No Jax, no Pytorch.

Caliban has Tensorflow version 2.1 hardcoded internally. Once the range of
possible values expands we'll make this customizable.

See `AI Platform's runtime version list
<https://cloud.google.com/ml-engine/docs/runtime-version-list>`_ for more
detail.


If you supply the ``--tpu_spec NUM_TPUSxTPU_TYPE`` argument to your ``caliban
cloud`` job, AI Platform will configure a worker node with that number of TPUs
and attach it to the master node where your code runs.

``--tpu_spec`` is compatible with ``--gpu_spec``\ ; the latter configures the master
node where your code lives, while the former sets up a separate worker instance.

CPU mode by Default
~~~~~~~~~~~~~~~~~~~

Normally, all jobs default to GPU mode unless you supply ``--nogpu`` explicitly.
This default flips when you supply a ``--tpu_spec`` and no explicit ``--gpu_spec``.
In that case, ``caliban cloud`` will NOT attach a default GPU to your master
instance. You have to ask for it explicitly.

A CPU mode default also means that by default Caliban will try to install the
``'cpu'`` extra dependency set in your ``setup.py``\ , as described in the
:doc:`../explore/declaring_requirements` guide.

Authorizing TPU Access
~~~~~~~~~~~~~~~~~~~~~~

Before you can pass ``--tpu_spec`` to a job you'll need to authorize your Cloud
TPU to access your service account. Check out `the AI Platform TPU tutorial
<https://cloud.google.com/ml-engine/docs/tensorflow/using-tpus#authorize-tpu>`_
for detailed steps on how to achieve this.

Example Workflows
~~~~~~~~~~~~~~~~~

Next you'll need to get the repository of TPU examples on your machine.

.. code-block:: bash

   mkdir tpu-demos && cd tpu-demos
   curl https://codeload.github.com/tensorflow/tpu/tar.gz/r1.14 -o r1.14.tar.gz
   tar -xzvf r1.14.tar.gz && rm r1.14.tar.gz

Check out the
`AI Platform TPU tutorial <https://cloud.google.com/ml-engine/docs/tensorflow/using-tpus#authorize-tpu>`_
for the next steps, and check back for more detail about how to use that
tutorial with Caliban.
Job Labels
^^^^^^^^^^

AI Platform provides you with the ability to label your jobs with key-value
pairs. Any arguments you provide using either :doc:`custom script arguments
<../explore/custom_script_args>` or an :doc:`experiment broadcast
<../explore/experiment_broadcasting>` will be added to your job as labels, like
this:

In addition to arguments Caliban will add these labels to each job:


* **job_name**: ``caliban_totoro`` by default, or the argument you pass
  using ``caliban cloud --name custom_name``
* **gpu_enabled**\ : ``true`` by default, or ``false`` if you ran your job with
  ``--nogpu``

Cloud has fairly strict requirements on the format of each label's key and
value; Caliban will transform your arguments into labels with the proper
formatting, so you don't have to think about these.

Additional Custom Labels
~~~~~~~~~~~~~~~~~~~~~~~~

You can also pass extra custom labels using ``-l`` or ``--label``\ :

.. code-block:: bash

   caliban cloud -l key:value --label another_k:my_value ...

These labels will be applied to every job if you're running an :doc:`experiment
broadcast <../explore/experiment_broadcasting>`, or to the single job you're
submitting otherwise.

If you provide a label that conflicts with a user argument or experiment flag,
your label will get knocked out.

.. NOTE:: periods aren't allowed in labels, but are often quite meaningful;
   because of this caliban replaces periods with underscores before stripping
   out any restricted characters.
Rate Limiting
^^^^^^^^^^^^^

``caliban cloud`` relies on AI Platform for rate limiting, so you can submit many,
many jobs using an ``--experiment_config`` (up to ~1500 total, I believe?) and AI
Platform will throttle submissions to the default limit of 60 submissions per
minute. If your project's been granted higher quotas, you won't be throttled
until you hit your project's rate limit.

Job submission on Cloud presents a nice progress bar, with terminal colors and
more. The log commands, URLs, jobIds and custom arguments are highlighted so
it's clear which jobs are going through. On a failure the error message prints
in red.
Creating a Bucket
^^^^^^^^^^^^^^^^^

If you need to store data that you generate during a :doc:`../cli/caliban_cloud`
run, storing data in a Cloud bucket is the easiest choice.

Your bucket is a reserved "folder" on the Cloud filesystem; you'll use this to
save models and measurements, and as a staging ground for model workflows you're
submitting to Cloud.

To create your bucket, add the following lines to your ``~/.bashrc`` file:

.. code-block:: bash

   export BUCKET_NAME="totoro_bucket"
   export REGION="us-central1"

Run ``source ~/.bashrc`` to pick up the changes, then run the following command
to create your new bucket:

.. code-block:: bash

   gsutil mb -l $REGION gs://$BUCKET_NAME

That's it.
Creating a Service Account Key
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This page describes how to generate and install a `Service Account Key
<https://www.google.com/search?q=service+account+key+google&oq=service+account+key+google&aqs=chrome..69i57j69i60l2.1592j0j4&sourceid=chrome&ie=UTF-8>`_.
A service account key is a sort of "passport" that your code can use to
authenticate itself during communication with Google's Cloud services.

You can also provide Caliban with a service account key via the ``--cloud_key``
flag. If you do, Caliban will use this service account to authenticate itself
with AI Platform when submitting jobs. (You would do this if you wanted to
submit to some project you didn't own, for example.)

To create a service account key, visit the `Service Accounts page
<https://console.cloud.google.com/iam-admin/serviceaccounts?_ga=2.94132893.1698699355.1592403366-805054138.1592403366>`_
and select the project you created earlier.

Click "Create Service Account" at the top of the page:

.. image:: /_static/img/cloud/activate.png
  :width: 600
  :align: center
  :alt: Activate Billing

At the next form, under **"Service Account Name"**, type something like
**totoro_key** and click **"Create"**.

This will bring up a page titled **"Service Account Permissions"**. Select
**Project > Owner** from the list:

.. image:: /_static/img/cloud/service_acct_permissions.png
  :width: 600
  :align: center
  :alt: Service Account Permissions

Then click **"Continue"** and **"Done"**. You now have a service account. You'll
need to download it to your machine for Caliban to use it.

Downloading the Service Account Key
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Click on the hyperlinked name of the key - something like
``totoro-key@totoro-lives.iam.gserviceaccount.com`` - in the service accounts
list.

Near the bottom of the page, click "Add Key" > "Create New Key":

.. image:: /_static/img/cloud/create_new_key.png
  :width: 600
  :align: center
  :alt: Create New Key

Select **"JSON"** for key type and click **"Create"**. This will download a file
with a name like ``totoro-lives-3df07b8c97a0.json`` to your machine.

Find the file in your terminal (probably in your Downloads folder) and run the
following command to move it to a nice, easy to read location:

.. code-block:: bash

   mv [NEW_FILENAME].json ~/.config/service_key.json

To make this key accessible to Caliban, you'll need to set a variable called
``GOOGLE_APPLICATION_CREDENTIALS`` in your shell to the path of your new service
account key. Add the following line to your `~/.bashrc`:

.. code-block:: bash

   export GOOGLE_APPLICATION_CREDENTIALS=$HOME/.config/service_key.json

If Caliban sees this environment variable set, it will go ahead and bake these
credentials into your container, making them accessible to your code even inside
the Docker environment.
Application Default Credentials
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of a service account key, you might also generate "Application Default
Credentials" on your machine.

To install these on your workstation, run

.. code-block:: bash

   gcloud auth application-default login

at your terminal, as described in `these gcloud docs
<https://cloud.google.com/sdk/gcloud/reference/auth/application-default/login>`_.
That's it!
Customizing Machines and GPUs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section discusses the default configurations for accelerators and machine
types that Caliban requests when it submits jobs to Cloud. You'll also find
instructions on how to request different GPUs or machine types for your job.

Default GPU and Machine Types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, if you don't supply ``--gpu_spec`` or ``--machine_type`` (both discussed
below), Caliban will configure your jobs on the following hardware for each
mode:


* GPU mode (default): a single P100 GPU on an ``n1-standard-8`` machine
* CPU mode: an ``n1-highcpu-32`` machine with no GPU attached

You can read more about the various machine types available on AI platform `here
<https://cloud.google.com/ml-engine/docs/machine-types>`_\ , or scan the
following sections.


Custom GPU Specs
~~~~~~~~~~~~~~~~

The optional ``--gpu_spec`` argument allows you to attach a custom number and type
of GPU to the Cloud node that will run your containerized job on AI Platform.
The required format is ``GPU_COUNTxGPU_TYPE``\ , as in this example:

.. code-block:: bash

   caliban cloud --gpu_spec 2xV100 trainer.train

This will submit your job to a node configured with 2 V100 GPUs to a machine in
the region you specify via:


* your ``$REGION`` environment variable,
* the ``--region`` CLI argument
* or, in the absence of either of those, the safe default of ``us-central1``.

When you run any ``caliban cloud`` command, the program will immediately validate
that the combination of GPU count, region, GPU type and machine type are
compatible and error quickly if they're not. If you make the impossible request
for 3 V100 GPUs:

.. code-block:: bash

   caliban cloud --gpu_spec 3xV100 trainer.train

you'll see this error message:

.. code-block::

   caliban cloud: error: argument --gpu_spec: 3 GPUs of type V100 aren't available
   for any machine type. Try one of the following counts: {1, 2, 4, 8}

   For more help, consult this page for valid combinations of GPU count, GPU type
   and machine type: https://cloud.google.com/ml-engine/docs/using-gpus

If you ask for a valid count, but a count that's not possible on the machine
type you specified - 2 V100s on an ``n1-standard-96`` machine, for example:

.. code-block:: bash

   caliban cloud --gpu_spec 2xV100 --machine_type n1-standard-96 trainer.train

You'll see this error:

.. code-block::

   'n1-standard-96' isn't a valid machine type for 2 V100 GPUs.

   Try one of these: ['n1-highcpu-16', 'n1-highmem-16', 'n1-highmem-2',
   'n1-highmem-4', 'n1-highmem-8', 'n1-standard-16', 'n1-standard-4', 'n1-standard-8']

   For more help, consult this page for valid combinations of GPU count, GPU type
   and machine type: https://cloud.google.com/ml-engine/docs/using-gpus

If you know that your combination is correct, but Caliban's internal
compatibility table hasn't been updated to support some new combination, you can
skip all of these validations by providing ``--force`` as an option.

Custom Machine Types
~~~~~~~~~~~~~~~~~~~~

The ``--machine_type`` option allows you to specify a custom node type for the
master node where your containerized job will run. ``caliban cloud --help`` will
show you all available choices.; You can also read about the various machine
types available on AI platform
`here <https://cloud.google.com/ml-engine/docs/machine-types>`_.

As an example, the following command will configure your job to run on an
``n1-highcpu-96`` instance with 8 V100 GPUs attached:

.. code-block:: bash

   caliban cloud --gpu_spec 8xV100 --machine_type n1-highcpu-96 trainer.train

As described above in :ref:`Custom GPU Specs`, ``--machine_type`` works with
``--gpu_spec`` to validate that the combination of GPU count, GPU type and
machine type are all valid, and returns an error immediately if the combination
is invalid.
GKE Concepts
^^^^^^^^^^^^

Caliban makes it easy to create your own GKE Cluster - similar to your own
personal copy of AI Platform - in your Cloud project, and submit jobs to that
cluster. The advantage over AI Platform currently is that you can get more
quota, often 10x what you have available in AI Platform, and many features are
supported in GKE much earlier than they are in AI Platform.

The quota disparity is particularly notable with TPUs. AI Platform currently
only allows 8 TPUs, while a GKE cluster lets you specify 32, 64, etc TPUs for a
given job.

A good collection of GKE documentation can be found
`here <https://cloud.google.com/kubernetes-engine/docs/concepts>`_

Cluster
~~~~~~~

A
`cluster <https://cloud.google.com/kubernetes-engine/docs/concepts/cluster-architecture>`_
is a collection of cloud machines, combining a set of *nodes* that run your
processing jobs, and *control plane* (also referred to as a *cluster master*\ )
that manages these worker nodes and handles scheduling your jobs and creating
worker nodes to run them.

Cluster Master
~~~~~~~~~~~~~~

A
`cluster master <https://cloud.google.com/kubernetes-engine/docs/concepts/cluster-architecture#master>`_
is the controller for the cluster and all its resources. It handles creating and
deleting worker nodes, and scheduling jobs submitted by users.

Nodes
~~~~~

A
`node <https://cloud.google.com/kubernetes-engine/docs/concepts/cluster-architecture#nodes>`_
is a worker machine (a cloud compute engine instance) that actually performs the
work your job requires. The cluster control plane creates and manages these
instances.

Node Pool
~~~~~~~~~

A
`node pool <https://cloud.google.com/kubernetes-engine/docs/concepts/node-pools>`_
is a collection of identical nodes (cpu, memory, gpu, tpu).

Job
~~~

A
`job <https://cloud.google.com/kubernetes-engine/docs/concepts/batch-reference#batchjobs>`_
is a task that is to be run to completion using cluster resources. The cluster
control plane manages the resources the job needs and handles restarting the job
in case of failure or preemption. A job probably matches the concept you have in
mind when you think of a job you submit to AI platform. A job is a top-level
task, which may be run on multiple machines/containers, which in GKE are
referred to as *pods*\ , described below.

Pod
~~~

A `pod <https://cloud.google.com/kubernetes-engine/docs/concepts/pod>`_ is a
single, ephemeral, running execution of your container. A job may run on several
pods.
GKE Prerequisites
^^^^^^^^^^^^^^^^^

There are a few prerequisites for creating and submitting jobs to a gke cluster.

Required Permissions
~~~~~~~~~~~~~~~~~~~~

To create and use a GKE cluster, you'll need to modify your service account key
to give it Account Owner permissions. Those instructions live at the
:doc:`/cloud/service_account` docs page. Note that this only applies if you are
using a service account key.
Cluster Management
^^^^^^^^^^^^^^^^^^

This section describes how to create and delete clusters. We'll add
documentation on other relevant cluster lifecycle tasks as we go.

Cluster Creation
~~~~~~~~~~~~~~~~

As described in the ``create`` section of :doc:`../cli/caliban_cluster`, you
will typically create a cluster once for a given project and leave it running.

You can create a cluster for your project as follows:

.. code-block:: bash

   totoro@totoro:$ caliban cluster create --cluster_name cluster_name --zone us-central1-a
   I0204 09:24:08.710866 139910209476416 cli.py:165] creating cluster cluster_name in project totoro-project in us-central1-a...
   I0204 09:24:08.711183 139910209476416 cli.py:166] please be patient, this may take several minutes
   I0204 09:24:08.711309 139910209476416 cli.py:167] visit https://console.cloud.google.com/kubernetes/clusters/details/us-central1-a/cluster_name?project=totoro-project to monitor cluster creation progress
   I0204 09:28:05.274621 139910209476416 cluster.py:1091] created cluster cluster_name successfully
   I0204 09:28:05.274888 139910209476416 cluster.py:1092] applying nvidia driver daemonset...

The command will typically take several minutes to complete. The command will
provide you with an url you can follow to monitor the creation process. The page
will look something like the following:

.. image:: /_static/img/gke/cluster_create_progress.png
  :width: 600
  :align: center
  :alt: Cluster creation progress

Once your cluster is created and running, you can view and inspect it from the
cloud dashboard from the ``Kuberenetes Engine > Clusters`` menu option:

.. image:: /_static/img/gke/cluster_dashboard.png
  :width: 600
  :align: center
  :alt: Cluster dashboard

Cluster Deletion
~~~~~~~~~~~~~~~~

In most cases you will bring up your cluster and leave it running. The cluster
master does consume resources, however, so if you know that you are not going to
be submitting jobs to your cluster for some length of time, you may want to
delete your cluster to save money. Before doing this, please make sure that all
of your jobs are complete, as deleting the cluster will also kill any running
jobs. Deleting the cluster is very straightforward, simply using the
:doc:`../cli/caliban_cluster` ``delete`` command.
Single Job Submission
^^^^^^^^^^^^^^^^^^^^^

This is a simple walkthrough for gke job submission from caliban.

Pre-submission Cluster Status
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, we have an existing cluster with no jobs currently running. You
can inspect the cluster from the GCP dashboard for your project under the
``Kubernetes Engine > Clusters`` menu.

.. image:: /_static/img/gke/pre_job_submission.png
  :width: 600
  :align: center
  :alt: Pre-submission

Selecting our ``foo`` cluster, we can see more details.

.. image:: /_static/img/gke/pre_job_details.png
  :width: 600
  :align: center
  :alt: Pre-submission details

Here we can see that our cluster has only a single node pool: the default pool
created when we started the cluster. We will submit a job that uses gpu
acceleration, so we will see how the cluster autoscaler will add a new node pool
for our job based on the gpu and machine specs we provide in the job submission.

We can also see here our cluster limits for autoscaling, which are derived from
our zone quota. These limits control how many instances of different accelerator
resources we can get via autoprovisioning. These limits are cluster-wide, so in
this example we can get at most eight K80 gpus, and at most four T4 gpus.

Submit the Job
~~~~~~~~~~~~~~

To submit a job to your cluster, use the ``caliban cluster job submit`` command.
(see :doc:`../cli/caliban_cluster` for additional examples and documentation.)

Here we create our cluster job (some of output elided):

.. code-block:: bash

   totoro@totoro:$ caliban cluster job submit --gpu_spec 1xK80 --name cifar10-test cifar10_resnet_train.sh --
   I0204 11:33:48.564418 139920906995520 core.py:386] Generating Docker image with parameters:
   I0204 11:33:48.565413 139920906995520 core.py:387] {'adc_path': '/usr/local/google/home/totoro/.config/gcloud/application_default_credentials.json',
    'credentials_path': '/usr/local/google/home/totoro/.config/service_keys/totoro_key.json',
    'extra_dirs': None,
    'job_mode': <JobMode.GPU: 2>,
    'package': Package(executable=['/bin/bash'], package_path='.', script_path='cifar10_resnet_train.sh', main_module=None),
    'requirements_path': 'requirements.txt',
    'setup_extras': None}
   I0204 11:33:48.566865 139920906995520 docker.py:497] Running command: docker build --rm -f- /usr/local/google/home/totoro/sw/tensorflow_models
   Sending build context to Docker daemon  1.058GB

   Step 1/15 : FROM gcr.io/blueshift-playground/blueshift:gpu
    ---> 74f198a8ba19

    ...

   6cebf3abed5f: Layer already exists
   latest: digest: sha256:99c759693d78c24d0b6441e70d5b5538541cccaa158142b5896fadebc30b7ab9 size: 6608
   I0204 11:35:12.189604 139920906995520 cli.py:431] submitted job:
   cifar10-test-tsnlf:
   https://console.cloud.google.com/kubernetes/job/us-central1-a/foo/default/cifar10-test-tsnlf

Our job has now been submitted to our cluster. Due to various factors, it will
take a short time before the job is actually running. We can use the link
provided by caliban to monitor the life cycle of our job.

Monitor Autoscaling/Job Placement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When we first submit the job, we will often see that the job shows what appears
to be an error with a big, ugly, red message saying something along the lines of
"unschedulable".

.. image:: /_static/img/gke/unschedulable.png
  :width: 600
  :align: center
  :alt: Unschedulable

We need to look at the 'details' on the right side to see how the Kubernetes pod
associated with this job is progressing. The job right now is unschedulable
because the cluster has not yet scaled up to accomodate our request. Choosing
the 'details' button, we see this.

.. image:: /_static/img/gke/unschedulable_details.png
  :width: 600
  :align: center
  :alt: Unschedulable details

This is the pod associated with our job. Clicking on this shows us details on
the pod, where we can watch its development. On the pod page, choose the
'Events' tab.

.. image:: /_static/img/gke/pod_events.png
  :width: 600
  :align: center
  :alt: Pod events

Here we can see the progression of the pod. (note that the events here are in
order of 'last seen', so they appear out-of-order when trying to divine the
logical progression of your job) The first event indicates that initially the
cluster does not have any resources to support the pod. The second event shows
that the cluster is scaling up to accomodate this job. This is often the crucial
step. The next relevant event (3) shows that our docker image is being pulled
for our new container. This is then followed by (4) container creation, and then
(5) container start. At this point our job is up and running. Note from the
timestamps that this process took (in this case) approximately ten minutes from
submission to container operation.

While this process is progressing, we can also monitor the cluster and its node
pools from the cluster page:

.. image:: /_static/img/gke/node_pool_autoprovision.png
  :width: 600
  :align: center
  :alt: Node pool autoprovisioning

Now we can see that the cluster has auto-provisioned a new node pool for us in
response to our job submission. Exploring this further you can find the new node
instance that was created and inspect its properties. Once your job has
completed, and if there are no more jobs pending, the cluster will scale down,
deleting the compute node and deleting the node pool.

Monitor Job Logs
~~~~~~~~~~~~~~~~

Now that our job is running, we can monitor the logs from the container from the
dashboard using stackdriver (Kubernetes Engine > Workloads > our-job):

.. image:: /_static/img/gke/job_logs.png
  :width: 600
  :align: center
  :alt: Job logs

This will take you to the stackdriver log viewer for the container:

.. image:: /_static/img/gke/stackdriver_logs.png
  :width: 600
  :align: center
  :alt: Stackdriver logs

Clean up Job
~~~~~~~~~~~~

Once our job has finished, its logs and other data will persist until we delete
it, even though the container has been stopped and no compute resources are
still active. This is quite useful of course, but at some point you will want to
delete the job (which will delete all of the logs and associated metadata, so
use caution)

.. image:: /_static/img/gke/cleanup_job.png
  :width: 600
  :align: center
  :alt: Cleanup job
GCloud and GSUtil Authentication
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Caliban supports authentication with GCloud and GSUtil via two methods:


* `Service Account Keys <https://cloud.google.com/iam/docs/creating-managing-service-account-keys>`_\ ,
  and
* `Application Default Credentials <https://cloud.google.com/sdk/gcloud/reference/auth/application-default/login>`_

Service accounts keys (described in :doc:`../getting_started/cloud`) are the
method of authentication you'll find recommended by most Cloud documentation for
authentication within Docker containers.

 You might also come across a different method of authentication called
"Application Default Credentials", or ADC Credentials. See :doc:`../cloud/adc`
for more information.

.. NOTE:: to set up service account keys, visit the :doc:`service
   account instructions </cloud/service_account>`. To generate application default
   credentials on your machine, simply run ``gcloud auth application-default
   login`` at your terminal, as described `in the Google Cloud docs
   <https://cloud.google.com/sdk/gcloud/reference/auth/application-default/login>`_.

If you've logged in to ``gcloud`` on your machine using application default
credentials, Caliban will copy your stored ADC credentials into your container.
If you DON'T have a service account, gcloud and the cloud python SDK will use
these ADC credentials inside the container and work just as they do on your
workstation.

If you've followed the service account key instructions above and declared a
``GOOGLE_APPLICATION_CREDENTIALS`` environment variable on your system pointing to
a Cloud JSON service account key, Caliban will copy that key into the container
that it builds and set up an environment variable in the container pointing to
the key copy.

You can set or override this variable for a specific caliban command by
supplying ``--cloud_key ~/path/to/my_key.json``\ , like so:

.. code-block:: bash

   caliban run --cloud_key ~/path/to/my_key.json trainer.train

.. WARNING:: If you supply this option to ``caliban shell`` or ``caliban
   notebook`` and have ``GOOGLE_APPLICATION_CREDENTIALS`` set in your
   ``.bashrc``, that variable will overwrite the key that the ``--cloud_key``
   option pushes into your container. To get around this, pass ``--bare`` to
   ``caliban shell`` or ``caliban notebook`` to prevent your home directory from
   mounting and, by extension, any of your environment variables from
   overwriting the environment variable set inside the container.

The environment variable and/or option aren't necessary, but if you don't have
either of them AND you don't have ADC credentials on your machine, you won't be
able to use the GCloud Python API or the ``gsutil`` or ``gcloud`` commands inside
the container.

As noted above, if you don't have this variable set up yet and want to get it
working, check out the :doc:`service account instructions
</cloud/service_account>`. To generate application default credentials on your
machine, simply run ``gcloud auth application-default login`` at your terminal,
as described `in the Cloud docs
<https://cloud.google.com/sdk/gcloud/reference/auth/application-default/login>`_.

GCloud SDK
~~~~~~~~~~

The `GCloud SDK <https://cloud.google.com/sdk/>`_ (\ ``gsutil``\ , ``gcloud`` and friends)
is also available inside of the containerized environment.

On your local machine, ``gsutil`` and ``gcloud`` are authorized using your Google
credentials and have full administrative access to anything in your project.
Inside of the container, these tools are authenticated using the JSON service
account key; this means that if your service account key is missing permissions,
you may see a mismatch in behavior inside the container vs on your workstation.

Shell Mode Caveats
~~~~~~~~~~~~~~~~~~

``caliban shell`` introduces one potentially confusing behavior with these Cloud
credentials. By default, ``caliban shell`` will mount your home directory inside
the container; it does this so that you have all of your bash aliases and your
familiar environment inside of the container. (You can disable this with the
``--bare`` option by running ``caliban shell --bare``\ ).

Mounting your ``$HOME`` directory will trigger an evaluation of your
``$HOME/.bashrc`` file, which will ``export GOOGLE_APPLICATION_CREDENTIALS`` and
overwrite the service key variable that Caliban has set up inside of the
container.

If you use a relative path for this variable on your workstation, like:

.. code-block:: bash

   export GOOGLE_APPLICATION_CREDENTIALS="$HOME/.config/devkey.json"

then everything will still work out wonderfully; inside of the container,
``$HOME`` will resolve to the in-container ``$HOME``\ , but because everything on your
workstation's ``$HOME`` is mounted the container environment will find the key.

If, instead, you use an absolute path, like:

.. code-block:: bash

   export GOOGLE_APPLICATION_CREDENTIALS="/usr/local/google/home/totoro/.config/devkey.json"

The key won't resolve inside the container. (This only applies in ``caliban
shell`` and ``caliban notebook``\ , not in ``caliban {cloud,run}``.)

To fix this, just change your absolute path to a relative path and everything
will work as expected:

.. code-block:: bash

   export GOOGLE_APPLICATION_CREDENTIALS="$HOME/.config/devkey.json"
Custom Docker Run Arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^

``caliban {shell, notebook, run}`` all perform some combination of ``docker build``
and ``docker run`` to provide their functionality. Each provides various sane
defaults that should be fine for most use cases; sometimes, however, you might
need to break through the ``caliban`` abstraction layer and pass arguments to
``docker run`` directly.

One example would be if you need to set environment variables inside the
container, or limit which GPUs are mounted into the container.

To pass custom options to ``docker run``\ , use ``--docker_run_args``\ , like this:

.. code-block:: bash

   caliban run --docker_run_args "--env MY_VARIABLE" trainer.train

This particular command will set ``MY_VARIABLE`` inside the container to its
current value in the shell where you run the above command, as described in the
`docker run <https://docs.docker.com/engine/reference/commandline/run/>`_
documentation. (The
`\ ``docker run`` <https://docs.docker.com/engine/reference/commandline/run/>`_ docs
have information on all possible options.)

This argument is available in ``caliban run``\ , ``caliban shell`` and ``caliban
notebook``.

You may see an error if you pass some flag or argument that ``caliban`` already
supplies. Caliban prints the ``docker run`` command it executes on each
invocation, so if you need full control you can always use ``docker run``
directly.
What's the Base Docker Image?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Caliban's modes build docker images using a dynamically generated ``Dockerfile``.
You'll see this ``Dockerfile`` stream to stdout when you run any of Caliban's
commands. You can specify the base image you wish to use in your ``.calibanconfig.json``
file as described :ref:`here<calibanconfig>`.

In addition to the isolation Docker provides, the images set up a Python virtual
environment inside of each container. This guarantees you a truly blank slate;
the dependencies you declare in your code directory are the only Python
libraries that will be present. No more version clashes or surprises.

Caliban uses a set of base images covering a set of common combinations of
python and cuda versions. You can find our base images
`here <https://gcr.io/blueshift-playground/blueshift>`_.
The format of our base image names is ``gcr.io/blueshift-playground/blueshift:TAG``,
where ``TAG`` describes the configuration of the base image.

For example, ``gcr.io/blueshift-playground/blueshift:gpu-ubuntu1804-py38-cuda101`` is a
gpu base image that uses Ubuntu 18.04, CUDA 10.1 and python 3.8, while
``gcr.io/blueshift-playground/blueshift:cpu-ubuntu2004-py38`` is a cpu-only Ubuntu 20.04
base image that has no CUDA support and uses python 3.8.

Our current supported combinations:

+--------------+------------+------------+
| Ubuntu 18.04 | python 3.7 | python 3.8 |
+--------------+------------+------------+
|   no cuda    |    yes     |    yes     |
+--------------+------------+------------+
|   cuda 10.0  |    yes     |    yes     |
+--------------+------------+------------+
|   cuda 10.1  |    yes     |    yes     |
+--------------+------------+------------+


+--------------+------------+------------+
| Ubuntu 20.04 | python 3.7 | python 3.8 |
+--------------+------------+------------+
|   no cuda    |    yes     |    yes     |
+--------------+------------+------------+
|   cuda 10.0  |    no      |     no     |
+--------------+------------+------------+
|   cuda 10.1  |    no      |     no     |
+--------------+------------+------------+


These images are automatically updated, and if you have an image combination that
we don't support, please file an issue and we'll consider adding it to our set
of supported images. We are planning to add support for custom base images so
you can build and use your own specialized image.

The dockerfiles we use to generate our supported images can be found
`here <https://github.com/google/caliban/tree/master/dockerfiles>`_. We create
base gpu images from the `Dockerfile.gpu <https://github.com/google/caliban/blob/master/dockerfiles/Dockerfile.gpu>`_
file, and then use these as base images for creating full GPU images with
support for specific python versions using this `Dockerfile <https://github.com/google/caliban/blob/master/dockerfiles/Dockerfile>`_.

We base our gpu base images on the `nvidia/cuda <https://hub.docker.com/r/nvidia/cuda/>`_
images, which contain the relevant CUDA drivers required for GPU use. The virtual
environment inside of the Caliban container isolates you from these low-level details,
so you can install any tensorflow version you like, or use Jax or Pytorch or any
other system.

Details for Maintainers
~~~~~~~~~~~~~~~~~~~~~~~

We utilize Google's `Cloud Build <http://cloud.google.com/cloud-build/docs>`_ service
to build Caliban's base images. Our Cloud Build configuration file that controls
our image generation can be found in the source repository
`here <https://github.com/google/caliban/blob/master/cloudbuild.json>`_.

This file can quickly get lengthy and difficult to maintain, so we generate this file
using `a script <https://github.com/google/caliban/blob/master/scripts/cloudbuild.py>`_
and `a configuration file <https://github.com/google/caliban/blob/master/scripts/cloudbuild_config.json>`_.
In the configuration file, we specify our supported CUDA versions, our supported
python versions, and a list of the combinations we use in our supported images.
For our CUDA and python versions, we specify a list of build-args that we then
pass to the docker build process for the Dockerfiles described above.

To generate a new ``cloudbuild.json`` file, invoke the ``cloudbuild.py`` utility with
your configuration file:

.. code-block:: bash

   python ./scripts/cloudbuild.py --config scripts/cloudbuild_config.json  --output cloudbuild.json

This will generate a new ``cloudbuild.json`` file which is used by the Cloud Build service
to generate our base docker images. For testing, you can set a different base image url for
the docker images by using the ``--base_url`` keyword argument.

To manually start a Cloud Build for these docker images, navigate to the top-level
of the caliban source repository, and use the
`gcloud builds <https://cloud.google.com/cloud-build/docs/running-builds/start-build-manually#gcloud>`_
command:

.. code-block:: bash

   gcloud builds submit --project=<destination project> --config=cloudbuild.json .

By default this uses your default project and the ``cloudbuild.json`` file in your current
directory. If you are pushing the images to a different project than your ``gcloud`` default,
then you may need to set the ``--project`` flag to the target project where you are pushing
your images. The logs from the build process will be streamed to your console, but they are
also available from the ``Cloud Build`` tab in the GCP dashboard for your project.

To automate the generation of these images, we utilize
`build triggers <https://cloud.google.com/cloud-build/docs/automating-builds/create-manage-triggers>`_
to start a new cloud build whenever the Caliban Dockerfiles are modified in the repository.
calibanconfig
^^^^^^^^^^^^^^^^^^^^^^

Caliban supports customization through a file called ``.calibanconfig.json``
that lives in your project's directory. Features are limited for now, but stay
tuned for more.

Custom Apt Packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Caliban provides support for custom aptitude packages inside your container. To
require custom apt packages, create a file called ``.calibanconfig.json`` inside
your project's directory.

The ``.calibanconfig.json`` should contain a single JSON dictionary with an
``"apt_packages"`` key. The value under this key can be either a list, or a
dictionary with ``"gpu"`` and ``"cpu"'`` keys. For example, any of the following are
valid:

.. code-block::

   # This is a list by itself. Comments are fine, by the way.
   {
        "apt_packages": ["libsm6", "libxext6", "libxrender-dev"]
   }

This works too:

.. code-block:: json

   # You can also include a dictionary with different deps
   # for gpu and cpu modes. It's fine to leave either of these blank,
   # or not include it.
   {
       "apt_packages": {
           "gpu": ["libsm6", "libxext6", "libxrender-dev"],
           "cpu": ["some_other_package"]
       }
   }

These values will do what you expect and run ``apt-get install <package_name>``
for each package. Packages are alphabetized, so changing the order won't
invalidate Docker's build cache.

Custom Base Images
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For details on Caliban's base images, see :ref:`What's the Base Docker Image?`.

You can specify a custom base image for Caliban to use in your ``.calibanconfig.json`` file
by adding an entry with the ``base_image`` key as follows:

.. code-block:: json

   {
       "base_image": "gcr.io/blueshift-playground/blueshift:gpu-ubuntu1804-py38-cuda101"
   }

You can also specify different base images for ``cpu`` and ``gpu`` modes as follows:

.. code-block:: json

   {
       "base_image": {
           "cpu": "gcr.io/blueshift-playground/blueshift:cpu-ubuntu1804-py38",
           "gpu": "gcr.io/blueshift-playground/blueshift:gpu-ubuntu1804-py38-cuda101"
       }
   }
Custom Script Arguments
^^^^^^^^^^^^^^^^^^^^^^^

In ``caliban run`` or ``caliban cloud`` modes, if you pass ``--`` to the CLI, Caliban
will stop parsing commands and pass everything after ``--`` through to your
script, untouched. If you run:

.. code-block:: bash

   caliban cloud trainer.train -- --epochs 2 --job_dir my_directory

Your script will execute inside the container environment with the following
command:

.. code-block:: bash

   python -m trainer.train --epochs 2 --job_dir my_directory

This feature is compatible with :doc:`experiment_broadcasting` in ``cloud``,
``run`` or ``cluster`` mode; arguments are prepended to the list generated by
the specific experiment being executed from your experiment config.
Experiment Groups
^^^^^^^^^^^^^^^^^

Caliban supports grouping experiments into a collection called an *experiment
group*. This allows you to do things like monitor all of the jobs in a given
group, stop all running jobs in a group, or re-run all of the jobs in a group.

Each of the caliban compute backends supports specifying an experiment group via
the ``--xgroup`` flag:

.. code-block::

   $ caliban run --xgroup my-xgroup ...
   $ caliban cloud --xgroup my-xgroup ...
   $ caliban cluster job submit --xgroup my-xgroup ...

If you don't specify an experiment group when submitting jobs via caliban, a new
experiment group will be generated for you, so you don't need to use them if you
don't want to. Also, the existence of this group should be transparent to you.

You can add new jobs to an existing experiment group simply by specifying the
same group on different caliban job submission calls:

.. code-block::

   caliban cloud --xgroup my-xgroup ... foo.py --
   ...
   (some time later...)
   caliban cloud --xgroup my-xgroup ... bar.py --

The experiment group ``my-xgroup`` will contain the experiments generated by both
of the caliban calls, and you can then perform different operations on these as
described in the sections below.
What can Caliban Execute?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Caliban's commands can run python files as modules or scripts. If you need more
customization, you can run arbitrary shell scripts with Caliban.

Script vs Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Inside the containerized environment, your Python script will run as a module or
a script, depending on the format of the argument you supply to caliban. If you
explicitly pass a python module, with components separated by dots:

.. code-block:: bash

   caliban cloud trainer.train -- --epochs 2 --job_dir my_directory

Your script will execute inside the container environment with the following
command:

.. code-block:: bash

   python -m trainer.train --epochs 2 --job_dir my_directory

If instead you supply a relative path to the python file, like this:

.. code-block:: bash

   caliban cloud trainer/train.py -- --epochs 2 --job_dir my_directory

Caliban will execute your code as a python *script* by passing it directly to
python without the ``-m`` flag, like this:

.. code-block:: bash

   python trainer/train.py --epochs 2 --job_dir my_directory

What does this mean for you? Concretely it means that if you execute your code
as a module, all imports inside of your script have to be declared relative to
the root directory, ie, the directory where you run the caliban command. If you
have other files inside of the ``trainer`` directory, you'll have to import them
from ``trainer/train.py`` like this:

.. code-block:: python

   import trainer.util
   from trainer.cloud import load_bucket

We do this because it enforces a common structure for all code. The reproducible
unit is the directory that holds all of the code. The script doesn't live in
isolation; it's part of a project, and depends on the other files in the code
tree as well as the dependencies declared in the root directory.

If you run your code as a script, imports will only work if they're relative to
the file itself, not to the running code.

I highly recommend running code as a module!

Using Caliban with Shell Scripts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Caliban can build containers for you that will execute arbitrary shell scripts,
in addition to python code.

If you pass a relative path that points to any file other other than:


* a python module, or
* an explicit path to a python file ending with ``.py``\ ,

to ``caliban cloud``\ , ``caliban run`` or one of the other modes that accepts
modules, caliban will execute the code as a bash script.

This feature is compatible with :doc:`custom script arguments
<custom_script_args>` or an :doc:`experiment broadcast
<experiment_broadcasting>`; your shell script will receive the same flags that
any python module would receive.
Experiment Config via stdin, pipes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to passing an explicit JSON file to ``caliban cloud
--experiment_config``\ , if you pass the string ``stdin`` as the flag's value
``caliban cloud`` will attempt to read the experiment config in off of ``stdin``.

As an example, this command pipes in a config and also passes ``--dry_run`` to
show the series of jobs that WILL be submitted when the ``--dry_run`` flag is
removed:

.. code-block:: bash

   cat experiment.json | caliban cloud --experiment_config stdin --dry_run trainer.train

Because ``experiment.json`` is a file on disk, the above command is not that
interesting, and equivalent to running:

.. code-block:: bash

   caliban cloud --experiment_config experiment.json --dry_run trainer.train

Things get more interesting when you need to dynamically generate an experiment
config.

Imagine you've written some python script ``generate_config.py`` that builds up a
list of complex, interdependent experiments. If you modify that script to print
a ``json`` list of ``json`` dicts when executed, you can pipe the results of the
script directly into ``caliban cloud``\ :

.. code-block:: bash

   python generate_config.py --turing_award 'winning' | \
     caliban cloud --experiment_config stdin --dry_run trainer.train

And see immediately (thanks to ``--dry_run``\ ) the list of jobs that would be
executed on AI Platform with a real run.


Experiment File Expansion and Pipes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :doc:`../cli/expansion` command described :doc:`above <../cli/expansion>`
allows you to expand an experiment config into its component JSON objects.
Because these are printed to ``stdout``\ , you can pipe them directly in to
Caliban's commands, like this:

.. code-block:: bash

   expansion experiment.json | caliban cloud --experiment_config stdin trainer.train

You can also insert your own script into the middle of this pipeline. Imagine a
script called ``my_script.py`` that:


* reads a JSON list of experiments in via ``stdin``
* modifies each entry by inserting a new key whose value is a function of one
  or more existing entries
* prints the resulting JSON list back out to ``stdout``

You could sequence these steps together like so:

.. code-block:: bash

   cat experiment.json | \
     expansion experiment.json | \
     my_script.py | \
     caliban cloud --experiment_config stdin --dry_run trainer.train

If you supply ``--dry_run`` to caliban, as in the example above, caliban will
print out all of the jobs that this particular command will kick off when you
remove ``--dry_run``. This is a great way to generate complex experiments and test
everything out before submitting your jobs.
Why Caliban and Docker?
^^^^^^^^^^^^^^^^^^^^^^^

Caliban uses Docker to build isolated environments for your research code. What
does this mean, and why would you want to do this?

One major source of friction in machine learning research is the potential
mismatch between the environment where your code runs during local development
and the environment in AI Platform or Cloud. Here's a typical situation:


* You run your code locally against some set of dependencies you installed
  months ago in the virtual environment you use for all your code.
* You get everything working and submit it to Cloud. Minutes later you see a
  failure - your specified Tensorflow version is wrong. You submit again,
  specifying the beta of TF 2.0 that you've been using... and the job fails.
  That version's not available in Cloud.
* Finally the submission works, but the job fails again. The ``gsutil`` command
  you've been shelling out to to save your models locally isn't available on
  AI Platform.
* You sigh and look at the clock. It's 4pm. Should I have another cup of
  coffee? What am I even doing? Is this what my life has become?

Each of these issues is small, but they stack up and turn you into a broken,
cautious person, afraid to flex the wings you've forgotten are attached to your
back.

Docker is the answer to this problem. `Docker <https://www.docker.com/>`_ is a
piece of software that allows you to build and run "containers"; you can think
of a container as a tiny Linux machine that you can run on your Mac or
workstation, or ship off to execute on AI platform. The container gets access to
the resources of the machine where it's running, but can't affect that machine
in any other way.

If you design your Python code to run inside of a container, you can move that
container between different environments and know that the code's behavior won't
change.

The Trouble with Bare Docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To build a Docker container for your code you need to write a ``Dockerfile``. If
you try this you'll realize that you actually need many ``Dockerfile`` copies...
one for GPU mode. One for CPU mode locally. Slight tweaks show up every time you
want to add some environment variable; locally, you don't want to copy your code
into the container, since you can live-mount the directory using ``docker run``\ ,
but on AI Platform you DO need a copy.

Soon your ``Dockerfile`` is infested with comments and instructions to a future,
less patient version of yourself, even less capable of remembering all of this
than you are now.

Caliban + Docker = <3
~~~~~~~~~~~~~~~~~~~~~

If you've felt this pain, you now understand the motivation for Caliban. Caliban
is a tool that dynamically builds docker images (by dynamically generating
``Dockerfile`` instances) for the various modes you rely on for machine learning
research:


* Jupyter notebook development
* Local, interactive development at the shell
* Local execution on your workstation on GPU
* AI platform execution of 100s of jobs for some experiment

By developing your research workflows inside of Docker containers (made easy by
Caliban) you're much closer to that noble goal of reproducible research.

Theoretically, you could publish the container that Caliban builds along with
the range of experiment parameters you used to produce your data.
Experiment Broadcasting
^^^^^^^^^^^^^^^^^^^^^^^

The ``--experiment_config`` keyword argument allows you to pass Caliban a config
that can run many instances of your containerized job by passing each job a
different combination of some set of parameters. These parameters are passed to
your job as ``--key value`` style flags that you can parse with
`\ ``abseil`` <https://abseil.io/docs/python/quickstart>`_ or
`\ ``argparse`` <https://docs.python.org/3/library/argparse.html>`_.

This keyword is accepted by the following subcommands:


* ``caliban cloud``\ , to submit experiments to AI Platform
* ``caliban run`` to run experiments in sequence on a local workstation
* ``caliban cluster`` to execute experiments on a GKE cluster

The documentation below will refer to ``caliban cloud``\ , but all commands will
work just as well with these other modes unless explicitly called out otherwise.

``--experiment_config`` accepts a path, local or absolute, to a JSON file on your
local machine. That JSON file defines a sweep of parameters that you'd like to
explore in an experiment. Let's look at the format, and what it means for job
submission.

Experiment.json Format
~~~~~~~~~~~~~~~~~~~~~~

You can name the file whatever you like, but we'll refer to it here as
``experiment.json`` always. Here's an example ``experiment.json`` file:

.. code-block:: json

   {
       # comments work inside the JSON file!
       "epochs": [2, 3],
       "batch_size": [64, 128], # end of line comments too.
       "constant_arg": "something"
       "important_toggle": [true, false]
   }

The following command will submit an experiment using the above experiment
definition:

.. code-block:: bash

   caliban cloud --experiment_config ~/path/to/experiment.json trainer.train

For this particular ``experiment.json`` file, Caliban will submit 8 different jobs
to AI Platform with the following combinations of flags, one combination for
each job:

.. code-block:: bash

   --epochs 2 --batch_size 64 --constant_arg 'something' --important_toggle
   --epochs 2 --batch_size 64 --constant_arg 'something'
   --epochs 2 --batch_size 128 --constant_arg 'something' --important_toggle
   --epochs 2 --batch_size 128 --constant_arg 'something'
   --epochs 3 --batch_size 64 --constant_arg 'something' --important_toggle
   --epochs 3 --batch_size 64 --constant_arg 'something'
   --epochs 3 --batch_size 128 --constant_arg 'something' --important_toggle
   --epochs 3 --batch_size 128 --constant_arg 'something'

As you can see, keys get expanded out into ``--key`` style flags by prepending a
``--`` onto the key string. Here are the rules for value expansion:


* ``int`` and ``string`` values are passed on to every job untouched.
* lists generate multiple jobs. ``caliban cloud`` takes the cartesian product of
  all list-type values and generates a job for each combination. Three lists
  of length 2 in the above example gives us 8 total jobs; one for each
  possible combination of items from each list.
* if a value equals ``true``\ , the key is passed through as ``--key``\ , with no
  value; it's treated as a boolean flag.
* a ``false`` boolean value means that the ``--key`` flag is ignored.

All arguments generated from the experiment config will create labels in the AI
Platform Job UI for each job as described in the :doc:`../cloud/labels` section.

Any :doc:`custom script arguments <custom_script_args>` you pass after the
module name, separated by ``--``\ , will be passed along to every job as if they
were static key-value pairs in the ``experiment.json`` file. As an example, the
following command:

.. code-block:: bash

   caliban cloud --experiment_config ~/path/to/experiment.json trainer.train -- --key value

would trigger the same jobs as before, with ``--key value`` appended BEFORE the
arguments broadcast out by the experiment config:

.. code-block:: bash

   --key value --epochs 2 --batch_size 64 --constant_arg 'something' --important_toggle
   --key value --epochs 2 --batch_size 64 --constant_arg 'something'
   # ....etc

Lists of Experiment Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can pass either an experiment config or a LIST of experiment configs in your
``experiment.json`` file; caliban will expand each entry in the list recursively.
This makes it possible to generate experiment configs that aren't strict
cartesian products.

For example you might add the following to a file called ``experiment.json``\ , and
pass it to your job with ``--experiment_config experiment.json``\ :

.. code-block:: json

   [
       {
           "epochs": [1,2, 3, 4],
           "batch_size": [64, 128],
           "constant_arg": "something",
           "important_toggle": [true, false]
       },
       {
           "epochs": [9, 10],
           "batch_size": [512, 1024],
           "constant_arg": "something"
       }
       {
           "epochs": 1000,
           "batch_size": 1
       }
   ]

This config will generate:


* 16 combinations for the first dictionary (every combination of 4 epoch
  entries, 2 batch sizes, and 2 ``"important_toggle"`` combos, with
  ``"constant_arg"`` appended to each)
* 4 combos for the second dictionary
* 1 static combo for the third entry.

for a total of 21 jobs. You can always pass ``--dry_run`` (see below) to ``caliban
cloud`` to see what jobs will be generated for some experiment config, or to
validate that it's well-formed at all.

Compound keys
~~~~~~~~~~~~~

By default, an experiment specification in which multiple values are lists will
be expanded using a Cartesian product, as described above. If you want multiple
arguments to vary in concert, you can use a compound key. For example, the
following (w/o compound keys) experiment config file will result in four jobs
total:

.. code-block:: json

   {
     "a": ["a1", "a2"],
     "b": ["b1", "b2"]
   }

Results in:

.. code-block:: bash

   --a a1 --b b1
   --a a1 --b b2
   --a a2 --b b1
   --a a2 --b b2

To tie the values of ``a`` and ``b`` together, specify them in a compound key:

.. code-block:: json

   {
     "[a,b]": [["a1", "b1"], ["a2", "b2"]]
   }

This will result in only two jobs: ``bash --a a1 --b b1 --a a2 --b b2``

``--dry_run``
~~~~~~~~~~~~~~~~~

Passing an ``--experiment_config`` to ``caliban cloud`` could potentially submit
many, many jobs. To verify that you have no errors and are submitting the number
of jobs you expect, you can add the ``--dry_run`` flag to your command, like this:

.. code-block:: bash

   caliban cloud --dry_run --experiment_config ~/path/to/experiment.json trainer.train

``--dry_run`` will trigger all of the logging side effects you'd see on job
submission, so you can verify that all of your settings are correct. This
command will skip any docker build and push phases, so it will return
immediately with no side effects other than logging.

Once you're sure that your jobs look good and you pass all validations, you can
remove ``--dry_run`` to submit all jobs.

Experiments and Custom Machine + GPUs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you supply a ``--gpu_spec`` or ``--machine_type`` in addition to
``--experiment_config``\ , every job in the experiment submission will be configured
with those options.
Declaring Requirements
^^^^^^^^^^^^^^^^^^^^^^

To use a Python library in your Caliban-based workflow you'll need to declare it
in either a


* ``requirements.txt`` file in the directory, or a
* ``setup.py`` file, or
* both of these together.

If you run any of the Caliban commands in a directory without these, your image
will have access to bare Python alone with no dependencies.

A ``requirements.txt`` file is the simplest way to get started. See the
`pip docs <https://pip.readthedocs.io/en/1.1/requirements.html>`_ for more
information on the structure here. You've got ``git`` inside the container, so
``git`` dependencies will work fine.

Setup.py and Extra Dependency Sets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Declaring your dependencies in a ``setup.py`` file gives you the ability to
declare different sets of dependencies for the different Caliban modes (CPU vs
GPU), in addition to your own custom dependency sets.

This solves the problem of depending on, say, ``tensorflow-gpu`` for a GPU job,
and ``tensorflow`` for normal, CPU-only jobs, without having to modify your
dependency file.

Here's an example ``setup.py`` file:

.. code-block:: python

   from setuptools import find_packages
   from setuptools import setup

   setup(
       name='hello-tensorflow',
       version='0.1',
       install_requires=['absl-py', 'google-cloud-storage'],
       extras_require={
           'cpu': ['tensorflow==2.0.*'],
           'gpu': ['tensorflow-gpu==2.0.*'],
       },
       packages=find_packages(),
       description='Hello Tensorflow setup file.')

This project has two normal dependencies - ``'absl-py'`` for flags, and
``'google-cloud-storage'`` to interact with Cloud buckets.

The ``setup.py`` file declares its Tensorflow dependencies in a dictionary under
the ``extras_require`` key. If you're using pip, you would install dependencies
from just ``install_requires`` by running

.. code-block:: bash

   pip install .

If you instead ran

.. code-block:: bash

   pip install .[gpu]

``pip`` would install


* the entries under ``install_requires``\ ,
* AND, additionally, the entries under the ``'gpu'`` key of the ``extras_require``
  dictionary.

By default, if you have a ``setup.py`` file in your directory, caliban will do the
latter and attempt to install a ``'gpu'`` set of extras, like

.. code-block::

   pip install .[gpu]

If you pass ``--nogpu`` to any of the commands, Caliban will similarly attempt to
run

.. code-block::

   pip install .[cpu]

If you don't declare these keys, don't worry. You'll see a warning that the
extras dependencies didn't exist, and everything will proceed, no problem.

If you have some other set of dependencies you want to install, you can pass
``--extras my_deps``\ , or ``-e my_deps``\ , to any of the caliban modes install those
in addition to the ``cpu`` or ``gpu`` dependency set.

You can provide many sets, like this:

.. code-block:: bash

   caliban cloud -e my_deps -e logging_extras <remaining args>

And Caliban will install the dependencies from all declared sets inside of the
containerized environment.
Caliban on a Mac
^^^^^^^^^^^^^^^^^^^^^^

If you're developing on your Macbook, you'll be able to build GPU containers,
but you won't be able to run them locally. You can still submit GPU jobs to AI
Platform!

To use Caliban's ``shell``\ , ``notebook`` and ``run``\ , you'll have to pass
``--nogpu`` as a keyword argument. If you don't do this you'll see the following
error:

.. code-block:: text

   [totoro@totoro-macbookpro hello-tensorflow (master)]$ caliban run trainer.train

   'caliban run' doesn't support GPU usage on Macs! Please pass --nogpu to use this command.

   (GPU mode is fine for 'caliban cloud' from a Mac; just nothing that runs locally.)

The :doc:`../getting_started/prerequisites` page covers Macbook installation of
Docker and other dependencies.
Mounting a Local Directory for Data Persistence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's say you're using ``caliban run`` with an experiment configuration to run
many experiments locally. Because ``caliban run`` attempts to look just like the
environment you'll see in the Cloud, the command doesn't mount any local
directories by default; the container is completely isolated, and you (usually)
have to persist data by writing it to a Cloud bucket.

It's possible to avoid this, however, and use Caliban to mount a local directory
into the Docker container. If you do this, you can take advantage of local
experiment broadcasting to loop through many experimental runs on your
workstation, and still persist all results and models to your local machine.

The answer comes from the :doc:`../explore/custom_docker_run` feature. If you
pass

.. code-block:: bash

   --docker_run_args "--volume workstation_dir:/foo"

to ``caliban run``\ , Caliban will mount the directory at ``workstation_dir`` into
your container at ``/foo``. (You can use any name or directory you choose instead
of ``/foo``\ , of course.)

Let's look at an example. The following command will mount a folder called
``data`` in your workstation's home directory into your container.

.. code-block:: bash

   caliban run \
     --docker_run_args "--volume /usr/local/google/home/totoro/data:/foo"
     --experiment_config exp_config.json \
     trainer.train

When you look at ``/foo`` inside the container, you'll see all of the files on
your workstation at ``/usr/local/google/home/totoro/data``. If you create or
edit any files, those changes will happen to the files on your workstation as
well.

.. WARNING:: For some reason I don't understand, if you pass ``-v`` instead of
   ``--volume``\ , as in ``--docker_run_args "-v mydir:containerdir"``\ , the
   argument parser in Caliban will break. Use ``--volume`` and you'll be set!

If you want to play around with volume mounting, you can pass the same argument
to ``caliban shell`` to get an interactive view of the filesystem your container
will have access to when you run the above command:

.. code-block:: bash

   # "--bare" prevents your home directory from mounting.
   caliban shell --bare \
   --docker_run_args "--volume /usr/local/google/home/totoro/data:/foo"

In the shell that launches you'll see the directory mirrored:

.. code-block::

   $ caliban shell --docker_run_args "--volume /usr/local/google/home/totoro/data:/foo" --nogpu --bare
   I0122 14:30:24.923780 4445842880 docker.py:438] Running command: docker build --rm -f- /Users/totoro/code/python/tutorials/hello-tensorflow
   Sending build context to Docker daemon  36.56MB
   <....lots of Docker output....>
   Successfully built f2ba6fb7b628
   I0122 14:30:33.125234 4445842880 docker.py:666] Running command: docker run --ipc host -w /usr/app -u 735994:89939 -v /Users/totoro/code/python/tutorials/hello-tensorflow:/usr/app -it --entrypoint /bin/bash --volume /usr/local/google/home/totoro/data:/foo f2ba6fb7b628
      _________    __    ________  ___    _   __  __  __
     / ____/   |  / /   /  _/ __ )/   |  / | / /  \ \ \ \
    / /   / /| | / /    / // __  / /| | /  |/ /    \ \ \ \
   / /___/ ___ |/ /____/ // /_/ / ___ |/ /|  /     / / / /
   \____/_/  |_/_____/___/_____/_/  |_/_/ |_/     /_/ /_/

   You are running caliban shell as user with ID 735994 and group 89939,
   which should map to the ID and group for your user on the Docker host. Great!

   caliban-shell /usr/app > ls -al /foo
   total 9788
   drwx------ 21 totoro 89939     672 Jan 22 20:35  .
   drwxr-xr-x  1 root       root     4096 Jan 22 21:30  ..
   -rw-r--r--  1 totoro 89939   41689 Jan 20 21:48  sets.png
   -rw-r--r--  1 totoro 89939   82811 Jan 20 21:48  tree.png
   caliban-shell /usr/app >
dockerignore speeds up builds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Many of Caliban's commands begin their work by triggering a ``docker build``
command; this command has a side effect of bundling up the entire directory
where you run the command into a "build context", which is zipped up and sent
off to the Docker build process on your machine.

In a directory containing machine learning code, it's not unusual that you might
also have subdirectories that contain, for example:


* large datasets that you've cached locally
* tensorboard output from local runs
* metrics

If you don't want to include any of these things in the Docker container that
caliban builds for you, you can significantly speed up your builds by creating a
file called ``.dockerignore`` in the directory of your project.

Here's an example ``.dockerignore`` file, with comments explaining each line:

.. code-block::

   # ignore the git repository info and the pip installation cache
   .git
   .cache

   # this is huge - ignore the virtualenv we've created inside the folder!
   env

   # tests don't belong inside the repo.
   tests

   # no need to package info about the packaged-up code in egg form.
   *.egg-info

   # These files are here for local development, but have nothing
   # to do with the code itself, and don't belong on the docker image.
   Makefile
   pylintrc
   setup.cfg
   __pycache__
   .coverage
   .pytest_cache

As a starting point, you might take your project's ``.gitignore`` file, copy
everything other to ``.dockerignore`` and then delete any entries that you
actually DO need inside your Docker container. An example might be some data you
don't control with ``git``\ , but that you do want to include in the container using
Caliban's ``-d`` flag.
Using a Single GPU
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, ``docker run`` will make all GPUs on your workstation available
inside of the container. This means that in ``caliban shell``\ , ``caliban
notebook`` or ``caliban run``\ , any jobs executed on your workstation will
attempt to use:


* your huge GPU, custom-built and installed for ML Supremacy
* the dinky GPU that exists solely to power your monitor, NOT to help train
  models

The second GPU will slow down everything.

To stop this from happening you need to set the ``CUDA_VISIBLE_DEVICES``
environment variable equal to ``0``\ , as described on this
`nvidia blog <https://devblogs.nvidia.com/cuda-pro-tip-control-gpu-visibility-cuda_visible_devices/>`_
about the issue.

You can set the environment variable inside your container by passing
``--docker_run_args`` to caliban, like this:

.. code-block:: bash

   caliban run --docker_run_args "--env CUDA_VISIBLE_DEVICES=0" trainer.train

.. NOTE:: you may have noticed that this problem doesn't happen when you run a
   job inside ``caliban shell``. Your local environment may have
   ``CUDA_VISIBLE_DEVICES`` set. ``caliban shell`` and ``caliban notebook``
   mount your home directory by default, which loads all of your local
   environment variables into the container and, if you've set this environment
   variable, modifies this setting inside your container. This doesn't happen
   with ``caliban run`` or ``caliban cloud``. You will always need to use this
   trick with those modes.

There are two other ways to solve this problem using the
`custom ``docker run`` arguments detailed here <https://docs.docker.com/engine/reference/commandline/run/>`_.
You can directly limit the GPUs that mount into the container using the ``--gpus``
argument:

.. code-block:: bash

   caliban run --docker_run_args "--gpus device=0" trainer.train

If you run ``nvidia-smi`` in the container after passing this argument you won't
see more than 1 GPU. This is useful if you know that some library you're using
doesn't respect the ``CUDA_VISIBLE_DEVICES`` environment variable for any reason.

You could also pass this and other environment variables using an env file.
Given some file, say, ``myvars.env``\ , whose contents look like this:

.. code-block:: text

   CUDA_VISIBLE_DEVICES=0
   IS_THIS_A_VARIABLE=yes

The ``--env-file`` argument will load all of the referenced variables into the
docker environment:

.. code-block:: bash

   caliban run --docker_run_args "--env-file myvars.env" trainer.train

Check out :doc:`../explore/custom_docker_run` for more information.
Passing Flags via --flagfile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you find yourself passing lots of flags in to some caliban subcommand, you
might consider Abseil's ``--flagfile`` feature.

.. NOTE:: `Abseil <https://abseil.io/docs/python>`_ is a Google library that we
   use to generate Caliban's CLI. You can see the options `Abseil
   <https://abseil.io/docs/python>`_ provides on top of Caliban's arguments by
   passing ``--helpfull`` to any command; ``caliban cloud --helpfull``\ , for
   example.

``--flagfile`` allows you to put any number of flags or arguments to caliban into
a file, one pair per line. Given some file like ``my_args.txt`` with the following
contents:

.. code-block::

   --docker_run_args "CUDA_VISIBLE_DEVICES=0"
   --experiment_config experiment_one.json
   --cloud_key my_key.json
   --extras extra_deps

You could run the following command:

.. code-block:: bash

   caliban run --flagfile my_args.txt trainer.train

All arguments expand in-line, so the above command would be equivalent to
running:

.. code-block:: bash

   caliban run --docker_run_args "CUDA_VISIBLE_DEVICES=0" \
               --experiment_config experiment_one.json \
               --cloud_key my_key.json \
               --extras extra_deps \
               trainer.train

One major benefit is that you can share groups of arguments between various
subcommand invocations, like ``caliban run`` and ``caliban cloud``\ , without having
to store large duplicated strings of arguments.

Nested Flagfiles
~~~~~~~~~~~~~~~~

You can supply ``--flagfile some_file`` arguments inside flag files! This allows
you to build up trees of arguments in a fine grained way. Imagine some flagfile
called ``v100_project.flags``\ :

.. code-block:: text

   # Definition for big iron GPUs.
   --gpu_spec 8xV100
   --machine_type n1-highcpu-64
   --cloud_key my_key.json

And then some further file called ``tpu_plus_gpu.flags``\ :

.. code-block:: text

   --flagfile v100_project.flags
   --tpu_spec 8xV3
   --region us-central1

The command:

.. code-block:: bash

   caliban cloud --flagfile tpu_plus_gpu.flags trainer.train

Would expand out **both** sets of flags, as expected. (I don't know what would
happen if each file referenced the other... feel free to try!)

For more information, check out the
`Abseil docs on ``--flagfile`` <https://abseil.io/docs/python/guides/flags#a-note-about---flagfile>`_.
Getting Caliban
---------------

.. warning:: If you're currently in a ``virtualenv``\ , please run ``deactivate``
   to disable it before proceeding.

We recommend installing ``caliban`` using `pipx
<https://pypi.org/project/pipx/>`_. `pipx <https://pypi.org/project/pipx/>`_ is
a tool that lets you install command line utilities written in Python into their
own virtual environments, completely isolated from your system python packages.

You don't HAVE to do this - you can install caliban in your global environment,
or in a virtualenv - but ``pipx`` is the sanest way we've found to install
Python CLI command tools.

.. NOTE:: Before you install Caliban, you'll need to visit the
          :doc:`prerequisites` page and make sure you have Docker installed and
          the correct version of Python 3.

Install ``pipx`` into your global python environment like this:

.. code-block:: bash

   python3 -m pip install --user pipx
   python3 -m pipx ensurepath

Once ``pipx`` is installed, use it to install ``caliban``:

.. code-block:: bash

   pipx install caliban

If you don't want to use `pipx`, install Caliban via pip:

.. code-block:: bash

   pip install -U caliban

Upgrading Caliban
^^^^^^^^^^^^^^^^^

With ``pipx``\ , upgrading Caliban is simple. The following command will do it:

.. code-block:: bash

   pipx upgrade caliban

If you've installed Caliban with pip:

.. code-block:: bash

   pip upgrade caliban

Check your Installation
^^^^^^^^^^^^^^^^^^^^^^^

To check if all is well, run

.. code-block:: bash

   caliban --help

To take Caliban through its paces, visit the `"Getting Started with Caliban"
<https://github.com/google/caliban#getting-started-with-caliban>`_ tutorial on
the main page of `Caliban's github repository
<https://github.com/google/caliban>`_.
Prerequisites
-------------

Before you can use Caliban, you'll need to install Docker and make sure your
Python 3 is up to date. Follow these steps to get set up.

Python 3
^^^^^^^^

Caliban requires Python >= 3.6. Check your current version at the terminal:

.. code-block:: bash

   $ python3 --version
   Python 3.6.9 # Or something above 3.6.0

If you need to upgrade:

- on MacOS, download `the latest Python from python.org
  <https://www.python.org/downloads/mac-osx>`_.
- On Linux, make sure your ``python3`` is up to date by running the following
  command at your terminal:

.. code-block:: bash

   sudo apt-get install python3 python3-venv python3-pip

Once that's all set, run ``python3 --help`` again to verify that you're running
python 3.6 or above.

Docker
^^^^^^

To use Caliban, you'll need a working Docker installation. If you have a GPU and
want to run jobs that use it, you'll have to install ``nvidia-docker2``, as
described below in :ref:`GPU Support on Linux Machines`

- On MacOS, install `Docker Desktop for Mac
  <https://hub.docker.com/editions/community/docker-ce-desktop-mac>`_. You'll
  only be able to run in CPU mode, as MacOS doesn't support Docker's nvidia
  runtime. You will, however, be able to build GPU containers and submit them to
  Google Cloud.
- On Linux, install Docker with `these instructions
  <https://docs.docker.com/install/linux/docker-ce/ubuntu/>`_.

Add your username to the docker group so that you can run Docker without using
``sudo``:

.. code-block:: bash

   sudo usermod -a -G docker ${USER}

GPU Support on Linux Machines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On Linux, Caliban can run jobs locally that take advantage of a GPU you may have installed.

To use this feature, install the ``nvidia-docker2`` runtime by following the
instructions at the `nvidia-docker2
<https://github.com/NVIDIA/nvidia-docker/wiki/Installation-(version-2.0)>`_
page.

.. NOTE:: It's important that you install ``nvidia-docker2``, not
          ``nvidia-docker``! The `nvidia-docker2
          <https://github.com/NVIDIA/nvidia-docker/wiki/Installation-(version-2.0)>`_
          instructions discuss how to upgrade if you accidentally install
          ``nvidia-docker``.

.. NOTE:: The most recent versions of docker don't need the ``nvidia-docker2``
          dependency. In a future version of Caliban we'll remove this
          dependency and upgrade the documentation.
Setting up Google Cloud
=======================

This document will guide you through the process of configuring your machine to
submit jobs to Google's Cloud AI Platform with :doc:`../cli/caliban_cloud`.

.. note:: Caliban will eventually support other Cloud providers like AWS, but
   Google Cloud will remain a great starting option. Google Cloud gives you $300
   of credit, so you can get started immediately with
   :doc:`../cli/caliban_cloud`.

By the end you'll be able to complete the `final step
<https://github.com/google/caliban#submitting-to-cloud-ai-platform>`_ of the
`"Getting Started with Caliban"
<https://github.com/google/caliban#getting-started-with-caliban>`_ tutorial and
train an ML model on AI Platform.

**These are the steps we'll cover**:

.. contents:: :local:
   :depth: 1

To make it through, you'll **need to know**:

- Know how to `save environment variables on your machine
  <https://scotch.io/tutorials/how-to-use-environment-variables>`_ by saving
  them in your ``~/.bashrc`` file.

Create a Cloud Account
----------------------

To submit jobs to AI Platform you need to create a Google Cloud account with an
active "`project
<https://cloud.google.com/resource-manager/docs/creating-managing-projects>`_".
Every job you submit to Cloud will be associated with this project.

Visit the `Google Cloud Console <https://console.cloud.google.com>`_ and click
"Select a Project":

.. image:: /_static/img/cloud/select_project.png
  :width: 600
  :align: center
  :alt: Select a Project

Click "New Project" in the dialogue that appears:

.. image:: /_static/img/cloud/new_project.png
  :width: 600
  :align: center
  :alt: New Project

Give your project a memorable name like "totoro-lives" and click Create. This
should take you to your new project's dashboard.

Note the **Project ID** in the "Project info" panel:

.. image:: /_static/img/cloud/project_id.png
  :width: 600
  :align: center
  :alt: Project ID

Caliban will use this project ID to submit jobs to the correct project in Cloud.

Add the following line to the file ``~/.bashrc`` or ``~/.bash_profile`` on your
machine, to make the ID available to Caliban:

.. code-block:: bash

   export PROJECT_ID=<your-project-id>

.. note:: If you don't know what this means, see `this page
          <https://scotch.io/tutorials/how-to-use-environment-variables>`_ for a
          tutorial on environment variables. We'll remind you at the end of the
          tutorial which variables you'll need.


Activate Free Cloud Trial
-------------------------

Every new Google Cloud project comes with a free $300 credit. To activate this,
click "Activate" at the top right of your new Cloud account's console and follow
the prompts.

.. image:: /_static/img/cloud/activate.png
  :width: 600
  :align: center
  :alt: Activate Billing

You'll have to set up a billing account as well.

The system will ask you for a credit card to verify your identity, but if you
use up the entire credit it won't automatically charge you. You can decide at
that point whether you'd like to continue or not.

Enable AI Platform and Container Registry
-----------------------------------------

Google Cloud has a `dizzying number of products
<https://cloud.google.com/products>`_. To submit jobs with
:doc:`/cli/caliban_cloud`, you'll need to activate just these two:

- `Cloud AI Platform <https://cloud.google.com/ai-platform/docs>`_
- `Container Registry <https://cloud.google.com/container-registry/docs/quickstart>`_

Follow the instructions at `this link to enable the Container Registry API
<https://console.cloud.google.com/flows/enableapi?apiid=containerregistry.googleapis.com&redirect=https://cloud.google.com/container-registry/docs/quickstart&_ga=2.204958805.498449691.1592416944-1401171737.1587152715>`_
by selecting the project you created above and clicking "Continue".

Click `this link to Enable the AI Platform Jobs API
<https://console.cloud.google.com/ai-platform/ml-enable-api/jobs>`_ by clicking
"Enable API" and waiting for the spinner to stop.

Install the Cloud SDK
---------------------

The final step is to install the `Google Cloud SDK
<https://cloud.google.com/sdk/install>`_ on your machine.

Visit the `Google Cloud SDK installation page
<https://cloud.google.com/sdk/install>`_ for a full set of installation
instructions. Here is the distilled version:

- For MacOS, run the `interactive installer
  <https://cloud.google.com/sdk/docs/downloads-interactive>`_.
- For Linux, use `apt-get
  <https://cloud.google.com/sdk/docs/downloads-apt-get>`_ to get the latest
  release.

If you didn't do it during installation, you'll need to initialize the SDK with this command:

.. code-block:: bash

   gcloud init

When you see this output:

.. code-block:: bash

   You are logged in as: [totoro@gmail.com].

   This account has a lot of projects! Listing them all can take a while.
    [1] Enter a project ID
    [2] Create a new project
    [3] List projects
   Please enter your numeric choice:

Enter ``1``, then type in your project ID you noted earlier. (You should have
saved it as an environment variable called ``$PROJECT_ID``).

If you'd like to set a default zone, anything beginning with ``us-central1`` is
a great choice. ``us-central1`` has `the most capability
<https://cloud.google.com/ml-engine/docs/regions>`_ of any region.

.. note:: the Cloud SDK is quite powerful, and gives you access to Cloud buckets
          and all sorts of Google services. You might want to peruse the full
          set of `SDK documentation
          <https://cloud.google.com/sdk/gcloud/reference/>`_ once you've got
          everything working.

Configure Docker Authentication
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To submit a job with :doc:`/cli/caliban_cloud`, Caliban needs to push Docker
images to the Container Registry service that you enabled earlier. To allow
Docker to push to the Container Registry, run this command at a terminal:

.. code-block:: bash

   gcloud auth configure-docker

You should see output that includes the text ``Adding credentials for all GCR
repositories.``.

Test your Environment
---------------------

To check if your SDK installation was successful, run ``gcloud auth list`` in
your terminal. You should see your email address listed as the active account:

.. code-block:: bash

   [totoro@totoro ~]$ gcloud auth list
       Credentialed Accounts
   ACTIVE  ACCOUNT
   *       totoro@google.com

   To set the active account, run:
       $ gcloud config set account `ACCOUNT`

As a final step, confirm that you've set the following environment variables.
(If you set a custom region above, add it here as a ``$REGION`` variable).

.. code-block:: bash

   export REGION="us-central1"
   export PROJECT_ID="research-3141"

If you have all of this, you're set!

Train a model in Cloud
----------------------

Now that you have a working Cloud configuration and a new project, you can use
:doc:`/cli/caliban_cloud` to submit jobs to Cloud AI platform.

The `"Getting Started with Caliban"
<https://github.com/google/caliban#getting-started-with-caliban>`_ tutorial ends
with a nice demo that has you training models in Cloud. Head over to `the
tutorial <https://github.com/google/caliban#getting-started-with-caliban>`_ and
complete the `final step
<https://github.com/google/caliban#submitting-to-cloud-ai-platform>`_ to train a
digit-classifying neural network on AI Platform.
