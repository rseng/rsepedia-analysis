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
