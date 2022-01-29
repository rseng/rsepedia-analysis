Dask is a community maintained project. We welcome contributions in the form of bug reports, documentation, code, design proposals, and more. 

For general information on how to contribute see https://docs.dask.org/en/latest/develop.html.
See [developer documentation](https://docs.dask.org/en/latest/develop.html)
for tips on how to get started.
- [ ] Closes #xxxx
- [ ] Tests added / passed
- [ ] Passes `pre-commit run --all-files`
## HDFS testing on Travis CI

Dask & HDFS testing relies on a docker container. The tests are setup to run on
Travis CI, but only under the following conditions:

- Merges to main
- PRs where the commit message contains the string `"test-hdfs"`

If you make a PR changing HDFS functionality it'd be good to have the HDFS
tests run, please add `"test-hdfs"` to your commit message.

## Setting up HDFS testing locally using Docker

Assumes docker is already installed and the docker-daemon is running.

From the root directory in the repo:

- First get the docker container:

```bash
# Either pull it from docker hub
docker pull daskdev/dask-hdfs-testing

# Or build it locally
docker build -t daskdev/dask-hdfs-testing continuous_integration/hdfs/
```

- Start the container and wait for it to be ready:

```bash
source continuous_integration/hdfs/startup_hdfs.sh
```

- Install dependencies and dask on the container

```bash
source continuous_integration/hdfs/install.sh
```

- Run the tests

```bash
source continuous_integration/hdfs/run_tests.sh
```

- Alternatively, you can start a terminal on the container and run the tests
  manually. This can be nicer for debugging:

```bash
# CONTAINER_ID should be defined from above, but if it isn't you can get it from
export CONTAINER_ID=$(docker ps -l -q)

# Start the bash session
docker exec -it $CONTAINER_ID bash

# Test just the hdfs tests
py.test dask/bytes/tests/test_hdfs.py -s -vv
```
There are a variety of other projects related to dask that are often
co-released.  We may want to check their status while releasing


Release per project:

*   Raise an issue in the https://github.com/dask/community issue tracker
    signaling your intent to release and the motivation.  Let that issue
    collect comments for a day to ensure that other maintainers are comfortable
    with releasing.

*   Update release notes in docs/source/changelog.rst

*   Commit

        git commit -a -m "bump version to x.x.x"

*   Tag commit

        git tag -a x.x.x -m 'Version x.x.x'

*   Push to GitHub

        git push dask main --tags

*   Upload to PyPI

        git clean -xfd
        python setup.py sdist bdist_wheel
        twine upload dist/*

*   Wait for [conda-forge](https://conda-forge.github.io) bots to track the
    change to PyPI

    This will typically happen in an hour or two.  There will be two PRs, one
    to `dask-core`, which you will likely be able to merge after tests pass,
    and another to `dask`, for which you might have to change version numbers
    if other packages (like `distributed`) are changing at the same time.

    In some cases you may also have to zero out build numbers.

    If for some reason you have to do this manually, then follow these steps:

    *  Update conda-smithy and run conda-smithy rerender

            git clone git@github.com:conda-forge/dask-core-feedstock
            cd dask-core-feedstock
            conda install conda-smithy
            conda-smithy rerender

    *  Get sha256 hash from pypi.org
    *  Update version number and hash in recipe
    *  Check dependencies
    *  Do the same for the dask-feedstock meta-package

*   Automated systems internal to Anaconda Inc then handle updating the
    Anaconda defaults channel
