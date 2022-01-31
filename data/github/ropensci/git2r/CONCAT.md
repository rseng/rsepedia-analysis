[![R-CMD-check](https://github.com/ropensci/git2r/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/git2r/actions)
[![CRAN status](http://www.r-pkg.org/badges/version/git2r)](http://cran.r-project.org/web/packages/git2r/index.html)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/last-month/git2r)](http://cran.r-project.org/web/packages/git2r/index.html)
[![Coverage Status](https://coveralls.io/repos/github/ropensci/git2r/badge.svg?branch=master)](https://coveralls.io/github/ropensci/git2r?branch=master)

# Introduction

The `git2r` package gives you programmatic access to Git repositories
from R. Internally the package uses the libgit2 library which is a
pure C implementation of the Git core methods. For more information
about libgit2, check out libgit2's website
[(http://libgit2.github.com)](http://libgit2.github.com).

Suggestions, bugs, forks and pull requests are appreciated. Get in
touch.

## Installation

To install the version available on CRAN:

```coffee
install.packages("git2r")
```

To install the development version of `git2r`, it's easiest to use the
devtools package:

```coffee
# install.packages("devtools")
library(devtools)
install_github("ropensci/git2r")
```

Another alternative is to use `git` and `make`

```coffee
$ git clone https://github.com/ropensci/git2r.git
$ cd git2r
$ make install
```

## Usage

### Repository

The central object in the `git2r` package is the S3 class
`git_repository`. The following three methods can instantiate a
repository; `init`, `repository` and `clone`.

#### Create a new repository

Create a new repository in a temporary directory using `init`

```coffee
library(git2r)
```

```
#> Loading required package: methods
```

```coffee

## Create a temporary directory to hold the repository
path <- tempfile(pattern="git2r-")
dir.create(path)

## Initialize the repository
repo <- init(path)

## Display a brief summary of the new repository
repo
```

```
#> Local:    /tmp/Rtmp7CXPlx/git2r-1ae2305c0e8d/
#> Head:     nothing commited (yet)
```

```coffee

## Check if repository is bare
is_bare(repo)
```

```
#> [1] FALSE
```

```coffee

## Check if repository is empty
is_empty(repo)
```

```
#> [1] TRUE
```

#### Create a new bare repository

```coffee
## Create a temporary directory to hold the repository
path <- tempfile(pattern="git2r-")
dir.create(path)

## Initialize the repository
repo <- init(path, bare=TRUE)

## Check if repository is bare
is_bare(repo)
```

```
#> [1] TRUE
```

#### Clone a repository

```coffee
## Create a temporary directory to hold the repository
path <- file.path(tempfile(pattern="git2r-"), "git2r")
dir.create(path, recursive=TRUE)

## Clone the git2r repository
repo <- clone("https://github.com/ropensci/git2r", path)
```

```
#> cloning into '/tmp/Rtmp7CXPlx/git2r-1ae27d811539/git2r'...
#> Receiving objects:   1% (24/2329),   12 kb
#> Receiving objects:  11% (257/2329),   60 kb
#> Receiving objects:  21% (490/2329),  100 kb
#> Receiving objects:  31% (722/2329),  125 kb
#> Receiving objects:  41% (955/2329),  237 kb
#> Receiving objects:  51% (1188/2329),  574 kb
#> Receiving objects:  61% (1421/2329), 1014 kb
#> Receiving objects:  71% (1654/2329), 1350 kb
#> Receiving objects:  81% (1887/2329), 1733 kb
#> Receiving objects:  91% (2120/2329), 2614 kb
#> Receiving objects: 100% (2329/2329), 2641 kb, done.
```

```coffee

## Summary of repository
summary(repo)
```

```
#> Remote:   @ origin (https://github.com/ropensci/git2r)
#> Local:    master /tmp/Rtmp7CXPlx/git2r-1ae27d811539/git2r/
#>
#> Branches:          1
#> Tags:              0
#> Commits:         320
#> Contributors:      3
#> Ignored files:     0
#> Untracked files:   0
#> Unstaged files:    0
#> Staged files:      0
```

```coffee

## List all references in repository
references(repo)
```

```
#> $`refs/heads/master`
#> [6fb440] master
#>
#> $`refs/remotes/origin/master`
#> [6fb440] origin/master
```

```coffee

## List all branches in repository
branches(repo)
```

```
#> [[1]]
#> [6fb440] (Local) (HEAD) master
#>
#> [[2]]
#> [6fb440] (origin @ https://github.com/ropensci/git2r) master
```

#### Open an existing repository

```coffee
## Open an existing repository
repo <- repository(path)

## Workdir of repository
workdir(repo)
```

```
#> [1] "/tmp/Rtmp7CXPlx/git2r-1ae27d811539/git2r/"
```

```coffee

## List all commits in repository
commits(repo)[[1]] # Truncated here for readability
```

```
#> Commit:  6fb440133765e80649de8d714eaea17b114bd0a7
#> Author:  Stefan Widgren <stefan.widgren@gmail.com>
#> When:    2014-04-22 21:43:19
#> Summary: Fixed clone progress to end line with newline
```

```coffee

## Get HEAD of repository
repository_head(repo)
```

```
#> [6fb440] (Local) (HEAD) master
```

```coffee

## Check if HEAD is head
is_head(repository_head(repo))
```

```
#> [1] TRUE
```

```coffee

## Check if HEAD is local
is_local(repository_head(repo))
```

```
#> [1] TRUE
```

```coffee

## List all tags in repository
tags(repo)
```

```
#> list()
```

### Configuration

```coffee
config(repo, user.name="Git2r Readme", user.email="git2r.readme@example.org")

## Display configuration
config(repo)
```

```
#> global:
#>         core.autocrlf=input
#> local:
#>         branch.master.merge=refs/heads/master
#>         branch.master.remote=origin
#>         core.bare=false
#>         core.filemode=true
#>         core.logallrefupdates=true
#>         core.repositoryformatversion=0
#>         remote.origin.fetch=+refs/heads/*:refs/remotes/origin/*
#>         remote.origin.url=https://github.com/ropensci/git2r
#>         user.email=git2r.readme@example.org
#>         user.name=Git2r Readme
```

### Commit

```coffee
## Create a new file
writeLines("Hello world!", file.path(path, "test.txt"))

## Add file and commit
add(repo, "test.txt")
commit(repo, "Commit message")
```

```
#> Commit:  0a6af48cedf43208bde34230662280514e0956eb
#> Author:  Git2r Readme <git2r.readme@example.org>
#> When:    2014-04-22 21:44:57
#> Summary: Commit message
```

# Included software

- The C library [libgit2](https://github.com/libgit2/libgit2). See
  `inst/AUTHORS` for the authors of libgit2.

- The libgit2 library has been modified, e.g. to use the R printing
  and error routines, and to use `runif` instead of `rand`.

# License

The `git2r` package is licensed under the GPLv2. See these files for additional details:

- LICENSE      - `git2r` package license (GPLv2)
- inst/COPYING - Copyright notices for additional included software

---

[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# git2r 0.29.0 (2021-11-18)

## CHANGES

* Added a 'branch' argument to the 'init' function to make it possible
  to specify the branch name.

* Updated the build configuration script on Windows and MacOS to use
  libgit2 version 1.3.0.

* Updated the bundled libgit2 source code to version 1.3.0.

* Renamed the NEWS file to NEWS.md and changed to use markdown format
  style.

# git2r 0.28.0 (2021-01-10)

## IMPROVEMENTS

* Updated to use libgit2 version 1.1.0 on Windows.

* Fix handling of a symbolic reference when checking out previous
  branch.

* Added a configure option '--without-libgit2' to ignore presence of a
  system libgit2 library and instead use the internal git2r libgit2
  library. Usage:
  R CMD INSTALL --configure-args='--without-libgit2' git2r_x.y.z.tar.gz

* Updated some tests to work with libgit2 version 1.1.0.

# git2r 0.27.1 (2020-05-03)

## CHANGES

* Fixed the CITATION file to pass 'R CMD check' without a NOTE.

# git2r 0.27.0 (2020-05-01)

## IMPROVEMENTS

* Updated the bundled libgit2 source code to version '0.28.5'.

* Updated the build configuration script to be able to build git2r
  with a system installation of libgit2 version >= 1.0.

* Updated to use libgit2 version 1.0.0 on Windows.

* The build configuration script checks for minimum required version
  of libssh2 (version >= 1.8). Issue #420.

* Updated to use roxygen2 version 7.1.0 to build the documentation.

* Make it easier to view and change the timezone (John Blischak in
  #408).

* Fixed 'ls_tree' to handle content in subfolder, see description in
  PR #402.

* The 'branch_create' function has been changed to use the
  'last_commit()' function as default to determine the commit to which
  the new branch should point.

# git2r 0.26.1 (2019-06-30)

## BUG FIXES

* Fixed the broken build on Solaris.

# git2r 0.26.0 (2019-06-29)

## IMPROVEMENTS

* Updated the bundled libgit2 source code to version '0.28.2'.

* Added the 'force' argument to the 'tag' function to overwrite an
  existing tag.

* Allow a zero length tag message.

* Make it possible to create a lighweight tag.

* Added the 'ref' argument to the 'commits' function to give a
  reference to list commits from.

* Added the utility function 'lookup_commit' to lookup a commit
  related to a git object.

* The 'path' argument was added to the 'commits' function to make it
  possible to specify that only commits modifying this file ('path')
  will be returned to reproduce 'git log' with '--no-follow', see the
  documentation. (Peter Carbonetto and John Blischak in PR #372)

## BUG FIXES

* Removed the timezone offset from the commit time to fix an incorrect
  time in GMT when reading information from a repository (Thierry
  Onkelinx in PR #393).

# git2r 0.25.2 (2019-03-20)

## CHANGES

* Improved the build configuration script: if the system installation
  of libgit2 is to old, use the bundled libgit2 instead of raising an
  error.

## BUG FIXES

* Fixed the broken build on Solaris.

# git2r 0.25.1 (2019-03-17)

## BUG FIXES

* Fixed significant warning from 'R CMD check'

# git2r 0.25.0 (2019-03-17)

## CHANGES

* Updated the bundled libgit2 source code to version '0.28.1'.

* Added additional parameters to the 'diff' function to control the
  output, see the documentation.

* Added getPass option to the password argument in 'cred_user_pass'
  (Annie Wang in PR #383)

* Changed the 'print' functions to return its argument invisibly.

* Changed the 'git_config_files' function to return a 'data.frame'

* Changed the 'ahead_behind' function to accept a tag or a branch for
  the local and upstrean commit.

* Changed the 'descendent_of' function to accept a tag or a branch for
  the 'commit' and 'ancestor' commit.

## BUG FIXES

* Fixed memory protection errors in the git2r C source code reported
  by the 'rchk' tool.

* Fixed listing of 'commits' from a shallow repository.

* Fixed the configuration script to include the missing macro
  'AM_ICONV'.

# git2r 0.24.0 (2019-01-07)

This is a bug-fix release.

## BUG FIXES

* Fixed memory protection errors in the git2r C source code reported
  by the 'rchk' tool.

* Raise an error if the path argument to the 'hashfile' function is
  NA.

# git2r 0.23.0 (2018-07-17)

## IMPROVEMENTS

* Updated the bundled libgit2 source code to v0.27.3 (504bd54).

## BREAKING CHANGE

* On macOS, git2r no longer enables SSH transport by default. This is
  due to the complexity to build the dependencies for SSH transport in
  an R package when macOS no longer ships the OpenSSL headers.
  However, you can install git2r from source on macOS (see the 'R
  Installation and Administration' manual) with SSH transport enabled
  if you first install the libgit2 library, for example, using the
  Homebrew package manager. Another possibility is to let the build
  configuration automatically download the libgit2 library from the
  Homebrew package manager with:

  install.packages('git2r', type='source', configure.vars='autobrew=yes')

# git2r 0.22.1 (2018-07-10)

## NEW FEATURES

* Added the 'git_config_files' method to locate configuration files.

* Added the 'stash_pop' method to apply a single stashed state from
  the stash list and remove it from the list if successful.

* Added the 'stash_apply' method to apply a single stashed state from
  the stash list.

## IMPROVEMENTS

* Updated the bundled libgit2 source code to v0.27.2 (8d36dc6).

* git2r can now build against a system installation of libgit2
  (Elliott Sales de Andrade in PR #345, #344 and #336).

* Refactoring of the configuration scripts to use a prebuilt libgit2
  on macOS and Windows (Thanks Jeroen).

* Ensure that git2r writes the config file to the correct location on
  Windows (John Blischak in PR #320).

* Better default location to find ssh keys in 'cred_ssh_key()' (Ian
  Lyttle in PR #317).

## BUG FIXES

* If a merge results in no change, the returned 'git_merge_result' now returns
  'FALSE' for 'fast_forward' and 'conflicts' and 'NA' for 'sha'. Previously it
  returned 'logical(0)' for 'fast_forward' and 'conflicts' and 'character(0)'
  for 'sha'.

## BREAKING CHANGES

* Changed from S4 classes to S3 classes to simplify the design and
  facilitate future development.

* Removed the trailing slash from the directory name when reporting
  repository path or workdir.

* Removed the 'libgit2_sha' method. Use the 'libgit2_version' method
  instead.

* Changed the 'stash_drop' argument 'index' from zero-based to
  one-based i.e. use index = 1 to drop the first stash.

# git2r 0.21.0 (2018-01-04)

* Added methods 'odb_blobs' and 'odb_objects' with missing repository
  signature. Internally, they use 'getwd' and 'discover_repository' to
  open a repository.

## BUG FIXES

* The bundled libgit2 source code has been reverted to libgit2 v0.26.0
  (15e1193) from 14 June 2017 (same as in git2r v0.19.0) to fix memory
  alignment errors.

# git2r 0.20.0 (2017-12-17)

## IMPROVEMENTS

* Updated the bundled libgit2 source code to commit (fa8cf14) from 16
  December 2017.

* Improvements to the build configuration script.

## BUG FIXES

* Fixed the internal callback for remote host authentication from
  hanging indefinitely when querying an ssh-agent for
  credentials. Now, the callback signals an error instead of trying
  again if the authentication failed the first time.

# git2r 0.19.0 (2017-07-19)

## IMPROVEMENTS

* Updated the bundled libgit2 source code to commit (15e1193)
  (v0.26.0) from 14 June 2017.

* Added 'checkout' argument to 'clone()'. Allows to control whether
  checkout of HEAD is performed after the clone is complete. Setting
  'checkout=FALSE' has similar effect as the git command line option
  '--no-checkout'. Andrzej K. Oles in #282.

## BUG FIXES

* Fixed memory protection errors in the git2r C source code reported
  by the 'rchk' tool.

* Added missing calls to 'R_forceSymbols' and 'R_useDynamicSymbols' in
  the C init function.

* Enable argument 'all' to commit multiple modified (or deleted)
  files. John Blischak in #283

* Changed the configure script to determine the architecture of the
  machine earlier in order to fix an unsupported architecture error
  encountered on CentOS (#268, #288).

# git2r 0.18.0 (2017-01-01)

## BUG FIXES

* This is a bug-fix release to solve an error introduced in the build
  configuration on mac in version 0.17.0. The package failed with
  'unable to load shared object', see issue #267.

# git2r 0.17.0 (2016-12-29)

## IMPROVEMENTS

* Updated the bundled libgit2 source code to commit (6b0510e) from 20
  December 2016.

* Static linking of LibSSH2 on mac to support redistributable binary
  package with SSH transport enabled. Version 1.8.0 of LibSSH2 is
  downloaded and built from 'https://www.libssh2.org/download/'.

# git2r 0.16.0 (2016-11-20)

## IMPROVEMENTS

* Updated libgit2 source code to commit (6b0510e) from
  17 November 2016.

* Add the option 'all_untracked' to the 'status' method to show
  individual files in untracked directories if the 'untracked' option
  is TRUE.

* Add the 'tag_delete' method to delete an existing tag reference.

* Update build configuration to support OpenSSL 1.1.0.

* If the the 'getPass' package is installed the 'cred_ssh_key' method
  to create a new passphrase-protected ssh key credential object will
  call the 'getPass::getPass()' method if the private key is
  passphrase protected to allow for interactive input of the
  passphrase. The 'getPass' package is a suggested package. (Peter
  Meissner in PR #254)

* Add 'path' argument to the 'reset' method to enable path-specific
  unstage, i.e. resets the index entries for all paths to their state
  at HEAD

## BUG FIXES

* Build configuration: use portable string equality comparison
  operator. This fixes the build e.g. for those without Bash as
  /bin/sh. (Sander Maijers in PR #243).

# git2r 0.15.0 (2016-05-11)

## IMPROVEMENTS

* Build configuration: 'pkg-config' is now used to find 'libssl', if
  possible (Elias Pipping in PR #231).

* Added a method to coerce a 'git_commit' object to a 'data.frame'.

* Added the method 'is_branch' to check if an object is a
  'git_branch'.

## BUG FIXES

* Build configuration: fixed installation with parallel make (Kirill
  Müller in PR #228).

# git2r 0.14.0 (2016-03-13)

## IMPROVEMENTS

* Updated libgit2 source code to commit (785d8c48) from
  2016-03-04. This is release v0.24.0 of libgit2.

* Refactoring of the build scripts.

* Added a check that the configuration key is valid when setting a
  configuration variable and output a warning if the key is invalid.

* The status method now prints "working directory clean" instead of
  nothing when the working directory is clean.

* Added the 'refspec' argument to the 'fetch' method to specify the
  refs to fetch and which local refs to update.

* Added a workaround to the 'commit' method to list commits in a
  shallow clone since the libgit2 library does not yet support this.

# git2r 0.13.1 (2015-12-10)

## BUG FIXES

* This is a bug-fix release to solve problems introduced in version
  0.12.1:

  - The bundled libgit2 source code has been reverted to commit
    (98f7bd2) from 2015-08-05 (same as in v0.11.0) to fix memory
    alignment errors (clang-UBSAN and gcc-UBSAN).

  - OpenSSL is now used again on OS X to provide the cryptographic
    support for HTTPS connections to fix a significant compiler
    warning (arithmetic on a pointer to void is a GNU extension
    [-Wpointer-arith]) on r-devel-osx-x86_64-clang.

  - Several fixes to the build configuration on non-Windows platforms.

# git2r 0.12.1 (2015-12-05)

## NEW FEATURES

* Add 'remote_ls' method to list references in a remote repository akin to the
  `git ls-remote` command.

* Add 'remote_set_url' method to set the remote's url in the
  configuration.

* Add 'cred_token' S4 class to hold the name of the environmental
  variable with the secret. Default value for the name is GITHUB_PAT.

* It is now possible to checkout a specific file with the 'checkout'
  method.

* Add 'ssl_cert_locations' method to set libgit2 global option
  'GIT_OPT_SET_SSL_CERT_LOCATIONS'

* Add 'ceiling' argument to 'discover_repository' method to prevent
  search from walking up the parent directories.

## CHANGES

* Improvments to the cred_* functions documentation.

* Add the following default arguments to the 'cred_ssh_key' method:
  publickey = '~/.ssh/id_rsa.pub' and privatekey = '~/.ssh/id_rsa'

* On OSX, cascade CPPFLAGS and LDFLAGS to libssh2 build to allow
  libssh2 to be built against a user-installed openssl, discovered by
  configure or from R's Makeconf. Necessary to build on OS X ≥ 10.11

* On OS X, SecureTransport is now used to provide the cryptographic
  support for HTTPS connections insead of OpenSSL.

* The search for libssh2 during configuration (non Windows) is now
  done via pkg-config.

* Update OpenSSL on Windows to v1.0.2d

* Update libgit2 source code to commit (3f5877d) from 2015-11-12.

## BUG FIXES

* Add missing credentials argument to pull method.

* Fix config failure when user.name or user.email are passed as
  variables.

* Include 'configure.ac' in the distribution.

# git2r 0.11.0 (2015-08-12)

## NEW FEATURES

* Add punch card plot.

* Add branch argument to clone with name of the branch to checkout.

* Add 'force' argument to 'add' method to add ignored files.

* The following methods can now be called without the repository
  argument: 'branches', 'references', 'remotes', 'tags' and 'workdir'.
  When these methods are called without the repository argument, the
  repository is searched for with 'discover_repository' in the current
  working directory.

* Add name of branch to each item in branch_list.

* Add name of tag to each item in tags list.

* Add S4 class 'cred_env' to pass credentials in environment
  variables.

* SSH transport on Windows. This requires 'LibSSH2' and
  'OpenSSL'. These two libraries are downloaded from
  'https://github.com/rwinlib' during configuration of the package.

* Static linking of LibSSH2 on OSX to support redistributable binary
  package with SSH transport enabled. Version 1.6.0 of LibSSH2 is
  downloaded and built from 'https://github.com/libssh2/libssh2'.

## IMPROVEMENTS

* Better summary output from S4 classes 'git_commit' and
  'git_repository'.

* Updated libgit2 source code to commit (98f7bd2) from 2015-08-05.

## BUG FIXES

* Add imports to DESCRIPTION to fix CRAN notes.

* Fix plot function to use the repository argument 'x'

* Update configuration to build on OpenBSD.

* Fix checkout branch in empty repository.

* Fix path argument in rm_file.

* Internal refactoring of C code that raise error to prevent segfault.

# git2r 0.10.1 (2015-05-07)

## CHANGES

* Rename 'bundle_repo' method to 'bundle_r_package'

# git2r 0.10.0 (2015-05-07)

## NEW FEATURES

* Added method libgit2_sha that returns the commit id of the libgit2
  library that the bundled source code is based on.

* Added the method in_repository to determine if a directory is in a
  git repository.

## CHANGES

* Add brief summary of the five latest commits when summarizing a
  git_respository.

* Added argument 'n' to the commits method to limit the number of
  commits in the output.

* Added the following methods with missing repository signature;
  commits, is_shallow, is_empty, is_detached, repository and
  status. Internally, these methods use getwd and discover_repository
  to open a repository.

* Changed configuration to raise error if the OpenSSL library is not
  found on non-Windows systems.

* Changed configuration to raise error if the iconv library is not
  found on OSX.

* Removed print of the configuration in the config method. Changed to
  return S3 class git_config.

* Removed print of the status in the status method. Changed to return
  S3 class git_status.

## BUG FIXES

* Use OPENSSL_INCLUDES variable to build on Solaris.

* Use bundled regex library on Solaris.

git2 0.9 (2015-04-25)

## CHANGES

* Single quote 'libgit2' and 'Git' in Description field

git2 0.8 (2015-04-24)

## CHANGES

* Added bare argument to clone method to create a bare repository

* Added force argument to push to force local revision to the remote
  repo

* Updated libgit2 source code (2a0f67f)

* Internal refactoring of push

## NEW FEATURES

* Added method rm_file to remove files

* Added 'all' argument to commit method to stage modified and deleted
  files

* Added shortcut to checkout previous branch with "-" which is
  synonymous with "@{-1}"

* Added session argument to commit method to add sessionInfo to commit
  message

* Added session argument to tag method to add sessionInfo to tag
  message

* Added  method to coerce POSIXlt to S4 class git_time

* Added method 'revparse_single' to find object specified by revision

* Added plot method

git2 0.7 (2015-02-23)

## CHANGES

* Update libgit2 source code to commit (366e53d)

* Fix configuration of compiler options when the OpenSSL library is
  found on non-Windows platforms

# git2r 0.6 (2015-02-18)

## CHANGES

* Update Title and Description field in DESCRIPTION file.

# git2r 0.5 (2015-02-17)

## CHANGES

* Update libgit2 source code to commit (a291790)

* Use Alice and Bob as placeholder names in examples.

* Add COPYRIGHTS file to list all copyright holders.

* Fix significant compiler warnings from R CMD check with pedantic
  flag.

# git2r 0.4 (2015-01-13)

## CHANGES

* Fix build on Windows

# git2r 0.3 (2015-01-13)

## CHANGES

* Internal refactoring of merge method and merge tests.

* Update libgit2 source code to version v0.22.0

## BUG FIXES

* Fix build on OSX.

# git2r 0.2 (2015-01-05)

## NEW FEATURES

* Add method 'odb_objects' to list all objects available in the
  database as a data.frame

* Add method 'odb_blobs' to list all blobs reachable from the commits
  in the object database.

## DOCUMENTATION

* Added examples to all exported methods.

## CHANGES

* Removed ggplot2 dependency. Moved plot functionality to the ggit
  package (https://github.com/ropensci/ggit).

* Renamed note_list method to notes

* Removed markdown_link method

* Renamed diff and merge arguments

## IMPROVEMENTS

* Better performance when summarizing contributions.

* Improved build of package.

## BUG FIXES

* Fixed memory leaks.

* Fixed use of allocVector without protection.

* Added unload hook to unload DLL.

* Fix tree and blob tests to use writeChar instead of writeLines to
  have more control over line endings.

# git2r 0.1 (2014-09-09)

## NEW FEATURES

* Many new features and methods added, see the documention for a
  description of the methods below:

  - Blob: content, blob_create, hash, hashfile, is_binary, is_blob,
          length, show, summary.
  - Branch: branch_create, branch_delete, branch_get_upstream,
            branch_remote_name, branch_remote_url, branch_rename,
            branch_set_upstream and branch_target.
  - Commit: is_commit and parents.
  - Diff: diff and diff_print.
  - Fetch: fetch and fetch_heads.
  - Libgit2: libgit2_features and libgit2_version.
  - Merge: merge.
  - Note: note_create, note_default_ref, note_list and note_remove.
  - Pull: pull.
  - Push: push.
  - Remote: remote_add, remote_remove, remote_rename and remote_url.
  - Repository: discover_repository and is_shallow
  - Reset: reset.
  - Stash: stash, stash_drop, stash_list, show and summary.

* Improved error messages to give more detailed information including
  which function raised the error.

## NEW S4 CLASSES

* The following new S4 classes to handle the libgit2 data structures:
  - cred_ssh_key
  - cred_user_pass
  - git_blame
  - git_blame_hunk
  - git_blob
  - git_diff
  - git_diff_file
  - git_diff_hunk
  - git_diff_line
  - git_fetch_head
  - git_merge_result
  - git_note
  - git_reflog_entry
  - git_stash
  - git_transfer_progress
  - git_tree

## CHANGES

* Renamed methods:
  - is.bare to is_bare
  - is.empty to is_empty
  - is.head to is_head
  - is.local to is_local

* Rename hex to sha for the 40-character SHA-1 hash in method
  arguments and S4 class slots.

# git2r 0.0.8 (2014-03-20)

## NEW FEATURES

* Added method to clone repository

* Added method config to set user.name and user.email in a repository

* Added method status to display state of a repository

# git2r 0.0.7 (2014-03-16)

## NEW FEATURES

* Added method to create a commit

## CHANGES

* Improved error checking

# git2r 0.0.6 (2014-02-21)

## NEW FEATURES

* Added method init to create a new Git repository

## CHANGES

* Removed usage of testthat package when testing the package

* Removed bundled zlib in src/zlib and instead link against zlib
  shipped with R.

* Dropped usage of external pointers, S4 git_repository now keeps
  track of the path of the repository.

# git2r 0.0.5 (2014-01-01)

## CHANGES

* Renamed S4 class repository to git_repository

## NEW FEATURES

* Added method commits to list all commits in repository

* Added S4 class git_commit to hold infformation of a commit

* Added S4 class git_time to hold time of an action

* Added slot walker to S4 class git_repository

# git2r 0.0.4 (2013-12-30)

## NEW FEATURES

* Added method remote_url to get the url a remote in a repository

* Added method workdir to get workdir of a repository

* Added method remotes to list remotes of a repository

* Added S4 class git_signature to hold information of an action
  signature (e.g. for committers, taggers, etc)

## CHANGES

* Renamed S4 class tag to git_tag

* Renamed S4 class branch to git_branch

* Renamed S4 class reference to git_reference

# git2r 0.0.3 (2013-12-29)

## NEW FEATURES

* Added method branches to list branches

* Added method head to retrieve head

* Added method is.head to check if a branch is head

* Added method is.local to check if a branch is local

* Added S4 class branch to hold information of a git branch

* Added method to show a reference

* Added method to list all references in a repository

* Added S4 class reference to hold information of a git reference

# git2r 0.0.2 (2013-12-28)

## NEW FEATURES

* Added is.bare method to check if a repository is bare

* Added is.empty method to check if a repository is empty

# git2r 0.0.1 (2013-12-28)

## NEW FEATURES

* Added S4 class repository to work with a git repository

* Initial package structure
# Welcome to git2r!

Contributions to the `git2r` project are welcome and appreciated. Here
follow some guidelines and description of the project to make it
easier for you to submit pull requests, report issues, and provide
feedback.

## Licensing

By contributing to git2r, you agree to release your contribution under
the terms of the [GPL v2 license](LICENSE).

## Sign your work

To improve tracking of who did what, we've borrowed the "sign-off"
procedure from the Linux kernel project -
[A Developer’s Certificate of Origin](http://elinux.org/Developer_Certificate_Of_Origin).
Although `git2r` is a lot smaller project it is a good discipline to
follow it.

The sign-off is a simple line at the end of the explanation for
the patch, which certifies that you wrote it or otherwise have
the right to pass it on as a open-source patch. The rules are
pretty simple: if you can certify the below:

```
Developer's Certificate of Origin 1.1

By making a contribution to this project, I certify that:

(a) The contribution was created in whole or in part by me and I
    have the right to submit it under the open source license
    indicated in the file; or

(b) The contribution is based upon previous work that, to the best
    of my knowledge, is covered under an appropriate open source
    license and I have the right under that license to submit that
    work with modifications, whether created in whole or in part
    by me, under the same open source license (unless I am
    permitted to submit under a different license), as indicated
    in the file; or

(c) The contribution was provided directly to me by some other
    person who certified (a), (b) or (c) and I have not modified
    it.

(d) I understand and agree that this project and the contribution
    are public and that a record of the contribution (including all
    personal information I submit with it, including my sign-off) is
    maintained indefinitely and may be redistributed consistent with
    this project or the open source license(s) involved.
```

Then you just add a line to every git commit message:

```
Signed-off-by: Random J Developer <random@developer.example.org>
```

Using your real name (sorry, no pseudonyms or anonymous contributions).

This line can be automatically added by Git if you
[configure](http://git-scm.com/book/en/Customizing-Git-Git-Configuration)
your `user.name` and `user.email` and run the git commit command with
the
[-s](https://www.kernel.org/pub/software/scm/git/docs/git-commit.html)
option:

```
git commit -s.
```

## Overall design `git2r <- libgit2 + R`

Internally the git2r package uses the
[libgit2](https://libgit2.github.com/) C library to interact with a Git
repository. C code can be called from R using the `.Call`
method. However, in order to map the R data structures to the C data
structures used by libgit2 an intermediate layer of C code is between
the R methods and libgit2, see figure below.

![Overall design](figure/git2r-design.png)

### R code

The package is based on S4 classes/methods. Each S4 class is prefixed
with `git_` to have the same name as the corresponding libgit2 data
structure. The naming strategy for methods is to use the data
structure (branch, note, stash) followed by the action (create, list,
remove), e.g. `note_create`, `note_list`, and `note_remove`, etc. However, the
naming is not completely consistent in the package, e.g. `init`
instead of `repository_init`.


#### Documentation

The package documentation is written with
[roxygen2](http://cran.r-project.org/web/packages/roxygen2/index.html)
version 4.1.1.  All contributed methods and classes must be
documented. Moreover, all exported methods should have examples of
their usage. Please file an issue if you spot a method in the
documentation without an example, or, better yet, a pull request. The
recommended way to generate man files from the roxygen documentation
is to run the `roxygen` target in the `Makefile`

```
make roxygen
```

### Tests

The tests for the package can be found in the `tests` folder. All
contributions of code should have corresponding tests.

### C code

The `git2r` C code is in files prefixed with `git2r_` in the `src`
folder. The `libgit2` C code is located in `src/libgit2` with
dependencies in `src/http-parser` and `src/regex`.

#### Naming

The git2r C functions are named with the prefix `git2r_module` where
`module` groups related functionality, e.g. `git2r_repository_open` and
`git2r_repository_init`. The source code for each module is in
`git2r_module.c` and `git2r_module.h`.

#### Code layout

1. **Argument checking:** All R arguments must be checked prior to
   their usage. There are several help functions in
   `src/git2r_arg.h` to facilitate argument checking. The `repo`
   argument is checked in the `git2r_repository_open` function and
   returns `NULL` on failure.

2. **Perform the actual work:** The error code (return value) from
   each call to a libit2 function must be checked. In case of an error,
   initiate cleanup.

3. **Cleanup:** All allocated resources must be freed. In case of an
   error, call `Rf_error` with the error message.

#### Documentation

Document functions and parameters with Doxygen to enable generation of
documentatation of the internal git2r C code.

### Code style

It's very important to have a consistent code style, so before
submitting, make sure the code matches surrounding code.

### Makefile

To facilitate development a [Makefile](Makefile) exists in the root
folder with targets for several common tasks. For an introduction to
make, see e.g. the
[minimal make](http://kbroman.github.io/minimal_make/) tutorial by
Karl Broman.

* **roxygen** Rebuild R documentation with
  [roxygen2](http://cran.r-project.org/web/packages/roxygen2/index.html). Checks
  that the correct version of *roxygen2* is installed.

* **check** Build package and then run `R CMD check` to
  [test](http://cran.r-project.org/doc/manuals/r-release/R-exts.html#Checking-packages)
  the package.

* **check_valgrind** Build package and then run `R CMD check` with
  [valgrind](http://valgrind.org/) to detect possible problems.

* **valgrind** Run all test code with valgrind.

* **sync_libgit2** Sync git2r with changes in the libgit2 C-library

  1. Clone or pull libgit2 to parent directory from
     https://github.com/libgit2/libgit2.git

  2. Run `make sync_libgit2`. It first removes files and then copy
     files from libgit2 directory. Next it runs an R script to build
     Makevars.in and Makevars.win based on source files. To
     [pass](http://cran.r-project.org/doc/manuals/r-release/R-exts.html#Writing-portable-packages)
     `R CMD check git2r` the printf calls must be changed to use the R
     [printing](http://cran.r-project.org/doc/manuals/r-release/R-exts.html#Printing)
     routine `Rprintf`. Therefore it runs a `patch` command to change
     some lines in the source code to pass `R CMD check git2r`.

  3. Build and check updated package `make check`

* **Makevars** The code in src is compiled with a
  [Makevars](http://cran.r-project.org/doc/manuals/r-release/R-exts.html#Using-Makevars)
  file ([Makevars.in](src/Makevars.in) and
  [Makevars.win](src/Makevars.win)). They are automatically generated
  from an R [script](tools/build_Makevars.r) when running this target.

* **configure** Generate a `configure` script from `configure.ac` with
  [autoconf](https://www.gnu.org/software/autoconf/). The
  [configure](http://cran.r-project.org/doc/manuals/r-release/R-exts.html#Configure-and-cleanup)
  script is used on non-Windows platforms to generate a
  system-dependent configuration and detect e.g. `OpenSSL` and
  `LibSSH2` before installation.

* **clean** Remove temporary and object files.

### If you have any questions, please don't hesitate to open an issue
# Scripts for configuring/building git2r

http://cran.r-project.org/doc/manuals/r-release/R-exts.html#Configure-and-cleanup

> If your configure script needs auxiliary files, it is recommended
> that you ship them in a tools directory (as R itself does).

The macro `AC_CANONICAL_HOST` in `configure.ac` determines the `host`
with the helper scripts `config.guess` and `config.sub`. This also
requires the script `install-sh`, a wrapper around the system's own
install utility.

## From

```
wget -O config.guess 'https://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.guess;hb=HEAD'
wget -O config.sub 'https://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.sub;hb=HEAD'
wget -O config.rpath 'https://git.savannah.gnu.org/gitweb/?p=gnulib.git;a=blob_plain;f=build-aux/config.rpath;hb=HEAD'
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch.R
\name{fetch}
\alias{fetch}
\title{Fetch new data and update tips}
\usage{
fetch(
  repo = ".",
  name = NULL,
  credentials = NULL,
  verbose = TRUE,
  refspec = NULL
)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{name}{the remote's name}

\item{credentials}{The credentials for remote repository
access. Default is NULL. To use and query an ssh-agent for the
ssh key credentials, let this parameter be NULL (the default).}

\item{verbose}{Print information each time a reference is updated
locally. Default is \code{TRUE}.}

\item{refspec}{The refs to fetch and which local refs to update,
see examples. Pass NULL to use the
\code{remote.<repository>.fetch} variable. Default is
\code{NULL}.}
}
\value{
invisible list of class \code{git_transfer_progress}
    with statistics from the fetch operation:
\describe{
  \item{total_objects}{
    Number of objects in the packfile being downloaded
  }
  \item{indexed_objects}{
    Received objects that have been hashed
  }
  \item{received_objects}{
    Objects which have been downloaded
  }
  \item{total_deltas}{
    Total number of deltas in the pack
  }
  \item{indexed_deltas}{
    Deltas which have been indexed
  }
  \item{local_objects}{
    Locally-available objects that have been injected in order to
    fix a thin pack
  }
  \item{received_bytes}{
    Size of the packfile received up to now
  }
}
}
\description{
Fetch new data and update tips
}
\examples{
\dontrun{
## Initialize three temporary repositories
path_bare <- tempfile(pattern="git2r-")
path_repo_1 <- tempfile(pattern="git2r-")
path_repo_2 <- tempfile(pattern="git2r-")

dir.create(path_bare)
dir.create(path_repo_1)
dir.create(path_repo_2)

bare_repo <- init(path_bare, bare = TRUE)
repo_1 <- clone(path_bare, path_repo_1)
repo_2 <- clone(path_bare, path_repo_2)

config(repo_1, user.name = "Alice", user.email = "alice@example.org")
config(repo_2, user.name = "Bob", user.email = "bob@example.org")

## Add changes to repo 1
writeLines("Lorem ipsum dolor sit amet",
           con = file.path(path_repo_1, "example.txt"))
add(repo_1, "example.txt")
commit(repo_1, "Commit message")

## Push changes from repo 1 to origin (bare_repo)
push(repo_1, "origin", "refs/heads/master")

## Fetch changes from origin (bare_repo) to repo 2
fetch(repo_2, "origin")

## List updated heads
fetch_heads(repo_2)

## Checking out GitHub pull requests locally
path <- tempfile(pattern="ghit-")
repo <- clone("https://github.com/leeper/ghit", path)
fetch(repo, "origin", refspec = "pull/13/head:refs/heads/BRANCHNAME")
checkout(repo, "BRANCHNAME")
summary(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{discover_repository}
\alias{discover_repository}
\title{Find path to repository for any file}
\usage{
discover_repository(path = ".", ceiling = NULL)
}
\arguments{
\item{path}{A character vector specifying the path to a file or
folder}

\item{ceiling}{The default is to not use the ceiling argument and
start the lookup from path and walk across parent
directories. When ceiling is 0, the lookup is only in
path. When ceiling is 1, the lookup is in both the path and
the parent to path.}
}
\value{
Character vector with path (terminated by a file
    separator) to repository or NULL if this cannot be
    established.
}
\description{
Find path to repository for any file
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "example-1.txt"))
add(repo, "example-1.txt")
commit(repo, "First commit message")

## Create a second file. The file is not added for version control
## in the repository.
dir.create(file.path(path, "example"))
file_2 <- file.path(path, "example/example-2.txt")
writeLines("Not under version control", file_2)

## Find the path to the repository using the path to the second file
discover_repository(file_2)

## Demonstrate the 'ceiling' argument
wd <- workdir(repo)
dir.create(file.path(wd, "temp"))

## Lookup repository in 'file.path(wd, "temp")'. Should return NULL
discover_repository(file.path(wd, "temp"), ceiling = 0)

## Lookup repository in parent to 'file.path(wd, "temp")'.
## Should not return NULL
discover_repository(file.path(wd, "temp"), ceiling = 1)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{repository_head}
\alias{repository_head}
\title{Get HEAD for a repository}
\usage{
repository_head(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
NULL if unborn branch or not found. A git_branch if not a
    detached head. A git_commit if detached head
}
\description{
Get HEAD for a repository
}
\examples{
\dontrun{
## Create and initialize a repository in a temporary directory
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "Commit message")

## Get HEAD of repository
repository_head(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stash.R
\name{summary.git_stash}
\alias{summary.git_stash}
\title{Summary of a stash}
\usage{
\method{summary}{git_stash}(object, ...)
}
\arguments{
\item{object}{The stash \code{object}}

\item{...}{Additional arguments affecting the summary produced.}
}
\value{
None (invisible 'NULL').
}
\description{
Summary of a stash
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

# Configure a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

# Create a file, add and commit
writeLines("Hello world!", file.path(path, "test.txt"))
add(repo, 'test.txt')
commit(repo, "Commit message")

# Change file
writeLines(c("Hello world!", "HELLO WORLD!"), file.path(path, "test.txt"))

# Create stash in repository
stash(repo, "Stash message")

# View summary of stash
summary(stash_list(repo)[[1]])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remote.R
\name{remote_set_url}
\alias{remote_set_url}
\title{Set the remote's url in the configuration}
\usage{
remote_set_url(repo = ".", name = NULL, url = NULL)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{name}{The name of the remote}

\item{url}{The \code{url} to set}
}
\value{
NULL, invisibly
}
\description{
This assumes the common case of a single-url remote and will
otherwise raise an error.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name="Alice", user.email="alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Add a remote
remote_add(repo, "playground", "https://example.org/git2r/playground")
remotes(repo)
remote_url(repo, "playground")

## Rename a remote
remote_rename(repo, "playground", "foobar")
remotes(repo)
remote_url(repo, "foobar")

## Set remote url
remote_set_url(repo, "foobar", "https://example.org/git2r/foobar")
remotes(repo)
remote_url(repo, "foobar")

## Remove a remote
remote_remove(repo, "foobar")
remotes(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/note.R
\name{note_create}
\alias{note_create}
\title{Add note for a object}
\usage{
note_create(
  object = NULL,
  message = NULL,
  ref = NULL,
  author = NULL,
  committer = NULL,
  force = FALSE
)
}
\arguments{
\item{object}{The object to annotate (git_blob, git_commit or
git_tree).}

\item{message}{Content of the note to add}

\item{ref}{Canonical name of the reference to use. Default is
\code{note_default_ref}.}

\item{author}{Signature of the notes note author}

\item{committer}{Signature of the notes note committer}

\item{force}{Overwrite existing note. Default is FALSE}
}
\value{
git_note
}
\description{
Add note for a object
}
\examples{
\dontrun{
## Create and initialize a repository in a temporary directory
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit_1 <- commit(repo, "Commit message 1")

## Create another commit
writeLines(c("Hello world!",
             "HELLO WORLD!"),
           file.path(path, "example.txt"))
add(repo, "example.txt")
commit_2 <- commit(repo, "Commit message 2")

## Check that notes is an empty list
notes(repo)

## Create note in default namespace
note_create(commit_1, "Note-1")

## Create note in named (review) namespace
note_create(commit_1, "Note-2", ref="refs/notes/review")
note_create(commit_2, "Note-3", ref="review")

## Create note on blob and tree
note_create(tree(commit_1), "Note-4")
note_create(tree(commit_1)["example.txt"], "Note-5")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge.R
\name{merge_base}
\alias{merge_base}
\title{Find a merge base between two commits}
\usage{
merge_base(one = NULL, two = NULL)
}
\arguments{
\item{one}{One of the commits}

\item{two}{The other commit}
}
\value{
git_commit
}
\description{
Find a merge base between two commits
}
\examples{
\dontrun{
## Create a directory in tempdir
path <- tempfile(pattern="git2r-")
dir.create(path)

## Initialize a repository
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
writeLines("Master branch", file.path(path, "master_branch.txt"))
add(repo, "master_branch.txt")
commit_1 <- commit(repo, "Commit message 1")

## Create first branch, checkout, add file and commit
branch_1 <- branch_create(commit_1, "branch_1")
checkout(branch_1)
writeLines("Branch 1", file.path(path, "branch_1.txt"))
add(repo, "branch_1.txt")
commit_2 <- commit(repo, "Commit message branch_1")

## Create second branch, checkout, add file and commit
branch_2 <- branch_create(commit_1, "branch_2")
checkout(branch_2)
writeLines("Branch 2", file.path(path, "branch_2.txt"))
add(repo, "branch_2.txt")
commit_3 <- commit(repo, "Commit message branch_2")

## Check that merge base equals commit_1
stopifnot(identical(merge_base(commit_2, commit_3), commit_1))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree.R
\name{is_tree}
\alias{is_tree}
\title{Check if object is S3 class git_tree}
\usage{
is_tree(object)
}
\arguments{
\item{object}{Check if object is S3 class git_tree}
}
\value{
TRUE if object is S3 class git_tree, else FALSE
}
\description{
Check if object is S3 class git_tree
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Commit a text file
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit_1 <- commit(repo, "First commit message")
tree_1 <- tree(commit_1)

## Check if tree
is_tree(commit_1)
is_tree(tree_1)
}
}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{summary.git_repository}
\alias{summary.git_repository}
\title{Summary of repository}
\usage{
\method{summary}{git_repository}(object, ...)
}
\arguments{
\item{object}{The repository \code{object}}

\item{...}{Additional arguments affecting the summary produced.}
}
\value{
None (invisible 'NULL').
}
\description{
Summary of repository
}
\examples{
\dontrun{
## Initialize a repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Config user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file
writeLines("Hello world!", file.path(path, "test.txt"))
summary(repo)

## Add file
add(repo, "test.txt")
summary(repo)

## Commit
commit(repo, "First commit message")
summary(repo)

## Change the file
writeLines(c("Hello again!", "Here is a second line", "And a third"),
           file.path(path, "test.txt"))
summary(repo)

## Add file and commit
add(repo, "test.txt")
commit(repo, "Second commit message")
summary(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blob.R
\name{hashfile}
\alias{hashfile}
\title{Determine the sha from a blob in a file}
\usage{
hashfile(path = NULL)
}
\arguments{
\item{path}{The path vector with files to hash.}
}
\value{
A vector with the sha for each file in path.
}
\description{
The blob is not written to the object database.
}
\examples{
\dontrun{
## Create a file. NOTE: The line endings from writeLines gives
## LF (line feed) on Unix/Linux and CRLF (carriage return, line feed)
## on Windows. The example use writeChar to have more control.
path <- tempfile()
f <- file(path, "wb")
writeChar("Hello, world!\n", f, eos = NULL)
close(f)

## Generate hash
hashfile(path)
identical(hashfile(path), hash("Hello, world!\n"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stash.R
\name{stash_list}
\alias{stash_list}
\title{List stashes in repository}
\usage{
stash_list(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
list of stashes in repository
}
\description{
List stashes in repository
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

# Configure a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

# Create a file, add and commit
writeLines("Hello world!", file.path(path, "test-1.txt"))
add(repo, 'test-1.txt')
commit(repo, "Commit message")

# Make one more commit
writeLines(c("Hello world!", "HELLO WORLD!"), file.path(path, "test-1.txt"))
add(repo, 'test-1.txt')
commit(repo, "Next commit message")

# Create one more file
writeLines("Hello world!", file.path(path, "test-2.txt"))

# Check that there are no stashes
stash_list(repo)

# Stash
stash(repo)

# Only untracked changes, therefore no stashes
stash_list(repo)

# Stash and include untracked changes
stash(repo, "Stash message", untracked=TRUE)

# View stash
stash_list(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{lookup_commit}
\alias{lookup_commit}
\alias{lookup_commit.git_branch}
\alias{lookup_commit.git_commit}
\alias{lookup_commit.git_tag}
\alias{lookup_commit.git_reference}
\title{Lookup the commit related to a git object}
\usage{
lookup_commit(object)

\method{lookup_commit}{git_branch}(object)

\method{lookup_commit}{git_commit}(object)

\method{lookup_commit}{git_tag}(object)

\method{lookup_commit}{git_reference}(object)
}
\arguments{
\item{object}{a git object to get the related commit from.}
}
\value{
A git commit object.
}
\description{
Lookup the commit related to a git_reference, git_tag or
git_branch object.
}
\examples{
\dontrun{
## Create a directory in tempdir
path <- tempfile(pattern="git2r-")
dir.create(path)

## Initialize a repository
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, con = file.path(path, "test.txt"))
add(repo, "test.txt")
commit(repo, "Commit message 1")

## Get the commit pointed to by the 'master' branch
lookup_commit(repository_head(repo))

## Create a tag
a_tag <- tag(repo, "Tagname", "Tag message")

## Get the commit pointed to by 'a_tag'
lookup_commit(a_tag)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reference.R
\name{references}
\alias{references}
\title{Get all references that can be found in a repository.}
\usage{
references(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
Character vector with references
}
\description{
Get all references that can be found in a repository.
}
\examples{
\dontrun{
## Initialize two temporary repositories
path_bare <- tempfile(pattern="git2r-")
path_repo <- tempfile(pattern="git2r-")
dir.create(path_bare)
dir.create(path_repo)
repo_bare <- init(path_bare, bare = TRUE)
repo <- clone(path_bare, path_repo)

## Config user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path_repo, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Push commits from repository to bare repository
## Adds an upstream tracking branch to branch 'master'
push(repo, "origin", "refs/heads/master")

## Add tag to HEAD
tag(repo, "v1.0", "First version")

## Create a note
note_create(commits(repo)[[1]], "My note")

## List all references in repository
references(repo)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree.R
\name{tree}
\alias{tree}
\title{Tree}
\usage{
tree(object = NULL)
}
\arguments{
\item{object}{the \code{commit} or \code{stash} object}
}
\value{
A S3 class git_tree object
}
\description{
Get the tree pointed to by a commit or stash.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a first user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

tree(last_commit(repo))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reset.R
\name{reset}
\alias{reset}
\title{Reset current HEAD to the specified state}
\usage{
reset(object, reset_type = c("soft", "mixed", "hard"), path = NULL)
}
\arguments{
\item{object}{Either a \code{git_commit}, a \code{git_repository}
or a character vector. If \code{object} is a
\code{git_commit}, HEAD is moved to the \code{git_commit}. If
\code{object} is a \code{git_repository}, resets the index
entries in the \code{path} argument to their state at HEAD. If
\code{object} is a character vector with paths, resets the
index entries in \code{object} to their state at HEAD if the
current working directory is in a repository.}

\item{reset_type}{If object is a 'git_commit', the kind of reset
operation to perform. 'soft' means the HEAD will be moved to
the commit. 'mixed' reset will trigger a 'soft' reset, plus
the index will be replaced with the content of the commit
tree. 'hard' reset will trigger a 'mixed' reset and the
working directory will be replaced with the content of the
index.}

\item{path}{If object is a 'git_repository', resets the index
entries for all paths to their state at HEAD.}
}
\value{
invisible NULL
}
\description{
Reset current HEAD to the specified state
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

# Configure a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
writeLines("Hello world!", file.path(path, "test-1.txt"))
add(repo, "test-1.txt")
commit_1 <- commit(repo, "Commit message")

## Change and stage the file
writeLines(c("Hello world!", "HELLO WORLD!"), file.path(path, "test-1.txt"))
add(repo, "test-1.txt")
status(repo)

## Unstage file
reset(repo, path = "test-1.txt")
status(repo)

## Make one more commit
add(repo, "test-1.txt")
commit(repo, "Next commit message")

## Create one more file
writeLines("Hello world!", file.path(path, "test-2.txt"))

## 'soft' reset to first commit and check status
reset(commit_1)
status(repo)

## 'mixed' reset to first commit and check status
commit(repo, "Next commit message")
reset(commit_1, "mixed")
status(repo)

## 'hard' reset to first commit and check status
add(repo, "test-1.txt")
commit(repo, "Next commit message")
reset(commit_1, "hard")
status(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blob.R
\name{is_blob}
\alias{is_blob}
\title{Check if object is S3 class git_blob}
\usage{
is_blob(object)
}
\arguments{
\item{object}{Check if object is S3 class git_blob}
}
\value{
TRUE if object is S3 class git_blob, else FALSE
}
\description{
Check if object is S3 class git_blob
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Commit a text file
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit_1 <- commit(repo, "First commit message")
blob_1 <- tree(commit_1)["example.txt"]

## Check if blob
is_blob(commit_1)
is_blob(blob_1)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stash.R
\name{stash_drop}
\alias{stash_drop}
\title{Drop stash}
\usage{
stash_drop(object = ".", index = 1)
}
\arguments{
\item{object}{path to a repository, or a \code{git_repository}
object, or the stash \code{object} to drop. Default is a
\code{path = '.'} to a reposiory.}

\item{index}{The index to the stash to drop. Only used when
\code{object} is a path to a repository or a
\code{git_repository} object. Default is \code{index = 1}.}
}
\value{
invisible NULL
}
\description{
Drop stash
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

# Configure a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

# Create a file, add and commit
writeLines("Hello world!", file.path(path, "test.txt"))
add(repo, 'test.txt')
commit(repo, "Commit message")

# Change file
writeLines(c("Hello world!", "HELLO WORLD!"), file.path(path, "test.txt"))

# Create stash in repository
stash(repo)

# Change file
writeLines(c("Hello world!", "HeLlO wOrLd!"), file.path(path, "test.txt"))

# Create stash in repository
stash(repo)

# View stashes
stash_list(repo)

# Drop git_stash object in repository
stash_drop(stash_list(repo)[[1]])

## Drop stash using an index to stash
stash_drop(repo, 1)

# View stashes
stash_list(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch.R
\name{is_head}
\alias{is_head}
\title{Check if branch is head}
\usage{
is_head(branch = NULL)
}
\arguments{
\item{branch}{The branch \code{object} to check if it's head.}
}
\value{
\code{TRUE} if branch is head, else \code{FALSE}.
}
\description{
Check if branch is head
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## List branches
branches(repo)

## Check that 'master' is_head
master <- branches(repo)[[1]]
is_head(master)

## Create and checkout 'dev' branch
checkout(repo, "dev", create = TRUE)

## List branches
branches(repo)

## Check that 'master' is no longer head
is_head(master)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree.R
\name{length.git_tree}
\alias{length.git_tree}
\title{Number of entries in tree}
\usage{
\method{length}{git_tree}(x)
}
\arguments{
\item{x}{The tree \code{object}}
}
\value{
a non-negative integer or double (which will be rounded
down)
}
\description{
Number of entries in tree
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/when.R
\name{when}
\alias{when}
\title{When}
\usage{
when(object, tz = "GMT", origin = "1970-01-01", usetz = TRUE)
}
\arguments{
\item{object}{the \code{object} to extract the time slot from.}

\item{tz}{time zone specification to be used for the conversion,
    \emph{if one is required}.  System-specific (see \link[base]{time zones}),
    but \code{""} is the current time zone, and \code{"GMT"} is UTC
    (Universal Time, Coordinated).  Invalid values are most commonly
    treated as UTC, on some platforms with a warning.}

\item{origin}{a date-time object, or something which can be coerced by
    \code{as.POSIXct(tz = "GMT")} to such an object.}

\item{usetz}{logical.  Should the time zone abbreviation be appended
    to the output?  This is used in printing times, and more reliable
    than using \code{"\%Z"}.}
}
\value{
A \code{character} vector of length one.
}
\description{
Help method to extract the time as a character string from a
git_commit, git_signature, git_tag and git_time object.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a first user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Create tag
tag(repo, "Tagname", "Tag message")

when(commits(repo)[[1]])
when(tags(repo)[[1]])
when(tags(repo)[[1]], tz = Sys.timezone())
}
}
\seealso{
\code{\link{git_time}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/credential.R
\name{cred_user_pass}
\alias{cred_user_pass}
\title{Create a new plain-text username and password credential object}
\usage{
cred_user_pass(username = NULL, password = NULL)
}
\arguments{
\item{username}{The username of the credential}

\item{password}{The password of the credential. If getPass is installed
and the only input is username, \code{getPass::getPass()} will be
called to allow for interactive and obfuscated interactive
input of the password.}
}
\value{
A list of class \code{cred_user_pass} with entries:
\describe{
  \item{username}{
    The username of the credential
  }
  \item{password}{
    The password of the credential
  }
}
}
\description{
Create a new plain-text username and password credential object
}
\examples{
\dontrun{
## Create a plain-text username and password credential object
cred_user_pass("Random Developer", "SecretPassword")
}
}
\seealso{
Other git credential functions: 
\code{\link{cred_env}()},
\code{\link{cred_ssh_key}()},
\code{\link{cred_token}()}
}
\concept{git credential functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remote.R
\name{remotes}
\alias{remotes}
\title{Get the configured remotes for a repo}
\usage{
remotes(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
Character vector with remotes
}
\description{
Get the configured remotes for a repo
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name="Alice", user.email="alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Add a remote
remote_add(repo, "playground", "https://example.org/git2r/playground")
remotes(repo)
remote_url(repo, "playground")

## Rename a remote
remote_rename(repo, "playground", "foobar")
remotes(repo)
remote_url(repo, "foobar")

## Set remote url
remote_set_url(repo, "foobar", "https://example.org/git2r/foobar")
remotes(repo)
remote_url(repo, "foobar")

## Remove a remote
remote_remove(repo, "foobar")
remotes(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{config}
\alias{config}
\title{Config}
\usage{
config(repo = NULL, global = FALSE, user.name, user.email, ...)
}
\arguments{
\item{repo}{The \code{repository}. Default is NULL.}

\item{global}{Write option(s) to global configuration
file. Default is FALSE.}

\item{user.name}{The user name. Use NULL to delete the entry}

\item{user.email}{The e-mail address. Use NULL to delete the entry}

\item{...}{Additional options to write or delete from the
configuration.}
}
\value{
S3 class \code{git_config}. When writing options, the
configuration is returned invisible.
}
\description{
Config file management. To display the configuration variables,
call method \code{config} without the \code{user.name},
\code{user.email} or \code{...} options.
}
\details{
There are two ways git2r can find the local repository when
writing local options (1) Use the \code{repo} argument. (2) If the
\code{repo} argument is \code{NULL} but the current working
directory is inside the local repository, then \code{git2r} uses
that repository.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern = "git2r-")
dir.create(path)
repo <- init(path)

## Set user name and email.
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Display configuration
config(repo)

## Delete user email.
config(repo, user.email = NULL)

## Display configuration
config(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{lookup}
\alias{lookup}
\title{Lookup}
\usage{
lookup(repo = ".", sha = NULL)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{sha}{The identity of the object to lookup. Must be 4 to 40
characters long.}
}
\value{
a \code{git_blob} or \code{git_commit} or \code{git_tag}
or \code{git_tree} object
}
\description{
Lookup one object in a repository.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit_1 <- commit(repo, "First commit message")

## Create tag
tag(repo, "Tagname", "Tag message")

## First, get SHAs to lookup in the repository
sha_commit <- sha(commit_1)
sha_tree <- sha(tree(commit_1))
sha_blob <- sha(tree(commit_1)["example.txt"])
sha_tag <- sha(tags(repo)[[1]])

## SHAs
sha_commit
sha_tree
sha_blob
sha_tag

## Lookup objects
lookup(repo, sha_commit)
lookup(repo, sha_tree)
lookup(repo, sha_blob)
lookup(repo, sha_tag)

## Lookup objects, using only the first seven characters
lookup(repo, substr(sha_commit, 1, 7))
lookup(repo, substr(sha_tree, 1, 7))
lookup(repo, substr(sha_blob, 1, 7))
lookup(repo, substr(sha_tag, 1, 7))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree.R
\name{as.list.git_tree}
\alias{as.list.git_tree}
\title{Coerce entries in a git_tree to a list of entry objects}
\usage{
\method{as.list}{git_tree}(x, ...)
}
\arguments{
\item{x}{The tree \code{object}}

\item{...}{Unused}
}
\value{
list of entry objects
}
\description{
Coerce entries in a git_tree to a list of entry objects
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
dir.create(file.path(path, "subfolder"))
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create three files and commit
writeLines("First file",  file.path(path, "example-1.txt"))
writeLines("Second file", file.path(path, "subfolder/example-2.txt"))
writeLines("Third file",  file.path(path, "example-3.txt"))
add(repo, c("example-1.txt", "subfolder/example-2.txt", "example-3.txt"))
commit(repo, "Commit message")

## Inspect size of each blob in tree
invisible(lapply(as(tree(last_commit(repo)), "list"),
  function(obj) {
    if (is_blob(obj))
      summary(obj)
    NULL
  }))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.git_repository}
\alias{plot.git_repository}
\title{Plot commits over time}
\usage{
\method{plot}{git_repository}(
  x,
  breaks = c("month", "year", "quarter", "week", "day"),
  main = NULL,
  ...
)
}
\arguments{
\item{x}{The repository to plot}

\item{breaks}{Default is \code{month}. Change to year, quarter,
week or day as necessary.}

\item{main}{Default title for the plot is "Commits on repo:" and
repository workdir basename. Supply a new title if you desire one.}

\item{...}{Additional arguments affecting the plot}
}
\description{
Plot commits over time
}
\examples{
\dontrun{
## Initialize repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- clone("https://github.com/ropensci/git2r.git", path)

## Plot commits
plot(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blob.R
\name{content}
\alias{content}
\title{Content of blob}
\usage{
content(blob = NULL, split = TRUE)
}
\arguments{
\item{blob}{The blob object.}

\item{split}{Split blob content to text lines. Default TRUE.}
}
\value{
The content of the blob. NA_character_ if the blob is binary.
}
\description{
Content of blob
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Display content of blob.
content(tree(commits(repo)[[1]])["example.txt"])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/libgit2.R
\name{libgit2_features}
\alias{libgit2_features}
\title{Compile time options for libgit2.}
\usage{
libgit2_features()
}
\value{
A list with threads, https and ssh set to TRUE/FALSE.
}
\description{
Compile time options for libgit2.
}
\examples{
libgit2_features()
}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/note.R
\name{note_remove}
\alias{note_remove}
\title{Remove the note for an object}
\usage{
note_remove(note = NULL, author = NULL, committer = NULL)
}
\arguments{
\item{note}{The note to remove}

\item{author}{Signature of the notes commit author.}

\item{committer}{Signature of the notes commit committer.}
}
\value{
invisible NULL
}
\description{
Remove the note for an object
}
\examples{
\dontrun{
## Create and initialize a repository in a temporary directory
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit_1 <- commit(repo, "Commit message 1")


## Create note in default namespace
note_1 <- note_create(commit_1, "Note-1")

## Create note in named (review) namespace
note_2 <- note_create(commit_1, "Note-2", ref="refs/notes/review")

## List notes in default namespace
notes(repo)

## List notes in 'review' namespace
notes(repo, "review")

## Remove notes
note_remove(note_1)
note_remove(note_2)

## List notes in default namespace
notes(repo)

## List notes in 'review' namespace
notes(repo, "review")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch.R
\name{branch_get_upstream}
\alias{branch_get_upstream}
\title{Get remote tracking branch}
\usage{
branch_get_upstream(branch = NULL)
}
\arguments{
\item{branch}{The branch}
}
\value{
\code{git_branch} object or NULL if no remote tracking
    branch.
}
\description{
Get remote tracking branch, given a local branch.
}
\examples{
\dontrun{
## Initialize two temporary repositories
path_bare <- tempfile(pattern="git2r-")
path_repo <- tempfile(pattern="git2r-")
dir.create(path_bare)
dir.create(path_repo)
repo_bare <- init(path_bare, bare = TRUE)
repo <- clone(path_bare, path_repo)

## Config user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path_repo, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Push commits from repository to bare repository
## Adds an upstream tracking branch to branch 'master'
push(repo, "origin", "refs/heads/master")

## Get remote tracking branch
branch_get_upstream(repository_head(repo))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blame.R
\name{blame}
\alias{blame}
\title{Get blame for file}
\usage{
blame(repo = ".", path = NULL)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{path}{Path to the file to consider}
}
\value{
git_blame object with the following entries:
\describe{
  \item{path}{
    The path to the file of the blame
  }
  \item{hunks}{
    List of blame hunks
  }
  \item{repo}{
    The git_repository that contains the file
  }
}
\describe{
  \item{lines_in_hunk}{
    The number of lines in this hunk
  }
  \item{final_commit_id}{
    The sha of the commit where this line was last changed
  }
  \item{final_start_line_number}{
    The 1-based line number where this hunk begins, in the final
    version of the file
  }
  \item{final_signature}{
    Final committer
  }
  \item{orig_commit_id}{
    The sha of the commit where this hunk was found. This will usually
    be the same as 'final_commit_id'.
  }
  \item{orig_start_line_number}{
     The 1-based line number where this hunk begins in the file
     named by 'orig_path' in the commit specified by 'orig_commit_id'.
  }
  \item{orig_signature}{
    Origin committer
  }
  \item{orig_path}{
    The path to the file where this hunk originated, as of the commit
    specified by 'orig_commit_id'
  }
  \item{boundary}{
    TRUE iff the hunk has been tracked to a boundary commit.
  }
  \item{repo}{
    The \code{git_repository} object that contains the blame hunk
  }
}
}
\description{
Get blame for file
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a first user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Create a second user and change the file
config(repo, user.name = "Bob", user.email = "bob@example.org")
writeLines(c("Hello world!", "HELLO WORLD!", "HOLA"),
           file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "Second commit message")

## Check blame
blame(repo, "example.txt")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commit.R
\name{commits}
\alias{commits}
\title{Commits}
\usage{
commits(
  repo = ".",
  topological = TRUE,
  time = TRUE,
  reverse = FALSE,
  n = NULL,
  ref = NULL,
  path = NULL
)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{topological}{Sort the commits in topological order (parents
before children); can be combined with time sorting. Default
is TRUE.}

\item{time}{Sort the commits by commit time; Can be combined with
topological sorting. Default is TRUE.}

\item{reverse}{Sort the commits in reverse order; can be combined
with topological and/or time sorting. Default is FALSE.}

\item{n}{The upper limit of the number of commits to output. The
default is NULL for unlimited number of commits.}

\item{ref}{The name of a reference to list commits from e.g. a tag
or a branch. The default is NULL for the current branch.}

\item{path}{The path to a file. If not NULL, only commits modifying
this file will be returned. Note that modifying commits that
occurred before the file was given its present name are not
returned; that is, the output of \code{git log} with
\code{--no-follow} is reproduced.}
}
\value{
list of commits in repository
}
\description{
Commits
}
\examples{
\dontrun{
## Initialize a repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Config user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Change file and commit
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "Second commit message")

## Create a tag
tag(repo, "Tagname", "Tag message")

## Change file again and commit
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad",
  "minim veniam, quis nostrud exercitation ullamco laboris nisi ut")
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "Third commit message")

## Create a new file containing R code, and commit.
writeLines(c("x <- seq(1,100)",
             "print(mean(x))"),
           file.path(path, "mean.R"))
add(repo, "mean.R")
commit(repo, "Fourth commit message")

## List the commits in the repository
commits(repo)

## List the commits starting from the tag
commits(repo, ref = "Tagname")

## List the commits modifying example.txt and mean.R.
commits(repo, path = "example.txt")
commits(repo, path = "mean.R")

## Create and checkout 'dev' branch in the repo
checkout(repo, "dev", create = TRUE)

## Add changes to the 'dev' branch
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "Commit message in dev branch")

## Checkout the 'master' branch again and list the commits
## starting from the 'dev' branch.
checkout(repo, "master")
commits(repo, ref = "dev")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/revparse.R
\name{revparse_single}
\alias{revparse_single}
\title{Revparse}
\usage{
revparse_single(repo = ".", revision = NULL)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{revision}{The revision string, see
http://git-scm.com/docs/git-rev-parse.html#_specifying_revisions}
}
\value{
a \code{git_commit} or \code{git_tag} or \code{git_tree}
object
}
\description{
Find object specified by revision.
}
\examples{
\dontrun{
## Create a directory in tempdir
path <- tempfile(pattern="git2r-")
dir.create(path)

## Initialize a repository
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "test.txt"))
add(repo, "test.txt")
commit(repo, "First commit message")

# Change file and commit
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path, "test.txt"))
add(repo, "test.txt")
commit(repo, "Second commit message")

revparse_single(repo, "HEAD^")
revparse_single(repo, "HEAD:test.txt")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch.R
\name{is_branch}
\alias{is_branch}
\title{Check if object is \code{git_branch}}
\usage{
is_branch(object)
}
\arguments{
\item{object}{Check if object is of class \code{git_branch}}
}
\value{
TRUE if object is class \code{git_branch}, else FALSE
}
\description{
Check if object is \code{git_branch}
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Commit a text file
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

branch <- branches(repo)[[1]]

## Check if branch
is_branch(branch)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/note.R
\name{notes}
\alias{notes}
\title{List notes}
\usage{
notes(repo = ".", ref = NULL)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{ref}{Reference to read from. Default (ref = NULL) is to call
\code{note_default_ref}.}
}
\value{
list with git_note objects
}
\description{
List all the notes within a specified namespace.
}
\examples{
\dontrun{
## Create and initialize a repository in a temporary directory
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit_1 <- commit(repo, "Commit message 1")

## Create another commit
writeLines(c("Hello world!",
             "HELLO WORLD!"),
           file.path(path, "example.txt"))
add(repo, "example.txt")
commit_2 <- commit(repo, "Commit message 2")

## Create note in default namespace
note_create(commit_1, "Note-1")
note_create(commit_1, "Note-2", force = TRUE)

## Create note in named (review) namespace
note_create(commit_1, "Note-3", ref="refs/notes/review")
note_create(commit_2, "Note-4", ref="review")

## Create note on blob and tree
note_create(tree(commit_1), "Note-5")
note_create(tree(commit_1)["example.txt"], "Note-6")

## List notes in default namespace
notes(repo)

## List notes in 'review' namespace
notes(repo, "review")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{repository}
\alias{repository}
\title{Open a repository}
\usage{
repository(path = ".", discover = TRUE)
}
\arguments{
\item{path}{A path to an existing local git repository.}

\item{discover}{Discover repository from path. Default is TRUE.}
}
\value{
A \code{git_repository} object with entries:
\describe{
  \item{path}{
    Path to a git repository
  }
}
}
\description{
Open a repository
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

# Configure a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
writeLines("Hello world!", file.path(path, "test-1.txt"))
add(repo, 'test-1.txt')
commit_1 <- commit(repo, "Commit message")

## Make one more commit
writeLines(c("Hello world!", "HELLO WORLD!"),
           file.path(path, "test-1.txt"))
add(repo, 'test-1.txt')
commit(repo, "Next commit message")

## Create one more file
writeLines("Hello world!",
           file.path(path, "test-2.txt"))

## Brief summary of repository
repo

## Summary of repository
summary(repo)

## Workdir of repository
workdir(repo)

## Check if repository is bare
is_bare(repo)

## Check if repository is empty
is_empty(repo)

## Check if repository is a shallow clone
is_shallow(repo)

## List all references in repository
references(repo)

## List all branches in repository
branches(repo)

## Get HEAD of repository
repository_head(repo)

## Check if HEAD is head
is_head(repository_head(repo))

## Check if HEAD is local
is_local(repository_head(repo))

## List all tags in repository
tags(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge.R
\name{merge.git_branch}
\alias{merge.git_branch}
\alias{merge.git_repository}
\alias{merge.character}
\title{Merge a branch into HEAD}
\usage{
\method{merge}{git_branch}(x, y = NULL, commit_on_success = TRUE, merger = NULL, fail = FALSE, ...)

\method{merge}{git_repository}(x, y = NULL, commit_on_success = TRUE, merger = NULL, fail = FALSE, ...)

\method{merge}{character}(
  x = ".",
  y = NULL,
  commit_on_success = TRUE,
  merger = NULL,
  fail = FALSE,
  ...
)
}
\arguments{
\item{x}{A path (default '.') to a repository, or a
\code{git_repository} object, or a \code{git_branch}.}

\item{y}{If \code{x} is a \code{git_repository}, the name of the
branch to merge into HEAD. Not used if \code{x} is a
\code{git_branch}.}

\item{commit_on_success}{If there are no conflicts written to the
index, the merge commit will be committed. Default is TRUE.}

\item{merger}{Who made the merge. The default (\code{NULL}) is to
use \code{default_signature} for the repository.}

\item{fail}{If a conflict occurs, exit immediately instead of
attempting to continue resolving conflicts. Default is
\code{FALSE}.}

\item{...}{Additional arguments (unused).}
}
\value{
A list of class \code{git_merge_result} with entries:
\describe{
  \item{up_to_date}{
    TRUE if the merge is already up-to-date, else FALSE.
  }
  \item{fast_forward}{
    TRUE if a fast-forward merge, else FALSE.
  }
  \item{conflicts}{
    TRUE if the index contain entries representing file conflicts,
    else FALSE.
  }
  \item{sha}{
    If the merge created a merge commit, the sha of the merge
    commit. NA if no merge commit created.
  }
}
}
\description{
Merge a branch into HEAD
}
\examples{
\dontrun{
## Create a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)
config(repo, user.name="Alice", user.email = "alice@example.org")

## Create a file, add and commit
writeLines("Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
           con = file.path(path, "test.txt"))
add(repo, "test.txt")
commit_1 <- commit(repo, "Commit message 1")

## Create first branch, checkout, add file and commit
checkout(repo, "branch1", create = TRUE)
writeLines("Branch 1", file.path(path, "branch-1.txt"))
add(repo, "branch-1.txt")
commit(repo, "Commit message branch 1")

## Create second branch, checkout, add file and commit
b_2 <- branch_create(commit_1, "branch2")
checkout(b_2)
writeLines("Branch 2", file.path(path, "branch-2.txt"))
add(repo, "branch-2.txt")
commit(repo, "Commit message branch 2")

## Make a change to 'test.txt'
writeLines(c("Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
             "eiusmod tempor incididunt ut labore et dolore magna aliqua."),
           con = file.path(path, "test.txt"))
add(repo, "test.txt")
commit(repo, "Second commit message branch 2")

## Checkout master
checkout(repo, "master", force = TRUE)

## Merge branch 1
merge(repo, "branch1")

## Merge branch 2
merge(repo, "branch2")

## Create third branch, checkout, change file and commit
checkout(repo, "branch3", create=TRUE)
writeLines(c("Lorem ipsum dolor amet sit, consectetur adipisicing elit, sed do",
             "eiusmod tempor incididunt ut labore et dolore magna aliqua."),
           con = file.path(path, "test.txt"))
add(repo, "test.txt")
commit(repo, "Commit message branch 3")

## Checkout master and create a change that creates a merge conflict
checkout(repo, "master", force=TRUE)
writeLines(c("Lorem ipsum dolor sit amet, adipisicing consectetur elit, sed do",
             "eiusmod tempor incididunt ut labore et dolore magna aliqua."),
           con = file.path(path, "test.txt"))
add(repo, "test.txt")
commit(repo, "Some commit message branch 1")

## Merge branch 3
merge(repo, "branch3")

## Check status; Expect to have one unstaged unmerged conflict.
status(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blob.R
\name{length.git_blob}
\alias{length.git_blob}
\title{Size in bytes of the contents of a blob}
\usage{
\method{length}{git_blob}(x)
}
\arguments{
\item{x}{The blob \code{object}}
}
\value{
a non-negative integer
}
\description{
Size in bytes of the contents of a blob
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Commit a text file
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit_1 <- commit(repo, "First commit message")
blob_1 <- tree(commit_1)["example.txt"]

## Get length in size of bytes of the content of the blob
length(blob_1)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remote.R
\name{remote_remove}
\alias{remote_remove}
\title{Remove a remote}
\usage{
remote_remove(repo = ".", name = NULL)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{name}{The name of the remote to remove}
}
\value{
NULL, invisibly
}
\description{
All remote-tracking branches and configuration settings for the
remote will be removed.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name="Alice", user.email="alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Add a remote
remote_add(repo, "playground", "https://example.org/git2r/playground")
remotes(repo)
remote_url(repo, "playground")

## Rename a remote
remote_rename(repo, "playground", "foobar")
remotes(repo)
remote_url(repo, "foobar")

## Set remote url
remote_set_url(repo, "foobar", "https://example.org/git2r/foobar")
remotes(repo)
remote_url(repo, "foobar")

## Remove a remote
remote_remove(repo, "foobar")
remotes(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{head.git_repository}
\alias{head.git_repository}
\title{Get HEAD for a repository}
\usage{
\method{head}{git_repository}(x, ...)
}
\arguments{
\item{x}{The repository \code{x} to check head}

\item{...}{Additional arguments. Unused.}
}
\value{
NULL if unborn branch or not found. A git_branch if not a
    detached head. A git_commit if detached head
}
\description{
Get HEAD for a repository
}
\examples{
\dontrun{
## Create and initialize a repository in a temporary directory
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "Commit message")

## Get HEAD of repository
repository_head(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remote.R
\name{remote_add}
\alias{remote_add}
\title{Add a remote to a repo}
\usage{
remote_add(repo = ".", name = NULL, url = NULL)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{name}{Short name of the remote repository}

\item{url}{URL of the remote repository}
}
\value{
NULL, invisibly
}
\description{
Add a remote to a repo
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name="Alice", user.email="alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Add a remote
remote_add(repo, "playground", "https://example.org/git2r/playground")
remotes(repo)
remote_url(repo, "playground")

## Rename a remote
remote_rename(repo, "playground", "foobar")
remotes(repo)
remote_url(repo, "foobar")

## Set remote url
remote_set_url(repo, "foobar", "https://example.org/git2r/foobar")
remotes(repo)
remote_url(repo, "foobar")

## Remove a remote
remote_remove(repo, "foobar")
remotes(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tag.R
\name{tags}
\alias{tags}
\title{Tags}
\usage{
tags(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
list of tags in repository
}
\description{
Tags
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Commit a text file
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Create tag
tag(repo, "Tagname", "Tag message")

## List tags
tags(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/time.R
\name{git_time}
\alias{git_time}
\alias{as.character.git_time}
\alias{format.git_time}
\alias{as.POSIXct.git_time}
\alias{print.git_time}
\title{Time}
\usage{
\method{as.character}{git_time}(x, tz = "GMT", origin = "1970-01-01", usetz = TRUE, ...)

\method{format}{git_time}(x, tz = "GMT", origin = "1970-01-01", usetz = TRUE, ...)

\method{as.POSIXct}{git_time}(x, tz = "GMT", origin = "1970-01-01", ...)

\method{print}{git_time}(x, tz = "GMT", origin = "1970-01-01", usetz = TRUE, ...)
}
\arguments{
\item{x}{\R  object to be converted.}

\item{tz}{time zone specification to be used for the conversion,
    \emph{if one is required}.  System-specific (see \link[base]{time zones}),
    but \code{""} is the current time zone, and \code{"GMT"} is UTC
    (Universal Time, Coordinated).  Invalid values are most commonly
    treated as UTC, on some platforms with a warning.}

\item{origin}{a date-time object, or something which can be coerced by
    \code{as.POSIXct(tz = "GMT")} to such an object.}

\item{usetz}{logical.  Should the time zone abbreviation be appended
    to the output?  This is used in printing times, and more reliable
    than using \code{"\%Z"}.}

\item{...}{further arguments to be passed to or from other methods.}
}
\description{
The class \code{git_time} stores the time a Git object was created.
}
\details{
The default is to use \code{tz = "GMT"} and \code{origin =
"1970-01-01"}. To use your local timezone, set \code{tz =
Sys.timezone()}.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a first user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Create tag
tag(repo, "Tagname", "Tag message")

as.POSIXct(commits(repo)[[1]]$author$when)
as.POSIXct(tags(repo)[[1]]$tagger$when)
as.POSIXct(tags(repo)[[1]]$tagger$when, tz = Sys.timezone())
}
}
\seealso{
\code{\link{when}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remote.R
\name{remote_url}
\alias{remote_url}
\title{Get the remote url for remotes in a repo}
\usage{
remote_url(repo = ".", remote = NULL)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{remote}{Character vector with the remotes to get the url
from. Default is the remotes of the repository.}
}
\value{
Character vector with remote_url for each of the remote
}
\description{
Get the remote url for remotes in a repo
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name="Alice", user.email="alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Add a remote
remote_add(repo, "playground", "https://example.org/git2r/playground")
remotes(repo)
remote_url(repo, "playground")

## Rename a remote
remote_rename(repo, "playground", "foobar")
remotes(repo)
remote_url(repo, "foobar")

## Set remote url
remote_set_url(repo, "foobar", "https://example.org/git2r/foobar")
remotes(repo)
remote_url(repo, "foobar")

## Remove a remote
remote_remove(repo, "foobar")
remotes(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{length.git_diff}
\alias{length.git_diff}
\title{Number of files in git_diff object}
\usage{
\method{length}{git_diff}(x)
}
\arguments{
\item{x}{The git_diff \code{object}}
}
\value{
a non-negative integer
}
\description{
Number of files in git_diff object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{clone}
\alias{clone}
\title{Clone a remote repository}
\usage{
clone(
  url = NULL,
  local_path = NULL,
  bare = FALSE,
  branch = NULL,
  checkout = TRUE,
  credentials = NULL,
  progress = TRUE
)
}
\arguments{
\item{url}{The remote repository to clone}

\item{local_path}{Local directory to clone to.}

\item{bare}{Create a bare repository. Default is FALSE.}

\item{branch}{The name of the branch to checkout. Default is NULL
which means to use the remote's default branch.}

\item{checkout}{Checkout HEAD after the clone is complete. Default
is TRUE.}

\item{credentials}{The credentials for remote repository
access. Default is NULL. To use and query an ssh-agent for the
ssh key credentials, let this parameter be NULL (the default).}

\item{progress}{Show progress. Default is TRUE.}
}
\value{
A \code{git_repository} object.
}
\description{
Clone a remote repository
}
\examples{
\dontrun{
## Initialize repository
path_repo_1 <- tempfile(pattern="git2r-")
path_repo_2 <- tempfile(pattern="git2r-")
dir.create(path_repo_1)
dir.create(path_repo_2)
repo_1 <- init(path_repo_1)

## Config user and commit a file
config(repo_1, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
writeLines(
    "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
    file.path(path_repo_1, "example.txt"))
add(repo_1, "example.txt")
commit(repo_1, "First commit message")

## Change file and commit
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path_repo_1, "example.txt"))
add(repo_1, "example.txt")
commit(repo_1, "Second commit message")

## Change file again and commit.
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad",
  "minim veniam, quis nostrud exercitation ullamco laboris nisi ut")
writeLines(lines, file.path(path_repo_1, "example.txt"))
add(repo_1, "example.txt")
commit(repo_1, "Third commit message")

## Clone to second repository
repo_2 <- clone(path_repo_1, path_repo_2)

## List commits in repositories
commits(repo_1)
commits(repo_2)
}
}
\seealso{
\link{repository}, \code{\link{cred_user_pass}},
    \code{\link{cred_ssh_key}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{diff.git_repository}
\alias{diff.git_repository}
\alias{diff.git_tree}
\title{Changes between commits, trees, working tree, etc.}
\usage{
\method{diff}{git_repository}(
  x,
  index = FALSE,
  as_char = FALSE,
  filename = NULL,
  context_lines = 3,
  interhunk_lines = 0,
  old_prefix = "a",
  new_prefix = "b",
  id_abbrev = NULL,
  path = NULL,
  max_size = NULL,
  ...
)

\method{diff}{git_tree}(
  x,
  new_tree = NULL,
  index = FALSE,
  as_char = FALSE,
  filename = NULL,
  context_lines = 3,
  interhunk_lines = 0,
  old_prefix = "a",
  new_prefix = "b",
  id_abbrev = NULL,
  path = NULL,
  max_size = NULL,
  ...
)
}
\arguments{
\item{x}{A \code{git_repository} object or the old \code{git_tree}
object to compare to.}

\item{index}{\describe{
  \item{\emph{When object equals a git_repository}}{
    Whether to compare the index to HEAD. If FALSE (the default),
    then the working tree is compared to the index.
  }
  \item{\emph{When object equals a git_tree}}{
    Whether to use the working directory (by default), or the index
    (if set to TRUE) in the comparison to \code{object}.
  }
}}

\item{as_char}{logical: should the result be converted to a
character string?. Default is FALSE.}

\item{filename}{If as_char is TRUE, then the diff can be written
to a file with name filename (the file is overwritten if it
exists). Default is NULL.}

\item{context_lines}{The number of unchanged lines that define the
boundary of a hunk (and to display before and after). Defaults
to 3.}

\item{interhunk_lines}{The maximum number of unchanged lines
between hunk boundaries before the hunks will be merged into
one. Defaults to 0.}

\item{old_prefix}{The virtual "directory" prefix for old file
names in hunk headers. Default is "a".}

\item{new_prefix}{The virtual "directory" prefix for new file
names in hunk headers. Defaults to "b".}

\item{id_abbrev}{The abbreviation length to use when formatting
object ids. Defaults to the value of 'core.abbrev' from the
config, or 7 if NULL.}

\item{path}{A character vector of paths / fnmatch patterns to
constrain diff. Default is NULL which include all paths.}

\item{max_size}{A size (in bytes) above which a blob will be
marked as binary automatically; pass a negative value to
disable. Defaults to 512MB when max_size is NULL.}

\item{...}{Not used.}

\item{new_tree}{The new git_tree object to compare, or NULL.  If
NULL, then we use the working directory or the index (see the
\code{index} argument).}
}
\value{
A \code{git_diff} object if as_char is FALSE. If as_char
    is TRUE and filename is NULL, a character string, else NULL.
}
\description{
Changes between commits, trees, working tree, etc.
}
\examples{
\dontrun{
## Initialize a repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Config user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add, commit
writeLines("Hello world!", file.path(path, "test.txt"))
add(repo, "test.txt")
commit(repo, "Commit message")

## Change the file
writeLines(c("Hello again!", "Here is a second line", "And a third"),
           file.path(path, "test.txt"))

## diff between index and workdir
diff_1 <- diff(repo)
summary(diff_1)
cat(diff(repo, as_char=TRUE))

## Diff between index and HEAD is empty
diff_2 <- diff(repo, index=TRUE)
summary(diff_2)
cat(diff(repo, index=TRUE, as_char=TRUE))

## Diff between tree and working dir, same as diff_1
diff_3 <- diff(tree(commits(repo)[[1]]))
summary(diff_3)
cat(diff(tree(commits(repo)[[1]]), as_char=TRUE))

## Add changes, diff between index and HEAD is the same as diff_1
add(repo, "test.txt")
diff_4 <- diff(repo, index=TRUE)
summary(diff_4)
cat(diff(repo, index=TRUE, as_char=TRUE))

## Diff between tree and index
diff_5 <- diff(tree(commits(repo)[[1]]), index=TRUE)
summary(diff_5)
cat(diff(tree(commits(repo)[[1]]), index=TRUE, as_char=TRUE))

## Diff between two trees
commit(repo, "Second commit")
tree_1 <- tree(commits(repo)[[2]])
tree_2 <- tree(commits(repo)[[1]])
diff_6 <- diff(tree_1, tree_2)
summary(diff_6)
cat(diff(tree_1, tree_2, as_char=TRUE))

## Binary files
set.seed(42)
writeBin(as.raw((sample(0:255, 1000, replace=TRUE))),
         con=file.path(path, "test.bin"))
add(repo, "test.bin")
diff_7 <- diff(repo, index=TRUE)
summary(diff_7)
cat(diff(repo, index=TRUE, as_char=TRUE))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/index.R
\name{index_remove_bypath}
\alias{index_remove_bypath}
\title{Remove an index entry corresponding to a file on disk}
\usage{
index_remove_bypath(repo = ".", path = NULL)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{path}{character vector with filenames to remove. The path
must be relative to the repository's working folder. It may
exist. If this file currently is the result of a merge
conflict, this file will no longer be marked as
conflicting. The data about the conflict will be moved to the
"resolve undo" (REUC) section.}
}
\value{
invisible(NULL)
}
\description{
Remove an index entry corresponding to a file on disk
}
\examples{
\dontrun{
## Initialize a repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file
writeLines("Hello world!", file.path(path, "file-to-remove.txt"))

## Add file to repository
add(repo, "file-to-remove.txt")

## View status of repository
status(repo)

## Remove file
index_remove_bypath(repo, "file-to-remove.txt")

## View status of repository
status(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/status.R
\name{status}
\alias{status}
\title{Status}
\usage{
status(
  repo = ".",
  staged = TRUE,
  unstaged = TRUE,
  untracked = TRUE,
  ignored = FALSE,
  all_untracked = FALSE
)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{staged}{Include staged files. Default TRUE.}

\item{unstaged}{Include unstaged files. Default TRUE.}

\item{untracked}{Include untracked files and directories. Default
TRUE.}

\item{ignored}{Include ignored files. Default FALSE.}

\item{all_untracked}{Shows individual files in untracked
directories if \code{untracked} is \code{TRUE}.}
}
\value{
\code{git_status} with repository status
}
\description{
Display state of the repository working directory and the staging
area.
}
\examples{
\dontrun{
## Initialize a repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Config user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file
writeLines("Hello world!", file.path(path, "test.txt"))

## Check status; untracked file
status(repo)

## Add file
add(repo, "test.txt")

## Check status; staged file
status(repo)

## Commit
commit(repo, "First commit message")

## Check status; clean
status(repo)

## Change the file
writeLines(c("Hello again!", "Here is a second line", "And a third"),
           file.path(path, "test.txt"))

## Check status; unstaged file
status(repo)

## Add file and commit
add(repo, "test.txt")
commit(repo, "Second commit message")

## Check status; clean
status(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reflog.R
\name{reflog}
\alias{reflog}
\title{List and view reflog information}
\usage{
reflog(repo = ".", refname = "HEAD")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{refname}{The name of the reference to list. 'HEAD' by
default.}
}
\value{
S3 class \code{git_reflog} with git_reflog_entry objects.
}
\description{
List and view reflog information
}
\examples{
\dontrun{
## Initialize a repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Config user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Change file and commit
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "Second commit message")

## Change file again and commit
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad",
  "minim veniam, quis nostrud exercitation ullamco laboris nisi ut")
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "Third commit message")

## View reflog
reflog(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/credential.R
\name{cred_token}
\alias{cred_token}
\title{Create a new personal access token credential object}
\usage{
cred_token(token = "GITHUB_PAT")
}
\arguments{
\item{token}{The name of the environmental variable that holds the
personal access token for the authentication. Default is
\code{GITHUB_PAT}.}
}
\value{
A list of class \code{cred_token} with entry:
\describe{
  \item{token}{
    The name of the environmental variable that holds
    the personal access token for the authentication.
  }
}
}
\description{
The personal access token is stored in an envrionmental variable.
Environmental variables can be written to the file
\code{.Renviron}. This file is read by \emph{R} during startup,
see \code{\link[base]{Startup}}. On GitHub, personal access tokens
function like ordinary OAuth access tokens. They can be used
instead of a password for Git over HTTPS, see
\url{https://help.github.com/articles/creating-an-access-token-for-command-line-use}
}
\examples{
\dontrun{
## Create a personal access token credential object.
## This example assumes that the token is stored in
## the 'GITHUB_PAT' environmental variable.
repo <- repository("git2r")
cred <- cred_token()
push(repo, credentials = cred)
}
}
\seealso{
Other git credential functions: 
\code{\link{cred_env}()},
\code{\link{cred_ssh_key}()},
\code{\link{cred_user_pass}()}
}
\concept{git credential functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch.R
\name{branch_remote_name}
\alias{branch_remote_name}
\title{Remote name of a branch}
\usage{
branch_remote_name(branch = NULL)
}
\arguments{
\item{branch}{The branch}
}
\value{
character string with remote name
}
\description{
The name of remote that the remote tracking branch belongs to
}
\examples{
\dontrun{
## Initialize two temporary repositories
path_bare <- tempfile(pattern="git2r-")
path_repo <- tempfile(pattern="git2r-")
dir.create(path_bare)
dir.create(path_repo)
repo_bare <- init(path_bare, bare = TRUE)
repo <- clone(path_bare, path_repo)

## Config user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path_repo, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Push commits from repository to bare repository
## Adds an upstream tracking branch to branch 'master'
push(repo, "origin", "refs/heads/master")

## Get remote name
branch_remote_name(branches(repo)[[2]])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/credential.R
\name{ssh_path}
\alias{ssh_path}
\title{Compose usual path to ssh keys}
\usage{
ssh_path(file = "")
}
\arguments{
\item{file}{basename of file for which path is requested}
}
\value{
Full path to the file
}
\description{
This function provides a consistent means across OS-types to access the
\code{.ssh} directory.
}
\details{
On Windows-based systems,
\code{path.expand("~")} returns \code{"C:/Users/username/Documents"},
whereas the usual path to the \code{.ssh} directory is
\code{"C:/Users/username"}.

On other operating systems, \code{path.expand("~")} returns the usual path
to the \code{.ssh} directory.

Calling \code{ssh_path()} with no arguments will return the usual path to
the \code{.ssh} directory.
}
\examples{
ssh_path()
ssh_path("is_rsa.pub")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{is_shallow}
\alias{is_shallow}
\title{Determine if the repository is a shallow clone}
\usage{
is_shallow(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
\code{TRUE} if shallow clone, else \code{FALSE}
}
\description{
Determine if the repository is a shallow clone
}
\examples{
\dontrun{
## Initialize repository
path_repo_1 <- tempfile(pattern="git2r-")
path_repo_2 <- tempfile(pattern="git2r-")
dir.create(path_repo_1)
dir.create(path_repo_2)
repo_1 <- init(path_repo_1)

## Config user and commit a file
config(repo_1, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path_repo_1, "example.txt"))
add(repo_1, "example.txt")
commit(repo_1, "First commit message")

## Change file and commit
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path_repo_1, "example.txt"))
add(repo_1, "example.txt")
commit(repo_1, "Second commit message")

## Change file again and commit.
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad",
  "minim veniam, quis nostrud exercitation ullamco laboris nisi ut")
writeLines(lines, file.path(path_repo_1, "example.txt"))
add(repo_1, "example.txt")
commit(repo_1, "Third commit message")

## Clone to second repository
repo_2 <- clone(path_repo_1, path_repo_2)

## Check if it's a shallow clone
is_shallow(repo_2)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch.R
\name{branch_remote_url}
\alias{branch_remote_url}
\title{Remote url of a branch}
\usage{
branch_remote_url(branch = NULL)
}
\arguments{
\item{branch}{The branch}
}
\value{
character string with remote url
}
\description{
Remote url of a branch
}
\examples{
\dontrun{
## Initialize two temporary repositories
path_bare <- tempfile(pattern="git2r-")
path_repo <- tempfile(pattern="git2r-")
dir.create(path_bare)
dir.create(path_repo)
repo_bare <- init(path_bare, bare = TRUE)
repo <- clone(path_bare, path_repo)

## Config user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path_repo, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Push commits from repository to bare repository
## Adds an upstream tracking branch to branch 'master'
push(repo, "origin", "refs/heads/master")

## Get remote url of tracking branch to branch 'master'
branch_remote_url(branch_get_upstream(repository_head(repo)))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{in_repository}
\alias{in_repository}
\title{Determine if a directory is in a git repository}
\usage{
in_repository(path = ".")
}
\arguments{
\item{path}{The path to the directory.}
}
\value{
TRUE if directory is in a git repository else FALSE
}
\description{
The lookup start from path and walk across parent directories if
nothing has been found.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Check if path is in a git repository
in_repository(path)

## Check if working directory is in a git repository
setwd(path)
in_repository()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/punch_card.R
\name{punch_card}
\alias{punch_card}
\title{Punch card}
\usage{
punch_card(repo = ".", main = NULL, ...)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{main}{Default title for the plot is "Punch card on repo:"
and repository workdir basename. Supply a new title if you
desire one.}

\item{...}{Additional arguments affecting the plot}
}
\value{
invisible NULL
}
\description{
Punch card
}
\examples{
\dontrun{
## Initialize repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- clone("https://github.com/ropensci/git2r.git", path)

## Plot
punch_card(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree.R
\name{[.git_tree}
\alias{[.git_tree}
\title{Extract object from tree}
\usage{
\method{[}{git_tree}(x, i)
}
\arguments{
\item{x}{The tree \code{object}}

\item{i}{The index (integer or logical) of the tree object to
extract. If negative values, all elements except those indicated
are selected. A character vector to match against the names of
objects to extract.}
}
\value{
Git object
}
\description{
Lookup a tree entry by its position in the tree
}
\examples{
\dontrun{
##' Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
dir.create(file.path(path, "subfolder"))
repo <- init(path)

##' Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

##' Create three files and commit
writeLines("First file",  file.path(path, "example-1.txt"))
writeLines("Second file", file.path(path, "subfolder/example-2.txt"))
writeLines("Third file",  file.path(path, "example-3.txt"))
add(repo, c("example-1.txt", "subfolder/example-2.txt", "example-3.txt"))
new_commit <- commit(repo, "Commit message")

##' Pick a tree in the repository
tree_object <- tree(new_commit)

##' Display tree
tree_object

##' Select item by name
tree_object["example-1.txt"]

##' Select first item in tree
tree_object[1]

##' Select first three items in tree
tree_object[1:3]

##' Select all blobs in tree
tree_object[vapply(as(tree_object, 'list'), is_blob, logical(1))]
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch.R
\name{is_local}
\alias{is_local}
\title{Check if branch is local}
\usage{
is_local(branch)
}
\arguments{
\item{branch}{The branch \code{object} to check if it's local}
}
\value{
\code{TRUE} if branch is local, else \code{FALSE}.
}
\description{
Check if branch is local
}
\examples{
\dontrun{
## Initialize repositories
path_bare <- tempfile(pattern="git2r-")
path_repo <- tempfile(pattern="git2r-")
dir.create(path_bare)
dir.create(path_repo)
repo_bare <- init(path_bare, bare = TRUE)
repo <- clone(path_bare, path_repo)

## Config first user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path_repo, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Push commits from repository to bare repository
## Adds an upstream tracking branch to branch 'master'
push(repo, "origin", "refs/heads/master")

## List branches
branches(repo)

## Check if first branch is_local
is_local(branches(repo)[[1]])

## Check if second branch is_local
is_local(branches(repo)[[2]])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commit.R
\name{is_merge}
\alias{is_merge}
\title{Is merge}
\usage{
is_merge(commit = NULL)
}
\arguments{
\item{commit}{a git_commit object.}
}
\value{
TRUE if commit has more than one parent, else FALSE
}
\description{
Determine if a commit is a merge commit, i.e. has more than one
parent.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
writeLines(c("First line in file 1.", "Second line in file 1."),
           file.path(path, "example-1.txt"))
add(repo, "example-1.txt")
commit(repo, "First commit message")

## Create and add one more file
writeLines(c("First line in file 2.", "Second line in file 2."),
           file.path(path, "example-2.txt"))
add(repo, "example-2.txt")
commit(repo, "Second commit message")

## Create a new branch 'fix'
checkout(repo, "fix", create = TRUE)

## Update 'example-1.txt' (swap words in first line) and commit
writeLines(c("line First in file 1.", "Second line in file 1."),
           file.path(path, "example-1.txt"))
add(repo, "example-1.txt")
commit(repo, "Third commit message")

checkout(repo, "master")

## Update 'example-2.txt' (swap words in second line) and commit
writeLines(c("First line in file 2.", "line Second in file 2."),
           file.path(path, "example-2.txt"))
add(repo, "example-2.txt")
commit(repo, "Fourth commit message")

## Merge 'fix'
merge(repo, "fix")

## Display parents of last commit
parents(lookup(repo, branch_target(repository_head(repo))))

## Check that last commit is a merge
is_merge(lookup(repo, branch_target(repository_head(repo))))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{is_bare}
\alias{is_bare}
\title{Check if repository is bare}
\usage{
is_bare(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
\code{TRUE} if bare repository, else \code{FALSE}
}
\description{
Check if repository is bare
}
\examples{
\dontrun{
## Initialize a repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)
is_bare(repo)

## Initialize a bare repository
path_bare <- tempfile(pattern="git2r-")
dir.create(path_bare)
repo_bare <- init(path_bare, bare = TRUE)
is_bare(repo_bare)
}
}
\seealso{
\link{init}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tag.R
\name{tag_delete}
\alias{tag_delete}
\title{Delete an existing tag reference}
\usage{
tag_delete(object = ".", name = NULL)
}
\arguments{
\item{object}{Can be either the path (default is ".") to a
repository, or a \code{git_repository} object, or a
\code{git_tag} object. or the tag name.}

\item{name}{If the \code{object} argument is a path to a
repository or a \code{git_repository}, the name of the tag to
delete.}
}
\value{
\code{invisible(NULL)}
}
\description{
Delete an existing tag reference
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Commit a text file
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Create two tags
tag(repo, "Tag1", "Tag message 1")
t2 <- tag(repo, "Tag2", "Tag message 2")

## List the two tags in the repository
tags(repo)

## Delete the two tags in the repository
tag_delete(repo, "Tag1")
tag_delete(t2)

## Show the empty list with tags in the repository
tags(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/libgit2.R
\name{ssl_cert_locations}
\alias{ssl_cert_locations}
\title{Set the SSL certificate-authority locations}
\usage{
ssl_cert_locations(filename = NULL, path = NULL)
}
\arguments{
\item{filename}{Location of a file containing several certificates
concatenated together. Default NULL.}

\item{path}{Location of a directory holding several certificates,
one per file. Default NULL.}
}
\value{
invisible(NULL)
}
\description{
Set the SSL certificate-authority locations
}
\note{
Either parameter may be 'NULL', but not both.
}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odb.R
\name{odb_blobs}
\alias{odb_blobs}
\title{Blobs in the object database}
\usage{
odb_blobs(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
A data.frame with the following columns:
\describe{
  \item{sha}{The sha of the blob}
  \item{path}{The path to the blob from the tree and sub-trees}
  \item{name}{The name of the blob from the tree that contains the blob}
  \item{len}{The length of the blob}
  \item{commit}{The sha of the commit}
  \item{author}{The author of the commit}
  \item{when}{The timestamp of the author signature in the commit}
}
}
\description{
List all blobs reachable from the commits in the object
database. For each commit, list blob's in the commit tree and
sub-trees.
}
\note{
A blob sha can have several entries
}
\examples{
\dontrun{
## Create a directory in tempdir
path <- tempfile(pattern="git2r-")
dir.create(path)

## Initialize a repository
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "test.txt"))
add(repo, "test.txt")
commit(repo, "Commit message 1")

## Change file and commit
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path, "test.txt"))
add(repo, "test.txt")
commit(repo, "Commit message 2")

## Commit same content under different name in a sub-directory
dir.create(file.path(path, "sub-directory"))
file.copy(file.path(path, "test.txt"),
          file.path(path, "sub-directory", "copy.txt"))
add(repo, "sub-directory/copy.txt")
commit(repo, "Commit message 3")

## List blobs
odb_blobs(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pull.R
\name{pull}
\alias{pull}
\title{Pull}
\usage{
pull(repo = ".", credentials = NULL, merger = NULL)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{credentials}{The credentials for remote repository
access. Default is NULL. To use and query an ssh-agent for the
ssh key credentials, let this parameter be NULL (the default).}

\item{merger}{Who made the merge, if the merge is non-fast forward
merge that creates a merge commit. The
\code{default_signature} for \code{repo} is used if this
parameter is \code{NULL}.}
}
\value{
A list of class \code{git_merge_result} with entries:
\describe{
  \item{up_to_date}{
    TRUE if the merge is already up-to-date, else FALSE.
  }
  \item{fast_forward}{
    TRUE if a fast-forward merge, else FALSE.
  }
  \item{conflicts}{
    TRUE if the index contain entries representing file conflicts,
    else FALSE.
  }
  \item{sha}{
    If the merge created a merge commit, the sha of the merge
    commit. NA if no merge commit created.
  }
}
}
\description{
Pull
}
\examples{
\dontrun{
## Initialize repositories
path_bare <- tempfile(pattern="git2r-")
path_repo_1 <- tempfile(pattern="git2r-")
path_repo_2 <- tempfile(pattern="git2r-")
dir.create(path_bare)
dir.create(path_repo_1)
dir.create(path_repo_2)
repo_bare <- init(path_bare, bare = TRUE)
repo_1 <- clone(path_bare, path_repo_1)

## Config first user and commit a file
config(repo_1, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path_repo_1, "example.txt"))
add(repo_1, "example.txt")
commit(repo_1, "First commit message")

## Push commits from first repository to bare repository
## Adds an upstream tracking branch to branch 'master'
push(repo_1, "origin", "refs/heads/master")

## Clone to second repository
repo_2 <- clone(path_bare, path_repo_2)
config(repo_2, user.name = "Bob", user.email = "bob@example.org")

## Change file and commit
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path_repo_1, "example.txt"))
add(repo_1, "example.txt")
commit(repo_1, "Second commit message")

## Push commits from first repository to bare repository
push(repo_1)

## Pull changes to repo_2
pull(repo_2)

## Change file again and commit. This time in repository 2
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad",
  "minim veniam, quis nostrud exercitation ullamco laboris nisi ut")
writeLines(lines, file.path(path_repo_2, "example.txt"))
add(repo_2, "example.txt")
commit(repo_2, "Third commit message")

## Push commits from second repository to bare repository
push(repo_2)

## Pull changes to repo_1
pull(repo_1)

## List commits in repositories
commits(repo_1)
commits(repo_2)
commits(repo_bare)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch.R
\name{branch_rename}
\alias{branch_rename}
\title{Rename a branch}
\usage{
branch_rename(branch = NULL, name = NULL, force = FALSE)
}
\arguments{
\item{branch}{Branch to rename}

\item{name}{The new name for the branch}

\item{force}{Overwrite existing branch. Default is FALSE}
}
\value{
invisible renamed \code{git_branch} object
}
\description{
Rename a branch
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Config user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Rename 'master' branch to 'dev'
branches(repo)
branch_rename(repository_head(repo), "dev")
branches(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/libgit2.R
\name{libgit2_version}
\alias{libgit2_version}
\title{Version of the libgit2 library}
\usage{
libgit2_version()
}
\value{
A list with major, minor and rev
}
\description{
Version of the libgit2 library that the bundled source code is
based on
}
\examples{
libgit2_version()
}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commit.R
\name{is_commit}
\alias{is_commit}
\title{Check if object is a git_commit object}
\usage{
is_commit(object)
}
\arguments{
\item{object}{Check if object is a git_commit object}
}
\value{
TRUE if object is a git_commit, else FALSE
}
\description{
Check if object is a git_commit object
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Commit a text file
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit_1 <- commit(repo, "First commit message")

## Check if commit
is_commit(commit_1)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commit.R
\name{commit}
\alias{commit}
\title{Commit}
\usage{
commit(
  repo = ".",
  message = NULL,
  all = FALSE,
  session = FALSE,
  author = NULL,
  committer = NULL
)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{message}{The commit message.}

\item{all}{Stage modified and deleted files. Files not added to
Git are not affected.}

\item{session}{Add sessionInfo to commit message. Default is
FALSE.}

\item{author}{Signature with author and author time of commit.}

\item{committer}{Signature with committer and commit time of
commit.}
}
\value{
A list of class \code{git_commit} with entries:
\describe{
  \item{sha}{
    The 40 character hexadecimal string of the SHA-1
  }
  \item{author}{
    An author signature
  }
  \item{committer}{
    The committer signature
  }
  \item{summary}{
    The short "summary" of a git commit message, comprising the first
    paragraph of the message with whitespace trimmed and squashed.
  }
  \item{message}{
    The message of a commit
  }
  \item{repo}{
    The \code{git_repository} object that contains the commit
  }
}
}
\description{
Commit
}
\examples{
\dontrun{
## Initialize a repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Config user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commit.R
\name{last_commit}
\alias{last_commit}
\title{Last commit}
\usage{
last_commit(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\description{
Get last commit in the current branch.
}
\examples{
\dontrun{
## Initialize a repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Config user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Get last commit
last_commit(repo)
last_commit(path)

## Coerce the last commit to a data.frame
as.data.frame(last_commit(path), "data.frame")

## Summary of last commit in repository
summary(last_commit(repo))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R, R/merge.R, R/repository.R, R/tree.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{diff}
\alias{merge}
\alias{head}
\alias{as.data.frame}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{base}{\code{\link[base]{as.data.frame}}, \code{\link[base]{diff}}, \code{\link[base]{merge}}}

  \item{utils}{\code{\link[utils]{head}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blob.R
\name{is_binary}
\alias{is_binary}
\title{Is blob binary}
\usage{
is_binary(blob = NULL)
}
\arguments{
\item{blob}{The blob \code{object}.}
}
\value{
TRUE if binary data, FALSE if not.
}
\description{
Is blob binary
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Commit a text file
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit_1 <- commit(repo, "First commit message")

## Check if binary
b_text <- tree(commit_1)["example.txt"]
is_binary(b_text)

## Commit plot file (binary)
x <- 1:100
y <- x^2
png(file.path(path, "plot.png"))
plot(y ~ x, type = "l")
dev.off()
add(repo, "plot.png")
commit_2 <- commit(repo, "Second commit message")

## Check if binary
b_png <- tree(commit_2)["plot.png"]
is_binary(b_png)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sha.R
\name{sha}
\alias{sha}
\alias{sha.git_blob}
\alias{sha.git_branch}
\alias{sha.git_commit}
\alias{sha.git_note}
\alias{sha.git_reference}
\alias{sha.git_reflog_entry}
\alias{sha.git_tag}
\alias{sha.git_tree}
\alias{sha.git_fetch_head}
\alias{sha.git_merge_result}
\title{Get the SHA-1 of a git object}
\usage{
sha(object)

\method{sha}{git_blob}(object)

\method{sha}{git_branch}(object)

\method{sha}{git_commit}(object)

\method{sha}{git_note}(object)

\method{sha}{git_reference}(object)

\method{sha}{git_reflog_entry}(object)

\method{sha}{git_tag}(object)

\method{sha}{git_tree}(object)

\method{sha}{git_fetch_head}(object)

\method{sha}{git_merge_result}(object)
}
\arguments{
\item{object}{a git object to get the SHA-1 from.}
}
\value{
The 40 character hexadecimal string of the SHA-1.
}
\description{
Get the 40 character hexadecimal string of the SHA-1.
}
\examples{
\dontrun{
## Create a directory in tempdir
path <- tempfile(pattern="git2r-")
dir.create(path)

## Initialize a repository
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "test.txt"))
add(repo, "test.txt")
commit(repo, "Commit message 1")

## Get the SHA-1 of the last commit
sha(last_commit(repo))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commit.R
\name{parents}
\alias{parents}
\title{Parents}
\usage{
parents(object = NULL)
}
\arguments{
\item{object}{a git_commit object.}
}
\value{
list of git_commit objects
}
\description{
Get parents of a commit.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
writeLines("First line.",
           file.path(path, "example.txt"))
add(repo, "example.txt")
commit_1 <- commit(repo, "First commit message")

## commit_1 has no parents
parents(commit_1)

## Update 'example.txt' and commit
writeLines(c("First line.", "Second line."),
           file.path(path, "example.txt"))
add(repo, "example.txt")
commit_2 <- commit(repo, "Second commit message")

## commit_2 has commit_1 as parent
parents(commit_2)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tag.R
\name{tag}
\alias{tag}
\title{Create tag targeting HEAD commit in repository}
\usage{
tag(
  object = ".",
  name = NULL,
  message = NULL,
  session = FALSE,
  tagger = NULL,
  force = FALSE
)
}
\arguments{
\item{object}{The repository \code{object}.}

\item{name}{Name for the tag.}

\item{message}{The tag message. Specify a tag message to create an
annotated tag. A lightweight tag is created if the message
parameter is \code{NULL}.}

\item{session}{Add sessionInfo to tag message. Default is FALSE.}

\item{tagger}{The tagger (author) of the tag}

\item{force}{Overwrite existing tag. Default = FALSE}
}
\value{
invisible(\code{git_tag}) object
}
\description{
Create tag targeting HEAD commit in repository
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Commit a text file
filename <- file.path(path, "example.txt")
writeLines("Hello world!", filename)
add(repo, "example.txt")
commit(repo, "First commit message")

## Create an annotated tag
tag(repo, "v1.0", "Tag message")

## List tags
tags(repo)

## Make a change to the text file and commit.
writeLines(c("Hello world!", "HELLO WORLD!"), filename)
add(repo, "example.txt")
commit(repo, "Second commit message")

## Create a lightweight tag
tag(repo, "v2.0")

## List tags
tags(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree.R
\name{ls_tree}
\alias{ls_tree}
\title{List the contents of a tree object}
\usage{
ls_tree(tree = NULL, repo = ".", recursive = TRUE)
}
\arguments{
\item{tree}{default (\code{NULL}) is the tree of the last commit
in \code{repo}. Can also be a \code{git_tree} object or a
character that identifies a tree in the repository (see
\sQuote{Examples}).}

\item{repo}{never used if \code{tree} is a \code{git_tree}
object. A \code{git_repository} object, or a path (default =
'.') to a repository.}

\item{recursive}{default is to recurse into sub-trees.}
}
\value{
A data.frame with the following columns: \describe{
    \item{mode}{UNIX file attribute of the tree entry}
    \item{type}{type of object} \item{sha}{sha of the object}
    \item{path}{path relative to the root tree}
    \item{name}{filename of the tree entry} \item{len}{object size
    of blob (file) entries. NA for other objects.}  }
}
\description{
Traverse the entries in a tree and its subtrees.  Akin to the 'git
ls-tree' command.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
dir.create(file.path(path, "subfolder"))
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create three files and commit
writeLines("First file",  file.path(path, "example-1.txt"))
writeLines("Second file", file.path(path, "subfolder/example-2.txt"))
writeLines("Third file",  file.path(path, "example-3.txt"))
add(repo, c("example-1.txt", "subfolder/example-2.txt", "example-3.txt"))
commit(repo, "Commit message")

## Traverse tree entries and its subtrees.
## Various approaches that give identical result.
ls_tree(tree = tree(last_commit(path)))
ls_tree(tree = tree(last_commit(repo)))
ls_tree(repo = path)
ls_tree(repo = repo)

## Skip content in subfolder
ls_tree(repo = repo, recursive = FALSE)

## Start in subfolder
ls_tree(tree = "HEAD:subfolder", repo = repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stash.R
\name{stash_apply}
\alias{stash_apply}
\title{Apply stash}
\usage{
stash_apply(object = ".", index = 1)
}
\arguments{
\item{object}{path to a repository, or a \code{git_repository}
object, or the stash \code{object} to pop. Default is a
\code{path = '.'} to a reposiory.}

\item{index}{The index to the stash to apply. Only used when
\code{object} is a path to a repository or a
\code{git_repository} object. Default is \code{index = 1}.}
}
\value{
invisible NULL
}
\description{
Apply a single stashed state from the stash list.
}
\details{
If local changes in the working directory conflict with changes in
the stash then an error will be raised. In this case, the index
will always remain unmodified and all files in the working
directory will remain unmodified. However, if you are restoring
untracked files or ignored files and there is a conflict when
applying the modified files, then those files will remain in the
working directory.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

# Configure a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

# Create a file, add and commit
writeLines("Hello world!", file.path(path, "test.txt"))
add(repo, 'test.txt')
commit(repo, "Commit message")

# Change file
writeLines(c("Hello world!", "HELLO WORLD!"), file.path(path, "test.txt"))

# Create stash in repository
stash(repo)

# Change file
writeLines(c("Hello world!", "HeLlO wOrLd!"), file.path(path, "test.txt"))

# Create stash in repository
stash(repo)

# View stashes
stash_list(repo)

# Read file
readLines(file.path(path, "test.txt"))

# Apply latest git_stash object in repository
stash_apply(stash_list(repo)[[1]])

# Read file
readLines(file.path(path, "test.txt"))

# View stashes
stash_list(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch.R
\name{fetch_heads}
\alias{fetch_heads}
\title{Get updated heads during the last fetch.}
\usage{
fetch_heads(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
list with \code{git_fetch_head} entries. NULL if there is
    no FETCH_HEAD file.
}
\description{
Get updated heads during the last fetch.
}
\examples{
\dontrun{
## Initialize three temporary repositories
path_bare <- tempfile(pattern="git2r-")
path_repo_1 <- tempfile(pattern="git2r-")
path_repo_2 <- tempfile(pattern="git2r-")

dir.create(path_bare)
dir.create(path_repo_1)
dir.create(path_repo_2)

bare_repo <- init(path_bare, bare = TRUE)
repo_1 <- clone(path_bare, path_repo_1)
repo_2 <- clone(path_bare, path_repo_2)

config(repo_1, user.name = "Alice", user.email = "alice@example.org")
config(repo_2, user.name = "Bob", user.email = "bob@example.org")

## Add changes to repo 1
writeLines("Lorem ipsum dolor sit amet",
           con = file.path(path_repo_1, "example.txt"))
add(repo_1, "example.txt")
commit(repo_1, "Commit message")

## Push changes from repo 1 to origin (bare_repo)
push(repo_1, "origin", "refs/heads/master")

## Fetch changes from origin (bare_repo) to repo 2
fetch(repo_2, "origin")

## List updated heads
fetch_heads(repo_2)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bundle_r_package.R
\name{bundle_r_package}
\alias{bundle_r_package}
\title{Bundle bare repo of package}
\usage{
bundle_r_package(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
Invisible bundled \code{git_repository} object
}
\description{
Clone the package git repository as a bare repository to
\code{pkg/inst/pkg.git}
}
\examples{
\dontrun{
## Initialize repository
path <- tempfile()
dir.create(path)
path <- file.path(path, "git2r")
repo <- clone("https://github.com/ropensci/git2r.git", path)

## Bundle bare repository in package
bundle_r_package(repo)

## Build and install bundled package
wd <- setwd(dirname(path))
system(sprintf("R CMD build \%s", path))
pkg <- list.files(".", pattern = "[.]tar[.]gz$")
system(sprintf("R CMD INSTALL \%s", pkg))
setwd(wd)

## Reload package
detach("package:git2r", unload = TRUE)
library(git2r)

## Summarize last five commits of bundled repo
repo <- repository(system.file("git2r.git", package = "git2r"))
invisible(lapply(commits(repo, n = 5), summary))

## Plot content of bundled repo
plot(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree.R
\name{as.data.frame.git_tree}
\alias{as.data.frame.git_tree}
\title{Coerce entries in a git_tree to a \code{data.frame}}
\usage{
\method{as.data.frame}{git_tree}(x, ...)
}
\arguments{
\item{x}{The tree \code{object}}

\item{...}{Additional arguments. Not used.}
}
\value{
\code{data.frame}
}
\description{
The entries in a tree are coerced to a \code{data.frame}
}
\details{
The \code{data.frame} have the following columns:
\describe{

  \item{filemode}{
    The UNIX file attributes of a tree entry
  }

  \item{type}{
    String representation of the tree entry type
  }

  \item{sha}{
    The sha of a tree entry
  }

  \item{name}{
    The filename of a tree entry
  }

}
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
dir.create(file.path(path, "subfolder"))
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create three files and commit
writeLines("First file",  file.path(path, "example-1.txt"))
writeLines("Second file", file.path(path, "subfolder/example-2.txt"))
writeLines("Third file",  file.path(path, "example-3.txt"))
add(repo, c("example-1.txt", "subfolder/example-2.txt", "example-3.txt"))
commit(repo, "Commit message")

## Display tree
tree(last_commit(repo))

## Coerce tree to a data.frame
df <- as.data.frame(tree(last_commit(repo)))
df
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/index.R
\name{rm_file}
\alias{rm_file}
\title{Remove files from the working tree and from the index}
\usage{
rm_file(repo = ".", path = NULL)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{path}{character vector with filenames to remove. Only files
known to Git are removed.}
}
\value{
invisible(NULL)
}
\description{
Remove files from the working tree and from the index
}
\examples{
\dontrun{
## Initialize a repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file
writeLines("Hello world!", file.path(path, "file-to-remove.txt"))

## Add file to repository
add(repo, "file-to-remove.txt")
commit(repo, "First commit message")

## Remove file
rm_file(repo, "file-to-remove.txt")

## View status of repository
status(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{init}
\alias{init}
\title{Init a repository}
\usage{
init(path = ".", bare = FALSE, branch = NULL)
}
\arguments{
\item{path}{A path to where to init a git repository}

\item{bare}{If TRUE, a Git repository without a working directory
is created at the pointed path. If FALSE, provided path will
be considered as the working directory into which the .git
directory will be created.}

\item{branch}{Use the specified name for the initial branch in the
newly created repository. If \code{branch=NULL}, fall back to
the default name.}
}
\value{
A \code{git_repository} object
}
\description{
Init a repository
}
\examples{
\dontrun{
## Initialize a repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)
is_bare(repo)

## Initialize a bare repository
path_bare <- tempfile(pattern="git2r-")
dir.create(path_bare)
repo_bare <- init(path_bare, bare = TRUE)
is_bare(repo_bare)
}
}
\seealso{
\link{repository}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stash.R
\name{stash_pop}
\alias{stash_pop}
\title{Pop stash}
\usage{
stash_pop(object = ".", index = 1)
}
\arguments{
\item{object}{path to a repository, or a \code{git_repository}
object, or the stash \code{object} to pop. Default is a
\code{path = '.'} to a reposiory.}

\item{index}{The index to the stash to pop. Only used when
\code{object} is a path to a repository or a
\code{git_repository} object. Default is \code{index = 1}.}
}
\value{
invisible NULL
}
\description{
Apply a single stashed state from the stash list and remove it
from the list if successful.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

# Configure a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

# Create a file, add and commit
writeLines("Hello world!", file.path(path, "test.txt"))
add(repo, 'test.txt')
commit(repo, "Commit message")

# Change file
writeLines(c("Hello world!", "HELLO WORLD!"), file.path(path, "test.txt"))

# Create stash in repository
stash(repo)

# Change file
writeLines(c("Hello world!", "HeLlO wOrLd!"), file.path(path, "test.txt"))

# Create stash in repository
stash(repo)

# View stashes
stash_list(repo)

# Read file
readLines(file.path(path, "test.txt"))

# Pop latest git_stash object in repository
stash_pop(stash_list(repo)[[1]])

# Read file
readLines(file.path(path, "test.txt"))

# View stashes
stash_list(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signature.R
\name{default_signature}
\alias{default_signature}
\title{Get the signature}
\usage{
default_signature(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
A \code{git_signature} object with entries:
}
\description{
Get the signature according to the repository's configuration
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Get the default signature
default_signature(repo)

## Change user
config(repo, user.name = "Bob", user.email = "bob@example.org")

## Get the default signature
default_signature(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stash.R
\name{stash}
\alias{stash}
\title{Stash}
\usage{
stash(
  repo = ".",
  message = as.character(Sys.time()),
  index = FALSE,
  untracked = FALSE,
  ignored = FALSE,
  stasher = NULL
)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{message}{Optional description. Defaults to current time.}

\item{index}{All changes already added to the index are left
intact in the working directory. Default is FALSE}

\item{untracked}{All untracked files are also stashed and then
cleaned up from the working directory. Default is FALSE}

\item{ignored}{All ignored files are also stashed and then cleaned
up from the working directory. Default is FALSE}

\item{stasher}{Signature with stasher and time of stash}
}
\value{
invisible \code{git_stash} object if anything to stash
    else NULL
}
\description{
Stash
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

# Configure a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

# Create a file, add and commit
writeLines("Hello world!", file.path(path, "test.txt"))
add(repo, 'test.txt')
commit(repo, "Commit message")

# Change file
writeLines(c("Hello world!", "HELLO WORLD!"), file.path(path, "test.txt"))

# Check status of repository
status(repo)

# Create stash in repository
stash(repo)

# Check status of repository
status(repo)

# View stash
stash_list(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tag.R
\name{is_tag}
\alias{is_tag}
\title{Check if object is a git_tag object}
\usage{
is_tag(object)
}
\arguments{
\item{object}{Check if object is a git_tag object}
}
\value{
TRUE if object is a git_tag, else FALSE
}
\description{
Check if object is a git_tag object
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Commit a text file
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Create tag
tag(repo, "Tagname", "Tag message")

is_tag(tags(repo)[[1]])
is_tag(last_commit(repo))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{is_detached}
\alias{is_detached}
\title{Check if HEAD of repository is detached}
\usage{
is_detached(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
\code{TRUE} if repository HEAD is detached, else
    \code{FALSE}.
}
\description{
Check if HEAD of repository is detached
}
\examples{
\dontrun{
## Create and initialize a repository in a temporary directory
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit_1 <- commit(repo, "Commit message 1")

## Change file, add and commit
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "Commit message 2")

## HEAD of repository is not detached
is_detached(repo)

## Checkout first commit
checkout(commit_1)

## HEAD of repository is detached
is_detached(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remote.R
\name{remote_rename}
\alias{remote_rename}
\title{Rename a remote}
\usage{
remote_rename(repo = ".", oldname = NULL, newname = NULL)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{oldname}{Old name of the remote}

\item{newname}{New name of the remote}
}
\value{
NULL, invisibly
}
\description{
Rename a remote
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name="Alice", user.email="alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Add a remote
remote_add(repo, "playground", "https://example.org/git2r/playground")
remotes(repo)
remote_url(repo, "playground")

## Rename a remote
remote_rename(repo, "playground", "foobar")
remotes(repo)
remote_url(repo, "foobar")

## Set remote url
remote_set_url(repo, "foobar", "https://example.org/git2r/foobar")
remotes(repo)
remote_url(repo, "foobar")

## Remove a remote
remote_remove(repo, "foobar")
remotes(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{workdir}
\alias{workdir}
\title{Workdir of repository}
\usage{
workdir(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
Character vector with the path of the workdir. If the
repository is bare, \code{NULL} will be returned.
}
\description{
Workdir of repository
}
\examples{
\dontrun{
## Create a directory in tempdir
path <- tempfile(pattern="git2r-")
dir.create(path)

## Initialize a repository
repo <- init(path)

## Get the path of the workdir for repository
workdir(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blob.R
\name{hash}
\alias{hash}
\title{Determine the sha from a blob string}
\usage{
hash(data = NULL)
}
\arguments{
\item{data}{The string vector to hash.}
}
\value{
A string vector with the sha for each string in data.
}
\description{
The blob is not written to the object database.
}
\examples{
\dontrun{
identical(hash(c("Hello, world!\n",
                 "test content\n")),
               c("af5626b4a114abcb82d63db7c8082c3c4756e51b",
                 "d670460b4b4aece5915caf5c68d12f560a9fe3e4"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blob.R
\name{blob_create}
\alias{blob_create}
\title{Create blob from file on disk}
\usage{
blob_create(repo = ".", path = NULL, relative = TRUE)
}
\arguments{
\item{repo}{The repository where the blob(s) will be written. Can
be a bare repository. A \code{git_repository} object, or a
path to a repository, or \code{NULL}.  If the \code{repo}
argument is \code{NULL}, the repository is searched for with
\code{\link{discover_repository}} in the current working
directory.}

\item{path}{The file(s) from which the blob will be created.}

\item{relative}{TRUE if the file(s) from which the blob will be
created is relative to the repository's working dir. Default
is TRUE.}
}
\value{
list of S3 class git_blob \code{objects}
}
\description{
Read a file from the filesystem and write its content to the
Object Database as a loose blob. The method is vectorized and
accepts a vector of files to create blobs from.
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create blobs from files relative to workdir
writeLines("Hello, world!", file.path(path, "example-1.txt"))
writeLines("test content", file.path(path, "example-2.txt"))
blob_list_1 <- blob_create(repo, c("example-1.txt",
                                   "example-2.txt"))

## Create blobs from files not relative to workdir
temp_file_1 <- tempfile()
temp_file_2 <- tempfile()
writeLines("Hello, world!", temp_file_1)
writeLines("test content", temp_file_2)
blob_list_2 <- blob_create(repo, c(temp_file_1, temp_file_2),
                           relative = FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{git_config_files}
\alias{git_config_files}
\title{Locate the path to configuration files}
\usage{
git_config_files(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
a \code{data.frame} with one row per potential
    configuration file where \code{NA} means not found.
}
\description{
Potential configuration files:
\describe{
  \item{system}{
    Locate the path to the system configuration file. If
    '/etc/gitconfig' doesn't exist, it will look for
    '\%PROGRAMFILES\%'.
  }
  \item{xdg}{
    Locate the path to the global xdg compatible configuration
    file. The xdg compatible configuration file is usually located
    in '$HOME/.config/git/config'. This method will try to guess
    the full path to that file, if the file exists.
  }
  \item{global}{
    The user or global configuration file is usually located in
    '$HOME/.gitconfig'. This method will try to guess the full
    path to that file, if the file exists.
  }
  \item{local}{
    Locate the path to the repository specific configuration file,
    if the file exists.
  }
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkout.R
\name{checkout}
\alias{checkout}
\title{Checkout}
\usage{
checkout(
  object = NULL,
  branch = NULL,
  create = FALSE,
  force = FALSE,
  path = NULL,
  ...
)
}
\arguments{
\item{object}{A path to a repository, or a \code{git_repository}
object, or a \code{git_commit} object, or a \code{git_tag}
object, or a \code{git_tree} object.}

\item{branch}{name of the branch to check out. Only used if object
is a path to a repository or a \code{git_repository} object.}

\item{create}{create branch if it doesn't exist. Only used if
object is a path to a repository or a \code{git_repository}
object.}

\item{force}{If \code{TRUE}, then make working directory match
target. This will throw away local changes. Default is
\code{FALSE}.}

\item{path}{Limit the checkout operation to only certain
paths. This argument is only used if branch is NULL. Default
is \code{NULL}.}

\item{...}{Additional arguments. Not used.}
}
\value{
invisible NULL
}
\description{
Update files in the index and working tree to match the content of
the tree pointed at by the treeish object (commit, tag or tree).
The default checkout strategy (\code{force = FALSE}) will only
make modifications that will not lose changes. Use \code{force =
TRUE} to force working directory to look like index.
}
\examples{
\dontrun{
## Create directories and initialize repositories
path_bare <- tempfile(pattern="git2r-")
path_repo_1 <- tempfile(pattern="git2r-")
path_repo_2 <- tempfile(pattern="git2r-")
dir.create(path_bare)
dir.create(path_repo_1)
dir.create(path_repo_2)
repo_bare <- init(path_bare, bare = TRUE)

## Clone to repo 1 and config user
repo_1 <- clone(path_bare, path_repo_1)
config(repo_1, user.name = "Alice", user.email = "alice@example.org")

## Add changes to repo 1 and push to bare
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path_repo_1, "test.txt"))
add(repo_1, "test.txt")
commit(repo_1, "First commit message")
push(repo_1, "origin", "refs/heads/master")

## Create and checkout 'dev' branch in repo 1
checkout(repo_1, "dev", create = TRUE)

## Add changes to 'dev' branch in repo 1 and push to bare
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path_repo_1, "test.txt"))
add(repo_1, "test.txt")
commit(repo_1, "Second commit message")
push(repo_1, "origin", "refs/heads/dev")

## Clone to repo 2
repo_2 <- clone(path_bare, path_repo_2)
config(repo_2, user.name = "Bob", user.email = "bob@example.org")

## Read content of 'test.txt'
readLines(file.path(path_repo_2, "test.txt"))

## Checkout dev branch
checkout(repo_2, "dev")

## Read content of 'test.txt'
readLines(file.path(path_repo_2, "test.txt"))

## Edit "test.txt" in repo_2
writeLines("Hello world!", con = file.path(path_repo_2, "test.txt"))

## Check status
status(repo_2)

## Checkout "test.txt"
checkout(repo_2, path = "test.txt")

## Check status
status(repo_2)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odb.R
\name{odb_objects}
\alias{odb_objects}
\title{List all objects available in the database}
\usage{
odb_objects(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
A data.frame with the following columns:
\describe{
  \item{sha}{The sha of the object}
  \item{type}{The type of the object}
  \item{len}{The length of the object}
}
}
\description{
List all objects available in the database
}
\examples{
\dontrun{
## Create a directory in tempdir
path <- tempfile(pattern="git2r-")
dir.create(path)

## Initialize a repository
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "test.txt"))
add(repo, "test.txt")
commit(repo, "Commit message 1")

## Create tag
tag(repo, "Tagname", "Tag message")

## List objects in repository
odb_objects(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/credential.R
\name{cred_env}
\alias{cred_env}
\title{Create a new environmental credential object}
\usage{
cred_env(username = NULL, password = NULL)
}
\arguments{
\item{username}{The name of the environmental variable that holds
the username for the authentication.}

\item{password}{The name of the environmental variable that holds
the password for the authentication.}
}
\value{
A list of class \code{cred_env} with entries:
\describe{
  \item{username}{
    The name of the environmental variable that holds
    the username for the authentication.
  }
  \item{password}{
    The name of the environmental variable that holds
    the password for the authentication.
  }
}
}
\description{
Environmental variables can be written to the file
\code{.Renviron}. This file is read by \emph{R} during startup,
see \code{\link[base]{Startup}}.
}
\examples{
\dontrun{
## Create an environmental credential object for the username and
## password.
cred <- cred_env("NAME_OF_ENV_VARIABLE_WITH_USERNAME",
                 "NAME_OF_ENV_VARIABLE_WITH_PASSWORD")
repo <- repository("git2r")
push(repo, credentials = cred)
}
}
\seealso{
Other git credential functions: 
\code{\link{cred_ssh_key}()},
\code{\link{cred_token}()},
\code{\link{cred_user_pass}()}
}
\concept{git credential functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch.R
\name{branch_delete}
\alias{branch_delete}
\title{Delete a branch}
\usage{
branch_delete(branch = NULL)
}
\arguments{
\item{branch}{The branch}
}
\value{
invisible NULL
}
\description{
Delete a branch
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit_1 <- commit(repo, "First commit message")

## Create a 'dev' branch
dev <- branch_create(commit_1, name = "dev")
branches(repo)

## Delete 'dev' branch
branch_delete(dev)
branches(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/index.R
\name{add}
\alias{add}
\title{Add file(s) to index}
\usage{
add(repo = ".", path = NULL, force = FALSE)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{path}{Character vector with file names or shell glob
patterns that will matched against files in the repository's
working directory. Each file that matches will be added to the
index (either updating an existing entry or adding a new
entry).}

\item{force}{Add ignored files. Default is FALSE.}
}
\value{
invisible(NULL)
}
\description{
Add file(s) to index
}
\examples{
\dontrun{
## Initialize a repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file
writeLines("a", file.path(path, "a.txt"))

## Add file to repository and view status
add(repo, "a.txt")
status(repo)

## Add file with a leading './' when the repository working
## directory is the current working directory
setwd(path)
writeLines("b", file.path(path, "b.txt"))
add(repo, "./b.txt")
status(repo)

## Add a file in a sub-folder with sub-folder as the working
## directory. Create a file in the root of the repository
## working directory that will remain untracked.
dir.create(file.path(path, "sub_dir"))
setwd("./sub_dir")
writeLines("c", file.path(path, "c.txt"))
writeLines("c", file.path(path, "sub_dir/c.txt"))
add(repo, "c.txt")
status(repo)

## Add files with glob expansion when the current working
## directory is outside the repository's working directory.
setwd(tempdir())
dir.create(file.path(path, "glob_dir"))
writeLines("d", file.path(path, "glob_dir/d.txt"))
writeLines("e", file.path(path, "glob_dir/e.txt"))
writeLines("f", file.path(path, "glob_dir/f.txt"))
writeLines("g", file.path(path, "glob_dir/g.md"))
add(repo, "glob_dir/*txt")
status(repo)

## Add file with glob expansion with a relative path when
## the current working directory is inside the repository's
## working directory.
setwd(path)
add(repo, "./glob_dir/*md")
status(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch.R
\name{branch_target}
\alias{branch_target}
\title{Get target (sha) pointed to by a branch}
\usage{
branch_target(branch = NULL)
}
\arguments{
\item{branch}{The branch}
}
\value{
sha or NA if not a direct reference
}
\description{
Get target (sha) pointed to by a branch
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Config user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Get target (sha) pointed to by 'master' branch
branch_target(repository_head(repo))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/credential.R
\name{cred_ssh_key}
\alias{cred_ssh_key}
\title{Create a new passphrase-protected ssh key credential object}
\usage{
cred_ssh_key(
  publickey = ssh_path("id_rsa.pub"),
  privatekey = ssh_path("id_rsa"),
  passphrase = character(0)
)
}
\arguments{
\item{publickey}{The path to the public key of the
credential. Default is \code{ssh_path("id_rsa.pub")}}

\item{privatekey}{The path to the private key of the
credential. Default is \code{ssh_path("id_rsa")}}

\item{passphrase}{The passphrase of the credential. Default is
\code{character(0)}. If getPass is installed and private key
is passphrase protected \code{getPass::getPass()} will be
called to allow for interactive and obfuscated interactive
input of the passphrase.}
}
\value{
A list of class \code{cred_ssh_key} with entries:
\describe{
  \item{publickey}{
    The path to the public key of the credential
  }
  \item{privatekey}{
    The path to the private key of the credential
  }
  \item{passphrase}{
    The passphrase of the credential
  }
}
}
\description{
Create a new passphrase-protected ssh key credential object
}
\examples{
\dontrun{
## Create a ssh key credential object. It can optionally be
## passphrase-protected
cred <- cred_ssh_key(ssh_path("id_rsa.pub"), ssh_path("id_rsa"))
repo <- repository("git2r")
push(repo, credentials = cred)
}
}
\seealso{
Other git credential functions: 
\code{\link{cred_env}()},
\code{\link{cred_token}()},
\code{\link{cred_user_pass}()}
}
\concept{git credential functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch.R
\name{branch_create}
\alias{branch_create}
\title{Create a branch}
\usage{
branch_create(commit = last_commit(), name = NULL, force = FALSE)
}
\arguments{
\item{commit}{Commit to which the branch should point. The default
is to use the \code{last_commit()} function to determine the
commit to which the branch should point.}

\item{name}{Name for the branch}

\item{force}{Overwrite existing branch. Default = FALSE}
}
\value{
invisible git_branch object
}
\description{
Create a branch
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
lines <- "Hello world!"
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit_1 <- commit(repo, "First commit message")

## Create a branch
branch_1 <- branch_create(commit_1, name = "test-branch")

## Add one more commit
lines <- c("Hello world!", "HELLO WORLD!")
writeLines(lines, file.path(path, "example.txt"))
add(repo, "example.txt")
commit_2 <- commit(repo, "Another commit message")

## Create a branch with the same name should fail
try(branch_create(commit_2, name = "test-branch"), TRUE)

## Force it
branch_2 <- branch_create(commit_2, name = "test-branch", force = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contributions.R
\name{contributions}
\alias{contributions}
\title{Contributions}
\usage{
contributions(
  repo = ".",
  breaks = c("month", "year", "quarter", "week", "day"),
  by = c("commits", "author")
)
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{breaks}{Default is \code{month}. Change to year, quarter,
week or day as necessary.}

\item{by}{Contributions by "commits" or "author". Default is "commits".}
}
\value{
A \code{data.frame} with contributions.
}
\description{
See contributions to a Git repo
}
\examples{
\dontrun{
## Create directories and initialize repositories
path_bare <- tempfile(pattern="git2r-")
path_repo_1 <- tempfile(pattern="git2r-")
path_repo_2 <- tempfile(pattern="git2r-")
dir.create(path_bare)
dir.create(path_repo_1)
dir.create(path_repo_2)
repo_bare <- init(path_bare, bare = TRUE)

## Clone to repo 1 and config user
repo_1 <- clone(path_bare, path_repo_1)
config(repo_1, user.name = "Alice", user.email = "alice@example.org")

## Add changes to repo 1 and push to bare
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path_repo_1, "test.txt"))
add(repo_1, "test.txt")
commit(repo_1, "First commit message")

## Add more changes to repo 1
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path_repo_1, "test.txt"))
add(repo_1, "test.txt")
commit(repo_1, "Second commit message")

## Push to bare
push(repo_1, "origin", "refs/heads/master")

## Clone to repo 2
repo_2 <- clone(path_bare, path_repo_2)
config(repo_2, user.name = "Bob", user.email = "bob@example.org")

## Add changes to repo 2
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad",
  "minim veniam, quis nostrud exercitation ullamco laboris nisi ut")
writeLines(lines, file.path(path_repo_2, "test.txt"))
add(repo_2, "test.txt")
commit(repo_2, "Third commit message")

## Push to bare
push(repo_2, "origin", "refs/heads/master")

## Pull changes to repo 1
pull(repo_1)

## View contributions by day
contributions(repo_1)

## View contributions by author and day
contributions(repo_1, by = "author")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch.R
\name{branches}
\alias{branches}
\title{Branches}
\usage{
branches(repo = ".", flags = c("all", "local", "remote"))
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}

\item{flags}{Filtering flags for the branch listing. Valid values
are 'all', 'local' or 'remote'}
}
\value{
list of branches in repository
}
\description{
List branches in repository
}
\examples{
\dontrun{
## Initialize repositories
path_bare <- tempfile(pattern="git2r-")
path_repo <- tempfile(pattern="git2r-")
dir.create(path_bare)
dir.create(path_repo)
repo_bare <- init(path_bare, bare = TRUE)
repo <- clone(path_bare, path_repo)

## Config first user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path_repo, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Push commits from repository to bare repository
## Adds an upstream tracking branch to branch 'master'
push(repo, "origin", "refs/heads/master")

## List branches
branches(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/push.R
\name{push}
\alias{push}
\title{Push}
\usage{
push(
  object = ".",
  name = NULL,
  refspec = NULL,
  force = FALSE,
  credentials = NULL,
  set_upstream = FALSE
)
}
\arguments{
\item{object}{path to repository, or a \code{git_repository} or
\code{git_branch}.}

\item{name}{The remote's name. Default is NULL.}

\item{refspec}{The refspec to be pushed. Default is NULL.}

\item{force}{Force your local revision to the remote repo. Use it
with care. Default is FALSE.}

\item{credentials}{The credentials for remote repository
access. Default is NULL. To use and query an ssh-agent for the
ssh key credentials, let this parameter be NULL (the default).}

\item{set_upstream}{Set the current local branch to track the
remote branch. Default is FALSE.}
}
\value{
invisible(NULL)
}
\description{
Push
}
\examples{
\dontrun{
## Initialize two temporary repositories
path_bare <- tempfile(pattern="git2r-")
path_repo <- tempfile(pattern="git2r-")
dir.create(path_bare)
dir.create(path_repo)
repo_bare <- init(path_bare, bare = TRUE)

## Clone the bare repository. This creates remote-tracking
## branches for each branch in the cloned repository.
repo <- clone(path_bare, path_repo)

## Config user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path_repo, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Push commits from repository to bare repository
push(repo, "origin", "refs/heads/master")

## Now, unset the remote-tracking branch to NULL to demonstrate
## the 'set_upstream' argument. Then push with 'set_upstream = TRUE'
## to add the upstream tracking branch to branch 'master' again.
branch_get_upstream(repository_head(repo))
branch_set_upstream(repository_head(repo), NULL)
branch_get_upstream(repository_head(repo))
push(repo, "origin", "refs/heads/master", set_upstream = TRUE)
branch_get_upstream(repository_head(repo))

## Change file and commit
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path_repo, "example.txt"))
add(repo, "example.txt")
commit(repo, "Second commit message")

## Push commits from repository to bare repository
push(repo)

## List commits in repository and bare repository
commits(repo)
commits(repo_bare)
}
}
\seealso{
\code{\link{cred_user_pass}}, \code{\link{cred_ssh_key}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/note.R
\name{note_default_ref}
\alias{note_default_ref}
\title{Default notes reference}
\usage{
note_default_ref(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
Character vector of length one with name of default notes
    reference
}
\description{
Get the default notes reference for a repository
}
\examples{
\dontrun{
## Create and initialize a repository in a temporary directory
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## View default notes reference
note_default_ref(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remote.R
\name{remote_ls}
\alias{remote_ls}
\title{List references in a remote repository}
\usage{
remote_ls(name = NULL, repo = NULL, credentials = NULL)
}
\arguments{
\item{name}{Character vector with the "remote" repository URL to
query or the name of the remote if a \code{repo} argument is
given.}

\item{repo}{an optional repository object used if remotes are
specified by name.}

\item{credentials}{The credentials for remote repository
access. Default is NULL. To use and query an ssh-agent for the
ssh key credentials, let this parameter be NULL (the default).}
}
\value{
Character vector for each reference with the associated
    commit IDs.
}
\description{
Displays references available in a remote repository along with
the associated commit IDs.  Akin to the 'git ls-remote' command.
}
\examples{
\dontrun{
remote_ls("https://github.com/ropensci/git2r")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reflog.R
\name{print.git_reflog_entry}
\alias{print.git_reflog_entry}
\title{Print a reflog entry}
\usage{
\method{print}{git_reflog_entry}(x, ...)
}
\arguments{
\item{x}{The reflog entry}

\item{...}{Unused}
}
\value{
None (invisible 'NULL').
}
\description{
Print a reflog entry
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## View repository HEAD reflog
reflog(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{is_empty}
\alias{is_empty}
\title{Check if repository is empty}
\usage{
is_empty(repo = ".")
}
\arguments{
\item{repo}{a path to a repository or a \code{git_repository}
object. Default is '.'}
}
\value{
\code{TRUE} if repository is empty else \code{FALSE}.
}
\description{
Check if repository is empty
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Check if it's an empty repository
is_empty(repo)

## Commit a file
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Check if it's an empty repository
is_empty(repo)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{as.data.frame.git_repository}
\alias{as.data.frame.git_repository}
\title{Coerce Git repository to a \code{data.frame}}
\usage{
\method{as.data.frame}{git_repository}(x, ...)
}
\arguments{
\item{x}{The repository \code{object}}

\item{...}{Additional arguments. Not used.}
}
\value{
\code{data.frame}
}
\description{
The commits in the repository are coerced to a \code{data.frame}
}
\details{
The \code{data.frame} have the following columns:
\describe{

  \item{sha}{
    The 40 character hexadecimal string of the SHA-1
  }

  \item{summary}{
    the short "summary" of the git commit message.
  }

  \item{message}{
    the full message of a commit
  }

  \item{author}{
    full name of the author
  }

  \item{email}{
    email of the author
  }

  \item{when}{
    time when the commit happened
  }

}
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create three files and commit
writeLines("First file",  file.path(path, "example-1.txt"))
writeLines("Second file", file.path(path, "example-2.txt"))
writeLines("Third file",  file.path(path, "example-3.txt"))
add(repo, "example-1.txt")
commit(repo, "Commit first file")
add(repo, "example-2.txt")
commit(repo, "Commit second file")
add(repo, "example-3.txt")
commit(repo, "Commit third file")

## Coerce commits to a data.frame
df <- as.data.frame(repo)
df
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commit.R
\name{ahead_behind}
\alias{ahead_behind}
\title{Ahead Behind}
\usage{
ahead_behind(local = NULL, upstream = NULL)
}
\arguments{
\item{local}{a git_commit object. Can also be a tag or a branch,
and in that case the commit will be the target of the tag or
branch.}

\item{upstream}{a git_commit object. Can also be a tag or a
branch, and in that case the commit will be the target of the
tag or branch.}
}
\value{
An integer vector of length 2 with number of commits that
    the upstream commit is ahead and behind the local commit
}
\description{
Count the number of unique commits between two commit objects.
}
\examples{
\dontrun{
## Create a directory in tempdir
path <- tempfile(pattern="git2r-")
dir.create(path)

## Initialize a repository
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "test.txt"))
add(repo, "test.txt")
commit_1 <- commit(repo, "Commit message 1")
tag_1 <- tag(repo, "Tagname1", "Tag message 1")

# Change file and commit
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path, "test.txt"))
add(repo, "test.txt")
commit_2 <- commit(repo, "Commit message 2")
tag_2 <- tag(repo, "Tagname2", "Tag message 2")

ahead_behind(commit_1, commit_2)
ahead_behind(tag_1, tag_2)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree.R
\name{summary.git_tree}
\alias{summary.git_tree}
\title{Summary of tree}
\usage{
\method{summary}{git_tree}(object, ...)
}
\arguments{
\item{object}{The tree \code{object}}

\item{...}{Additional arguments affecting the summary produced.}
}
\value{
None (invisible 'NULL').
}
\description{
Summary of tree
}
\examples{
\dontrun{
## Initialize a temporary repository
path <- tempfile(pattern="git2r-")
dir.create(path)
repo <- init(path)

## Create a user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")
writeLines("Hello world!", file.path(path, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

summary(tree(last_commit(repo)))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commit.R
\name{descendant_of}
\alias{descendant_of}
\title{Descendant}
\usage{
descendant_of(commit = NULL, ancestor = NULL)
}
\arguments{
\item{commit}{a git_commit object. Can also be a tag or a branch,
and in that case the commit will be the target of the tag or
branch.}

\item{ancestor}{a git_commit object to check if ancestor to
\code{commit}. Can also be a tag or a branch, and in that case
the commit will be the target of the tag or branch.}
}
\value{
TRUE if \code{commit} is descendant of \code{ancestor},
    else FALSE
}
\description{
Determine if a commit is the descendant of another commit
}
\examples{
\dontrun{
## Create a directory in tempdir
path <- tempfile(pattern="git2r-")
dir.create(path)

## Initialize a repository
repo <- init(path)
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Create a file, add and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path, "test.txt"))
add(repo, "test.txt")
commit_1 <- commit(repo, "Commit message 1")
tag_1 <- tag(repo, "Tagname1", "Tag message 1")

# Change file and commit
lines <- c(
  "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do",
  "eiusmod tempor incididunt ut labore et dolore magna aliqua.")
writeLines(lines, file.path(path, "test.txt"))
add(repo, "test.txt")
commit_2 <- commit(repo, "Commit message 2")
tag_2 <- tag(repo, "Tagname2", "Tag message 2")

descendant_of(commit_1, commit_2)
descendant_of(commit_2, commit_1)
descendant_of(tag_1, tag_2)
descendant_of(tag_2, tag_1)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/git2r.R
\docType{package}
\name{git2r}
\alias{git2r}
\title{git2r: R bindings to the libgit2 library}
\description{
git2r: R bindings to the libgit2 library.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch.R
\name{branch_set_upstream}
\alias{branch_set_upstream}
\title{Set remote tracking branch}
\usage{
branch_set_upstream(branch = NULL, name)
}
\arguments{
\item{branch}{The branch to configure}

\item{name}{remote-tracking or local branch to set as
upstream. Pass NULL to unset.}
}
\value{
invisible NULL
}
\description{
Set the upstream configuration for a given local branch
}
\examples{
\dontrun{
## Initialize two temporary repositories
path_bare <- tempfile(pattern="git2r-")
path_repo <- tempfile(pattern="git2r-")
dir.create(path_bare)
dir.create(path_repo)
repo_bare <- init(path_bare, bare = TRUE)
repo <- clone(path_bare, path_repo)

## Config user and commit a file
config(repo, user.name = "Alice", user.email = "alice@example.org")

## Write to a file and commit
lines <- "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do"
writeLines(lines, file.path(path_repo, "example.txt"))
add(repo, "example.txt")
commit(repo, "First commit message")

## Push commits from repository to bare repository
## Adds an upstream tracking branch to branch 'master'
push(repo, "origin", "refs/heads/master")

## Unset remote remote tracking branch
branch_get_upstream(repository_head(repo))
branch_set_upstream(repository_head(repo), NULL)
branch_get_upstream(repository_head(repo))

## Set remote tracking branch
branch_set_upstream(repository_head(repo), "origin/master")
branch_get_upstream(repository_head(repo))
}
}
