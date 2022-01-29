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
