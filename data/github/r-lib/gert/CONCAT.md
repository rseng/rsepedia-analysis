# gert <img src="man/figures/logo.png" align="right" alt="logo" width="120" height = "139" style = "border: none; float: right;">

*This package is a joint effort from [rOpenSci](https://ropensci.org/) and the [Tidyverse](https://www.tidyverse.org/) team.*

> Simple Git Client for R

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![runiverse-name](https://ropensci.r-universe.dev/badges/:name)
![runiverse-package](https://ropensci.r-universe.dev/badges/gert)
![cran-badge](http://www.r-pkg.org/badges/version/gert)
<!-- badges: end -->

Simple git client for R based on 'libgit2' with support for SSH and 
HTTPS remotes. All functions in gert use basic R data types (such as vectors
and data-frames) for their arguments and return values. User credentials are
shared with command line 'git' through the git-credential store and ssh keys
stored on disk or ssh-agent. On Linux, a somewhat recent version of 'libgit2'
is required; we provide a PPA for older Ubuntu LTS versions.

## Installation

```r
# To install the latest version
install.packages("gert", repos = c(
    ropensci = 'https://ropensci.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))
    
# To install the CRAN release:
install.packages("gert")
```

On Linux you need to install libgit2:

 - Debian: [libgit2-dev](https://packages.debian.org/buster/libgit2-dev)
 - Fedora / CentOS: [libgit2-devel](https://src.fedoraproject.org/rpms/libgit2)
 - Arch Linux [libgit2](https://archlinux.org/packages/extra/x86_64/libgit2)
 
If no suitable version of libgit2 is found, the package automatically tries to download a static build.


## Documentation:

 - homepage: https://docs.ropensci.org/gert
 - slides: https://jeroen.github.io/gert2019/#1

## Hello world

Some basic commands to get started with gert:

``` r
library(gert)
repo <- git_clone("https://github.com/r-lib/gert")
setwd("gert")

# Show some info
git_log(max = 10)

# Create a branch
git_branch_create("mybranch", checkout = TRUE)

# Commit things
writeLines("Lorem ipsum dolor sit amet", 'test.txt')
git_add('test.txt')
git_commit("Adding a file", author = "jerry <jerry@gmail.com>")
git_log(max = 10)

# Merge it in master
git_branch_checkout("master")
git_merge("mybranch")
git_branch_delete("mybranch")

# Remove the commit
git_reset_hard("HEAD^")
```

## Should I use HTTPS or SSH remotes?

On most platforms, gert supports both HTTPS or SSH remotes. If you don't have any preference, the safest choice is  __HTTPS remotes using a PAT as the password__. This is what I use myself as well. HTTPS remotes have the following benefits:

  - Your credentials are safely stored by your OS, accessible both to gert and command line `git`.
  - Https works on any network. However the ssh protocol requires port 22, which is often blocked on public wifi networks.
  - You can authenticate over https using the same GITHUB_PAT that you use for the GitHub API.
  - libgit2 supports https on all platforms (SSH support depends on libssh2 availability).
  
Again: no need to use your Github master password in gert/git. Instead [generate a personal access token](https://github.com/settings/tokens/new) and enter this as the password when pushing/pulling from https remotes. This works both with gert and with the git command line, even when you have 2FA enabled (which you should).

Ninja tip: use `credentials::set_github_pat()` to automatically set the `GITHUB_PAT` environment variable in your R session using the value stored in your git credential store. This is a safer way to store your PAT than hardcoding it in your `.Renviron`.

## Differences with `git2r`

Gert is based on [libgit2](https://libgit2.org/), just like the rOpenSci package [git2r](https://docs.ropensci.org/git2r/). Both are good packages. The well established git2r has been on CRAN since 2015, is actively maintained by Stefan Widgren, and is widely used. Gert was started in 2019, and takes a fresh approach based on more recent APIs in libgit2 and lessons learned from using git2r. Some of the main differences:

### Simplicity 

Gert is focused on high-level functions that shield the end-user from the complexity of libgit2. Functions in gert use standard R data types (such as vectors and data-frames) for their arguments and return values, which should be easy to work with for R users/packages. The target repository is either inferred from current working directory or is specified as a filepath. Branches and remotes are referred to by name, much like command line `git`. None of the functions in gert expose any externalptr types to the user.

```
> gert::git_log(max=6)
# A tibble: 6 x 6
  commit                        author                    time                files merge message             
* <chr>                         <chr>                     <dttm>              <int> <lgl> <chr>               
1 6f39ba6dae890d679970c0f8bf03… Jeroen Ooms <jeroenooms@… 2020-06-16 01:16:33    17 FALSE "Add some family ta…
2 c023c407a0f0bfa3955576bc3551… Jeroen Ooms <jeroenooms@… 2020-06-16 01:06:38     1 FALSE "Check for matching…
3 24234060ea8e54c73ddd0bce90ff… Jeroen Ooms <jeroenooms@… 2020-06-15 13:17:57     1 FALSE "Update fedora link…
4 e60b0fbad129f470a2f7065063fa… Jeroen Ooms <jeroenooms@… 2020-06-15 13:05:45     4 FALSE "Tweak docs and rea…
5 629420ddccbab51c1e78f472bf06… Jeroen Ooms <jeroenooms@… 2020-06-15 12:14:25     1 FALSE "More tests\n"      
6 a62ce14eb887e183ad0a3cf0e22c… Jeroen Ooms <jeroenooms@… 2020-06-15 12:06:41     1 FALSE "Fix unit test\n"   
```

For R users who are familiar with the `git` command line, gert should be mostly self-explanatory, and generally "just work".

### Automatic authentication

The overall goal for auth is that gert should successfully discover credentials whenever that would also be true for command line `git`. And, should that fail, there is a way to debug it.

To authenticate with a remote in git2r, you often need to manually pass your credentials in every call to, e.g., `git2r::clone()`. This is always the case for an https remote and is often the case even for an ssh remote. This creates special challenges for those new to `git` or for indirect use of git2r.

In gert, authentication is done automatically using the [credentials](https://docs.ropensci.org/credentials/articles/intro.html) package. This package calls out to the local OS credential store which is also used by the `git` command line. Therefore gert will automatically pick up on https credentials that are safely stored in your OS keychain. 

If no credentials are available from the store, gert will try to authenticate using your `GITHUB_PAT` (if set) for GitHub https remotes. If none of that works, it safely prompts the user for credentials using [askpass](https://github.com/jeroen/askpass#readme). Together, these methods should make https authentication "just work" in any scenario, without having to manually provide passwords in R.

Authentication with ssh remotes is a bit more complicated, but gert will again try to make this as smooth as possible. First of all, gert will tell you if SSH is supported when attaching the package (this will be the case on all modern systems):

```r
> library(gert)
Linking to libgit2 v1.0.0, ssh support: YES
Global config: /Users/jeroen/.gitconfig
Default user: Jeroen Ooms <jeroenooms@gmail.com
```

On Mac/Linux, gert first tries to authenticate using credentials from your `ssh-agent`. If that doesn't work it will look for a suitable ssh key on your system (usually `id_rsa`), and if it is protected with a passphrase, gert will safely prompt the user for the passphrase using [askpass](https://github.com/jeroen/askpass#readme).
If the user does not have an SSH key yet, the [credentials](https://docs.ropensci.org/credentials/articles/intro.html) package makes it easy to set that up.

```r
> library(credentials)
Found git version 2.24.3 (Apple Git-128)
Supported HTTPS credential helpers: cache, store
Found OpenSSH_8.1p1, LibreSSL 2.7.3
Default SSH key: /Users/jeroen/.ssh/id_rsa
```

One limitation that remains is that libgit2 does not support `ssh-agent` on Windows. This is [unlikely to change](https://github.com/libgit2/libgit2/issues/4958) because ssh-agent uses unix-sockets which do not exist in native Windows software.

### The libgit2 dependency

If you use Windows or macOS and you install gert from CRAN, it comes with "batteries included". Gert brings in prebuilt versions of external dependencies, like libgit2 and the 3rd party libraries needed to support SSH and TLS (for HTTPS). This approach guarantees that gert uses libraries that are properly configured for your operating system.

The git2r package takes another approach by bundling the libgit2 source code in the R package, and automatically building libgit2 on-the-fly when the R package is compiled. This is mostly for historical reasons, because until recently, libgit2 was not available on every Linux system.

However the problem is that configuring and building libgit2 is complicated (like most system libraries) and requires several platform-specific flags and system dependencies. As a result, git2r is sometimes installed with missing functionality, depending on what was detected during compilation. On macOS for example, some git2r users have SSH support but others do not. Weird problems due to missing libgit2 features turn out to be very persistent, and have caused a lot of frustration. For this reason, gert does not bundle and compile the libgit2 source, but instead always links to system libraries.

As usual, those who install gert as a source package, by choice on Windows and macOS or by necessity on Linux, **do** need to ensure the necessary system libraries are present, e.g.:

  * [libgit2-dev](https://packages.ubuntu.com/focal/libgit2-dev) on Debian/Ubuntu
  * [libgit2-devel](https://src.fedoraproject.org/rpms/libgit2) on Fedora
  * [libgit2](https://archlinux.org/packages/extra/x86_64/libgit2) on Arch Linux
  * [Homebrew libgit2](https://github.com/Homebrew/homebrew-core/blob/master/Formula/libgit2.rb) on macOS
  * [rtools40 libgit2](https://github.com/r-windows/rtools-packages/blob/master/mingw-w64-libgit2/PKGBUILD) on Windows

One disadvantage of this approach is that on very old versions of Ubuntu, the system-provided version of libgit2 is out of date, and we need to enable a PPA with more recent libgit2 backports. This is the case for Ubuntu Xenial (16.04) which is a system from 2016 that will be EOL in April 2021.

```sh
# Needed on Ubuntu 16.04
sudo add-apt-repository ppa:cran/libgit2
sudo apt-get install libgit2-dev
```

CI users do not need to worry about this, because we automatically enable this PPA on Travis and GitHub Actions. Outside of CI systems, very few people are running Ubuntu 16 anymore, most production servers have updated to Ubuntu 18 or 20 by now, so this is rarely an issue in practice.

---
title: "gert"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{gert}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  error = TRUE
)
```

gert is a simple git client based on 'libgit2' ([libgit2.org](https://libgit2.org)):

> libgit2 is a portable, pure C implementation of the Git core methods provided as a re-entrant linkable library with a solid API, allowing you to write native speed custom Git applications in any language which supports C bindings.

What this means for R users is that we can work with local and remote Git repositories from the comfort of R!

User-friendly authentication is a high priority for gert, which supports both SSH and HTTPS remotes on all platforms. User credentials are shared with command line Git through the git-credential store and ssh keys stored on disk or ssh-agent.

Let's attach gert.

```{r setup}
library(gert)
```

## Introduce yourself to Git

Before you can do anything with Git, you must first configure your user name and email. When you attach gert, it actually reveals whether you've already done this and, above, you can see that we have.

But what if you have not already configured your user name and email? Do this with `git_config_global_set()`:

```{r eval = FALSE}
git_config_global_set("user.name", "Jerry Johnson")
git_config_global_set("user.email", "jerry@gmail.com")
```

We can verify our success (and see all global options) with `git_config_global()`.

```{r}
git_config_global()
```

This is equivalent to these commands in command line Git:

```
git config --global user.name 'Jerry Johnson'
git config --global user.email 'jerry@gmail.com'
git config --global --list
```

To inspect and change local Git config, i.e. options specific to one repository, use `git_config()` and `git_config_set()`.

## Local repository basics

`gert::git_init()` is essentially `git init`; it's how we create a new local repository. You provide the path to the repository you want to create.

```{r}
(path <- file.path(tempdir(), "aaa", "bbb", "repo_ccc"))
dir.exists(path)

(r <- git_init(path))
dir.exists(path)
```

Note that all non-existing parts of the path are created: `aaa`, `bbb`, and `repo_ccc` (the actual git repository).

`git_find()` finds a git repository at or above the path you provide and errors otherwise.

```{r}
git_find(r)

dir.create(file.path(r, "child_dir"))
git_find(file.path(r, "child_dir"))

git_find(file.path(tempdir(), "aaa", "bbb"))
```

`git_init()` can also create a repository in a pre-existing directory, as long as it is empty.

```{r}
r2 <- file.path(tempdir(), "repo_ddd")
dir.create(r2)

git_init(r2)
```

Cleanup.

```{r}
unlink(r, recursive = TRUE)
unlink(r2, recursive = TRUE)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commit.R
\name{git_commit}
\alias{git_commit}
\alias{git_commit_all}
\alias{git_commit_info}
\alias{git_commit_id}
\alias{git_commit_descendant_of}
\alias{git_add}
\alias{git_rm}
\alias{git_status}
\alias{git_conflicts}
\alias{git_ls}
\alias{git_log}
\alias{git_stat_files}
\title{Stage and commit changes}
\usage{
git_commit(message, author = NULL, committer = NULL, repo = ".")

git_commit_all(message, author = NULL, committer = NULL, repo = ".")

git_commit_info(ref = "HEAD", repo = ".")

git_commit_id(ref = "HEAD", repo = ".")

git_commit_descendant_of(ancestor, ref = "HEAD", repo = ".")

git_add(files, force = FALSE, repo = ".")

git_rm(files, repo = ".")

git_status(staged = NULL, repo = ".")

git_conflicts(repo = ".")

git_ls(repo = ".")

git_log(ref = "HEAD", max = 100, repo = ".")

git_stat_files(files, ref = "HEAD", repo = ".")
}
\arguments{
\item{message}{a commit message}

\item{author}{A \link{git_signature} value, default is \code{\link[=git_signature_default]{git_signature_default()}}.}

\item{committer}{A \link{git_signature} value, default is same as \code{author}}

\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}

\item{ref}{revision string with a branch/tag/commit value}

\item{ancestor}{a reference to a potential ancestor commit}

\item{files}{vector of paths relative to the git root directory.
Use \code{"."} to stage all changed files.}

\item{force}{add files even if in gitignore}

\item{staged}{return only staged (TRUE) or unstaged files (FALSE).
Use \code{NULL} or \code{NA} to show both (default).}

\item{max}{lookup at most latest n parent commits}
}
\value{
\itemize{
\item \code{git_status()}, \code{git_ls()}: A data frame with one row per file
\item \code{git_log()}: A data frame with one row per commit
\item \code{git_commit()}, \code{git_commit_all()}: A SHA
}
}
\description{
To commit changes, start by \emph{staging} the files to be included in the commit
using \code{git_add()} or \code{git_rm()}. Use \code{git_status()} to see an overview of
staged and unstaged changes, and finally \code{git_commit()} creates a new commit
with currently staged files.

\code{git_commit_all()} is a convenience function that automatically stages and
commits all modified files. Note that \code{git_commit_all()} does \strong{not} add
new, untracked files to the repository. You need to make an explicit call to
\code{git_add()} to start tracking new files.

\code{git_log()} shows the most recent commits and \code{git_ls()} lists all the files
that are being tracked in the repository. \code{git_stat_files()}
}
\examples{
oldwd <- getwd()
repo <- file.path(tempdir(), "myrepo")
git_init(repo)
setwd(repo)

# Set a user if no default
if(!user_is_configured()){
  git_config_set("user.name", "Jerry")
  git_config_set("user.email", "jerry@gmail.com")
}

writeLines(letters[1:6], "alphabet.txt")
git_status()

git_add("alphabet.txt")
git_status()

git_commit("Start alphabet file")
git_status()

git_ls()

git_log()

cat(letters[7:9], file = "alphabet.txt", sep = "\n", append = TRUE)
git_status()

git_commit_all("Add more letters")

# cleanup
setwd(oldwd)
unlink(repo, recursive = TRUE)
}
\seealso{
Other git: 
\code{\link{git_archive}},
\code{\link{git_branch}()},
\code{\link{git_config}()},
\code{\link{git_diff}()},
\code{\link{git_fetch}()},
\code{\link{git_merge}()},
\code{\link{git_rebase}()},
\code{\link{git_remote}},
\code{\link{git_repo}},
\code{\link{git_signature}()},
\code{\link{git_stash}},
\code{\link{git_tag}}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stash.R
\name{git_stash}
\alias{git_stash}
\alias{git_stash_save}
\alias{git_stash_pop}
\alias{git_stash_drop}
\alias{git_stash_list}
\title{Stashing changes}
\usage{
git_stash_save(
  message = "",
  keep_index = FALSE,
  include_untracked = FALSE,
  include_ignored = FALSE,
  repo = "."
)

git_stash_pop(index = 0, repo = ".")

git_stash_drop(index = 0, repo = ".")

git_stash_list(repo = ".")
}
\arguments{
\item{message}{optional message to store the stash}

\item{keep_index}{changes already added to the index are left intact in
the working directory}

\item{include_untracked}{untracked files are also stashed and then
cleaned up from the working directory}

\item{include_ignored}{ignored files are also stashed and then cleaned
up from the working directory}

\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}

\item{index}{The position within the stash list. 0 points to the
most recent stashed state.}
}
\description{
Temporary stash away changed from the working directory.
}
\seealso{
Other git: 
\code{\link{git_archive}},
\code{\link{git_branch}()},
\code{\link{git_commit}()},
\code{\link{git_config}()},
\code{\link{git_diff}()},
\code{\link{git_fetch}()},
\code{\link{git_merge}()},
\code{\link{git_rebase}()},
\code{\link{git_remote}},
\code{\link{git_repo}},
\code{\link{git_signature}()},
\code{\link{git_tag}}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/open.R
\name{git_open}
\alias{git_open}
\title{Open local repository}
\usage{
git_open(repo = ".")
}
\arguments{
\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}
}
\value{
an pointer to the libgit2 repository
}
\description{
Returns a pointer to a libgit2 repository object.This function is mainly
for internal use; users should simply reference a repository in gert by
by the path to the directory.
}
\examples{
r <- tempfile(pattern = "gert")
git_init(r)
r_ptr <- git_open(r)
r_ptr
git_open(r_ptr)
git_info(r)

# cleanup
unlink(r, recursive = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repository.R
\name{git_repo}
\alias{git_repo}
\alias{git_init}
\alias{git_find}
\alias{git_info}
\title{Create or discover a local Git repository}
\usage{
git_init(path = ".", bare = FALSE)

git_find(path = ".")

git_info(repo = ".")
}
\arguments{
\item{path}{the location of the git repository, see details.}

\item{bare}{if true, a Git repository without a working directory is created}

\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}
}
\value{
The path to the Git repository.
}
\description{
Use \code{git_init()} to create a new repository or \code{git_find()} to discover an
existing local repository. \code{git_info()} shows basic information about a
repository, such as the SHA and branch of the current HEAD.
}
\details{
For \code{git_init()} the \code{path} parameter sets the directory of the git repository
to create. If this directory already exists, it must be empty. If it does
not exist, it is created, along with any intermediate directories that don't
yet exist. For \code{git_find()} the \code{path} arguments specifies the directory at
which to start the search for a git repository. If it is not a git repository
itself, then its parent directory is consulted, then the parent's parent, and
so on.
}
\examples{
# directory does not yet exist
r <- tempfile(pattern = "gert")
git_init(r)
git_find(r)

# create a child directory, then a grandchild, then search
r_grandchild_dir <- file.path(r, "aaa", "bbb")
dir.create(r_grandchild_dir, recursive = TRUE)
git_find(r_grandchild_dir)

# cleanup
unlink(r, recursive = TRUE)

# directory exists but is empty
r <- tempfile(pattern = "gert")
dir.create(r)
git_init(r)
git_find(r)

# cleanup
unlink(r, recursive = TRUE)
}
\seealso{
Other git: 
\code{\link{git_archive}},
\code{\link{git_branch}()},
\code{\link{git_commit}()},
\code{\link{git_config}()},
\code{\link{git_diff}()},
\code{\link{git_fetch}()},
\code{\link{git_merge}()},
\code{\link{git_rebase}()},
\code{\link{git_remote}},
\code{\link{git_signature}()},
\code{\link{git_stash}},
\code{\link{git_tag}}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pr.R
\name{git_checkout_pull_request}
\alias{git_checkout_pull_request}
\alias{git_fetch_pull_requests}
\title{GitHub Wrappers}
\usage{
git_checkout_pull_request(pr = 1, remote = NULL, repo = ".")

git_fetch_pull_requests(pr = "*", remote = NULL, repo = ".")
}
\arguments{
\item{pr}{number with PR to fetch or check out. Use \code{"*"} to fetch all
pull requests.}

\item{remote}{Optional. Name of a remote listed in \code{\link[=git_remote_list]{git_remote_list()}}. If
unspecified and the current branch is already tracking branch a remote
branch, that remote is honored. Otherwise, defaults to \code{origin}.}

\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}
}
\description{
Fetch and checkout pull requests.
}
\details{
By default \code{git_fetch_pull_requests} will download all PR branches. To
remove these again simply use \code{git_fetch(prune = TRUE)}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/submodules.R
\name{git_submodule_list}
\alias{git_submodule_list}
\alias{git_submodule_info}
\alias{git_submodule_init}
\alias{git_submodule_set_to}
\alias{git_submodule_add}
\alias{git_submodule_fetch}
\title{Submodules}
\usage{
git_submodule_list(repo = ".")

git_submodule_info(submodule, repo = ".")

git_submodule_init(submodule, overwrite = FALSE, repo = ".")

git_submodule_set_to(submodule, ref, checkout = TRUE, repo = ".")

git_submodule_add(url, path = basename(url), ref = "HEAD", ..., repo = ".")

git_submodule_fetch(submodule, ..., repo = ".")
}
\arguments{
\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}

\item{submodule}{name of the submodule}

\item{overwrite}{overwrite existing entries}

\item{ref}{a branch or tag or hash with}

\item{checkout}{actually switch the contents of the directory to this commit}

\item{url}{full git url of the submodule}

\item{path}{relative of the submodule}

\item{...}{extra arguments for \link{git_fetch} for authentication things}
}
\description{
Interact with submodules in the repository.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gert-package.R
\docType{package}
\name{gert-package}
\alias{gert}
\alias{gert-package}
\title{gert: Simple Git Client for R}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

Simple git client for R based on 'libgit2' with support for SSH and HTTPS remotes. All functions in 'gert' use basic R data types (such as vectors and data-frames) for their arguments and return values. User credentials are shared with command line 'git' through the git-credential store and ssh keys stored on disk or ssh-agent.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/gert/ (website)}
  \item \url{https://github.com/r-lib/gert (devel)}
  \item \url{https://libgit2.org (upstream)}
  \item Report bugs at \url{https://github.com/r-lib/gert/issues}
}

}
\author{
\strong{Maintainer}: Jeroen Ooms \email{jeroen@berkeley.edu} (\href{https://orcid.org/0000-0002-4035-0289}{ORCID})

Other contributors:
\itemize{
  \item Jennifer Bryan \email{jenny@rstudio.com} (\href{https://orcid.org/0000-0002-6983-2759}{ORCID}) [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remotes.R
\name{git_remote}
\alias{git_remote}
\alias{git_remote_list}
\alias{git_remote_add}
\alias{git_remote_remove}
\alias{git_remote_info}
\alias{git_remote_set_url}
\alias{git_remote_set_pushurl}
\alias{git_remote_refspecs}
\title{Git Remotes}
\usage{
git_remote_list(repo = ".")

git_remote_add(url, name = "origin", refspec = NULL, repo = ".")

git_remote_remove(remote, repo = ".")

git_remote_info(remote = NULL, repo = ".")

git_remote_set_url(url, remote = NULL, repo = ".")

git_remote_set_pushurl(url, remote = NULL, repo = ".")

git_remote_refspecs(remote = NULL, repo = ".")
}
\arguments{
\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}

\item{url}{server url (https or ssh)}

\item{name}{unique name for the new remote}

\item{refspec}{optional string with the remote fetch value}

\item{remote}{name of an existing remote. Default \code{NULL} means the remote
from the upstream of the current branch.}
}
\description{
List, add, configure, or remove remotes.
}
\seealso{
Other git: 
\code{\link{git_archive}},
\code{\link{git_branch}()},
\code{\link{git_commit}()},
\code{\link{git_config}()},
\code{\link{git_diff}()},
\code{\link{git_fetch}()},
\code{\link{git_merge}()},
\code{\link{git_rebase}()},
\code{\link{git_repo}},
\code{\link{git_signature}()},
\code{\link{git_stash}},
\code{\link{git_tag}}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch.R
\name{git_branch}
\alias{git_branch}
\alias{git_branch_list}
\alias{git_branch_checkout}
\alias{git_branch_create}
\alias{git_branch_delete}
\alias{git_branch_move}
\alias{git_branch_fast_forward}
\alias{git_branch_set_upstream}
\alias{git_branch_exists}
\title{Git Branch}
\usage{
git_branch(repo = ".")

git_branch_list(local = NULL, repo = ".")

git_branch_checkout(branch, force = FALSE, orphan = FALSE, repo = ".")

git_branch_create(branch, ref = "HEAD", checkout = TRUE, repo = ".")

git_branch_delete(branch, repo = ".")

git_branch_move(branch, new_branch, force = FALSE, repo = ".")

git_branch_fast_forward(ref, repo = ".")

git_branch_set_upstream(upstream, branch = git_branch(repo), repo = ".")

git_branch_exists(branch, local = TRUE, repo = ".")
}
\arguments{
\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}

\item{local}{set TRUE to only check for local branches, FALSE to check for remote
branches. Use NULL to return all branches.}

\item{branch}{name of branch to check out}

\item{force}{ignore conflicts and overwrite modified files}

\item{orphan}{if branch does not exist, checkout unborn branch}

\item{ref}{string with a branch/tag/commit}

\item{checkout}{move HEAD to the newly created branch}

\item{new_branch}{target name of the branch once the move is performed; this name is validated for consistency.}

\item{upstream}{remote branch from \link{git_branch_list}, for example \code{"origin/master"}}
}
\description{
Create, list, and checkout branches.
}
\seealso{
Other git: 
\code{\link{git_archive}},
\code{\link{git_commit}()},
\code{\link{git_config}()},
\code{\link{git_diff}()},
\code{\link{git_fetch}()},
\code{\link{git_merge}()},
\code{\link{git_rebase}()},
\code{\link{git_remote}},
\code{\link{git_repo}},
\code{\link{git_signature}()},
\code{\link{git_stash}},
\code{\link{git_tag}}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{libgit2_config}
\alias{libgit2_config}
\title{Show libgit2 version and capabilities}
\usage{
libgit2_config()
}
\description{
\code{libgit2_config()} reveals which version of libgit2 gert is using and which
features are supported, such whether you are able to use ssh remotes.
}
\examples{
libgit2_config()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{git_config}
\alias{git_config}
\alias{git_config_global}
\alias{git_config_set}
\alias{git_config_global_set}
\title{Get or set Git configuration}
\usage{
git_config(repo = ".")

git_config_global()

git_config_set(name, value, repo = ".")

git_config_global_set(name, value)
}
\arguments{
\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}

\item{name}{Name of the option to set}

\item{value}{Value to set. Must be a string, logical, number or \code{NULL} (to
unset).}
}
\value{
\itemize{
\item \code{git_config()}: a \code{data.frame} of the Git options "in force" in the context
of \code{repo}, one row per option. The \code{level} column reveals whether the
option is determined from global or local config.
\item \code{git_config_global()}: a \code{data.frame}, as for \code{git_config()}, except only
for global Git options.
\item \code{git_config_set()}, \code{git_config_global_set()}: The previous value of
\code{name} in local or global config, respectively. If this option was
previously unset, returns \code{NULL}. Returns invisibly.
}
}
\description{
Get or set Git options, as \verb{git config} does on the command line. \strong{Global}
settings affect all of a user's Git operations (\verb{git config --global}),
whereas \strong{local} settings are scoped to a specific repository (\verb{git config --local}). When both exist, local options always win. Four functions address
the four possible combinations of getting vs setting and global vs. local.\tabular{lll}{
    \tab \strong{local} \tab \strong{global} \cr
   get \tab \code{git_config()} \tab \code{git_config_global()} \cr
   set \tab \code{git_config_set()} \tab \code{git_config_global_set()} \cr
}
}
\note{
All entries in the \code{name} column are automatically normalised to
lowercase (see
\url{https://libgit2.org/libgit2/#HEAD/type/git_config_entry} for details).
}
\examples{
# Set and inspect a local, custom Git option
r <- file.path(tempdir(), "gert-demo")
git_init(r)

previous <- git_config_set("aaa.bbb", "ccc", repo = r)
previous
cfg <- git_config(repo = r)
subset(cfg, level == "local")
cfg$value[cfg$name == "aaa.bbb"]

previous <- git_config_set("aaa.bbb", NULL, repo = r)
previous
cfg <- git_config(repo = r)
subset(cfg, level == "local")
cfg$value[cfg$name == "aaa.bbb"]

unlink(r, recursive = TRUE)

\dontrun{
# Set global Git options
git_config_global_set("user.name", "Your Name")
git_config_global_set("user.email", "your@email.com")
git_config_global()
}
}
\seealso{
Other git: 
\code{\link{git_archive}},
\code{\link{git_branch}()},
\code{\link{git_commit}()},
\code{\link{git_diff}()},
\code{\link{git_fetch}()},
\code{\link{git_merge}()},
\code{\link{git_rebase}()},
\code{\link{git_remote}},
\code{\link{git_repo}},
\code{\link{git_signature}()},
\code{\link{git_stash}},
\code{\link{git_tag}}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rebase.R
\name{git_rebase}
\alias{git_rebase}
\alias{git_rebase_list}
\alias{git_rebase_commit}
\alias{git_reset_hard}
\alias{git_reset_soft}
\alias{git_reset_mixed}
\alias{git_cherry_pick}
\alias{git_ahead_behind}
\title{Cherry-Pick and Rebase}
\usage{
git_rebase_list(upstream = NULL, repo = ".")

git_rebase_commit(upstream = NULL, repo = ".")

git_reset_hard(ref = "HEAD", repo = ".")

git_reset_soft(ref = "HEAD", repo = ".")

git_reset_mixed(ref = "HEAD", repo = ".")

git_cherry_pick(commit, repo = ".")

git_ahead_behind(upstream = NULL, ref = "HEAD", repo = ".")
}
\arguments{
\item{upstream}{branch to which you want to rewind and re-apply your
local commits. The default uses the remote upstream branch with the
current state on the git server, simulating \link{git_pull}.}

\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}

\item{ref}{string with a branch/tag/commit}

\item{commit}{id of the commit to cherry pick}
}
\description{
A cherry-pick applies the changes from a given commit (from another branch)
onto the current branch. A rebase resets the branch to the state of another
branch (upstream) and then re-applies your local changes by cherry-picking
each of your local commits onto the upstream commit history.
}
\details{
\code{git_rebase_list} shows your local commits that are missing from the \code{upstream}
history, and if they conflict with upstream changes. It does so by performing
a rebase dry-run, without committing anything. If there are no conflicts, you
can use \code{git_rebase_commit} to rewind and rebase your branch onto \code{upstream}.
Gert only support a clean rebase; it never leaves the repository in unfinished
"rebasing" state. If conflicts arise, \code{git_rebase_commit} will raise an error
without making changes.
}
\seealso{
Other git: 
\code{\link{git_archive}},
\code{\link{git_branch}()},
\code{\link{git_commit}()},
\code{\link{git_config}()},
\code{\link{git_diff}()},
\code{\link{git_fetch}()},
\code{\link{git_merge}()},
\code{\link{git_remote}},
\code{\link{git_repo}},
\code{\link{git_signature}()},
\code{\link{git_stash}},
\code{\link{git_tag}}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge.R
\name{git_merge}
\alias{git_merge}
\alias{git_merge_stage_only}
\alias{git_merge_find_base}
\alias{git_merge_analysis}
\alias{git_merge_abort}
\title{Merging tools}
\usage{
git_merge(ref, commit = TRUE, squash = FALSE, repo = ".")

git_merge_stage_only(ref, squash = FALSE, repo = ".")

git_merge_find_base(ref, target = "HEAD", repo = ".")

git_merge_analysis(ref, repo = ".")

git_merge_abort(repo = ".")
}
\arguments{
\item{ref}{branch or commit that you want to merge}

\item{commit}{automatically create a merge commit if the merge succeeds without
conflicts. Set this to \code{FALSE} if you want to customize your commit message/author.}

\item{squash}{omits the second parent from the commit, which make the merge a regular
single-parent commit.}

\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}

\item{target}{the branch where you want to merge into. Defaults to current \code{HEAD}.}
}
\description{
Use \code{git_merge} to merge a branch into the current head. Based on how the branches
have diverged, the function will select a fast-forward or merge-commit strategy.
}
\details{
By default \code{git_merge} automatically commits the merge commit upon success.
However if the merge fails with merge-conflicts, or if \code{commit} is set to
\code{FALSE}, the changes are staged and the repository is put in merging state,
and you have to manually run \code{git_commit} or \code{git_merge_abort} to proceed.

Other functions are more low-level tools that are used by \code{git_merge}.
\code{git_merge_find_base} looks up the commit where two branches have diverged
(i.e. the youngest common ancestor). The \code{git_merge_analysis} is used to
test if a merge can simply be fast forwarded or not.

The \code{git_merge_stage_only} function applies and stages changes, without
committing or fast-forwarding.
}
\seealso{
Other git: 
\code{\link{git_archive}},
\code{\link{git_branch}()},
\code{\link{git_commit}()},
\code{\link{git_config}()},
\code{\link{git_diff}()},
\code{\link{git_fetch}()},
\code{\link{git_rebase}()},
\code{\link{git_remote}},
\code{\link{git_repo}},
\code{\link{git_signature}()},
\code{\link{git_stash}},
\code{\link{git_tag}}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/archive.R
\name{git_archive}
\alias{git_archive}
\alias{git_archive_zip}
\title{Git Archive}
\usage{
git_archive_zip(file = NULL, repo = ".")
}
\arguments{
\item{file}{name of the output zip file. Default is returned
by the function}

\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}
}
\value{
path to the zip file that was created
}
\description{
Exports the files in your repository to a zip file that
is returned by the function.
}
\seealso{
Other git: 
\code{\link{git_branch}()},
\code{\link{git_commit}()},
\code{\link{git_config}()},
\code{\link{git_diff}()},
\code{\link{git_fetch}()},
\code{\link{git_merge}()},
\code{\link{git_rebase}()},
\code{\link{git_remote}},
\code{\link{git_repo}},
\code{\link{git_signature}()},
\code{\link{git_stash}},
\code{\link{git_tag}}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{user_is_configured}
\alias{user_is_configured}
\title{Test if a Git user is configured}
\usage{
user_is_configured(repo = ".")
}
\arguments{
\item{repo}{An optional \code{repo}, in the sense of \code{\link[=git_open]{git_open()}}.}
}
\value{
\code{TRUE} if \code{user.name} and \code{user.email} are set locally or globally,
\code{FALSE} otherwise.
}
\description{
This function exists mostly to guard examples that rely on having a user
configured, in order to make commits. \code{user_is_configured()} makes no
distinction between local or global user config.
}
\examples{
user_is_configured()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tag.R
\name{git_tag}
\alias{git_tag}
\alias{git_tag_list}
\alias{git_tag_create}
\alias{git_tag_delete}
\alias{git_tag_push}
\title{Git Tag}
\usage{
git_tag_list(match = "*", repo = ".")

git_tag_create(name, message, ref = "HEAD", repo = ".")

git_tag_delete(name, repo = ".")

git_tag_push(name, ..., repo = ".")
}
\arguments{
\item{match}{pattern to filter tags (use \code{*} for wildcard)}

\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}

\item{name}{tag name}

\item{message}{tag message}

\item{ref}{target reference to tag}

\item{...}{other arguments passed to \link{git_push}}
}
\description{
Create and list tags.
}
\seealso{
Other git: 
\code{\link{git_archive}},
\code{\link{git_branch}()},
\code{\link{git_commit}()},
\code{\link{git_config}()},
\code{\link{git_diff}()},
\code{\link{git_fetch}()},
\code{\link{git_merge}()},
\code{\link{git_rebase}()},
\code{\link{git_remote}},
\code{\link{git_repo}},
\code{\link{git_signature}()},
\code{\link{git_stash}}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signature.R
\name{git_signature}
\alias{git_signature}
\alias{git_signature_default}
\alias{git_signature_parse}
\title{Author Signature}
\usage{
git_signature_default(repo = ".")

git_signature(name, email, time = NULL)

git_signature_parse(sig)
}
\arguments{
\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}

\item{name}{Real name of the committer}

\item{email}{Email address of the committer}

\item{time}{timestamp of class POSIXt or NULL}

\item{sig}{string in proper \code{"First Last <your@email.com>"} format, see details.}
}
\description{
A signature contains the author and timestamp of a commit. Each commit
includes a signature of the author and committer (which can be identical).
}
\details{
A signature string has format \code{"Real Name <email> timestamp tzoffset"}. The
\verb{timestamp tzoffset} piece can be omitted in which case the current local
time is used. If not omitted, \code{timestamp} must contain the number
of seconds since the Unix epoch and \code{tzoffset} is the timezone offset in
\code{hhmm} format (note the lack of a colon separator)
}
\examples{
# Your default user
try(git_signature_default())

# Specify explicit name and email
git_signature("Some committer", "sarah@gmail.com")

# Create signature for an hour ago
(sig <- git_signature("Han", "han@company.com", Sys.time() - 3600))

# Parse a signature
git_signature_parse(sig)
git_signature_parse("Emma <emma@mu.edu>")
}
\seealso{
Other git: 
\code{\link{git_archive}},
\code{\link{git_branch}()},
\code{\link{git_commit}()},
\code{\link{git_config}()},
\code{\link{git_diff}()},
\code{\link{git_fetch}()},
\code{\link{git_merge}()},
\code{\link{git_rebase}()},
\code{\link{git_remote}},
\code{\link{git_repo}},
\code{\link{git_stash}},
\code{\link{git_tag}}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch.R
\name{git_fetch}
\alias{git_fetch}
\alias{git_remote_ls}
\alias{git_push}
\alias{git_clone}
\alias{git_pull}
\title{Push and pull}
\usage{
git_fetch(
  remote = NULL,
  refspec = NULL,
  password = askpass,
  ssh_key = NULL,
  prune = FALSE,
  verbose = interactive(),
  repo = "."
)

git_remote_ls(
  remote = NULL,
  password = askpass,
  ssh_key = NULL,
  verbose = interactive(),
  repo = "."
)

git_push(
  remote = NULL,
  refspec = NULL,
  set_upstream = NULL,
  password = askpass,
  ssh_key = NULL,
  mirror = FALSE,
  force = FALSE,
  verbose = interactive(),
  repo = "."
)

git_clone(
  url,
  path = NULL,
  branch = NULL,
  password = askpass,
  ssh_key = NULL,
  bare = FALSE,
  mirror = FALSE,
  verbose = interactive()
)

git_pull(remote = NULL, rebase = FALSE, ..., repo = ".")
}
\arguments{
\item{remote}{Optional. Name of a remote listed in \code{\link[=git_remote_list]{git_remote_list()}}. If
unspecified and the current branch is already tracking branch a remote
branch, that remote is honored. Otherwise, defaults to \code{origin}.}

\item{refspec}{string with mapping between remote and local refs. Default
uses the default refspec from the remote, which usually fetches all branches.}

\item{password}{a string or a callback function to get passwords for authentication
or password protected ssh keys. Defaults to \link[askpass:askpass]{askpass} which
checks \code{getOption('askpass')}.}

\item{ssh_key}{path or object containing your ssh private key. By default we
look for keys in \code{ssh-agent} and \link[credentials:ssh_credentials]{credentials::ssh_key_info}.}

\item{prune}{delete tracking branches that no longer exist on the remote, or
are not in the refspec (such as pull requests).}

\item{verbose}{display some progress info while downloading}

\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}

\item{set_upstream}{change the branch default upstream to \code{remote}.
If \code{NULL}, this will set the branch upstream only if the push was
successful and if the branch does not have an upstream set yet.}

\item{mirror}{use the \code{--mirror} flag}

\item{force}{use the \code{--force} flag}

\item{url}{remote url. Typically starts with \verb{https://github.com/} for public
repositories, and \verb{https://yourname@github.com/} or \verb{git@github.com/} for
private repos. You will be prompted for a password or pat when needed.}

\item{path}{Directory of the Git repository to create.}

\item{branch}{name of branch to check out locally}

\item{bare}{use the \code{--bare} flag}

\item{rebase}{if TRUE we try to rebase instead of merge local changes. This
is not possible in case of conflicts (you will get an error).}

\item{...}{arguments passed to \link{git_fetch}}
}
\description{
Functions to connect with a git server (remote) to fetch or push changes.
The 'credentials' package is used to handle authentication, the
\href{https://docs.ropensci.org/credentials/articles/intro.html}{credentials vignette}
explains the various authentication methods for SSH and HTTPS remotes.
}
\details{
Use \code{\link[=git_fetch]{git_fetch()}} and \code{\link[=git_push]{git_push()}} to sync a local branch with a remote
branch. Here \code{\link[=git_pull]{git_pull()}} is a wrapper for \code{\link[=git_fetch]{git_fetch()}} which then tries to
\link[=git_branch_fast_forward]{fast-forward} the local branch after fetching.
}
\examples{
{# Clone a small repository
git_dir <- file.path(tempdir(), 'antiword')
git_clone('https://github.com/ropensci/antiword', git_dir)

# Change into the repo directory
olddir <- getwd()
setwd(git_dir)

# Show some stuff
git_log()
git_branch_list()
git_remote_list()

# Add a file
write.csv(iris, 'iris.csv')
git_add('iris.csv')

# Commit the change
jerry <- git_signature("Jerry", "jerry@hotmail.com")
git_commit('added the iris file', author = jerry)

# Now in the log:
git_log()

# Cleanup
setwd(olddir)
unlink(git_dir, recursive = TRUE)
}
}
\seealso{
Other git: 
\code{\link{git_archive}},
\code{\link{git_branch}()},
\code{\link{git_commit}()},
\code{\link{git_config}()},
\code{\link{git_diff}()},
\code{\link{git_merge}()},
\code{\link{git_rebase}()},
\code{\link{git_remote}},
\code{\link{git_repo}},
\code{\link{git_signature}()},
\code{\link{git_stash}},
\code{\link{git_tag}}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{git_diff}
\alias{git_diff}
\alias{git_diff_patch}
\title{Git Diff}
\usage{
git_diff(ref = NULL, repo = ".")

git_diff_patch(ref = NULL, repo = ".")
}
\arguments{
\item{ref}{a reference such as \code{"HEAD"}, or a commit id, or \code{NULL}
to the diff the working directory against the repository index.}

\item{repo}{The path to the git repository. If the directory is not a
repository, parent directories are considered (see \link{git_find}). To disable
this search, provide the filepath protected with \code{\link[=I]{I()}}. When using this
parameter, always explicitly call by name (i.e. \verb{repo = }) because future
versions of gert may have additional parameters.}
}
\description{
View changes in a commit or in the current working directory.
}
\seealso{
Other git: 
\code{\link{git_archive}},
\code{\link{git_branch}()},
\code{\link{git_commit}()},
\code{\link{git_config}()},
\code{\link{git_fetch}()},
\code{\link{git_merge}()},
\code{\link{git_rebase}()},
\code{\link{git_remote}},
\code{\link{git_repo}},
\code{\link{git_signature}()},
\code{\link{git_stash}},
\code{\link{git_tag}}
}
\concept{git}
