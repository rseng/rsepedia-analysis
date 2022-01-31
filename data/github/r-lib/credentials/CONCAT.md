# credentials <img src="man/figures/logo.png" align="right" alt="logo" width="120" height = "139" style = "border: none; float: right;">

*This package is a joint effort from [rOpenSci](https://ropensci.org/) and the [Tidyverse](https://www.tidyverse.org/) team.*

> Tools for Managing SSH and Git Credentials

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/r-lib/credentials?branch=master&svg=true)](https://ci.appveyor.com/project/jeroen/credentials)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/credentials)](http://cran.r-project.org/package=credentials)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/credentials)](http://cran.r-project.org/web/packages/credentials/index.html)


Setup and retrieve HTTPS and SSH credentials for use with 'git' and 
other services. For HTTPS remotes the package interfaces the 'git-credential' 
utility which 'git' uses to store HTTP usernames and passwords. For SSH 
remotes we provide convenient functions to find or generate appropriate SSH 
keys. The package both helps the user to setup a local git installation, and
also provides a back-end for git/ssh client libraries to authenticate with 
existing user credentials.

### Setting your GITHUB_PAT

Automatically populate your `GITHUB_PAT` environment variable from the native git credential store. The credential manager will safely prompt the user for credentials when needed.

```r
credentials::set_github_pat()
```

Use this function in your `.Rprofile` if you want to automatically set `GITHUB_PAT` for each R session, without hardcoding your secret in plain text.

### Manage HTTPS credentials

Load or prompt the user for GitHub username and password:

```r
library(credentials)
git_credential_ask('https://github.com')
```

See which credential helper back-end your `git-credential` store is using:

```r
credentials::credential_helper_get()
```

### Manage SSH keys

Lookup the appropriate key, or prompt the user to generate one:

```r
library(credentials)
ssh_key_info()
```

You can copy-paste the public key directly to your [GitHub profile](https://github.com/settings/ssh/new)!

## For developers

Use the openssl package to read the user private key in R for encryption / signatures: 

```r
user <- ssh_key_info()
key <- ssh_read_key(user$key)
openssl::write_pem(key)
```

---
title: "Managing SSH and Git Credentials in R"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Managing SSH and Git Credentials in R}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The `credentials` package contains tools for configuring and retrieving SSH and HTTPS credentials for use with `git` or other services. It helps users to setup their git installation, and also provides a back-end for packages to authenticate with existing user credentials.

## Two types of remotes

```{r, echo=FALSE}
has_git <- credentials:::has_git_cmd() && .Platform$OS.type != 'windows'
delete_git_config_on_exit <- !file.exists('~/.gitconfig')
credentials:::set_default_cred_helper()

library <- function(package){
  withCallingHandlers(base::library(credentials), 
                      packageStartupMessage = function(e) {
                        cat(e$message)
                        invokeRestart("muffleMessage")
                      })
}  
```


```{r}
library(credentials)
```

Git supports two types of remotes: SSH and HTTPS. These two use completely distinct authentication systems. 


|         | SSH REMOTES      | HTTPS REMOTES          |
|---------|------------------|------------------------|
| __url__     | `git@server.com` | `https://server.com`   |
| __authentication__    | Personal SSH key          | Password or PAT         |
| __stored in__ | file `id_rsa` or `ssh-agent` | `git credential` store |

For HTTPS remotes, git authenticates with a username + password. With GitHub, instead of a password, you can also use a [Personal Access Token](https://github.com/settings/tokens) (PAT). PATs are preferred because they provide more granular control of permissions (via the PAT's scopes), you can have many of them for different purposes, you can give them informative names, and they can be selectively revoked. Note that, if you use 2FA with GitHub, you must authenticate with a PAT if you use the HTTPS protocol.

For SSH remotes, git shells out to `ssh` on your system, and lets `ssh` take care of authentication. This means you have to setup an ssh key (usually `~/.ssh/id_rsa`) which you then [add to your git profile](https://github.com/settings/ssh/new).

### Special note for Windows

Windows does not include a native `git` installation by default. We recommended to use the latest version of [Git for Windows](https://git-scm.com/download/win). This bundle also includes `ssh` and [git credential manager for windows](https://github.com/Microsoft/Git-Credential-Manager-for-Windows) which is all you need.

Important: ssh keys are stored in your home directory for example: `C:\Users\Jeroen\.ssh\id_rsa`, and __not in the Documents folder__ (which is what R treats as the home sometimes). The `ssh_home()` function shows the correct `.ssh` directory on all platforms.


## Part 1: Storing HTTPS credentials

HTTPS remotes do not always require authentication. You can clone from a public repository without providing any credentials. But for pushing, or private repositories, `git` will prompt for a username/password.

```
git clone https://github.com/jeroen/jsonlite
```

To save you from entering your password over and over, git includes a [credential helper](https://git-scm.com/docs/gitcredentials). It has two modes:

 - `cache`: Cache credentials in memory for a short period of time.
 - `store`: Store credentials permanently in your operating system password manager.
 
To see which helper is configured for a given repo, run:

```{r, eval=has_git}
credential_helper_get()
```

Most `git` installations default to `store` if supported because it is more convenient and secure. However the look and policy of the git credential store for entering and retrieving passwords can vary by system, because it uses the OS native password manager.

### Accessing the HTTPS Credential Store from R

The `credentials` R package provides a wrapper around the `git credential` command line API for reading and saving credentials. The `git_credential_ask()` function looks up suitable credentials for a given URL from the store. If no credentials are available, it will attempt to prompt the user for credentials and return those instead.

```{r, echo=FALSE, eval=has_git}
# This hack may not work on MacOS server where cred helper is osxkeychain 
# which always requires user interaction. Hence error=TRUE in the next block.
example <- list(protocol = "https", host = "example.com",
  username = "jeroen", password = "supersecret")
credential_approve(example)
```

```{r error=TRUE, eval=has_git}
library(credentials)
git_credential_ask('https://example.com')
```

The function `git_credential_update()` looks similar but it behaves slightly different: it first removes existing credentials from the store (if any), then prompts the user for a new username/password, and saves these to the store.

```r
# This should always prompt for new credentials
git_credential_update('https://example.com')
```

In a terminal window this will result in an interactive password prompt. In Windows the user might see something like this (depending on the version of Windows and git configuration):

<img style="width:690px; padding:0;" src="wincred.png">

```{r, echo=FALSE, eval=has_git}
credential_reject(list(protocol = "https", host = "example.com"))
```

### Setting your GITHUB_PAT

Automatically populate your `GITHUB_PAT` environment variable from the native git credential store. The credential manager will safely prompt the user for credentials when needed.

```r
credentials::set_github_pat()
## Using GITHUB_PAT from Jeroen Ooms (credential helper: osxkeychain)
```

Use this function in your `.Rprofile` if you want to automatically set `GITHUB_PAT` for each R session, without hardcoding your secrets in plain text, such as in your `.Renviron` file.

### Non-interactive use

Retrieving credentials is by definition interactive, because the user may be required to enter a password or unlock the system keychain. However, saving or deleting credentials can sometimes be done non-interactively, but this depends on which credential helper is used.

The manual page for `credential_approve` and `credential_reject` has more details about how to call the basic git credential api. 


## Part 2: Managing SSH Keys

```{r, echo = FALSE}
ssh_key_info <- function(){
  try({
    out <- credentials:::ssh_key_info(auto_keygen = FALSE)
    out$pubkey = paste(substring(out$pubkey, 1, 80), "...")
    out
  })
}
```

For SSH remotes, git does not handle authentication itself. Git simply shells out to `ssh` on your system and uses your standard user ssh configuration. Hence authenticating with SSH git remotes comes down to setting up your ssh keys and copying these to your profile.

The `credentials` package provides a few utility functions to make this easier. The `ssh_key_info()` function calls out to look up which key `ssh` uses to connect to a given server. This is usually `~/.ssh/id_rsa` unless you have a fancy custom ssh configuration.

```{r}
ssh_key_info()
```

The output shows both the path to your (private) key as well as the ssh pubkey string. The latter is what you have to enter in your [GitHub profile](https://github.com/settings/ssh/new) to associate this key with your user account. You will then be automatically authenticated when using GitHub SSH remotes.

### Generating a key

To use SSH you need a personal key, which is usually stored in `~/.ssh/id_rsa`. If you do not have a key yet, the `ssh_key_info()` function will automatically ask if you want to generate one.

<img style="width:690px; padding:0;" src="keygen.png">

You can also generate a key manually elsewhere using the `ssh_keygen()` function.

```{r echo=FALSE}
if(isTRUE(delete_git_config_on_exit))
  unlink("~/.gitconfig")
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/credential-helpers.R
\name{credential_helper}
\alias{credential_helper}
\alias{credential_helper_list}
\alias{credential_helper_get}
\alias{credential_helper_set}
\title{Credential Helpers}
\usage{
credential_helper_list()

credential_helper_get(global = FALSE)

credential_helper_set(helper, global = FALSE)
}
\arguments{
\item{global}{if FALSE the setting is done per git repository, if
TRUE it is in your global user git configuration.}

\item{helper}{string with one of the supported helpers from \link{credential_helper_list}}
}
\description{
Git supports several back-end stores for HTTPS credentials called
helpers. Default helpers include \code{cache} and \code{store}, see the
\href{https://git-scm.com/docs/gitcredentials}{git-credentials} manual
page for details.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssh-keys.R
\name{ssh_credentials}
\alias{ssh_credentials}
\alias{ssh_key_info}
\alias{ssh_keygen}
\alias{ssh_setup_github}
\alias{ssh_home}
\alias{ssh_agent_add}
\alias{ssh_update_passphrase}
\alias{ssh_read_key}
\title{Managing Your SSH Key}
\usage{
ssh_key_info(host = NULL, auto_keygen = NA)

ssh_keygen(file = ssh_home("id_rsa"))

ssh_setup_github()

ssh_home(file = NULL)

ssh_agent_add(file = NULL)

ssh_update_passphrase(file = ssh_home("id_rsa"))

ssh_read_key(file = ssh_home("id_rsa"), password = askpass)
}
\arguments{
\item{host}{target host (only matters if you have configured specific keys per host)}

\item{auto_keygen}{if \code{TRUE} automatically generates a key if none exists yet.
Default \code{NA} is to prompt the user what to.}

\item{file}{destination path of the private key. For the public key, \code{.pub}
is appended to the filename.}

\item{password}{a passphrase or callback function}
}
\description{
Utility functions to find or generate your SSH key for use with git remotes
or other ssh servers.
}
\details{
Use \code{\link[=ssh_key_info]{ssh_key_info()}} to find the appropriate key file on your system to connect with a
given target host. In most cases this will simply be \code{ssh_home('id_rsa')} unless
you have configured ssh to use specific keys for specific hosts.

To use your key to authenticate with GitHub, copy the pubkey from \code{ssh_key_info()} to
your profile: \url{https://github.com/settings/ssh/new}.

If this is the first time you use ssh, \link{ssh_keygen} can help generate a key and
save it in the default location. This will also automatically opens the above Github
page in your browser where you can add the key to your profile.

\code{ssh_read_key} reads a private key and caches the result (in memory) for the
duration of the R session. This prevents having to enter the key passphrase many
times. Only use this if \code{ssh-agent} is not available (i.e. Windows)
}
\seealso{
Other credentials: 
\code{\link{http_credentials}}
}
\concept{credentials}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/credential-api.R
\name{credential_api}
\alias{credential_api}
\alias{credential_fill}
\alias{credential_approve}
\alias{credential_reject}
\title{Retrieve and store git HTTPS credentials}
\usage{
credential_fill(cred, verbose = TRUE)

credential_approve(cred, verbose = TRUE)

credential_reject(cred, verbose = TRUE)
}
\arguments{
\item{cred}{named list with at least fields \code{protocol} and \code{host} and
optionally also \code{path}, \code{username} ,\code{password}.}

\item{verbose}{emit some useful output about what is happening}
}
\description{
Low-level wrappers for the \href{https://git-scm.com/docs/git-credential}{git-credential}
command line tool. Try the user-friendly \link{git_credential_ask}
and \link{git_credential_update} functions first.
}
\details{
The \link{credential_fill} function looks up credentials for a given host, and
if none exists it will attempt to prompt the user for new credentials. Upon
success it returns a list with the same \code{protocol} and \code{host} fields as the
\code{cred} input, and additional \code{username} and \code{password} fields.

After you have tried to authenticate the provided credentials, you can report
back if the credentials were valid or not. Call \link{credential_approve} and
\link{credential_reject} with the \code{cred} that was returned by \link{credential_fill}
in order to validate or invalidate a credential from the store.

Because git credential interacts with the system password manager, the appearance
of the prompts vary by OS and R frontend.  Note that \link{credential_fill} should
only be used interactively, because it may require the user to enter credentials
or unlock the system keychain. On the other hand \link{credential_approve} and
\link{credential_reject} are non-interactive and could be used to save or delete
credentials in a scripted program. However note that some credential helpers
(e.g. on Windows) have additional security restrictions that limit use of
\link{credential_approve} and \link{credential_reject} to credentials that were actually
entered by the user via \link{credential_fill}. Here it is not possible at all to
update the credential store without user interaction.
}
\examples{
\donttest{
# Insert example cred
example <- list(protocol = "https", host = "example.org",
  username = "test", password = "secret")
credential_approve(example)

# Retrieve it from the store
cred <- credential_fill(list(protocol = "https", host = "example.org", path = "/foo"))
print(cred)

# Delete it
credential_reject(cred)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/github-pat.R
\name{set_github_pat}
\alias{set_github_pat}
\title{Set your Github Personal Access Token}
\usage{
set_github_pat(force_new = FALSE, validate = interactive(), verbose = validate)
}
\arguments{
\item{force_new}{forget existing pat, always ask for new one.}

\item{validate}{checks with the github API that this token works. Defaults to
\code{TRUE} only in an interactive R session (not when running e.g. CMD check).}

\item{verbose}{prints a message showing the credential helper and PAT owner.}
}
\value{
Returns \code{TRUE} if a valid GITHUB_PAT was set, and FALSE if not.
}
\description{
Populates the \code{GITHUB_PAT} environment variable using the \link[=http_credentials]{git_credential}
manager, which \code{git} itself uses for storing passwords. The credential manager
returns stored credentials if available, and securely prompt the user for
credentials when needed.
}
\details{
Packages that require a \code{GITHUB_PAT} can call this function to automatically
set the \code{GITHUB_PAT} when needed. Users may call this function in their
\link[=Startup]{.Rprofile} script to automatically set \code{GITHUB_PAT} for each R
session without hardcoding any tokens on disk in plain-text.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/http-credential.R
\name{http_credentials}
\alias{http_credentials}
\alias{git_credential_ask}
\alias{credentials}
\alias{git_credential_update}
\alias{git_credential_forget}
\title{Load and store git HTTPS credentials}
\usage{
git_credential_ask(url = "https://github.com", save = TRUE, verbose = TRUE)

git_credential_update(url = "https://github.com", verbose = TRUE)

git_credential_forget(url = "https://github.com", verbose = TRUE)
}
\arguments{
\item{url}{target url, possibly including username or path}

\item{save}{in case the user is prompted for credentials, attempt to
remember them.}

\item{verbose}{print errors from \verb{git credential} to stdout}
}
\description{
This requires you have the \code{git} command line program installed.The
\link{git_credential_ask} function looks up a suitable username/password
from the \href{https://git-scm.com/docs/gitcredentials}{\code{git-credential} store}.
If none are available it will prompt the user for credentials which
may be saved the store. On subsequent calls for the same URL, the
function will then return the stored credentials without prompting
the user.
}
\details{
The appearance and security policy of the credential store depends
on your version of git, your operating system, your R frontend and
which \link{credential_helper} is used. On Windows and MacOS the credentials
are stored in the system password manager by default.

It should be assumed that reading credentials always involves user
interaction. The user may be asked to unlock the system keychain or
enter new credentials. In reality, user interaction is usually only
required on the first authentication attempt, but the security policy
of most credential helpers prevent you from programmatically testing
if the credentials are already unlocked.
}
\seealso{
Other credentials: 
\code{\link{ssh_credentials}}
}
\concept{credentials}
