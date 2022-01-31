# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who contribute through reporting issues, posting feature requests, updating documentation, submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for everyone, regardless of level of experience, gender, gender identity and expression, sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or imagery, derogatory comments or personal attacks, trolling, public or private harassment, insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant (https://contributor-covenant.org), version 1.0.0, available at https://contributor-covenant.org/version/1/0/0/


# cyphr

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R build status](https://github.com/ropensci/cyphr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/cyphr/actions)
[![codecov.io](https://codecov.io/github/ropensci/cyphr/coverage.svg?branch=master)](https://codecov.io/github/ropensci/cyphr?branch=master)
[![](https://www.r-pkg.org/badges/version/cyphr)](https://cran.r-project.org/package=cyphr)
[![](https://badges.ropensci.org/114_status.svg)](https://github.com/ropensci/onboarding/issues/114)
<!-- badges: end -->

High-level functions for supporting encryption and decryption of data from R.  This allows secure storage and exchange of information, while trying to keep the encryption/decryption code from taking over your analyses.  `cyphr` wraps the lower-level support from [`sodium`](https://github.com/jeroenooms/sodium) and [`openssl`](https://github.com/jeroenooms/openssl).  This package is designed to be easy to use, rather than the most secure thing (you're using R, remember - for examples of what `cyphr` can't protect against see [`jammr`](https://github.com/Ironholds/jammr), [`rpwnd`](https://github.com/hrbrmstr/rpwnd) and [`evil.R`](https://github.com/romainfrancois/evil.R).)

`cyphr` provides high level functions to:

* Encrypt and decrypt
  * **strings**: `encrypt_string` / `decrypt_string`
  * **objects**: `encrypt_object` / `decrypt_object`
  * **raw data**: `encrypt_data` / `decrypt_data`
  * **files**: `encrypt_file` / `decrypt_file`
* User-friendly wrappers (`encrypt` and `decrypt`) around R's file reading and writing functions that enable transparent encryption (support included for `readRDS`/`writeRDS`, `read.csv`/`write.csv`, etc).

The package aims to make encrypting and decrypting as easy as


```r
cyphr::encrypt(save.csv(dat, "file.csv"), key)
```

and


```r
dat <- cyphr::decrypt(read.csv("file.csv", stringsAsFactors = FALSE), key)
```

In addition, the package implements a workflow that allows a group to securely share data by encrypting it with a shared ("symmetric") key that is in turn encrypted with each users ssh keys.  The use case is a group of researchers who are collaborating on a dataset that cannot be made public, for example containing sensitive data.  However, they have decided or need to store it in a setting that they are not 100% confident about the security of the data.  So encrypt the data at each read/write.

## Installation

Install `cyphr` from CRAN with


```r
install.packages("cyphr")
```

To install a development version from github, you can use `remotes`

To install `cyphr` from github:

```r
remotes::install_github("ropensci/cyphr", upgrade = FALSE)
```

## Scope

The scope of the package is to protect data that has been saved to disk.  It is not designed to stop an attacker targeting the R process itself to determine the contents of sensitive data.  The package does try to prevent you accidentally saving to disk the contents of sensitive information, including the keys that could decrypt such information.

## Objects to handle keys:

Decide on a style of encryption and create a key object

* `key_sodium`: Symmetric encryption, using [sodium](https://cran.r-project.org/package=sodium) -- everyone shares the same key (which must be kept secret!) and can encrypt and decrypt data with it.  This is used as a building block but is inflexible because of the need to keep the key secret.
* `key_openssl`: Symmetric encryption using [openssl](https://cran.r-project.org/package=openssl)
* `keypair_sodium`: Public key encryption with [sodium](https://cran.r-project.org/package=sodium) -- this lets people encrypt messages using your public key that only you can read using your private key.
* `keypair_openssl`: Public key encryption, using [openssl](https://cran.r-project.org/package=openssl), which has the big advantage that many people already have compatible (ssh) keys in standard places with standard file formats (see `?encrypt_envelope` in the the `openssl` package).

`cyphr` does not include wrappers for key generation for sodium - sodium keys do not have a file format:  So a secret symmetric key in `sodium` might be:


```r
k <- sodium::keygen()
k
```

```
##  [1] 48 35 a2 6c 05 27 65 75 cb 08 01 de 76 8f 71 fe 3f d7 e4 7a df bf d8 e7 08
## [26] d5 fb e9 61 c8 5f d1
```

With this key we can create the `key_sodium` object:


```r
key <- cyphr::key_sodium(k)
key
```

```
## <cyphr_key: sodium>
```

If the key was saved to file that would work too:

If you load a password protected ssh key you will be prompted for your passphrase.  `cyphr` will ensure that this is not echoed onto the console.


```r
key <- cyphr::key_openssl()
## Please enter private key passphrase:
key
```

## Encrypt / decrypt a file

If you have files that already exist and you want to encrypt or decrypt, the functions `cyphr::encrypt_file` and `cyphr::decrypt_file` will do that (these are workhorse functions that are used internally throughout the package)


```r
saveRDS(iris, "myfile")
cyphr::encrypt_file("myfile", key, "myfile.encrypted")
```

The file is encrypted now:


```r
readRDS("myfile.encrypted")
```

```
## Error in readRDS("myfile.encrypted"): unknown input format
```

Decrypt the file and read it:


```r
cyphr::decrypt_file("myfile.encrypted", key, "myfile.clear")
identical(readRDS("myfile.clear"), iris)
```

```
## [1] TRUE
```



## Wrappers around R's file functions

Encrypting files like the above risks leaving a cleartext (i.e., unencrypted) version around.  If you want to wrap the output of something like `write.csv` or `saveRDS` you really have no choice but to write out the file first, encrypt it, and delete the clear version.  Making sure that this happens even if a step fails is error prone and takes a surprising number of repetitive lines of code.

Alternatively, to encrypt the output of a file producing command, just wrap it in `cyphr::encrypt`


```r
cyphr::encrypt(saveRDS(iris, "myfile.rds"), key)
```

Then to decrypt the a file to feed into a file consuming command, wrap it in `cyphr::decrypt`


```r
dat <- cyphr::decrypt(readRDS("myfile.rds"), key)
```

The round-trip preserves the data:

```r
identical(dat, iris) # yay
```

```
## [1] TRUE
```

But without the key, it cannot be read:

```r
readRDS("myfile.rds")
```

```
## Error in readRDS("myfile.rds"): unknown input format
```



The above commands work through computing on the language, rewriting the `readRDS` and `saveRDS` commands.  Commands for reading and writing tabular and plain text files (`read.csv`, `readLines`, etc) are also supported, and the way the rewriting is done is designed to be extensible.

The argument to the wrapped functions can be connection objects.  In this case the *actual* command is written to a file and the contents of that file are encrypted and written to the connection.  When reading/writing multiple objects from/to a single connection though, this is likely to go very badly.

### Supporting additional functions

The functions supported so far are:

* `readLines` / `writeLines`
* `readRDS` / `writeRDS`
* `read` / `save`
* `read.table` / `write.table`
* `read.csv` / `read.csv2` / `write.csv`
* `read.delim` / `read.delim2`

However, there are bound to be more functions that could be useful to add here (e.g., `readxl::read_excel`).  Either pass the name of the file argument to `cyphr::encrypt` / `cyphr::decrypt` as

```r
cyphr::decrypt(readxl::read_excel("myfile.xlsx"), key, file_arg = "path")
```

or *register* the function with the package using `rewrite_register`:

```r
cyphr::rewrite_register("readxl", "read_excel", "path")
```

Then you can use

```r
cyphr::decrypt(readxl::read_excel("myfile.xlsx"), key)
```

to decrypt the file (these are equivalent, but the former will likely be more convenient if you're only dealing with a couple of files, the latter will be more convenient if you are dealing with many).

## Collaborating with encrypted data

Even with high-level functions to ease encrypting and decrypting things given a key, there is some work to be done to distribute a set of keys across a group of people who are working together so that everyone can encrypt and decrypt the data but so that the keys themselves are not compromised.

The package contains support for a group of people are working on a sensitive data set.  The data will be stored with a symmetric key.  However, we never actually store the key directly, instead we'll store a copy for each user that is encrypted with the user's key.  Any user with access to the data can authorise another user to access the data.  This is described in more detail in the [vignette](https://docs.ropensci.org/cyphr/articles/data.html) (in R: `vignette("data", package = "cyphr")`).

## Why are wrappers needed?

The low level functions in `sodium` and `openssl` work with raw data, for generality.  Few users encounter raw vectors in their typical use of R, so these require serialisation.  Most of the encryption involves a little extra random data (the "nonce" in `sodium` and similar additional pieces with `openssl`).  These need storing with the data, and then separating from the data when decryption happens.

## Licence

MIT © [Rich FitzJohn](https://github.com/richfitz).

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/cyphr/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# cyphr 1.0.6

* Added wrappers for `readxl::read_excel` and `writexl::write_xlsx` (reside-109)
* Support for custom messages during requesting and authorising access to data (reside-108)

# cyphr 1.0.1

* First CRAN release

# cyphr 0.2.0

* Authenticated encryption (with signed messages) is now supported for openssl, and is enabled by default.  Along with this the pack format for openssl has changed which will break *all* existing uses of the package (I don't believe there are any)

# cyphr 0.1.0

* Initial prototype, as sent to rOpenSci onboarding for review
## cyphr internal structures

This directory is used by `cyphr` and contains:

* `test` - a small file used to test whether encryption/decrpytion is possible with your key
* Files in `keys/` are encrypted copies of the (symmetric) data key, encrypted with different users' public keys.  The filename is based on the fingerprint of the public key (see `?openssl::fingerprint`)
* Files in `requests/` which are pending requests for access

### Templates

You can edit the files `template/request` and `template/authorise` and they will be used to provide feedback when requesting access and authorising keys.  Because this step requires some out-of-band communication this can be useful.

Within the `template/request` template the string `$HASH` will be substituted for your key's hash, and within `template/authorise` the string `$USERS` will contain the usernames of added users.
```{r, setup, echo = FALSE, message = FALSE}
knitr::opts_chunk$set(
  tidy = FALSE,
  error = FALSE,
  fig.width = 8,
  fig.height = 8)
```

# cyphr

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R build status](https://github.com/ropensci/cyphr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/cyphr/actions)
[![codecov.io](https://codecov.io/github/ropensci/cyphr/coverage.svg?branch=master)](https://codecov.io/github/ropensci/cyphr?branch=master)
[![](https://www.r-pkg.org/badges/version/cyphr)](https://cran.r-project.org/package=cyphr)
[![](https://badges.ropensci.org/114_status.svg)](https://github.com/ropensci/onboarding/issues/114)
<!-- badges: end -->

High-level functions for supporting encryption and decryption of data from R.  This allows secure storage and exchange of information, while trying to keep the encryption/decryption code from taking over your analyses.  `cyphr` wraps the lower-level support from [`sodium`](https://github.com/jeroenooms/sodium) and [`openssl`](https://github.com/jeroenooms/openssl).  This package is designed to be easy to use, rather than the most secure thing (you're using R, remember - for examples of what `cyphr` can't protect against see [`jammr`](https://github.com/Ironholds/jammr), [`rpwnd`](https://github.com/hrbrmstr/rpwnd) and [`evil.R`](https://github.com/romainfrancois/evil.R).)

`cyphr` provides high level functions to:

* Encrypt and decrypt
  * **strings**: `encrypt_string` / `decrypt_string`
  * **objects**: `encrypt_object` / `decrypt_object`
  * **raw data**: `encrypt_data` / `decrypt_data`
  * **files**: `encrypt_file` / `decrypt_file`
* User-friendly wrappers (`encrypt` and `decrypt`) around R's file reading and writing functions that enable transparent encryption (support included for `readRDS`/`writeRDS`, `read.csv`/`write.csv`, etc).

The package aims to make encrypting and decrypting as easy as

```{r, eval = FALSE}
cyphr::encrypt(save.csv(dat, "file.csv"), key)
```

and

```{r, eval = FALSE}
dat <- cyphr::decrypt(read.csv("file.csv", stringsAsFactors = FALSE), key)
```

In addition, the package implements a workflow that allows a group to securely share data by encrypting it with a shared ("symmetric") key that is in turn encrypted with each users ssh keys.  The use case is a group of researchers who are collaborating on a dataset that cannot be made public, for example containing sensitive data.  However, they have decided or need to store it in a setting that they are not 100% confident about the security of the data.  So encrypt the data at each read/write.

## Installation

Install `cyphr` from CRAN with

```{r, eval = FALSE}
install.packages("cyphr")
```

To install a development version from github, you can use `remotes`

To install `cyphr` from github:
```{r, eval = FALSE}
remotes::install_github("ropensci/cyphr", upgrade = FALSE)
```

## Scope

The scope of the package is to protect data that has been saved to disk.  It is not designed to stop an attacker targeting the R process itself to determine the contents of sensitive data.  The package does try to prevent you accidentally saving to disk the contents of sensitive information, including the keys that could decrypt such information.

## Objects to handle keys:

Decide on a style of encryption and create a key object

* `key_sodium`: Symmetric encryption, using [sodium](https://cran.r-project.org/package=sodium) -- everyone shares the same key (which must be kept secret!) and can encrypt and decrypt data with it.  This is used as a building block but is inflexible because of the need to keep the key secret.
* `key_openssl`: Symmetric encryption using [openssl](https://cran.r-project.org/package=openssl)
* `keypair_sodium`: Public key encryption with [sodium](https://cran.r-project.org/package=sodium) -- this lets people encrypt messages using your public key that only you can read using your private key.
* `keypair_openssl`: Public key encryption, using [openssl](https://cran.r-project.org/package=openssl), which has the big advantage that many people already have compatible (ssh) keys in standard places with standard file formats (see `?encrypt_envelope` in the the `openssl` package).

`cyphr` does not include wrappers for key generation for sodium - sodium keys do not have a file format:  So a secret symmetric key in `sodium` might be:

```{r}
k <- sodium::keygen()
k
```

With this key we can create the `key_sodium` object:

```{r}
key <- cyphr::key_sodium(k)
key
```

If the key was saved to file that would work too:

If you load a password protected ssh key you will be prompted for your passphrase.  `cyphr` will ensure that this is not echoed onto the console.

```{r, eval = FALSE}
key <- cyphr::key_openssl()
## Please enter private key passphrase:
key
```

## Encrypt / decrypt a file

If you have files that already exist and you want to encrypt or decrypt, the functions `cyphr::encrypt_file` and `cyphr::decrypt_file` will do that (these are workhorse functions that are used internally throughout the package)

```{r}
saveRDS(iris, "myfile")
cyphr::encrypt_file("myfile", key, "myfile.encrypted")
```

The file is encrypted now:

```{r, error = TRUE}
readRDS("myfile.encrypted")
```

Decrypt the file and read it:

```{r}
cyphr::decrypt_file("myfile.encrypted", key, "myfile.clear")
identical(readRDS("myfile.clear"), iris)
```

```{r, echo = FALSE, results = "hide"}
file.remove("myfile", "myfile.encrypted", "myfile.clear")
```

## Wrappers around R's file functions

Encrypting files like the above risks leaving a cleartext (i.e., unencrypted) version around.  If you want to wrap the output of something like `write.csv` or `saveRDS` you really have no choice but to write out the file first, encrypt it, and delete the clear version.  Making sure that this happens even if a step fails is error prone and takes a surprising number of repetitive lines of code.

Alternatively, to encrypt the output of a file producing command, just wrap it in `cyphr::encrypt`

```{r}
cyphr::encrypt(saveRDS(iris, "myfile.rds"), key)
```

Then to decrypt the a file to feed into a file consuming command, wrap it in `cyphr::decrypt`

```{r}
dat <- cyphr::decrypt(readRDS("myfile.rds"), key)
```

The round-trip preserves the data:
```{r}
identical(dat, iris) # yay
```

But without the key, it cannot be read:
```{r, error = TRUE}
readRDS("myfile.rds")
```

```{r, echo = FALSE, results = "hide"}
file.remove("myfile.rds")
```

The above commands work through computing on the language, rewriting the `readRDS` and `saveRDS` commands.  Commands for reading and writing tabular and plain text files (`read.csv`, `readLines`, etc) are also supported, and the way the rewriting is done is designed to be extensible.

The argument to the wrapped functions can be connection objects.  In this case the *actual* command is written to a file and the contents of that file are encrypted and written to the connection.  When reading/writing multiple objects from/to a single connection though, this is likely to go very badly.

### Supporting additional functions

The functions supported so far are:

* `readLines` / `writeLines`
* `readRDS` / `writeRDS`
* `read` / `save`
* `read.table` / `write.table`
* `read.csv` / `read.csv2` / `write.csv`
* `read.delim` / `read.delim2`

However, there are bound to be more functions that could be useful to add here (e.g., `readxl::read_excel`).  Either pass the name of the file argument to `cyphr::encrypt` / `cyphr::decrypt` as

```r
cyphr::decrypt(readxl::read_excel("myfile.xlsx"), key, file_arg = "path")
```

or *register* the function with the package using `rewrite_register`:

```r
cyphr::rewrite_register("readxl", "read_excel", "path")
```

Then you can use

```r
cyphr::decrypt(readxl::read_excel("myfile.xlsx"), key)
```

to decrypt the file (these are equivalent, but the former will likely be more convenient if you're only dealing with a couple of files, the latter will be more convenient if you are dealing with many).

## Collaborating with encrypted data

Even with high-level functions to ease encrypting and decrypting things given a key, there is some work to be done to distribute a set of keys across a group of people who are working together so that everyone can encrypt and decrypt the data but so that the keys themselves are not compromised.

The package contains support for a group of people are working on a sensitive data set.  The data will be stored with a symmetric key.  However, we never actually store the key directly, instead we'll store a copy for each user that is encrypted with the user's key.  Any user with access to the data can authorise another user to access the data.  This is described in more detail in the [vignette](https://docs.ropensci.org/cyphr/articles/data.html) (in R: `vignette("data", package = "cyphr")`).

## Why are wrappers needed?

The low level functions in `sodium` and `openssl` work with raw data, for generality.  Few users encounter raw vectors in their typical use of R, so these require serialisation.  Most of the encryption involves a little extra random data (the "nonce" in `sodium` and similar additional pieces with `openssl`).  These need storing with the data, and then separating from the data when decryption happens.

## Licence

MIT © [Rich FitzJohn](https://github.com/richfitz).

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/cyphr/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Data Encryption"
author: "Rich FitzJohn"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Data Encryption}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

**The scenario:**

A group of people are working on a sensitive data set that for
practical reasons needs to be stored in a place that we're not 100%
happy with the security (e.g., Dropbox), or we're concerned that
files stored in plain text on users computers (e.g. laptops) may
lead to the data being compromised.

If the data can be stored encrypted but everyone in the group can
still read and write the data then we've improved the situation
somewhat.  But organising for everyone to get a copy of the key to
decrypt the data files is non-trivial.  The workflow described here
aims to simplify this procedure using lower-level functions in the
`cyphr` package.

The general procedure is this:

1. A person will set up a set of personal keys and a key for the
data.  The data key will be encrypted with their personal key so
they have access to the data but nobody else does.  At this point
the data can be encrypted.

2. Additional users set up personal keys and request access to the
data.  Anyone with access to the data can grant access to anyone
else.

Before doing any of this, everyone needs to have ssh keys set up.
By default the package will use your ssh keys found at "~/.ssh";
see the main package vignette for how to use this.

For clarity here we will generate two sets of key pairs for two
actors Alice and Bob:
``` {r }
path_key_alice <- cyphr::ssh_keygen(password = FALSE)
path_key_bob <- cyphr::ssh_keygen(password = FALSE)
```

These would ordinarily be on different machines (nobody has access
to anyone else's private key) and they would be password protected.
In the function calls below, all the `path_user` arguments would be
omitted.

We'll store data in the directory `data`; at present there is
nothing there (this is in a temporary directory for compliance with
CRAN policies but would ordinarily be somewhere persistent and
under version control ideally).
``` {r }
data_dir <- file.path(tempdir(), "data")
dir.create(data_dir)
dir(data_dir)
```

**First**, create a personal set of keys.  These will be shared
across all projects and stored away from the data.  Ideally one
would do this with `ssh-keygen` at the command line, following one
of the many guides available.  A utility function `ssh_keygen`
(which simply calls `ssh-keygen` for you) is available in this
package though.  You will need to generate a key on each computer
you want access from.  Don't copy the key around.  If you lose your
user key you will lose access to the data!

**Second**, create a key for the data and encrypt that key with
your personal key.  Note that the data key is never stored directly
- it is always stored encrypted by a personal key.
``` {r }
cyphr::data_admin_init(data_dir, path_user = path_key_alice)
```

The data key is very important.  If it is deleted, then the data
cannot be decrypted.  So do not delete the directory
`data_dir/.cyphr`!  Ideally add it to your version control
system so that it cannot be lost.  Of course, if you're working in
a group, there are multiple copies of the data key (each encrypted
with a different person's personal key) which reduces the chance of
total loss.

This command can be run multiple times safely; if it detects it has
been rerun and the data key will not be regenerated.
``` {r }
cyphr::data_admin_init(data_dir, path_user = path_key_alice)
```

**Third**, you can add encrypted data to the directory (or to
anywhere really).  When run, `cyphr::config_data` will verify
that it can actually decrypt things.
``` {r }
key <- cyphr::data_key(data_dir, path_user = path_key_alice)
```

This object can be used with all the `cyphr` functions (see the
"cyphr" vignette; `vignette("cyphr")`)
``` {r }
filename <- file.path(data_dir, "iris.rds")
cyphr::encrypt(saveRDS(iris, filename), key)
dir(data_dir)
```

The file is encrypted and so cannot be read with `readRDS`:
``` {r error = TRUE}
readRDS(filename)
```

But we can decrypt and read it:
``` {r }
head(cyphr::decrypt(readRDS(filename), key))
```

**Fourth**, have someone else join in.  Recall that to simulate
another person here, I'm going to pass an argument `path_user =
path_key_bob` though to the functions.  This contains the path to
"Bob"'s ssh keypair.  If run on an actually different computer this
would not be needed; this is just to simulate two users in a single
session for this vignette (see minimal example below where this is
simulated).  Again, typically this user would also not use the
`cyphr::ssh_keygen` function but use the `ssh-keygen` command from
their shell.

We're going to assume that the user can read and write to the data.
This is the case for my use case where the data are stored on
dropbox and will be the case with GitHub based distribution, though
there would be a pull request step in here.

This user cannot read the data, though trying to will print a
message explaining how you might request access:
``` {r error = TRUE}
key_bob <- cyphr::data_key(data_dir, path_user = path_key_bob)
```

But `bob` is your collaborator and needs access!  What they need
to do is run:
``` {r }
cyphr::data_request_access(data_dir, path_user = path_key_bob)
```

(again, ordinarily you would not need the `bob` bit here)

The user should the send an email to someone with access and quote
the hash in the message above.

**Fifth**, back on the first computer we can authorise the second
user.  First, see who has requested access:
``` {r }
req <- cyphr::data_admin_list_requests(data_dir)
req
```

We can see the same hash here as above (``r names(req)[[1]]``)

...and then grant access to them with the
`cyphr::data_admin_authorise` function.
``` {r }
cyphr::data_admin_authorise(data_dir, yes = TRUE, path_user = path_key_alice)
```

If you do not specify `yes = TRUE` will prompt for confirmation at
each key added.

This has cleared the request queue:
``` {r }
cyphr::data_admin_list_requests(data_dir)
```

and added it to our set of keys:
``` {r }
cyphr::data_admin_list_keys(data_dir)
```

**Finally**, as soon as the authorisation has happened, the user
can encrypt and decrypt files:
``` {r }
key_bob <- cyphr::data_key(data_dir, path_user = path_key_bob)
head(cyphr::decrypt(readRDS(filename), key_bob))
```

## Minimal example

As above, but with less discussion:

``` {r echo = FALSE, results = "hide"}
unlink(data_dir, recursive = TRUE)
dir.create(data_dir)
```

Setup, on Alice's computer:
``` {r }
cyphr::data_admin_init(data_dir, path_user = path_key_alice)
```

Get the data key key:
``` {r }
key <- cyphr::data_key(data_dir, path_user = path_key_alice)
```

Encrypt a file:
``` {r }
cyphr::encrypt(saveRDS(iris, filename), key)
```

Request access, on Bob's computer:
``` {r }
hash <- cyphr::data_request_access(data_dir, path_user = path_key_bob)
```

Alice authorises this request::
``` {r }
cyphr::data_admin_authorise(data_dir, yes = TRUE, path_user = path_key_alice)
```

Bob can get the data key:
``` {r }
key <- cyphr::data_key(data_dir, path_user = path_key_bob)
```

Bob can read the secret data:
``` {r }
head(cyphr::decrypt(readRDS(filename), key))
```

## Details & disclosure

Encryption does not work through security through obscurity; it
works because we can rely on the underlying maths enough to be open
about how things are stored and where.

Most encryption libraries require some degree of security in
the underlying software.  Because of the way R works this is very
difficult to guarantee; it is trivial to rewrite code in running
packages to skip past verification checks.  So this package is
_not_ designed to (or able to) avoid exploits in your running code;
an attacker could intercept your private keys, the private key to
the data, or skip the verification checks that are used to make
sure that the keys you load are what they say they are.  However,
the _data_ are safe; only people who have keys to the data will be
able to read it.

`cyphr` uses two different encryption algorithms; it uses RSA
encryption via the `openssl` package for user keys, because there
is a common file format for these keys so it makes user
configuration easier.  It uses the modern sodium package (and
through that the libsodium library) for data encryption because it
is very fast and simple to work with.  This does leave two possible
points of weakness as a vulnerability in either of these libraries
could lead to an exploit that could allow decryption of your data.

Each user has a public/private key pair.  Typically this is in
`~/.ssh/id_rsa.pub` and `~/.ssh/id_rsa`, and if found these will be
used.  Alternatively the location of the keypair can be stored
elsewhere and pointed at with the `USER_KEY` or `USER_PUBKEY`
environment variables.  The key may be password protected (and this
is recommended!) and the password will be requested without ever
echoing it to the terminal.

The data directory has a hidden directory `.cyphr` in it.
``` {r }
dir(data_dir, all.files = TRUE, no.. = TRUE)
```

This does not actually need to be stored with the data but it makes
sense to (there are workflows where data is stored remotely where
storing this directory might make sense).  The "keys" directory
contains a number of files; one for each person who has access to
the data.
``` {r }
dir(file.path(data_dir, ".cyphr", "keys"))
names(cyphr::data_admin_list_keys(data_dir))
```

(the file `test` is a small file encrypted with the data key used
to verify everything is working OK).

Each file is stored in RDS format and is a list with elements:

* user: the reported user name of the person who created request for data
* host: the reported computer name
* date: the time the request was generated
* pub: the RSA public key of the user
* key: the data key, encrypted with the user key.  Without the
  private key, this cannot be used.  With the user's private key
  this can be used to generate the symmetric key to the data.

``` {r }
h <- names(cyphr::data_admin_list_keys(data_dir))[[1]]
readRDS(file.path(data_dir, ".cyphr", "keys", h))
```

You can see that the hash of the public key is the same as name of
the stored file here (which is used to prevent collisions when
multiple people request access at the same time).
``` {r }
h
```

When a request is posted it is an RDS file with all of the above
except for the `key` element, which is added during authorisation.

(Note that the verification relies on the package code not being
attacked, and given R's highly dynamic nature an attacker could
easily swap out the definition for the verification function with
something that always returns `TRUE`.)

When an authorised user creates the `data_key` object (which
allows decryption of the data) `secret` will:

* read their private user key (probably from `~/.ssh/id_rsa`)
* read the encrypted data key from the data directory (the `$key`
  element from the list above).
* decrypt this data key using their user key to yield the the data
  symmetric key.

## Limitations

In the Dropbox scenario, non-password protected keys will afford
only limited protection.  This is because even though the keys and
data are stored separately on Dropbox, they will be in the same
place on a local computer; if that computer is lost then the only
thing preventing an attacker recovering the data is security
through obscurity (the data would appear to be random junk but
they will be able to run your analysis scripts as easily as you
can).  Password protected keys will improve this situation
considerably as without a password the data cannot be recovered.

The data is not encrypted during a running R session.  R allows
arbitrary modification of code at runtime so this package provides
no security from the point where the data can be decrypted.  If
your computer was compromised then stealing the data while you are
running R should be assumed to be straightforward.

``` {r echo = FALSE, results = "hide"}
unlink(c(data_dir, path_key_alice, path_key_bob), recursive = TRUE)
```
---
title: "Introduction"
author: "Rich FitzJohn"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

This package tries to smooth over some of the differences in
encryption approaches (symmetric vs. asymmetric, sodium vs. openssl)
to provide a simple interface for users who just want to encrypt or
decrypt things.

The scope of the package is to protect data that has been saved to
disk.  It is not designed to stop an attacker targeting the R
process itself to determine the contents of sensitive data.  The
package does try to prevent you accidentally saving to disk the
contents of sensitive information, including the keys that could
decrypt such information.

This vignette works through the basic functionality of the package.
It does not offer much in the way of an introduction to encryption
itself; for that see the excellent vignettes in the `openssl` and
`sodium` packages (see `vignette("crypto101")` and
`vignette("bignum")` for information about how encryption works).
This package is a wrapper around those packages in order to make
them more accessible.

# Keys and the like

To encrypt anything we need a key.  There are two sorts of key
"types" we will concern ourselves with here "symmetric" and
"asymmetric".

* "symmetric" keys are used for storing secrets that multiple
  people need to access.  Everyone has the same key (which is just
  a bunch of bytes) and with that we can either encrypt data or
  decrypt it.

* a "key pair" is a public and a private key; this is used in
  communication.  You hold a private key that nobody else ever sees
  and a public key that you can copy around all over the show.
  These can be used for a couple of different patterns of
  communication (see below).

We support symmetric keys and asymmetric key pairs from the
`openssl` and `sodium` packages (which wrap around
industry-standard cryptographic libraries) - this vignette will
show how to create and load keys of different types as they're
used.

The `openssl` keys have the advantage of a standard key format, and
that many people (especially on Linux and macOS) have a keypair
already (see below if you're not sure if you do).  The `sodium`
keys have the advantage of being a new library, starting from a
clean slate rather than carrying with it accumulated ideas from the
last 20 years of development.

The idea in `cyphr` is that we can abstract away some differences
in the types of keys and the functions that go with them to create
a standardised interface to encrypting and decrypting strings, R
objects, files and raw vectors.  With that, we can then create
wrappers around functions that create files and simplify the
process of adding encryption into a data workflow.

Below, I'll describe the sorts of keys that `cyphr` supports and in
the sections following describe how these can be used to actually
do some encryption.

## Symmetric encryption

![Illustration of Symmetric Encryption](symmetric.png)

This is the simplest form of encryption because everyone has the
same key (like a key to your house or a single password).  This
raises issues (like how do you *store* the key without other people
reading it) but we can deal with that below.

### `openssl`

To generate a key with `openssl`, you can use:
``` {r }
k <- openssl::aes_keygen()
```

which generates a raw vector
``` {r }
k
```

(this prints nicely but it really is stored as a 16 byte raw
vector).

The encryption functions that this key supports are
`openssl::aes_cbc_encrypt`, `openssl::aes_ctr_encrypt` and
`openssl::aes_gcm_encrypt` (along with the corresponding decryption
functions).  The `cyphr` package tries to abstract this away by
using a wrapper `cyphr::key_openssl
``` {r }
key <- cyphr::key_openssl(k)
key
```

With this key, one can encrypt a string with `cyphr::encrypt_string`:
``` {r }
secret <- cyphr::encrypt_string("my secret string", key)
```

and decrypt it again with `cyphr::decrypt_string`:
``` {r }
cyphr::decrypt_string(secret, key)
```

See below for more functions that use these key objects.

### `sodium`

The interface is almost identical using sodium symmetric keys.  To
generate a symmetric key with libsodium you would use
`sodium::keygen`
``` {r }
k <- sodium::keygen()
```

This is really just a raw vector of length 32, without even any
class attribute!

The encryption functions that this key supports are
`sodium::data_encrypt` and `sodium::data_decrypt`.  To create a key
for use with `cyphr` that knows this, use:
``` {r }
key <- cyphr::key_sodium(k)
key
```

This key can then be used with the high-level cyphr encryption
functions described below.

## Asymmetric encryption ("key pairs")

![Illustration of Asymmetric Encryption](asymmetric.png)

With asymmetric encryption everybody has two keys that differ from
everyone else's key.  One key is public and can be shared freely
with anyone you would like to communicate with and the other is
private and must never be disclosed.

In the `sodium` package there is a vignette
(`vignette("crypto101")`) that gives a gentle introduction to how
this all works.  In practice, you end up creating a pair of keys
for yourself.  Then to encrypt or decrypt something you encrypt
messages with the recipient's *public key* and they (and only they)
can decrypt it with their *private key*.

One use for asymmetric encryption is to encrypt a shared secret
(such as a symmetric key) - with this you can then safely store or
communicate a symmetric key without disclosing it.

### `openssl`

Let's suppose that we have two parties "Alice" and "Bob" who want
to talk with one another.  For demonstration purposes we need to
generate SSH keys (with no password) in temporary directories (to
comply with CRAN policies).  In a real situation these would be on
different machines (Alice has no access to Bob's key!) and these
keys would be password protected.
``` {r }
path_key_alice <- cyphr::ssh_keygen(password = FALSE)
path_key_bob <- cyphr::ssh_keygen(password = FALSE)
```

Note that each directory contains a public key (`id_rsa.pub`) and a
private key (`id_rsa`).
``` {r }
dir(path_key_alice)
dir(path_key_bob)
```

Below, the full path to the key (e.g., `.../id_rsa`) could be
used in place of the directory name if you prefer.

If Alice wants to send a message to Bob she needs to use her
private key and his public key
``` {r }
pair_a <- cyphr::keypair_openssl(path_key_bob, path_key_alice)
pair_a
```

with this pair she can write a message to "bob":
``` {r }
secret <- cyphr::encrypt_string("secret message", pair_a)
```

The secret is now just a big pile of bytes
``` {r }
secret
```

Note that unlike symmetric encryption above, Alice cannot decrypt
her own message:
``` {r error = TRUE}
cyphr::decrypt_string(secret, pair_a)
```

For Bob to read the message, he uses his private key and Alice's
public key (which she has transmitted to him previously).
``` {r }
pair_b <- cyphr::keypair_openssl(path_key_alice, path_key_bob)
```

With this keypair, Bob can decrypt Alice's message
``` {r }
cyphr::decrypt_string(secret, pair_b)
```

And send one back of his own:
``` {r }
secret2 <- cyphr::encrypt_string("another message", pair_b)
secret2
```

which she can decrypt
``` {r }
cyphr::decrypt_string(secret2, pair_a)
```

Chances are, you have an openssl keypair in your `.ssh/` directory.
If so, you would pass `NULL` as the path for the private (or less
usefully, the public) key pair part.  So to send a message to Bob,
we'd include the path to Bob's public key.
```r
pair_us <- cyphr::keypair_openssl(path_key_bob, NULL)
```

This all skips over how Alice and Bob will exchange this secret
information.  Because the secret is bytes, it's a bit odd to work
with.  Alice could save the secret to disk with
``` {r }
secret <- cyphr::encrypt_string("secret message", pair_a)
path_for_bob <- file.path(tempdir(), "for_bob_only")
writeBin(secret, path_for_bob)
```

And then send Bob the file `for_bob_only` (over email or any other
insecure medium).

and bob could read the secret in with:
``` {r }
secret <- readBin(path_for_bob, raw(), file.size(path_for_bob))
cyphr::decrypt_string(secret, pair_b)
```

As an alternative, you can "base64 encode" the bytes into something
that you can just email around:
``` {r }
secret_base64 <- openssl::base64_encode(secret)
secret_base64
```

This can be converted back with `openssl::base64_decode`:
``` {r }
identical(openssl::base64_decode(secret_base64), secret)
```

Or, less compactly but also suitable for email, you might just
convert the bytes into their hex representation:
``` {r }
secret_hex <- sodium::bin2hex(secret)
secret_hex
```

and the reverse with `sodium::hex2bin`:
``` {r }
identical(sodium::hex2bin(secret_hex), secret)
```

(this is somewhat less space efficient than base64 encoding.

As a final option, you can just save the secret with `saveRDS` and
read it in with `readRDS` like any other option.  This will be the
best route if the secret is saved into a more complicated R object
(e.g., a list or `data.frame`).

See the other cyphr vignette (`vignette("data", package =
"cyphr")`) for a suggested workflow for exchanging secrets within a
team, and the wrapper functions below for more convenient ways of
working with encrypted data.

**Do you already have an ssh keypair?** To find out, run

```r
cyphr::keypair_openssl(NULL, NULL)
```

One of three things will happen:

1. you will be prompted for your password to decrypt your private
key, and then after entering it an object `<cyphr_keypair:
openssl>` will be returned - you're good to go!

2. you were _not_ prompted for your password, but got a
`<cyphr_keypair: openssl>` object.  You should consider whether
this is appropriate and consider generating a new keypair with the
private key encrypted.  If you don't then anyone who can read your
private key can decrypt any message intended for you.

3. you get an error like `Did not find default ssh public key at
~/.ssh/id_rsa.pub`.  You need to create a keypair.

To create a keypair, you can use the `cyphr::ssh_keygen()` function as

```r
cyphr::ssh_keygen("~/.ssh")
```

This will create the keypair as `~/.ssh/id_rsa` and
`~/.ssh/id_rsa.pub`, which is where `cyphr` will look for your keys
by default.  See `?ssh_keygen` for more information.  (On Linux and
macOS you might use the `ssh-keygen` command line utility.  On
windows, PuTTY` has a utility for creating keys.)

### `sodium`

With `sodium`, things are largely the same with the exception that
there is no standard format for saving sodium keys.  The bits below
use an in-memory key (which is just a collection of bytes) but
these can also be filenames, each of which contains the contents of
the key written out with `writeBin`.

First, generate keys for Alice:
``` {r }
key_a <- sodium::keygen()
pub_a <- sodium::pubkey(key_a)
```

the public key is derived from the private key, and Alice can share
that with Bob.  We next generate Bob's keys
``` {r }
key_b <- sodium::keygen()
pub_b <- sodium::pubkey(key_b)
```

Bob would now share is public key with Alice.

If Alice wants to send a message to Bob she again uses her private
key and Bob's public key:
``` {r }
pair_a <- cyphr::keypair_sodium(pub_b, key_a)
```

As above, she can now send a message:
``` {r }
secret <- cyphr::encrypt_string("secret message", pair_a)
secret
```

Note how this line is identical to the one in the `openssl` section.

To decrypt this message, Bob would use Alice's public key and his
private key:
``` {r }
pair_b <- cyphr::keypair_sodium(pub_a, key_b)
cyphr::decrypt_string(secret, pair_b)
```

# Encrypting things

Above, we used `cyphr::encrypt_string` and `cyphr::decrypt_string`
to encrypt and decrypt a string.  There are several such functions
in the package that encrypt and decrypt

* R objects `encrypt_object` / `decrypt_object` (using serialization
  and deserialization)
* strings: `encrypt_string` / `decrypt_string`
* raw vectors: `encrypt_data` / `decrypt_data`
* files: `encrypt_file` / `decrypt_file`

For this section we will just use a sodium symmetric encryption
key
``` {r }
key <- cyphr::key_sodium(sodium::keygen())
```

For the examples below, in the case of asymmetric encryption (using
either `cyphr::keypair_openssl` or `cyphr::keypair_sodium`) the
sender would use their private key and the recipient's public key
and the recipient would use the complementary key pair.

## Objects

Here's an object to encrypt:
``` {r }
obj <- list(x = 1:10, y = "secret")
```

This creates a bunch of raw bytes corresponding to the data (it's
not really possible to print this as anything nicer than bytes).
``` {r }
secret <- cyphr::encrypt_object(obj, key)
secret
```

The data can be decrypted with the `decrypt_object` function:
``` {r }
cyphr::decrypt_object(secret, key)
```

Optionally, this process can go via a file, using a third argument
to the functions (note that temporary files are used here for
compliance with CRAN policies - any path may be used in practice).
``` {r }
path_secret <- file.path(tempdir(), "secret.rds")
cyphr::encrypt_object(obj, key, path_secret)
```

There is now a file called `secret.rds` in the temporary directory:
``` {r }
file.exists(path_secret)
```

though it is not actually an rds file:
``` {r error = TRUE}
readRDS(path_secret)
```

When passed a filename (as opposed to a raw vector),
`cyphr::decrypt_object` will read the object in before decrypting
it
``` {r }
cyphr::decrypt_object(path_secret, key)
```

## Strings

For the case of strings we can do this in a slightly more
lightweight way (the above function routes through `serialize` /
`deserialize` which can be slow and will create larger objects than
using `charToRaw` / `rawToChar`)
``` {r }
secret <- cyphr::encrypt_string("secret", key)
secret
```

and decrypt:
``` {r }
cyphr::decrypt_string(secret, key)
```

## Plain raw data

If these are not enough for you, you can work directly with raw
objects (bunches of bytes) by using `encrypt_data`:
``` {r }
dat <- sodium::random(100)
dat # some random bytes

secret <- cyphr::encrypt_data(dat, key)
secret
```

Decrypted data is the same as a the original data
``` {r }
identical(cyphr::decrypt_data(secret, key), dat)
```

## Files

Suppose we have written a file that we want to encrypt to send to
someone (in a temporary directory for compliance with CRAN
policies)
``` {r }
path_data_csv <- file.path(tempdir(), "iris.csv")
write.csv(iris, path_data_csv, row.names = FALSE)
```

You can encrypt that file with
``` {r }
path_data_enc <- file.path(tempdir(), "iris.csv.enc")
cyphr::encrypt_file(path_data_csv, key, path_data_enc)
```

This encrypted file can then be decrypted with
``` {r }
path_data_decrypted <- file.path(tempdir(), "idis2.csv")
cyphr::decrypt_file(path_data_enc, key, path_data_decrypted)
```

Which is identical to the original:
``` {r }
tools::md5sum(c(path_data_csv, path_data_decrypted))
```

# An even higher level interface for files

This is the most user-friendly way of using the package when the
aim is to encrypt and decrypt files.  The package provides a pair
of functions `cyphr::encrypt` and `cyphr::decrypt` that wrap file
writing and file reading functions.  In general you would use
`encrypt` when writing a file and `decrypt` when reading one.
They're designed to be used like so:

Suppose you have a super-secret object that you want to share privately
``` {r }
key <- cyphr::key_sodium(sodium::keygen())
x <- list(a = 1:10, b = "don't tell anyone else")
```

If you save `x` to disk with `saveRDS` it will be readable by
everyone until it is deleted.  But if you encrypted the file that
`saveRDS` produced it would be protected and only people with the
key can read it:
``` {r }
path_object <- file.path(tempdir(), "secret.rds")
cyphr::encrypt(saveRDS(x, path_object), key)
```

(see below for some more details on how this works).

This file cannot be read with `readRDS`:

``` {r error = TRUE}
readRDS(path_object)
```

but if we wrap the call with `decrypt` and pass in the config
object it can be decrypted and read:
``` {r }
cyphr::decrypt(readRDS(path_object), key)
```

What happens in the call above is `cyphr` uses "non standard
evaluation" to rewrite the call above so that it becomes
(approximately)

1. use `cyphr::decrypt_file` to decrypt "secret.rds" as a temporary file
2. call `readRDS` on that temporary file
3. delete the temporary file (even if there is an error in the above calls)

This non-standard evaluation breaks referential integrity (so may
not be suitable for programming).  You can always do this manually
with `encrypt_file` / `decrypt_file` so long as you make sure to
clean up after yourself.

The `encrypt` function inspects the call in the first argument
passed to it and works out for the function provided (`saveRDS`)
which argument corresponds to the filename (here `"secret.rds"`).
It then rewrites the call to write out to a temporary file (using
`tempfile()`).  Then it calls `encrypt_file` (see below) on this
temporary file to create the file asked for (`"secret.rds"`).  Then
it deletes the temporary file, though this will also happen in case
of an error in any of the above.

The `decrypt` function works similarly.  It inspects the call and
detects that the first argument represents the filename.  It
decrypts that file to create a temporary file, and then runs
`readRDS` on that file.  Again it will delete the temporary file on
exit.

The functions supported via this interface are:

* `readLines` / `writeLines`
* `readRDS` / `writeRDS`
* `read` / `save`
* `read.table` / `write.table`
* `read.csv` / `read.csv2` / `write.csv`
* `read.delim` / `read.delim2`

But new functions can be added with the `rewrite_register`
function.  For example, to support the excellent
[rio](https://cran.r-project.org/package=rio) package, whose
`import` and `export` functions take the filename `file` you could
use:

```r
cyphr::rewrite_register("rio", "import", "file")
cyphr::rewrite_register("rio", "export", "file")
```

now you can read and write tabular data into and out of a great
many different file formats with encryption with calls like

```r
cyphr::encrypt(rio::export(mtcars, "file.json"), key)
cyphr::decrypt(rio::import("file.json"), key)
```

The functions above use [non standard evaluation](http://adv-r.had.co.nz/Computing-on-the-language.html)
and so may not be suitable for programming or use in packages.  An
"escape hatch" is provided via `encrypt_` and `decrypt_` where the
first argument is a quoted expression.
``` {r }
cyphr::encrypt_(quote(saveRDS(x, path_object)), key)
cyphr::decrypt_(quote(readRDS(path_object)), key)
```

# Session keys

When using `key_openssl`, `keypair_openssl`, `key_sodium`, or
`keypair_sodium` we generate something that can decrypt data.  The
objects that are returned by these functions can encrypt and
decrypt data and so it is reasonable to be concerned that if these
objects were themselves saved to disk your data would be
compromised.

To avoid this, `cyphr` does not store private or symmetric keys
directly in these objects but instead encrypts the sensitive keys
with a `cyphr`-specific session key that is regenerated each time
the package is loaded.  This means that the objects are practically
only useful within one session, and if saved with `save.image`
(perhaps automatically at the end of a session) the keys cannot be
used to decrypt data.

To manually invalidate all keys you can use the
`cyphr::session_key_refresh` function.  For example, here is a
symmetric key:
``` {r }
key <- cyphr::key_sodium(sodium::keygen())
```

which we can use to encrypt a secret string
``` {r }
secret <- cyphr::encrypt_string("my secret", key)
```

and decrypt it:
``` {r }
cyphr::decrypt_string(secret, key)
```

If we refresh the session key we invalidate the `key` object
``` {r }
cyphr::session_key_refresh()
```

and after this point the key cannot be used any further
``` {r error = TRUE}
cyphr::decrypt_string(secret, key)
```

This approach works because the package holds the session key
within its environment (in `cyphr:::session$key`) which R will not
serialize.  As noted above - this approach does not prevent an
attacker with the ability to snoop on your R session from
discovering your private keys or sensitive data but it does prevent
accidentally saving keys in a way that would be useful for an
attacker to use in a subsequent session.

``` {r include = FALSE}
unlink(c(path_secret, path_object, path_data_csv, path_data_enc,
         path_data_decrypted, path_for_bob,
         path_key_alice, path_key_bob),
       recursive = TRUE)
```

# Further reading

* The wikipedia page on Public Key cryptography has some nice
  diagrams that explain how key and data interact
  https://en.wikipedia.org/wiki/Public-key_cryptography
* The vignettes in the `openssl` (`vignette(package = "openssl")`)
  and `sodium` (`vignette(package = "openssl")`) packages have
  explanations of how the tools used in `cyphr` work and interface
  with R.

Confused?  Need help?  Found a bug?

* Post an issue on the [`cyphr` issue
  tracker](https://github.com/ropensci/cyphr/issues)
* Start a discussion on the [rOpenSci discussion
  forum](https://discuss.ropensci.org/)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/encrypt_wrapper.R
\name{encrypt}
\alias{encrypt}
\alias{decrypt}
\alias{encrypt_}
\alias{decrypt_}
\title{Easy encryption and decryption}
\usage{
encrypt(expr, key, file_arg = NULL, envir = parent.frame())

decrypt(expr, key, file_arg = NULL, envir = parent.frame())

encrypt_(expr, key, file_arg = NULL, envir = parent.frame())

decrypt_(expr, key, file_arg = NULL, envir = parent.frame())
}
\arguments{
\item{expr}{A single expression representing a function call that
would be called for the side effect of creating or reading a
file.}

\item{key}{A \code{cyphr_key} object describing the
encryption approach to use.}

\item{file_arg}{Optional hint indicating which argument to
\code{expr} is the filename.  This is done automatically for
some built-in functions.}

\item{envir}{Environment in which \code{expr} is to be evaluated.}
}
\description{
Wrapper functions for encryption.  These functions wrap
expressions that produce or consume a file and arrange to encrypt
(for producing functions) or decrypt (for consuming functions).
The forms with a trailing underscore (\code{encrypt_},
\code{decrypt_}) do not use any non-standard evaluation and may be
more useful for programming.
}
\details{
These functions will not work for all functions.  For example
\code{pdf}/\code{dev.off} will create a file but we can't wrap
those up (yet!).  Functions that \emph{modify} a file (e.g.,
appending) also will not work and may cause data loss.
}
\examples{
# To do anything we first need a key:
key <- cyphr::key_sodium(sodium::keygen())

# Encrypted write.csv - note how any number of arguments to
# write.csv will be passed along
path <- tempfile(fileext = ".csv")
cyphr::encrypt(write.csv(iris, path, row.names = FALSE), key)

# The new file now exists, but you would not be able to read it
# with read.csv because it is now binary data.
file.exists(path)

# Wrap the read.csv call with cyphr::decrypt()
dat <- cyphr::decrypt(read.csv(path, stringsAsFactors = FALSE), key)
head(dat)

file.remove(path)

# If you have a function that is not supported you can specify the
# filename argument directly.  For example, with "write.dcf" the
# filename argument is called "file"; we can pass that along
path <- tempfile()
cyphr::encrypt(write.dcf(list(a = 1), path), key, file_arg = "file")

# Similarly for decryption:
cyphr::decrypt(read.dcf(path), key, file_arg = "file")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/encrypt_data.R
\name{encrypt_data}
\alias{encrypt_data}
\alias{encrypt_object}
\alias{encrypt_string}
\alias{encrypt_file}
\alias{decrypt_data}
\alias{decrypt_object}
\alias{decrypt_string}
\alias{decrypt_file}
\title{Encrypt and decrypt data and other things}
\usage{
encrypt_data(data, key, dest = NULL)

encrypt_object(object, key, dest = NULL, rds_version = NULL)

encrypt_string(string, key, dest = NULL)

encrypt_file(path, key, dest = NULL)

decrypt_data(data, key, dest = NULL)

decrypt_object(data, key)

decrypt_string(data, key)

decrypt_file(path, key, dest = NULL)
}
\arguments{
\item{data}{(for \code{encrypt_data}, \code{decrypt_data},
\code{decrypt_object}, \code{decrypt_string}) a raw vector with
the data to be encrypted or decrypted.  For the decryption
functions this must be data derived by encrypting something or
you will get an error.}

\item{key}{A \code{cyphr_key} object describing the encryption approach
to use.}

\item{dest}{The destination filename for the encrypted or
decrypted data, or \code{NULL} to return a raw vector.  This is
not used by \code{decrypt_object} or \code{decrypt_string} which
always return an object or string.}

\item{object}{(for \code{encrypt_object}) an arbitrary R object to
encrypt.  It will be serialised to raw first (see
\link{serialize}).}

\item{rds_version}{RDS serialisation version to use (see
\link{serialize}.  The default in R version 3.3 and below is version
2 - in the R 3.4 series version 3 was introduced and is becoming
the default.  Version 3 format serialisation is not understood
by older versions so if you need to exchange data with older R
versions, you will need to use \code{rds_version = 2}.  The default
argument here (\code{NULL}) will ensure the same serialisation is
used as R would use by default.}

\item{string}{(for \code{encrypt_string}) a scalar character
vector to encrypt.  It will be converted to raw first with
\link{charToRaw}.}

\item{path}{(for \code{encrypt_file}) the name of a file to
encrypt.  It will first be read into R as binary (see
\link{readBin}).}
}
\description{
Encrypt and decrypt raw data, objects, strings and files.  The
core functions here are \code{encrypt_data} and
\code{decrypt_data} which take raw data and decrypt it, writing
either to file or returning a raw vector.  The other functions
encrypt and decrypt arbitrary R objects (\code{encrypt_object},
\code{decrypt_object}), strings (\code{encrypt_string},
\code{decrypt_string}) and files (\code{encrypt_file},
\code{decrypt_file}).
}
\examples{
key <- key_sodium(sodium::keygen())
# Some super secret data we want to encrypt:
x <- runif(10)
# Convert the data into a raw vector:
data <- serialize(x, NULL)
data
# Encrypt the data; without the key above we will never be able to
# decrypt this.
data_enc <- encrypt_data(data, key)
data_enc
# Our random numbers:
unserialize(decrypt_data(data_enc, key))
# Same as the never-encrypted version:
x

# This can be achieved more easily using `encrypt_object`:
data_enc <- encrypt_object(x, key)
identical(decrypt_object(data_enc, key), x)

# Encrypt strings easily:
str_enc <- encrypt_string("secret message", key)
str_enc
decrypt_string(str_enc, key)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{data_request_access}
\alias{data_request_access}
\alias{data_key}
\title{User commands}
\usage{
data_request_access(path_data = NULL, path_user = NULL, quiet = FALSE)

data_key(
  path_data = NULL,
  path_user = NULL,
  test = TRUE,
  quiet = FALSE,
  cache = TRUE
)
}
\arguments{
\item{path_data}{Path to the data.  If not given, then we look
recursively down below the working directory for a ".cyphr"
directory, and use that as the data directory.}

\item{path_user}{Path to the directory with your user key.
Usually this can be omitted.  This argument is passed in as both
\code{pub} and \code{key} to \code{\link[=keypair_openssl]{keypair_openssl()}}.
Briefly, if this argument is not given we look at the
environment variables \code{USER_PUBKEY} and \code{USER_KEY} -
if set then these must refer to path of your public and private
keys.  If these environment variables are not set then we fall
back on \verb{~/.ssh/id_rsa.pub} and \verb{~/.ssh/id_rsa},
which should work in most environments.  Alternatively, provide
a path to a directory where the file \code{id_rsa.pub} and
\code{id_rsa} can be found.}

\item{quiet}{Suppress printing of informative messages.}

\item{test}{Test that the encryption is working?  (Recommended)}

\item{cache}{Cache the key within the session.  This will be
useful if you are using ssh keys that have passwords, as if the
key is found within the cache, then you will not have to
re-enter your password.  Using \code{cache = FALSE} neither
looks for the key in the cache, nor saves it.}
}
\description{
User commands
}
\examples{

# The workflow here does not really lend itself to an example,
# please see the vignette.

# Suppose that Alice has created a data directory:
path_alice <- tempfile()
cyphr::ssh_keygen(path_alice, password = FALSE)
path_data <- tempfile()
dir.create(path_data, FALSE, TRUE)
cyphr::data_admin_init(path_data, path_user = path_alice)

# If Bob can also write to the data directory (e.g., it is a
# shared git repo, on a shared drive, etc), then he can request
# access
path_bob <- tempfile()
cyphr::ssh_keygen(path_bob, password = FALSE)
hash <- cyphr::data_request_access(path_data, path_user = path_bob)

# Alice can authorise Bob
cyphr::data_admin_authorise(path_data, path_user = path_alice, yes = TRUE)

# After which Bob can get the data key
cyphr::data_key(path_data, path_user = path_bob)

# See the vignette for more details.  This is not the best medium
# to explore this.

# Cleanup
unlink(path_alice, recursive = TRUE)
unlink(path_bob, recursive = TRUE)
unlink(path_data, recursive = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sodium.R
\name{key_sodium}
\alias{key_sodium}
\title{Symmetric encryption with sodium}
\usage{
key_sodium(key)
}
\arguments{
\item{key}{A sodium key (i.e., generated with \code{\link[sodium:keygen]{sodium::keygen()}}}
}
\description{
Wrap a sodium symmetric key.  This can be used with the functions
\code{\link[=encrypt_data]{encrypt_data()}} and \code{\link[=decrypt_data]{decrypt_data()}}, along
with the higher level wrappers \code{\link[=encrypt]{encrypt()}} and
\code{\link[=decrypt]{decrypt()}}.  With a symmetric key, everybody uses the
same key for encryption and decryption.
}
\examples{
# Create a new key
key <- cyphr::key_sodium(sodium::keygen())
key

# With this key encrypt a string
secret <- cyphr::encrypt_string("my secret string", key)
# And decrypt it again:
cyphr::decrypt_string(secret, key)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{data_admin_init}
\alias{data_admin_init}
\alias{data_admin_authorise}
\alias{data_admin_list_requests}
\alias{data_admin_list_keys}
\title{Encrypted data administration}
\usage{
data_admin_init(path_data, path_user = NULL, quiet = FALSE)

data_admin_authorise(
  path_data = NULL,
  hash = NULL,
  path_user = NULL,
  yes = FALSE,
  quiet = FALSE
)

data_admin_list_requests(path_data = NULL)

data_admin_list_keys(path_data = NULL)
}
\arguments{
\item{path_data}{Path to the data set.  We will store a bunch of
things in a hidden directory within this path.  By default in
most functions we will search down the tree until we find the
.cyphr directory}

\item{path_user}{Path to the directory with your ssh key.
Usually this can be omitted.}

\item{quiet}{Suppress printing of informative messages.}

\item{hash}{A vector of hashes to add.  If provided, each hash can
be the binary or string representation of the hash to add.  Or
omit to add each request.}

\item{yes}{Skip the confirmation prompt?  If any request is
declined then the function will throw an error on exit.}
}
\description{
Encrypted data administration; functions for setting up, adding
users, etc.
}
\details{
\code{data_admin_init} initialises the system; it will create a
data key if it does not exist and authorise you.  If it already
exists and you do not have access it will throw an error.

\code{data_admin_authorise} authorises a key by creating a key to
the data that the user can use in conjunction with their personal
key.

\code{data_admin_list_requests} lists current requests.

\code{data_admin_list_keys} lists known keys that can access the
data.  Note that this is \emph{not secure}; keys not listed here
may still be able to access the data (if a key was authorised and
moved elsewhere for example).  Conversely, if the user has deleted
or changed their key they will not be able to access the data
despite the key being listed here.
}
\examples{

# The workflow here does not really lend itself to an example,
# please see the vignette instead.

# First we need a set of user ssh keys.  In a non example
# environment your personal ssh keys will probably work well, but
# hopefully they are password protected so cannot be used in
# examples.  The password = FALSE argument is only for testing,
# and should not be used for data that you care about.
path_ssh_key <- tempfile()
cyphr::ssh_keygen(path_ssh_key, password = FALSE)

# Initialise the data directory, using this key path.  Ordinarily
# the path_user argument would not be needed because we would be
# using your user ssh keys:
path_data <- tempfile()
dir.create(path_data, FALSE, TRUE)
cyphr::data_admin_init(path_data, path_user = path_ssh_key)

# Now you can get the data key
key <- cyphr::data_key(path_data, path_user = path_ssh_key)

# And encrypt things with it
cyphr::encrypt_string("hello", key)

# See the vignette for more details.  This is not the best medium
# to explore this.

# Cleanup
unlink(path_ssh_key, recursive = TRUE)
unlink(path_data, recursive = TRUE)
}
\seealso{
\code{\link[=data_request_access]{data_request_access()}} for requesting access
to the data, and and \code{data_key} for using the data
itself.  But for a much more thorough overview, see the vignette
(\code{vignette("data", package = "cyphr")}).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/session.R
\name{session_key_refresh}
\alias{session_key_refresh}
\title{Refresh the session key}
\usage{
session_key_refresh()
}
\description{
Refresh the session key, invalidating all keys created by
\code{\link[=key_openssl]{key_openssl()}}, \code{\link[=keypair_openssl]{keypair_openssl()}},
\code{\link[=key_sodium]{key_sodium()}} and \code{\link[=keypair_sodium]{keypair_sodium()}}.
}
\details{
Running this function will invalidate \emph{all} keys loaded with
the above functions.  It should not be needed very often.
}
\examples{

# Be careful - if you run this then all keys loaded from file will
# no longer work until reloaded
if (FALSE) {
  cyphr::session_key_refresh()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/openssl.R
\name{key_openssl}
\alias{key_openssl}
\title{Symmetric encryption with openssl}
\usage{
key_openssl(key, mode = "cbc")
}
\arguments{
\item{key}{An openssl aes key (i.e., an object of class \code{aes}).}

\item{mode}{The encryption mode to use.  Options are \code{cbc},
\code{ctr} and \code{gcm} (see the \code{openssl} package for
more details)}
}
\description{
Wrap an openssl symmetric (aes) key.  This can be used with the
functions \code{\link[=encrypt_data]{encrypt_data()}} and
\code{\link[=decrypt_data]{decrypt_data()}}, along with the higher level wrappers
\code{\link[=encrypt]{encrypt()}} and \code{\link[=decrypt]{decrypt()}}.  With a symmetric
key, everybody uses the same key for encryption and decryption.
}
\examples{
# Create a new key
key <- cyphr::key_openssl(openssl::aes_keygen())
key

# With this key encrypt a string
secret <- cyphr::encrypt_string("my secret string", key)
# And decrypt it again:
cyphr::decrypt_string(secret, key)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sodium.R
\name{keypair_sodium}
\alias{keypair_sodium}
\title{Asymmetric encryption with sodium}
\usage{
keypair_sodium(pub, key, authenticated = TRUE)
}
\arguments{
\item{pub}{A sodium public key.  This is either a raw vector of
length 32 or a path to file containing the contents of the key
(written by \code{\link[=writeBin]{writeBin()}}).}

\item{key}{A sodium private key.  This is either a raw vector of
length 32 or a path to file containing the contents of the key
(written by \code{\link[=writeBin]{writeBin()}}).}

\item{authenticated}{Logical, indicating if authenticated
encryption (via \code{\link[sodium:messaging]{sodium::auth_encrypt()}} /
\code{\link[sodium:messaging]{sodium::auth_decrypt()}}) should be used.  If \code{FALSE}
then \code{\link[sodium:simple]{sodium::simple_encrypt()}} /
\code{\link[sodium:simple]{sodium::simple_decrypt()}} will be used.  The difference is
that with \code{authenticated = TRUE} the message is signed with
your private key so that tampering with the message will be
detected.}
}
\description{
Wrap a pair of sodium keys for asymmetric encryption.  You should
pass your private key and the public key of the person that you
are communicating with.
}
\details{
\emph{NOTE}: the order here (pub, key) is very important; if the
wrong order is used you cannot decrypt things.  Unfortunately
because sodium keys are just byte sequences there is nothing to
distinguish the public and private keys so this is a pretty easy
mistake to make.
}
\examples{

# Generate two keypairs, one for Alice, and one for Bob
key_alice <- sodium::keygen()
pub_alice <- sodium::pubkey(key_alice)
key_bob <- sodium::keygen()
pub_bob <- sodium::pubkey(key_bob)

# Alice wants to send Bob a message so she creates a key pair with
# her private key and bob's public key (she does not have bob's
# private key).
pair_alice <- cyphr::keypair_sodium(pub = pub_bob, key = key_alice)

# She can then encrypt a secret message:
secret <- cyphr::encrypt_string("hi bob", pair_alice)
secret

# Bob wants to read the message so he creates a key pair using
# Alice's public key and his private key:
pair_bob <- cyphr::keypair_sodium(pub = pub_alice, key = key_bob)

cyphr::decrypt_string(secret, pair_bob)
}
\seealso{
\code{\link[=keypair_openssl]{keypair_openssl()}} for a similar function using
openssl keypairs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/openssl_keygen.R
\name{ssh_keygen}
\alias{ssh_keygen}
\title{Create ssh keypairs}
\usage{
ssh_keygen(path = tempfile(), password = TRUE, use_shell = FALSE)
}
\arguments{
\item{path}{A directory in which to create a keypair.  If the path
does not exist it will be created.}

\item{password}{The password for the key.  The default will prompt
interactively (but without echoing the password).  Other valid
options are \code{FALSE} (no password) or a string.}

\item{use_shell}{Try to use \code{ssh-keygen} (the shell utility)
rather than functions in the \code{openssl} package.  This will
be necessary on at least very old versions of OS/X (Yosemite and
older at least) where the keys generated by the \code{openssl}
package cannot be read by the system ssh commands (e.g.,
\code{ssh-add}).}
}
\value{
The \code{path}, invisibly.  This is useful in the case
where \code{path} is \code{\link[=tempfile]{tempfile()}}.
}
\description{
Create openssl key pairs in the manner of \code{ssh-keygen}(1).
In general this should not be used (generate keys yourself with
\code{ssh-keygen} at the command line.  However this is useful for
testing and demonstration so I have included it to make that
easier.  Once a keypair has been generated it can be used with
\code{\link[=keypair_openssl]{keypair_openssl()}}.
}
\examples{
# Generate a new key in a temporary directory:
path <- cyphr::ssh_keygen(password = FALSE)
dir(path) # will contain id_rsa and id_rsa.pub

# This key can now be used via keypair_openssl:
key <- cyphr::keypair_openssl(path, path)
secret <- cyphr::encrypt_string("hello", key)
cyphr::decrypt_string(secret, key)

# Cleanup
unlink(path, recursive = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/openssl.R
\name{keypair_openssl}
\alias{keypair_openssl}
\title{Asymmetric encryption with openssl}
\usage{
keypair_openssl(
  pub,
  key,
  envelope = TRUE,
  password = NULL,
  authenticated = TRUE
)
}
\arguments{
\item{pub}{An openssl public key.  Usually this will be the path
to the key, in which case it may either the path to a public key
or be the path to a directory containing a file
\code{id_rsa.pub}.  If \code{NULL}, then your public key will be
used (found via the environment variable \code{USER_PUBKEY},
then \verb{~/.ssh/id_rsa.pub}).  However, it is not that common
to use your own public key - typically you want either the
sender of a message you are going to decrypt, or the recipient
of a message you want to send.}

\item{key}{An openssl private key.  Usually this will be the path
to the key, in which case it may either the path to a private
key or be the path to a directory containing a file.  You may
specify \code{NULL} here, in which case the environment variable
\code{USER_KEY} is checked and if that is not defined then
\verb{~/.ssh/id_rsa} will be used.}

\item{envelope}{A logical indicating if "envelope" encryption
functions should be used.  If so, then we use
\code{\link[openssl:encrypt_envelope]{openssl::encrypt_envelope()}} and
\code{\link[openssl:encrypt_envelope]{openssl::decrypt_envelope()}}.  If \code{FALSE} then we use
\code{\link[openssl:rsa_encrypt]{openssl::rsa_encrypt()}} and \code{\link[openssl:rsa_encrypt]{openssl::rsa_decrypt()}}.
See the openssl docs for further details.  The main effect of
this is that using \code{envelope = TRUE} will allow you to
encrypt much larger data than \code{envelope = FALSE}; this is
because openssl asymmetric encryption can only encrypt data up
to the size of the key itself.}

\item{password}{A password for the private key.  If \code{NULL}
then you will be prompted interactively for your password, and
if a string then that string will be used as the password (but
be careful in scripts!)}

\item{authenticated}{Logical, indicating if the result should be
signed with your public key.  If \code{TRUE} then your key will
be verified on decryption.  This provides tampering detection.}
}
\description{
Wrap a pair of openssl keys.  You should pass your private key and
the public key of the person that you are communicating with.
}
\examples{

# Note this uses password = FALSE for use in examples only, but
# this should not be done for any data you actually care about.

# Note that the vignette contains much more information than this
# short example and should be referred to before using these
# functions.

# Generate two keypairs, one for Alice, and one for Bob
path_alice <- tempfile()
path_bob <- tempfile()
cyphr::ssh_keygen(path_alice, password = FALSE)
cyphr::ssh_keygen(path_bob, password = FALSE)

# Alice wants to send Bob a message so she creates a key pair with
# her private key and bob's public key (she does not have bob's
# private key).
pair_alice <- cyphr::keypair_openssl(pub = path_bob, key = path_alice)

# She can then encrypt a secret message:
secret <- cyphr::encrypt_string("hi bob", pair_alice)
secret

# Bob wants to read the message so he creates a key pair using
# Alice's public key and his private key:
pair_bob <- cyphr::keypair_openssl(pub = path_alice, key = path_bob)

cyphr::decrypt_string(secret, pair_bob)

# Clean up
unlink(path_alice, recursive = TRUE)
unlink(path_bob, recursive = TRUE)
}
\seealso{
\code{\link[=keypair_sodium]{keypair_sodium()}} for a similar function using
sodium keypairs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/encrypt_wrapper.R
\name{rewrite_register}
\alias{rewrite_register}
\title{Register functions to work with encrypt/decrypt}
\usage{
rewrite_register(package, name, arg, fn = NULL)
}
\arguments{
\item{package}{The name of the package with the function to
support (as a scalar character).  If your function has no
package (e.g., a function you are working on outside of a
package, use "" as the name).}

\item{name}{The name of the function to support.}

\item{arg}{The name of the argument in the target function that
refers to the file that should be encrypted or decrypted.  This
is the value you would pass through to \code{file_arg} in
\link{encrypt}.}

\item{fn}{Optional (and should be rare) argument used to work
around functions that pass all their arguments through to a
second function as dots.  This is how \code{read.csv} works.  If
needed this function is a length-2 character vector in the form
"package", "name" with the actual function that is used.  But
this should be very rare!}
}
\description{
Add information about argument rewriting so that they can be used
with \link{encrypt} and \link{decrypt}.
}
\details{
If your package uses cyphr, it might be useful to add this as
an \code{.onLoad()} hook.
}
\examples{
# The saveRDS function is already supported.  But if we wanted to
# support it we could look at the arguments for the function:
args(saveRDS)
# The 'file' argument is the one that refers to the filename, so
# we'd write:
cyphr::rewrite_register("base", "saveRDS", "file")
# It's non-API but you can see what is supported in the package by
# looking at
ls(cyphr:::db)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cyphr.R
\docType{package}
\name{cyphr}
\alias{cyphr}
\title{High Level Encryption Wrappers}
\description{
Encryption wrappers, using low-level support from sodium and openssl.
}
\details{
It is \emph{strongly} recommended that you read \emph{both} vignettes before
attempting to use \code{cyphr}.
\itemize{
\item \href{https://docs.ropensci.org/cyphr/articles/cyphr.html}{introduction};
in R: \code{vignette("cyphr", package = "cyphr")}
\item \href{https://docs.ropensci.org/cyphr/articles/data.html}{data vignette};
in R: \code{vignette("data", package = "cyphr")}
}
}
\author{
Rich FitzJohn (rich.fitzjohn@gmail.com)
}
