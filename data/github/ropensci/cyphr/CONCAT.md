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
