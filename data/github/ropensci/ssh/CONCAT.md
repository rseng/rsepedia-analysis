# SSH

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/ropensci/ssh.svg?branch=master)](https://travis-ci.org/ropensci/ssh)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/ssh?branch=master&svg=true)](https://ci.appveyor.com/project/jeroen/ssh)
[![Coverage Status](https://codecov.io/github/ropensci/ssh/coverage.svg?branch=master)](https://codecov.io/github/ropensci/ssh?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ssh)](http://cran.r-project.org/package=ssh)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/ssh)](http://cran.r-project.org/web/packages/ssh/index.html)

> Secure Shell (SSH) Client for R

## Installation

This package is available on CRAN and can be installed via:

```r
install.packages('ssh')
```

Alternatively it can be installed from source using `devtools`:

```r
remotes::install_github('ropensci/ssh')
```

Installation from source on MacOS or Linux requires [`libssh`](https://www.libssh.org/) (the original `libssh`, __not__ the unrelated `libssh2` library). On __Debian__ or __Ubuntu__ use [libssh-dev](https://packages.debian.org/testing/libssh-dev):

```
sudo apt-get install -y libssh-dev
```

On __Fedora__ we need [libssh-devel](https://src.fedoraproject.org/rpms/libssh):

```
sudo yum install libssh-devel
````

On __CentOS / RHEL__ we install [libssh-devel](https://src.fedoraproject.org/rpms/libssh) via EPEL:

```
sudo yum install epel-release
sudo yum install libssh-devel
```

On __OS-X__ use [libssh](https://github.com/Homebrew/homebrew-core/blob/master/Formula/libssh.rb) from Homebrew:

```
brew install libssh
```

Using __conda__ (need a conda R environment `conda create -n Renv r-base r-essentials`)

```
conda install --channel conda-forge r-ssh
```

If you have issues with the conda installation please submit an issue in [`conda-forge/r-ssh-feedstock`](https://github.com/conda-forge/r-ssh-feedstock/issues)
## Getting Started

First create an ssh session by connecting to an SSH server. You can either use private key or passphrase authentication: 

```{r}
session <- ssh_connect("jeroen@dev.opencpu.org")
```

You can use the session in subsequent ssh functions below.

### Run a command

Run a command or script on the host while streaming stdout and stderr directly to the client.

```{r}
ssh_exec_wait(session, command = c(
  'curl -fOL https://cloud.r-project.org/src/contrib/Archive/jsonlite/jsonlite_1.5.tar.gz',
  'R CMD check jsonlite_1.5.tar.gz',
  'rm -f jsonlite_1.5.tar.gz'
))
```

If you want to capture the stdout/stderr:

```r
out <- ssh_exec_internal(session, "R -e 'rnorm(100)'")
cat(rawToChar(out$stdout))
```
#### Using 'sudo'

Note that the exec functions are non interactive so they cannot prompt for a sudo password. A trick is to use `-S` which reads the password from stdin:

```
out <- ssh_exec_wait(session, 'echo "mypassword!" | sudo -s -S apt-get update -y')
```

Be very careful with hardcoding passwords!

### Uploading and Downloading via SCP

Upload and download files via SCP. Directories are automatically traversed as in `scp -r`.

```r
# Upload a file to the server
file_path <- R.home("COPYING")
scp_upload(session, file_path)
```

```r
# Download the file back and verify it is the same
scp_download(session, "COPYING", to = tempdir())
tools::md5sum(file_path)
tools::md5sum(file.path(tempdir(), "COPYING"))
```

### Create a Tunnel

Opens a port on your machine and tunnel all traffic to a custom target host via the SSH server.

```r
ssh_tunnel(session, port = 5555,target = "ds043942.mongolab.com:43942")
```

This function blocks while the tunnel is active. Use the tunnel by connecting to `localhost:5555` from a separate process. The tunnel can only be used once and will automatically be closed when the client disconnects.

### Disconnect


When you are done with the session you should disconnect:

```r
ssh_disconnect(session)
```

If you forgot to disconnect, the garbage collector will do so for you (with a warning).
---
title: "Secure Shell (SSH) Client for R"
date: "2019-03-26"
output:
  html_document:
    fig_caption: false
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    toc_depth: 2
vignette: >
  %\VignetteIndexEntry{Secure Shell (SSH) Client for R}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




## Installation

On Windows or MacOS you can install the binary package directly from CRAN:

```r
install.packages("ssh")
```

Installation from source requires [`libssh`](https://www.libssh.org/) (which is __not__ `libssh2`). On __Debian__ or __Ubuntu__ use [libssh-dev](https://packages.debian.org/testing/libssh-dev):

```
sudo apt-get install -y libssh-dev
```

On __Fedora__ we need [libssh-devel](https://src.fedoraproject.org/rpms/libssh):

```
sudo yum install libssh-devel
````

On __CentOS / RHEL__ we install [libssh-devel](https://src.fedoraproject.org/rpms/libssh) via EPEL:

```
sudo yum install epel-release
sudo yum install libssh-devel
```

On __OS-X__ use [libssh](https://github.com/Homebrew/homebrew-core/blob/master/Formula/libssh.rb) from Homebrew:

```
brew install libssh
```

## Connecting to an SSH server

First create an ssh session by connecting to an SSH server.



```r
session <- ssh_connect("jeroen@dev.opencpu.org")
print(session)
```

```
<ssh session>
connected: jeroen@dev.opencpu.org:22
server: 1c:86:dc:b9:77:2c:74:ee:2d:00:49:6f:2e:c8:b2:61:5e:14:53:73
```

Once established, a session is closed automatically by the garbage collector when the object goes out of scope or when R quits. You can also manually close it using `ssh_disconnect()` but this is not strictly needed.


## Authentication

The client attempts to use the following authentication methods (in this order) until one succeeds:

 1. try key from `privkey` argument in `ssh_connect()` if specified
 2. if ssh-agent is available, try private key from ssh-agent
 3. try user key specified in `~/.ssh/config` or any of the default locations: `~/.ssh/id_ed25519`, `~/.ssh/id_ecdsa`, `~/.ssh/id_rsa`, or `.ssh/id_dsa`.
 4. Try challenge-response password authentication (if permitted by the server)
 5. Try plain password authentication (if permitted by the server)

To debug authentication set verbosity to at least level 2 or 3:

```r
session <- ssh_connect("jeroen@dev.opencpu.org", verbose = 2)
```
```
ssh_socket_connect: Nonblocking connection socket: 7
ssh_connect: Socket connecting, now waiting for the callbacks to work
socket_callback_connected: Socket connection callback: 1 (0)
ssh_client_connection_callback: SSH server banner: SSH-2.0-OpenSSH_7.2p2 Ubuntu-4ubuntu2.4
ssh_analyze_banner: Analyzing banner: SSH-2.0-OpenSSH_7.2p2 Ubuntu-4ubuntu2.4
ssh_analyze_banner: We are talking to an OpenSSH client version: 7.2 (70200)
ssh_packet_dh_reply: Received SSH_KEXDH_REPLY
ssh_client_curve25519_reply: SSH_MSG_NEWKEYS sent
ssh_packet_newkeys: Received SSH_MSG_NEWKEYS
ssh_packet_newkeys: Signature verified and valid
ssh_packet_userauth_failure: Access denied. Authentication that can continue: publickey
ssh_packet_userauth_failure: Access denied. Authentication that can continue: publickey
ssh_agent_get_ident_count: Answer type: 12, expected answer: 12
ssh_userauth_publickey_auto: Successfully authenticated using /Users/jeroen/.ssh/id_rsa
```

Tools for setting up and debugging your ssh keys are provided via the [`credentials` package](https://cran.r-project.org/package=credentials).


```r
ssh_key_info()
```

```
$key
[1] "/Users/jeroen/.ssh/id_rsa"

$pubkey
[1] "ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQDwI/G2v42n/fVFbuyOL/7eh06MOK71djbe7wr8op0LcMbCiMOHOFdMRpb+UcwlfnXDYWwNZsfr3l4Hvfbdi5CesDwZmx4FQ4BC4e+2IRqCbwknVwdIMfBwLWd9HNs8hJ8M7YW4lt0TucMtVkzGN4wFR+onvzyU0VSiSiT1JQUY+7g4YMJlmhKN3lES+as1RFMExLosFlzydD7o35nMuKp1VzdkUaWm3xjFCTl7MVl6cbCWh7m3IDl95tTT0jH+pxvK2Vgb3u+r+tH7dWbbjm9tkwH7NS5HwEjy08cKUv+8BUnDLC27EuK2/k+Ro6BsEDHmq0z0nf6FgpTXQ0iJxnUZ"
```

## Execute Script or Command

Run a command or script on the host and block while it runs. By default stdout and stderr are steamed directly back to the client. This function returns the exit status of the remote command (hence it does not automatically error for an unsuccessful exit status). 


```r
out <- ssh_exec_wait(session, command = 'whoami')
```

```
jeroen
```

```r
print(out)
```

```
[1] 0
```

You can also run a script that consists of multiple commands.


```r
ssh_exec_wait(session, command = c(
  'curl -O https://cran.r-project.org/src/contrib/Archive/jsonlite/jsonlite_1.4.tar.gz',
  'R CMD check jsonlite_1.4.tar.gz',
  'rm -f jsonlite_1.4.tar.gz'
))
```
```
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 1071k  100 1071k    0     0   654k      0  0:00:01  0:00:01 --:--:--  654k
* using log directory '/home/jeroen/jsonlite.Rcheck'
* using R version 3.4.3 (2017-11-30)
* using platform: x86_64-pc-linux-gnu (64-bit)
* using session charset: ASCII
* checking for file 'jsonlite/DESCRIPTION' ... OK
* this is package 'jsonlite' version '1.4'
* checking package namespace information ... OK
* checking package dependencies ...
...
```

### Capturing output

The `ssh_exec_internal()` is a convenient wrapper for `ssh_exec_wait()` which buffers the output steams and returns them as a raw vector. Also it raises an error by default when the remote command was not successful.


```r
out <- ssh_exec_internal(session, "R -e 'rnorm(10)'")
print(out$status)
```

```
[1] 0
```

```r
cat(rawToChar(out$stdout))
```

```

R version 3.5.3 (2019-03-11) -- "Great Truth"
Copyright (C) 2019 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> rnorm(10)
 [1] -0.4657579  0.1628843 -1.8308908 -1.4333284 -1.6689843 -0.3603318
 [7] -0.3299621 -0.2162147  1.0788035 -0.9536821
> 
> 
```

This function is very useful if you are running a remote command and want to use it's output as if you had executed it locally.

### Using sudo

Note that the exec functions are non interactive so they cannot prompt for a sudo password. A trick is to use `-S` which reads the password from stdin:

```r
command <- 'echo "mypassword!" | sudo -s -S apt-get update -y'
out <- ssh_exec_wait(session, command)
```

Be very careful with hardcoding passwords!

## Transfer Files via SCP

Upload and download files via SCP. Directories are automatically traversed as in `scp -r`.

```r
# Upload a file to the server
file_path <- R.home("COPYING")
scp_upload(session, file_path)
```

```r
# Download the file back and verify it is the same
scp_download(session, "COPYING", to = tempdir())
tools::md5sum(file_path)
tools::md5sum(file.path(tempdir(), "COPYING"))
```

## Hosting a Tunnel

Opens a port on your machine and tunnel all traffic to a custom target host via the SSH server.

```r
ssh_tunnel(session, port = 5555, target = "ds043942.mongolab.com:43942")
```

This function blocks while the tunnel is active. Use the tunnel by connecting to `localhost:5555` from a separate process. The tunnel can only be used once and will automatically be closed when the client disconnects.

## Disconnecting

When you are done with the session you should disconnect:


```r
ssh_disconnect(session)
```

If you forgot to disconnect, the garbage collector will do so for you (with a warning).

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connect.R
\name{ssh_connect}
\alias{ssh_connect}
\alias{ssh}
\alias{ssh_session_info}
\alias{ssh_info}
\alias{ssh_disconnect}
\alias{libssh_version}
\title{SSH Client}
\usage{
ssh_connect(host, keyfile = NULL, passwd = askpass, verbose = FALSE)

ssh_session_info(session)

ssh_disconnect(session)

libssh_version()
}
\arguments{
\item{host}{an ssh server string of the form \verb{[user@]hostname[:@port]}. An ipv6
hostname should be wrapped in brackets like this: \verb{[2001:db8::1]:80}.}

\item{keyfile}{path to private key file. Must be in OpenSSH format (see details)}

\item{passwd}{either a string or a callback function for password prompt}

\item{verbose}{either TRUE/FALSE or a value between 0 and 4 indicating log level:
0: no logging, 1: only warnings, 2: protocol, 3: packets or 4: full stack trace.}

\item{session}{ssh connection created with \code{\link[=ssh_connect]{ssh_connect()}}}
}
\description{
Create an ssh session using \code{ssh_connect()}. The session can be used to execute
commands, scp files or setup a tunnel.
}
\details{
The client first tries to authenticate using a private key, either from ssh-agent
or \verb{/.ssh/id_rsa} in the user home directory. If this fails it falls back on
challenge-response (interactive) and password auth if allowed by the server. The
\code{passwd} parameter can be used to provide a passphrase or a callback function to
ask prompt the user for the passphrase when needed.

The session will automatically be disconnected when the session object is removed
or when R exits but you can also use \code{\link[=ssh_disconnect]{ssh_disconnect()}}.

\strong{Windows users:} the private key must be in OpenSSH PEM format. If you open it in
a text editor the first line must be: \verb{-----BEGIN RSA PRIVATE KEY-----}.
To convert a Putty PKK key, open it in the \emph{PuttyGen} utility and go to
\emph{Conversions -> Export OpenSSH}.
}
\examples{
\dontrun{
session <- ssh_connect("dev.opencpu.org")
ssh_exec_wait(session, command = "whoami")
ssh_disconnect(session)
}
}
\seealso{
Other ssh: 
\code{\link{scp}},
\code{\link{ssh_credentials}},
\code{\link{ssh_exec}},
\code{\link{ssh_tunnel}()}
}
\concept{ssh}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tunnel.R
\name{ssh_tunnel}
\alias{ssh_tunnel}
\title{Create SSH tunnel}
\usage{
ssh_tunnel(session, port = 5555, target = "rainmaker.wunderground.com:23")
}
\arguments{
\item{session}{ssh connection created with \code{\link[=ssh_connect]{ssh_connect()}}}

\item{port}{integer of local port on which to listen for incoming connections}

\item{target}{string with target host and port to connect to via ssh tunnel}
}
\description{
Opens a port on your machine and tunnel all traffic to a custom target host via the
SSH server, for example to connect with a database server behind a firewall.
}
\details{
This function blocks while the tunnel is active. Use the tunnel by connecting to
\code{localhost:5555} from a separate process. Each tunnel can only be used once and will
automatically be closed when the client disconnects. It is intended to tunnel a single
connection, not as a long running proxy server.
}
\seealso{
Other ssh: 
\code{\link{scp}},
\code{\link{ssh_connect}()},
\code{\link{ssh_credentials}},
\code{\link{ssh_exec}}
}
\concept{ssh}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scp.R
\name{scp}
\alias{scp}
\alias{scp_download}
\alias{scp_upload}
\title{SCP (Secure Copy)}
\usage{
scp_download(session, files, to = ".", verbose = TRUE)

scp_upload(session, files, to = ".", verbose = TRUE)
}
\arguments{
\item{session}{ssh connection created with \code{\link[=ssh_connect]{ssh_connect()}}}

\item{files}{path to files or directory to transfer}

\item{to}{existing directory on the destination where \code{files} will be copied into}

\item{verbose}{print progress while copying files}
}
\description{
Upload and download files to/from the SSH server via the scp protocol.
Directories in the \code{files} argument are automatically traversed and
uploaded / downloaded recursively.
}
\details{
Note that the syntax is slightly different from the \code{scp} command line
tool because the \code{to} parameter is always a target \emph{directory} where
all \code{files} will be copied \strong{into}. If \code{to} does not exist, it will be
created.

The \code{files} parameter in \code{\link[=scp_upload]{scp_upload()}} is vectorised hence all files
and directories will be recursively uploaded \strong{into} the \code{to} directory.
For \code{\link[=scp_download]{scp_download()}} the \code{files} parameter must be a single string which
may contain wildcards.

The default path \code{to = "."} means that files get downloaded to the current
working directory and uploaded to the user home directory on the server.
}
\examples{
\dontrun{
# recursively upload files and directories
session <- ssh_connect("dev.opencpu.org")
files <- c(R.home("doc"), R.home("COPYING"))
scp_upload(session, files, to = "~/target")

# download it back
scp_download(session, "~/target/*", to = tempdir())

# delete it from the server
ssh_exec_wait(session, command = "rm -Rf ~/target")
ssh_disconnect(session)
}
}
\seealso{
Other ssh: 
\code{\link{ssh_connect}()},
\code{\link{ssh_credentials}},
\code{\link{ssh_exec}},
\code{\link{ssh_tunnel}()}
}
\concept{ssh}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/credentials.R
\docType{import}
\name{ssh_credentials}
\alias{ssh_credentials}
\alias{ssh_key_info}
\alias{reexports}
\alias{ssh_keygen}
\alias{ssh_read_key}
\alias{ssh_home}
\alias{ssh_agent_add}
\title{Managing Your SSH Key}
\seealso{
Other ssh: 
\code{\link{scp}},
\code{\link{ssh_connect}()},
\code{\link{ssh_exec}},
\code{\link{ssh_tunnel}()}
}
\concept{ssh}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{credentials}{\code{\link[credentials:ssh_credentials]{ssh_agent_add}}, \code{\link[credentials:ssh_credentials]{ssh_home}}, \code{\link[credentials:ssh_credentials]{ssh_key_info}}, \code{\link[credentials:ssh_credentials]{ssh_keygen}}, \code{\link[credentials:ssh_credentials]{ssh_read_key}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exec.R
\name{ssh_exec}
\alias{ssh_exec}
\alias{ssh_exec_wait}
\alias{ssh_exec_internal}
\title{Execute Remote Command}
\usage{
ssh_exec_wait(
  session,
  command = "whoami",
  std_out = stdout(),
  std_err = stderr()
)

ssh_exec_internal(session, command = "whoami", error = TRUE)
}
\arguments{
\item{session}{ssh connection created with \code{\link[=ssh_connect]{ssh_connect()}}}

\item{command}{The command or script to execute}

\item{std_out}{callback function, filename, or connection object to handle stdout stream}

\item{std_err}{callback function, filename, or connection object to handle stderr stream}

\item{error}{automatically raise an error if the exit status is non-zero}
}
\description{
Run a command or script on the host while streaming stdout and stderr directly
to the client.
}
\details{
The \code{\link[=ssh_exec_wait]{ssh_exec_wait()}} function is the remote equivalent of the local \code{\link[sys:exec]{sys::exec_wait()}}.
It runs a command or script on the ssh server and streams stdout and stderr to the client
to a file or connection. When done it returns the exit status for the remotely executed command.

Similarly \code{\link[=ssh_exec_internal]{ssh_exec_internal()}} is a small wrapper analogous to \code{\link[sys:exec]{sys::exec_internal()}}.
It buffers all stdout and stderr output into a raw vector and returns it in a list along with
the exit status. By default this function raises an error if the remote command was unsuccessful.
}
\examples{
\dontrun{
session <- ssh_connect("dev.opencpu.org")
ssh_exec_wait(session, command = c(
  'curl -O https://cran.r-project.org/src/contrib/jsonlite_1.5.tar.gz',
  'R CMD check jsonlite_1.5.tar.gz',
  'rm -f jsonlite_1.5.tar.gz'
))
ssh_disconnect(session)}
}
\seealso{
Other ssh: 
\code{\link{scp}},
\code{\link{ssh_connect}()},
\code{\link{ssh_credentials}},
\code{\link{ssh_tunnel}()}
}
\concept{ssh}
