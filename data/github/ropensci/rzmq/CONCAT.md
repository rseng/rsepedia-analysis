# rzmq

> R Bindings for 'ZeroMQ'

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/rzmq?branch=master)](https://ci.appveyor.com/project/jeroen/rzmq)
[![Package-License](http://img.shields.io/badge/license-GPL--3-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html) [![CRAN](http://www.r-pkg.org/badges/version/rzmq)](https://cran.r-project.org/package=rzmq) [![Downloads](http://cranlogs.r-pkg.org/badges/rzmq?color=brightgreen)](http://www.r-pkg.org/pkg/rzmq)

Interface to the 'ZeroMQ' lightweight messaging kernel (see <http://www.zeromq.org/> for more information).


## Features

rzmq is a message queue for serialized R objects.
* rzmq implements most the standard socket pairs that ZMQ offers.
* ZMQ devices are not implemented yet, nor is zmq_poll.
* Look for more features shortly.

## Installation

Binary packages for __OS-X__ or __Windows__ can be installed directly from CRAN:

```r
install.packages("rzmq")
```

## Build from source

Installation from source requires [`ZeroMQ`](http://zeromq.org/area:download). On __Debian__ or __Ubuntu__ use [libzmq3-dev](https://packages.debian.org/testing/libzmq3-dev):

```
sudo apt-get install -y libzmq3-dev
```

On __Fedora__ we need [zeromq-devel](https://apps.fedoraproject.org/packages/zeromq-devel):

```
sudo yum install zeromq-devel
````

On __CentOS / RHEL__ we install [zeromq3-devel](https://apps.fedoraproject.org/packages/zeromq3-devel) via EPEL:

```
sudo yum install epel-release
sudo yum install zeromq3-devel
```

On __OS-X__ use [zeromq](https://github.com/Homebrew/homebrew-core/blob/master/Formula/zeromq.rb) from Homebrew:

```
brew install zeromq
```


## Usage

A minimal example of remote execution.

execute this R script on the remote server:
```{.r}
#!/usr/bin/env Rscript
library(rzmq)
context = init.context()
socket = init.socket(context,"ZMQ_REP")
bind.socket(socket,"tcp://*:5555")
while(1) {
    msg = receive.socket(socket);
    fun <- msg$fun
    args <- msg$args
    print(args)
    ans <- do.call(fun,args)
    send.socket(socket,ans);
}
```	

and execute this bit locally:
```{.r}
library(rzmq)

remote.exec <- function(socket,fun,...) {
    send.socket(socket,data=list(fun=fun,args=list(...)))
    receive.socket(socket)
}

substitute(expr)
context = init.context()
socket = init.socket(context,"ZMQ_REQ")
connect.socket(socket,"tcp://localhost:5555")

ans <- remote.exec(socket,sqrt,10000)
```
\name{receive.multipart}
\alias{receive.multipart}
\title{Receive multipart ZMQ message}
\usage{
receive.multipart(socket)
}
\arguments{
  \item{socket}{The ZMQ socket from which to receive data}
}
\description{
  Returns a list of raw vectors for the parts of a multipart message.
}

\name{socket.options}
\alias{set.hwm}
\alias{set.swap}
\alias{set.affinity}
\alias{set.identity}
\alias{subscribe}
\alias{unsubscribe}
\alias{set.rate}
\alias{set.recovery.ivl}
\alias{set.recovery.ivl.msec}
\alias{set.mcast.loop}
\alias{set.sndbuf}
\alias{set.rcvbuf}
\alias{set.linger}
\alias{set.reconnect.ivl}
\alias{set.zmq.backlog}
\alias{set.reconnect.ivl.max}
\alias{get.rcvmore}
\alias{get.last.endpoint}
\alias{get.send.timeout}
\alias{set.send.timeout}

\title{
  set a socket option.
}
\description{

The zmq_setsockopt() function shall set the option specified by the
option_name argument to the value pointed to by the option_value
argument for the ZMQ socket pointed to by the socket argument.
}
\usage{
set.hwm(socket, option.value)
set.swap(socket, option.value)
set.affinity(socket, option.value)
set.identity(socket, option.value)
subscribe(socket, option.value)
unsubscribe(socket, option.value)
set.rate(socket, option.value)
set.recovery.ivl(socket, option.value)
set.recovery.ivl.msec(socket, option.value)
set.mcast.loop(socket, option.value)
set.sndbuf(socket, option.value)
set.rcvbuf(socket, option.value)
set.linger(socket, option.value)
set.reconnect.ivl(socket, option.value)
set.zmq.backlog(socket, option.value)
set.reconnect.ivl.max(socket, option.value)
get.rcvmore(socket)
get.last.endpoint(socket)
get.send.timeout(socket)
set.send.timeout(socket, option.value)

}

\arguments{
  \item{socket}{a zmq socket object}
  \item{option.value}{the new option value to bet set}
}
\value{
  a boolean indicating success or failure of the operation or in the
  case of getsocketoptions, the value of the requsted option.
}
\references{
  http://www.zeromq.org
  http://api.zeromq.org
  http://zguide.zeromq.org/page:all
}
\author{
  ZMQ was written by Martin Sustrik <sustrik@250bpm.com> and Martin Lucina <mato@kotelna.sk>.
  rzmq was written by Whit Armstrong.
}

\seealso{
  \code{\link{connect.socket},\link{bind.socket},\link{receive.socket},\link{send.socket},\link{poll.socket}}
}
\examples{\dontrun{

library(rzmq)
context = init.context()
socket = init.socket(context,"ZMQ_REQ")

set.hwm(socket, 1L)
set.swap(socket, 100L)
set.identity(socket, "big.ass.socket")
}}
\keyword{utilities}
\name{receive.socket}
\alias{receive.socket}
\alias{receive.null.msg}
\alias{receive.string}
\alias{receive.int}
\alias{receive.double}
\title{
Receive a message from the socket referenced by the socket argument.
}
\description{
The zmq_recv() function shall receive a message from the socket
referenced by the socket argument. If there are no messages available
on the specified socket, by default the function shall block until the request
can be satisfied.
A non-blocking receive can be obtained by setting dont.wait to TRUE
If there are no messages available on the specified socket, the
receive.socket() call will return NULL immediately.

}
\usage{
receive.socket(socket, unserialize=TRUE, dont.wait=FALSE)
receive.null.msg(socket)
receive.string(socket)
receive.int(socket)
receive.double(socket)
}

\arguments{
\item{socket}{a zmq socket object}
\item{unserialize}{whether to call unserialize on the received data}
\item{dont.wait}{defaults to false, for blocking receive. Set to TRUE for non-blocking receive.}
}
\value{
  the value sent from the remote server or NULL on failure.
  If dont.wait was TRUE and a message was not immediately
  available for receipt, NULL is returned and get.zmq.errno() is set to 11
  or get.zmq.strerror() is set to EAGAIN.
}
\references{
http://www.zeromq.org
http://api.zeromq.org
http://zguide.zeromq.org/page:all
}
\author{
ZMQ was written by Martin Sustrik <sustrik@250bpm.com> and Martin Lucina <mato@kotelna.sk>.
rzmq was written by Whit Armstrong.
}
\seealso{
  \code{\link{connect.socket},\link{bind.socket},\link{receive.socket},\link{send.socket},\link{poll.socket}}
}
\examples{\dontrun{
library(rzmq)

remote.exec <- function(out.socket,in.socket,fun,...) {
    send.socket(out.socket,data=list(fun=fun,args=list(...)))
    receive.socket(in.socket)
}

context = init.context()
out.socket = init.socket(context,"ZMQ_PUSH")
bind.socket(out.socket,"tcp://*:5557")

in.socket = init.socket(context,"ZMQ_PULL")
bind.socket(in.socket,"tcp://*:5558")


myfun <- function(x) {
    sum(abs(x))
}

remote.exec(out.socket,in.socket,myfun,rnorm(1e3))

}}

\keyword{utilities}
\name{bind.socket}
\alias{bind.socket}
\title{
Create an endpoint for accepting connections and bind it to the socket referenced by the socket argument.
}
\description{
The zmq_bind() function shall create an endpoint for accepting connections and bind it to the socket referenced by the socket argument.

The endpoint argument is a string consisting of two parts as follows: transport ://address. The transport part specifies the underlying transport protocol to use. The meaning of the address part is specific to the underlying transport protocol selected.

The following transports are defined:

inproc
local in-process (inter-thread) communication transport, see zmq_inproc(7)
ipc
local inter-process communication transport, see zmq_ipc(7)
tcp
unicast transport using TCP, see zmq_tcp(7)
pgm, epgm
reliable multicast transport using PGM, see zmq_pgm(7)
With the exception of ZMQ_PAIR sockets, a single socket may be connected to multiple endpoints using zmq_connect(), while simultaneously accepting incoming connections from multiple endpoints bound to the socket using zmq_bind(). Refer to zmq_socket(3) for a description of the exact semantics involved when connecting or binding a socket to multiple endpoints.
}
\usage{
bind.socket(socket, address)
}
\arguments{
  \item{socket}{a zmq socket object.}
  \item{address}{a transport as described above.}
}
\value{TRUE if operation succeeds or FALSE if the operation fails}
\references{
http://www.zeromq.org
http://api.zeromq.org
http://zguide.zeromq.org/page:all
}
\author{
ZMQ was written by Martin Sustrik <sustrik@250bpm.com> and Martin Lucina <mato@kotelna.sk>.
rzmq was written by Whit Armstrong.
}


\seealso{
  \code{\link{connect.socket},\link{bind.socket},\link{receive.socket},\link{send.socket},\link{poll.socket}}
}
\examples{\dontrun{

library(rzmq)
context = init.context()
in.socket = init.socket(context,"ZMQ_PULL")
bind.socket(in.socket,"tcp://*:5557")

out.socket = init.socket(context,"ZMQ_PUSH")
bind.socket(out.socket,"tcp://*:5558")
}}

\keyword{utilities}
\name{connect.socket}
\alias{connect.socket}
\alias{disconnect.socket}
\title{
Connect the socket referenced by the socket argument to the endpoint specified by the endpoint argument.
}
\description{
  The zmq_connect() function shall connect the socket referenced by the socket argument to the endpoint specified by the endpoint argument.

The endpoint argument is a string consisting of two parts as follows: transport ://address. The transport part specifies the underlying transport protocol to use. The meaning of the address part is specific to the underlying transport protocol selected.

The following transports are defined:

inproc
local in-process (inter-thread) communication transport, see zmq_inproc(7)
ipc
local inter-process communication transport, see zmq_ipc(7)
tcp
unicast transport using TCP, see zmq_tcp(7)
pgm, epgm
reliable multicast transport using PGM, see zmq_pgm(7)
With the exception of ZMQ_PAIR sockets, a single socket may be connected to multiple endpoints using zmq_connect(), while simultaneously accepting incoming connections from multiple endpoints bound to the socket using zmq_bind(). Refer to zmq_socket(3) for a description of the exact semantics involved when connecting or binding a socket to multiple endpoints.
}
\usage{
connect.socket(socket, address)
disconnect.socket(socket, address)
}

\arguments{
  \item{socket}{a zmq socket object.}
  \item{address}{a transport as described above.}
}
\value{TRUE if operation succeeds or FALSE if the operation fails}
\references{
http://www.zeromq.org
http://api.zeromq.org
http://zguide.zeromq.org/page:all
}
\author{
ZMQ was written by Martin Sustrik <sustrik@250bpm.com> and Martin Lucina <mato@kotelna.sk>.
rzmq was written by Whit Armstrong.
}
\seealso{
  \code{\link{connect.socket},\link{bind.socket},\link{receive.socket},\link{send.socket},\link{poll.socket}}
}
\examples{\dontrun{
library(rzmq)
context = init.context()
in.socket = init.socket(context,"ZMQ_PULL")
bind.socket(in.socket,"tcp://*:5557")

out.socket = init.socket(context,"ZMQ_PUSH")
bind.socket(out.socket,"tcp://*:5558")
}}
\keyword{utilities}
\name{init.context}
\alias{init.context}
\alias{init.socket}
\title{
  initailize zmq context and zmq socket
}
\description{
  initialize zmq context and zmq socket for to be used for further zmq operations.
}
\usage{
init.context(threads=1L)
init.socket(context, socket.type)
}

\arguments{
  \item{threads}{number of threads for the context to use}
  \item{context}{a zmq context object.}
  \item{socket.type}{ The ZMQ socket type requested
    e.g. ZMQ_REQ,ZMQ_REP,ZMQ_PULL,ZMQ_PUSH, etc.}
}
\value{
  \code{init.context}{ returns a zmq context object.}
  \code{init.socket}{returns a zmq socket object.}
}
\references{
  http://www.zeromq.org
  http://api.zeromq.org
  http://zguide.zeromq.org/page:all
}
\author{
  ZMQ was written by Martin Sustrik <sustrik@250bpm.com> and Martin Lucina <mato@kotelna.sk>.
  rzmq was written by Whit Armstrong.
}

\seealso{
  \code{\link{connect.socket},\link{bind.socket},\link{receive.socket},\link{send.socket},\link{poll.socket}}
}

\examples{\dontrun{

library(rzmq)
context = init.context()
in.socket = init.socket(context,"ZMQ_PULL")
}}
\keyword{utilities}
\name{zmq.version}
\alias{zmq.version}
\title{
  get version of libzmq
}
\description{
  return the version string of the system zmq library
}
\usage{
zmq.version()
}

\value{
  a string of the following format: major.minor.patch
}
\references{
  http://www.zeromq.org
  http://api.zeromq.org
  http://zguide.zeromq.org/page:all
}
\author{
  ZMQ was written by Martin Sustrik <sustrik@250bpm.com> and Martin Lucina <mato@kotelna.sk>.
  rzmq was written by Whit Armstrong.
}

\seealso{
  \code{\link{connect.socket},\link{bind.socket},\link{receive.socket},\link{send.socket}}
}

\examples{\dontrun{

library(rzmq)
zmq.version()
}}
\keyword{utilities}
\name{send.multipart}
\alias{send.multipart}
\title{Send multipart ZMQ message.}
\usage{
send.multipart(socket, parts)
}
\arguments{
  \item{socket}{The ZMQ socket on which to send data}

  \item{parts}{A list of raw vectors; each component will be sent
  as one part of the message, in the order of the list}
}
\description{
  Queue a list of raw vectors to be sent as a series of ZMQ message parts. Each
  part before the last will be sent with the SNDMORE flag.
}
\name{zmq.error}
\alias{zmq.errno}
\alias{zmq.strerror}
\title{
  get libzmq error numbers and error strings
}
\description{
  return the error number or error description after a zmq call
}
\usage{
zmq.errno()
zmq.strerror()
}

\value{
  an integer for zmq.errno or a string for zmq.strerror
}
\references{
  http://www.zeromq.org
  http://api.zeromq.org
  http://zguide.zeromq.org/page:all
}
\author{
  ZMQ was written by Martin Sustrik <sustrik@250bpm.com> and Martin Lucina <mato@kotelna.sk>.
  rzmq was written by Whit Armstrong.
}

\seealso{
  \code{\link{connect.socket},\link{bind.socket},\link{receive.socket},\link{send.socket}}
}

\examples{\dontrun{

library(rzmq)
zmq.errno()
zmq.strerror()
}}
\keyword{utilities}
\name{send.socket}
\alias{send.socket}
\alias{send.null.msg}
\alias{send.raw.string}
\alias{send.message.object}
\title{
  send a message.
}
\description{
  Queue the message referenced by the msg argument to be sent to the socket referenced by the socket argument. 

  A successful invocation of send.socket does not indicate that the message has been transmitted to the network, only that it has been queued on the socket and ZMQ has assumed responsibility for the message.
}
\usage{
send.socket(socket, data, send.more=FALSE, serialize=TRUE, xdr=.Platform$endian=="big")
send.null.msg(socket, send.more=FALSE)
send.raw.string(socket,data,send.more=FALSE)
}

\arguments{
  \item{socket}{a zmq socket object}
  \item{data}{the R object to be sent}
  \item{send.more}{whether this message has more frames to be sent}
  \item{serialize}{whether to call serialize before sending the data}
  \item{xdr}{passed directly to serialize command if serialize is requested}
}
\value{
  a boolean indicating success or failure of the operation.
}
\references{
  http://www.zeromq.org
  http://api.zeromq.org
  http://zguide.zeromq.org/page:all
}
\author{
  ZMQ was written by Martin Sustrik <sustrik@250bpm.com> and Martin Lucina <mato@kotelna.sk>.
  rzmq was written by Whit Armstrong.
}

\seealso{
  \code{\link{connect.socket},\link{bind.socket},\link{receive.socket},\link{send.socket},\link{poll.socket}}
}
\examples{\dontrun{

## remote execution server in rzmq
library(rzmq)
context = init.context()
in.socket = init.socket(context,"ZMQ_PULL")
bind.socket(in.socket,"tcp://*:5557")

out.socket = init.socket(context,"ZMQ_PUSH")
bind.socket(out.socket,"tcp://*:5558")

while(1) {
   msg = receive.socket(in.socket)
   fun <- msg$fun
   args <- msg$args
   print(args)
   ans <- do.call(fun,args)
   send.socket(out.socket,ans)
}
}}
\keyword{utilities}
\name{poll.socket}
\alias{poll.socket}
\title{Polls a list of sockets, waiting for the presence of a nonblocking read, write, or error event.}
\description{The zmq_poll() function shall poll a list of a sockets for either read, write, or error conditions subject to a millisecond resolution timeout.}
\usage{
    poll.socket(sockets, events, timeout=0L)
}

\arguments{
    \item{sockets}{a list of zmq socket objects.}
    \item{events}{a list of character vectors containing one or more events in \{read, write, error\}. The first element in the list corresponds to the first zmq socket, and so on...}
    \item{timeout}{the numbers of seconds to wait for events. Fractional seconds are supported. ZeroMQ guarantees at most millisecond resolution. A timeout of -1L blocks until an event occurs; a timeout of 0L is non-blocking.}
}
\value{A list of pairlists corresponding to the polled zmq sockets. Each list has one of more tags from \{read, write, error\} with logical values indicating the results of the poll operation.}

\references{
http://www.zeromq.org
http://api.zeromq.org
http://zguide.zeromq.org/page:all
}

\author{
ZMQ was written by Martin Sustrik <sustrik@250bpm.com> and Martin Lucina <mato@kotelna.sk>.
rzmq was written by Whit Armstrong.
}
\seealso{
  \code{\link{connect.socket},\link{bind.socket},\link{receive.socket},\link{send.socket},\link{poll.socket}}
}
\examples{\dontrun{
library(rzmq)

# Create a set of REP-REQ sockets that
# have a Send, Receive, Send, Receive, ...
# pattern.
context = init.context()
in.socket = init.socket(context,"ZMQ_REP")
bind.socket(in.socket,"tcp://*:5557")

out.socket = init.socket(context,"ZMQ_REQ")
connect.socket(out.socket,"tcp://*:5557")

# Poll the REP and REQ sockets for all events.
events <- poll.socket(list(in.socket, out.socket),
                      list(c("read", "write", "error"),
                           c("read", "write", "error")),
                      timeout=0L)

# The REQ socket is writable without blocking.
paste("Is upstream REP socket readable without blocking?", events[[1]]$read)
paste("Is upstream REP socket writable without blocking?", events[[1]]$write)
paste("Is downstream REQ socket readable without blocking?", events[[2]]$read)
paste("Is downstream REQ socket writable without blocking?", events[[2]]$write)

# Send a message to the REP socket from the REQ socket. The
# REQ socket must respond before the REP socket can send
# another message.
send.socket(out.socket, "Hello World")

events <- poll.socket(list(in.socket, out.socket),
                      list(c("read", "write", "error"),
                           c("read", "write", "error")),
                      timeout=0L)

# The incoming message is readable on the REP socket.
paste("Is upstream REP socket readable without blocking?", events[[1]]$read)
paste("Is upstream REP socket writable without blocking?", events[[1]]$write)
paste("Is downstream REQ socket readable without blocking?", events[[2]]$read)
paste("Is downstream REQ socket writable without blocking?", events[[2]]$write)

receive.socket(in.socket)

events <- poll.socket(list(in.socket, out.socket),
                      list(c("read", "write", "error"),
                           c("read", "write", "error")),
                      timeout=0L)

# The REQ socket is waiting for a response from the REP socket. 
paste("Is upstream REP socket readable without blocking?", events[[1]]$read)
paste("Is upstream REP socket writable without blocking?", events[[1]]$write)
paste("Is downstream REQ socket readable without blocking?", events[[2]]$read)
paste("Is downstream REQ socket writable without blocking?", events[[2]]$write)

send.socket(in.socket, "Greetings")

events <- poll.socket(list(in.socket, out.socket),
                      list(c("read", "write", "error"),
                           c("read", "write", "error")),
                      timeout=0L)

# The REP response is waiting to be read on the REQ socket.
paste("Is upstream REP socket readable without blocking?", events[[1]]$read)
paste("Is upstream REP socket writable without blocking?", events[[1]]$write)
paste("Is downstream REQ socket readable without blocking?", events[[2]]$read)
paste("Is downstream REQ socket writable without blocking?", events[[2]]$write)

# Complete the REP-REQ transaction cycle by reading
# the REP response.
receive.socket(out.socket)
}}
\keyword{utilities}
\name{init.message}
\alias{init.message}
\title{
  create a message object.
}
\description{
  Create a ZeroMQ message object that can be sent multiple times
}
\usage{
init.message(data, serialize=TRUE, xdr=.Platform$endian=="big")
}

\arguments{
  \item{data}{the R object to be sent}
  \item{serialize}{whether to call serialize before sending the data}
  \item{xdr}{passed directly to serialize command if serialize is requested}
}
\value{
  a ZeroMQ message object as external pointer
}
\references{
  http://www.zeromq.org
  http://api.zeromq.org
  http://zguide.zeromq.org/page:all
}
\author{
  ZMQ was written by Martin Sustrik <sustrik@250bpm.com> and Martin Lucina <mato@kotelna.sk>.
  rzmq was written by Whit Armstrong.
}

\seealso{
  \code{\link{send.message.object}}
}
\examples{\dontrun{

## remote execution server in rzmq
library(rzmq)
data = list(x=5)
msg = init.message(data)
}}
\keyword{utilities}
