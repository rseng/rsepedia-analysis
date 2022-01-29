![curl logo](https://curl.se/logo/curl-logo.svg)

[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/63/badge)](https://bestpractices.coreinfrastructure.org/projects/63)
[![Coverity passed](https://scan.coverity.com/projects/curl/badge.svg)](https://scan.coverity.com/projects/curl)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/l1vv31029huhf4g4?svg=true)](https://ci.appveyor.com/project/curlorg/curl)
[![Azure DevOps Build Status](https://dev.azure.com/daniel0244/curl/_apis/build/status/curl.curl?branchName=master)](https://dev.azure.com/daniel0244/curl/_build/latest?definitionId=1&branchName=master)
[![Cirrus Build Status](https://api.cirrus-ci.com/github/curl/curl.svg?branch=master)](https://cirrus-ci.com/github/curl/curl)
[![Backers on Open Collective](https://opencollective.com/curl/backers/badge.svg)](#backers)
[![Sponsors on Open Collective](https://opencollective.com/curl/sponsors/badge.svg)](#sponsors)
[![Language Grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/curl/curl.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/curl/curl/context:cpp)
[![Fuzzing Status](https://oss-fuzz-build-logs.storage.googleapis.com/badges/curl.svg)](https://bugs.chromium.org/p/oss-fuzz/issues/list?sort=-opened&can=1&q=proj:curl)

Curl is a command-line tool for transferring data specified with URL
syntax. Find out how to use curl by reading [the curl.1 man
page](https://curl.se/docs/manpage.html) or [the MANUAL
document](https://curl.se/docs/manual.html). Find out how to install Curl
by reading [the INSTALL document](https://curl.se/docs/install.html).

libcurl is the library curl is using to do its job. It is readily available to
be used by your software. Read [the libcurl.3 man
page](https://curl.se/libcurl/c/libcurl.html) to learn how.

You can find answers to the most frequent questions we get in [the FAQ
document](https://curl.se/docs/faq.html).

Study [the COPYING file](https://curl.se/docs/copyright.html) for
distribution terms.

## Contact

If you have problems, questions, ideas or suggestions, please contact us by
posting to a suitable [mailing list](https://curl.se/mail/).

All contributors to the project are listed in [the THANKS
document](https://curl.se/docs/thanks.html).

## Commercial support

For commercial support, maybe private and dedicated help with your problems or
applications using (lib)curl visit [the support page](https://curl.se/support.html).

## Website

Visit the [curl website](https://curl.se/) for the latest news and
downloads.

## Git

To download the latest source from the Git server do this:

    git clone https://github.com/curl/curl.git

(you will get a directory named curl created, filled with the source code)

## Security problems

Report suspected security problems via [our HackerOne
page](https://hackerone.com/curl) and not in public!

## Notice

Curl contains pieces of source code that is Copyright (c) 1998, 1999 Kungliga
Tekniska Högskolan. This notice is included here to comply with the
distribution terms.

## Backers

Thank you to all our backers! 🙏 [[Become a backer](https://opencollective.com/curl#backer)]

<a href="https://opencollective.com/curl#backers" target="_blank"><img src="https://opencollective.com/curl/backers.svg?width=890"></a>

## Sponsors

Support this project by becoming a sponsor. Your logo will show up here with a
link to your website. [[Become a
sponsor](https://opencollective.com/curl#sponsor)]

<a href="https://opencollective.com/curl/sponsor/0/website" target="_blank"><img src="https://opencollective.com/curl/sponsor/0/avatar.svg"></a>
<a href="https://opencollective.com/curl/sponsor/1/website" target="_blank"><img src="https://opencollective.com/curl/sponsor/1/avatar.svg"></a>
<a href="https://opencollective.com/curl/sponsor/2/website" target="_blank"><img src="https://opencollective.com/curl/sponsor/2/avatar.svg"></a>
<a href="https://opencollective.com/curl/sponsor/3/website" target="_blank"><img src="https://opencollective.com/curl/sponsor/3/avatar.svg"></a>
<a href="https://opencollective.com/curl/sponsor/4/website" target="_blank"><img src="https://opencollective.com/curl/sponsor/4/avatar.svg"></a>
<a href="https://opencollective.com/curl/sponsor/5/website" target="_blank"><img src="https://opencollective.com/curl/sponsor/5/avatar.svg"></a>
<a href="https://opencollective.com/curl/sponsor/6/website" target="_blank"><img src="https://opencollective.com/curl/sponsor/6/avatar.svg"></a>
<a href="https://opencollective.com/curl/sponsor/7/website" target="_blank"><img src="https://opencollective.com/curl/sponsor/7/avatar.svg"></a>
<a href="https://opencollective.com/curl/sponsor/8/website" target="_blank"><img src="https://opencollective.com/curl/sponsor/8/avatar.svg"></a>
<a href="https://opencollective.com/curl/sponsor/9/website" target="_blank"><img src="https://opencollective.com/curl/sponsor/9/avatar.svg"></a>
# Security Policy

See [docs/SECURITY-PROCESS.md](docs/SECURITY-PROCESS.md) for full details.

## Reporting a Vulnerability

If you have found or just suspect a security problem somewhere in curl or libcurl,
report it on [https://hackerone.com/curl](https://hackerone.com/curl).

We treat security issues with confidentiality until controlled and disclosed responsibly.
# include

Public include files for libcurl, external users.

They're all placed in the curl subdirectory here for better fit in any kind of
environment. You must include files from here using...

    #include <curl/curl.h>

... style and point the compiler's include path to the directory holding the
curl subdirectory. It makes it more likely to survive future modifications.

The public curl include files can be shared freely between different platforms
and different architectures.
# Continuous Integration for curl

Curl runs in many different environments, so every change is run against a large
number of test suites.

Every pull request is verified for each of the following:

 - ... it still builds, warning-free, on Linux and macOS, with both
   clang and gcc
 - ... it still builds fine on Windows with several MSVC versions
 - ... it still builds with cmake on Linux, with gcc and clang
 - ... it follows rudimentary code style rules
 - ... the test suite still runs 100% fine
 - ... the release tarball (the "dist") still works
 - ... it builds fine in-tree as well as out-of-tree
 - ... code coverage does not shrink drastically
 - ... different TLS backends still compile and pass tests

If the pull-request fails one of these tests, it will show up as a red X and
you are expected to fix the problem. If you do not understand when the issue is
or have other problems to fix the complaint, just ask and other project
members will likely be able to help out.

Consider the following table while looking at pull request failures:

 | CI platform as shown in PR          | State  | What to look at next       |
 | ----------------------------------- | ------ | -------------------------- |
 | CI / codeql                         | stable | quality check results      |
 | CI / fuzzing                        | stable | fuzzing results            |
 | CI / macos ...                      | stable | all errors and failures    |
 | Code scanning results / CodeQL      | stable | quality check results      |
 | FreeBSD FreeBSD: ...                | stable | all errors and failures    |
 | LGTM analysis: Python               | stable | new findings               |
 | LGTM analysis:  C/C++               | stable | new findings               |
 | buildbot/curl_winssl_ ...           | stable | all errors and failures    |
 | continuous-integration/appveyor/pr  | stable | all errors and failures    |
 | curl.curl (linux ...)               | stable | all errors and failures    |
 | curl.curl (windows ...)             | flaky  | repetitive errors/failures |
 | deepcode-ci-bot                     | stable | new findings               |
 | musedev                             | stable | new findings               |

Sometimes the tests fail due to a dependency service temporarily being offline
or otherwise unavailable, eg. package downloads. In this case you can just
try to update your pull requests to rerun the tests later as described below.

## CI servers

Here are the different CI environments that are currently in use, and how they
are configured:

### Github Actions

Github Actions runs the following tests:

- Mac OS tests with a variety of different compilation options
- Fuzz tests ([see tests/fuzz/README for
    more info](https://github.com/curl/curl/blob/master/tests/fuzz/README)).
- Curl compiled using the Rust TLS backend with Hyper
- CodeQL static analysis

These are each configured in different files in `.github/workflows`.

### Azure

The following tests are run in Microsoft Azure CI environment:

- Ubuntu tests with a variety of different compilation options.
- Windows tests with a variety of different compilation options.

These are all configured in `.azure-pipelines.yml`.

As of November 2021 @bagder and @mback2k are the only people with administrator
access to the Azure CI environment. Additional admins/group members can be added
on request.

### Appveyor

Appveyor runs a variety of different Windows builds, with different compilation
options.

As of November 2021 @bagder, @mback2k, @jay, @vszakats, @dfandrich and
@danielgustafsson have administrator access to the Appveyor CI environment.
Additional admins/group members can be added on request.

The tests are configured in appveyor.yml.

### Zuul

[Zuul](https://zuul-ci.org/) is an open source CI tool. A number of Curl tests
are run at [curl.zuul.vexxhost.dev](https://curl.zuul.vexxhost.dev/builds):

- Source code is formatted according to expectations (`make checksrc`).
- Curl compiled with a number of different TLS configurations (WolfSSL, rustls,
BoringSSL, etc).
- Curl compiled with different C compilers.

As of November 2021, the tests run (sometimes) but do not run consistently and
do not report results to the Github checks runner - you need to manually check
for failures. See [#7522](https://github.com/curl/curl/issues/7522) for more
information.

As of November 2021 Daniel Stenberg is the only person with administrator access
to the Zuul CI environment.

These are configured in `zuul.d` and have test runners in `scripts/zuul`.

### CircleCI

CircleCI runs a basic Linux test suite on Ubuntu for both x86 and ARM
processors. This is configured in `.circleci/config.yml`.

You can [view the full list of CI jobs on CircleCI's
website](https://app.circleci.com/pipelines/github/curl/curl).

@bagder has access to edit the "Project Settings" on that page.
Additional admins/group members can be added on request.

### Cirrus CI

Cirrus CI runs a basic test suite on FreeBSD and Windows. This is configured in
`.cirrus.yml`.

You can [view the full list of CI jobs on Cirrus CI's
website](https://cirrus-ci.com/github/curl/curl).

@bagder has access to edit the "Project Settings" on that page.
Additional admins/group members can be added on request.
# curl test suite file format

The curl test suite's file format is very simple and extensible, closely
resembling XML. All data for a single test case resides in a single ASCII
file. Labels mark the beginning and the end of all sections, and each label
must be written in its own line.  Comments are either XML-style (enclosed with
`<!--` and `-->`) or shell script style (beginning with `#`) and must appear
on their own lines and not alongside actual test data.  Most test data files
are syntactically valid XML, although a few files are not (lack of support for
character entities and the preservation of CR/LF characters at the end of
lines are the biggest differences).

Each test case source exists as a file matching the format
`tests/data/testNUM`, where NUM is the unique test number, and must begin with
a 'testcase' tag, which encompasses the remainder of the file.

# Preprocessing

When a test is to be executed, the source file is first preprocessed and
variables are substituted by the their respective contents and the output
version of the test file is stored as `log/testNUM`. That version is what will
be read and used by the test servers.

## Base64 Encoding

In the preprocess stage, a special instruction can be used to have runtests.pl
base64 encode a certain section and insert in the generated output file. This
is in particular good for test cases where the test tool is expected to pass
in base64 encoded content that might use dynamic information that is unique
for this particular test invocation, like the server port number.

To insert a base64 encoded string into the output, use this syntax:

    %b64[ data to encode ]b64%

The data to encode can then use any of the existing variables mentioned below,
or even percent-encoded individual bytes. As an example, insert the HTTP
server's port number (in ASCII) followed by a space and the hexadecimal byte
9a:

    %b64[%HTTPPORT %9a]b64%

## Hexadecimal decoding

In the preprocess stage, a special instruction can be used to have runtests.pl
generate a sequence of binary bytes.

To insert a sequence of bytes from a hex encoded string, use this syntax:

    %hex[ %XX-encoded data to decode ]hex%

For example, to insert the binary octets 0, 1 and 255 into the test file:

    %hex[ %00%01%FF ]hex%

## Repeat content

In the preprocess stage, a special instruction can be used to have runtests.pl
generate a repetetive sequence of bytes.

To insert a sequence of repeat bytes, use this syntax to make the `<string>`
get repeated `<number>` of times. The number has to be 1 or larger and the
string may contain `%HH` hexadecimal codes:

    %repeat[<number> x <string>]%

For example, to insert the word hello a 100 times:

    %repeat[100 x hello]%

## Conditional lines

Lines in the test file can be made to appear conditionally on a specific
feature (see the "features" section below) being set or not set. If the
specific feature is present, the following lines will be output, otherwise it
outputs nothing, until a following else or endif clause. Like this:

    %if brotli
    Accept-Encoding
    %endif

It can also check for the inversed condition, so if the feature us *not* set by
the use of an exclamation mark:

    %if !brotli
    Accept-Encoding: not-brotli
    %endif

You can also make an "else" clause to get output for the opposite condition,
like:

    %if brotli
    Accept-Encoding: brotli
    %else
    Accept-Encoding: nothing
    %endif

**Note** that there can be no nested conditions. You can only do one
conditional at a time and you can only check for a single feature in it.

# Variables

When the test is preprocessed, a range of "variables" in the test file will be
replaced by their content at that time.

Available substitute variables include:

- `%CLIENT6IP` - IPv6 address of the client running curl
- `%CLIENTIP` - IPv4 address of the client running curl
- `%CURL` - Path to the curl executable
- `%FILE_PWD` - Current directory, on windows prefixed with a slash
- `%FTP6PORT` - IPv6 port number of the FTP server
- `%FTPPORT` - Port number of the FTP server
- `%FTPSPORT` - Port number of the FTPS server
- `%FTPTIME2` - Timeout in seconds that should be just sufficient to receive a
  response from the test FTP server
- `%FTPTIME3` - Even longer than %FTPTIME2
- `%GOPHER6PORT` - IPv6 port number of the Gopher server
- `%GOPHERPORT` - Port number of the Gopher server
- `%GOPHERSPORT` - Port number of the Gophers server
- `%HOST6IP` - IPv6 address of the host running this test
- `%HOSTIP` - IPv4 address of the host running this test
- `%HTTP6PORT` - IPv6 port number of the HTTP server
- `%HTTPPORT` - Port number of the HTTP server
- `%HTTP2PORT` - Port number of the HTTP/2 server
- `%HTTPSPORT` - Port number of the HTTPS server
- `%HTTPSPROXYPORT` - Port number of the HTTPS-proxy
- `%HTTPTLS6PORT` - IPv6 port number of the HTTP TLS server
- `%HTTPTLSPORT` - Port number of the HTTP TLS server
- `%HTTPUNIXPATH` - Path to the Unix socket of the HTTP server
- `%IMAP6PORT` - IPv6 port number of the IMAP server
- `%IMAPPORT` - Port number of the IMAP server
- `%MQTTPORT` - Port number of the MQTT server
- `%TELNETPORT` - Port number of the telnet server
- `%NOLISTENPORT` - Port number where no service is listening
- `%POP36PORT` - IPv6 port number of the POP3 server
- `%POP3PORT` - Port number of the POP3 server
- `%POSIX_PWD` - Current directory somewhat mingw friendly
- `%PROXYPORT` - Port number of the HTTP proxy
- `%PWD` - Current directory
- `%RTSP6PORT` - IPv6 port number of the RTSP server
- `%RTSPPORT` - Port number of the RTSP server
- `%SMBPORT` - Port number of the SMB server
- `%SMBSPORT` - Port number of the SMBS server
- `%SMTP6PORT` - IPv6 port number of the SMTP server
- `%SMTPPORT` - Port number of the SMTP server
- `%SOCKSPORT` - Port number of the SOCKS4/5 server
- `%SRCDIR` - Full path to the source dir
- `%SSHPORT` - Port number of the SCP/SFTP server
- `%SSHSRVMD5` - MD5 of SSH server's public key
- `%SSHSRVSHA256` - SHA256 of SSH server's public key
- `%SSH_PWD` - Current directory friendly for the SSH server
- `%TESTNUMBER` - Number of the test case
- `%TFTP6PORT` - IPv6 port number of the TFTP server
- `%TFTPPORT` - Port number of the TFTP server
- `%USER` - Login ID of the user running the test
- `%VERSION` - the full version number of the tested curl

# `<testcase>`

Each test is always specified entirely within the testcase tag. Each test case
is split up in four main sections: `info`, `reply`, `client` and `verify`.

- **info** provides information about the test case

- **reply** is used for the server to know what to send as a reply for the
requests curl sends

- **client** defines how the client should behave

- **verify** defines how to verify that the data stored after a command has
been run ended up correctly

Each main section has a number of available subsections that can be specified,
that will be checked/used if specified.

## `<info>`

### `<keywords>`
A newline-separated list of keywords describing what this test case uses and
tests. Try to use already used keywords.  These keywords will be used for
statistical/informational purposes and for choosing or skipping classes of
tests.  "Keywords" must begin with an alphabetic character, "-", "[" or "{"
and may actually consist of multiple words separated by spaces which are
treated together as a single identifier.

When using curl built with Hyper, the keywords must include HTTP or HTTPS for
'hyper mode' to kick in and make line ending checks work for tests.
## `<reply>`

### `<data [nocheck="yes"] [sendzero="yes"] [base64="yes"] [hex="yes"] [nonewline="yes"]>`

data to be sent to the client on its request and later verified that it
arrived safely. Set `nocheck="yes"` to prevent the test script from verifying
the arrival of this data.

If the data contains `swsclose` anywhere within the start and end tag, and
this is a HTTP test, then the connection will be closed by the server after
this response is sent. If not, the connection will be kept persistent.

If the data contains `swsbounce` anywhere within the start and end tag, the
HTTP server will detect if this is a second request using the same test and
part number and will then increase the part number with one. This is useful
for auth tests and similar.

`sendzero=yes` means that the (FTP) server will "send" the data even if the
size is zero bytes. Used to verify curl's behavior on zero bytes transfers.

`base64=yes` means that the data provided in the test-file is a chunk of data
encoded with base64. It is the only way a test case can contain binary
data. (This attribute can in fact be used on any section, but it doesn't make
much sense for other sections than "data").

`hex=yes` means that the data is a sequence of hex pairs. It will get decoded
and used as "raw" data.

`nonewline=yes` means that the last byte (the trailing newline character)
should be cut off from the data before sending or comparing it.

For FTP file listings, the `<data>` section will be used *only* if you make
sure that there has been a CWD done first to a directory named `test-[num]`
where [num] is the test case number. Otherwise the ftp server can't know from
which test file to load the list content.

### `<dataNUM>`

Send back this contents instead of the <data> one. The num is set by:

 - The test number in the request line is >10000 and this is the remainder
   of [test case number]%10000.
 - The request was HTTP and included digest details, which adds 1000 to NUM
 - If a HTTP request is NTLM type-1, it adds 1001 to num
 - If a HTTP request is NTLM type-3, it adds 1002 to num
 - If a HTTP request is Basic and num is already >=1000, it adds 1 to num
 - If a HTTP request is Negotiate, num gets incremented by one for each
   request with Negotiate authorization header on the same test case.

Dynamically changing num in this way allows the test harness to be used to
test authentication negotiation where several different requests must be sent
to complete a transfer. The response to each request is found in its own data
section.  Validating the entire negotiation sequence can be done by specifying
a datacheck section.

### `<connect>`
The connect section is used instead of the 'data' for all CONNECT
requests. The remainder of the rules for the data section then apply but with
a connect prefix.

### `<socks>`
Address type and address details as logged by the SOCKS proxy.

### `<datacheck [mode="text"] [nonewline="yes"]>`
if the data is sent but this is what should be checked afterwards. If
`nonewline=yes` is set, runtests will cut off the trailing newline from the
data before comparing with the one actually received by the client.

Use the `mode="text"` attribute if the output is in text mode on platforms
that have a text/binary difference.

### `<datacheckNUM [nonewline="yes"] [mode="text"]>`
The contents of numbered datacheck sections are appended to the non-numbered
one.

### `<size>`
number to return on a ftp SIZE command (set to -1 to make this command fail)

### `<mdtm>`
what to send back if the client sends a (FTP) MDTM command, set to -1 to
have it return that the file doesn't exist

### `<postcmd>`
special purpose server-command to control its behavior *after* the
reply is sent
For HTTP/HTTPS, these are supported:

`wait [secs]` - Pause for the given time

### `<servercmd>`
Special-commands for the server.

The first line of this file will always be set to `Testnum [number]` by the
test script, to allow servers to read that to know what test the client is
about to issue.

#### For FTP/SMTP/POP/IMAP

- `REPLY [command] [return value] [response string]` - Changes how the server
  responds to the [command]. [response string] is evaluated as a perl string,
  so it can contain embedded \r\n, for example. There's a special [command]
  named "welcome" (without quotes) which is the string sent immediately on
  connect as a welcome.
- `REPLYLF` (like above but sends the response terminated with LF-only and not
   CRLF)
- `COUNT [command] [num]` - Do the `REPLY` change for `[command]` only `[num]`
  times and then go back to the built-in approach
- `DELAY [command] [secs]` - Delay responding to this command for the given
  time
- `RETRWEIRDO` - Enable the "weirdo" RETR case when multiple response lines
   appear at once when a file is transferred
- `RETRNOSIZE` - Make sure the RETR response doesn't contain the size of the
  file
- `NOSAVE` - Don't actually save what is received
- `SLOWDOWN` - Send FTP responses with 0.01 sec delay between each byte
- `PASVBADIP` - makes PASV send back an illegal IP in its 227 response
- `CAPA [capabilities]` - Enables support for and specifies a list of space
   separated capabilities to return to the client for the IMAP `CAPABILITY`,
   POP3 `CAPA` and SMTP `EHLO` commands
- `AUTH [mechanisms]` - Enables support for SASL authentication and specifies
   a list of space separated mechanisms for IMAP, POP3 and SMTP
- `STOR [msg]` respond with this instead of default after `STOR`

#### For HTTP/HTTPS

- `auth_required` if this is set and a POST/PUT is made without auth, the
  server will NOT wait for the full request body to get sent
- `idle` - do nothing after receiving the request, just "sit idle"
- `stream` - continuously send data to the client, never-ending
- `writedelay: [secs]` delay this amount between reply packets
- `skip: [num]` - instructs the server to ignore reading this many bytes from
  a PUT or POST request
- `rtp: part [num] channel [num] size [num]` - stream a fake RTP packet for
  the given part on a chosen channel with the given payload size
- `connection-monitor` - When used, this will log `[DISCONNECT]` to the
  `server.input` log when the connection is disconnected.
- `upgrade` - when an HTTP upgrade header is found, the server will upgrade to
  http2
- `swsclose` - instruct server to close connection after response
- `no-expect` - don't read the request body if Expect: is present

#### For TFTP
`writedelay: [secs]` delay this amount between reply packets (each packet
  being 512 bytes payload)

## `<client>`

### `<server>`
What server(s) this test case requires/uses. Available servers:

- `file`
- `ftp-ipv6`
- `ftp`
- `ftps`
- `gopher`
- `gophers`
- `http-ipv6`
- `http-proxy`
- `http-unix`
- `http/2`
- `http`
- `https`
- `httptls+srp-ipv6`
- `httptls+srp`
- `imap`
- `mqtt`
- `none`
- `pop3`
- `rtsp-ipv6`
- `rtsp`
- `scp`
- `sftp`
- `smtp`
- `socks4`
- `socks5`

Give only one per line.  This subsection is mandatory.

### `<features>`
A list of features that MUST be present in the client/library for this test to
be able to run. If a required feature is not present then the test will be
SKIPPED.

Alternatively a feature can be prefixed with an exclamation mark to indicate a
feature is NOT required. If the feature is present then the test will be
SKIPPED.

Features testable here are:

- `alt-svc`
- `c-ares`
- `cookies`
- `crypto`
- `debug`
- `DoH`
- `getrlimit`
- `GnuTLS`
- `GSS-API`
- `HSTS`
- `HTTP-auth`
- `http/2`
- `hyper`
- `idn`
- `ipv6`
- `Kerberos`
- `large_file`
- `ld_preload`
- `libssh2`
- `libssh`
- `libz`
- `manual`
- `Mime`
- `netrc`
- `NSS`
- `NTLM`
- `OpenSSL`
- `parsedate`
- `proxy`
- `PSL`
- `rustls`
- `Schannel`
- `sectransp`
- `shuffle-dns`
- `socks`
- `SPNEGO`
- `SSL`
- `SSLpinning`
- `SSPI`
- `threaded-resolver`
- `TLS-SRP`
- `TrackMemory`
- `typecheck`
- `Unicode`
- `unittest`
- `unix-sockets`
- `verbose-strings`
- `wakeup`
- `win32`
- `wolfssh`
- `wolfssl`

as well as each protocol that curl supports.  A protocol only needs to be
specified if it is different from the server (useful when the server
is `none`).

### `<killserver>`
Using the same syntax as in `<server>` but when mentioned here these servers
are explicitly KILLED when this test case is completed. Only use this if there
is no other alternatives. Using this of course requires subsequent tests to
restart servers.

### `<precheck>`
A command line that if set gets run by the test script before the test. If an
output is displayed by the command or if the return code is non-zero, the test
will be skipped and the (single-line) output will be displayed as reason for
not running the test.

### `<postcheck>`
A command line that if set gets run by the test script after the test. If
the command exists with a non-zero status code, the test will be considered
to have failed.

### `<tool>`
Name of tool to invoke instead of "curl". This tool must be built and exist
either in the libtest/ directory (if the tool name starts with 'lib') or in
the unit/ directory (if the tool name starts with 'unit').

### `<name>`
Brief test case description, shown when the test runs.

### `<setenv>`
    variable1=contents1
    variable2=contents2

Set the given environment variables to the specified value before the actual
command is run. They are cleared again after the command has been run.

### `<command [option="no-output/no-include/force-output/binary-trace"] [timeout="secs"][delay="secs"][type="perl/shell"]>`
Command line to run.

Note that the URL that gets passed to the server actually controls what data
that is returned. The last slash in the URL must be followed by a number. That
number (N) will be used by the test-server to load test case N and return the
data that is defined within the `<reply><data></data></reply>` section.

If there's no test number found above, the HTTP test server will use the
number following the last dot in the given hostname (made so that a CONNECT
can still pass on test number) so that "foo.bar.123" gets treated as test case
123. Alternatively, if an IPv6 address is provided to CONNECT, the last
hexadecimal group in the address will be used as the test number! For example
the address "[1234::ff]" would be treated as test case 255.

Set `type="perl"` to write the test case as a perl script. It implies that
there's no memory debugging and valgrind gets shut off for this test.

Set `type="shell"` to write the test case as a shell script. It implies that
there's no memory debugging and valgrind gets shut off for this test.

Set `option="no-output"` to prevent the test script to slap on the `--output`
argument that directs the output to a file. The `--output` is also not added
if the verify/stdout section is used.

Set `option="force-output"` to make use of `--output` even when the test is
otherwise written to verify stdout.

Set `option="no-include"` to prevent the test script to slap on the
`--include` argument.

Set `option="binary-trace"` to use `--trace` instead of `--trace-ascii` for
tracing.  Suitable for binary-oriented protocols such as MQTT.

Set `timeout="secs"` to override default server logs advisor read lock
timeout.  This timeout is used by the test harness, once that the command has
completed execution, to wait for the test server to write out server side log
files and remove the lock that advised not to read them. The "secs" parameter
is the not negative integer number of seconds for the timeout. This `timeout`
attribute is documented for completeness sake, but is deep test harness stuff
and only needed for very singular and specific test cases. Avoid using it.

Set `delay="secs"` to introduce a time delay once that the command has
completed execution and before the `<postcheck>` section runs. The "secs"
parameter is the not negative integer number of seconds for the delay. This
'delay' attribute is intended for very specific test cases, and normally not
needed.

### `<file name="log/filename" [nonewline="yes"]>`
This creates the named file with this content before the test case is run,
which is useful if the test case needs a file to act on.

If 'nonewline="yes"` is used, the created file will have the final newline
stripped off.

### `<stdin [nonewline="yes"]>`
Pass this given data on stdin to the tool.

If 'nonewline' is set, we will cut off the trailing newline of this given data
before comparing with the one actually received by the client

## `<verify>`
### `<errorcode>`
numerical error code curl is supposed to return. Specify a list of accepted
error codes by separating multiple numbers with comma. See test 237 for an
example.

### `<strip>`
One regex per line that is removed from the protocol dumps before the
comparison is made. This is very useful to remove dependencies on dynamically
changing protocol data such as port numbers or user-agent strings.

### `<strippart>`
One perl op per line that operates on the protocol dump. This is pretty
advanced. Example: `s/^EPRT .*/EPRT stripped/`.

### `<protocol [nonewline="yes"]>`

the protocol dump curl should transmit, if 'nonewline' is set, we will cut off
the trailing newline of this given data before comparing with the one actually
sent by the client The `<strip>` and `<strippart>` rules are applied before
comparisons are made.

### `<proxy [nonewline="yes"]>`

The protocol dump curl should transmit to a HTTP proxy (when the http-proxy
server is used), if 'nonewline' is set, we will cut off the trailing newline
of this given data before comparing with the one actually sent by the client
The `<strip>` and `<strippart>` rules are applied before comparisons are made.

### `<stderr [mode="text"] [nonewline="yes"]>`
This verifies that this data was passed to stderr.

Use the mode="text" attribute if the output is in text mode on platforms that
have a text/binary difference.

If 'nonewline' is set, we will cut off the trailing newline of this given data
before comparing with the one actually received by the client

### `<stdout [mode="text"] [nonewline="yes"]>`
This verifies that this data was passed to stdout.

Use the mode="text" attribute if the output is in text mode on platforms that
have a text/binary difference.

If 'nonewline' is set, we will cut off the trailing newline of this given data
before comparing with the one actually received by the client

### `<file name="log/filename" [mode="text"]>`
The file's contents must be identical to this after the test is complete.  Use
the mode="text" attribute if the output is in text mode on platforms that have
a text/binary difference.

### `<file1>`
1 to 4 can be appended to 'file' to compare more files.

### `<file2>`

### `<file3>`

### `<file4>`

### `<stripfile>`
One perl op per line that operates on the output file or stdout before being
compared with what is stored in the test file. This is pretty
advanced. Example: "s/^EPRT .*/EPRT stripped/"

### `<stripfile1>`
1 to 4 can be appended to 'stripfile' to strip the corresponding <fileN>
content

### `<stripfile2>`

### `<stripfile3>`

### `<stripfile4>`

### `<upload>`
the contents of the upload data curl should have sent

### `<valgrind>`
disable - disables the valgrind log check for this test
# The curl Test Suite

# Running

## Requires to run

  - perl (and a unix-style shell)
  - python (and a unix-style shell, for SMB and TELNET tests)
  - python-impacket (for SMB tests)
  - diff (when a test fails, a diff is shown)
  - stunnel (for HTTPS and FTPS tests)
  - OpenSSH or SunSSH (for SCP, SFTP and SOCKS4/5 tests)
  - nghttpx (for HTTP/2 tests)
  - nroff (for --manual tests)
  - An available `en_US.UTF-8` locale

### Installation of python-impacket

  The Python-based test servers support both recent Python 2 and 3.
  You can figure out your default Python interpreter with python -V

  Please install python-impacket in the correct Python environment.
  You can use pip or your OS' package manager to install 'impacket'.

  On Debian/Ubuntu the package names are:

  -  Python 2: 'python-impacket'
  -  Python 3: 'python3-impacket'

  On FreeBSD the package names are:

  -  Python 2: 'py27-impacket'
  -  Python 3: 'py37-impacket'

  On any system where pip is available:

  -  Python 2: 'pip2 install impacket'
  -  Python 3: 'pip3 install impacket'

  You may also need to manually install the Python package 'six'
  as that may be a missing requirement for impacket on Python 3.

### Port numbers used by test servers

  All test servers run on "random" port numbers. All tests should be written
  to use suitable variables instead of fixed port numbers so that test cases
  continue to work independent on what port numbers the test servers actually
  use.

  See [FILEFORMAT](FILEFORMAT.md) for the port number variables.

### Test servers

  The test suite runs stand-alone servers on random ports to which it makes
  requests. For SSL tests, it runs stunnel to handle encryption to the regular
  servers. For SSH, it runs a standard OpenSSH server. For SOCKS4/5 tests SSH
  is used to perform the SOCKS functionality and requires a SSH client and
  server.

  The listen port numbers for the test servers are picked randomly to allow
  users to run multiple test cases concurrently and to not collide with other
  existing services that might listen to ports on the machine.

  The HTTP server supports listening on a Unix domain socket, the default
  location is 'http.sock'.

### Run

  `./configure && make && make test`. This builds the test suite support code
  and invokes the 'runtests.pl' perl script to run all the tests. Edit the top
  variables of that script in case you have some specific needs, or run the
  script manually (after the support code has been built).

  The script breaks on the first test that doesn't do OK. Use `-a` to prevent
  the script from aborting on the first error. Run the script with `-v` for
  more verbose output. Use `-d` to run the test servers with debug output
  enabled as well. Specifying `-k` keeps all the log files generated by the
  test intact.

  Use `-s` for shorter output, or pass test numbers to run specific tests only
  (like `./runtests.pl 3 4` to test 3 and 4 only). It also supports test case
  ranges with 'to', as in `./runtests.pl 3 to 9` which runs the seven tests
  from 3 to 9. Any test numbers starting with ! are disabled, as are any test
  numbers found in the files `data/DISABLED` or `data/DISABLED.local` (one per
  line). The latter is meant for local temporary disables and will be ignored
  by git.

  Test cases mentioned in `DISABLED` can still be run if `-f` is provided.

  When `-s` is not present, each successful test will display on one line the
  test number and description and on the next line a set of flags, the test
  result, current test sequence, total number of tests to be run and an
  estimated amount of time to complete the test run. The flags consist of
  these letters describing what is checked in this test:

    s stdout
    d data
    u upload
    p protocol
    o output
    e exit code
    m memory
    v valgrind

### Shell startup scripts

  Tests which use the ssh test server, SCP/SFTP/SOCKS tests, might be badly
  influenced by the output of system wide or user specific shell startup
  scripts, .bashrc, .profile, /etc/csh.cshrc, .login, /etc/bashrc, etc. which
  output text messages or escape sequences on user login.  When these shell
  startup messages or escape sequences are output they might corrupt the
  expected stream of data which flows to the sftp-server or from the ssh
  client which can result in bad test behavior or even prevent the test
  server from running.

  If the test suite ssh or sftp server fails to start up and logs the message
  'Received message too long' then you are certainly suffering the unwanted
  output of a shell startup script.  Locate, cleanup or adjust the shell
  script.

### Memory test

  The test script will check that all allocated memory is freed properly IF
  curl has been built with the `CURLDEBUG` define set. The script will
  automatically detect if that is the case, and it will use the
  'memanalyze.pl' script to analyze the memory debugging output.

  Also, if you run tests on a machine where valgrind is found, the script will
  use valgrind to run the test with (unless you use `-n`) to further verify
  correctness.

  runtests.pl's `-t` option will enable torture testing mode, which runs each
  test many times and makes each different memory allocation fail on each
  successive run.  This tests the out of memory error handling code to ensure
  that memory leaks do not occur even in those situations. It can help to
  compile curl with `CPPFLAGS=-DMEMDEBUG_LOG_SYNC` when using this option, to
  ensure that the memory log file is properly written even if curl crashes.

### Debug

  If a test case fails, you can conveniently get the script to invoke the
  debugger (gdb) for you with the server running and the exact same command
  line parameters that failed. Just invoke `runtests.pl <test number> -g` and
  then just type 'run' in the debugger to perform the command through the
  debugger.

### Logs

  All logs are generated in the log/ subdirectory (it is emptied first in the
  runtests.pl script). They remain in there after a test run.

### Test input files

  All test cases are put in the `data/` subdirectory. Each test is stored in
  the file named according to the test number.

  See [FILEFORMAT.md](FILEFORMAT.md) for a description of the test case file
  format.

### Code coverage

  gcc provides a tool that can determine the code coverage figures for the
  test suite.  To use it, configure curl with `CFLAGS='-fprofile-arcs
  -ftest-coverage -g -O0'`.  Make sure you run the normal and torture tests to
  get more full coverage, i.e. do:

    make test
    make test-torture

  The graphical tool ggcov can be used to browse the source and create
  coverage reports on *NIX hosts:

    ggcov -r lib src

  The text mode tool gcov may also be used, but it doesn't handle object files
  in more than one directory very well.

### Remote testing

  The runtests.pl script provides some hooks to allow curl to be tested on a
  machine where perl can not be run.  The test framework in this case runs on
  a workstation where perl is available, while curl itself is run on a remote
  system using ssh or some other remote execution method.  See the comments at
  the beginning of runtests.pl for details.

## Test case numbering

  Test cases used to be numbered by category ranges, but the ranges filled
  up. Subsets of tests can now be selected by passing keywords to the
  runtests.pl script via the make `TFLAGS` variable.

  New tests are added by finding a free number in `tests/data/Makefile.inc`.

## Write tests

  Here's a quick description on writing test cases. We basically have three
  kinds of tests: the ones that test the curl tool, the ones that build small
  applications and test libcurl directly and the unit tests that test
  individual (possibly internal) functions.

### test data

  Each test has a master file that controls all the test data. What to read,
  what the protocol exchange should look like, what exit code to expect and
  what command line arguments to use etc.

  These files are `tests/data/test[num]` where `[num]` is just a unique
  identifier described above, and the XML-like file format of them is
  described in the separate [FILEFORMAT.md](FILEFORMAT.md) document.

### curl tests

  A test case that runs the curl tool and verifies that it gets the correct
  data, it sends the correct data, it uses the correct protocol primitives
  etc.

### libcurl tests

  The libcurl tests are identical to the curl ones, except that they use a
  specific and dedicated custom-built program to run instead of "curl". This
  tool is built from source code placed in `tests/libtest` and if you want to
  make a new libcurl test that is where you add your code.

### unit tests

  Unit tests are placed in `tests/unit`. There's a tests/unit/README
  describing the specific set of checks and macros that may be used when
  writing tests that verify behaviors of specific individual functions.

  The unit tests depend on curl being built with debug enabled.
# Unit tests

The goal is to add tests for *all* functions in libcurl. If functions are too
big and complicated, we should split them into smaller and testable ones.

## Build Unit Tests

`./configure --enable-debug` is required for the unit tests to build. To
enable unit tests, there will be a separate static libcurl built that will be
used exclusively for linking unit test programs. Just build everything as
normal, and then you can run the unit test cases as well.

## Run Unit Tests

Unit tests are run as part of the regular test suite. If you have built
everything to run unit tests, to can do 'make test' at the root level. Or you
can `cd tests` and `make` and then invoke individual unit tests with
`./runtests.pl NNNN` where `NNNN` is the specific test number.

## Debug Unit Tests

If a specific test fails you will get told. The test case then has output left
in the log/ subdirectory, but most importantly you can re-run the test again
using gdb by doing `./runtests.pl -g NNNN`. That is, add a `-g` to make it
start up gdb and run the same case using that.

## Write Unit Tests

We put tests that focus on an area or a specific function into a single C
source file. The source file should be named 'unitNNNN.c' where NNNN is a
previously unused number.

Add your test to `tests/unit/Makefile.inc` (if it is a unit test).  Add your
test data file name to `tests/data/Makefile.inc`

You also need a separate file called `tests/data/testNNNN` (using the same
number) that describes your test case. See the test1300 file for inspiration
and the `tests/FILEFORMAT.md` documentation.

For the actual C file, here's a very simple example:
~~~c
#include "curlcheck.h"

#include "a libcurl header.h" /* from the lib dir */

static CURLcode unit_setup( void )
{
  /* whatever you want done first */
  return CURLE_OK;
}

static void unit_stop( void )
{
  /* done before shutting down and exiting */
}

UNITTEST_START

  /* here you start doing things and checking that the results are good */

  fail_unless( size == 0 , "initial size should be zero" );
  fail_if( head == NULL , "head should not be initiated to NULL" );

  /* you end the test code like this: */

UNITTEST_STOP
How to contribute to curl
=========================

Join the community
------------------

 1. Click 'watch' on the GitHub repo

 2. Subscribe to the suitable [mailing lists](https://curl.se/mail/)

Read [CONTRIBUTE](../docs/CONTRIBUTE.md)
---------------------------------------

Send your suggestions using one of these methods:
-------------------------------------------------

 1. in a mail to the mailing list

 2. as a [pull request](https://github.com/curl/curl/pulls)

 3. as an [issue](https://github.com/curl/curl/issues)

/ The curl team!
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

<!-- Only file bugs here! Ask questions on the mailing lists https://curl.se/mail/

     SECURITY RELATED? Post it here: https://hackerone.com/curl

     There are collections of known issues to be aware of:
     https://curl.se/docs/knownbugs.html
     https://curl.se/docs/todo.html       -->

### I did this

### I expected the following

### curl/libcurl version

[curl -V output]

### operating system

<!-- On Unix please post the output of "uname -a" -->
# Building curl with Visual C++

 This document describes how to compile, build and install curl and libcurl
 from sources using the Visual C++ build tool. To build with VC++, you will of
 course have to first install VC++. The minimum required version of VC is 6
 (part of Visual Studio 6). However using a more recent version is strongly
 recommended.

 VC++ is also part of the Windows Platform SDK. You do not have to install the
 full Visual Studio or Visual C++ if all you want is to build curl.

 The latest Platform SDK can be downloaded freely from [Windows SDK and
 emulator
 archive](https://developer.microsoft.com/en-us/windows/downloads/sdk-archive)

## Prerequisites

 If you wish to support zlib, openssl, c-ares, ssh2, you will have to download
 them separately and copy them to the deps directory as shown below:

    somedirectory\
     |_curl-src
     | |_winbuild
     |
     |_deps
       |_ lib
       |_ include
       |_ bin

 It is also possible to create the deps directory in some other random places
 and tell the Makefile its location using the WITH_DEVEL option.

## Building straight from git

 When you check out code git and build it, as opposed from a released source
 code archive, you need to first run the `buildconf.bat` batch file (present
 in the source code root directory) to set things up.

## Open a command prompt

Open a Visual Studio Command prompt:

 Using the **'Developer Command Prompt for VS [version]'** menu entry: where
 [version} is the Visual Studio version. The developer prompt at default uses
 the x86 mode. It is required to call `Vcvarsall.bat` to setup the prompt for
 the machine type you want. This type of command prompt may not exist in all
 Visual Studio versions.

 See also: [Developer Command Prompt for Visual
 Studio](https://docs.microsoft.com/en-us/dotnet/framework/tools/developer-command-prompt-for-vs)
 and [How to: Enable a 64-Bit, x64 hosted MSVC toolset on the command
 line](https://docs.microsoft.com/en-us/cpp/build/how-to-enable-a-64-bit-visual-cpp-toolset-on-the-command-line)

 Using the **'VS [version] [platform] [type] Command Prompt'** menu entry:
 where [version] is the Visual Studio version, [platform] is e.g. x64 and
 [type] Native of Cross platform build.  This type of command prompt may not
 exist in all Visual Studio versions.

 See also: [Set the Path and Environment Variables for Command-Line Builds](https://msdn.microsoft.com/en-us/library/f2ccy3wt.aspx)

## Build in the console

 Once you are in the console, go to the winbuild directory in the Curl
 sources:

    cd curl-src\winbuild

 Then you can call `nmake /f Makefile.vc` with the desired options (see
 below). The builds will be in the top src directory, `builds\` directory, in
 a directory named using the options given to the nmake call.

    nmake /f Makefile.vc mode=<static or dll> <options>

where `<options>` is one or many of:

 - `VC=<num>`                    - VC version. 6 or later.
 - `WITH_DEVEL=<path>`           - Paths for the development files (SSL, zlib, etc.)
                                   Defaults to sibbling directory deps: ../deps
                                   Libraries can be fetched at https://windows.php.net/downloads/php-sdk/deps/
                                   Uncompress them into the deps folder.
 - `WITH_SSL=<dll/static>`       - Enable OpenSSL support, DLL or static
 - `WITH_NGHTTP2=<dll/static>`   - Enable HTTP/2 support, DLL or static
 - `WITH_MBEDTLS=<dll/static>`   - Enable mbedTLS support, DLL or static
 - `WITH_CARES=<dll/static>`     - Enable c-ares support, DLL or static
 - `WITH_ZLIB=<dll/static>`      - Enable zlib support, DLL or static
 - `WITH_SSH2=<dll/static>`      - Enable libSSH2 support, DLL or static
 - `WITH_PREFIX=<dir>`           - Where to install the build
 - `ENABLE_SSPI=<yes/no>`        - Enable SSPI support, defaults to yes
 - `ENABLE_IPV6=<yes/no>`        - Enable IPv6, defaults to yes
 - `ENABLE_IDN=<yes or no>`      - Enable use of Windows IDN APIs, defaults to yes
                                   Requires Windows Vista or later
 - `ENABLE_SCHANNEL=<yes/no>`    - Enable native Windows SSL support, defaults
                                   to yes if SSPI and no other SSL library
 - `ENABLE_OPENSSL_AUTO_LOAD_CONFIG=<yes/no>`
                                 - Enable loading OpenSSL configuration
                                   automatically, defaults to yes
 - `ENABLE_UNICODE=<yes/no>`     - Enable UNICODE support, defaults to no
 - `GEN_PDB=<yes/no>`            - Generate External Program Database
                                   (debug symbols for release build)
 - `DEBUG=<yes/no>`              - Debug builds
 - `MACHINE=<x86/x64>`           - Target architecture (default is x86)
 - `CARES_PATH=<path>`           - Custom path for c-ares
 - `MBEDTLS_PATH=<path>`         - Custom path for mbedTLS
 - `NGHTTP2_PATH=<path>`         - Custom path for nghttp2
 - `SSH2_PATH=<path>`            - Custom path for libSSH2
 - `SSL_PATH=<path>`             - Custom path for OpenSSL
 - `ZLIB_PATH=<path>`            - Custom path for zlib

## Static linking of Microsoft's C RunTime (CRT):

 If you are using mode=static nmake will create and link to the static build
 of libcurl but *not* the static CRT. If you must you can force nmake to link
 in the static CRT by passing RTLIBCFG=static. Typically you shouldn't use
 that option, and nmake will default to the DLL CRT. RTLIBCFG is rarely used
 and therefore rarely tested. When passing RTLIBCFG for a configuration that
 was already built but not with that option, or if the option was specified
 differently, you must destroy the build directory containing the
 configuration so that nmake can build it from scratch.

## Building your own application with a static libcurl

 When building an application that uses the static libcurl library on Windows,
 you must define CURL_STATICLIB. Otherwise the linker will look for dynamic
 import symbols.

## Legacy Windows and SSL

 When you build curl using the build files in this directory the default SSL
 backend will be Schannel (Windows SSPI), the native SSL library that comes
 with the Windows OS. Schannel in Windows <= XP is not able to connect to
 servers that no longer support the legacy handshakes and algorithms used by
 those versions. If you will be using curl in one of those earlier versions of
 Windows you should choose another SSL backend like OpenSSL.
# checksrc

This is the tool we use within the curl project to scan C source code and
check that it adheres to our [Source Code Style guide](CODE_STYLE.md).

## Usage

    checksrc.pl [options] [file1] [file2] ...

## Command line options

`-W[file]` skip that file and excludes it from being checked. Helpful
when, for example, one of the files is generated.

`-D[dir]` directory name to prepend to file names when accessing them.

`-h` shows the help output, that also lists all recognized warnings

## What does checksrc warn for?

checksrc does not check and verify the code against the entire style guide.
The script is an effort to detect the most common mistakes and syntax mistakes
that contributors make before they get accustomed to our code style. Heck,
many of us regulars do the mistakes too and this script helps us keep the code
in shape.

    checksrc.pl -h

Lists how to use the script and it lists all existing warnings it has and
problems it detects. At the time of this writing, the existing checksrc
warnings are:

- `ASSIGNWITHINCONDITION`: Assignment within a conditional expression. The
  code style mandates the assignment to be done outside of it.

- `ASTERISKNOSPACE`: A pointer was declared like `char* name` instead of the
   more appropriate `char *name` style. The asterisk should sit next to the
   name.

- `ASTERISKSPACE`: A pointer was declared like `char * name` instead of the
   more appropriate `char *name` style. The asterisk should sit right next to
   the name without a space in between.

- `BADCOMMAND`: There's a bad `!checksrc!` instruction in the code. See the
   **Ignore certain warnings** section below for details.

- `BANNEDFUNC`: A banned function was used. The functions sprintf, vsprintf,
   strcat, strncat, gets are **never** allowed in curl source code.

- `BRACEELSE`: '} else' on the same line. The else is supposed to be on the
   following line.

- `BRACEPOS`: wrong position for an open brace (`{`).

- `BRACEWHILE`: more than once space between end brace and while keyword

- `COMMANOSPACE`: a comma without following space

- `COPYRIGHT`: the file is missing a copyright statement!

- `CPPCOMMENTS`: `//` comment detected, that is not C89 compliant

- `DOBRACE`: only use one space after do before open brace

- `EMPTYLINEBRACE`: found empty line before open brace

- `EQUALSNOSPACE`: no space after `=` sign

- `EQUALSNULL`: comparison with `== NULL` used in if/while. We use `!var`.

- `EXCLAMATIONSPACE`: space found after exclamations mark

- `FOPENMODE`: `fopen()` needs a macro for the mode string, use it

- `INDENTATION`: detected a wrong start column for code. Note that this
   warning only checks some specific places and will certainly miss many bad
   indentations.

- `LONGLINE`: A line is longer than 79 columns.

- `MULTISPACE`: Multiple spaces were found where only one should be used.

- `NOSPACEEQUALS`: An equals sign was found without preceding space. We prefer
  `a = 2` and *not* `a=2`.

- `NOTEQUALSZERO`: check found using `!= 0`. We use plain `if(var)`.

- `ONELINECONDITION`: do not put the conditional block on the same line as `if()`

- `OPENCOMMENT`: File ended with a comment (`/*`) still "open".

- `PARENBRACE`: `){` was used without sufficient space in between.

- `RETURNNOSPACE`: `return` was used without space between the keyword and the
   following value.

- `SEMINOSPACE`: There was no space (or newline) following a semicolon.

- `SIZEOFNOPAREN`: Found use of sizeof without parentheses. We prefer
  `sizeof(int)` style.

- `SNPRINTF` - Found use of `snprintf()`. Since we use an internal replacement
   with a different return code etc, we prefer `msnprintf()`.

- `SPACEAFTERPAREN`: there was a space after open parenthesis, `( text`.

- `SPACEBEFORECLOSE`: there was a space before a close parenthesis, `text )`.

- `SPACEBEFORECOMMA`: there was a space before a comma, `one , two`.

- `SPACEBEFOREPAREN`: there was a space before an open parenthesis, `if (`,
   where one was not expected

- `SPACESEMICOLON`: there was a space before semicolon, ` ;`.

- `TABS`: TAB characters are not allowed!

- `TRAILINGSPACE`: Trailing whitespace on the line

- `TYPEDEFSTRUCT`: we frown upon (most) typedefed structs

- `UNUSEDIGNORE`: a checksrc inlined warning ignore was asked for but not used,
   that is an ignore that should be removed or changed to get used.

### Extended warnings

Some warnings are quite computationally expensive to perform, so they are
turned off by default. To enable these warnings, place a `.checksrc` file in
the directory where they should be activated with commands to enable the
warnings you are interested in. The format of the file is to enable one
warning per line like so: `enable <EXTENDEDWARNING>`

Currently these are the extended warnings which can be enabled:

- `COPYRIGHTYEAR`: the current changeset has not updated the copyright year in
   the source file

- `STRERROR`: use of banned function strerror()

## Ignore certain warnings

Due to the nature of the source code and the flaws of the checksrc tool, there
is sometimes a need to ignore specific warnings. checksrc allows a few
different ways to do this.

### Inline ignore

You can control what to ignore within a specific source file by providing
instructions to checksrc in the source code itself. You need a magic marker
that is `!checksrc!` followed by the instruction. The instruction can ask to
ignore a specific warning N number of times or you ignore all of them until
you mark the end of the ignored section.

Inline ignores are only done for that single specific source code file.

Example

    /* !checksrc! disable LONGLINE all */

This will ignore the warning for overly long lines until it is re-enabled with:

    /* !checksrc! enable LONGLINE */

If the enabling is not performed before the end of the file, it will be enabled
automatically for the next file.

You can also opt to ignore just N violations so that if you have a single long
line you just cannot shorten and is agreed to be fine anyway:

    /* !checksrc! disable LONGLINE 1 */

... and the warning for long lines will be enabled again automatically after
it has ignored that single warning. The number `1` can of course be changed to
any other integer number. It can be used to make sure only the exact intended
instances are ignored and nothing extra.

### Directory wide ignore patterns

This is a method we have transitioned away from. Use inline ignores as far as
possible.

Make a `checksrc.skip` file in the directory of the source code with the
false positive, and include the full offending line into this file.
How curl Became Like This
=========================

Towards the end of 1996, Daniel Stenberg was spending time writing an IRC bot
for an Amiga related channel on EFnet. He then came up with the idea to make
currency-exchange calculations available to Internet Relay Chat (IRC)
users. All the necessary data were published on the Web; he just needed to
automate their retrieval.

1996
----

On November 11, 1996 the Brazilian developer Rafael Sagula wrote and released
HttpGet version 0.1.

Daniel extended this existing command-line open-source tool. After a few minor
adjustments, it did just what he needed. The first release with Daniel's
additions was 0.2, released on December 17, 1996. Daniel quickly became the
new maintainer of the project.

1997
----

HttpGet 0.3 was released in January 1997 and now it accepted HTTP URLs on the
command line.

HttpGet 1.0 was released on April 8th 1997 with brand new HTTP proxy support.

We soon found and fixed support for getting currencies over GOPHER. Once FTP
download support was added, the name of the project was changed and urlget 2.0
was released in August 1997. The http-only days were already passed.

Version 2.2 was released on August 14 1997 and introduced support to build for
and run on Windows and Solaris.

November 24 1997: Version 3.1 added FTP upload support.

Version 3.5 added support for HTTP POST.

1998
----

February 4: urlget 3.10

February 9: urlget 3.11

March 14: urlget 3.12 added proxy authentication.

The project slowly grew bigger. With upload capabilities, the name was once
again misleading and a second name change was made. On March 20, 1998 curl 4
was released. (The version numbering from the previous names was kept.)

(Unrelated to this project a company called Curl Corporation registered a US
trademark on the name "CURL" on May 18 1998. That company had then already
registered the curl.com domain back in November of the previous year. All this
was revealed to us much later.)

SSL support was added, powered by the SSLeay library.

August: first announcement of curl on freshmeat.net.

October: with the curl 4.9 release and the introduction of cookie support,
curl was no longer released under the GPL license. Now we are at 4000 lines of
code, we switched over to the MPL license to restrict the effects of
"copyleft".

November: configure script and reported successful compiles on several
major operating systems. The never-quite-understood -F option was added and
curl could now simulate quite a lot of a browser. TELNET support was added.

Curl 5 was released in December 1998 and introduced the first ever curl man
page. People started making Linux RPM packages out of it.

1999
----

January: DICT support added.

OpenSSL took over and SSLeay was abandoned.

May: first Debian package.

August: LDAP:// and FILE:// support added. The curl website gets 1300 visits
weekly. Moved site to curl.haxx.nu.

September: Released curl 6.0. 15000 lines of code.

December 28: added the project on Sourceforge and started using its services
for managing the project.

2000
----

Spring: major internal overhaul to provide a suitable library interface.
The first non-beta release was named 7.1 and arrived in August. This offered
the easy interface and turned out to be the beginning of actually getting
other software and programs to be based on and powered by libcurl. Almost
20000 lines of code.

June: the curl site moves to "curl.haxx.se"

August, the curl website gets 4000 visits weekly.

The PHP guys adopted libcurl already the same month, when the first ever third
party libcurl binding showed up. CURL has been a supported module in PHP since
the release of PHP 4.0.2. This would soon get followers. More than 16
different bindings exist at the time of this writing.

September: kerberos4 support was added.

November: started the work on a test suite for curl. It was later re-written
from scratch again. The libcurl major SONAME number was set to 1.

2001
----

January: Daniel released curl 7.5.2 under a new license again: MIT (or
MPL). The MIT license is extremely liberal and can be combined with GPL
in other projects. This would finally put an end to the "complaints" from
people involved in GPLed projects that previously were prohibited from using
libcurl while it was released under MPL only. (Due to the fact that MPL is
deemed "GPL incompatible".)

March 22: curl supports HTTP 1.1 starting with the release of 7.7. This
also introduced libcurl's ability to do persistent connections. 24000 lines of
code. The libcurl major SONAME number was bumped to 2 due to this overhaul.
The first experimental ftps:// support was added.

August: The curl website gets 8000 visits weekly. Curl Corporation contacted
Daniel to discuss "the name issue". After Daniel's reply, they have never
since got back in touch again.

September: libcurl 7.9 introduces cookie jar and curl_formadd(). During the
forthcoming 7.9.x releases, we introduced the multi interface slowly and
without many whistles.

September 25: curl (7.7.2) is bundled in Mac OS X (10.1) for the first time. It was
already becoming more and more of a standard utility of Linux distributions
and a regular in the BSD ports collections.

2002
----

June: the curl website gets 13000 visits weekly. curl and libcurl is
35000 lines of code. Reported successful compiles on more than 40 combinations
of CPUs and operating systems.

To estimate number of users of the curl tool or libcurl library is next to
impossible. Around 5000 downloaded packages each week from the main site gives
a hint, but the packages are mirrored extensively, bundled with numerous OS
distributions and otherwise retrieved as part of other software.

October 1: with the release of curl 7.10 it is released under the MIT license
only.

Starting with 7.10, curl verifies SSL server certificates by default.

2003
----

January: Started working on the distributed curl tests. The autobuilds.

February: the curl site averages at 20000 visits weekly. At any given moment,
there's an average of 3 people browsing the website.

Multiple new authentication schemes are supported: Digest (May), NTLM (June)
and Negotiate (June).

November: curl 7.10.8 is released. 45000 lines of code. ~55000 unique visitors
to the website. Five official web mirrors.

December: full-fledged SSL for FTP is supported.

2004
----

January: curl 7.11.0 introduced large file support.

June: curl 7.12.0 introduced IDN support. 10 official web mirrors.

This release bumped the major SONAME to 3 due to the removal of the
curl_formparse() function

August: Curl and libcurl 7.12.1

    Public curl release number:                82
    Releases counted from the very beginning: 109
    Available command line options:            96
    Available curl_easy_setopt() options:     120
    Number of public functions in libcurl:     36
    Amount of public website mirrors:         12
    Number of known libcurl bindings:          26

2005
----

April: GnuTLS can now optionally be used for the secure layer when curl is
built.

April: Added the multi_socket() API

September: TFTP support was added.

More than 100,000 unique visitors of the curl website. 25 mirrors.

December: security vulnerability: libcurl URL Buffer Overflow

2006
----

January: We dropped support for Gopher. We found bugs in the implementation
that turned out to have been introduced years ago, so with the conclusion that
nobody had found out in all this time we removed it instead of fixing it.

March: security vulnerability: libcurl TFTP Packet Buffer Overflow

September: The major SONAME number for libcurl was bumped to 4 due to the
removal of ftp third party transfer support.

November: Added SCP and SFTP support

2007
----

February: Added support for the Mozilla NSS library to do the SSL/TLS stuff

July: security vulnerability: libcurl GnuTLS insufficient cert verification

2008
----

November:

    Command line options:         128
    curl_easy_setopt() options:   158
    Public functions in libcurl:   58
    Known libcurl bindings:        37
    Contributors:                 683

 145,000 unique visitors. >100 GB downloaded.

2009
----

March: security vulnerability: libcurl Arbitrary File Access

April: added CMake support

August: security vulnerability: libcurl embedded zero in cert name

December: Added support for IMAP, POP3 and SMTP

2010
----

January: Added support for RTSP

February: security vulnerability: libcurl data callback excessive length

March: The project switched over to use git (hosted by GitHub) instead of CVS
for source code control

May: Added support for RTMP

Added support for PolarSSL to do the SSL/TLS stuff

August:

    Public curl releases:         117
    Command line options:         138
    curl_easy_setopt() options:   180
    Public functions in libcurl:   58
    Known libcurl bindings:        39
    Contributors:                 808

 Gopher support added (re-added actually, see January 2006)

2011
----

February: added support for the axTLS backend

April: added the cyassl backend (later renamed to WolfSSL)

2012
----

 July: Added support for Schannel (native Windows TLS backend) and Darwin SSL
 (Native Mac OS X and iOS TLS backend).

 Supports metalink

 October: SSH-agent support.

2013
----

 February: Cleaned up internals to always uses the "multi" non-blocking
 approach internally and only expose the blocking API with a wrapper.

 September: First small steps on supporting HTTP/2 with nghttp2.

 October: Removed krb4 support.

 December: Happy eyeballs.

2014
----

 March: first real release supporting HTTP/2

 September: Website had 245,000 unique visitors and served 236GB data

 SMB and SMBS support

2015
----

 June: support for multiplexing with HTTP/2

 August: support for HTTP/2 server push

 December: Public Suffix List

2016
----

 January: the curl tool defaults to HTTP/2 for HTTPS URLs

 December: curl 7.52.0 introduced support for HTTPS-proxy!

 First TLS 1.3 support

2017
----

 July: OSS-Fuzz started fuzzing libcurl

 September: Added Multi-SSL support

 The website serves 3100 GB/month

    Public curl releases:         169
    Command line options:         211
    curl_easy_setopt() options:   249
    Public functions in libcurl:  74
    Contributors:                 1609

 October: SSLKEYLOGFILE support, new MIME API

 October: Daniel received the Polhem Prize for his work on curl

 November: brotli

2018
----

 January: new SSH backend powered by libssh

 March: starting with the 1803 release of Windows 10, curl is shipped bundled
 with Microsoft's operating system.

 July: curl shows headers using bold type face

 October: added DNS-over-HTTPS (DoH) and the URL API

 MesaLink is a new supported TLS backend

 libcurl now does HTTP/2 (and multiplexing) by default on HTTPS URLs

 curl and libcurl are installed in an estimated 5 *billion* instances
 world-wide.

 October 31: Curl and libcurl 7.62.0

    Public curl releases:         177
    Command line options:         219
    curl_easy_setopt() options:   261
    Public functions in libcurl:  80
    Contributors:                 1808

 December: removed axTLS support

2019
----

 March: added experimental alt-svc support

 August: the first HTTP/3 requests with curl.

 September: 7.66.0 is released and the tool offers parallel downloads

2020
----

 curl and libcurl are installed in an estimated 10 *billion* instances
 world-wide.

 January: added BearSSL support

 March: removed support for PolarSSL, added wolfSSH support

 April: experimental MQTT support

 August: zstd support

 November: the website moves to curl.se. The website serves 10TB data monthly.

 December: alt-svc support

2021
----

 February 3: curl 7.75.0 ships with support for Hyper is a HTTP backend

 March 31: curl 7.76.0 ships with support for rustls

 July: HSTS is supported
# Alt-Svc

curl features support for the Alt-Svc: HTTP header.

## Enable Alt-Svc in build

`./configure --enable-alt-svc`

(enabled by default since 7.73.0)

## Standard

[RFC 7838](https://datatracker.ietf.org/doc/html/rfc7838)

# Alt-Svc cache file format

This a text based file with one line per entry and each line consists of nine
space separated fields.

## Example

    h2 quic.tech 8443 h3-22 quic.tech 8443 "20190808 06:18:37" 0 0

## Fields

1. The ALPN id for the source origin
2. The host name for the source origin
3. The port number for the source origin
4. The ALPN id for the destination host
5. The host name for the destination host
6. The host number for the destination host
7. The expiration date and time of this entry within double quotes. The date format is "YYYYMMDD HH:MM:SS" and the time zone is GMT.
8. Boolean (1 or 0) if "persist" was set for this entry
9. Integer priority value (not currently used)

# TODO

- handle multiple response headers, when one of them says `clear` (should
  override them all)
- using `Age:` value for caching age as per spec
- `CURLALTSVC_IMMEDIATELY` support
# The curl bug bounty

The curl project runs a bug bounty program in association with
[HackerOne](https://www.hackerone.com) and the [Internet Bug
Bounty](https://internetbugbounty.org).

# How does it work?

Start out by posting your suspected security vulnerability directly to [curl's
HackerOne program](https://hackerone.com/curl).

After you have reported a security issue, it has been deemed credible, and a
patch and advisory has been made public, you may be eligible for a bounty from
this program.

See all details at [https://hackerone.com/curl](https://hackerone.com/curl)

This bounty is relying on funds from sponsors. If you use curl professionally,
consider help funding this! See
[https://opencollective.com/curl](https://opencollective.com/curl) for
details.

# What are the reward amounts?

The curl project offers monetary compensation for reported and published
security vulnerabilities. The amount of money that is rewarded depends on how
serious the flaw is determined to be.

We offer reward money *up to* a certain amount per severity. The curl security
team determines the severity of each reported flaw on a case by case basis and
the exact amount rewarded to the reporter is then decided.

Check out the current award amounts at [https://hackerone.com/curl](https://hackerone.com/curl)

# Who is eligible for a reward?

Everyone and anyone who reports a security problem in a released curl version
that has not already been reported can ask for a bounty.

Vulnerabilities in features that are off by default and documented as
experimental are not eligible for a reward.

The vulnerability has to be fixed and publicly announced (by the curl project)
before a bug bounty will be considered.

Bounties need to be requested within twelve months from the publication of the
vulnerability.

# Product vulnerabilities only

This bug bounty only concerns the curl and libcurl products and thus their
respective source codes - when running on existing hardware. It does not
include documentation, websites, or other infrastructure.

The curl security team is the sole arbiter if a reported flaw is subject to a
bounty or not.

# How are vulnerabilities graded?

The grading of each reported vulnerability that makes a reward claim will be
performed by the curl security team. The grading will be based on the CVSS
(Common Vulnerability Scoring System) 3.0.

# How are reward amounts determined?

The curl security team first gives the vulnerability a score, as mentioned
above, and based on that level we set an amount depending on the specifics of
the individual case. Other sponsors of the program might also get involved and
can raise the amounts depending on the particular issue.

# What happens if the bounty fund is drained?

The bounty fund depends on sponsors. If we pay out more bounties than we add,
the fund will eventually drain. If that end up happening, we will simply not
be able to pay out as high bounties as we would like and hope that we can
convince new sponsors to help us top up the fund again.

# Regarding taxes, etc. on the bounties

In the event that the individual receiving a curl bug bounty needs to pay
taxes on the reward money, the responsibility lies with the receiver. The
curl project or its security team never actually receive any of this money,
hold the money, or pay out the money.
SSL Certificate Verification
============================

SSL is TLS
----------

SSL is the old name. It is called TLS these days.


Native SSL
----------

If libcurl was built with Schannel or Secure Transport support (the native SSL
libraries included in Windows and Mac OS X), then this does not apply to
you. Scroll down for details on how the OS-native engines handle SSL
certificates. If you are not sure, then run "curl -V" and read the results. If
the version string says `Schannel` in it, then it was built with Schannel
support.

It is about trust
-----------------

This system is about trust. In your local CA certificate store you have certs
from *trusted* Certificate Authorities that you then can use to verify that the
server certificates you see are valid. they are signed by one of the CAs you
trust.

Which CAs do you trust? You can decide to trust the same set of companies your
operating system trusts, or the set one of the known browsers trust. That is
basically trust via someone else you trust. You should just be aware that
modern operating systems and browsers are setup to trust *hundreds* of
companies and recent years several such CAs have been found untrustworthy.

Certificate Verification
------------------------

libcurl performs peer SSL certificate verification by default. This is done
by using a CA certificate store that the SSL library can use to make sure the
peer's server certificate is valid.

If you communicate with HTTPS, FTPS or other TLS-using servers using
certificates that are signed by CAs present in the store, you can be sure
that the remote server really is the one it claims to be.

If the remote server uses a self-signed certificate, if you do not install a CA
cert store, if the server uses a certificate signed by a CA that is not
included in the store you use or if the remote host is an impostor
impersonating your favorite site, and you want to transfer files from this
server, do one of the following:

 1. Tell libcurl to *not* verify the peer. With libcurl you disable this with
    `curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, FALSE);`

    With the curl command line tool, you disable this with -k/--insecure.

 2. Get a CA certificate that can verify the remote server and use the proper
    option to point out this CA cert for verification when connecting. For
    libcurl hackers: `curl_easy_setopt(curl, CURLOPT_CAINFO, cacert);`

    With the curl command line tool: --cacert [file]

 3. Add the CA cert for your server to the existing default CA certificate
    store. The default CA certificate store can be changed at compile time with
    the following configure options:

    --with-ca-bundle=FILE: use the specified file as CA certificate store. CA
    certificates need to be concatenated in PEM format into this file.

    --with-ca-path=PATH: use the specified path as CA certificate store. CA
    certificates need to be stored as individual PEM files in this directory.
    You may need to run c_rehash after adding files there.

    If neither of the two options is specified, configure will try to auto-detect
    a setting. It's also possible to explicitly not hardcode any default store
    but rely on the built in default the crypto library may provide instead.
    You can achieve that by passing both --without-ca-bundle and
    --without-ca-path to the configure script.

    If you use Internet Explorer, this is one way to get extract the CA cert
    for a particular server:

     - View the certificate by double-clicking the padlock
     - Find out where the CA certificate is kept (Certificate>
       Authority Information Access>URL)
     - Get a copy of the crt file using curl
     - Convert it from crt to PEM using the openssl tool:
       openssl x509 -inform DES -in yourdownloaded.crt \
       -out outcert.pem -text
     - Add the 'outcert.pem' to the CA certificate store or use it stand-alone
       as described below.

    If you use the 'openssl' tool, this is one way to get extract the CA cert
    for a particular server:

     - `openssl s_client -showcerts -servername server -connect server:443 > cacert.pem`
     - type "quit", followed by the "ENTER" key
     - The certificate will have "BEGIN CERTIFICATE" and "END CERTIFICATE"
       markers.
     - If you want to see the data in the certificate, you can do: "openssl
       x509 -inform PEM -in certfile -text -out certdata" where certfile is
       the cert you extracted from logfile. Look in certdata.
     - If you want to trust the certificate, you can add it to your CA
       certificate store or use it stand-alone as described. Just remember that
       the security is no better than the way you obtained the certificate.

 4. If you are using the curl command line tool, you can specify your own CA
    cert file by setting the environment variable `CURL_CA_BUNDLE` to the path
    of your choice.

    If you are using the curl command line tool on Windows, curl will search
    for a CA cert file named "curl-ca-bundle.crt" in these directories and in
    this order:
      1. application's directory
      2. current working directory
      3. Windows System directory (e.g. C:\windows\system32)
      4. Windows Directory (e.g. C:\windows)
      5. all directories along %PATH%

 5. Get a better/different/newer CA cert bundle! One option is to extract the
    one a recent Firefox browser uses by running 'make ca-bundle' in the curl
    build tree root, or possibly download a version that was generated this
    way for you: [CA Extract](https://curl.se/docs/caextract.html)

Neglecting to use one of the above methods when dealing with a server using a
certificate that is not signed by one of the certificates in the installed CA
certificate store, will cause SSL to report an error ("certificate verify
failed") during the handshake and SSL will then refuse further communication
with that server.

Certificate Verification with NSS
---------------------------------

If libcurl was built with NSS support, then depending on the OS distribution,
it is probably required to take some additional steps to use the system-wide
CA cert db. RedHat ships with an additional module, libnsspem.so, which
enables NSS to read the OpenSSL PEM CA bundle. On openSUSE you can install
p11-kit-nss-trust which makes NSS use the system wide CA certificate store. NSS
also has a new [database format](https://wiki.mozilla.org/NSS_Shared_DB).

Starting with version 7.19.7, libcurl automatically adds the 'sql:' prefix to
the certdb directory (either the hardcoded default /etc/pki/nssdb or the
directory configured with SSL_DIR environment variable). To check which certdb
format your distribution provides, examine the default certdb location:
/etc/pki/nssdb; the new certdb format can be identified by the filenames
cert9.db, key4.db, pkcs11.txt; filenames of older versions are cert8.db,
key3.db, secmod.db.

Certificate Verification with Schannel and Secure Transport
-----------------------------------------------------------

If libcurl was built with Schannel (Microsoft's native TLS engine) or Secure
Transport (Apple's native TLS engine) support, then libcurl will still perform
peer certificate verification, but instead of using a CA cert bundle, it will
use the certificates that are built into the OS. These are the same
certificates that appear in the Internet Options control panel (under Windows)
or Keychain Access application (under OS X). Any custom security rules for
certificates will be honored.

Schannel will run CRL checks on certificates unless peer verification is
disabled. Secure Transport on iOS will run OCSP checks on certificates unless
peer verification is disabled. Secure Transport on OS X will run either OCSP
or CRL checks on certificates if those features are enabled, and this behavior
can be adjusted in the preferences of Keychain Access.

HTTPS proxy
-----------

Since version 7.52.0, curl can do HTTPS to the proxy separately from the
connection to the server. This TLS connection is handled separately from the
server connection so instead of `--insecure` and `--cacert` to control the
certificate verification, you use `--proxy-insecure` and `--proxy-cacert`.
With these options, you make sure that the TLS connection and the trust of the
proxy can be kept totally separate from the TLS connection to the server.
# Hyper

Hyper is a separate HTTP library written in Rust. curl can be told to use this
library as a backend to deal with HTTP.

## Experimental!

Hyper support in curl is considered **EXPERIMENTAL** until further notice. It
needs to be explicitly enabled at build-time.

Further development and tweaking of the Hyper backend support in curl will
happen in in the master branch using pull-requests, just like ordinary
changes.

## Hyper version

The C API for Hyper is brand new and is still under development.

## build curl with hyper

Build hyper and enable the C API:

     % git clone https://github.com/hyperium/hyper
     % cd hyper
     % RUSTFLAGS="--cfg hyper_unstable_ffi" cargo build --features client,http1,http2,ffi

Build curl to use hyper's C API:

     % git clone https://github.com/curl/curl
     % cd curl
     % ./buildconf
     % ./configure --with-hyper=<hyper dir>
     % make

# using Hyper internally

Hyper is a low level HTTP transport library. curl itself provides all HTTP
headers and Hyper provides all received headers back to curl.

Therefore, most of the "header logic" in curl as in responding to and acting
on specific input and output headers are done the same way in curl code.

The API in Hyper delivers received HTTP headers as (cleaned up) name=value
pairs, making it impossible for curl to know the exact byte representation
over the wire with Hyper.

## Limitations

The hyper backend does not support

- `CURLOPT_IGNORE_CONTENT_LENGTH`
- `--raw` and disabling `CURLOPT_HTTP_TRANSFER_DECODING`
- RTSP
- hyper is much stricter about what HTTP header contents it allow in requests
- HTTP/0.9

## Remaining issues

This backend is still not feature complete with the native backend. Areas that
still need attention and verification include:

- multiplexed HTTP/2
- h2 Upgrade:
- pausing transfers
- receiving HTTP/1 trailers
- sending HTTP/1 trailers

# curl the next few years - perhaps

Roadmap of things Daniel Stenberg wants to work on next. It is intended to
serve as a guideline for others for information, feedback and possible
participation.

## "Complete" the HTTP/3 support

curl has experimental support for HTTP/3 since a good while back. There are
some functionality missing and once the final specs are published we want to
eventually remove the "experimental" label from this functionality.

## HTTPS DNS records

As a DNS version of alt-svc and also a pre-requisite for ECH (see below).

See: https://datatracker.ietf.org/doc/html/draft-ietf-dnsop-svcb-https-02

## ECH (Encrypted Client Hello - formerly known as ESNI)

 See Daniel's post on [Support of Encrypted
 SNI](https://curl.se/mail/lib-2019-03/0000.html) on the mailing list.

 Initial work exists in https://github.com/curl/curl/pull/4011
# How to get started helping out in the curl project

We are always in need of more help. If you are new to the project and are
looking for ways to contribute and help out, this document aims to give a few
good starting points.

A good idea is to start by subscribing to the [curl-library mailing
list](https://lists.haxx.se/listinfo/curl-library) to keep track of the
current discussion topics.

## Scratch your own itch

One of the best ways is to start working on any problems or issues you have
found yourself or perhaps got annoyed at in the past. It can be a spelling
error in an error text or a weirdly phrased section in a man page. Hunt it
down and report the bug. Or make your first pull request with a fix for that.

## Smaller tasks

Some projects mark small issues as "beginner friendly", "bite-sized" or
similar. We do not do that in curl since such issues never linger around long
enough. Simple issues get handled fast.

If you are looking for a smaller or simpler task in the project to help out
with as an entry-point into the project, perhaps because you are a newcomer or
even maybe not a terribly experienced developer, here's our advice:

 - Read through this document to get a grasp on a general approach to use
 - Consider adding a test case for something not currently tested (correctly)
 - Consider updating or adding documentation
 - One way to get started gently in the project, is to participate in an
   existing issue/PR and help out by reproducing the issue, review the code in
   the PR etc.

## Help wanted

In the issue tracker we occasionally mark bugs with [help
wanted](https://github.com/curl/curl/labels/help%20wanted), as a sign that the
bug is acknowledged to exist and that there's nobody known to work on this
issue for the moment. Those are bugs that are fine to "grab" and provide a
pull request for. The complexity level of these will of course vary, so pick
one that piques your interest.

## Work on known bugs

Some bugs are known and have not yet received attention and work enough to get
fixed. We collect such known existing flaws in the
[KNOWN_BUGS](https://curl.se/docs/knownbugs.html) page. Many of them link
to the original bug report with some additional details, but some may also
have aged a bit and may require some verification that the bug still exists in
the same way and that what was said about it in the past is still valid.

## Fix autobuild problems

On the [autobuilds page](https://curl.se/dev/builds.html) we show a
collection of test results from the automatic curl build and tests that are
performed by volunteers. Fixing compiler warnings and errors shown there is
something we value greatly. Also, if you own or run systems or architectures
that are not already tested in the autobuilds, we also appreciate more
volunteers running builds automatically to help us keep curl portable.

## TODO items

Ideas for features and functions that we have considered worthwhile to
implement and provide are kept in the
[TODO](https://curl.se/docs/todo.html) file. Some of the ideas are
rough. Some are well thought out. Some probably are not really suitable
anymore.

Before you invest a lot of time on a TODO item, do bring it up for discussion
on the mailing list. For discussion on applicability but also for ideas and
brainstorming on specific ways to do the implementation etc.

## You decide

You can also come up with a completely new thing you think we should do. Or
not do. Or fix. Or add to the project. You then either bring it to the mailing
list first to see if people will shoot down the idea at once, or you bring a
first draft of the idea as a pull request and take the discussion there around
the specific implementation. Either way is fine.

## CONTRIBUTE

We offer [guidelines](https://curl.se/dev/contribute.html) that are
suitable to be familiar with before you decide to contribute to curl. If
you are used to open source development, you will probably not find many
surprises in there.
Version Numbers and Releases
============================

 Curl is not only curl. Curl is also libcurl. they are actually individually
 versioned, but they usually follow each other closely.

 The version numbering is always built up using the same system:

        X.Y.Z

  - X is main version number
  - Y is release number
  - Z is patch number

## Bumping numbers

 One of these numbers will get bumped in each new release. The numbers to the
 right of a bumped number will be reset to zero.

 The main version number will get bumped when *really* big, world colliding
 changes are made. The release number is bumped when changes are performed or
 things/features are added. The patch number is bumped when the changes are
 mere bugfixes.

 It means that after release 1.2.3, we can release 2.0.0 if something really
 big has been made, 1.3.0 if not that big changes were made or 1.2.4 if only
 bugs were fixed.

 Bumping, as in increasing the number with 1, is unconditionally only
 affecting one of the numbers (except the ones to the right of it, that may be
 set to zero). 1 becomes 2, 3 becomes 4, 9 becomes 10, 88 becomes 89 and 99
 becomes 100. So, after 1.2.9 comes 1.2.10. After 3.99.3, 3.100.0 might come.

 All original curl source release archives are named according to the libcurl
 version (not according to the curl client version that, as said before, might
 differ).

 As a service to any application that might want to support new libcurl
 features while still being able to build with older versions, all releases
 have the libcurl version stored in the curl/curlver.h file using a static
 numbering scheme that can be used for comparison. The version number is
 defined as:

```c
#define LIBCURL_VERSION_NUM 0xXXYYZZ
```

 Where XX, YY and ZZ are the main version, release and patch numbers in
 hexadecimal. All three number fields are always represented using two digits
 (eight bits each). 1.2 would appear as "0x010200" while version 9.11.7
 appears as "0x090b07".

 This 6-digit hexadecimal number is always a greater number in a more recent
 release. It makes comparisons with greater than and less than work.

 This number is also available as three separate defines:
 `LIBCURL_VERSION_MAJOR`, `LIBCURL_VERSION_MINOR` and `LIBCURL_VERSION_PATCH`.
# dynbuf

This is the internal module for creating and handling "dynamic buffers". This
means buffers that can be appended to, dynamically and grow to adapt.

There will always be a terminating zero put at the end of the dynamic buffer.

The `struct dynbuf` is used to hold data for each instance of a dynamic
buffer. The members of that struct **MUST NOT** be accessed or modified
without using the dedicated dynbuf API.

## init

```c
void Curl_dyn_init(struct dynbuf *s, size_t toobig);
```

This inits a struct to use for dynbuf and it cannot fail. The `toobig` value
**must** be set to the maximum size we allow this buffer instance to grow to.
The functions below will return `CURLE_OUT_OF_MEMORY` when hitting this limit.

## free

```c
void Curl_dyn_free(struct dynbuf *s);
```

Free the associated memory and clean up. After a free, the `dynbuf` struct can
be re-used to start appending new data to.

## addn

```c
CURLcode Curl_dyn_addn(struct dynbuf *s, const void *mem, size_t len);
```

Append arbitrary data of a given length to the end of the buffer.

## add

```c
CURLcode Curl_dyn_add(struct dynbuf *s, const char *str);
```

Append a C string to the end of the buffer.

## addf

```c
CURLcode Curl_dyn_addf(struct dynbuf *s, const char *fmt, ...);
```

Append a `printf()`-style string to the end of the buffer.

## vaddf

```c
CURLcode Curl_dyn_vaddf(struct dynbuf *s, const char *fmt, va_list ap);
```

Append a `vprintf()`-style string to the end of the buffer.

## reset

```c
void Curl_dyn_reset(struct dynbuf *s);
```

Reset the buffer length, but leave the allocation.

## tail

```c
CURLcode Curl_dyn_tail(struct dynbuf *s, size_t length);
```

Keep `length` bytes of the buffer tail (the last `length` bytes of the
buffer). The rest of the buffer is dropped. The specified `length` must not be
larger than the buffer length.

## ptr

```c
char *Curl_dyn_ptr(const struct dynbuf *s);
```

Returns a `char *` to the buffer if it has a length, otherwise a NULL. Since
the buffer may be reallocated, this pointer should not be trusted or used
anymore after the next buffer manipulation call.

## uptr

```c
unsigned char *Curl_dyn_uptr(const struct dynbuf *s);
```

Returns an `unsigned char *` to the buffer if it has a length, otherwise a
NULL. Since the buffer may be reallocated, this pointer should not be trusted
or used anymore after the next buffer manipulation call.

## len

```c
size_t Curl_dyn_len(const struct dynbuf *s);
```

Returns the length of the buffer in bytes. Does not include the terminating
zero byte.
# Items to be removed from future curl releases

If any of these deprecated features is a cause for concern for you, please
email the
[curl-library mailing list](https://lists.haxx.se/listinfo/curl-library)
as soon as possible and explain to us why this is a problem for you and
how your use case cannot be satisfied properly using a workaround.

## Past removals

 - Pipelining
 - axTLS
 - PolarSSL
# How to do code reviews for curl

Anyone and everyone is encouraged and welcome to review code submissions in
curl. This is a guide on what to check for and how to perform a successful
code review.

## All submissions should get reviewed

All pull requests and patches submitted to the project should be reviewed by
at least one experienced curl maintainer before that code is accepted and
merged.

## Let the tools and tests take the first rounds

On initial pull requests, let the tools and tests do their job first and then
start out by helping the submitter understand the test failures and tool
alerts.

## How to provide feedback to author

Be nice. Ask questions. Provide examples or suggestions of improvements.
Assume the best intentions. Remember language barriers.

All first-time contributors can become regulars. Let's help them go there.

## Is this a change we want?

If this is not a change that seems to be aligned with the project's path
forward and as such cannot be accepted, inform the author about this sooner
rather than later. Do it gently and explain why and possibly what could be
done to make it more acceptable.

## API/ABI stability or changed behavior

Changing the API and the ABI may be fine in a change but it needs to be done
deliberately and carefully. If not, a reviewer must help the author to realize
the mistake.

curl and libcurl are similarly strict on not modifying existing behavior. API
and ABI stability is not enough, the behavior should also remain intact as far
as possible.

## Code style

Most code style nits are detected by checksrc but not all. Only leave remarks
on style deviation once checksrc does not find anymore.

Minor nits from fresh submitters can also be handled by the maintainer when
merging, in case it seems like the submitter is not clear on what to do. We
want to make the process fun and exciting for new contributors.

## Encourage consistency

Make sure new code is written in a similar style as existing code. Naming,
logic, conditions, etc.

## Are pointers always non-NULL?

If a function or code rely on pointers being non-NULL, take an extra look if
that seems to be a fair assessment.

## Asserts

Conditions that should never be false can be verified with `DEBUGASSERT()`
calls to get caught in tests and debugging easier, while not having an impact
on final or release builds.

## Memory allocation

Can the mallocs be avoided? Do not introduce mallocs in any hot paths. If
there are (new) mallocs, can they be combined into fewer calls?

Are all allocations handled in errorpaths to avoid leaks and crashes?

## Thread-safety

We do not like static variables as they break thread-safety and prevent
functions from being reentrant.

## Should features be `#ifdef`ed?

Features and functionality may not be present everywhere and should therefore
be `#ifdef`ed. Additionally, some features should be possible to switch on/off
in the build.

Write `#ifdef`s to be as little of a "maze" as possible.

## Does it look portable enough?

curl runs "everywhere". Does the code take a reasonable stance and enough
precautions to be possible to build and run on most platforms?

Remember that we live by C89 restrictions.

## Tests and testability

New features should be added in conjunction with one or more test cases.
Ideally, functions should also be written so that unit tests can be done to
test individual functions.

## Documentation

New features or changes to existing functionality **must** be accompanied by
updated documentation. Submitting that in a separate follow-up pull request is
not OK. A code review must also verify that the submitted documentation update
matches the code submission.

English is not everyone's first language, be mindful of this and help the
submitter improve the text if it needs a rewrite to read better.

## Code should not be hard to understand

Source code should be written to maximize readability and be easy to
understand.

## Functions should not be large

A single function should never be large as that makes it hard to follow and
understand all the exit points and state changes. Some existing functions in
curl certainly violate this ground rule but when reviewing new code we should
propose splitting into smaller functions.

## Duplication is evil

Anything that looks like duplicated code is a red flag. Anything that seems to
introduce code that we *should* already have or provide needs a closer check.

## Sensitive data

When credentials are involved, take an extra look at what happens with this
data. Where it comes from and where it goes.

## Variable types differ

`size_t` is not a fixed size. `time_t` can be signed or unsigned and have
different sizes. Relying on variable sizes is a red flag.

Also remember that endianness and >= 32 bit accesses to unaligned addresses
are problematic areas.

## Integer overflows

Be careful about integer overflows. Some variable types can be either 32 bit
or 64 bit. Integer overflows must be detected and acted on *before* they
happen.

## Dangerous use of functions

Maybe use of `realloc()` should rather use the dynbuf functions?

Do not allow new code that grows buffers without using dynbuf.

Use of C functions that rely on a terminating zero must only be used on data
that really do have a zero terminating zero.

## Dangerous "data styles"

Make extra precautions and verify that memory buffers that need a terminating
zero always have exactly that. Buffers *without* a zero terminator must not be
used as input to string functions.

# Commit messages

Tightly coupled with a code review is making sure that the commit message is
good. It is the responsibility of the person who merges the code to make sure
that the commit message follows our standard (detailed in the
[CONTRIBUTE.md](CONTRIBUTE.md) document). This includes making sure the PR
identifies related issues and giving credit to reporters and helpers.
# The Art Of Scripting HTTP Requests Using Curl

## Background

 This document assumes that you are familiar with HTML and general networking.

 The increasing amount of applications moving to the web has made "HTTP
 Scripting" more frequently requested and wanted. To be able to automatically
 extract information from the web, to fake users, to post or upload data to
 web servers are all important tasks today.

 Curl is a command line tool for doing all sorts of URL manipulations and
 transfers, but this particular document will focus on how to use it when
 doing HTTP requests for fun and profit. I will assume that you know how to
 invoke `curl --help` or `curl --manual` to get basic information about it.

 Curl is not written to do everything for you. It makes the requests, it gets
 the data, it sends data and it retrieves the information. You probably need
 to glue everything together using some kind of script language or repeated
 manual invokes.

## The HTTP Protocol

 HTTP is the protocol used to fetch data from web servers. It is a simple
 protocol that is built upon TCP/IP. The protocol also allows information to
 get sent to the server from the client using a few different methods, as will
 be shown here.

 HTTP is plain ASCII text lines being sent by the client to a server to
 request a particular action, and then the server replies a few text lines
 before the actual requested content is sent to the client.

 The client, curl, sends a HTTP request. The request contains a method (like
 GET, POST, HEAD etc), a number of request headers and sometimes a request
 body. The HTTP server responds with a status line (indicating if things went
 well), response headers and most often also a response body. The "body" part
 is the plain data you requested, like the actual HTML or the image etc.

## See the Protocol

 Using curl's option [`--verbose`](https://curl.se/docs/manpage.html#-v)
 (`-v` as a short option) will display what kind of commands curl sends to the
 server, as well as a few other informational texts.

 `--verbose` is the single most useful option when it comes to debug or even
 understand the curl<->server interaction.

 Sometimes even `--verbose` is not enough. Then
 [`--trace`](https://curl.se/docs/manpage.html#-trace) and
 [`--trace-ascii`](https://curl.se/docs/manpage.html#--trace-ascii)
 offer even more details as they show **everything** curl sends and
 receives. Use it like this:

    curl --trace-ascii debugdump.txt http://www.example.com/

## See the Timing

 Many times you may wonder what exactly is taking all the time, or you just
 want to know the amount of milliseconds between two points in a transfer. For
 those, and other similar situations, the
 [`--trace-time`](https://curl.se/docs/manpage.html#--trace-time) option
 is what you need. It will prepend the time to each trace output line:

    curl --trace-ascii d.txt --trace-time http://example.com/

## See the Response

 By default curl sends the response to stdout. You need to redirect it
 somewhere to avoid that, most often that is done with ` -o` or `-O`.

# URL

## Spec

 The Uniform Resource Locator format is how you specify the address of a
 particular resource on the Internet. You know these, you have seen URLs like
 https://curl.se or https://yourbank.com a million times. RFC 3986 is the
 canonical spec. And yeah, the formal name is not URL, it is URI.

## Host

 The host name is usually resolved using DNS or your /etc/hosts file to an IP
 address and that is what curl will communicate with. Alternatively you specify
 the IP address directly in the URL instead of a name.

 For development and other trying out situations, you can point to a different
 IP address for a host name than what would otherwise be used, by using curl's
 [`--resolve`](https://curl.se/docs/manpage.html#--resolve) option:

    curl --resolve www.example.org:80:127.0.0.1 http://www.example.org/

## Port number

 Each protocol curl supports operates on a default port number, be it over TCP
 or in some cases UDP. Normally you do not have to take that into
 consideration, but at times you run test servers on other ports or
 similar. Then you can specify the port number in the URL with a colon and a
 number immediately following the host name. Like when doing HTTP to port
 1234:

    curl http://www.example.org:1234/

 The port number you specify in the URL is the number that the server uses to
 offer its services. Sometimes you may use a proxy, and then you may
 need to specify that proxy's port number separately from what curl needs to
 connect to the server. Like when using a HTTP proxy on port 4321:

    curl --proxy http://proxy.example.org:4321 http://remote.example.org/

## User name and password

 Some services are setup to require HTTP authentication and then you need to
 provide name and password which is then transferred to the remote site in
 various ways depending on the exact authentication protocol used.

 You can opt to either insert the user and password in the URL or you can
 provide them separately:

    curl http://user:password@example.org/

 or

    curl -u user:password http://example.org/

 You need to pay attention that this kind of HTTP authentication is not what
 is usually done and requested by user-oriented websites these days. They tend
 to use forms and cookies instead.

## Path part

 The path part is just sent off to the server to request that it sends back
 the associated response. The path is what is to the right side of the slash
 that follows the host name and possibly port number.

# Fetch a page

## GET

 The simplest and most common request/operation made using HTTP is to GET a
 URL. The URL could itself refer to a web page, an image or a file. The client
 issues a GET request to the server and receives the document it asked for.
 If you issue the command line

    curl https://curl.se

 you get a web page returned in your terminal window. The entire HTML document
 that that URL holds.

 All HTTP replies contain a set of response headers that are normally hidden,
 use curl's [`--include`](https://curl.se/docs/manpage.html#-i) (`-i`)
 option to display them as well as the rest of the document.

## HEAD

 You can ask the remote server for ONLY the headers by using the
 [`--head`](https://curl.se/docs/manpage.html#-I) (`-I`) option which
 will make curl issue a HEAD request. In some special cases servers deny the
 HEAD method while others still work, which is a particular kind of annoyance.

 The HEAD method is defined and made so that the server returns the headers
 exactly the way it would do for a GET, but without a body. It means that you
 may see a `Content-Length:` in the response headers, but there must not be an
 actual body in the HEAD response.

## Multiple URLs in a single command line

 A single curl command line may involve one or many URLs. The most common case
 is probably to just use one, but you can specify any amount of URLs. Yes
 any. No limits. You will then get requests repeated over and over for all the
 given URLs.

 Example, send two GETs:

    curl http://url1.example.com http://url2.example.com

 If you use [`--data`](https://curl.se/docs/manpage.html#-d) to POST to
 the URL, using multiple URLs means that you send that same POST to all the
 given URLs.

 Example, send two POSTs:

    curl --data name=curl http://url1.example.com http://url2.example.com


## Multiple HTTP methods in a single command line

 Sometimes you need to operate on several URLs in a single command line and do
 different HTTP methods on each. For this, you will enjoy the
 [`--next`](https://curl.se/docs/manpage.html#-:) option. It is basically
 a separator that separates a bunch of options from the next. All the URLs
 before `--next` will get the same method and will get all the POST data
 merged into one.

 When curl reaches the `--next` on the command line, it will sort of reset the
 method and the POST data and allow a new set.

 Perhaps this is best shown with a few examples. To send first a HEAD and then
 a GET:

    curl -I http://example.com --next http://example.com

 To first send a POST and then a GET:

    curl -d score=10 http://example.com/post.cgi --next http://example.com/results.html

# HTML forms

## Forms explained

 Forms are the general way a website can present a HTML page with fields for
 the user to enter data in, and then press some kind of 'OK' or 'Submit'
 button to get that data sent to the server. The server then typically uses
 the posted data to decide how to act. Like using the entered words to search
 in a database, or to add the info in a bug tracking system, display the
 entered address on a map or using the info as a login-prompt verifying that
 the user is allowed to see what it is about to see.

 Of course there has to be some kind of program on the server end to receive
 the data you send. You cannot just invent something out of the air.

## GET

 A GET-form uses the method GET, as specified in HTML like:

```html
<form method="GET" action="junk.cgi">
  <input type=text name="birthyear">
  <input type=submit name=press value="OK">
</form>
```

 In your favorite browser, this form will appear with a text box to fill in
 and a press-button labeled "OK". If you fill in '1905' and press the OK
 button, your browser will then create a new URL to get for you. The URL will
 get `junk.cgi?birthyear=1905&press=OK` appended to the path part of the
 previous URL.

 If the original form was seen on the page `www.example.com/when/birth.html`,
 the second page you will get will become
 `www.example.com/when/junk.cgi?birthyear=1905&press=OK`.

 Most search engines work this way.

 To make curl do the GET form post for you, just enter the expected created
 URL:

    curl "http://www.example.com/when/junk.cgi?birthyear=1905&press=OK"

## POST

 The GET method makes all input field names get displayed in the URL field of
 your browser. That is generally a good thing when you want to be able to
 bookmark that page with your given data, but it is an obvious disadvantage if
 you entered secret information in one of the fields or if there are a large
 amount of fields creating a long and unreadable URL.

 The HTTP protocol then offers the POST method. This way the client sends the
 data separated from the URL and thus you will not see any of it in the URL
 address field.

 The form would look similar to the previous one:

```html
<form method="POST" action="junk.cgi">
  <input type=text name="birthyear">
  <input type=submit name=press value=" OK ">
</form>
```

 And to use curl to post this form with the same data filled in as before, we
 could do it like:

    curl --data "birthyear=1905&press=%20OK%20" http://www.example.com/when.cgi

 This kind of POST will use the Content-Type
 `application/x-www-form-urlencoded` and is the most widely used POST kind.

 The data you send to the server MUST already be properly encoded, curl will
 not do that for you. For example, if you want the data to contain a space,
 you need to replace that space with `%20`, etc. Failing to comply with this will
 most likely cause your data to be received wrongly and messed up.

 Recent curl versions can in fact url-encode POST data for you, like this:

    curl --data-urlencode "name=I am Daniel" http://www.example.com

 If you repeat `--data` several times on the command line, curl will
 concatenate all the given data pieces - and put a `&` symbol between each
 data segment.

## File Upload POST

 Back in late 1995 they defined an additional way to post data over HTTP. It
 is documented in the RFC 1867, why this method sometimes is referred to as
 RFC1867-posting.

 This method is mainly designed to better support file uploads. A form that
 allows a user to upload a file could be written like this in HTML:

```html
<form method="POST" enctype='multipart/form-data' action="upload.cgi">
  <input type=file name=upload>
  <input type=submit name=press value="OK">
</form>
```

 This clearly shows that the Content-Type about to be sent is
 `multipart/form-data`.

 To post to a form like this with curl, you enter a command line like:

    curl --form upload=@localfilename --form press=OK [URL]

## Hidden Fields

 A common way for HTML based applications to pass state information between
 pages is to add hidden fields to the forms. Hidden fields are already filled
 in, they are not displayed to the user and they get passed along just as all
 the other fields.

 A similar example form with one visible field, one hidden field and one
 submit button could look like:

```html
<form method="POST" action="foobar.cgi">
  <input type=text name="birthyear">
  <input type=hidden name="person" value="daniel">
  <input type=submit name="press" value="OK">
</form>
```

 To POST this with curl, you will not have to think about if the fields are
 hidden or not. To curl they are all the same:

    curl --data "birthyear=1905&press=OK&person=daniel" [URL]

## Figure Out What A POST Looks Like

 When you are about fill in a form and send to a server by using curl instead
 of a browser, you are of course interested in sending a POST exactly the way
 your browser does.

 An easy way to get to see this, is to save the HTML page with the form on
 your local disk, modify the 'method' to a GET, and press the submit button
 (you could also change the action URL if you want to).

 You will then clearly see the data get appended to the URL, separated with a
 `?`-letter as GET forms are supposed to.

# HTTP upload

## PUT

 Perhaps the best way to upload data to a HTTP server is to use PUT. Then
 again, this of course requires that someone put a program or script on the
 server end that knows how to receive a HTTP PUT stream.

 Put a file to a HTTP server with curl:

    curl --upload-file uploadfile http://www.example.com/receive.cgi

# HTTP Authentication

## Basic Authentication

 HTTP Authentication is the ability to tell the server your username and
 password so that it can verify that you are allowed to do the request you are
 doing. The Basic authentication used in HTTP (which is the type curl uses by
 default) is **plain text** based, which means it sends username and password
 only slightly obfuscated, but still fully readable by anyone that sniffs on
 the network between you and the remote server.

 To tell curl to use a user and password for authentication:

    curl --user name:password http://www.example.com

## Other Authentication

 The site might require a different authentication method (check the headers
 returned by the server), and then
 [`--ntlm`](https://curl.se/docs/manpage.html#--ntlm),
 [`--digest`](https://curl.se/docs/manpage.html#--digest),
 [`--negotiate`](https://curl.se/docs/manpage.html#--negotiate) or even
 [`--anyauth`](https://curl.se/docs/manpage.html#--anyauth) might be
 options that suit you.

## Proxy Authentication

 Sometimes your HTTP access is only available through the use of a HTTP
 proxy. This seems to be especially common at various companies. A HTTP proxy
 may require its own user and password to allow the client to get through to
 the Internet. To specify those with curl, run something like:

    curl --proxy-user proxyuser:proxypassword curl.se

 If your proxy requires the authentication to be done using the NTLM method,
 use [`--proxy-ntlm`](https://curl.se/docs/manpage.html#--proxy-ntlm), if
 it requires Digest use
 [`--proxy-digest`](https://curl.se/docs/manpage.html#--proxy-digest).

 If you use any one of these user+password options but leave out the password
 part, curl will prompt for the password interactively.

## Hiding credentials

 Do note that when a program is run, its parameters might be possible to see
 when listing the running processes of the system. Thus, other users may be
 able to watch your passwords if you pass them as plain command line
 options. There are ways to circumvent this.

 It is worth noting that while this is how HTTP Authentication works, many
 websites will not use this concept when they provide logins etc. See the Web
 Login chapter further below for more details on that.

# More HTTP Headers

## Referer

 A HTTP request may include a 'referer' field (yes it is misspelled), which
 can be used to tell from which URL the client got to this particular
 resource. Some programs/scripts check the referer field of requests to verify
 that this was not arriving from an external site or an unknown page. While
 this is a stupid way to check something so easily forged, many scripts still
 do it. Using curl, you can put anything you want in the referer-field and
 thus more easily be able to fool the server into serving your request.

 Use curl to set the referer field with:

    curl --referer http://www.example.come http://www.example.com

## User Agent

 Similar to the referer field, all HTTP requests may set the User-Agent
 field. It names what user agent (client) that is being used. Many
 applications use this information to decide how to display pages. Silly web
 programmers try to make different pages for users of different browsers to
 make them look the best possible for their particular browsers. They usually
 also do different kinds of javascript, vbscript etc.

 At times, you will see that getting a page with curl will not return the same
 page that you see when getting the page with your browser. Then you know it
 is time to set the User Agent field to fool the server into thinking you are
 one of those browsers.

 To make curl look like Internet Explorer 5 on a Windows 2000 box:

    curl --user-agent "Mozilla/4.0 (compatible; MSIE 5.01; Windows NT 5.0)" [URL]

 Or why not look like you are using Netscape 4.73 on an old Linux box:

    curl --user-agent "Mozilla/4.73 [en] (X11; U; Linux 2.2.15 i686)" [URL]

## Redirects

## Location header

 When a resource is requested from a server, the reply from the server may
 include a hint about where the browser should go next to find this page, or a
 new page keeping newly generated output. The header that tells the browser to
 redirect is `Location:`.

 Curl does not follow `Location:` headers by default, but will simply display
 such pages in the same manner it displays all HTTP replies. It does however
 feature an option that will make it attempt to follow the `Location:`
 pointers.

 To tell curl to follow a Location:

    curl --location http://www.example.com

 If you use curl to POST to a site that immediately redirects you to another
 page, you can safely use
 [`--location`](https://curl.se/docs/manpage.html#-L) (`-L`) and
 `--data`/`--form` together. Curl will only use POST in the first request, and
 then revert to GET in the following operations.

## Other redirects

 Browser typically support at least two other ways of redirects that curl
 does not: first the html may contain a meta refresh tag that asks the browser
 to load a specific URL after a set number of seconds, or it may use
 javascript to do it.

# Cookies

## Cookie Basics

 The way the web browsers do "client side state control" is by using
 cookies. Cookies are just names with associated contents. The cookies are
 sent to the client by the server. The server tells the client for what path
 and host name it wants the cookie sent back, and it also sends an expiration
 date and a few more properties.

 When a client communicates with a server with a name and path as previously
 specified in a received cookie, the client sends back the cookies and their
 contents to the server, unless of course they are expired.

 Many applications and servers use this method to connect a series of requests
 into a single logical session. To be able to use curl in such occasions, we
 must be able to record and send back cookies the way the web application
 expects them. The same way browsers deal with them.

## Cookie options

 The simplest way to send a few cookies to the server when getting a page with
 curl is to add them on the command line like:

    curl --cookie "name=Daniel" http://www.example.com

 Cookies are sent as common HTTP headers. This is practical as it allows curl
 to record cookies simply by recording headers. Record cookies with curl by
 using the [`--dump-header`](https://curl.se/docs/manpage.html#-D) (`-D`)
 option like:

    curl --dump-header headers_and_cookies http://www.example.com

 (Take note that the
 [`--cookie-jar`](https://curl.se/docs/manpage.html#-c) option described
 below is a better way to store cookies.)

 Curl has a full blown cookie parsing engine built-in that comes in use if you
 want to reconnect to a server and use cookies that were stored from a
 previous connection (or hand-crafted manually to fool the server into
 believing you had a previous connection). To use previously stored cookies,
 you run curl like:

    curl --cookie stored_cookies_in_file http://www.example.com

 Curl's "cookie engine" gets enabled when you use the
 [`--cookie`](https://curl.se/docs/manpage.html#-b) option. If you only
 want curl to understand received cookies, use `--cookie` with a file that
 does not exist. Example, if you want to let curl understand cookies from a
 page and follow a location (and thus possibly send back cookies it received),
 you can invoke it like:

    curl --cookie nada --location http://www.example.com

 Curl has the ability to read and write cookie files that use the same file
 format that Netscape and Mozilla once used. It is a convenient way to share
 cookies between scripts or invokes. The `--cookie` (`-b`) switch
 automatically detects if a given file is such a cookie file and parses it,
 and by using the `--cookie-jar` (`-c`) option you will make curl write a new
 cookie file at the end of an operation:

    curl --cookie cookies.txt --cookie-jar newcookies.txt \
    http://www.example.com

# HTTPS

## HTTPS is HTTP secure

 There are a few ways to do secure HTTP transfers. By far the most common
 protocol for doing this is what is generally known as HTTPS, HTTP over
 SSL. SSL encrypts all the data that is sent and received over the network and
 thus makes it harder for attackers to spy on sensitive information.

 SSL (or TLS as the latest version of the standard is called) offers a
 truckload of advanced features to allow all those encryptions and key
 infrastructure mechanisms encrypted HTTP requires.

 Curl supports encrypted fetches when built to use a TLS library and it can be
 built to use one out of a fairly large set of libraries - `curl -V` will show
 which one your curl was built to use (if any!). To get a page from a HTTPS
 server, simply run curl like:

    curl https://secure.example.com

## Certificates

 In the HTTPS world, you use certificates to validate that you are the one
 you claim to be, as an addition to normal passwords. Curl supports client-
 side certificates. All certificates are locked with a pass phrase, which you
 need to enter before the certificate can be used by curl. The pass phrase
 can be specified on the command line or if not, entered interactively when
 curl queries for it. Use a certificate with curl on a HTTPS server like:

    curl --cert mycert.pem https://secure.example.com

 curl also tries to verify that the server is who it claims to be, by
 verifying the server's certificate against a locally stored CA cert
 bundle. Failing the verification will cause curl to deny the connection. You
 must then use [`--insecure`](https://curl.se/docs/manpage.html#-k)
 (`-k`) in case you want to tell curl to ignore that the server cannot be
 verified.

 More about server certificate verification and ca cert bundles can be read in
 the [SSLCERTS document](https://curl.se/docs/sslcerts.html).

 At times you may end up with your own CA cert store and then you can tell
 curl to use that to verify the server's certificate:

    curl --cacert ca-bundle.pem https://example.com/

# Custom Request Elements

## Modify method and headers

 Doing fancy stuff, you may need to add or change elements of a single curl
 request.

 For example, you can change the POST request to a PROPFIND and send the data
 as `Content-Type: text/xml` (instead of the default Content-Type) like this:

    curl --data "<xml>" --header "Content-Type: text/xml" \
      --request PROPFIND example.com

 You can delete a default header by providing one without content. Like you
 can ruin the request by chopping off the Host: header:

    curl --header "Host:" http://www.example.com

 You can add headers the same way. Your server may want a `Destination:`
 header, and you can add it:

    curl --header "Destination: http://nowhere" http://example.com

## More on changed methods

 It should be noted that curl selects which methods to use on its own
 depending on what action to ask for. `-d` will do POST, `-I` will do HEAD and
 so on. If you use the
 [`--request`](https://curl.se/docs/manpage.html#-X) / `-X` option you
 can change the method keyword curl selects, but you will not modify curl's
 behavior. This means that if you for example use -d "data" to do a POST, you
 can modify the method to a `PROPFIND` with `-X` and curl will still think it
 sends a POST . You can change the normal GET to a POST method by simply
 adding `-X POST` in a command line like:

    curl -X POST http://example.org/

 ... but curl will still think and act as if it sent a GET so it will not send
 any request body etc.

# Web Login

## Some login tricks

 While not strictly just HTTP related, it still causes a lot of people
 problems so here's the executive run-down of how the vast majority of all
 login forms work and how to login to them using curl.

 It can also be noted that to do this properly in an automated fashion, you
 will most certainly need to script things and do multiple curl invokes etc.

 First, servers mostly use cookies to track the logged-in status of the
 client, so you will need to capture the cookies you receive in the
 responses. Then, many sites also set a special cookie on the login page (to
 make sure you got there through their login page) so you should make a habit
 of first getting the login-form page to capture the cookies set there.

 Some web-based login systems feature various amounts of javascript, and
 sometimes they use such code to set or modify cookie contents. Possibly they
 do that to prevent programmed logins, like this manual describes how to...
 Anyway, if reading the code is not enough to let you repeat the behavior
 manually, capturing the HTTP requests done by your browsers and analyzing the
 sent cookies is usually a working method to work out how to shortcut the
 javascript need.

 In the actual `<form>` tag for the login, lots of sites fill-in
 random/session or otherwise secretly generated hidden tags and you may need
 to first capture the HTML code for the login form and extract all the hidden
 fields to be able to do a proper login POST. Remember that the contents need
 to be URL encoded when sent in a normal POST.

# Debug

## Some debug tricks

 Many times when you run curl on a site, you will notice that the site does not
 seem to respond the same way to your curl requests as it does to your
 browser's.

 Then you need to start making your curl requests more similar to your
 browser's requests:

 - Use the `--trace-ascii` option to store fully detailed logs of the requests
   for easier analyzing and better understanding

 - Make sure you check for and use cookies when needed (both reading with
   `--cookie` and writing with `--cookie-jar`)

 - Set user-agent (with [`-A`](https://curl.se/docs/manpage.html#-A)) to
   one like a recent popular browser does

 - Set referer (with [`-E`](https://curl.se/docs/manpage.html#-E)) like
   it is set by the browser

 - If you use POST, make sure you send all the fields and in the same order as
   the browser does it.

## Check what the browsers do

 A good helper to make sure you do this right, is the web browsers' developers
 tools that let you view all headers you send and receive (even when using
 HTTPS).

 A more raw approach is to capture the HTTP traffic on the network with tools
 such as Wireshark or tcpdump and check what headers that were sent and
 received by the browser. (HTTPS forces you to use `SSLKEYLOGFILE` to do
 that.)
HTTP/2 with curl
================

[HTTP/2 Spec](https://www.rfc-editor.org/rfc/rfc7540.txt)
[http2 explained](https://daniel.haxx.se/http2/)

Build prerequisites
-------------------
  - nghttp2
  - OpenSSL, libressl, BoringSSL, NSS, GnuTLS, mbedTLS, wolfSSL or Schannel
    with a new enough version.

[nghttp2](https://nghttp2.org/)
-------------------------------

libcurl uses this 3rd party library for the low level protocol handling
parts. The reason for this is that HTTP/2 is much more complex at that layer
than HTTP/1.1 (which we implement on our own) and that nghttp2 is an already
existing and well functional library.

We require at least version 1.12.0.

Over an http:// URL
-------------------

If `CURLOPT_HTTP_VERSION` is set to `CURL_HTTP_VERSION_2_0`, libcurl will
include an upgrade header in the initial request to the host to allow
upgrading to HTTP/2.

Possibly we can later introduce an option that will cause libcurl to fail if
not possible to upgrade. Possibly we introduce an option that makes libcurl
use HTTP/2 at once over http://

Over an https:// URL
--------------------

If `CURLOPT_HTTP_VERSION` is set to `CURL_HTTP_VERSION_2_0`, libcurl will use
ALPN (or NPN) to negotiate which protocol to continue with. Possibly introduce
an option that will cause libcurl to fail if not possible to use HTTP/2.

`CURL_HTTP_VERSION_2TLS` was added in 7.47.0 as a way to ask libcurl to prefer
HTTP/2 for HTTPS but stick to 1.1 by default for plain old HTTP connections.

ALPN is the TLS extension that HTTP/2 is expected to use. The NPN extension is
for a similar purpose, was made prior to ALPN and is used for SPDY so early
HTTP/2 servers are implemented using NPN before ALPN support is widespread.

`CURLOPT_SSL_ENABLE_ALPN` and `CURLOPT_SSL_ENABLE_NPN` are offered to allow
applications to explicitly disable ALPN or NPN.

SSL libs
--------

The challenge is the ALPN and NPN support and all our different SSL
backends. You may need a fairly updated SSL library version for it to provide
the necessary TLS features. Right now we support:

  - OpenSSL:          ALPN and NPN
  - libressl:         ALPN and NPN
  - BoringSSL:        ALPN and NPN
  - NSS:              ALPN and NPN
  - GnuTLS:           ALPN
  - mbedTLS:          ALPN
  - Schannel:         ALPN
  - wolfSSL:          ALPN
  - Secure Transport: ALPN

Multiplexing
------------

Starting in 7.43.0, libcurl fully supports HTTP/2 multiplexing, which is the
term for doing multiple independent transfers over the same physical TCP
connection.

To take advantage of multiplexing, you need to use the multi interface and set
`CURLMOPT_PIPELINING` to `CURLPIPE_MULTIPLEX`. With that bit set, libcurl will
attempt to re-use existing HTTP/2 connections and just add a new stream over
that when doing subsequent parallel requests.

While libcurl sets up a connection to a HTTP server there is a period during
which it does not know if it can pipeline or do multiplexing and if you add new
transfers in that period, libcurl will default to start new connections for
those transfers. With the new option `CURLOPT_PIPEWAIT` (added in 7.43.0), you
can ask that a transfer should rather wait and see in case there's a
connection for the same host in progress that might end up being possible to
multiplex on. It favours keeping the number of connections low to the cost of
slightly longer time to first byte transferred.

Applications
------------

We hide HTTP/2's binary nature and convert received HTTP/2 traffic to headers
in HTTP 1.1 style. This allows applications to work unmodified.

curl tool
---------

curl offers the `--http2` command line option to enable use of HTTP/2.

curl offers the `--http2-prior-knowledge` command line option to enable use of
HTTP/2 without HTTP/1.1 Upgrade.

Since 7.47.0, the curl tool enables HTTP/2 by default for HTTPS connections.

curl tool limitations
---------------------

The command line tool does not support HTTP/2 server push. It supports
multiplexing when the parallel transfer option is used.

HTTP Alternative Services
-------------------------

Alt-Svc is an extension with a corresponding frame (ALTSVC) in HTTP/2 that
tells the client about an alternative "route" to the same content for the same
origin server that you get the response from. A browser or long-living client
can use that hint to create a new connection asynchronously. For libcurl, we
may introduce a way to bring such clues to the application and/or let a
subsequent request use the alternate route automatically.

[Detailed in RFC 7838](https://datatracker.ietf.org/doc/html/rfc7838)
                                  _   _ ____  _
                              ___| | | |  _ \| |
                             / __| | | | |_) | |
                            | (__| |_| |  _ <| |___
                             \___|\___/|_| \_\_____|

# SSL problems

  First, let's establish that we often refer to TLS and SSL interchangeably as
  SSL here. The current protocol is called TLS, it was called SSL a long time
  ago.

  There are several known reasons why a connection that involves SSL might
  fail. This is a document that attempts to details the most common ones and
  how to mitigate them.

## CA certs

  CA certs are used to digitally verify the server's certificate. You need a
  "ca bundle" for this. See lots of more details on this in the SSLCERTS
  document.

## CA bundle missing intermediate certificates

  When using said CA bundle to verify a server cert, you will experience
  problems if your CA store does not contain the certificates for the
  intermediates if the server does not provide them.

  The TLS protocol mandates that the intermediate certificates are sent in the
  handshake, but as browsers have ways to survive or work around such
  omissions, missing intermediates in TLS handshakes still happen that
  browser-users will not notice.

  Browsers work around this problem in two ways: they cache intermediate
  certificates from previous transfers and some implement the TLS "AIA"
  extension that lets the client explicitly download such certificates on
  demand.

## Protocol version

  Some broken servers fail to support the protocol negotiation properly that
  SSL servers are supposed to handle. This may cause the connection to fail
  completely. Sometimes you may need to explicitly select a SSL version to use
  when connecting to make the connection succeed.

  An additional complication can be that modern SSL libraries sometimes are
  built with support for older SSL and TLS versions disabled!

  All versions of SSL and the TLS versions before 1.2 are considered insecure
  and should be avoided. Use TLS 1.2 or later.

## Ciphers

  Clients give servers a list of ciphers to select from. If the list does not
  include any ciphers the server wants/can use, the connection handshake
  fails.

  curl has recently disabled the user of a whole bunch of seriously insecure
  ciphers from its default set (slightly depending on SSL backend in use).

  You may have to explicitly provide an alternative list of ciphers for curl
  to use to allow the server to use a WEAK cipher for you.

  Note that these weak ciphers are identified as flawed. For example, this
  includes symmetric ciphers with less than 128 bit keys and RC4.

  Schannel in Windows XP is not able to connect to servers that no longer
  support the legacy handshakes and algorithms used by those versions, so we
  advice against building curl to use Schannel on really old Windows versions.

  References:

  https://datatracker.ietf.org/doc/html/draft-popov-tls-prohibiting-rc4-01

## Allow BEAST

  BEAST is the name of a TLS 1.0 attack that surfaced 2011. When adding means
  to mitigate this attack, it turned out that some broken servers out there in
  the wild did not work properly with the BEAST mitigation in place.

  To make such broken servers work, the --ssl-allow-beast option was
  introduced. Exactly as it sounds, it re-introduces the BEAST vulnerability
  but on the other hand it allows curl to connect to that kind of strange
  servers.

## Disabling certificate revocation checks

  Some SSL backends may do certificate revocation checks (CRL, OCSP, etc)
  depending on the OS or build configuration. The --ssl-no-revoke option was
  introduced in 7.44.0 to disable revocation checking but currently is only
  supported for Schannel (the native Windows SSL library), with an exception
  in the case of Windows' Untrusted Publishers block list which it seems cannot
  be bypassed. This option may have broader support to accommodate other SSL
  backends in the future.

  References:

  https://curl.se/docs/ssl-compared.html
# Adding a new protocol?

Every once in a while someone comes up with the idea of adding support for yet
another protocol to curl. After all, curl already supports 25 something
protocols and it is the Internet transfer machine for the world.

In the curl project we love protocols and we love supporting many protocols
and do it well.

So how do you proceed to add a new protocol and what are the requirements?

## No fixed set of requirements

This document is an attempt to describe things to consider. There is no
checklist of the twenty-seven things you need to cross off. We view the entire
effort as a whole and then judge if it seems to be the right thing - for
now. The more things that look right, fit our patterns and are done in ways
that align with our thinking, the better are the chances that we will agree
that supporting this protocol is a grand idea.

## Mutual benefit is preferred

curl is not here for your protocol. Your protocol is not here for curl. The
best cooperation and end result occur when all involved parties mutually see
and agree that supporting this protocol in curl would be good for everyone.
Heck, for the world.

Consider "selling us" the idea that we need an implementation merged in curl,
to be fairly important. *Why* do we want curl to support this new protocol?

## Protocol requirements

### Client-side

The protocol implementation is for a client's side of a "communication
session".

### Transfer oriented

The protocol itself should be focused on *transfers*. Be it uploads or
downloads or both. It should at least be possible to view the transfers as
such, like we can view reading emails over POP3 as a download and sending
emails over SMTP as an upload.

If you cannot even shoehorn the protocol into a transfer focused view, then
you are up for a tough argument.

### URL

There should be a documented URL format. If there is an RFC for it there is no
question about it but the syntax does not have to be a published RFC. It could
be enough if it is already in use by other implementations.

If you make up the syntax just in order to be able to propose it to curl, then
you are in a bad place. URLs are designed and defined for interoperability.
There should at least be a good chance that other clients and servers can be
implemented supporting the same URL syntax and work the same or similar way.

URLs work on registered 'schemes'. There is a register of [all officially
recognized
schemes](https://www.iana.org/assignments/uri-schemes/uri-schemes.xhtml). If
your protocol is not in there, is it really a protocol we want?

### Wide and public use

The protocol shall already be used or have an expectation of getting used
widely. Experimental protocols are better off worked on in experiments first,
to prove themselves before they are adopted by curl.

## Code

Of course the code needs to be written, provided, licensed agreeably and it
should follow our code guidelines and review comments have to be dealt with.
If the implementation needs third party code, that third party code should not
have noticeably lesser standards than the curl project itself.

## Tests

As much of the protocol implementation as possible needs to be verified by
curl test cases. We must have the implementation get tested by CI jobs,
torture tests and more.

We have experienced many times in the past how new implementations were brought
to curl and immediately once the code had been merged, the originator vanished
from the face of the earth. That is fine, but we need to take the necessary
precautions so when it happens we are still fine.

Our test infrastructure is powerful enough to test just about every possible
protocol - but it might require a bit of an effort to make it happen.

## Documentation

We cannot assume that users are particularly familiar with details and
peculiarities of the protocol. It needs documentation.

Maybe it even needs some internal documentation so that the developers who
will try to debug something five years from now can figure out functionality a
little easier!

The protocol specification itself should be freely available without requiring
any NDA or similar.

## Do not compare

We are constantly raising the bar and we are constantly improving the
project. A lot of things we did in the past would not be acceptable if done
today. Therefore, you might be tempted to use shortcuts or "hacks" you can
spot other - existing - protocol implementations have used, but there is
nothing to gain from that. The bar has been raised. Former "cheats" will not be
tolerated anymore.
# BUGS

## There are still bugs

 Curl and libcurl keep being developed. Adding features and changing code
 means that bugs will sneak in, no matter how hard we try not to.

 Of course there are lots of bugs left. And lots of misfeatures.

 To help us make curl the stable and solid product we want it to be, we need
 bug reports and bug fixes.

## Where to report

 If you cannot fix a bug yourself and submit a fix for it, try to report an as
 detailed report as possible to a curl mailing list to allow one of us to have
 a go at a solution. You can optionally also submit your problem in [curl's
 bug tracking system](https://github.com/curl/curl/issues).

 Please read the rest of this document below first before doing that.

 If you feel you need to ask around first, find a suitable [mailing list](
 https://curl.se/mail/) and post your questions there.

## Security bugs

 If you find a bug or problem in curl or libcurl that you think has a security
 impact, for example a bug that can put users in danger or make them
 vulnerable if the bug becomes public knowledge, then please report that bug
 using our security development process.

 Security related bugs or bugs that are suspected to have a security impact,
 should be reported on the [curl security tracker at
 HackerOne](https://hackerone.com/curl).

 This ensures that the report reaches the curl security team so that they
 first can deal with the report away from the public to minimize the harm
 and impact it will have on existing users out there who might be using the
 vulnerable versions.

 The curl project's process for handling security related issues is
 [documented separately](https://curl.se/dev/secprocess.html).

## What to report

 When reporting a bug, you should include all information that will help us
 understand what's wrong, what you expected to happen and how to repeat the
 bad behavior. You therefore need to tell us:

 - your operating system's name and version number

 - what version of curl you are using (`curl -V` is fine)

 - versions of the used libraries that libcurl is built to use

 - what URL you were working with (if possible), at least which protocol

 and anything and everything else you think matters. Tell us what you expected
 to happen, tell use what did happen, tell us how you could make it work
 another way. Dig around, try out, test. Then include all the tiny bits and
 pieces in your report. You will benefit from this yourself, as it will enable
 us to help you quicker and more accurately.

 Since curl deals with networks, it often helps us if you include a protocol
 debug dump with your bug report. The output you get by using the `-v` or
 `--trace` options.

 If curl crashed, causing a core dump (in unix), there is hardly any use to
 send that huge file to anyone of us. Unless we have the same system setup as
 you, we cannot do much with it. Instead, we ask you to get a stack trace and
 send that (much smaller) output to us instead.

 The address and how to subscribe to the mailing lists are detailed in the
 `MANUAL.md` file.

## libcurl problems

 When you have written your own application with libcurl to perform transfers,
 it is even more important to be specific and detailed when reporting bugs.

 Tell us the libcurl version and your operating system. Tell us the name and
 version of all relevant sub-components like for example the SSL library
 you are using and what name resolving your libcurl uses. If you use SFTP or
 SCP, the libssh2 version is relevant etc.

 Showing us a real source code example repeating your problem is the best way
 to get our attention and it will greatly increase our chances to understand
 your problem and to work on a fix (if we agree it truly is a problem).

 Lots of problems that appear to be libcurl problems are actually just abuses
 of the libcurl API or other malfunctions in your applications. It is advised
 that you run your problematic program using a memory debug tool like valgrind
 or similar before you post memory-related or "crashing" problems to us.

## Who will fix the problems

 If the problems or bugs you describe are considered to be bugs, we want to
 have the problems fixed.

 There are no developers in the curl project that are paid to work on bugs.
 All developers that take on reported bugs do this on a voluntary basis. We do
 it out of an ambition to keep curl and libcurl excellent products and out of
 pride.

 Please do not assume that you can just lump over something to us and it will
 then magically be fixed after some given time. Most often we need feedback
 and help to understand what you have experienced and how to repeat a
 problem. Then we may only be able to assist YOU to debug the problem and to
 track down the proper fix.

 We get reports from many people every month and each report can take a
 considerable amount of time to really go to the bottom with.

## How to get a stack trace

 First, you must make sure that you compile all sources with `-g` and that you
 do not 'strip' the final executable. Try to avoid optimizing the code as well,
 remove `-O`, `-O2` etc from the compiler options.

 Run the program until it cores.

 Run your debugger on the core file, like `<debugger> curl
 core`. `<debugger>` should be replaced with the name of your debugger, in
 most cases that will be `gdb`, but `dbx` and others also occur.

 When the debugger has finished loading the core file and presents you a
 prompt, enter `where` (without quotes) and press return.

 The list that is presented is the stack trace. If everything worked, it is
 supposed to contain the chain of functions that were called when curl
 crashed. Include the stack trace with your detailed bug report. it will help a
 lot.

## Bugs in libcurl bindings

 There will of course pop up bugs in libcurl bindings. You should then
 primarily approach the team that works on that particular binding and see
 what you can do to help them fix the problem.

 If you suspect that the problem exists in the underlying libcurl, then please
 convert your program over to plain C and follow the steps outlined above.

## Bugs in old versions

 The curl project typically releases new versions every other month, and we
 fix several hundred bugs per year. For a huge table of releases, number of
 bug fixes and more, see: https://curl.se/docs/releases.html

 The developers in the curl project do not have bandwidth or energy enough to
 maintain several branches or to spend much time on hunting down problems in
 old versions when chances are we already fixed them or at least that they have
 changed nature and appearance in later versions.

 When you experience a problem and want to report it, you really SHOULD
 include the version number of the curl you are using when you experience the
 issue. If that version number shows us that you are using an out-of-date curl,
 you should also try out a modern curl version to see if the problem persists
 or how/if it has changed in appearance.

 Even if you cannot immediately upgrade your application/system to run the
 latest curl version, you can most often at least run a test version or
 experimental build or similar, to get this confirmed or not.

 At times people insist that they cannot upgrade to a modern curl version, but
 instead they "just want the bug fixed". That is fine, just do not count on us
 spending many cycles on trying to identify which single commit, if that is
 even possible, that at some point in the past fixed the problem you are now
 experiencing.

 Security wise, it is almost always a bad idea to lag behind the current curl
 versions by a lot. We keep discovering and reporting security problems
 over time see you can see in [this
 table](https://curl.se/docs/vulnerabilities.html)

# Bug fixing procedure

## What happens on first filing

 When a new issue is posted in the issue tracker or on the mailing list, the
 team of developers first need to see the report. Maybe they took the day off,
 maybe they are off in the woods hunting. Have patience. Allow at least a few
 days before expecting someone to have responded.

 In the issue tracker you can expect that some labels will be set on the issue
 to help categorize it.

## First response

 If your issue/bug report was not perfect at once (and few are), chances are
 that someone will ask follow-up questions. Which version did you use? Which
 options did you use? How often does the problem occur? How can we reproduce
 this problem? Which protocols does it involve? Or perhaps much more specific
 and deep diving questions. It all depends on your specific issue.

 You should then respond to these follow-up questions and provide more info
 about the problem, so that we can help you figure it out. Or maybe you can
 help us figure it out. An active back-and-forth communication is important
 and the key for finding a cure and landing a fix.

## Not reproducible

 For problems that we cannot reproduce and cannot understand even after having
 gotten all the info we need and having studied the source code over again,
 are really hard to solve so then we may require further work from you who
 actually see or experience the problem.

## Unresponsive

 If the problem have not been understood or reproduced, and there's nobody
 responding to follow-up questions or questions asking for clarifications or
 for discussing possible ways to move forward with the task, we take that as a
 strong suggestion that the bug is unimportant.

 Unimportant issues will be closed as inactive sooner or later as they cannot
 be fixed. The inactivity period (waiting for responses) should not be shorter
 than two weeks but may extend months.

## Lack of time/interest

 Bugs that are filed and are understood can unfortunately end up in the
 "nobody cares enough about it to work on it" category. Such bugs are
 perfectly valid problems that *should* get fixed but apparently are not. We
 try to mark such bugs as `KNOWN_BUGS material` after a time of inactivity and
 if no activity is noticed after yet some time those bugs are added to the
 `KNOWN_BUGS` document and are closed in the issue tracker.

## `KNOWN_BUGS`

 This is a list of known bugs. Bugs we know exist and that have been pointed
 out but that have not yet been fixed. The reasons for why they have not been
 fixed can involve anything really, but the primary reason is that nobody has
 considered these problems to be important enough to spend the necessary time
 and effort to have them fixed.

 The `KNOWN_BUGS` items are always up for grabs and we love the ones who bring
 one of them back to life and offer solutions to them.

 The `KNOWN_BUGS` document has a sibling document known as `TODO`.

## `TODO`

 Issues that are filed or reported that are not really bugs but more missing
 features or ideas for future improvements and so on are marked as
 'enhancement' or 'feature-request' and will be added to the `TODO` document
 and the issues are closed. We do not keep TODO items open in the issue
 tracker.

 The `TODO` document is full of ideas and suggestions of what we can add or
 fix one day. you are always encouraged and free to grab one of those items and
 take up a discussion with the curl development team on how that could be
 implemented or provided in the project so that you can work on ticking it odd
 that document.

 If an issue is rather a bug and not a missing feature or functionality, it is
 listed in `KNOWN_BUGS` instead.

## Closing off stalled bugs

 The [issue and pull request trackers](https://github.com/curl/curl) only
 holds "active" entries open (using a non-precise definition of what active
 actually is, but they are at least not completely dead). Those that are
 abandoned or in other ways dormant will be closed and sometimes added to
 `TODO` and `KNOWN_BUGS` instead.

 This way, we only have "active" issues open on GitHub. Irrelevant issues and
 pull requests will not distract developers or casual visitors.
# MQTT in curl

## Usage

A plain "GET" subscribes to the topic and prints all published messages.
Doing a "POST" publishes the post data to the topic and exits.

Example subscribe:

    curl mqtt://host/home/bedroom/temp

Example publish:

    curl -d 75 mqtt://host/home/bedroom/dimmer

## What does curl deliver as a response to a subscribe

It outputs two bytes topic length (MSB | LSB), the topic followed by the
payload.

## Caveats

Remaining limitations:
 - Only QoS level 0 is implemented for publish
 - No way to set retain flag for publish
 - No TLS (mqtts) support
 - Naive EAGAIN handling will not handle split messages
# Parallel transfers

curl 7.66.0 introduces support for doing multiple transfers simultaneously; in
parallel.

## -Z, --parallel

When this command line option is used, curl will perform the transfers given
to it at the same time. It will do up to `--parallel-max` concurrent
transfers, with a default value of 50.

## Progress meter

The progress meter that is displayed when doing parallel transfers is
completely different than the regular one used for each single transfer.

  It shows:

 o percent download (if known, which means *all* transfers need to have a
   known size)
 o percent upload (if known, with the same caveat as for download)
 o total amount of downloaded data
 o total amount of uploaded data
 o number of transfers to perform
 o number of concurrent transfers being transferred right now
 o number of transfers queued up waiting to start
 o total time all transfers are expected to take (if sizes are known)
 o current time the transfers have spent so far
 o estimated time left (if sizes are known)
 o current transfer speed (the faster of UL/DL speeds measured over the last
   few seconds)

Example:

    DL% UL%  Dled  Uled  Xfers  Live   Qd Total     Current  Left    Speed
    72  --  37.9G     0   101    30    23  0:00:55  0:00:34  0:00:22 2752M

## Behavior differences

Connections are shared fine between different easy handles, but the
"authentication contexts" are not. So for example doing HTTP Digest auth with
one handle for a particular transfer and then continue on with another handle
that reuses the same connection, the second handle cannot send the necessary
Authorization header at once since the context is only kept in the original
easy handle.

To fix this, the authorization state could be made possible to share with the
share API as well, as a context per origin + path (realm?) basically.

Visible in test 153, 1412 and more.

## Feedback

This is early days for parallel transfer support. Keep your eyes open for
unintended side effects or downright bugs.

Tell us what you think and how you think we could improve this feature!

# URL syntax and their use in curl

## Specifications

The official "URL syntax" is primarily defined in these two different
specifications:

 - [RFC 3986](https://datatracker.ietf.org/doc/html/rfc3986) (although URL is called
   "URI" in there)
 - [The WHATWG URL Specification](https://url.spec.whatwg.org/)

RFC 3986 is the earlier one, and curl has always tried to adhere to that one
(since it shipped in January 2005).

The WHATWG URL spec was written later, is incompatible with the RFC 3986 and
changes over time.

## Variations

URL parsers as implemented in browsers, libraries and tools usually opt to
support one of the mentioned specifications. Bugs, differences in
interpretations and the moving nature of the WHATWG spec does however make it
unlikely that multiple parsers treat URLs the same way.

## Security

Due to the inherent differences between URL parser implementations, it is
considered a security risk to mix different implementations and assume the
same behavior!

For example, if you use one parser to check if a URL uses a good host name or
the correct auth field, and then pass on that same URL to a *second* parser,
there will always be a risk it treats the same URL differently. There is no
right and wrong in URL land, only differences of opinions.

libcurl offers a separate API to its URL parser for this reason, among others.

Applications may at times find it convenient to allow users to specify URLs
for various purposes and that string would then end up fed to curl. Getting a
URL from an external untrusted party and using it with curl brings several
security concerns:

1. If you have an application that runs as or in a server application, getting
   an unfiltered URL can trick your application to access a local resource
   instead of a remote resource. Protecting yourself against localhost accesses
   is hard when accepting user provided URLs.

2. Such custom URLs can access other ports than you planned as port numbers
   are part of the regular URL format. The combination of a local host and a
   custom port number can allow external users to play tricks with your local
   services.

3. Such a URL might use other schemes than you thought of or planned for.

## "RFC3986 plus"

curl recognizes a URL syntax that we call "RFC 3986 plus". It is grounded on
the well established RFC 3986 to make sure previously written command lines and
curl using scripts will remain working.

curl's URL parser allows a few deviations from the spec in order to
inter-operate better with URLs that appear in the wild.

### spaces

A URL provided to curl cannot contain spaces. They need to be provided URL
encoded to be accepted in a URL by curl.

An exception to this rule: `Location:` response headers that indicate to a
client where a resource has been redirected to, sometimes contain spaces. This
is a violation of RFC 3986 but is fine in the WHATWG spec. curl handles these
by re-encoding them to `%20`.

### non-ASCII

Byte values in a provided URL that are outside of the printable ASCII range
are percent-encoded by curl.

### multiple slashes

An absolute URL always starts with a "scheme" followed by a colon. For all the
schemes curl supports, the colon must be followed by two slashes according to
RFC 3986 but not according to the WHATWG spec - which allows one to infinity
amount.

curl allows one, two or three slashes after the colon to still be considered a
valid URL.

### "scheme-less"

curl supports "URLs" that do not start with a scheme. This is not supported by
any of the specifications. This is a shortcut to entering URLs that was
supported by browsers early on and has been mimicked by curl.

Based on what the host name starts with, curl will "guess" what protocol to
use:

 - `ftp.` means FTP
 - `dict.` means DICT
 - `ldap.` means LDAP
 - `imap.` means IMAP
 - `smtp.` means SMTP
 - `pop3.` means POP3
 - all other means HTTP

### globbing letters

The curl command line tool supports "globbing" of URLs. It means that you can
create ranges and lists using `[N-M]` and `{one,two,three}` sequences. The
letters used for this (`[]{}`) are reserved in RFC 3986 and can therefore not
legitimately be part of such a URL.

They are however not reserved or special in the WHATWG specification, so
globbing can mess up such URLs. Globbing can be turned off for such occasions
(using `--globoff`).

# URL syntax details

A URL may consist of the following components - many of them are optional:

    [scheme][divider][userinfo][hostname][port number][path][query][fragment]

Each component is separated from the following component with a divider
character or string.

For example, this could look like:

    http://user:password@www.example.com:80/index.hmtl?foo=bar#top

## Scheme

The scheme specifies the protocol to use. A curl build can support a few or
many different schemes. You can limit what schemes curl should accept.

curl supports the following schemes on URLs specified to transfer. They are
matched case insensitively:

`dict`, `file`, `ftp`, `ftps`, `gopher`, `gophers`, `http`, `https`, `imap`,
`imaps`, `ldap`, `ldaps`, `mqtt`, `pop3`, `pop3s`, `rtmp`, `rtmpe`, `rtmps`,
`rtmpt`, `rtmpte`, `rtmpts`, `rtsp`, `smb`, `smbs`, `smtp`, `smtps`, `telnet`,
`tftp`

When the URL is specified to identify a proxy, curl recognizes the following
schemes:

`http`, `https`, `socks4`, `socks4a`, `socks5`, `socks5h`, `socks`

## Userinfo

The userinfo field can be used to set user name and password for
authentication purposes in this transfer. The use of this field is discouraged
since it often means passing around the password in plain text and is thus a
security risk.

URLs for IMAP, POP3 and SMTP also support *login options* as part of the
userinfo field. they are provided as a semicolon after the password and then
the options.

## Hostname

The hostname part of the URL contains the address of the server that you want
to connect to. This can be the fully qualified domain name of the server, the
local network name of the machine on your network or the IP address of the
server or machine represented by either an IPv4 or IPv6 address (within
brackets). For example:

    http://www.example.com/

    http://hostname/

    http://192.168.0.1/

    http://[2001:1890:1112:1::20]/

### "localhost"

Starting in curl 7.77.0, curl uses loopback IP addresses for the name
`localhost`: `127.0.0.1` and `::1`. It does not resolve the name using the
resolver functions.

This is done to make sure the host accessed is truly the localhost - the local
machine.

### IDNA

If curl was built with International Domain Name (IDN) support, it can also
handle host names using non-ASCII characters.

When built with libidn2, curl uses the IDNA 2008 standard. This is equivalent
to the WHATWG URL spec, but differs from certain browsers that use IDNA 2003
Transitional Processing. The two standards have a huge overlap but differ
slightly, perhaps most famously in how they deal with the German "double s"
(`ß`).

When winidn is used, curl uses IDNA 2003 Transitional Processing, like the rest
of Windows.

## Port number

If there's a colon after the hostname, that should be followed by the port
number to use. 1 - 65535. curl also supports a blank port number field - but
only if the URL starts with a scheme.

If the port number is not specified in the URL, curl will used a default port
based on the provide scheme:

DICT 2628, FTP 21, FTPS 990, GOPHER 70, GOPHERS 70, HTTP 80, HTTPS 443,
IMAP 132, IMAPS 993, LDAP 369, LDAPS 636, MQTT 1883, POP3 110, POP3S 995,
RTMP 1935, RTMPS 443, RTMPT 80, RTSP 554, SCP 22, SFTP 22, SMB 445, SMBS 445,
SMTP 25, SMTPS 465, TELNET 23, TFTP 69

# Scheme specific behaviors

## FTP

The path part of an FTP request specifies the file to retrieve and from which
directory. If the file part is omitted then libcurl downloads the directory
listing for the directory specified. If the directory is omitted then the
directory listing for the root / home directory will be returned.

FTP servers typically put the user in its "home directory" after login, which
then differs between users. To explicitly specify the root directory of an FTP
server start the path with double slash `//` or `/%2f` (2F is the hexadecimal
value of the ascii code for the slash).

## FILE

When a `FILE://` URL is accessed on Windows systems, it can be crafted in a
way so that Windows attempts to connect to a (remote) machine when curl wants
to read or write such a path.

curl only allows the hostname part of a FILE URL to be one out of these three
alternatives: `localhost`, `127.0.0.1` or blank ("", zero characters).
Anything else will make curl fail to parse the URL.

### Windows-specific FILE details

curl accepts that the FILE URL's path starts with a "drive letter". That is a
single letter `a` to `z` followed by a colon or a pipe character (`|`).

The Windows operating system itself will convert some file accesses to perform
network accesses over SMB/CIFS, through several different file path patterns.
This way, a `file://` URL passed to curl *might* be converted into a network
access inadvertently and unknowingly to curl. This is a Windows feature curl
cannot control or disable.

## IMAP

The path part of an IMAP request not only specifies the mailbox to list or
select, but can also be used to check the `UIDVALIDITY` of the mailbox, to
specify the `UID`, `SECTION` and `PARTIAL` octets of the message to fetch and
to specify what messages to search for.

A top level folder list:

    imap://user:password@mail.example.com

A folder list on the user's inbox:

    imap://user:password@mail.example.com/INBOX

Select the user's inbox and fetch message with uid = 1:

    imap://user:password@mail.example.com/INBOX/;UID=1

Select the user's inbox and fetch the first message in the mail box:

    imap://user:password@mail.example.com/INBOX/;MAILINDEX=1

Select the user's inbox, check the `UIDVALIDITY` of the mailbox is 50 and
fetch message 2 if it is:

    imap://user:password@mail.example.com/INBOX;UIDVALIDITY=50/;UID=2

Select the user's inbox and fetch the text portion of message 3:

    imap://user:password@mail.example.com/INBOX/;UID=3/;SECTION=TEXT

Select the user's inbox and fetch the first 1024 octets of message 4:

    imap://user:password@mail.example.com/INBOX/;UID=4/;PARTIAL=0.1024

Select the user's inbox and check for NEW messages:

    imap://user:password@mail.example.com/INBOX?NEW

Select the user's inbox and search for messages containing "shadows" in the
subject line:

    imap://user:password@mail.example.com/INBOX?SUBJECT%20shadows

Searching via the query part of the URL `?` is a search request for the results
to be returned as message sequence numbers (MAILINDEX). It is possible to make
a search request for results to be returned as unique ID numbers (UID) by using
a custom curl request via `-X`. UID numbers are unique per session (and
multiple sessions when UIDVALIDITY is the same). For example, if you are
searching for `"foo bar"` in header+body (TEXT) and you want the matching
MAILINDEX numbers returned then you could search via URL:

    imap://user:password@mail.example.com/INBOX?TEXT%20%22foo%20bar%22

.. but if you wanted matching UID numbers you would have to use a custom request:

    imap://user:password@mail.example.com/INBOX -X "UID SEARCH TEXT \"foo bar\""

For more information about IMAP commands please see RFC 9051. For more
information about the individual components of an IMAP URL please see RFC 5092.

* Note old curl versions would FETCH by message sequence number when UID was
specified in the URL. That was a bug fixed in 7.62.0, which added MAILINDEX to
FETCH by mail sequence number.

## LDAP

The path part of a LDAP request can be used to specify the: Distinguished
Name, Attributes, Scope, Filter and Extension for a LDAP search. Each field is
separated by a question mark and when that field is not required an empty
string with the question mark separator should be included.

Search for the DN as `My Organisation`:

    ldap://ldap.example.com/o=My%20Organisation

the same search but will only return postalAddress attributes:

    ldap://ldap.example.com/o=My%20Organisation?postalAddress

Search for an empty DN and request information about the
`rootDomainNamingContext` attribute for an Active Directory server:

    ldap://ldap.example.com/?rootDomainNamingContext

For more information about the individual components of a LDAP URL please
see [RFC 4516](https://datatracker.ietf.org/doc/html/rfc4516).

## POP3

The path part of a POP3 request specifies the message ID to retrieve. If the
ID is not specified then a list of waiting messages is returned instead.

## SCP

The path part of an SCP URL specifies the path and file to retrieve or
upload. The file is taken as an absolute path from the root directory on the
server.

To specify a path relative to the user's home directory on the server, prepend
`~/` to the path portion.

## SFTP

The path part of an SFTP URL specifies the file to retrieve or upload. If the
path ends with a slash (`/`) then a directory listing is returned instead of a
file. If the path is omitted entirely then the directory listing for the root
/ home directory will be returned.

## SMB
The path part of a SMB request specifies the file to retrieve and from what
share and directory or the share to upload to and as such, may not be omitted.
If the user name is embedded in the URL then it must contain the domain name
and as such, the backslash must be URL encoded as %2f.

curl supports SMB version 1 (only)

## SMTP

The path part of a SMTP request specifies the host name to present during
communication with the mail server. If the path is omitted, then libcurl will
attempt to resolve the local computer's host name. However, this may not
return the fully qualified domain name that is required by some mail servers
and specifying this path allows you to set an alternative name, such as your
machine's fully qualified domain name, which you might have obtained from an
external function such as gethostname or getaddrinfo.

The default smtp port is 25. Some servers use port 587 as an alternative.

## RTMP

There's no official URL spec for RTMP so libcurl uses the URL syntax supported
by the underlying librtmp library. It has a syntax where it wants a
traditional URL, followed by a space and a series of space-separated
`name=value` pairs.

While space is not typically a "legal" letter, libcurl accepts them. When a
user wants to pass in a `#` (hash) character it will be treated as a fragment
and get cut off by libcurl if provided literally. You will instead have to
escape it by providing it as backslash and its ASCII value in hexadecimal:
`\23`.
curl internals
==============

 - [Intro](#intro)
 - [git](#git)
 - [Portability](#Portability)
 - [Windows vs Unix](#winvsunix)
 - [Library](#Library)
   - [`Curl_connect`](#Curl_connect)
   - [`multi_do`](#multi_do)
   - [`Curl_readwrite`](#Curl_readwrite)
   - [`multi_done`](#multi_done)
   - [`Curl_disconnect`](#Curl_disconnect)
 - [HTTP(S)](#http)
 - [FTP](#ftp)
 - [Kerberos](#kerberos)
 - [TELNET](#telnet)
 - [FILE](#file)
 - [SMB](#smb)
 - [LDAP](#ldap)
 - [Email](#email)
 - [General](#general)
 - [Persistent Connections](#persistent)
 - [multi interface/non-blocking](#multi)
 - [SSL libraries](#ssl)
 - [Library Symbols](#symbols)
 - [Return Codes and Informationals](#returncodes)
 - [AP/ABI](#abi)
 - [Client](#client)
 - [Memory Debugging](#memorydebug)
 - [Test Suite](#test)
 - [Asynchronous name resolves](#asyncdns)
   - [c-ares](#cares)
 - [`curl_off_t`](#curl_off_t)
 - [curlx](#curlx)
 - [Content Encoding](#contentencoding)
 - [`hostip.c` explained](#hostip)
 - [Track Down Memory Leaks](#memoryleak)
 - [`multi_socket`](#multi_socket)
 - [Structs in libcurl](#structs)
   - [Curl_easy](#Curl_easy)
   - [connectdata](#connectdata)
   - [Curl_multi](#Curl_multi)
   - [Curl_handler](#Curl_handler)
   - [conncache](#conncache)
   - [Curl_share](#Curl_share)
   - [CookieInfo](#CookieInfo)

<a name="intro"></a>
Intro
=====

 This project is split in two. The library and the client. The client part
 uses the library, but the library is designed to allow other applications to
 use it.

 The largest amount of code and complexity is in the library part.


<a name="git"></a>
git
===

 All changes to the sources are committed to the git repository as soon as
 they are somewhat verified to work. Changes shall be committed as independently
 as possible so that individual changes can be easily spotted and tracked
 afterwards.

 Tagging shall be used extensively, and by the time we release new archives we
 should tag the sources with a name similar to the released version number.

<a name="Portability"></a>
Portability
===========

 We write curl and libcurl to compile with C89 compilers. On 32-bit and up
 machines. Most of libcurl assumes more or less POSIX compliance but that is
 not a requirement.

 We write libcurl to build and work with lots of third party tools, and we
 want it to remain functional and buildable with these and later versions
 (older versions may still work but is not what we work hard to maintain):

Dependencies
------------

 - OpenSSL      0.9.7
 - GnuTLS       3.1.10
 - zlib         1.1.4
 - libssh2      1.0
 - c-ares       1.16.0
 - libidn2      2.0.0
 - wolfSSL      2.0.0
 - openldap     2.0
 - MIT Kerberos 1.2.4
 - GSKit        V5R3M0
 - NSS          3.14.x
 - Heimdal      ?
 - nghttp2      1.12.0
 - WinSock      2.2 (on Windows 95+ and Windows CE .NET 4.1+)

Operating Systems
-----------------

 On systems where configure runs, we aim at working on them all - if they have
 a suitable C compiler. On systems that do not run configure, we strive to keep
 curl running correctly on:

 - Windows      98
 - AS/400       V5R3M0
 - Symbian      9.1
 - Windows CE   ?
 - TPF          ?

Build tools
-----------

 When writing code (mostly for generating stuff included in release tarballs)
 we use a few "build tools" and we make sure that we remain functional with
 these versions:

 - GNU Libtool  1.4.2
 - GNU Autoconf 2.57
 - GNU Automake 1.7
 - GNU M4       1.4
 - perl         5.004
 - roffit       0.5
 - groff        ? (any version that supports `groff -Tps -man [in] [out]`)
 - ps2pdf (gs)  ?

<a name="winvsunix"></a>
Windows vs Unix
===============

 There are a few differences in how to program curl the Unix way compared to
 the Windows way. Perhaps the four most notable details are:

 1. Different function names for socket operations.

   In curl, this is solved with defines and macros, so that the source looks
   the same in all places except for the header file that defines them. The
   macros in use are `sclose()`, `sread()` and `swrite()`.

 2. Windows requires a couple of init calls for the socket stuff.

   That is taken care of by the `curl_global_init()` call, but if other libs
   also do it etc there might be reasons for applications to alter that
   behavior.

   We require WinSock version 2.2 and load this version during global init.

 3. The file descriptors for network communication and file operations are
    not as easily interchangeable as in Unix.

   We avoid this by not trying any funny tricks on file descriptors.

 4. When writing data to stdout, Windows makes end-of-lines the DOS way, thus
    destroying binary data, although you do want that conversion if it is
    text coming through... (sigh)

   We set stdout to binary under windows

 Inside the source code, We make an effort to avoid `#ifdef [Your OS]`. All
 conditionals that deal with features *should* instead be in the format
 `#ifdef HAVE_THAT_WEIRD_FUNCTION`. Since Windows cannot run configure scripts,
 we maintain a `curl_config-win32.h` file in lib directory that is supposed to
 look exactly like a `curl_config.h` file would have looked like on a Windows
 machine.

 Generally speaking: always remember that this will be compiled on dozens of
 operating systems. Do not walk on the edge.

<a name="Library"></a>
Library
=======

 (See [Structs in libcurl](#structs) for the separate section describing all
 major internal structs and their purposes.)

 There are plenty of entry points to the library, namely each publicly defined
 function that libcurl offers to applications. All of those functions are
 rather small and easy-to-follow. All the ones prefixed with `curl_easy` are
 put in the `lib/easy.c` file.

 `curl_global_init()` and `curl_global_cleanup()` should be called by the
 application to initialize and clean up global stuff in the library. As of
 today, it can handle the global SSL initialization if SSL is enabled and it
 can initialize the socket layer on Windows machines. libcurl itself has no
 "global" scope.

 All printf()-style functions use the supplied clones in `lib/mprintf.c`. This
 makes sure we stay absolutely platform independent.

 [ `curl_easy_init()`][2] allocates an internal struct and makes some
 initializations. The returned handle does not reveal internals. This is the
 `Curl_easy` struct which works as an "anchor" struct for all `curl_easy`
 functions. All connections performed will get connect-specific data allocated
 that should be used for things related to particular connections/requests.

 [`curl_easy_setopt()`][1] takes three arguments, where the option stuff must
 be passed in pairs: the parameter-ID and the parameter-value. The list of
 options is documented in the man page. This function mainly sets things in
 the `Curl_easy` struct.

 `curl_easy_perform()` is just a wrapper function that makes use of the multi
 API. It basically calls `curl_multi_init()`, `curl_multi_add_handle()`,
 `curl_multi_wait()`, and `curl_multi_perform()` until the transfer is done
 and then returns.

 Some of the most important key functions in `url.c` are called from
 `multi.c` when certain key steps are to be made in the transfer operation.

<a name="Curl_connect"></a>
Curl_connect()
--------------

   Analyzes the URL, it separates the different components and connects to the
   remote host. This may involve using a proxy and/or using SSL. The
   `Curl_resolv()` function in `lib/hostip.c` is used for looking up host
   names (it does then use the proper underlying method, which may vary
   between platforms and builds).

   When `Curl_connect` is done, we are connected to the remote site. Then it
   is time to tell the server to get a document/file. `Curl_do()` arranges
   this.

   This function makes sure there's an allocated and initiated `connectdata`
   struct that is used for this particular connection only (although there may
   be several requests performed on the same connect). A bunch of things are
   initialized/inherited from the `Curl_easy` struct.

<a name="multi_do"></a>
multi_do()
---------

   `multi_do()` makes sure the proper protocol-specific function is called.
   The functions are named after the protocols they handle.

   The protocol-specific functions of course deal with protocol-specific
   negotiations and setup. When they are ready to start the actual file
   transfer they call the `Curl_setup_transfer()` function (in
   `lib/transfer.c`) to setup the transfer and returns.

   If this DO function fails and the connection is being re-used, libcurl will
   then close this connection, setup a new connection and re-issue the DO
   request on that. This is because there is no way to be perfectly sure that
   we have discovered a dead connection before the DO function and thus we
   might wrongly be re-using a connection that was closed by the remote peer.

<a name="Curl_readwrite"></a>
Curl_readwrite()
----------------

   Called during the transfer of the actual protocol payload.

   During transfer, the progress functions in `lib/progress.c` are called at
   frequent intervals (or at the user's choice, a specified callback might get
   called). The speedcheck functions in `lib/speedcheck.c` are also used to
   verify that the transfer is as fast as required.

<a name="multi_done"></a>
multi_done()
-----------

   Called after a transfer is done. This function takes care of everything
   that has to be done after a transfer. This function attempts to leave
   matters in a state so that `multi_do()` should be possible to call again on
   the same connection (in a persistent connection case). It might also soon
   be closed with `Curl_disconnect()`.

<a name="Curl_disconnect"></a>
Curl_disconnect()
-----------------

   When doing normal connections and transfers, no one ever tries to close any
   connections so this is not normally called when `curl_easy_perform()` is
   used. This function is only used when we are certain that no more transfers
   are going to be made on the connection. It can be also closed by force, or
   it can be called to make sure that libcurl does not keep too many
   connections alive at the same time.

   This function cleans up all resources that are associated with a single
   connection.

<a name="http"></a>
HTTP(S)
=======

 HTTP offers a lot and is the protocol in curl that uses the most lines of
 code. There is a special file `lib/formdata.c` that offers all the
 multipart post functions.

 base64-functions for user+password stuff (and more) is in `lib/base64.c`
 and all functions for parsing and sending cookies are found in
 `lib/cookie.c`.

 HTTPS uses in almost every case the same procedure as HTTP, with only two
 exceptions: the connect procedure is different and the function used to read
 or write from the socket is different, although the latter fact is hidden in
 the source by the use of `Curl_read()` for reading and `Curl_write()` for
 writing data to the remote server.

 `http_chunks.c` contains functions that understands HTTP 1.1 chunked transfer
 encoding.

 An interesting detail with the HTTP(S) request, is the `Curl_add_buffer()`
 series of functions we use. They append data to one single buffer, and when
 the building is finished the entire request is sent off in one single write.
 This is done this way to overcome problems with flawed firewalls and lame
 servers.

<a name="ftp"></a>
FTP
===

 The `Curl_if2ip()` function can be used for getting the IP number of a
 specified network interface, and it resides in `lib/if2ip.c`.

 `Curl_ftpsendf()` is used for sending FTP commands to the remote server. It
 was made a separate function to prevent us programmers from forgetting that
 they must be CRLF terminated. They must also be sent in one single `write()`
 to make firewalls and similar happy.

<a name="kerberos"></a>
Kerberos
========

 Kerberos support is mainly in `lib/krb5.c` but also `curl_sasl_sspi.c` and
 `curl_sasl_gssapi.c` for the email protocols and `socks_gssapi.c` and
 `socks_sspi.c` for SOCKS5 proxy specifics.

<a name="telnet"></a>
TELNET
======

 Telnet is implemented in `lib/telnet.c`.

<a name="file"></a>
FILE
====

 The `file://` protocol is dealt with in `lib/file.c`.

<a name="smb"></a>
SMB
===

 The `smb://` protocol is dealt with in `lib/smb.c`.

<a name="ldap"></a>
LDAP
====

 Everything LDAP is in `lib/ldap.c` and `lib/openldap.c`.

<a name="email"></a>
Email
======

 The email related source code is in `lib/imap.c`, `lib/pop3.c` and
 `lib/smtp.c`.

<a name="general"></a>
General
=======

 URL encoding and decoding, called escaping and unescaping in the source code,
 is found in `lib/escape.c`.

 While transferring data in `Transfer()` a few functions might get used.
 `curl_getdate()` in `lib/parsedate.c` is for HTTP date comparisons (and
 more).

 `lib/getenv.c` offers `curl_getenv()` which is for reading environment
 variables in a neat platform independent way. That is used in the client, but
 also in `lib/url.c` when checking the proxy environment variables. Note that
 contrary to the normal unix `getenv()`, this returns an allocated buffer that
 must be `free()`ed after use.

 `lib/netrc.c` holds the `.netrc` parser.

 `lib/timeval.c` features replacement functions for systems that do not have
 `gettimeofday()` and a few support functions for timeval conversions.

 A function named `curl_version()` that returns the full curl version string
 is found in `lib/version.c`.

<a name="persistent"></a>
Persistent Connections
======================

 The persistent connection support in libcurl requires some considerations on
 how to do things inside of the library.

 - The `Curl_easy` struct returned in the [`curl_easy_init()`][2] call
   must never hold connection-oriented data. It is meant to hold the root data
   as well as all the options etc that the library-user may choose.

 - The `Curl_easy` struct holds the "connection cache" (an array of
   pointers to `connectdata` structs).

 - This enables the 'curl handle' to be reused on subsequent transfers.

 - When libcurl is told to perform a transfer, it first checks for an already
   existing connection in the cache that we can use. Otherwise it creates a
   new one and adds that to the cache. If the cache is full already when a new
   connection is added, it will first close the oldest unused one.

 - When the transfer operation is complete, the connection is left
   open. Particular options may tell libcurl not to, and protocols may signal
   closure on connections and then they will not be kept open, of course.

 - When `curl_easy_cleanup()` is called, we close all still opened connections,
   unless of course the multi interface "owns" the connections.

 The curl handle must be re-used in order for the persistent connections to
 work.

<a name="multi"></a>
multi interface/non-blocking
============================

 The multi interface is a non-blocking interface to the library. To make that
 interface work as well as possible, no low-level functions within libcurl
 must be written to work in a blocking manner. (There are still a few spots
 violating this rule.)

 One of the primary reasons we introduced c-ares support was to allow the name
 resolve phase to be perfectly non-blocking as well.

 The FTP and the SFTP/SCP protocols are examples of how we adapt and adjust
 the code to allow non-blocking operations even on multi-stage command-
 response protocols. They are built around state machines that return when
 they would otherwise block waiting for data. The DICT, LDAP and TELNET
 protocols are crappy examples and they are subject for rewrite in the future
 to better fit the libcurl protocol family.

<a name="ssl"></a>
SSL libraries
=============

 Originally libcurl supported SSLeay for SSL/TLS transports, but that was then
 extended to its successor OpenSSL but has since also been extended to several
 other SSL/TLS libraries and we expect and hope to further extend the support
 in future libcurl versions.

 To deal with this internally in the best way possible, we have a generic SSL
 function API as provided by the `vtls/vtls.[ch]` system, and they are the only
 SSL functions we must use from within libcurl. vtls is then crafted to use
 the appropriate lower-level function calls to whatever SSL library that is in
 use. For example `vtls/openssl.[ch]` for the OpenSSL library.

<a name="symbols"></a>
Library Symbols
===============

 All symbols used internally in libcurl must use a `Curl_` prefix if they are
 used in more than a single file. Single-file symbols must be made static.
 Public ("exported") symbols must use a `curl_` prefix. (There are exceptions,
 but they are to be changed to follow this pattern in future versions.) Public
 API functions are marked with `CURL_EXTERN` in the public header files so
 that all others can be hidden on platforms where this is possible.

<a name="returncodes"></a>
Return Codes and Informationals
===============================

 I have made things simple. Almost every function in libcurl returns a CURLcode,
 that must be `CURLE_OK` if everything is OK or otherwise a suitable error
 code as the `curl/curl.h` include file defines. The place that detects an
 error must use the `Curl_failf()` function to set the human-readable error
 description.

 In aiding the user to understand what's happening and to debug curl usage, we
 must supply a fair number of informational messages by using the
 `Curl_infof()` function. Those messages are only displayed when the user
 explicitly asks for them. They are best used when revealing information that
 is not otherwise obvious.

<a name="abi"></a>
API/ABI
=======

 We make an effort to not export or show internals or how internals work, as
 that makes it easier to keep a solid API/ABI over time. See docs/libcurl/ABI
 for our promise to users.

<a name="client"></a>
Client
======

 `main()` resides in `src/tool_main.c`.

 `src/tool_hugehelp.c` is automatically generated by the `mkhelp.pl` perl
 script to display the complete "manual" and the `src/tool_urlglob.c` file
 holds the functions used for the URL-"globbing" support. Globbing in the
 sense that the `{}` and `[]` expansion stuff is there.

 The client mostly sets up its `config` struct properly, then
 it calls the `curl_easy_*()` functions of the library and when it gets back
 control after the `curl_easy_perform()` it cleans up the library, checks
 status and exits.

 When the operation is done, the `ourWriteOut()` function in `src/writeout.c`
 may be called to report about the operation. That function is mostly using the
 `curl_easy_getinfo()` function to extract useful information from the curl
 session.

 It may loop and do all this several times if many URLs were specified on the
 command line or config file.

<a name="memorydebug"></a>
Memory Debugging
================

 The file `lib/memdebug.c` contains debug-versions of a few functions.
 Functions such as `malloc()`, `free()`, `fopen()`, `fclose()`, etc that
 somehow deal with resources that might give us problems if we "leak" them.
 The functions in the memdebug system do nothing fancy, they do their normal
 function and then log information about what they just did. The logged data
 can then be analyzed after a complete session,

 `memanalyze.pl` is the perl script present in `tests/` that analyzes a log
 file generated by the memory tracking system. It detects if resources are
 allocated but never freed and other kinds of errors related to resource
 management.

 Internally, definition of preprocessor symbol `DEBUGBUILD` restricts code
 which is only compiled for debug enabled builds. And symbol `CURLDEBUG` is
 used to differentiate code which is _only_ used for memory
 tracking/debugging.

 Use `-DCURLDEBUG` when compiling to enable memory debugging, this is also
 switched on by running configure with `--enable-curldebug`. Use
 `-DDEBUGBUILD` when compiling to enable a debug build or run configure with
 `--enable-debug`.

 `curl --version` will list 'Debug' feature for debug enabled builds, and
 will list 'TrackMemory' feature for curl debug memory tracking capable
 builds. These features are independent and can be controlled when running
 the configure script. When `--enable-debug` is given both features will be
 enabled, unless some restriction prevents memory tracking from being used.

<a name="test"></a>
Test Suite
==========

 The test suite is placed in its own subdirectory directly off the root in the
 curl archive tree, and it contains a bunch of scripts and a lot of test case
 data.

 The main test script is `runtests.pl` that will invoke test servers like
 `httpserver.pl` and `ftpserver.pl` before all the test cases are performed.
 The test suite currently only runs on Unix-like platforms.

 you will find a description of the test suite in the `tests/README` file, and
 the test case data files in the `tests/FILEFORMAT` file.

 The test suite automatically detects if curl was built with the memory
 debugging enabled, and if it was, it will detect memory leaks, too.

<a name="asyncdns"></a>
Asynchronous name resolves
==========================

 libcurl can be built to do name resolves asynchronously, using either the
 normal resolver in a threaded manner or by using c-ares.

<a name="cares"></a>
[c-ares][3]
------

### Build libcurl to use a c-ares

1. ./configure --enable-ares=/path/to/ares/install
2. make

### c-ares on win32

 First I compiled c-ares. I changed the default C runtime library to be the
 single-threaded rather than the multi-threaded (this seems to be required to
 prevent linking errors later on). Then I simply build the areslib project
 (the other projects adig/ahost seem to fail under MSVC).

 Next was libcurl. I opened `lib/config-win32.h` and I added a:
 `#define USE_ARES 1`

 Next thing I did was I added the path for the ares includes to the include
 path, and the libares.lib to the libraries.

 Lastly, I also changed libcurl to be single-threaded rather than
 multi-threaded, again this was to prevent some duplicate symbol errors. I'm
 not sure why I needed to change everything to single-threaded, but when I
 did not I got redefinition errors for several CRT functions (`malloc()`,
 `stricmp()`, etc.)

<a name="curl_off_t"></a>
`curl_off_t`
==========

 `curl_off_t` is a data type provided by the external libcurl include
 headers. It is the type meant to be used for the [`curl_easy_setopt()`][1]
 options that end with LARGE. The type is 64-bit large on most modern
 platforms.

<a name="curlx"></a>
curlx
=====

 The libcurl source code offers a few functions by source only. They are not
 part of the official libcurl API, but the source files might be useful for
 others so apps can optionally compile/build with these sources to gain
 additional functions.

 We provide them through a single header file for easy access for apps:
 `curlx.h`

`curlx_strtoofft()`
-------------------
   A macro that converts a string containing a number to a `curl_off_t` number.
   This might use the `curlx_strtoll()` function which is provided as source
   code in strtoofft.c. Note that the function is only provided if no
   `strtoll()` (or equivalent) function exist on your platform. If `curl_off_t`
   is only a 32-bit number on your platform, this macro uses `strtol()`.

Future
------

 Several functions will be removed from the public `curl_` name space in a
 future libcurl release. They will then only become available as `curlx_`
 functions instead. To make the transition easier, we already today provide
 these functions with the `curlx_` prefix to allow sources to be built
 properly with the new function names. The concerned functions are:

 - `curlx_getenv`
 - `curlx_strequal`
 - `curlx_strnequal`
 - `curlx_mvsnprintf`
 - `curlx_msnprintf`
 - `curlx_maprintf`
 - `curlx_mvaprintf`
 - `curlx_msprintf`
 - `curlx_mprintf`
 - `curlx_mfprintf`
 - `curlx_mvsprintf`
 - `curlx_mvprintf`
 - `curlx_mvfprintf`

<a name="contentencoding"></a>
Content Encoding
================

## About content encodings

 [HTTP/1.1][4] specifies that a client may request that a server encode its
 response. This is usually used to compress a response using one (or more)
 encodings from a set of commonly available compression techniques. These
 schemes include `deflate` (the zlib algorithm), `gzip`, `br` (brotli) and
 `compress`. A client requests that the server perform an encoding by including
 an `Accept-Encoding` header in the request document. The value of the header
 should be one of the recognized tokens `deflate`, ... (there's a way to
 register new schemes/tokens, see sec 3.5 of the spec). A server MAY honor
 the client's encoding request. When a response is encoded, the server
 includes a `Content-Encoding` header in the response. The value of the
 `Content-Encoding` header indicates which encodings were used to encode the
 data, in the order in which they were applied.

 It's also possible for a client to attach priorities to different schemes so
 that the server knows which it prefers. See sec 14.3 of RFC 2616 for more
 information on the `Accept-Encoding` header. See sec
 [3.1.2.2 of RFC 7231][15] for more information on the `Content-Encoding`
 header.

## Supported content encodings

 The `deflate`, `gzip` and `br` content encodings are supported by libcurl.
 Both regular and chunked transfers work fine. The zlib library is required
 for the `deflate` and `gzip` encodings, while the brotli decoding library is
 for the `br` encoding.

## The libcurl interface

 To cause libcurl to request a content encoding use:

  [`curl_easy_setopt`][1](curl, [`CURLOPT_ACCEPT_ENCODING`][5], string)

 where string is the intended value of the `Accept-Encoding` header.

 Currently, libcurl does support multiple encodings but only
 understands how to process responses that use the `deflate`, `gzip` and/or
 `br` content encodings, so the only values for [`CURLOPT_ACCEPT_ENCODING`][5]
 that will work (besides `identity`, which does nothing) are `deflate`,
 `gzip` and `br`. If a response is encoded using the `compress` or methods,
 libcurl will return an error indicating that the response could
 not be decoded. If `<string>` is NULL no `Accept-Encoding` header is
 generated. If `<string>` is a zero-length string, then an `Accept-Encoding`
 header containing all supported encodings will be generated.

 The [`CURLOPT_ACCEPT_ENCODING`][5] must be set to any non-NULL value for
 content to be automatically decoded. If it is not set and the server still
 sends encoded content (despite not having been asked), the data is returned
 in its raw form and the `Content-Encoding` type is not checked.

## The curl interface

 Use the [`--compressed`][6] option with curl to cause it to ask servers to
 compress responses using any format supported by curl.

<a name="hostip"></a>
`hostip.c` explained
====================

 The main compile-time defines to keep in mind when reading the `host*.c`
 source file are these:

## `CURLRES_IPV6`

 this host has `getaddrinfo()` and family, and thus we use that. The host may
 not be able to resolve IPv6, but we do not really have to take that into
 account. Hosts that are not IPv6-enabled have `CURLRES_IPV4` defined.

## `CURLRES_ARES`

 is defined if libcurl is built to use c-ares for asynchronous name
 resolves. This can be Windows or \*nix.

## `CURLRES_THREADED`

 is defined if libcurl is built to use threading for asynchronous name
 resolves. The name resolve will be done in a new thread, and the supported
 asynch API will be the same as for ares-builds. This is the default under
 (native) Windows.

 If any of the two previous are defined, `CURLRES_ASYNCH` is defined too. If
 libcurl is not built to use an asynchronous resolver, `CURLRES_SYNCH` is
 defined.

## `host*.c` sources

 The `host*.c` sources files are split up like this:

 - `hostip.c`      - method-independent resolver functions and utility functions
 - `hostasyn.c`    - functions for asynchronous name resolves
 - `hostsyn.c`     - functions for synchronous name resolves
 - `asyn-ares.c`   - functions for asynchronous name resolves using c-ares
 - `asyn-thread.c` - functions for asynchronous name resolves using threads
 - `hostip4.c`     - IPv4 specific functions
 - `hostip6.c`     - IPv6 specific functions

 The `hostip.h` is the single united header file for all this. It defines the
 `CURLRES_*` defines based on the `config*.h` and `curl_setup.h` defines.

<a name="memoryleak"></a>
Track Down Memory Leaks
=======================

## Single-threaded

  Please note that this memory leak system is not adjusted to work in more
  than one thread. If you want/need to use it in a multi-threaded app. Please
  adjust accordingly.

## Build

  Rebuild libcurl with `-DCURLDEBUG` (usually, rerunning configure with
  `--enable-debug` fixes this). `make clean` first, then `make` so that all
  files are actually rebuilt properly. It will also make sense to build
  libcurl with the debug option (usually `-g` to the compiler) so that
  debugging it will be easier if you actually do find a leak in the library.

  This will create a library that has memory debugging enabled.

## Modify Your Application

  Add a line in your application code:

```c
  curl_dbg_memdebug("dump");
```

  This will make the malloc debug system output a full trace of all resource
  using functions to the given file name. Make sure you rebuild your program
  and that you link with the same libcurl you built for this purpose as
  described above.

## Run Your Application

  Run your program as usual. Watch the specified memory trace file grow.

  Make your program exit and use the proper libcurl cleanup functions etc. So
  that all non-leaks are returned/freed properly.

## Analyze the Flow

  Use the `tests/memanalyze.pl` perl script to analyze the dump file:

    tests/memanalyze.pl dump

  This now outputs a report on what resources that were allocated but never
  freed etc. This report is fine for posting to the list.

  If this does not produce any output, no leak was detected in libcurl. Then
  the leak is mostly likely to be in your code.

<a name="multi_socket"></a>
`multi_socket`
==============

 Implementation of the `curl_multi_socket` API

 The main ideas of this API are simply:

 1. The application can use whatever event system it likes as it gets info
    from libcurl about what file descriptors libcurl waits for what action
    on. (The previous API returns `fd_sets` which is `select()`-centric).

 2. When the application discovers action on a single socket, it calls
    libcurl and informs that there was action on this particular socket and
    libcurl can then act on that socket/transfer only and not care about
    any other transfers. (The previous API always had to scan through all
    the existing transfers.)

 The idea is that [`curl_multi_socket_action()`][7] calls a given callback
 with information about what socket to wait for what action on, and the
 callback only gets called if the status of that socket has changed.

 We also added a timer callback that makes libcurl call the application when
 the timeout value changes, and you set that with [`curl_multi_setopt()`][9]
 and the [`CURLMOPT_TIMERFUNCTION`][10] option. To get this to work,
 Internally, there's an added struct to each easy handle in which we store
 an "expire time" (if any). The structs are then "splay sorted" so that we
 can add and remove times from the linked list and yet somewhat swiftly
 figure out both how long there is until the next nearest timer expires
 and which timer (handle) we should take care of now. Of course, the upside
 of all this is that we get a [`curl_multi_timeout()`][8] that should also
 work with old-style applications that use [`curl_multi_perform()`][11].

 We created an internal "socket to easy handles" hash table that given
 a socket (file descriptor) returns the easy handle that waits for action on
 that socket. This hash is made using the already existing hash code
 (previously only used for the DNS cache).

 To make libcurl able to report plain sockets in the socket callback, we had
 to re-organize the internals of the [`curl_multi_fdset()`][12] etc so that
 the conversion from sockets to `fd_sets` for that function is only done in
 the last step before the data is returned. I also had to extend c-ares to
 get a function that can return plain sockets, as that library too returned
 only `fd_sets` and that is no longer good enough. The changes done to c-ares
 are available in c-ares 1.3.1 and later.

<a name="structs"></a>
Structs in libcurl
==================

This section should cover 7.32.0 pretty accurately, but will make sense even
for older and later versions as things do not change drastically that often.

<a name="Curl_easy"></a>
## Curl_easy

  The `Curl_easy` struct is the one returned to the outside in the external API
  as a `CURL *`. This is usually known as an easy handle in API documentations
  and examples.

  Information and state that is related to the actual connection is in the
  `connectdata` struct. When a transfer is about to be made, libcurl will
  either create a new connection or re-use an existing one. The particular
  connectdata that is used by this handle is pointed out by
  `Curl_easy->easy_conn`.

  Data and information that regard this particular single transfer is put in
  the `SingleRequest` sub-struct.

  When the `Curl_easy` struct is added to a multi handle, as it must be in
  order to do any transfer, the `->multi` member will point to the `Curl_multi`
  struct it belongs to. The `->prev` and `->next` members will then be used by
  the multi code to keep a linked list of `Curl_easy` structs that are added to
  that same multi handle. libcurl always uses multi so `->multi` *will* point
  to a `Curl_multi` when a transfer is in progress.

  `->mstate` is the multi state of this particular `Curl_easy`. When
  `multi_runsingle()` is called, it will act on this handle according to which
  state it is in. The mstate is also what tells which sockets to return for a
  specific `Curl_easy` when [`curl_multi_fdset()`][12] is called etc.

  The libcurl source code generally use the name `data` for the variable that
  points to the `Curl_easy`.

  When doing multiplexed HTTP/2 transfers, each `Curl_easy` is associated with
  an individual stream, sharing the same connectdata struct. Multiplexing
  makes it even more important to keep things associated with the right thing!

<a name="connectdata"></a>
## connectdata

  A general idea in libcurl is to keep connections around in a connection
  "cache" after they have been used in case they will be used again and then
  re-use an existing one instead of creating a new as it creates a significant
  performance boost.

  Each `connectdata` identifies a single physical connection to a server. If
  the connection cannot be kept alive, the connection will be closed after use
  and then this struct can be removed from the cache and freed.

  Thus, the same `Curl_easy` can be used multiple times and each time select
  another `connectdata` struct to use for the connection. Keep this in mind,
  as it is then important to consider if options or choices are based on the
  connection or the `Curl_easy`.

  Functions in libcurl will assume that `connectdata->data` points to the
  `Curl_easy` that uses this connection (for the moment).

  As a special complexity, some protocols supported by libcurl require a
  special disconnect procedure that is more than just shutting down the
  socket. It can involve sending one or more commands to the server before
  doing so. Since connections are kept in the connection cache after use, the
  original `Curl_easy` may no longer be around when the time comes to shut down
  a particular connection. For this purpose, libcurl holds a special dummy
  `closure_handle` `Curl_easy` in the `Curl_multi` struct to use when needed.

  FTP uses two TCP connections for a typical transfer but it keeps both in
  this single struct and thus can be considered a single connection for most
  internal concerns.

  The libcurl source code generally use the name `conn` for the variable that
  points to the connectdata.

<a name="Curl_multi"></a>
## Curl_multi

  Internally, the easy interface is implemented as a wrapper around multi
  interface functions. This makes everything multi interface.

  `Curl_multi` is the multi handle struct exposed as `CURLM *` in external
  APIs.

  This struct holds a list of `Curl_easy` structs that have been added to this
  handle with [`curl_multi_add_handle()`][13]. The start of the list is
  `->easyp` and `->num_easy` is a counter of added `Curl_easy`s.

  `->msglist` is a linked list of messages to send back when
  [`curl_multi_info_read()`][14] is called. Basically a node is added to that
  list when an individual `Curl_easy`'s transfer has completed.

  `->hostcache` points to the name cache. It is a hash table for looking up
  name to IP. The nodes have a limited life time in there and this cache is
  meant to reduce the time for when the same name is wanted within a short
  period of time.

  `->timetree` points to a tree of `Curl_easy`s, sorted by the remaining time
  until it should be checked - normally some sort of timeout. Each `Curl_easy`
  has one node in the tree.

  `->sockhash` is a hash table to allow fast lookups of socket descriptor for
  which `Curl_easy` uses that descriptor. This is necessary for the
  `multi_socket` API.

  `->conn_cache` points to the connection cache. It keeps track of all
  connections that are kept after use. The cache has a maximum size.

  `->closure_handle` is described in the `connectdata` section.

  The libcurl source code generally use the name `multi` for the variable that
  points to the `Curl_multi` struct.

<a name="Curl_handler"></a>
## Curl_handler

  Each unique protocol that is supported by libcurl needs to provide at least
  one `Curl_handler` struct. It defines what the protocol is called and what
  functions the main code should call to deal with protocol specific issues.
  In general, there's a source file named `[protocol].c` in which there's a
  `struct Curl_handler Curl_handler_[protocol]` declared. In `url.c` there's
  then the main array with all individual `Curl_handler` structs pointed to
  from a single array which is scanned through when a URL is given to libcurl
  to work with.

  The concrete function pointer prototypes can be found in `lib/urldata.h`.

  `->scheme` is the URL scheme name, usually spelled out in uppercase. That is
  "HTTP" or "FTP" etc. SSL versions of the protocol need their own
  `Curl_handler` setup so HTTPS separate from HTTP.

  `->setup_connection` is called to allow the protocol code to allocate
  protocol specific data that then gets associated with that `Curl_easy` for
  the rest of this transfer. It gets freed again at the end of the transfer.
  It will be called before the `connectdata` for the transfer has been
  selected/created. Most protocols will allocate its private `struct
  [PROTOCOL]` here and assign `Curl_easy->req.p.[protocol]` to it.

  `->connect_it` allows a protocol to do some specific actions after the TCP
  connect is done, that can still be considered part of the connection phase.

  Some protocols will alter the `connectdata->recv[]` and
  `connectdata->send[]` function pointers in this function.

  `->connecting` is similarly a function that keeps getting called as long as
  the protocol considers itself still in the connecting phase.

  `->do_it` is the function called to issue the transfer request. What we call
  the DO action internally. If the DO is not enough and things need to be kept
  getting done for the entire DO sequence to complete, `->doing` is then
  usually also provided. Each protocol that needs to do multiple commands or
  similar for do/doing need to implement their own state machines (see SCP,
  SFTP, FTP). Some protocols (only FTP and only due to historical reasons) has
  a separate piece of the DO state called `DO_MORE`.

  `->doing` keeps getting called while issuing the transfer request command(s)

  `->done` gets called when the transfer is complete and DONE. That is after the
  main data has been transferred.

  `->do_more` gets called during the `DO_MORE` state. The FTP protocol uses
  this state when setting up the second connection.

  `->proto_getsock`
  `->doing_getsock`
  `->domore_getsock`
  `->perform_getsock`
  Functions that return socket information. Which socket(s) to wait for which
  I/O action(s) during the particular multi state.

  `->disconnect` is called immediately before the TCP connection is shutdown.

  `->readwrite` gets called during transfer to allow the protocol to do extra
  reads/writes

  `->attach` attaches a transfer to the connection.

  `->defport` is the default report TCP or UDP port this protocol uses

  `->protocol` is one or more bits in the `CURLPROTO_*` set. The SSL versions
  have their "base" protocol set and then the SSL variation. Like
  "HTTP|HTTPS".

  `->flags` is a bitmask with additional information about the protocol that will
  make it get treated differently by the generic engine:

  - `PROTOPT_SSL` - will make it connect and negotiate SSL

  - `PROTOPT_DUAL` - this protocol uses two connections

  - `PROTOPT_CLOSEACTION` - this protocol has actions to do before closing the
    connection. This flag is no longer used by code, yet still set for a bunch
    of protocol handlers.

  - `PROTOPT_DIRLOCK` - "direction lock". The SSH protocols set this bit to
    limit which "direction" of socket actions that the main engine will
    concern itself with.

  - `PROTOPT_NONETWORK` - a protocol that does not use network (read `file:`)

  - `PROTOPT_NEEDSPWD` - this protocol needs a password and will use a default
    one unless one is provided

  - `PROTOPT_NOURLQUERY` - this protocol cannot handle a query part on the URL
    (?foo=bar)

<a name="conncache"></a>
## conncache

  Is a hash table with connections for later re-use. Each `Curl_easy` has a
  pointer to its connection cache. Each multi handle sets up a connection
  cache that all added `Curl_easy`s share by default.

<a name="Curl_share"></a>
## Curl_share

  The libcurl share API allocates a `Curl_share` struct, exposed to the
  external API as `CURLSH *`.

  The idea is that the struct can have a set of its own versions of caches and
  pools and then by providing this struct in the `CURLOPT_SHARE` option, those
  specific `Curl_easy`s will use the caches/pools that this share handle
  holds.

  Then individual `Curl_easy` structs can be made to share specific things
  that they otherwise would not, such as cookies.

  The `Curl_share` struct can currently hold cookies, DNS cache and the SSL
  session cache.

<a name="CookieInfo"></a>
## CookieInfo

  This is the main cookie struct. It holds all known cookies and related
  information. Each `Curl_easy` has its own private `CookieInfo` even when
  they are added to a multi handle. They can be made to share cookies by using
  the share API.


[1]: https://curl.se/libcurl/c/curl_easy_setopt.html
[2]: https://curl.se/libcurl/c/curl_easy_init.html
[3]: https://c-ares.org/
[4]: https://datatracker.ietf.org/doc/html/rfc7230 "RFC 7230"
[5]: https://curl.se/libcurl/c/CURLOPT_ACCEPT_ENCODING.html
[6]: https://curl.se/docs/manpage.html#--compressed
[7]: https://curl.se/libcurl/c/curl_multi_socket_action.html
[8]: https://curl.se/libcurl/c/curl_multi_timeout.html
[9]: https://curl.se/libcurl/c/curl_multi_setopt.html
[10]: https://curl.se/libcurl/c/CURLMOPT_TIMERFUNCTION.html
[11]: https://curl.se/libcurl/c/curl_multi_perform.html
[12]: https://curl.se/libcurl/c/curl_multi_fdset.html
[13]: https://curl.se/libcurl/c/curl_multi_add_handle.html
[14]: https://curl.se/libcurl/c/curl_multi_info_read.html
[15]: https://datatracker.ietf.org/doc/html/rfc7231#section-3.1.2.2
![curl logo](https://curl.se/logo/curl-logo.svg)

# Documentation

you will find a mix of various documentation in this directory and
subdirectories, using several different formats. Some of them are not ideal
for reading directly in your browser.

If you would rather see the rendered version of the documentation, check out the
curl website's [documentation section](https://curl.se/docs/) for
general curl stuff or the [libcurl section](https://curl.se/libcurl/) for
libcurl related documentation.
# HTTP Cookies

## Cookie overview

  Cookies are `name=contents` pairs that a HTTP server tells the client to
  hold and then the client sends back those to the server on subsequent
  requests to the same domains and paths for which the cookies were set.

  Cookies are either "session cookies" which typically are forgotten when the
  session is over which is often translated to equal when browser quits, or
  the cookies are not session cookies they have expiration dates after which
  the client will throw them away.

  Cookies are set to the client with the Set-Cookie: header and are sent to
  servers with the Cookie: header.

  For a long time, the only spec explaining how to use cookies was the
  original [Netscape spec from 1994](https://curl.se/rfc/cookie_spec.html).

  In 2011, [RFC6265](https://www.ietf.org/rfc/rfc6265.txt) was finally
  published and details how cookies work within HTTP. In 2016, an update which
  added support for prefixes was
  [proposed](https://datatracker.ietf.org/doc/html/draft-ietf-httpbis-cookie-prefixes-00),
  and in 2017, another update was
  [drafted](https://datatracker.ietf.org/doc/html/draft-ietf-httpbis-cookie-alone-01)
  to deprecate modification of 'secure' cookies from non-secure origins. Both
  of these drafts have been incorporated into a proposal to
  [replace](https://datatracker.ietf.org/doc/html/draft-ietf-httpbis-rfc6265bis-02)
  RFC6265. Cookie prefixes and secure cookie modification protection has been
  implemented by curl.

## Cookies saved to disk

  Netscape once created a file format for storing cookies on disk so that they
  would survive browser restarts. curl adopted that file format to allow
  sharing the cookies with browsers, only to see browsers move away from that
  format. Modern browsers no longer use it, while curl still does.

  The netscape cookie file format stores one cookie per physical line in the
  file with a bunch of associated meta data, each field separated with
  TAB. That file is called the cookiejar in curl terminology.

  When libcurl saves a cookiejar, it creates a file header of its own in which
  there is a URL mention that will link to the web version of this document.

## Cookie file format

  The cookie file format is text based and stores one cookie per line. Lines
  that start with `#` are treated as comments.

  Each line that specifies a single cookie consists of seven text fields
  separated with TAB characters. A valid line must end with a newline
  character.

### Fields in the file

  Field number, what type and example data and the meaning of it:

  0. string `example.com` - the domain name
  1. boolean `FALSE` - include subdomains
  2. string `/foobar/` - path
  3. boolean `TRUE` - send/receive over HTTPS only
  4. number `1462299217` - expires at - seconds since Jan 1st 1970, or 0
  5. string `person` - name of the cookie
  6. string `daniel` - value of the cookie

## Cookies with curl the command line tool

  curl has a full cookie "engine" built in. If you just activate it, you can
  have curl receive and send cookies exactly as mandated in the specs.

  Command line options:

  `-b, --cookie`

  tell curl a file to read cookies from and start the cookie engine, or if it
  is not a file it will pass on the given string. -b name=var works and so does
  -b cookiefile.

  `-j, --junk-session-cookies`

  when used in combination with -b, it will skip all "session cookies" on load
  so as to appear to start a new cookie session.

  `-c, --cookie-jar`

  tell curl to start the cookie engine and write cookies to the given file
  after the request(s)

## Cookies with libcurl

  libcurl offers several ways to enable and interface the cookie engine. These
  options are the ones provided by the native API. libcurl bindings may offer
  access to them using other means.

  `CURLOPT_COOKIE`

  Is used when you want to specify the exact contents of a cookie header to
  send to the server.

  `CURLOPT_COOKIEFILE`

  Tell libcurl to activate the cookie engine, and to read the initial set of
  cookies from the given file. Read-only.

  `CURLOPT_COOKIEJAR`

  Tell libcurl to activate the cookie engine, and when the easy handle is
  closed save all known cookies to the given cookiejar file. Write-only.

  `CURLOPT_COOKIELIST`

  Provide detailed information about a single cookie to add to the internal
  storage of cookies. Pass in the cookie as a HTTP header with all the details
  set, or pass in a line from a netscape cookie file. This option can also be
  used to flush the cookies etc.

  `CURLINFO_COOKIELIST`

  Extract cookie information from the internal cookie storage as a linked
  list.

## Cookies with javascript

  These days a lot of the web is built up by javascript. The webbrowser loads
  complete programs that render the page you see. These javascript programs
  can also set and access cookies.

  Since curl and libcurl are plain HTTP clients without any knowledge of or
  capability to handle javascript, such cookies will not be detected or used.

  Often, if you want to mimic what a browser does on such websites, you can
  record web browser HTTP traffic when using such a site and then repeat the
  cookie operations using curl or libcurl.
libcurl bindings
================

 Creative people have written bindings or interfaces for various environments
 and programming languages. Using one of these allows you to take advantage of
 curl powers from within your favourite language or system.

 This is a list of all known interfaces as of this writing.

 The bindings listed below are not part of the curl/libcurl distribution
 archives, but must be downloaded and installed separately.

<!-- markdown-link-check-disable -->

[Ada95](https://web.archive.org/web/20070403105909/www.almroth.com/adacurl/index.html) Written by Andreas Almroth

[Basic](https://scriptbasic.com/) ScriptBasic bindings written by Peter Verhas

C++: [curlpp](https://github.com/jpbarrette/curlpp/) Written by Jean-Philippe Barrette-LaPierre,
[curlcpp](https://github.com/JosephP91/curlcpp) by Giuseppe Persico and [C++
Requests](https://github.com/libcpr/cpr) by Huu Nguyen

[Ch](https://chcurl.sourceforge.io/) Written by Stephen Nestinger and Jonathan Rogado

Cocoa: [BBHTTP](https://github.com/biasedbit/BBHTTP) written by Bruno de Carvalho
[curlhandle](https://github.com/karelia/curlhandle) Written by Dan Wood

Clojure: [clj-curl](https://github.com/lsevero/clj-curl) by Lucas Severo

[D](https://dlang.org/library/std/net/curl.html) Written by Kenneth Bogert

[Delphi](https://github.com/Mercury13/curl4delphi) Written by Mikhail Merkuryev

[Dylan](https://dylanlibs.sourceforge.io/) Written by Chris Double

[Eiffel](https://iron.eiffel.com/repository/20.11/package/ABEF6975-37AC-45FD-9C67-52D10BA0669B) Written by Eiffel Software

[Euphoria](https://web.archive.org/web/20050204080544/rays-web.com/eulibcurl.htm) Written by Ray Smith

[Falcon](http://www.falconpl.org/project_docs/curl/)

[Ferite](https://web.archive.org/web/20150102192018/ferite.org/) Written by Paul Querna

[Gambas](https://gambas.sourceforge.io/)

[glib/GTK+](https://web.archive.org/web/20100526203452/atterer.net/glibcurl) Written by Richard Atterer

Go: [go-curl](https://github.com/andelf/go-curl) by ShuYu Wang

[Guile](https://github.com/spk121/guile-curl) Written by Michael L. Gran

[Harbour](https://github.com/vszakats/hb/tree/main/contrib/hbcurl) Written by Viktor Szakats

[Haskell](https://hackage.haskell.org/package/curl) Written by Galois, Inc

[Java](https://github.com/pjlegato/curl-java)

[Julia](https://github.com/JuliaWeb/LibCURL.jl) Written by Amit Murthy

[Kapito](https://github.com/puzza007/katipo) is an Erlang HTTP library around libcurl.

[Lisp](https://common-lisp.net/project/cl-curl/) Written by Liam Healy

Lua: [luacurl](https://web.archive.org/web/20201205052437/luacurl.luaforge.net/) by Alexander Marinov, [Lua-cURL](https://github.com/Lua-cURL) by Jürgen Hötzel

[Mono](https://web.archive.org/web/20070606064500/https://forge.novell.com/modules/xfmod/project/?libcurl-mono) Written by Jeffrey Phillips

[.NET](https://sourceforge.net/projects/libcurl-net/) libcurl-net by Jeffrey Phillips

[Nim](https://nimble.directory/pkg/libcurl) wrapper for libcurl

[node.js](https://github.com/JCMais/node-libcurl) node-libcurl by Jonathan Cardoso Machado

[Object-Pascal](https://web.archive.org/web/20020610214926/www.tekool.com/opcurl) Free Pascal, Delphi and Kylix binding written by Christophe Espern.

[OCaml](https://opam.ocaml.org/packages/ocurl/) Written by Lars Nilsson and ygrek

[Pascal](https://web.archive.org/web/20030804091414/houston.quik.com/jkp/curlpas/) Free Pascal, Delphi and Kylix binding written by Jeffrey Pohlmeyer.

Perl: [WWW::Curl](https://github.com/szbalint/WWW--Curl) Maintained by Cris
Bailiff and Bálint Szilakszi,
[perl6-net-curl](https://github.com/azawawi/perl6-net-curl) by Ahmad M. Zawawi
[NET::Curl](https://metacpan.org/pod/Net::Curl) by Przemyslaw Iskra

[PHP](https://php.net/curl) Originally written by Sterling Hughes

[PostgreSQL](https://github.com/pramsey/pgsql-http) - HTTP client for PostgreSQL

[PostgreSQL](https://github.com/RekGRpth/pg_curl) - cURL client for PostgreSQL

[PureBasic](https://www.purebasic.com/documentation/http/index.html) uses libcurl in its "native" HTTP subsystem

[Python](http://pycurl.io/) PycURL by Kjetil Jacobsen

[Q](https://q-lang.sourceforge.io/) The libcurl module is part of the default install

[R](https://cran.r-project.org/package=curl)

[Rexx](https://rexxcurl.sourceforge.io/) Written Mark Hessling

[Ring](https://ring-lang.sourceforge.io/doc1.3/libcurl.html) RingLibCurl by Mahmoud Fayed

RPG, support for ILE/RPG on OS/400 is included in source distribution

Ruby: [curb](https://github.com/taf2/curb) written by Ross Bamford,
[ruby-curl-multi](https://github.com/kball/curl_multi.rb) by Kristjan Petursson and Keith Rarick

[Rust](https://github.com/alexcrichton/curl-rust) curl-rust - by Carl Lerche

[Scheme](http://www.metapaper.net/lisovsky/web/curl/) Bigloo binding by Kirill Lisovsky

[Scilab](https://help.scilab.org/docs/current/fr_FR/getURL.html) binding by Sylvestre Ledru

[S-Lang](https://www.jedsoft.org/slang/modules/curl.html) by John E Davis

[Smalltalk](https://www.squeaksource.com/CurlPlugin/) Written by Danil Osipchuk

[SP-Forth](https://sourceforge.net/p/spf/spf/ci/master/tree/devel/~ac/lib/lin/curl/) Written by Andrey Cherezov

[SPL](https://web.archive.org/web/20210203022158/http://www.clifford.at/spl/spldoc/curl.html) Written by Clifford Wolf

[Tcl](https://web.archive.org/web/20160826011806/mirror.yellow5.com/tclcurl/) Tclcurl by Andrés García

[Visual Basic](https://sourceforge.net/projects/libcurl-vb/) libcurl-vb by Jeffrey Phillips

[Visual Foxpro](https://web.archive.org/web/20130730181523/www.ctl32.com.ar/libcurl.asp) by Carlos Alloatti

[wxWidgets](https://wxcode.sourceforge.io/components/wxcurl/) Written by Casey O'Donnell

[XBLite](https://web.archive.org/web/20060426150418/perso.wanadoo.fr/xblite/libraries.html) Written by David Szafranski

[Xojo](https://github.com/charonn0/RB-libcURL) Written by Andrew Lambert
# Ciphers

With curl's options
[`CURLOPT_SSL_CIPHER_LIST`](https://curl.se/libcurl/c/CURLOPT_SSL_CIPHER_LIST.html)
and
[`--ciphers`](https://curl.se/docs/manpage.html#--ciphers)
users can control which ciphers to consider when negotiating TLS connections.

TLS 1.3 ciphers are supported since curl 7.61 for OpenSSL 1.1.1+ with options
[`CURLOPT_TLS13_CIPHERS`](https://curl.se/libcurl/c/CURLOPT_TLS13_CIPHERS.html)
and
[`--tls13-ciphers`](https://curl.se/docs/manpage.html#--tls13-ciphers)
. If you are using a different SSL backend you can try setting TLS 1.3 cipher
suites by using the respective regular cipher option.

The names of the known ciphers differ depending on which TLS backend that
libcurl was built to use. This is an attempt to list known cipher names.

## OpenSSL

(based on [OpenSSL docs](https://www.openssl.org/docs/manmaster/man1/openssl-ciphers.html))

When specifying multiple cipher names, separate them with colon (`:`).

### SSL3 cipher suites

`NULL-MD5`
`NULL-SHA`
`RC4-MD5`
`RC4-SHA`
`IDEA-CBC-SHA`
`DES-CBC3-SHA`
`DH-DSS-DES-CBC3-SHA`
`DH-RSA-DES-CBC3-SHA`
`DHE-DSS-DES-CBC3-SHA`
`DHE-RSA-DES-CBC3-SHA`
`ADH-RC4-MD5`
`ADH-DES-CBC3-SHA`

### TLS v1.0 cipher suites

`NULL-MD5`
`NULL-SHA`
`RC4-MD5`
`RC4-SHA`
`IDEA-CBC-SHA`
`DES-CBC3-SHA`
`DHE-DSS-DES-CBC3-SHA`
`DHE-RSA-DES-CBC3-SHA`
`ADH-RC4-MD5`
`ADH-DES-CBC3-SHA`

### AES ciphersuites from RFC3268, extending TLS v1.0

`AES128-SHA`
`AES256-SHA`
`DH-DSS-AES128-SHA`
`DH-DSS-AES256-SHA`
`DH-RSA-AES128-SHA`
`DH-RSA-AES256-SHA`
`DHE-DSS-AES128-SHA`
`DHE-DSS-AES256-SHA`
`DHE-RSA-AES128-SHA`
`DHE-RSA-AES256-SHA`
`ADH-AES128-SHA`
`ADH-AES256-SHA`

### SEED ciphersuites from RFC4162, extending TLS v1.0

`SEED-SHA`
`DH-DSS-SEED-SHA`
`DH-RSA-SEED-SHA`
`DHE-DSS-SEED-SHA`
`DHE-RSA-SEED-SHA`
`ADH-SEED-SHA`

### GOST ciphersuites, extending TLS v1.0

`GOST94-GOST89-GOST89`
`GOST2001-GOST89-GOST89`
`GOST94-NULL-GOST94`
`GOST2001-NULL-GOST94`

### Elliptic curve cipher suites

`ECDHE-RSA-NULL-SHA`
`ECDHE-RSA-RC4-SHA`
`ECDHE-RSA-DES-CBC3-SHA`
`ECDHE-RSA-AES128-SHA`
`ECDHE-RSA-AES256-SHA`
`ECDHE-ECDSA-NULL-SHA`
`ECDHE-ECDSA-RC4-SHA`
`ECDHE-ECDSA-DES-CBC3-SHA`
`ECDHE-ECDSA-AES128-SHA`
`ECDHE-ECDSA-AES256-SHA`
`AECDH-NULL-SHA`
`AECDH-RC4-SHA`
`AECDH-DES-CBC3-SHA`
`AECDH-AES128-SHA`
`AECDH-AES256-SHA`

### TLS v1.2 cipher suites

`NULL-SHA256`
`AES128-SHA256`
`AES256-SHA256`
`AES128-GCM-SHA256`
`AES256-GCM-SHA384`
`DH-RSA-AES128-SHA256`
`DH-RSA-AES256-SHA256`
`DH-RSA-AES128-GCM-SHA256`
`DH-RSA-AES256-GCM-SHA384`
`DH-DSS-AES128-SHA256`
`DH-DSS-AES256-SHA256`
`DH-DSS-AES128-GCM-SHA256`
`DH-DSS-AES256-GCM-SHA384`
`DHE-RSA-AES128-SHA256`
`DHE-RSA-AES256-SHA256`
`DHE-RSA-AES128-GCM-SHA256`
`DHE-RSA-AES256-GCM-SHA384`
`DHE-DSS-AES128-SHA256`
`DHE-DSS-AES256-SHA256`
`DHE-DSS-AES128-GCM-SHA256`
`DHE-DSS-AES256-GCM-SHA384`
`ECDHE-RSA-AES128-SHA256`
`ECDHE-RSA-AES256-SHA384`
`ECDHE-RSA-AES128-GCM-SHA256`
`ECDHE-RSA-AES256-GCM-SHA384`
`ECDHE-ECDSA-AES128-SHA256`
`ECDHE-ECDSA-AES256-SHA384`
`ECDHE-ECDSA-AES128-GCM-SHA256`
`ECDHE-ECDSA-AES256-GCM-SHA384`
`ADH-AES128-SHA256`
`ADH-AES256-SHA256`
`ADH-AES128-GCM-SHA256`
`ADH-AES256-GCM-SHA384`
`AES128-CCM`
`AES256-CCM`
`DHE-RSA-AES128-CCM`
`DHE-RSA-AES256-CCM`
`AES128-CCM8`
`AES256-CCM8`
`DHE-RSA-AES128-CCM8`
`DHE-RSA-AES256-CCM8`
`ECDHE-ECDSA-AES128-CCM`
`ECDHE-ECDSA-AES256-CCM`
`ECDHE-ECDSA-AES128-CCM8`
`ECDHE-ECDSA-AES256-CCM8`

### Camellia HMAC-Based ciphersuites from RFC6367, extending TLS v1.2

`ECDHE-ECDSA-CAMELLIA128-SHA256`
`ECDHE-ECDSA-CAMELLIA256-SHA384`
`ECDHE-RSA-CAMELLIA128-SHA256`
`ECDHE-RSA-CAMELLIA256-SHA384`

### TLS 1.3 cipher suites

(Note these ciphers are set with `CURLOPT_TLS13_CIPHERS` and `--tls13-ciphers`)

`TLS_AES_256_GCM_SHA384`
`TLS_CHACHA20_POLY1305_SHA256`
`TLS_AES_128_GCM_SHA256`
`TLS_AES_128_CCM_8_SHA256`
`TLS_AES_128_CCM_SHA256`

## NSS

### Totally insecure

`rc4`
`rc4-md5`
`rc4export`
`rc2`
`rc2export`
`des`
`desede3`

###  SSL3/TLS cipher suites

`rsa_rc4_128_md5`
`rsa_rc4_128_sha`
`rsa_3des_sha`
`rsa_des_sha`
`rsa_rc4_40_md5`
`rsa_rc2_40_md5`
`rsa_null_md5`
`rsa_null_sha`
`fips_3des_sha`
`fips_des_sha`
`fortezza`
`fortezza_rc4_128_sha`
`fortezza_null`

### TLS 1.0 Exportable 56-bit Cipher Suites

`rsa_des_56_sha`
`rsa_rc4_56_sha`

### AES ciphers

`dhe_dss_aes_128_cbc_sha`
`dhe_dss_aes_256_cbc_sha`
`dhe_rsa_aes_128_cbc_sha`
`dhe_rsa_aes_256_cbc_sha`
`rsa_aes_128_sha`
`rsa_aes_256_sha`

### ECC ciphers

`ecdh_ecdsa_null_sha`
`ecdh_ecdsa_rc4_128_sha`
`ecdh_ecdsa_3des_sha`
`ecdh_ecdsa_aes_128_sha`
`ecdh_ecdsa_aes_256_sha`
`ecdhe_ecdsa_null_sha`
`ecdhe_ecdsa_rc4_128_sha`
`ecdhe_ecdsa_3des_sha`
`ecdhe_ecdsa_aes_128_sha`
`ecdhe_ecdsa_aes_256_sha`
`ecdh_rsa_null_sha`
`ecdh_rsa_128_sha`
`ecdh_rsa_3des_sha`
`ecdh_rsa_aes_128_sha`
`ecdh_rsa_aes_256_sha`
`ecdhe_rsa_null`
`ecdhe_rsa_rc4_128_sha`
`ecdhe_rsa_3des_sha`
`ecdhe_rsa_aes_128_sha`
`ecdhe_rsa_aes_256_sha`
`ecdh_anon_null_sha`
`ecdh_anon_rc4_128sha`
`ecdh_anon_3des_sha`
`ecdh_anon_aes_128_sha`
`ecdh_anon_aes_256_sha`

### HMAC-SHA256 cipher suites

`rsa_null_sha_256`
`rsa_aes_128_cbc_sha_256`
`rsa_aes_256_cbc_sha_256`
`dhe_rsa_aes_128_cbc_sha_256`
`dhe_rsa_aes_256_cbc_sha_256`
`ecdhe_ecdsa_aes_128_cbc_sha_256`
`ecdhe_rsa_aes_128_cbc_sha_256`

### AES GCM cipher suites in RFC 5288 and RFC 5289

`rsa_aes_128_gcm_sha_256`
`dhe_rsa_aes_128_gcm_sha_256`
`dhe_dss_aes_128_gcm_sha_256`
`ecdhe_ecdsa_aes_128_gcm_sha_256`
`ecdh_ecdsa_aes_128_gcm_sha_256`
`ecdhe_rsa_aes_128_gcm_sha_256`
`ecdh_rsa_aes_128_gcm_sha_256`

### cipher suites using SHA384

`rsa_aes_256_gcm_sha_384`
`dhe_rsa_aes_256_gcm_sha_384`
`dhe_dss_aes_256_gcm_sha_384`
`ecdhe_ecdsa_aes_256_sha_384`
`ecdhe_rsa_aes_256_sha_384`
`ecdhe_ecdsa_aes_256_gcm_sha_384`
`ecdhe_rsa_aes_256_gcm_sha_384`

### chacha20-poly1305 cipher suites

`ecdhe_rsa_chacha20_poly1305_sha_256`
`ecdhe_ecdsa_chacha20_poly1305_sha_256`
`dhe_rsa_chacha20_poly1305_sha_256`

### TLS 1.3 cipher suites

`aes_128_gcm_sha_256`
`aes_256_gcm_sha_384`
`chacha20_poly1305_sha_256`

## GSKit

Ciphers are internally defined as [numeric
codes](https://www.ibm.com/support/knowledgecenter/ssw_ibm_i_73/apis/gsk_attribute_set_buffer.htm). libcurl
maps them to the following case-insensitive names.

### SSL2 cipher suites (insecure: disabled by default)

`rc2-md5`
`rc4-md5`
`exp-rc2-md5`
`exp-rc4-md5`
`des-cbc-md5`
`des-cbc3-md5`

### SSL3 cipher suites

`null-md5`
`null-sha`
`rc4-md5`
`rc4-sha`
`exp-rc2-cbc-md5`
`exp-rc4-md5`
`exp-des-cbc-sha`
`des-cbc3-sha`

### TLS v1.0 cipher suites

`null-md5`
`null-sha`
`rc4-md5`
`rc4-sha`
`exp-rc2-cbc-md5`
`exp-rc4-md5`
`exp-des-cbc-sha`
`des-cbc3-sha`
`aes128-sha`
`aes256-sha`

### TLS v1.1 cipher suites

`null-md5`
`null-sha`
`rc4-md5`
`rc4-sha`
`exp-des-cbc-sha`
`des-cbc3-sha`
`aes128-sha`
`aes256-sha`

### TLS v1.2 cipher suites

`null-md5`
`null-sha`
`null-sha256`
`rc4-md5`
`rc4-sha`
`des-cbc3-sha`
`aes128-sha`
`aes256-sha`
`aes128-sha256`
`aes256-sha256`
`aes128-gcm-sha256`
`aes256-gcm-sha384`

## WolfSSL

`RC4-SHA`,
`RC4-MD5`,
`DES-CBC3-SHA`,
`AES128-SHA`,
`AES256-SHA`,
`NULL-SHA`,
`NULL-SHA256`,
`DHE-RSA-AES128-SHA`,
`DHE-RSA-AES256-SHA`,
`DHE-PSK-AES256-GCM-SHA384`,
`DHE-PSK-AES128-GCM-SHA256`,
`PSK-AES256-GCM-SHA384`,
`PSK-AES128-GCM-SHA256`,
`DHE-PSK-AES256-CBC-SHA384`,
`DHE-PSK-AES128-CBC-SHA256`,
`PSK-AES256-CBC-SHA384`,
`PSK-AES128-CBC-SHA256`,
`PSK-AES128-CBC-SHA`,
`PSK-AES256-CBC-SHA`,
`DHE-PSK-AES128-CCM`,
`DHE-PSK-AES256-CCM`,
`PSK-AES128-CCM`,
`PSK-AES256-CCM`,
`PSK-AES128-CCM-8`,
`PSK-AES256-CCM-8`,
`DHE-PSK-NULL-SHA384`,
`DHE-PSK-NULL-SHA256`,
`PSK-NULL-SHA384`,
`PSK-NULL-SHA256`,
`PSK-NULL-SHA`,
`HC128-MD5`,
`HC128-SHA`,
`HC128-B2B256`,
`AES128-B2B256`,
`AES256-B2B256`,
`RABBIT-SHA`,
`NTRU-RC4-SHA`,
`NTRU-DES-CBC3-SHA`,
`NTRU-AES128-SHA`,
`NTRU-AES256-SHA`,
`AES128-CCM-8`,
`AES256-CCM-8`,
`ECDHE-ECDSA-AES128-CCM`,
`ECDHE-ECDSA-AES128-CCM-8`,
`ECDHE-ECDSA-AES256-CCM-8`,
`ECDHE-RSA-AES128-SHA`,
`ECDHE-RSA-AES256-SHA`,
`ECDHE-ECDSA-AES128-SHA`,
`ECDHE-ECDSA-AES256-SHA`,
`ECDHE-RSA-RC4-SHA`,
`ECDHE-RSA-DES-CBC3-SHA`,
`ECDHE-ECDSA-RC4-SHA`,
`ECDHE-ECDSA-DES-CBC3-SHA`,
`AES128-SHA256`,
`AES256-SHA256`,
`DHE-RSA-AES128-SHA256`,
`DHE-RSA-AES256-SHA256`,
`ECDH-RSA-AES128-SHA`,
`ECDH-RSA-AES256-SHA`,
`ECDH-ECDSA-AES128-SHA`,
`ECDH-ECDSA-AES256-SHA`,
`ECDH-RSA-RC4-SHA`,
`ECDH-RSA-DES-CBC3-SHA`,
`ECDH-ECDSA-RC4-SHA`,
`ECDH-ECDSA-DES-CBC3-SHA`,
`AES128-GCM-SHA256`,
`AES256-GCM-SHA384`,
`DHE-RSA-AES128-GCM-SHA256`,
`DHE-RSA-AES256-GCM-SHA384`,
`ECDHE-RSA-AES128-GCM-SHA256`,
`ECDHE-RSA-AES256-GCM-SHA384`,
`ECDHE-ECDSA-AES128-GCM-SHA256`,
`ECDHE-ECDSA-AES256-GCM-SHA384`,
`ECDH-RSA-AES128-GCM-SHA256`,
`ECDH-RSA-AES256-GCM-SHA384`,
`ECDH-ECDSA-AES128-GCM-SHA256`,
`ECDH-ECDSA-AES256-GCM-SHA384`,
`CAMELLIA128-SHA`,
`DHE-RSA-CAMELLIA128-SHA`,
`CAMELLIA256-SHA`,
`DHE-RSA-CAMELLIA256-SHA`,
`CAMELLIA128-SHA256`,
`DHE-RSA-CAMELLIA128-SHA256`,
`CAMELLIA256-SHA256`,
`DHE-RSA-CAMELLIA256-SHA256`,
`ECDHE-RSA-AES128-SHA256`,
`ECDHE-ECDSA-AES128-SHA256`,
`ECDH-RSA-AES128-SHA256`,
`ECDH-ECDSA-AES128-SHA256`,
`ECDHE-RSA-AES256-SHA384`,
`ECDHE-ECDSA-AES256-SHA384`,
`ECDH-RSA-AES256-SHA384`,
`ECDH-ECDSA-AES256-SHA384`,
`ECDHE-RSA-CHACHA20-POLY1305`,
`ECDHE-ECDSA-CHACHA20-POLY1305`,
`DHE-RSA-CHACHA20-POLY1305`,
`ECDHE-RSA-CHACHA20-POLY1305-OLD`,
`ECDHE-ECDSA-CHACHA20-POLY1305-OLD`,
`DHE-RSA-CHACHA20-POLY1305-OLD`,
`ADH-AES128-SHA`,
`QSH`,
`RENEGOTIATION-INFO`,
`IDEA-CBC-SHA`,
`ECDHE-ECDSA-NULL-SHA`,
`ECDHE-PSK-NULL-SHA256`,
`ECDHE-PSK-AES128-CBC-SHA256`,
`PSK-CHACHA20-POLY1305`,
`ECDHE-PSK-CHACHA20-POLY1305`,
`DHE-PSK-CHACHA20-POLY1305`,
`EDH-RSA-DES-CBC3-SHA`,

## Schannel

Schannel allows the enabling and disabling of encryption algorithms, but not
specific ciphersuites. They are
[defined](https://docs.microsoft.com/windows/desktop/SecCrypto/alg-id) by
Microsoft.

There is also the case that the selected algorithm is not supported by the
protocol or does not match the ciphers offered by the server during the SSL
negotiation. In this case curl will return error
`CURLE_SSL_CONNECT_ERROR (35) SEC_E_ALGORITHM_MISMATCH`
and the request will fail.

`CALG_MD2`,
`CALG_MD4`,
`CALG_MD5`,
`CALG_SHA`,
`CALG_SHA1`,
`CALG_MAC`,
`CALG_RSA_SIGN`,
`CALG_DSS_SIGN`,
`CALG_NO_SIGN`,
`CALG_RSA_KEYX`,
`CALG_DES`,
`CALG_3DES_112`,
`CALG_3DES`,
`CALG_DESX`,
`CALG_RC2`,
`CALG_RC4`,
`CALG_SEAL`,
`CALG_DH_SF`,
`CALG_DH_EPHEM`,
`CALG_AGREEDKEY_ANY`,
`CALG_HUGHES_MD5`,
`CALG_SKIPJACK`,
`CALG_TEK`,
`CALG_CYLINK_MEK`,
`CALG_SSL3_SHAMD5`,
`CALG_SSL3_MASTER`,
`CALG_SCHANNEL_MASTER_HASH`,
`CALG_SCHANNEL_MAC_KEY`,
`CALG_SCHANNEL_ENC_KEY`,
`CALG_PCT1_MASTER`,
`CALG_SSL2_MASTER`,
`CALG_TLS1_MASTER`,
`CALG_RC5`,
`CALG_HMAC`,
`CALG_TLS1PRF`,
`CALG_HASH_REPLACE_OWF`,
`CALG_AES_128`,
`CALG_AES_192`,
`CALG_AES_256`,
`CALG_AES`,
`CALG_SHA_256`,
`CALG_SHA_384`,
`CALG_SHA_512`,
`CALG_ECDH`,
`CALG_ECMQV`,
`CALG_ECDSA`,
`CALG_ECDH_EPHEM`,

As of curl 7.77.0, you can also pass `SCH_USE_STRONG_CRYPTO` as a cipher name
to [constrain the set of available ciphers as specified in the schannel
documentation](https://docs.microsoft.com/en-us/windows/win32/secauthn/tls-cipher-suites-in-windows-server-2022).
Note that the supported ciphers in this case follows the OS version, so if you
are running an outdated OS you might still be supporting weak ciphers.
# Decision making in the curl project

A rough guide to how we make decisions and who does what.

## BDFL

This project was started by and has to some extent been pushed forward over
the years with Daniel Stenberg as the driving force. It matches a standard
BDFL (Benevolent Dictator For Life) style project.

This setup has been used due to convenience and the fact that is has worked
fine this far. It is not because someone thinks of it as a superior project
leadership model. It will also only continue working as long as Daniel manages
to listen in to what the project and the general user population wants and
expects from us.

## Legal entity

There is no legal entity. The curl project is just a bunch of people scattered
around the globe with the common goal to produce source code that creates
great products. We are not part of any umbrella organization and we are not
located in any specific country. We are totally independent.

The copyrights in the project are owned by the individuals and organizations
that wrote those parts of the code.

## Decisions

The curl project is not a democracy, but everyone is entitled to state their
opinion and may argue for their sake within the community.

All and any changes that have been done or will be done are eligible to bring
up for discussion, to object to or to praise. Ideally, we find consensus for
the appropriate way forward in any given situation or challenge.

If there is no obvious consensus, a maintainer who's knowledgeable in the
specific area will take an "executive" decision that they think is the right
for the project.

## Donations

Donating plain money to curl is best done to curl's [Open Collective
fund](https://opencollective.com/curl). Open Collective is a US based
non-profit organization that holds on to funds for us. This fund is then used
for paying the curl security bug bounties, to reimburse project related
expenses etc.

Donations to the project can also come in form of server hosting, providing
services and paying for people to work on curl related code etc. Usually, such
donations are services paid for directly by the sponsors.

We grade sponsors in a few different levels and if they meet the criteria,
they can be mentioned on the Sponsors page on the curl website.

## Commercial Support

The curl project does not do or offer commercial support. It only hosts
mailing lists, runs bug trackers etc to facilitate communication and work.

However, Daniel works for wolfSSL and we offer commercial curl support there.

# Key roles

## User

Someone who uses or has used curl or libcurl.

## Contributor

Someone who has helped the curl project, who has contributed to bring it
forward. Contributing could be to provide advice, debug a problem, file a bug
report, run test infrastructure or writing code etc.

## Commit author

Sometimes also called 'committer'. Someone who has authored a commit in the
curl source code repository. Committers are recorded as `Author` in git.

## Maintainers

A maintainer in the curl project is an individual who has been given
permissions to push commits to one of the git repositories.

Maintainers are free to push commits to the repositories at their own will.
Maintainers are however expected to listen to feedback from users and any
change that is non-trivial in size or nature *should* be brought to the
project as a Pull-Request (PR) to allow others to comment/object before merge.

## Former maintainers

A maintainer who stops being active in the project will at some point get
their push permissions removed. We do this for security reasons but also to
make sure that we always have the list of maintainers as "the team that push
stuff to curl".

Getting push permissions removed is not a punishment. Everyone who ever worked
on maintaining curl is considered a hero, for all time hereafter.

## Security team members

We have a security team. That is the team of people who are subscribed to the
curl-security mailing list; the receivers of security reports from users and
developers. This list of people will vary over time but should be skilled
developers familiar with the curl project.

The security team works best when it consists of a small set of active
persons. We invite new members when the team seems to need it, and we also
expect to retire security team members as they "drift off" from the project or
just find themselves unable to perform their duties there.

## Server admins

We run a web server, a mailing list and more on the curl project's primary
server. That physical machine is owned and run by Haxx. Daniel is the primary
admin of all things curl related server stuff, but Björn Stenberg and Linus
Feltzing serve as backup admins for when Daniel is gone or unable.

The primary server is paid for by Haxx. The machine is physically located in a
server bunker in Stockholm Sweden, operated by the company Portlane.

The website contents are served to the web via Fastly and Daniel is the
primary curl contact with Fastly.

## BDFL

That is Daniel.

# Maintainers

A curl maintainer is a project volunteer who has the authority and rights to
merge changes into a git repository in the curl project.

Anyone can aspire to become a curl maintainer.

### Duties

There are no mandatory duties. We hope and wish that maintainers consider
reviewing patches and help merging them, especially when the changes are
within the area of personal expertise and experience.

### Requirements

- only merge code that meets our quality and style guide requirements.
- *never* merge code without doing a PR first, unless the change is "trivial"
- if in doubt, ask for input/feedback from others

### Recommendations

- we require two-factor authentication enabled on your GitHub account to
  reduce risk of malicious source code tampering
- consider enabling signed git commits for additional verification of changes

### Merge advice

When you are merging patches/PRs...

- make sure the commit messages follow our template
- squash patch sets into a few logical commits even if the PR did not, if
  necessary
- avoid the "merge" button on GitHub, do it "manually" instead to get full
  control and full audit trail (github leaves out you as "Committer:")
- remember to credit the reporter and the helpers.

## Who are maintainers?

The [list of maintainers](https://github.com/orgs/curl/people). Be aware that
the level of presence and activity in the project vary greatly between
different individuals and over time.

### Become a maintainer?

If you think you can help making the project better by shouldering some
maintaining responsibilities, then please get in touch.

You will be expected to be familiar with the curl project and its ways of
working. You need to have gotten a few quality patches merged as a proof of
this.

### Stop being a maintainer

If you (appear to) not be active in the project anymore, you may be removed as
a maintainer. Thank you for your service!
# Features -- what curl can do

## curl tool

 - config file support
 - multiple URLs in a single command line
 - range "globbing" support: [0-13], {one,two,three}
 - multiple file upload on a single command line
 - custom maximum transfer rate
 - redirectable stderr
 - parallel transfers

## libcurl

 - full URL syntax with no length limit
 - custom maximum download time
 - custom least download speed acceptable
 - custom output result after completion
 - guesses protocol from host name unless specified
 - uses .netrc
 - progress bar with time statistics while downloading
 - "standard" proxy environment variables support
 - compiles on win32 (reported builds on 70+ operating systems)
 - selectable network interface for outgoing traffic
 - IPv6 support on unix and Windows
 - happy eyeballs dual-stack connects
 - persistent connections
 - SOCKS 4 + 5 support, with or without local name resolving
 - supports user name and password in proxy environment variables
 - operations through HTTP proxy "tunnel" (using CONNECT)
 - replaceable memory functions (malloc, free, realloc, etc)
 - asynchronous name resolving (6)
 - both a push and a pull style interface
 - international domain names (11)

## HTTP

 - HTTP/0.9 responses are optionally accepted
 - HTTP/1.0
 - HTTP/1.1
 - HTTP/2, including multiplexing and server push (5)
 - GET
 - PUT
 - HEAD
 - POST
 - multipart formpost (RFC1867-style)
 - authentication: Basic, Digest, NTLM (9) and Negotiate (SPNEGO) (3)
   to server and proxy
 - resume (both GET and PUT)
 - follow redirects
 - maximum amount of redirects to follow
 - custom HTTP request
 - cookie get/send fully parsed
 - reads/writes the netscape cookie file format
 - custom headers (replace/remove internally generated headers)
 - custom user-agent string
 - custom referrer string
 - range
 - proxy authentication
 - time conditions
 - via HTTP proxy, HTTPS proxy or SOCKS proxy
 - retrieve file modification date
 - Content-Encoding support for deflate and gzip
 - "Transfer-Encoding: chunked" support in uploads
 - automatic data compression (12)

## HTTPS (1)

 - (all the HTTP features)
 - HTTP/3 experimental support
 - using client certificates
 - verify server certificate
 - via HTTP proxy, HTTPS proxy or SOCKS proxy
 - select desired encryption
 - select usage of a specific SSL version

## FTP

 - download
 - authentication
 - Kerberos 5 (13)
 - active/passive using PORT, EPRT, PASV or EPSV
 - single file size information (compare to HTTP HEAD)
 - 'type=' URL support
 - dir listing
 - dir listing names-only
 - upload
 - upload append
 - upload via http-proxy as HTTP PUT
 - download resume
 - upload resume
 - custom ftp commands (before and/or after the transfer)
 - simple "range" support
 - via HTTP proxy, HTTPS proxy or SOCKS proxy
 - all operations can be tunneled through proxy
 - customizable to retrieve file modification date
 - no dir depth limit

## FTPS (1)

 - implicit `ftps://` support that use SSL on both connections
 - explicit "AUTH TLS" and "AUTH SSL" usage to "upgrade" plain `ftp://`
   connection to use SSL for both or one of the connections

## SCP (8)

 - both password and public key auth

## SFTP (7)

 - both password and public key auth
 - with custom commands sent before/after the transfer

## TFTP

 - download
 - upload

## TELNET

 - connection negotiation
 - custom telnet options
 - stdin/stdout I/O

## LDAP (2)

 - full LDAP URL support

## DICT

 - extended DICT URL support

## FILE

 - URL support
 - upload
 - resume

## SMB

 - SMBv1 over TCP and SSL
 - download
 - upload
 - authentication with NTLMv1

## SMTP

 - authentication: Plain, Login, CRAM-MD5, Digest-MD5, NTLM (9), Kerberos 5
   (4) and External.
 - send emails
 - mail from support
 - mail size support
 - mail auth support for trusted server-to-server relaying
 - multiple recipients
 - via http-proxy

## SMTPS (1)

 - implicit `smtps://` support
 - explicit "STARTTLS" usage to "upgrade" plain `smtp://` connections to use SSL
 - via http-proxy

## POP3

 - authentication: Clear Text, APOP and SASL
 - SASL based authentication: Plain, Login, CRAM-MD5, Digest-MD5, NTLM (9),
   Kerberos 5 (4) and External.
 - list emails
 - retrieve emails
 - enhanced command support for: CAPA, DELE, TOP, STAT, UIDL and NOOP via
   custom requests
 - via http-proxy

## POP3S (1)

 - implicit `pop3s://` support
 - explicit "STLS" usage to "upgrade" plain `pop3://` connections to use SSL
 - via http-proxy

## IMAP

 - authentication: Clear Text and SASL
 - SASL based authentication: Plain, Login, CRAM-MD5, Digest-MD5, NTLM (9),
   Kerberos 5 (4) and External.
 - list the folders of a mailbox
 - select a mailbox with support for verifying the UIDVALIDITY
 - fetch emails with support for specifying the UID and SECTION
 - upload emails via the append command
 - enhanced command support for: EXAMINE, CREATE, DELETE, RENAME, STATUS,
   STORE, COPY and UID via custom requests
 - via http-proxy

## IMAPS (1)

 - implicit `imaps://` support
 - explicit "STARTTLS" usage to "upgrade" plain `imap://` connections to use SSL
 - via http-proxy

## MQTT

 - Subscribe to and publish topics using url scheme `mqtt://broker/topic`

## Footnotes

  1. requires a TLS library
  2. requires OpenLDAP or WinLDAP
  3. requires a GSS-API implementation (such as Heimdal or MIT Kerberos) or
     SSPI (native Windows)
  4. requires a GSS-API implementation, however, only Windows SSPI is
     currently supported
  5. requires nghttp2
  6. requires c-ares
  7. requires libssh2, libssh or wolfSSH
  8. requires libssh2 or libssh
  9. requires OpenSSL, GnuTLS, mbedTLS, NSS, yassl, Secure Transport or SSPI
     (native Windows)
  10. -
  11. requires libidn2 or Windows
  12. requires libz, brotli and/or zstd
  13. requires a GSS-API implementation (such as Heimdal or MIT Kerberos)
Contributor Code of Conduct
===========================

As contributors and maintainers of this project, we pledge to respect all
people who contribute through reporting issues, posting feature requests,
updating documentation, submitting pull requests or patches, and other
activities.

We are committed to making participation in this project a harassment-free
experience for everyone, regardless of level of experience, gender, gender
identity and expression, sexual orientation, disability, personal appearance,
body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual
language or imagery, derogatory comments or personal attacks, trolling, public
or private harassment, insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct. Project maintainers who do not
follow the Code of Conduct may be removed from the project team.

This code of conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community.

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by opening an issue or contacting one or more of the project
maintainers.

This Code of Conduct is adapted from the [Contributor
Covenant](https://contributor-covenant.org/), version 1.1.0, available at
[https://contributor-covenant.org/version/1/1/0/](https://contributor-covenant.org/version/1/1/0/)
# Contributing to the curl project

This document is intended to offer guidelines on how to best contribute to the
curl project. This concerns new features as well as corrections to existing
flaws or bugs.

## Learning curl

### Join the Community

Skip over to [https://curl.se/mail/](https://curl.se/mail/) and join
the appropriate mailing list(s). Read up on details before you post
questions. Read this file before you start sending patches. We prefer
questions sent to and discussions being held on the mailing list(s), not sent
to individuals.

Before posting to one of the curl mailing lists, please read up on the
[mailing list etiquette](https://curl.se/mail/etiquette.html).

We also hang out on IRC in #curl on libera.chat

If you are at all interested in the code side of things, consider clicking
'watch' on the [curl repo on GitHub](https://github.com/curl/curl) to be
notified of pull requests and new issues posted there.

### License and copyright

When contributing with code, you agree to put your changes and new code under
the same license curl and libcurl is already using unless stated and agreed
otherwise.

If you add a larger piece of code, you can opt to make that file or set of
files to use a different license as long as they do not enforce any changes to
the rest of the package and they make sense. Such "separate parts" can not be
GPL licensed (as we do not want copyleft to affect users of libcurl) but they
must use "GPL compatible" licenses (as we want to allow users to use libcurl
properly in GPL licensed environments).

When changing existing source code, you do not alter the copyright of the
original file(s). The copyright will still be owned by the original creator(s)
or those who have been assigned copyright by the original author(s).

By submitting a patch to the curl project, you are assumed to have the right
to the code and to be allowed by your employer or whatever to hand over that
patch/code to us. We will credit you for your changes as far as possible, to
give credit but also to keep a trace back to who made what changes. Please
always provide us with your full real name when contributing,

### What To Read

Source code, the man pages, the [INTERNALS
document](https://curl.se/dev/internals.html),
[TODO](https://curl.se/docs/todo.html),
[KNOWN_BUGS](https://curl.se/docs/knownbugs.html) and the [most recent
changes](https://curl.se/dev/sourceactivity.html) in git. Just lurking on
the [curl-library mailing
list](https://curl.se/mail/list.cgi?list=curl-library) will give you a
lot of insights on what's going on right now. Asking there is a good idea too.

## Write a good patch

### Follow code style

When writing C code, follow the
[CODE_STYLE](https://curl.se/dev/code-style.html) already established in
the project. Consistent style makes code easier to read and mistakes less
likely to happen. Run `make checksrc` before you submit anything, to make sure
you follow the basic style. That script does not verify everything, but if it
complains you know you have work to do.

### Non-clobbering All Over

When you write new functionality or fix bugs, it is important that you do not
fiddle all over the source files and functions. Remember that it is likely
that other people have done changes in the same source files as you have and
possibly even in the same functions. If you bring completely new
functionality, try writing it in a new source file. If you fix bugs, try to
fix one bug at a time and send them as separate patches.

### Write Separate Changes

It is annoying when you get a huge patch from someone that is said to fix 511
odd problems, but discussions and opinions do not agree with 510 of them - or
509 of them were already fixed in a different way. Then the person merging
this change needs to extract the single interesting patch from somewhere
within the huge pile of source, and that creates a lot of extra work.

Preferably, each fix that corrects a problem should be in its own patch/commit
with its own description/commit message stating exactly what they correct so
that all changes can be selectively applied by the maintainer or other
interested parties.

Also, separate changes enable bisecting much better for tracking problems
and regression in the future.

### Patch Against Recent Sources

Please try to get the latest available sources to make your patches against.
It makes the lives of the developers so much easier. The best is if you get
the most up-to-date sources from the git repository, but the latest release
archive is quite OK as well.

### Documentation

Writing docs is dead boring and one of the big problems with many open source
projects. But someone's gotta do it. It makes things a lot easier if you
submit a small description of your fix or your new features with every
contribution so that it can be swiftly added to the package documentation.

The documentation is always made in man pages (nroff formatted) or plain
ASCII files. All HTML files on the website and in the release archives are
generated from the nroff/ASCII versions.

### Test Cases

Since the introduction of the test suite, we can quickly verify that the main
features are working as they are supposed to. To maintain this situation and
improve it, all new features and functions that are added need to be tested
in the test suite. Every feature that is added should get at least one valid
test case that verifies that it works as documented. If every submitter also
posts a few test cases, it will not end up as a heavy burden on a single person!

If you do not have test cases or perhaps you have done something that is hard
to write tests for, do explain exactly how you have otherwise tested and
verified your changes.

## Sharing Your Changes

### How to get your changes into the main sources

Ideally you file a [pull request on
GitHub](https://github.com/curl/curl/pulls), but you can also send your plain
patch to [the curl-library mailing
list](https://curl.se/mail/list.cgi?list=curl-library).

Either way, your change will be reviewed and discussed there and you will be
expected to correct flaws pointed out and update accordingly, or the change
risks stalling and eventually just getting deleted without action. As a
submitter of a change, you are the owner of that change until it has been merged.

Respond on the list or on github about the change and answer questions and/or
fix nits/flaws. This is important. We will take lack of replies as a sign that
you are not anxious to get your patch accepted and we tend to simply drop such
changes.

### About pull requests

With github it is easy to send a [pull
request](https://github.com/curl/curl/pulls) to the curl project to have
changes merged.

We strongly prefer pull requests to mailed patches, as it makes it a proper
git commit that is easy to merge and they are easy to track and not that easy
to lose in the flood of many emails, like they sometimes do on the mailing
lists.

Every pull request submitted will automatically be
tested in several different ways. [See CI.md for more
information](https://github.com/curl/curl/blob/master/tests/CI.md).

Sometimes the tests fail due to a dependency service temporarily being offline
or otherwise unavailable, eg. package downloads. In this case you can just
try to update your pull requests to rerun the tests later as described below.

You can update your pull requests by pushing new commits or force-pushing
changes to existing commits. Force-pushing an amended commit without any
actual content changed also allows you to retrigger the tests for that commit.

When you adjust your pull requests after review, consider squashing the
commits so that we can review the full updated version more easily.

### Making quality patches

Make the patch against as recent source versions as possible.

If you have followed the tips in this document and your patch still has not been
incorporated or responded to after some weeks, consider resubmitting it to the
list or better yet: change it to a pull request.

### Write good commit messages

A short guide to how to write commit messages in the curl project.

    ---- start ----
    [area]: [short line describing the main effect]
           -- empty line --
    [full description, no wider than 72 columns that describe as much as
    possible as to why this change is made, and possibly what things
    it fixes and everything else that is related]
           -- empty line --
    [Closes/Fixes #1234 - if this closes or fixes a github issue]
    [Bug: URL to source of the report or more related discussion]
    [Reported-by: John Doe - credit the reporter]
    [whatever-else-by: credit all helpers, finders, doers]
    ---- stop ----

The first line is a succinct description of the change:

 - use the imperative, present tense: "change" not "changed" nor "changes"
 - do not capitalize first letter
 - no dot (.) at the end

The `[area]` in the first line can be `http2`, `cookies`, `openssl` or
similar. There's no fixed list to select from but using the same "area" as
other related changes could make sense.

Do not forget to use commit --author="" if you commit someone else's work, and
make sure that you have your own user and email setup correctly in git before
you commit

### Write Access to git Repository

If you are a frequent contributor, you may be given push access to the git
repository and then you will be able to push your changes straight into the git
repo instead of sending changes as pull requests or by mail as patches.

Just ask if this is what you would want. You will be required to have posted
several high quality patches first, before you can be granted push access.

### How To Make a Patch with git

You need to first checkout the repository:

    git clone https://github.com/curl/curl.git

You then proceed and edit all the files you like and you commit them to your
local repository:

    git commit [file]

As usual, group your commits so that you commit all changes at once that
constitute a logical change.

Once you have done all your commits and you are happy with what you see, you
can make patches out of your changes that are suitable for mailing:

    git format-patch remotes/origin/master

This creates files in your local directory named NNNN-[name].patch for each
commit.

Now send those patches off to the curl-library list. You can of course opt to
do that with the 'git send-email' command.

### How To Make a Patch without git

Keep a copy of the unmodified curl sources. Make your changes in a separate
source tree. When you think you have something that you want to offer the
curl community, use GNU diff to generate patches.

If you have modified a single file, try something like:

    diff -u unmodified-file.c my-changed-one.c > my-fixes.diff

If you have modified several files, possibly in different directories, you
can use diff recursively:

    diff -ur curl-original-dir curl-modified-sources-dir > my-fixes.diff

The GNU diff and GNU patch tools exist for virtually all platforms, including
all kinds of Unixes and Windows:

For unix-like operating systems:

 - [https://savannah.gnu.org/projects/patch/](https://savannah.gnu.org/projects/patch/)
 - [https://www.gnu.org/software/diffutils/](https://www.gnu.org/software/diffutils/)

For Windows:

 - [https://gnuwin32.sourceforge.io/packages/patch.htm](https://gnuwin32.sourceforge.io/packages/patch.htm)
 - [https://gnuwin32.sourceforge.io/packages/diffutils.htm](https://gnuwin32.sourceforge.io/packages/diffutils.htm)

### Useful resources
 - [Webinar on getting code into cURL](https://www.youtube.com/watch?v=QmZ3W1d6LQI)
# Code defines to disable features and protocols

## CURL_DISABLE_ALTSVC

Disable support for Alt-Svc: HTTP headers.

## CURL_DISABLE_COOKIES

Disable support for HTTP cookies.

## CURL_DISABLE_CRYPTO_AUTH

Disable support for authentication methods using crypto.

## CURL_DISABLE_DICT

Disable the DICT protocol

## CURL_DISABLE_DOH

Disable DNS-over-HTTPS

## CURL_DISABLE_FILE

Disable the FILE protocol

## CURL_DISABLE_FTP

Disable the FTP (and FTPS) protocol

## CURL_DISABLE_GETOPTIONS

Disable the `curl_easy_options` API calls that lets users get information
about existing options to `curl_easy_setopt`.

## CURL_DISABLE_GOPHER

Disable the GOPHER protocol.

## CURL_DISABLE_HSTS

Disable the HTTP Strict Transport Security support.

## CURL_DISABLE_HTTP

Disable the HTTP(S) protocols. Note that this then also disable HTTP proxy
support.

## CURL_DISABLE_HTTP_AUTH

Disable support for all HTTP authentication methods.

## CURL_DISABLE_IMAP

Disable the IMAP(S) protocols.

## CURL_DISABLE_LDAP

Disable the LDAP(S) protocols.

## CURL_DISABLE_LDAPS

Disable the LDAPS protocol.

## CURL_DISABLE_LIBCURL_OPTION

Disable the --libcurl option from the curl tool.

## CURL_DISABLE_MIME

Disable MIME support.

## CURL_DISABLE_MQTT

Disable MQTT support.

## CURL_DISABLE_NETRC

Disable the netrc parser.

## CURL_DISABLE_NTLM

Disable support for NTLM.

## CURL_DISABLE_OPENSSL_AUTO_LOAD_CONFIG

Disable the auto load config support in the OpenSSL backend.

## CURL_DISABLE_PARSEDATE

Disable date parsing

## CURL_DISABLE_POP3

Disable the POP3 protocol

## CURL_DISABLE_PROGRESS_METER

Disable the built-in progress meter

## CURL_DISABLE_PROXY

Disable support for proxies

## CURL_DISABLE_RTSP

Disable the RTSP protocol.

## CURL_DISABLE_SHUFFLE_DNS

Disable the shuffle DNS feature

## CURL_DISABLE_SMB

Disable the SMB(S) protocols

## CURL_DISABLE_SMTP

Disable the SMTP(S) protocols

## CURL_DISABLE_SOCKETPAIR

Disable the use of socketpair internally to allow waking up and canceling
curl_multi_poll().

## CURL_DISABLE_TELNET

Disable the TELNET protocol

## CURL_DISABLE_TFTP

Disable the TFTP protocol

## CURL_DISABLE_VERBOSE_STRINGS

Disable verbose strings and error messages.
# HTTP3 (and QUIC)

## Resources

[HTTP/3 Explained](https://http3-explained.haxx.se/en/) - the online free
book describing the protocols involved.

[QUIC implementation](https://github.com/curl/curl/wiki/QUIC-implementation) -
the wiki page describing the plan for how to support QUIC and HTTP/3 in curl
and libcurl.

[quicwg.org](https://quicwg.org/) - home of the official protocol drafts

## QUIC libraries

QUIC libraries we are experimenting with:

[ngtcp2](https://github.com/ngtcp2/ngtcp2)

[quiche](https://github.com/cloudflare/quiche)

## Experimental

HTTP/3 and QUIC support in curl is considered **EXPERIMENTAL** until further
notice. It needs to be enabled at build-time.

Further development and tweaking of the HTTP/3 support in curl will happen in
in the master branch using pull-requests, just like ordinary changes.

# ngtcp2 version

## Build with OpenSSL

Build (patched) OpenSSL

     % git clone --depth 1 -b openssl-3.0.0+quic https://github.com/quictls/openssl
     % cd openssl
     % ./config enable-tls1_3 --prefix=<somewhere1>
     % make
     % make install

Build nghttp3

     % cd ..
     % git clone https://github.com/ngtcp2/nghttp3
     % cd nghttp3
     % autoreconf -fi
     % ./configure --prefix=<somewhere2> --enable-lib-only
     % make
     % make install

Build ngtcp2

     % cd ..
     % git clone https://github.com/ngtcp2/ngtcp2
     % cd ngtcp2
     % autoreconf -fi
     % ./configure PKG_CONFIG_PATH=<somewhere1>/lib/pkgconfig:<somewhere2>/lib/pkgconfig LDFLAGS="-Wl,-rpath,<somewhere1>/lib" --prefix=<somewhere3> --enable-lib-only
     % make
     % make install

Build curl

     % cd ..
     % git clone https://github.com/curl/curl
     % cd curl
     % autoreconf -fi
     % LDFLAGS="-Wl,-rpath,<somewhere1>/lib" ./configure --with-openssl=<somewhere1> --with-nghttp3=<somewhere2> --with-ngtcp2=<somewhere3>
     % make
     % make install

For OpenSSL 3.0.0 or later builds on Linux for x86_64 architecture, substitute all occurrences of "/lib" with "/lib64"

## Build with GnuTLS

Build GnuTLS

     % git clone --depth 1 https://gitlab.com/gnutls/gnutls.git
     % cd gnutls
     % ./bootstrap
     % ./configure --prefix=<somewhere1>
     % make
     % make install

Build nghttp3

     % cd ..
     % git clone https://github.com/ngtcp2/nghttp3
     % cd nghttp3
     % autoreconf -fi
     % ./configure --prefix=<somewhere2> --enable-lib-only
     % make
     % make install

Build ngtcp2

     % cd ..
     % git clone https://github.com/ngtcp2/ngtcp2
     % cd ngtcp2
     % autoreconf -fi
     % ./configure PKG_CONFIG_PATH=<somewhere1>/lib/pkgconfig:<somewhere2>/lib/pkgconfig LDFLAGS="-Wl,-rpath,<somewhere1>/lib" --prefix=<somewhere3> --enable-lib-only --with-gnutls
     % make
     % make install

Build curl

     % cd ..
     % git clone https://github.com/curl/curl
     % cd curl
     % autoreconf -fi
     % ./configure --without-openssl --with-gnutls=<somewhere1> --with-nghttp3=<somewhere2> --with-ngtcp2=<somewhere3>
     % make
     % make install

# quiche version

## build

Build quiche and BoringSSL:

     % git clone --recursive https://github.com/cloudflare/quiche
     % cd quiche
     % cargo build --package quiche --release --features ffi,pkg-config-meta,qlog
     % mkdir quiche/deps/boringssl/src/lib
     % ln -vnf $(find target/release -name libcrypto.a -o -name libssl.a) quiche/deps/boringssl/src/lib/

Build curl:

     % cd ..
     % git clone https://github.com/curl/curl
     % cd curl
     % autoreconf -fi
     % ./configure LDFLAGS="-Wl,-rpath,$PWD/../quiche/target/release" --with-openssl=$PWD/../quiche/quiche/deps/boringssl/src --with-quiche=$PWD/../quiche/target/release
     % make
     % make install

 If `make install` results in `Permission denied` error, you will need to prepend it with `sudo`.

# `--http3`

Use HTTP/3 directly:

    curl --http3 https://nghttp2.org:4433/

Upgrade via Alt-Svc:

    curl --alt-svc altsvc.cache https://quic.aiortc.org/

See this [list of public HTTP/3 servers](https://bagder.github.io/HTTP3-test/)

## Known Bugs

Check out the [list of known HTTP3 bugs](https://curl.se/docs/knownbugs.html#HTTP3).

# HTTP/3 Test server

This is not advice on how to run anything in production. This is for
development and experimenting.

## Preqreqs

An existing local HTTP/1.1 server that hosts files. Preferably also a few huge
ones.  You can easily create huge local files like `truncate -s=8G 8GB` - they
are huge but do not occupy that much space on disk since they're just a big
hole.

In my Debian setup I just installed **apache2**. It runs on port 80 and has a
document root in `/var/www/html`. I can get the 8GB file from it with `curl
localhost/8GB -o dev/null`

In this description we setup and run a HTTP/3 reverse-proxy in front of the
HTTP/1 server.

## Setup

You can select either or both of these server solutions.

### nghttpx

Get, build and install **quictls**, **nghttp3** and **ngtcp2** as described
above.

Get, build and install **nghttp2**:

    git clone https://github.com/nghttp2/nghttp2.git
    cd nghttp2
    autoreconf -fi
    PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/home/daniel/build-quictls/lib/pkgconfig:/home/daniel/build-nghttp3/lib/pkgconfig:/home/daniel/build-ngtcp2/lib/pkgconfig  LDFLAGS=-L/home/daniel/build-quictls/lib CFLAGS=-I/home/daniel/build-quictls/include ./configure --enable-maintainer-mode --prefix=/home/daniel/build-nghttp2 --disable-shared --enable-app --enable-http3 --without-jemalloc --without-libxml2 --without-systemd
    make && make install

Run the local h3 server on port 9443, make it proxy all traffic through to
HTTP/1 on localhost port 80. For local toying, we can just use the test cert
that exists in curl's test dir.

    CERT=$CURLSRC/tests/stunnel.pem
    $HOME/bin/nghttpx $CERT $CERT --backend=localhost,80 \
      --frontend="localhost,9443;quic"

### Caddy

[Install caddy](https://caddyserver.com/docs/install), you can even put the
single binary in a separate directory if you prefer.

In the same directory you put caddy, create a `Caddyfile` with the following
content to run a HTTP/3 reverse-proxy on port 7443:
~~~
{
    auto_https disable_redirects
	servers :7443 {
		protocol {
			experimental_http3
		}
	}
}

localhost:7443 {
	reverse_proxy localhost:80
}
~~~

Then run caddy:

    ./caddy start
# Rustls

[Rustls is a TLS backend written in Rust.](https://docs.rs/rustls/). Curl can
be built to use it as an alternative to OpenSSL or other TLS backends. We use
the [rustls-ffi C bindings](https://github.com/rustls/rustls-ffi/). This
version of curl depends on version v0.8.2 of rustls-ffi.

# Building with rustls

First, [install Rust](https://rustup.rs/).

Next, check out, build, and install the appropriate version of rustls-ffi:

    % cargo install cbindgen
    % git clone https://github.com/rustls/rustls-ffi -b v0.8.2
    % cd rustls-ffi
    % make
    % make DESTDIR=${HOME}/rustls-ffi-built/ install

Now configure and build curl with rustls:

    % git clone https://github.com/curl/curl
    % cd curl
    % ./buildconf
    % ./configure --with-rustls=${HOME}/rustls-ffi-built
    % make
# curl tutorial

## Simple Usage

Get the main page from a web-server:

    curl https://www.example.com/

Get the README file the user's home directory at funet's ftp-server:

    curl ftp://ftp.funet.fi/README

Get a web page from a server using port 8000:

    curl http://www.weirdserver.com:8000/

Get a directory listing of an FTP site:

    curl ftp://ftp.funet.fi

Get the definition of curl from a dictionary:

    curl dict://dict.org/m:curl

Fetch two documents at once:

    curl ftp://ftp.funet.fi/ http://www.weirdserver.com:8000/

Get a file off an FTPS server:

    curl ftps://files.are.secure.com/secrets.txt

or use the more appropriate FTPS way to get the same file:

    curl --ftp-ssl ftp://files.are.secure.com/secrets.txt

Get a file from an SSH server using SFTP:

    curl -u username sftp://example.com/etc/issue

Get a file from an SSH server using SCP using a private key (not
password-protected) to authenticate:

    curl -u username: --key ~/.ssh/id_rsa scp://example.com/~/file.txt

Get a file from an SSH server using SCP using a private key
(password-protected) to authenticate:

    curl -u username: --key ~/.ssh/id_rsa --pass private_key_password
    scp://example.com/~/file.txt

Get the main page from an IPv6 web server:

    curl "http://[2001:1890:1112:1::20]/"

Get a file from an SMB server:

    curl -u "domain\username:passwd" smb://server.example.com/share/file.txt

## Download to a File

Get a web page and store in a local file with a specific name:

    curl -o thatpage.html http://www.example.com/

Get a web page and store in a local file, make the local file get the name of
the remote document (if no file name part is specified in the URL, this will
fail):

    curl -O http://www.example.com/index.html

Fetch two files and store them with their remote names:

    curl -O www.haxx.se/index.html -O curl.se/download.html

## Using Passwords

### FTP

To ftp files using name+passwd, include them in the URL like:

    curl ftp://name:passwd@machine.domain:port/full/path/to/file

or specify them with the -u flag like

    curl -u name:passwd ftp://machine.domain:port/full/path/to/file

### FTPS

It is just like for FTP, but you may also want to specify and use SSL-specific
options for certificates etc.

Note that using `FTPS://` as prefix is the "implicit" way as described in the
standards while the recommended "explicit" way is done by using FTP:// and the
`--ftp-ssl` option.

### SFTP / SCP

This is similar to FTP, but you can use the `--key` option to specify a
private key to use instead of a password. Note that the private key may itself
be protected by a password that is unrelated to the login password of the
remote system; this password is specified using the `--pass` option.
Typically, curl will automatically extract the public key from the private key
file, but in cases where curl does not have the proper library support, a
matching public key file must be specified using the `--pubkey` option.

### HTTP

Curl also supports user and password in HTTP URLs, thus you can pick a file
like:

    curl http://name:passwd@machine.domain/full/path/to/file

or specify user and password separately like in

    curl -u name:passwd http://machine.domain/full/path/to/file

HTTP offers many different methods of authentication and curl supports
several: Basic, Digest, NTLM and Negotiate (SPNEGO). Without telling which
method to use, curl defaults to Basic. You can also ask curl to pick the most
secure ones out of the ones that the server accepts for the given URL, by
using `--anyauth`.

**Note**! According to the URL specification, HTTP URLs can not contain a user
and password, so that style will not work when using curl via a proxy, even
though curl allows it at other times. When using a proxy, you _must_ use the
`-u` style for user and password.

### HTTPS

Probably most commonly used with private certificates, as explained below.

## Proxy

curl supports both HTTP and SOCKS proxy servers, with optional authentication.
It does not have special support for FTP proxy servers since there are no
standards for those, but it can still be made to work with many of them. You
can also use both HTTP and SOCKS proxies to transfer files to and from FTP
servers.

Get an ftp file using an HTTP proxy named my-proxy that uses port 888:

    curl -x my-proxy:888 ftp://ftp.leachsite.com/README

Get a file from an HTTP server that requires user and password, using the
same proxy as above:

    curl -u user:passwd -x my-proxy:888 http://www.get.this/

Some proxies require special authentication. Specify by using -U as above:

    curl -U user:passwd -x my-proxy:888 http://www.get.this/

A comma-separated list of hosts and domains which do not use the proxy can be
specified as:

    curl --noproxy localhost,get.this -x my-proxy:888 http://www.get.this/

If the proxy is specified with `--proxy1.0` instead of `--proxy` or `-x`, then
curl will use HTTP/1.0 instead of HTTP/1.1 for any `CONNECT` attempts.

curl also supports SOCKS4 and SOCKS5 proxies with `--socks4` and `--socks5`.

See also the environment variables Curl supports that offer further proxy
control.

Most FTP proxy servers are set up to appear as a normal FTP server from the
client's perspective, with special commands to select the remote FTP server.
curl supports the `-u`, `-Q` and `--ftp-account` options that can be used to
set up transfers through many FTP proxies. For example, a file can be uploaded
to a remote FTP server using a Blue Coat FTP proxy with the options:

    curl -u "username@ftp.server Proxy-Username:Remote-Pass"
      --ftp-account Proxy-Password --upload-file local-file
      ftp://my-ftp.proxy.server:21/remote/upload/path/

See the manual for your FTP proxy to determine the form it expects to set up
transfers, and curl's `-v` option to see exactly what curl is sending.

## Ranges

HTTP 1.1 introduced byte-ranges. Using this, a client can request to get only
one or more subparts of a specified document. Curl supports this with the `-r`
flag.

Get the first 100 bytes of a document:

    curl -r 0-99 http://www.get.this/

Get the last 500 bytes of a document:

    curl -r -500 http://www.get.this/

Curl also supports simple ranges for FTP files as well. Then you can only
specify start and stop position.

Get the first 100 bytes of a document using FTP:

    curl -r 0-99 ftp://www.get.this/README

## Uploading

### FTP / FTPS / SFTP / SCP

Upload all data on stdin to a specified server:

    curl -T - ftp://ftp.upload.com/myfile

Upload data from a specified file, login with user and password:

    curl -T uploadfile -u user:passwd ftp://ftp.upload.com/myfile

Upload a local file to the remote site, and use the local file name at the
remote site too:

    curl -T uploadfile -u user:passwd ftp://ftp.upload.com/

Upload a local file to get appended to the remote file:

    curl -T localfile -a ftp://ftp.upload.com/remotefile

Curl also supports ftp upload through a proxy, but only if the proxy is
configured to allow that kind of tunneling. If it does, you can run curl in a
fashion similar to:

    curl --proxytunnel -x proxy:port -T localfile ftp.upload.com

### SMB / SMBS

    curl -T file.txt -u "domain\username:passwd"
      smb://server.example.com/share/

### HTTP

Upload all data on stdin to a specified HTTP site:

    curl -T - http://www.upload.com/myfile

Note that the HTTP server must have been configured to accept PUT before this
can be done successfully.

For other ways to do HTTP data upload, see the POST section below.

## Verbose / Debug

If curl fails where it is not supposed to, if the servers do not let you in, if
you cannot understand the responses: use the `-v` flag to get verbose
fetching. Curl will output lots of info and what it sends and receives in
order to let the user see all client-server interaction (but it will not show you
the actual data).

    curl -v ftp://ftp.upload.com/

To get even more details and information on what curl does, try using the
`--trace` or `--trace-ascii` options with a given file name to log to, like
this:

    curl --trace trace.txt www.haxx.se


## Detailed Information

Different protocols provide different ways of getting detailed information
about specific files/documents. To get curl to show detailed information about
a single file, you should use `-I`/`--head` option. It displays all available
info on a single file for HTTP and FTP. The HTTP information is a lot more
extensive.

For HTTP, you can get the header information (the same as `-I` would show)
shown before the data by using `-i`/`--include`. Curl understands the
`-D`/`--dump-header` option when getting files from both FTP and HTTP, and it
will then store the headers in the specified file.

Store the HTTP headers in a separate file (headers.txt in the example):

      curl --dump-header headers.txt curl.se

Note that headers stored in a separate file can be useful at a later time if
you want curl to use cookies sent by the server. More about that in the
cookies section.

## POST (HTTP)

It's easy to post data using curl. This is done using the `-d <data>` option.
The post data must be urlencoded.

Post a simple "name" and "phone" guestbook.

    curl -d "name=Rafael%20Sagula&phone=3320780" http://www.where.com/guest.cgi

How to post a form with curl, lesson #1:

Dig out all the `<input>` tags in the form that you want to fill in.

If there's a "normal" post, you use `-d` to post. `-d` takes a full "post
string", which is in the format

    <variable1>=<data1>&<variable2>=<data2>&...

The 'variable' names are the names set with `"name="` in the `<input>` tags,
and the data is the contents you want to fill in for the inputs. The data
*must* be properly URL encoded. That means you replace space with + and that
you replace weird letters with %XX where XX is the hexadecimal representation
of the letter's ASCII code.

Example:

(page located at `http://www.formpost.com/getthis/`)

```html
<form action="post.cgi" method="post">
<input name=user size=10>
<input name=pass type=password size=10>
<input name=id type=hidden value="blablabla">
<input name=ding value="submit">
</form>
```

We want to enter user 'foobar' with password '12345'.

To post to this, you enter a curl command line like:

    curl -d "user=foobar&pass=12345&id=blablabla&ding=submit"
      http://www.formpost.com/getthis/post.cgi

While `-d` uses the application/x-www-form-urlencoded mime-type, generally
understood by CGI's and similar, curl also supports the more capable
multipart/form-data type. This latter type supports things like file upload.

`-F` accepts parameters like `-F "name=contents"`. If you want the contents to
be read from a file, use `@filename` as contents. When specifying a file, you
can also specify the file content type by appending `;type=<mime type>` to the
file name. You can also post the contents of several files in one field.  For
example, the field name 'coolfiles' is used to send three files, with
different content types using the following syntax:

    curl -F "coolfiles=@fil1.gif;type=image/gif,fil2.txt,fil3.html"
      http://www.post.com/postit.cgi

If the content-type is not specified, curl will try to guess from the file
extension (it only knows a few), or use the previously specified type (from an
earlier file if several files are specified in a list) or else it will use the
default type 'application/octet-stream'.

Emulate a fill-in form with `-F`. Let's say you fill in three fields in a
form. One field is a file name which to post, one field is your name and one
field is a file description. We want to post the file we have written named
"cooltext.txt". To let curl do the posting of this data instead of your
favourite browser, you have to read the HTML source of the form page and find
the names of the input fields. In our example, the input field names are
'file', 'yourname' and 'filedescription'.

    curl -F "file=@cooltext.txt" -F "yourname=Daniel"
      -F "filedescription=Cool text file with cool text inside"
      http://www.post.com/postit.cgi

To send two files in one post you can do it in two ways:

Send multiple files in a single "field" with a single field name:

    curl -F "pictures=@dog.gif,cat.gif" $URL

Send two fields with two field names

    curl -F "docpicture=@dog.gif" -F "catpicture=@cat.gif" $URL

To send a field value literally without interpreting a leading `@` or `<`, or
an embedded `;type=`, use `--form-string` instead of `-F`. This is recommended
when the value is obtained from a user or some other unpredictable
source. Under these circumstances, using `-F` instead of `--form-string` could
allow a user to trick curl into uploading a file.

## Referrer

An HTTP request has the option to include information about which address
referred it to the actual page. curl allows you to specify the referrer to be
used on the command line. It is especially useful to fool or trick stupid
servers or CGI scripts that rely on that information being available or
contain certain data.

    curl -e www.coolsite.com http://www.showme.com/

## User Agent

An HTTP request has the option to include information about the browser that
generated the request. Curl allows it to be specified on the command line. It
is especially useful to fool or trick stupid servers or CGI scripts that only
accept certain browsers.

Example:

    curl -A 'Mozilla/3.0 (Win95; I)' http://www.nationsbank.com/

Other common strings:

- `Mozilla/3.0 (Win95; I)` - Netscape Version 3 for Windows 95
- `Mozilla/3.04 (Win95; U)` - Netscape Version 3 for Windows 95
- `Mozilla/2.02 (OS/2; U)` - Netscape Version 2 for OS/2
- `Mozilla/4.04 [en] (X11; U; AIX 4.2; Nav)` - Netscape for AIX
- `Mozilla/4.05 [en] (X11; U; Linux 2.0.32 i586)` - Netscape for Linux

Note that Internet Explorer tries hard to be compatible in every way:

- `Mozilla/4.0 (compatible; MSIE 4.01; Windows 95)` - MSIE for W95

Mozilla is not the only possible User-Agent name:

- `Konqueror/1.0` - KDE File Manager desktop client
- `Lynx/2.7.1 libwww-FM/2.14` - Lynx command line browser

## Cookies

Cookies are generally used by web servers to keep state information at the
client's side. The server sets cookies by sending a response line in the
headers that looks like `Set-Cookie: <data>` where the data part then
typically contains a set of `NAME=VALUE` pairs (separated by semicolons `;`
like `NAME1=VALUE1; NAME2=VALUE2;`). The server can also specify for what path
the "cookie" should be used for (by specifying `path=value`), when the cookie
should expire (`expire=DATE`), for what domain to use it (`domain=NAME`) and
if it should be used on secure connections only (`secure`).

If you have received a page from a server that contains a header like:

```http
Set-Cookie: sessionid=boo123; path="/foo";
```

it means the server wants that first pair passed on when we get anything in a
path beginning with "/foo".

Example, get a page that wants my name passed in a cookie:

    curl -b "name=Daniel" www.sillypage.com

Curl also has the ability to use previously received cookies in following
sessions. If you get cookies from a server and store them in a file in a
manner similar to:

    curl --dump-header headers www.example.com

... you can then in a second connect to that (or another) site, use the
cookies from the 'headers' file like:

    curl -b headers www.example.com

While saving headers to a file is a working way to store cookies, it is
however error-prone and not the preferred way to do this. Instead, make curl
save the incoming cookies using the well-known netscape cookie format like
this:

    curl -c cookies.txt www.example.com

Note that by specifying `-b` you enable the "cookie awareness" and with `-L`
you can make curl follow a location: (which often is used in combination with
cookies). So that if a site sends cookies and a location, you can use a
non-existing file to trigger the cookie awareness like:

    curl -L -b empty.txt www.example.com

The file to read cookies from must be formatted using plain HTTP headers OR as
netscape's cookie file. Curl will determine what kind it is based on the file
contents. In the above command, curl will parse the header and store the
cookies received from www.example.com. curl will send to the server the
stored cookies which match the request as it follows the location. The file
"empty.txt" may be a nonexistent file.

To read and write cookies from a netscape cookie file, you can set both `-b`
and `-c` to use the same file:

    curl -b cookies.txt -c cookies.txt www.example.com

## Progress Meter

The progress meter exists to show a user that something actually is
happening. The different fields in the output have the following meaning:

    % Total    % Received % Xferd  Average Speed          Time             Curr.
                                   Dload  Upload Total    Current  Left    Speed
    0  151M    0 38608    0     0   9406      0  4:41:43  0:00:04  4:41:39  9287

From left-to-right:

 - %             - percentage completed of the whole transfer
 - Total         - total size of the whole expected transfer
 - %             - percentage completed of the download
 - Received      - currently downloaded amount of bytes
 - %             - percentage completed of the upload
 - Xferd         - currently uploaded amount of bytes
 - Average Speed Dload - the average transfer speed of the download
 - Average Speed Upload - the average transfer speed of the upload
 - Time Total    - expected time to complete the operation
 - Time Current  - time passed since the invoke
 - Time Left     - expected time left to completion
 - Curr.Speed    - the average transfer speed the last 5 seconds (the first
                   5 seconds of a transfer is based on less time of course.)

The `-#` option will display a totally different progress bar that does not
need much explanation!

## Speed Limit

Curl allows the user to set the transfer speed conditions that must be met to
let the transfer keep going. By using the switch `-y` and `-Y` you can make
curl abort transfers if the transfer speed is below the specified lowest limit
for a specified time.

To have curl abort the download if the speed is slower than 3000 bytes per
second for 1 minute, run:

    curl -Y 3000 -y 60 www.far-away-site.com

This can be used in combination with the overall time limit, so that the above
operation must be completed in whole within 30 minutes:

    curl -m 1800 -Y 3000 -y 60 www.far-away-site.com

Forcing curl not to transfer data faster than a given rate is also possible,
which might be useful if you are using a limited bandwidth connection and you
do not want your transfer to use all of it (sometimes referred to as
"bandwidth throttle").

Make curl transfer data no faster than 10 kilobytes per second:

    curl --limit-rate 10K www.far-away-site.com

or

    curl --limit-rate 10240 www.far-away-site.com

Or prevent curl from uploading data faster than 1 megabyte per second:

    curl -T upload --limit-rate 1M ftp://uploadshereplease.com

When using the `--limit-rate` option, the transfer rate is regulated on a
per-second basis, which will cause the total transfer speed to become lower
than the given number. Sometimes of course substantially lower, if your
transfer stalls during periods.

## Config File

Curl automatically tries to read the `.curlrc` file (or `_curlrc` file on
Microsoft Windows systems) from the user's home dir on startup.

The config file could be made up with normal command line switches, but you
can also specify the long options without the dashes to make it more
readable. You can separate the options and the parameter with spaces, or with
`=` or `:`. Comments can be used within the file. If the first letter on a
line is a `#`-symbol the rest of the line is treated as a comment.

If you want the parameter to contain spaces, you must enclose the entire
parameter within double quotes (`"`). Within those quotes, you specify a quote
as `\"`.

NOTE: You must specify options and their arguments on the same line.

Example, set default time out and proxy in a config file:

    # We want a 30 minute timeout:
    -m 1800
    # ... and we use a proxy for all accesses:
    proxy = proxy.our.domain.com:8080

Whitespaces ARE significant at the end of lines, but all whitespace leading
up to the first characters of each line are ignored.

Prevent curl from reading the default file by using -q as the first command
line parameter, like:

    curl -q www.thatsite.com

Force curl to get and display a local help page in case it is invoked without
URL by making a config file similar to:

    # default url to get
    url = "http://help.with.curl.com/curlhelp.html"

You can specify another config file to be read by using the `-K`/`--config`
flag. If you set config file name to `-` it will read the config from stdin,
which can be handy if you want to hide options from being visible in process
tables etc:

    echo "user = user:passwd" | curl -K - http://that.secret.site.com

## Extra Headers

When using curl in your own programs, you may end up needing to pass on your
own custom headers when getting a web page. You can do this by using the `-H`
flag.

Example, send the header `X-you-and-me: yes` to the server when getting a
page:

    curl -H "X-you-and-me: yes" www.love.com

This can also be useful in case you want curl to send a different text in a
header than it normally does. The `-H` header you specify then replaces the
header curl would normally send. If you replace an internal header with an
empty one, you prevent that header from being sent. To prevent the `Host:`
header from being used:

    curl -H "Host:" www.server.com

## FTP and Path Names

Do note that when getting files with a `ftp://` URL, the given path is
relative the directory you enter. To get the file `README` from your home
directory at your ftp site, do:

    curl ftp://user:passwd@my.site.com/README

If you want the README file from the root directory of that same site, you
need to specify the absolute file name:

    curl ftp://user:passwd@my.site.com//README

(I.e with an extra slash in front of the file name.)

## SFTP and SCP and Path Names

With sftp: and scp: URLs, the path name given is the absolute name on the
server. To access a file relative to the remote user's home directory, prefix
the file with `/~/` , such as:

    curl -u $USER sftp://home.example.com/~/.bashrc

## FTP and Firewalls

The FTP protocol requires one of the involved parties to open a second
connection as soon as data is about to get transferred. There are two ways to
do this.

The default way for curl is to issue the PASV command which causes the server
to open another port and await another connection performed by the
client. This is good if the client is behind a firewall that does not allow
incoming connections.

    curl ftp.download.com

If the server, for example, is behind a firewall that does not allow
connections on ports other than 21 (or if it just does not support the `PASV`
command), the other way to do it is to use the `PORT` command and instruct the
server to connect to the client on the given IP number and port (as parameters
to the PORT command).

The `-P` flag to curl supports a few different options. Your machine may have
several IP-addresses and/or network interfaces and curl allows you to select
which of them to use. Default address can also be used:

    curl -P - ftp.download.com

Download with `PORT` but use the IP address of our `le0` interface (this does
not work on windows):

    curl -P le0 ftp.download.com

Download with `PORT` but use 192.168.0.10 as our IP address to use:

    curl -P 192.168.0.10 ftp.download.com

## Network Interface

Get a web page from a server using a specified port for the interface:

    curl --interface eth0:1 http://www.example.com/

or

    curl --interface 192.168.1.10 http://www.example.com/

## HTTPS

Secure HTTP requires a TLS library to be installed and used when curl is
built. If that is done, curl is capable of retrieving and posting documents
using the HTTPS protocol.

Example:

    curl https://www.secure-site.com

curl is also capable of using client certificates to get/post files from sites
that require valid certificates. The only drawback is that the certificate
needs to be in PEM-format. PEM is a standard and open format to store
certificates with, but it is not used by the most commonly used browsers. If
you want curl to use the certificates you use with your (favourite) browser,
you may need to download/compile a converter that can convert your browser's
formatted certificates to PEM formatted ones.

Example on how to automatically retrieve a document using a certificate with a
personal password:

    curl -E /path/to/cert.pem:password https://secure.site.com/

If you neglect to specify the password on the command line, you will be
prompted for the correct password before any data can be received.

Many older HTTPS servers have problems with specific SSL or TLS versions,
which newer versions of OpenSSL etc use, therefore it is sometimes useful to
specify what SSL-version curl should use. Use -3, -2 or -1 to specify that
exact SSL version to use (for SSLv3, SSLv2 or TLSv1 respectively):

    curl -2 https://secure.site.com/

Otherwise, curl will attempt to use a sensible TLS default version.

## Resuming File Transfers

To continue a file transfer where it was previously aborted, curl supports
resume on HTTP(S) downloads as well as FTP uploads and downloads.

Continue downloading a document:

    curl -C - -o file ftp://ftp.server.com/path/file

Continue uploading a document:

    curl -C - -T file ftp://ftp.server.com/path/file

Continue downloading a document from a web server

    curl -C - -o file http://www.server.com/

## Time Conditions

HTTP allows a client to specify a time condition for the document it requests.
It is `If-Modified-Since` or `If-Unmodified-Since`. curl allows you to specify
them with the `-z`/`--time-cond` flag.

For example, you can easily make a download that only gets performed if the
remote file is newer than a local copy. It would be made like:

    curl -z local.html http://remote.server.com/remote.html

Or you can download a file only if the local file is newer than the remote
one. Do this by prepending the date string with a `-`, as in:

    curl -z -local.html http://remote.server.com/remote.html

You can specify a "free text" date as condition. Tell curl to only download
the file if it was updated since January 12, 2012:

    curl -z "Jan 12 2012" http://remote.server.com/remote.html

Curl will then accept a wide range of date formats. You always make the date
check the other way around by prepending it with a dash (`-`).

## DICT

For fun try

    curl dict://dict.org/m:curl
    curl dict://dict.org/d:heisenbug:jargon
    curl dict://dict.org/d:daniel:gcide

Aliases for 'm' are 'match' and 'find', and aliases for 'd' are 'define' and
'lookup'. For example,

    curl dict://dict.org/find:curl

Commands that break the URL description of the RFC (but not the DICT
protocol) are

    curl dict://dict.org/show:db
    curl dict://dict.org/show:strat

Authentication support is still missing

## LDAP

If you have installed the OpenLDAP library, curl can take advantage of it and
offer `ldap://` support. On Windows, curl will use WinLDAP from Platform SDK
by default.

Default protocol version used by curl is LDAPv3. LDAPv2 will be used as
fallback mechanism in case if LDAPv3 will fail to connect.

LDAP is a complex thing and writing an LDAP query is not an easy task. I do
advise you to dig up the syntax description for that elsewhere. One such place
might be: [RFC 2255, The LDAP URL
Format](https://curl.se/rfc/rfc2255.txt)

To show you an example, this is how I can get all people from my local LDAP
server that has a certain sub-domain in their email address:

    curl -B "ldap://ldap.frontec.se/o=frontec??sub?mail=*sth.frontec.se"

If I want the same info in HTML format, I can get it by not using the `-B`
(enforce ASCII) flag.

You also can use authentication when accessing LDAP catalog:

    curl -u user:passwd "ldap://ldap.frontec.se/o=frontec??sub?mail=*"
    curl "ldap://user:passwd@ldap.frontec.se/o=frontec??sub?mail=*"

By default, if user and password provided, OpenLDAP/WinLDAP will use basic
authentication. On Windows you can control this behavior by providing one of
`--basic`, `--ntlm` or `--digest` option in curl command line

    curl --ntlm "ldap://user:passwd@ldap.frontec.se/o=frontec??sub?mail=*"

On Windows, if no user/password specified, auto-negotiation mechanism will be
used with current logon credentials (SSPI/SPNEGO).

## Environment Variables

Curl reads and understands the following environment variables:

    http_proxy, HTTPS_PROXY, FTP_PROXY

They should be set for protocol-specific proxies. General proxy should be set
with

    ALL_PROXY

A comma-separated list of host names that should not go through any proxy is
set in (only an asterisk, `*` matches all hosts)

    NO_PROXY

If the host name matches one of these strings, or the host is within the
domain of one of these strings, transactions with that node will not be
proxied. When a domain is used, it needs to start with a period. A user can
specify that both www.example.com and foo.example.com should not use a proxy
by setting `NO_PROXY` to `.example.com`. By including the full name you can
exclude specific host names, so to make `www.example.com` not use a proxy but
still have `foo.example.com` do it, set `NO_PROXY` to `www.example.com`.

The usage of the `-x`/`--proxy` flag overrides the environment variables.

## Netrc

Unix introduced the `.netrc` concept a long time ago. It is a way for a user
to specify name and password for commonly visited FTP sites in a file so that
you do not have to type them in each time you visit those sites. You realize
this is a big security risk if someone else gets hold of your passwords, so
therefore most unix programs will not read this file unless it is only readable
by yourself (curl does not care though).

Curl supports `.netrc` files if told to (using the `-n`/`--netrc` and
`--netrc-optional` options). This is not restricted to just FTP, so curl can
use it for all protocols where authentication is used.

A simple `.netrc` file could look something like:

    machine curl.se login iamdaniel password mysecret

## Custom Output

To better allow script programmers to get to know about the progress of curl,
the `-w`/`--write-out` option was introduced. Using this, you can specify what
information from the previous transfer you want to extract.

To display the amount of bytes downloaded together with some text and an
ending newline:

    curl -w 'We downloaded %{size_download} bytes\n' www.download.com

## Kerberos FTP Transfer

Curl supports kerberos4 and kerberos5/GSSAPI for FTP transfers. You need the
kerberos package installed and used at curl build time for it to be available.

First, get the krb-ticket the normal way, like with the kinit/kauth tool.
Then use curl in way similar to:

    curl --krb private ftp://krb4site.com -u username:fakepwd

There's no use for a password on the `-u` switch, but a blank one will make
curl ask for one and you already entered the real password to kinit/kauth.

## TELNET

The curl telnet support is basic and easy to use. Curl passes all data passed
to it on stdin to the remote server. Connect to a remote telnet server using a
command line similar to:

    curl telnet://remote.server.com

And enter the data to pass to the server on stdin. The result will be sent to
stdout or to the file you specify with `-o`.

You might want the `-N`/`--no-buffer` option to switch off the buffered output
for slow connections or similar.

Pass options to the telnet protocol negotiation, by using the `-t` option. To
tell the server we use a vt100 terminal, try something like:

    curl -tTTYPE=vt100 telnet://remote.server.com

Other interesting options for it `-t` include:

 - `XDISPLOC=<X display>` Sets the X display location.
 - `NEW_ENV=<var,val>` Sets an environment variable.

NOTE: The telnet protocol does not specify any way to login with a specified
user and password so curl cannot do that automatically. To do that, you need to
track when the login prompt is received and send the username and password
accordingly.

## Persistent Connections

Specifying multiple files on a single command line will make curl transfer all
of them, one after the other in the specified order.

libcurl will attempt to use persistent connections for the transfers so that
the second transfer to the same host can use the same connection that was
already initiated and was left open in the previous transfer. This greatly
decreases connection time for all but the first transfer and it makes a far
better use of the network.

Note that curl cannot use persistent connections for transfers that are used
in subsequence curl invokes. Try to stuff as many URLs as possible on the same
command line if they are using the same host, as that will make the transfers
faster. If you use an HTTP proxy for file transfers, practically all transfers
will be persistent.

## Multiple Transfers With A Single Command Line

As is mentioned above, you can download multiple files with one command line
by simply adding more URLs. If you want those to get saved to a local file
instead of just printed to stdout, you need to add one save option for each
URL you specify. Note that this also goes for the `-O` option (but not
`--remote-name-all`).

For example: get two files and use `-O` for the first and a custom file
name for the second:

    curl -O http://url.com/file.txt ftp://ftp.com/moo.exe -o moo.jpg

You can also upload multiple files in a similar fashion:

    curl -T local1 ftp://ftp.com/moo.exe -T local2 ftp://ftp.com/moo2.txt

## IPv6

curl will connect to a server with IPv6 when a host lookup returns an IPv6
address and fall back to IPv4 if the connection fails. The `--ipv4` and
`--ipv6` options can specify which address to use when both are
available. IPv6 addresses can also be specified directly in URLs using the
syntax:

    http://[2001:1890:1112:1::20]/overview.html

When this style is used, the `-g` option must be given to stop curl from
interpreting the square brackets as special globbing characters. Link local
and site local addresses including a scope identifier, such as `fe80::1234%1`,
may also be used, but the scope portion must be numeric or match an existing
network interface on Linux and the percent character must be URL escaped. The
previous example in an SFTP URL might look like:

    sftp://[fe80::1234%251]/

IPv6 addresses provided other than in URLs (e.g. to the `--proxy`,
`--interface` or `--ftp-port` options) should not be URL encoded.

## Mailing Lists

For your convenience, we have several open mailing lists to discuss curl, its
development and things relevant to this. Get all info at
https://curl.se/mail/.

Please direct curl questions, feature requests and trouble reports to one of
these mailing lists instead of mailing any individual.

Available lists include:

### curl-users

Users of the command line tool. How to use it, what does not work, new
features, related tools, questions, news, installations, compilations,
running, porting etc.

### curl-library

Developers using or developing libcurl. Bugs, extensions, improvements.

### curl-announce

Low-traffic. Only receives announcements of new public versions. At worst,
that makes something like one or two mails per month, but usually only one
mail every second month.

### curl-and-php

Using the curl functions in PHP. Everything curl with a PHP angle. Or PHP with
a curl angle.

### curl-and-python

Python hackers using curl with or without the python binding pycurl.

curl release procedure - how to do a release
============================================

in the source code repo
-----------------------

- run `./scripts/copyright.pl` and correct possible omissions

- edit `RELEASE-NOTES` to be accurate

- update `docs/THANKS`

- make sure all relevant changes are committed on the master branch

- tag the git repo in this style: `git tag -a curl-7_34_0`. -a annotates the
  tag and we use underscores instead of dots in the version number. Make sure
  the tag is GPG signed (using -s).

- run "./maketgz 7.34.0" to build the release tarballs. It is important that
  you run this on a machine with the correct set of autotools etc installed
  as this is what then will be shipped and used by most users on \*nix like
  systems.

- push the git commits and the new tag

- gpg sign the 4 tarballs as maketgz suggests

- upload the 8 resulting files to the primary download directory

in the curl-www repo
--------------------

- edit `Makefile` (version number and date),

- edit `_newslog.html` (announce the new release) and

- edit `_changes.html` (insert changes+bugfixes from RELEASE-NOTES)

- commit all local changes

- tag the repo with the same name as used for the source repo.

- make sure all relevant changes are committed and pushed on the master branch

  (the website then updates its contents automatically)

on GitHub
---------

- edit the newly made release tag so that it is listed as the latest release

inform
------

- send an email to curl-users, curl-announce and curl-library. Insert the
  RELEASE-NOTES into the mail.

celebrate
---------

- suitable beverage intake is encouraged for the festivities

curl release scheduling
=======================

Release Cycle
-------------

We do releases every 8 weeks on Wednesdays. If critical problems arise, we can
insert releases outside of the schedule or we can move the release date - but
this is rare.

Each 8 week release cycle is split in two 4-week periods.

- During the first 4 weeks after a release, we allow new features and changes
  to curl and libcurl. If we accept any such changes, we bump the minor number
  used for the next release.

- During the second 4-week period we do not merge any features or changes, we
  then only focus on fixing bugs and polishing things to make a solid coming
  release.

- After a regular procedure-following release (made on Wednesdays), the
  feature window remains closed until the following Monday in case of special
  actions or patch releases etc.

If a future release date happens to end up on a "bad date", like in the middle
of common public holidays or when the lead release manager is away traveling,
the release date can be moved forwards or backwards a full week. This is then
advertised well in advance.

Coming dates
------------

Based on the description above, here are some planned release dates (at the
time of this writing):

- January 5, 2022 (7.81.0)
- March 2, 2022
- April 27, 2022
- June 22, 2022
- August 17, 2022
- October 12, 2022
- December 7, 2022
- February 1, 2023
- March 20, 2023 (8.0.0)
# curl C code style

Source code that has a common style is easier to read than code that uses
different styles in different places. It helps making the code feel like one
single code base. Easy-to-read is an important property of code and helps
making it easier to review when new things are added and it helps debugging
code when developers are trying to figure out why things go wrong. A unified
style is more important than individual contributors having their own personal
tastes satisfied.

Our C code has a few style rules. Most of them are verified and upheld by the
`lib/checksrc.pl` script. Invoked with `make checksrc` or even by default by
the build system when built after `./configure --enable-debug` has been used.

It is normally not a problem for anyone to follow the guidelines, as you just
need to copy the style already used in the source code and there are no
particularly unusual rules in our set of rules.

We also work hard on writing code that are warning-free on all the major
platforms and in general on as many platforms as possible. Code that obviously
will cause warnings will not be accepted as-is.

## Naming

Try using a non-confusing naming scheme for your new functions and variable
names. It does not necessarily have to mean that you should use the same as in
other places of the code, just that the names should be logical,
understandable and be named according to what they are used for. File-local
functions should be made static. We like lower case names.

See the [INTERNALS](https://curl.se/dev/internals.html#symbols) document on
how we name non-exported library-global symbols.

## Indenting

We use only spaces for indentation, never TABs. We use two spaces for each new
open brace.

```c
if(something_is_true) {
  while(second_statement == fine) {
    moo();
  }
}
```

## Comments

Since we write C89 code, **//** comments are not allowed. They were not
introduced in the C standard until C99. We use only __/* comments */__.

```c
/* this is a comment */
```

## Long lines

Source code in curl may never be wider than 79 columns and there are two
reasons for maintaining this even in the modern era of large and high
resolution screens:

1. Narrower columns are easier to read than wide ones. There's a reason
   newspapers have used columns for decades or centuries.

2. Narrower columns allow developers to easier show multiple pieces of code
   next to each other in different windows. I often have two or three source
   code windows next to each other on the same screen - as well as multiple
   terminal and debugging windows.

## Braces

In if/while/do/for expressions, we write the open brace on the same line as
the keyword and we then set the closing brace on the same indentation level as
the initial keyword. Like this:

```c
if(age < 40) {
  /* clearly a youngster */
}
```

You may omit the braces if they would contain only a one-line statement:

```c
if(!x)
  continue;
```

For functions the opening brace should be on a separate line:

```c
int main(int argc, char **argv)
{
  return 1;
}
```

## 'else' on the following line

When adding an **else** clause to a conditional expression using braces, we
add it on a new line after the closing brace. Like this:

```c
if(age < 40) {
  /* clearly a youngster */
}
else {
  /* probably grumpy */
}
```

## No space before parentheses

When writing expressions using if/while/do/for, there shall be no space
between the keyword and the open parenthesis. Like this:

```c
while(1) {
  /* loop forever */
}
```

## Use boolean conditions

Rather than test a conditional value such as a bool against TRUE or FALSE, a
pointer against NULL or != NULL and an int against zero or not zero in
if/while conditions we prefer:

```c
result = do_something();
if(!result) {
  /* something went wrong */
  return result;
}
```

## No assignments in conditions

To increase readability and reduce complexity of conditionals, we avoid
assigning variables within if/while conditions. We frown upon this style:

```c
if((ptr = malloc(100)) == NULL)
  return NULL;
```

and instead we encourage the above version to be spelled out more clearly:

```c
ptr = malloc(100);
if(!ptr)
  return NULL;
```

## New block on a new line

We never write multiple statements on the same source line, even for short
if() conditions.

```c
if(a)
  return TRUE;
else if(b)
  return FALSE;
```

and NEVER:

```c
if(a) return TRUE;
else if(b) return FALSE;
```

## Space around operators

Please use spaces on both sides of operators in C expressions. Postfix **(),
[], ->, ., ++, --** and Unary **+, -, !, ~, &** operators excluded they should
have no space.

Examples:

```c
bla = func();
who = name[0];
age += 1;
true = !false;
size += -2 + 3 * (a + b);
ptr->member = a++;
struct.field = b--;
ptr = &address;
contents = *pointer;
complement = ~bits;
empty = (!*string) ? TRUE : FALSE;
```

## No parentheses for return values

We use the 'return' statement without extra parentheses around the value:

```c
int works(void)
{
  return TRUE;
}
```

## Parentheses for sizeof arguments

When using the sizeof operator in code, we prefer it to be written with
parentheses around its argument:

```c
int size = sizeof(int);
```

## Column alignment

Some statements cannot be completed on a single line because the line would be
too long, the statement too hard to read, or due to other style guidelines
above. In such a case the statement will span multiple lines.

If a continuation line is part of an expression or sub-expression then you
should align on the appropriate column so that it's easy to tell what part of
the statement it is. Operators should not start continuation lines. In other
cases follow the 2-space indent guideline. Here are some examples from
libcurl:

```c
if(Curl_pipeline_wanted(handle->multi, CURLPIPE_HTTP1) &&
   (handle->set.httpversion != CURL_HTTP_VERSION_1_0) &&
   (handle->set.httpreq == HTTPREQ_GET ||
    handle->set.httpreq == HTTPREQ_HEAD))
  /* did not ask for HTTP/1.0 and a GET or HEAD */
  return TRUE;
```

If no parenthesis, use the default indent:

```c
data->set.http_disable_hostname_check_before_authentication =
  (0 != va_arg(param, long)) ? TRUE : FALSE;
```

Function invoke with an open parenthesis:

```c
if(option) {
  result = parse_login_details(option, strlen(option),
                               (userp ? &user : NULL),
                               (passwdp ? &passwd : NULL),
                               NULL);
}
```

Align with the "current open" parenthesis:

```c
DEBUGF(infof(data, "Curl_pp_readresp_ %d bytes of trailing "
             "server response left\n",
             (int)clipamount));
```

## Platform dependent code

Use **#ifdef HAVE_FEATURE** to do conditional code. We avoid checking for
particular operating systems or hardware in the #ifdef lines. The HAVE_FEATURE
shall be generated by the configure script for unix-like systems and they are
hard-coded in the `config-[system].h` files for the others.

We also encourage use of macros/functions that possibly are empty or defined
to constants when libcurl is built without that feature, to make the code
seamless. Like this example where the **magic()** function works differently
depending on a build-time conditional:

```c
#ifdef HAVE_MAGIC
void magic(int a)
{
  return a + 2;
}
#else
#define magic(x) 1
#endif

int content = magic(3);
```

## No typedefed structs

Use structs by all means, but do not typedef them. Use the `struct name` way
of identifying them:

```c
struct something {
   void *valid;
   size_t way_to_write;
};
struct something instance;
```

**Not okay**:

```c
typedef struct {
   void *wrong;
   size_t way_to_write;
} something;
something instance;
```
# Experimental

Some features and functionality in curl and libcurl are considered
**EXPERIMENTAL**.

Experimental support in curl means:

1. Experimental features are provided to allow users to try them out and
   provide feedback on functionality and API etc before they ship and get
   "carved in stone".
2. You must enable the feature when invoking configure as otherwise curl will
   not be built with the feature present.
3. We strongly advice against using this feature in production.
4. **We reserve the right to change behavior** of the feature without sticking
   to our API/ABI rules as we do for regular features, as long as it is marked
   experimental.
5. Experimental features are clearly marked so in documentation. Beware.

## Experimental features right now

 - The Hyper HTTP backend
 - HTTP/3 support and options
 - CURLSSLOPT_NATIVE_CA (No configure option, feature built in when supported)
# bufref

This is an internal module for handling buffer references. A referenced
buffer is associated with its destructor function that is implicitly called
when the reference is invalidated. Once referenced, a buffer cannot be
reallocated.

A data length is stored within the reference for binary data handling
purpose; it is not used by the bufref API.

The `struct bufref` is used to hold data referencing a buffer. The members of
that structure **MUST NOT** be accessed or modified without using the dedicated
bufref API.

## init

```c
void Curl_bufref_init(struct bufref *br);
```

Initialises a `bufref` structure. This function **MUST** be called before any
other operation is performed on the structure.

Upon completion, the referenced buffer is `NULL` and length is zero.

This function may also be called to bypass referenced buffer destruction while
invalidating the current reference.

## free

```c
void Curl_bufref_free(struct bufref *br);
```

Destroys the previously referenced buffer using its destructor and
reinitialises the structure for a possible subsequent reuse.

## set

```c
void Curl_bufref_set(struct bufref *br, const void *buffer, size_t length,
                     void (*destructor)(void *));
```

Releases the previously referenced buffer, then assigns the new `buffer` to
the structure, associated with its `destructor` function. The later can be
specified as `NULL`: this will be the case when the referenced buffer is
static.

if `buffer` is NULL, `length`must be zero.

## memdup

```c
CURLcode Curl_bufref_memdup(struct bufref *br, const void *data, size_t length);
```

Releases the previously referenced buffer, then duplicates the `length`-byte
`data` into a buffer allocated via `malloc()` and references the latter
associated with destructor `curl_free()`.

An additional trailing byte is allocated and set to zero as a possible
string zero-terminator; it is not counted in the stored length.

Returns `CURLE_OK` if successful, else `CURLE_OUT_OF_MEMORY`.

## ptr

```c
const unsigned char *Curl_bufref_ptr(const struct bufref *br);
```

Returns a `const unsigned char *` to the referenced buffer.

## len

```c
size_t Curl_bufref_len(const struct bufref *br);
```

Returns the stored length of the referenced buffer.
# HSTS support

HTTP Strict-Transport-Security. Added as experimental in curl
7.74.0. Supported "for real" since 7.77.0.

## Standard

[HTTP Strict Transport Security](https://datatracker.ietf.org/doc/html/rfc6797)

## Behavior

libcurl features an in-memory cache for HSTS hosts, so that subsequent
HTTP-only requests to a host name present in the cache will get internally
"redirected" to the HTTPS version.

## `curl_easy_setopt()` options:

 - `CURLOPT_HSTS_CTRL` - enable HSTS for this easy handle
 - `CURLOPT_HSTS` - specify file name where to store the HSTS cache on close
  (and possibly read from at startup)

## curl cmdline options

 - `--hsts [filename]` - enable HSTS, use the file as HSTS cache. If filename
   is `""` (no length) then no file will be used, only in-memory cache.

## HSTS cache file format

Lines starting with `#` are ignored.

For each hsts entry:

    [host name] "YYYYMMDD HH:MM:SS"

The `[host name]` is dot-prefixed if it is a includeSubDomain.

The time stamp is when the entry expires.

I considered using wget's file format for the HSTS cache. However, they store the time stamp as the epoch (number of seconds since 1970) and I strongly disagree with using that format. Instead I opted to use a format similar to the curl alt-svc cache file format.

## Possible future additions

 - `CURLOPT_HSTS_PRELOAD` - provide a set of preloaded HSTS host names
 - ability to save to something else than a file
# how to install curl and libcurl

## Installing Binary Packages

Lots of people download binary distributions of curl and libcurl. This
document does not describe how to install curl or libcurl using such a binary
package. This document describes how to compile, build and install curl and
libcurl from source code.

## Building using vcpkg

You can download and install curl and libcurl using the [vcpkg](https://github.com/Microsoft/vcpkg/) dependency manager:

    git clone https://github.com/Microsoft/vcpkg.git
    cd vcpkg
    ./bootstrap-vcpkg.sh
    ./vcpkg integrate install
    vcpkg install curl[tool]

The curl port in vcpkg is kept up to date by Microsoft team members and community contributors. If the version is out of date, please [create an issue or pull request](https://github.com/Microsoft/vcpkg) on the vcpkg repository.

## Building from git

If you get your code off a git repository instead of a release tarball, see
the `GIT-INFO` file in the root directory for specific instructions on how to
proceed.

# Unix

A normal Unix installation is made in three or four steps (after you have
unpacked the source archive):

    ./configure --with-openssl [--with-gnutls --with-wolfssl]
    make
    make test (optional)
    make install

(Adjust the configure line accordingly to use the TLS library you want.)

You probably need to be root when doing the last command.

Get a full listing of all available configure options by invoking it like:

    ./configure --help

If you want to install curl in a different file hierarchy than `/usr/local`,
specify that when running configure:

    ./configure --prefix=/path/to/curl/tree

If you have write permission in that directory, you can do 'make install'
without being root. An example of this would be to make a local install in
your own home directory:

    ./configure --prefix=$HOME
    make
    make install

The configure script always tries to find a working SSL library unless
explicitly told not to. If you have OpenSSL installed in the default search
path for your compiler/linker, you do not need to do anything special. If you
have OpenSSL installed in `/usr/local/ssl`, you can run configure like:

    ./configure --with-openssl

If you have OpenSSL installed somewhere else (for example, `/opt/OpenSSL`) and
you have pkg-config installed, set the pkg-config path first, like this:

    env PKG_CONFIG_PATH=/opt/OpenSSL/lib/pkgconfig ./configure --with-openssl

Without pkg-config installed, use this:

    ./configure --with-openssl=/opt/OpenSSL

If you insist on forcing a build without SSL support, even though you may
have OpenSSL installed in your system, you can run configure like this:

    ./configure --without-ssl

If you have OpenSSL installed, but with the libraries in one place and the
header files somewhere else, you have to set the `LDFLAGS` and `CPPFLAGS`
environment variables prior to running configure. Something like this should
work:

    CPPFLAGS="-I/path/to/ssl/include" LDFLAGS="-L/path/to/ssl/lib" ./configure

If you have shared SSL libs installed in a directory where your run-time
linker does not find them (which usually causes configure failures), you can
provide this option to gcc to set a hard-coded path to the run-time linker:

    LDFLAGS=-Wl,-R/usr/local/ssl/lib ./configure --with-openssl

## More Options

To force a static library compile, disable the shared library creation by
running configure like:

    ./configure --disable-shared

To tell the configure script to skip searching for thread-safe functions, add
an option like:

    ./configure --disable-thread

If you are a curl developer and use gcc, you might want to enable more debug
options with the `--enable-debug` option.

curl can be built to use a whole range of libraries to provide various useful
services, and configure will try to auto-detect a decent default. But if you
want to alter it, you can select how to deal with each individual library.

## Select TLS backend

These options are provided to select TLS backend to use.

 - AmiSSL: `--with-amissl`
 - BearSSL: `--with-bearssl`
 - GnuTLS: `--with-gnutls`.
 - mbedTLS: `--with-mbedtls`
 - NSS: `--with-nss`
 - OpenSSL: `--with-openssl` (also for BoringSSL and libressl)
 - rustls: `--with-rustls`
 - schannel: `--with-schannel`
 - secure transport: `--with-secure-transport`
 - wolfSSL: `--with-wolfssl`

# Windows

## Building Windows DLLs and C run-time (CRT) linkage issues

 As a general rule, building a DLL with static CRT linkage is highly
 discouraged, and intermixing CRTs in the same app is something to avoid at
 any cost.

 Reading and comprehending Microsoft Knowledge Base articles KB94248 and
 KB140584 is a must for any Windows developer. Especially important is full
 understanding if you are not going to follow the advice given above.

 - [How To Use the C Run-Time](https://support.microsoft.com/help/94248/how-to-use-the-c-run-time)
 - [Run-Time Library Compiler Options](https://docs.microsoft.com/cpp/build/reference/md-mt-ld-use-run-time-library)
 - [Potential Errors Passing CRT Objects Across DLL Boundaries](https://docs.microsoft.com/cpp/c-runtime-library/potential-errors-passing-crt-objects-across-dll-boundaries)

If your app is misbehaving in some strange way, or it is suffering from
memory corruption, before asking for further help, please try first to
rebuild every single library your app uses as well as your app using the
debug multithreaded dynamic C runtime.

 If you get linkage errors read section 5.7 of the FAQ document.

## MingW32

Make sure that MinGW32's bin dir is in the search path, for example:

```cmd
set PATH=c:\mingw32\bin;%PATH%
```

then run `mingw32-make mingw32` in the root dir. There are other
make targets available to build libcurl with more features, use:

 - `mingw32-make mingw32-zlib` to build with Zlib support;
 - `mingw32-make mingw32-ssl-zlib` to build with SSL and Zlib enabled;
 - `mingw32-make mingw32-ssh2-ssl-zlib` to build with SSH2, SSL, Zlib;
 - `mingw32-make mingw32-ssh2-ssl-sspi-zlib` to build with SSH2, SSL, Zlib
   and SSPI support.

If you have any problems linking libraries or finding header files, be sure
to verify that the provided `Makefile.m32` files use the proper paths, and
adjust as necessary. It is also possible to override these paths with
environment variables, for example:

```cmd
set ZLIB_PATH=c:\zlib-1.2.8
set OPENSSL_PATH=c:\openssl-1.0.2c
set LIBSSH2_PATH=c:\libssh2-1.6.0
```

It is also possible to build with other LDAP SDKs than MS LDAP; currently
it is possible to build with native Win32 OpenLDAP, or with the Novell CLDAP
SDK. If you want to use these you need to set these vars:

```cmd
set LDAP_SDK=c:\openldap
set USE_LDAP_OPENLDAP=1
```

or for using the Novell SDK:

```cmd
set USE_LDAP_NOVELL=1
```

If you want to enable LDAPS support then set LDAPS=1.

## Cygwin

Almost identical to the unix installation. Run the configure script in the
curl source tree root with `sh configure`. Make sure you have the `sh`
executable in `/bin/` or you will see the configure fail toward the end.

Run `make`

## Disabling Specific Protocols in Windows builds

The configure utility, unfortunately, is not available for the Windows
environment, therefore, you cannot use the various disable-protocol options of
the configure utility on this platform.

You can use specific defines to disable specific protocols and features. See
[CURL-DISABLE.md](CURL-DISABLE.md) for the full list.

If you want to set any of these defines you have the following options:

 - Modify `lib/config-win32.h`
 - Modify `lib/curl_setup.h`
 - Modify `winbuild/Makefile.vc`
 - Modify the "Preprocessor Definitions" in the libcurl project

Note: The pre-processor settings can be found using the Visual Studio IDE
under "Project -> Settings -> C/C++ -> General" in VC6 and "Project ->
Properties -> Configuration Properties -> C/C++ -> Preprocessor" in later
versions.

## Using BSD-style lwIP instead of Winsock TCP/IP stack in Win32 builds

In order to compile libcurl and curl using BSD-style lwIP TCP/IP stack it is
necessary to make definition of preprocessor symbol `USE_LWIPSOCK` visible to
libcurl and curl compilation processes. To set this definition you have the
following alternatives:

 - Modify `lib/config-win32.h` and `src/config-win32.h`
 - Modify `winbuild/Makefile.vc`
 - Modify the "Preprocessor Definitions" in the libcurl project

Note: The pre-processor settings can be found using the Visual Studio IDE
under "Project -> Settings -> C/C++ -> General" in VC6 and "Project ->
Properties -> Configuration Properties -> C/C++ -> Preprocessor" in later
versions.

Once that libcurl has been built with BSD-style lwIP TCP/IP stack support, in
order to use it with your program it is mandatory that your program includes
lwIP header file `<lwip/opt.h>` (or another lwIP header that includes this)
before including any libcurl header. Your program does not need the
`USE_LWIPSOCK` preprocessor definition which is for libcurl internals only.

Compilation has been verified with [lwIP
1.4.0](https://download.savannah.gnu.org/releases/lwip/lwip-1.4.0.zip) and
[contrib-1.4.0](https://download.savannah.gnu.org/releases/lwip/contrib-1.4.0.zip).

This BSD-style lwIP TCP/IP stack support must be considered experimental given
that it has been verified that lwIP 1.4.0 still needs some polish, and libcurl
might yet need some additional adjustment, caveat emptor.

## Important static libcurl usage note

When building an application that uses the static libcurl library on Windows,
you must add `-DCURL_STATICLIB` to your `CFLAGS`. Otherwise the linker will
look for dynamic import symbols.

## Legacy Windows and SSL

Schannel (from Windows SSPI), is the native SSL library in Windows. However,
Schannel in Windows <= XP is unable to connect to servers that
no longer support the legacy handshakes and algorithms used by those
versions. If you will be using curl in one of those earlier versions of
Windows you should choose another SSL backend such as OpenSSL.

# Apple Platforms (macOS, iOS, tvOS, watchOS, and their simulator counterparts)

On modern Apple operating systems, curl can be built to use Apple's SSL/TLS
implementation, Secure Transport, instead of OpenSSL. To build with Secure
Transport for SSL/TLS, use the configure option `--with-secure-transport`. (It
is not necessary to use the option `--without-openssl`.)

When Secure Transport is in use, the curl options `--cacert` and `--capath`
and their libcurl equivalents, will be ignored, because Secure Transport uses
the certificates stored in the Keychain to evaluate whether or not to trust
the server. This, of course, includes the root certificates that ship with the
OS. The `--cert` and `--engine` options, and their libcurl equivalents, are
currently unimplemented in curl with Secure Transport.

In general, a curl build for an Apple `ARCH/SDK/DEPLOYMENT_TARGET` combination
can be taken by providing appropriate values for `ARCH`, `SDK`, `DEPLOYMENT_TARGET`
below and running the commands:

```bash
# Set these three according to your needs
export ARCH=x86_64
export SDK=macosx
export DEPLOYMENT_TARGET=10.8

export CFLAGS="-arch $ARCH -isysroot $(xcrun -sdk $SDK --show-sdk-path) -m$SDK-version-min=$DEPLOYMENT_TARGET"
./configure --host=$ARCH-apple-darwin --prefix $(pwd)/artifacts --with-secure-transport
make -j8
make install
```

Above will build curl for macOS platform with `x86_64` architecture and `10.8` as deployment target.

Here is an example for iOS device:

```bash
export ARCH=arm64
export SDK=iphoneos
export DEPLOYMENT_TARGET=11.0

export CFLAGS="-arch $ARCH -isysroot $(xcrun -sdk $SDK --show-sdk-path) -m$SDK-version-min=$DEPLOYMENT_TARGET"
./configure --host=$ARCH-apple-darwin --prefix $(pwd)/artifacts --with-secure-transport
make -j8
make install
```

Another example for watchOS simulator for macs with Apple Silicon:

```bash
export ARCH=arm64
export SDK=watchsimulator
export DEPLOYMENT_TARGET=5.0

export CFLAGS="-arch $ARCH -isysroot $(xcrun -sdk $SDK --show-sdk-path) -m$SDK-version-min=$DEPLOYMENT_TARGET"
./configure --host=$ARCH-apple-darwin --prefix $(pwd)/artifacts --with-secure-transport
make -j8
make install
```

In all above, the built libraries and executables can be found in `artifacts` folder.

# Android

When building curl for Android it's recommended to use a Linux environment
since using curl's `configure` script is the easiest way to build curl
for Android. Before you can build curl for Android, you need to install the
Android NDK first. This can be done using the SDK Manager that is part of
Android Studio. Once you have installed the Android NDK, you need to figure out
where it has been installed and then set up some environment variables before
launching `configure`. On macOS, those variables could look like this to compile
for `aarch64` and API level 29:

```bash
export NDK=~/Library/Android/sdk/ndk/20.1.5948944
export HOST_TAG=darwin-x86_64
export TOOLCHAIN=$NDK/toolchains/llvm/prebuilt/$HOST_TAG
export AR=$TOOLCHAIN/bin/aarch64-linux-android-ar
export AS=$TOOLCHAIN/bin/aarch64-linux-android-as
export CC=$TOOLCHAIN/bin/aarch64-linux-android29-clang
export CXX=$TOOLCHAIN/bin/aarch64-linux-android29-clang++
export LD=$TOOLCHAIN/bin/aarch64-linux-android-ld
export RANLIB=$TOOLCHAIN/bin/aarch64-linux-android-ranlib
export STRIP=$TOOLCHAIN/bin/aarch64-linux-android-strip
```

When building on Linux or targeting other API levels or architectures, you need
to adjust those variables accordingly. After that you can build curl like this:

    ./configure --host aarch64-linux-android --with-pic --disable-shared

Note that this will not give you SSL/TLS support. If you need SSL/TLS, you have
to build curl against a SSL/TLS layer, e.g. OpenSSL, because it's impossible for
curl to access Android's native SSL/TLS layer. To build curl for Android using
OpenSSL, follow the OpenSSL build instructions and then install `libssl.a` and
`libcrypto.a` to `$TOOLCHAIN/sysroot/usr/lib` and copy `include/openssl` to
`$TOOLCHAIN/sysroot/usr/include`. Now you can build curl for Android using
OpenSSL like this:

    ./configure --host aarch64-linux-android --with-pic --disable-shared --with-openssl="$TOOLCHAIN/sysroot/usr"

Note, however, that you must target at least Android M (API level 23) or `configure`
will not be able to detect OpenSSL since `stderr` (and the like) were not defined
before Android M.

# IBM i

For IBM i (formerly OS/400), you can use curl in two different ways:

- Natively, running in the **ILE**. The obvious use is being able to call curl
  from ILE C or RPG applications.
  - You will need to build this from source. See `packages/OS400/README` for
    the ILE specific build instructions.
- In the **PASE** environment, which runs AIX programs. curl will be built as
  it would be on AIX.
  - IBM provides builds of curl in their Yum repository for PASE software.
  - To build from source, follow the Unix instructions.

There are some additional limitations and quirks with curl on this platform;
they affect both environments.

## Multithreading notes

By default, jobs in IBM i will not start with threading enabled. (Exceptions
include interactive PASE sessions started by `QP2TERM` or SSH.) If you use
curl in an environment without threading when options like async DNS were
enabled, you will get messages like:

```
getaddrinfo() thread failed to start
```

Do not panic. curl and your program are not broken. You can fix this by:

- Set the environment variable `QIBM_MULTI_THREADED` to `Y` before starting
  your program. This can be done at whatever scope you feel is appropriate.
- Alternatively, start the job with the `ALWMLTTHD` parameter set to `*YES`.

# Cross compile

Download and unpack the curl package.

`cd` to the new directory. (e.g. `cd curl-7.12.3`)

Set environment variables to point to the cross-compile toolchain and call
configure with any options you need. Be sure and specify the `--host` and
`--build` parameters at configuration time. The following script is an
example of cross-compiling for the IBM 405GP PowerPC processor using the
toolchain from MonteVista for Hardhat Linux.

```bash
#! /bin/sh

export PATH=$PATH:/opt/hardhat/devkit/ppc/405/bin
export CPPFLAGS="-I/opt/hardhat/devkit/ppc/405/target/usr/include"
export AR=ppc_405-ar
export AS=ppc_405-as
export LD=ppc_405-ld
export RANLIB=ppc_405-ranlib
export CC=ppc_405-gcc
export NM=ppc_405-nm

./configure --target=powerpc-hardhat-linux
    --host=powerpc-hardhat-linux
    --build=i586-pc-linux-gnu
    --prefix=/opt/hardhat/devkit/ppc/405/target/usr/local
    --exec-prefix=/usr/local
```

You may also need to provide a parameter like `--with-random=/dev/urandom` to
configure as it cannot detect the presence of a random number generating
device for a target system. The `--prefix` parameter specifies where curl
will be installed. If `configure` completes successfully, do `make` and `make
install` as usual.

In some cases, you may be able to simplify the above commands to as little as:

    ./configure --host=ARCH-OS

# REDUCING SIZE

There are a number of configure options that can be used to reduce the size of
libcurl for embedded applications where binary size is an important factor.
First, be sure to set the `CFLAGS` variable when configuring with any relevant
compiler optimization flags to reduce the size of the binary. For gcc, this
would mean at minimum the -Os option, and potentially the `-march=X`,
`-mdynamic-no-pic` and `-flto` options as well, e.g.

    ./configure CFLAGS='-Os' LDFLAGS='-Wl,-Bsymbolic'...

Note that newer compilers often produce smaller code than older versions
due to improved optimization.

Be sure to specify as many `--disable-` and `--without-` flags on the
configure command-line as you can to disable all the libcurl features that you
know your application is not going to need. Besides specifying the
`--disable-PROTOCOL` flags for all the types of URLs your application will not
use, here are some other flags that can reduce the size of the library by
disabling support for some feature:

 - `--disable-alt-svc` (HTTP Alt-Srv)
 - `--disable-ares` (the C-ARES DNS library)
 - `--disable-cookies` (HTTP cookies)
 - `--disable-crypto-auth` (cryptographic authentication)
 - `--disable-dateparse` (date parsing for time conditionals)
 - `--disable-dnsshuffle` (internal server load spreading)
 - `--disable-doh` (DNS-over-HTTP)
 - `--disable-get-easy-options` (lookup easy options at runtime)
 - `--disable-hsts` (HTTP Strict Transport Security)
 - `--disable-http-auth` (all HTTP authentication)
 - `--disable-ipv6` (IPv6)
 - `--disable-libcurl-option` (--libcurl C code generation support)
 - `--disable-manual` (built-in documentation)
 - `--disable-netrc`  (.netrc file)
 - `--disable-ntlm-wb` (NTLM WinBind)
 - `--disable-progress-meter` (graphical progress meter in library)
 - `--disable-proxy` (HTTP and SOCKS proxies)
 - `--disable-pthreads` (multithreading)
 - `--disable-socketpair` (socketpair for async name resolving)
 - `--disable-threaded-resolver`  (threaded name resolver)
 - `--disable-tls-srp` (Secure Remote Password authentication for TLS)
 - `--disable-unix-sockets` (UNIX sockets)
 - `--disable-verbose` (eliminates debugging strings and error code strings)
 - `--disable-versioned-symbols` (versioned symbols)
 - `--enable-symbol-hiding` (eliminates unneeded symbols in the shared library)
 - `--without-brotli` (Brotli on-the-fly decompression)
 - `--without-libpsl` (Public Suffix List in cookies)
 - `--without-nghttp2` (HTTP/2 using nghttp2)
 - `--without-ngtcp2` (HTTP/2 using ngtcp2)
 - `--without-zstd` (Zstd on-the-fly decompression)
 - `--without-libidn2` (internationalized domain names)
 - `--without-librtmp` (RTMP)
 - `--without-ssl` (SSL/TLS)
 - `--without-zlib` (on-the-fly decompression)

The GNU compiler and linker have a number of options that can reduce the
size of the libcurl dynamic libraries on some platforms even further.
Specify them by providing appropriate `CFLAGS` and `LDFLAGS` variables on
the configure command-line, e.g.

    CFLAGS="-Os -ffunction-sections -fdata-sections
            -fno-unwind-tables -fno-asynchronous-unwind-tables -flto"
    LDFLAGS="-Wl,-s -Wl,-Bsymbolic -Wl,--gc-sections"

Be sure also to strip debugging symbols from your binaries after compiling
using 'strip' (or the appropriate variant if cross-compiling). If space is
really tight, you may be able to remove some unneeded sections of the shared
library using the -R option to objcopy (e.g. the .comment section).

Using these techniques it is possible to create a basic HTTP-only libcurl
shared library for i386 Linux platforms that is only 133 KiB in size
(as of libcurl version 7.80.0, using gcc 11.2.0).

You may find that statically linking libcurl to your application will result
in a lower total size than dynamically linking.

Note that the curl test harness can detect the use of some, but not all, of
the `--disable` statements suggested above. Use will cause tests relying on
those features to fail. The test harness can be manually forced to skip the
relevant tests by specifying certain key words on the `runtests.pl` command
line. Following is a list of appropriate key words for those configure options
that aren't automatically detected:

 - `--disable-cookies`          !cookies
 - `--disable-dateparse`        !RETRY-AFTER !CURLOPT_TIMECONDITION !CURLINFO_FILETIME !If-Modified-Since !getdate !-z
 - `--disable-libcurl-option`   !--libcurl
 - `--disable-verbose`          !verbose\ logs

# PORTS

This is a probably incomplete list of known CPU architectures and operating
systems that curl has been compiled for. If you know a system curl compiles
and runs on, that is not listed, please let us know!

## 85 Operating Systems

AIX, AmigaOS, Android, Aros, BeOS, Blackberry 10, Blackberry Tablet OS, Cell
OS, ChromeOS, Cisco IOS, Cygwin, Dragonfly BSD, eCOS, FreeBSD, FreeDOS,
FreeRTOS, Fuchsia, Garmin OS, Genode, Haiku, HardenedBSD, HP-UX, Hurd,
Illumos, Integrity, iOS, ipadOS, IRIX, LineageOS, Linux, Lua RTOS, Mac OS 9,
macOS, Mbed, Micrium, MINIX, MorphOS, MPE/iX, MS-DOS, NCR MP-RAS, NetBSD,
Netware, Nintendo Switch, NonStop OS, NuttX, OpenBSD, OpenStep, Orbis OS,
OS/2, OS/400, OS21, Plan 9, PlayStation Portable, QNX, Qubes OS, ReactOS,
Redox, RICS OS, Sailfish OS, SCO Unix, Serenity, SINIX-Z, Solaris, SunOS,
Syllable OS, Symbian, Tizen, TPF, Tru64, tvOS, ucLinux, Ultrix, UNICOS,
UnixWare, VMS, vxWorks, WebOS, Wii system software, Windows, Windows CE, Xbox
System, z/OS, z/TPF, z/VM, z/VSE

## 22 CPU Architectures

Alpha, ARC, ARM, AVR32, Cell, HP-PA, Itanium, m68k, MicroBlaze, MIPS, Nios,
OpenRISC, POWER, PowerPC, RISC-V, s390, SH4, SPARC, VAX, x86, x86-64, Xtensa
curl security process
=====================

This document describes how security vulnerabilities should be handled in the
curl project.

Publishing Information
----------------------

All known and public curl or libcurl related vulnerabilities are listed on
[the curl website security page](https://curl.se/docs/security.html).

Security vulnerabilities **should not** be entered in the project's public bug
tracker.

Vulnerability Handling
----------------------

The typical process for handling a new security vulnerability is as follows.

No information should be made public about a vulnerability until it is
formally announced at the end of this process. That means, for example that a
bug tracker entry must NOT be created to track the issue since that will make
the issue public and it should not be discussed on any of the project's public
mailing lists. Also messages associated with any commits should not make any
reference to the security nature of the commit if done prior to the public
announcement.

- The person discovering the issue, the reporter, reports the vulnerability on
  [https://hackerone.com/curl](https://hackerone.com/curl). Issues filed there
  reach a handful of selected and trusted people.

- Messages that do not relate to the reporting or managing of an undisclosed
  security vulnerability in curl or libcurl are ignored and no further action
  is required.

- A person in the security team responds to the original report to acknowledge
  that a human has seen the report.

- The security team investigates the report and either rejects it or accepts
  it.

- If the report is rejected, the team writes to the reporter to explain why.

- If the report is accepted, the team writes to the reporter to let him/her
  know it is accepted and that they are working on a fix.

- The security team discusses the problem, works out a fix, considers the
  impact of the problem and suggests a release schedule. This discussion
  should involve the reporter as much as possible.

- The release of the information should be "as soon as possible" and is most
  often synchronized with an upcoming release that contains the fix. If the
  reporter, or anyone else involved, thinks the next planned release is too
  far away, then a separate earlier release should be considered.

- Write a security advisory draft about the problem that explains what the
  problem is, its impact, which versions it affects, solutions or workarounds,
  when the release is out and make sure to credit all contributors properly.
  Figure out the CWE (Common Weakness Enumeration) number for the flaw.

- Request a CVE number from
  [HackerOne](https://docs.hackerone.com/programs/cve-requests.html)

- Update the "security advisory" with the CVE number.

- The security team commits the fix in a private branch. The commit message
  should ideally contain the CVE number.

- The security team also decides on and delivers a monetary reward to the
  reporter as per the bug-bounty polices.

- No more than 10 days before release, inform
  [distros@openwall](https://oss-security.openwall.org/wiki/mailing-lists/distros)
  to prepare them about the upcoming public security vulnerability
  announcement - attach the advisory draft for information with CVE and
  current patch. 'distros' does not accept an embargo longer than 14 days and
  they do not care for Windows-specific flaws.

- No more than 48 hours before the release, the private branch is merged into
  the master branch and pushed. Once pushed, the information is accessible to
  the public and the actual release should follow suit immediately afterwards.
  The time between the push and the release is used for final tests and
  reviews.

- The project team creates a release that includes the fix.

- The project team announces the release and the vulnerability to the world in
  the same manner we always announce releases. It gets sent to the
  curl-announce, curl-library and curl-users mailing lists.

- The security web page on the website should get the new vulnerability
  mentioned.

security (at curl dot se)
------------------------------

This is a private mailing list for discussions on and about curl security
issues.

Who is on this list? There are a couple of criteria you must meet, and then we
might ask you to join the list or you can ask to join it. It really is not a
formal process. We basically only require that you have a long-term presence
in the curl project and you have shown an understanding for the project and
its way of working. You must have been around for a good while and you should
have no plans in vanishing in the near future.

We do not make the list of participants public mostly because it tends to vary
somewhat over time and a list somewhere will only risk getting outdated.

Publishing Security Advisories
------------------------------

1. Write up the security advisory, using markdown syntax. Use the same
   subtitles as last time to maintain consistency.

2. Name the advisory file after the allocated CVE id.

3. Add a line on the top of the array in `curl-www/docs/vuln.pm'.

4. Put the new advisory markdown file in the curl-www/docs/ directory. Add it
   to the git repo.

5. Run `make` in your local web checkout and verify that things look fine.

6. On security advisory release day, push the changes on the curl-www
   repository's remote master branch.

Hackerone
---------

Request the issue to be disclosed. If there are sensitive details present in
the report and discussion, those should be redacted from the disclosure. The
default policy is to disclose as much as possible as soon as the vulnerability
has been published.

Bug Bounty
----------

See [BUG-BOUNTY](https://curl.se/docs/bugbounty.html) for details on the
bug bounty program.
ABI - Application Binary Interface
==================================

 "ABI" describes the low-level interface between an application program and a
 library. Calling conventions, function arguments, return values, struct
 sizes/defines and more.

 [Wikipedia has a longer description](https://en.wikipedia.org/wiki/Application_binary_interface)

## Upgrades

 A libcurl upgrade does not break the ABI or change established and documented
 behavior. Your application can remain using libcurl just as before, only with
 fewer bugs and possibly with added new features.

## Version Numbers

 In libcurl land, you cannot tell by the libcurl version number if that
 libcurl is binary compatible or not with another libcurl version. As a rule,
 we do not break the ABI so you can *always* upgrade to a later version without
 any loss or change in functionality.

## Soname Bumps

 Whenever there are changes done to the library that will cause an ABI
 breakage, that may require your application to get attention or possibly be
 changed to adhere to new things, we will bump the soname. Then the library
 will get a different output name and thus can in fact be installed in
 parallel with an older installed lib (on most systems). Thus, old
 applications built against the previous ABI version will remain working and
 using the older lib, while newer applications build and use the newer one.

 During the first seven years of libcurl releases, there have only been four
 ABI breakages.

 We are determined to bump the SONAME as rarely as possible.  Ideally, we
 never do it again.

## Downgrades

 Going to an older libcurl version from one you are currently using can be a
 tricky thing. Mostly we add features and options to newer libcurls as that
 will not break ABI or hamper existing applications. This has the implication
 that going backwards may get you in a situation where you pick a libcurl that
 does not support the options your application needs. Or possibly you even
 downgrade so far so you cross an ABI break border and thus a different
 soname, and then your application may need to adapt to the modified ABI.

## History

 The previous major library soname number bumps (breaking backwards
 compatibility) happened the following times:

 0 - libcurl 7.1,   August 2000

 1 - libcurl 7.5    December 2000

 2 - libcurl 7.7    March 2001

 3 - libcurl 7.12.0 June 2004

 4 - libcurl 7.16.0 October 2006
# curl man page generator

This is the curl man page generator. It generates a single nroff man page
output from the set of sources files in this directory.

There is one source file for each supported command line option. The output
gets `page-header` prepended and `page-footer` appended. The format is
described below.

## Option files

Each command line option is described in a file named `<long name>.d`, where
option name is written without any prefixing dashes. Like the file name for
the -v, --verbose option is named `verbose.d`.

Each file has a set of meta-data and a body of text.

### Meta-data

    Short: (single letter, without dash)
    Long: (long form name, without dashes)
    Arg: (the argument the option takes)
    Magic: (description of "magic" options)
    Tags: (space separated list)
    Protocols: (space separated list for which protocols this option works)
    Added: (version number in which this was added)
    Mutexed: (space separated list of options this overrides, no dashes)
    Requires: (space separated list of features this requires, no dashes)
    See-also: (space separated list of related options, no dashes)
    Help: (short text for the --help output for this option)
    Example: (example command line, without "curl" and can use `$URL`)
    --- (end of meta-data)

### Body

The body of the description. Only refer to options with their long form option
version, like `--verbose`. The output generator will replace such with the
correct markup that shows both short and long version.

Text written within `*asterisks*` will get shown using italics. Text within
two `**asterisks**` will get shown using bold.

Text that is prefixed with a space will be treated like an "example" and will
be output in monospace.

## Header and footer

`page-header` is the file that will be output before the generated options
output for the master man page.

`page-footer` is appended after all the individual options.

## Generate

`./gen.pl mainpage`

This command outputs a single huge nroff file, meant to become `curl.1`. The
full curl man page.

`./gen.pl listhelp`

Generates a full `curl --help` output for all known command line options.
# libcurl examples

This directory is for libcurl programming examples. They are meant to show
some simple steps on how you can build your own application to take full
advantage of libcurl.

If you end up with other small but still useful example sources, please mail
them for submission in future packages and on the website.

## Building

The Makefile.example is an example makefile that could be used to build these
examples. Just edit the file according to your system and requirements first.

Most examples should build fine using a command line like this:

    `curl-config --cc --cflags --libs` -o example example.c

Some compilers do not like having the arguments in this order but instead
want you do reorganize them like:

    `curl-config --cc` -o example example.c `curl-config --cflags --libs`

**Please** do not use the `curl.se` site as a test target for your
libcurl applications/experiments. Even if some of the examples use that site
as a URL at some places, it does not mean that the URLs work or that we expect
you to actually torture our website with your tests!  Thanks.

## Examples

Each example source code file is designed to be and work stand-alone and
rather self-explanatory. The examples may at times lack the level of error
checks you need in a real world, but that is then only for the sake of
readability: to make the code smaller and easier to follow.
