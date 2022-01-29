# MathJax
## Beautiful math in all browsers

![GitHub release version](https://img.shields.io/github/v/release/mathjax/MathJax-src.svg?sort=semver)
![GitHub release version (v2)](https://img.shields.io/github/package-json/v/mathjax/MathJax/legacy-v2.svg?label=release-v2)
![NPM version](https://img.shields.io/npm/v/mathjax.svg?style=flat)
<a href="http://www.numfocus.org">![powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat)</a>  
![jsdelivr rank](https://flat.badgen.net/jsdelivr/rank/npm/mathjax?color=green)
![jsDelivr hits (npm)](https://img.shields.io/jsdelivr/npm/hm/mathjax)
![npm monthly downloads (full)](https://img.shields.io/npm/dm/mathjax?label=npm)
![npm total downloads](https://img.shields.io/npm/dt/mathjax.svg?style=flat&label=npm%20total)

MathJax is an open-source JavaScript display engine for LaTeX, MathML,
and AsciiMath notation that works in all modern browsers.  It was
designed with the goal of consolidating the recent advances in web
technologies into a single, definitive, math-on-the-web platform
supporting the major browsers and operating systems.  It requires no
setup on the part of the user (no plugins to download or software to
install), so the page author can write web documents that include
mathematics and be confident that users will be able to view it
naturally and easily.  Simply include MathJax and some mathematics in
a web page, and MathJax does the rest.

Some of the main features of MathJax include:

- High-quality display of LaTeX, MathML, and AsciiMath notation in HTML pages

- Supported in most browsers with no plug-ins, extra fonts, or special
  setup for the reader

- Easy for authors, flexible for publishers, extensible for developers

- Supports math accessibility, cut-and-paste interoperability, and other
  advanced functionality

- Powerful API for integration with other web applications

See <http://www.mathjax.org/> for additional details about MathJax,
and <https://docs.mathjax.org> for the MathJax documentation.

## MathJax Components

MathJax version 3 uses files called *components* that contain the
various MathJax modules that you can include in your web pages or
access on a server through NodeJS.  Some components combine all the
pieces you need to run MathJax with one or more input formats and a
particular output format, while other components are pieces that can
be loaded on demand when needed, or by a configuration that specifies
the pieces you want to combine in a custom way.  For usage
instructions, see the [MathJax documentation](https://docs.mathjax.org).

Components provide a convenient packaging of MathJax's modules, but it
is possible for you to form your own custom components, or to use
MathJax's modules directly in a node application on a server.  There
are [web examples](https://github.com/mathjax/MathJax-demos-web)
showing how to use MathJax in web pages and how to build your own
components, and [node
examples](https://github.com/mathjax/MathJax-demos-node) illustrating
how to use components in node applications or call MathJax modules
directly.

## What's in this Repository

This repository contains only the component files for MathJax, not the
source code for MathJax (which are available in a separate [MathJax
source repository](https://github.com/mathjax/MathJax-src/)).  These
component files are the ones served by the CDNs that offer MathJax to
the web.  In version 2, the files used on the web were also the source
files for MathJax, but in version 3, the source files are no longer on
the CDN, as they are not what are run in the browser.

The components are stored in the `es5` directory, and are in ES5 format
for the widest possible compatibility.  In the future, we may make an
`es6` directory containing ES6 versions of the components.

## Installation and Use

### Using MathJax components from a CDN on the web

If you are loading MathJax from a CDN into a web page, there is no
need to install anything.  Simply use a `script` tag that loads
MathJax from the CDN.  E.g.,

``` html
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
```

See the [MathJax
documentation](https://docs.mathjax.org/en/latest/index.html#browser-components),
the [MathJax Web Demos](https://github.com/mathjax/MathJax-demos-web),
and the [MathJax Component
Repository](https://github.com/mathjax/MathJax-demos-web) for more information.

### Hosting your own copy of the MathJax Components

If you want to host MathJax from your own server, you can do so by
installing the `mathjax` package using `npm` and moving the `es5`
directory to an appropriate location on your server:

``` bash
npm install mathjax@3
mv node_modules/mathjax/es5 <path-to-server-location>/mathjax
```

Note that we are still making updates to version 2, so include `@3`
when you install, since the latest chronological version may not be
version 3.

Alternatively, you can get the files via GitHub:

``` bash
git clone https://github.com/mathjax/MathJax.git mj-tmp
mv mj-tmp/es5 <path-to-server-location>/mathjax
rm -rf mj-tmp
```

Then (in either case) you can use a script tag like the following:

``` html
<script id="MathJax-script" async src="<url-to-your-site>/mathjax/tex-chtml.js"></script>
```

where `<url-to-your-site>` is replaced by the URL to the location
where you moved the MathJax files above.

See the
[documentation](https://docs.mathjax.org/en/latest/web/hosting.html)
for details.

### Using MathJax components in a node application

To use MathJax components in a node application, install the `mathjax` package:

``` bash
npm install mathjax@3
```

(we are still making updates to version 2, so you should include `@3`
since the latest chronological version may not be version 3).

Then require `mathjax` within your application:

```js
require('mathjax').init({ ... }).then((MathJax) => { ... });
```
    
where the first `{ ... }` is a MathJax configuration, and the second
`{ ... }` is the code to run after MathJax has been loaded.  E.g.

```js
require('mathjax').init({
  loader: {load: ['input/tex', 'output/svg']}
}).then((MathJax) => {
  const svg = MathJax.tex2svg('\\frac{1}{x^2-1}', {display: true});
  console.log(MathJax.startup.adaptor.outerHTML(svg));
}).catch((err) => console.log(err.message));
```

**Note:** this technique is for node-based application only, not for
browser applications.  This method sets up an alternative DOM
implementation, which you don't need in the browser, and tells MathJax
to use node's `require()` command to load external modules.  This
setup will not work properly in the browser, even if you webpack it or
bundle it in other ways.

See the
[documentation](https://docs.mathjax.org/en/latest/index.html#server-nodejs)
and the [MathJax Node
Repository](https://github.com/mathjax/MathJax-demos-node) for more details.

## Reducing the Size of the Components Directory

Since the `es5` directory contains *all* the component files, so if
you are only planning one use one configuration, you can reduce the
size of the MathJax directory by removing unused components. For
example, if you are using the `tex-chtml.js` component, then you can
remove the `tex-mml-chtml.js`, `tex-svg.js`, `tex-mml-svg.js`,
`tex-chtml-full.js`, and `tex-svg-full.js` configurations, which will
save considerable space.  Indeed, you should be able to remove
everything other than `tex-chtml.js`, and the `input/tex/extensions`,
`output/chtml/fonts/woff-v2`, `adaptors`, `a11y`, and `sre`
directories.  If you are using the results only on the web, you can
remove `adaptors` as well.

If you are not using A11Y support (e.g., speech generation, or
semantic enrichment), then you can remove `a11y` and `sre` as well
(though in this case you may need to disable the assistive tools in
the MathJax contextual menu in order to avoid MathJax trying to load
them when they aren't there).

If you are using SVG rather than CommonHTML output (e.g., `tex-svg.js`
rather than `tex-chtml.js`), you can remove the
`output/chtml/fonts/woff-v2` directory.  If you are using MathML input
rather than TeX (e.g., `mml-chtml.js` rather than `tex-chtml.js`),
then you can remove `input/tex/extensions` as well.


## The Component Files and Pull Requests

The `es5` directory is generated automatically from the contents of the
MathJax source repository.  You can rebuild the components using the
command

``` bash
npm run make-es5 --silent
```

Note that since the contents of this repository are generated
automatically, you should not submit pull requests that modify the
contents of the `es5` directory.  If you wish to submit a modification
to MathJax, you should make a pull request in the [MathJax source
repository](https://github.com/mathjax/MathJax-src).

## MathJax Community

The main MathJax website is <http://www.mathjax.org>, and it includes
announcements and other important information.  A [MathJax user
forum](http://groups.google.com/group/mathjax-users) for asking
questions and getting assistance is hosted at Google, and the [MathJax
bug tracker](https://github.com/mathjax/MathJax/issues) is hosted
at GitHub.

Before reporting a bug, please check that it has not already been
reported.  Also, please use the bug tracker (rather than the help
forum) for reporting bugs, and use the user's forum (rather than the
bug tracker) for questions about how to use MathJax.

## MathJax Resources

* [MathJax Documentation](https://docs.mathjax.org)
* [MathJax Components](https://github.com/mathjax/MathJax)
* [MathJax Source Code](https://github.com/mathjax/MathJax-src)
* [MathJax Web Examples](https://github.com/mathjax/MathJax-demos-web)
* [MathJax Node Examples](https://github.com/mathjax/MathJax-demos-node)
* [MathJax Bug Tracker](https://github.com/mathjax/MathJax/issues)
* [MathJax Users' Group](http://groups.google.com/group/mathjax-users)

# Contributing to MathJax

You are interested in giving us a hand? That's awesome! We've put
together some brief guidelines that should help you get started
quickly and easily.

There are lots and lots of ways to get involved, this document covers:

* [reporting an issue](#reporting-an-issue)
    * [bug reports](#bug-reports)
    * [feature requests](#feature-requests)
    * [change requests](#change-requests)
* [working on MathJax core](#working-on-mathjax-core)
    * [key branches and tags](#key-branches--tags)
    * [submitting pull requests](#submitting-pull-requests)
    * [testing and quality assurance](#testing-and-quality-assurance)
    * [writing documentation](#writing-documentation)
    * [translation](#translation)
* [Conduct](#conduct)


## Reporting An Issue

If you're about to raise an issue because you think you've found a
problem with MathJax, or you'd like to make a request for a new
feature in the codebase, or any other reason, please read this first.

The [MathJax issue
tracker](https://github.com/mathjax/MathJax/issues) is the
preferred channel for [bug reports](#bug-reports), [feature
requests](#feature-requests), [change requests](#change-requests), and
[submitting pull requests](#submitting-pull-requests), but please
respect the following restrictions:

* Please **search for existing issues**. Help us keep duplicate issues
  to a minimum by checking to see if someone has already reported your
  problem or requested your idea.

* Please **do not** use the issue tracker for personal support
  requests (use [the MathJax User Group](https://groups.google.com/forum/#!forum/mathjax-users)).

* Please **be civil**. Keep the discussion on topic and respect the
  opinions of others. See also our [Conduct Guidelines](#conduct)

### Bug Reports

A bug is a *demonstrable problem* that is caused by the code in the
repository. Good bug reports are extremely helpful &mdash; thank you!

Guidelines for bug reports:

1. **Use the GitHub issue search** &mdash; check if the issue has already been
   reported.

2. **Check if the issue has been fixed** &mdash; look for [closed
   issues in the current
   milestone](https://github.com/mathjax/MathJax/issues?q=is%3Aclosed)
   or try to reproduce it using the latest `develop` branch. Please
   note that you will need to
   [compile MathJax and make the components](https://docs.mathjax.org/en/latest/web/hosting.html#getting-mathjax-via-git)
   in order to test MathJax from the source repository.

3. **Share a live sample of the problem** &mdash; without a live page
   it is usually impossible to debug problems; see also the [Bug Report
   Template](#template) below.

4. **Isolate the problem** &mdash; a live sample is a starting point
   but if you want to speed things up, create a [reduced test
   case](https://css-tricks.com/reduced-test-cases/). Be specific
   about your setup (browser, OS versions, etc). Use services like
   [jsbin](http://jsbin.com), [CodePen](http://codepen.io), or
   [jsFiddle](http://jsfiddle.com) to make collaboration on minimal
   test cases easier for everyone.

5. **Include a screenshot/cast as a last resort** &mdash; Is your
   issue about a layout or design feature or bug that is hard to reproduce
   or isolate? Then please provide a screenshot or screencast. Tools
   like [LICEcap](http://www.cockos.com/licecap/) or
   [SauceLabs](http://www.saucelabs.com) allow you to quickly and
   easily record a screencasts. If you make it an animated gif, you can
   embed it directly into your GitHub issue.

6. Use the [Bug Report Template](#template) below or [click this
   link](https://github.com/MathJax/MathJax/issues/new?title=Bug%3A&body=%23%23%23%20Issue%20Summary%0A%0A%23%23%23%20Steps%20to%20Reproduce%0A%0A1.%20This%20is%20the%20first%20step%0A%0AThis%20is%20a%20bug%20because...%0A%0A%23%23%23%20Technical%20details%0A%0A*%20MathJax%20Version%3A%20master%20-%20latest%20commit%3A%20%20INSERT%20COMMIT%20REF%0A*%20Client%20OS%3A%20%0A*%20Browser%3A%20%0A*%20)
   to start creating a bug report with the template automatically.

A good bug report shouldn't leave others needing to request
more information from you. Be sure to include the details of your environment.

<a id="template"></a>

Template Example ([click to use](https://github.com/MathJax/MathJax/issues/new?title=Bug%3A&body=%23%23%23%20Issue%20Summary%0A%0A%23%23%23%20Steps%20to%20Reproduce%0A%0A1.%20This%20is%20the%20first%20step%0A%0AThis%20is%20a%20bug%20because...%0A%0A%23%23%23%20Technical%20details%0A%0A*%20MathJax%20Version%3A%20master%20-%20latest%20commit%3A%20%20INSERT%20COMMIT%20REF%0A*%20Client%20OS%3A%20%0A*%20Browser%3A%20%0A*%20)):

```
Short and descriptive example bug report title

### Issue Summary

A summary of the issue and the browser/OS environment in which it occurs. If
suitable, include the steps required to reproduce the bug.

### Steps to Reproduce

1. This is the first step
2. This is the second step
3. Further steps, etc.

Any other information you want to share that is relevant to the issue
being reported. Especially, why do you consider this to be a bug? What
do you expect to happen instead?

### Technical details:

* MathJax Version: 2.3 (latest commit: f3aaf3a2a3e964df2770dc4aaaa9c87ce5f47e2c)
* Client OS: Mac OS X 10.8.4
* Browser: Chrome 29.0.1547.57
```


### Feature Requests

Feature requests are welcome. Before you submit one, be sure to have:

1. **Used the GitHub search** to check that the feature hasn't already
   been requested.
2. Take a moment to think about whether your idea fits with the scope
   and aims of the project, or if it might better fit being a [custom
   extension](https://github.com/mathjax/MathJax-third-party-extensions).
3. Remember, it's up to *you* to make a strong case to convince the
   project's leaders of the merits of this feature. Please provide as
   much detail and context as possible, this means explaining the use
   case and why it is likely to be common.

### Change Requests

Change requests cover both architectural and functional changes to how
MathJax works. If you have an idea for a new or different dependency,
a refactor, or an improvement to a feature, etc., please be sure to:

1. **Use the GitHub search** to check that someone else didn't get there first.
2. Take a moment to think about the best way to make a case for, and
   explain what you're thinking. Are you sure this shouldn't really be
   a [bug report](#bug-reports) or a [feature
   request](#feature-requests)?  Is it really one idea or is it many?
   What's the context? What problem are you solving? Why is what you
   are suggesting better than what's already there?

## Working on MathJax core

You want to contribute code? We describe how below.  First, note that
the MathJax source code is in the
<https://github.com/mathjax/MathJax-src> repository, not the
<https://github.com/mathjax/MathJax> repository, which contains the
packaged component files for distribution on CDNs and the [mathjax npm
package](https://www.npmjs.com/package/mathjax) (the source code is
included in the [mathjax-full npm
package](https://www.npmjs.com/package/mathjax-full)).

### Key Branches & Tags

MathJax uses several permanent branches in the [MathJax source repository](https://github.com/mathjax/MathJax-src):

- **[develop](https://github.com/mathjax/MathJax-src/tree/develop)**
  is the development branch. All work on the next release happens here
  so you should generally branch off `develop` if you are going to
  submit a pull request. Do **NOT** use this branch for a production
  site.

- **[master](https://github.com/mathjax/MathJax-src)** contains the latest
  release of MathJax. This branch may be used in production. Do 
  **NOT** use this branch to work on MathJax's source.

These branches reflect version 3 of MathJax, which is substantially
different from the version 2 codebase.  Version 2 will continue to be
maintained while web sites transition to version 3, with work being
done using the following branches in the [MathJax distribution
repository](https://github.com/mathjax/MathJax):

- **[legacy-v2-develop](https://github.com/mathjax/MathJax/tree/legacy-v2-develop)**
  is the development branch for changes to the legacy version 2 code.
  Any pull requests for version 2 should be branched from here.  Do
  **NOT** use this branch for a production site.

- **[legacy-v2](https://github.com/mathjax/MathJax/tree/legacy-v2)**
  is the branch that contains any updates to version 2 following
  the release of version 3.  Do **NOT** use this branch to work on
  MathJax's source.

In addition to these branches, MathJax uses tags to identify the
various versions.  These can be checked out to obtain the specified
release; for example, `git checkout 2.7.5` would get you the files for
version 2.7.5 of MathJax.

Note that version 3 is written in Typescript, and so must be compiled
to obtain usable javascript files, and that the components need to be
built once that is done.  See the
[documentation](https://docs.mathjax.org/en/latest/web/hosting.html#getting-mathjax-via-git)
for details. For version 2, the source javascript files are not
compressed until a release is made, so you should use the copies in
the `unpacked` directory during development.
  

### Submitting Pull Requests

Pull requests are welcome. If you're looking to submit a PR for
something that doesn't have an open issue, please consider [raising an
issue](#reporting-an-issue) that your PR can close, especially if
you're fixing a bug. This makes it more likely that there will be
enough information available for your PR to be properly tested and
merged.

##### Need Help?

If you're not completely clear on how to submit/update/*do* Pull
Requests, please check out our [source control
policies](https://github.com/mathjax/MathJax/wiki/Source-control-policies). For
more insights, check the excellent in depth [Git Workflow
guide](https://github.com/TryGhost/Ghost/wiki/Git-Workflow) from
Ghost, in particular

* [Ghost Workflow guide: commit messages](https://github.com/TryGhost/Ghost/wiki/Git-workflow#commit-messages)

### Testing and Quality Assurance

If you're
looking to get involved with the code base and don't know where to
start, checking out and testing a pull request is one of the most
useful things you could do.

These are some [excellent
instructions](https://gist.github.com/piscisaureus/3342247) on
configuring your GitHub repository to allow you to checkout pull
requests in the same way as branches.


### Writing documentation

MathJax's main documentation can be found at [docs.mathjax.org](http://docs.mathjax.org).
The source of the docs is hosted in the
[mathjax/MathJax-docs](http://github.com/mathjax/MathJax-docs) repo here on GitHub.

The documentation is generated using
[Sphinx-doc](http://sphinx-doc.org/) and hosted on [Read the
docs](http://readthedocs.org).  You can clone the repo and submit pull
requests following the [pull-request](#submitting-pull-requests)
guidelines.


### Translation

If you wish to add or update translations of MathJax, please do it on
[TranslateWiki.net](https://translatewiki.net/w/i.php?title=Special:Translate&group=out-mathjax-0-all)
(and while you're there you can help other open source projects,
too).

For bug reports and other questions that don't fit on
TranslateWiki.net, head over to the
[mathjax/mathjax-i18n](https://github.com/mathjax/MathJax-i18n)
repository.

The translation files currently are for version 2, as localization
hasn't been added to version 3 yet.

## Conduct

As a NumFOCUS fiscally sponsored project, MathJax is governed by the
[NumFOCUS code of conduct](https://numfocus.org/code-of-conduct),
which we summarize as follows:

We are committed to providing a friendly, safe and welcoming
environment for all, regardless of gender, sexual orientation,
disability, ethnicity, religion, or similar personal characteristic.

Please be kind and courteous. There's no need to be mean or rude.
Respect that people have differences of opinion and that every design
or implementation choice carries a trade-off and numerous costs. There
is seldom a right answer, merely an optimal answer given a set of
values and circumstances.

Please keep unstructured critique to a minimum. If you have solid
ideas you want to experiment with, make a fork and see how it works.

We will exclude you from interaction if you insult, demean or harass
anyone.  That is not welcome behaviour. We interpret the term
"harassment" as including the definition in the [Unacceptable
Behavior](https://numfocus.org/code-of-conduct#unacceptable-behavior)
section of the [NumFOCUS code of
conduct](https://numfocus.org/code-of-conduct); if you have any lack
of clarity about what might be included in that concept, please read
that definition. In particular, we don't tolerate behavior that
excludes people in socially marginalized groups.

Private harassment is also unacceptable. No matter who you are, if you
feel you have been or are being harassed or made uncomfortable by a
community member, please contact one of the channel ops or any of the
[MathJax](https://github.com/MathJax/MathJax) core team
immediately. Whether you're a regular contributor or a newcomer, we
care about making this community a safe place for you and we've got
your back.

Likewise any spamming, trolling, flaming, baiting, or other
attention-stealing behaviour is not welcome.

We also recommend that you read [discourse's
rules](http://blog.discourse.org/2013/03/the-universal-rules-of-civilized-discourse/)
for further suggestions on appropriate behavior.

## References

* We heavily borrowed from Mozilla and Ghost -- thank you!
  * <https://github.com/TryGhost/Ghost/blob/master/CONTRIBUTING.md>
  * <https://github.com/mozilla/rust/wiki/Note-development-policy>
* <https://github.com/jden/CONTRIBUTING.md/blob/master/CONTRIBUTING.md>
* <http://blog.discourse.org/2013/03/the-universal-rules-of-civilized-discourse/>
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

### Issue Summary

A summary of the issue and the browser/OS environment in which it occurs. If
suitable, include the steps required to reproduce the bug.

### Steps to Reproduce:

1. This is the first step
2. This is the second step
3. Further steps, etc.

Any other information you want to share that is relevant to the issue
being reported. Especially, why do you consider this to be a bug? What
do you expect to happen instead?

### Technical details:

* MathJax Version: 2.3 (latest commit: f3aaf3a2a3e964df2770dc4aaaa9c87ce5f47e2c)
* Client OS: Mac OS X 10.8.4
* Browser: Chrome 29.0.1547.57

### Supporting information:

 * Please supply a link to a (live) minimal example page, when possible.
 * If your issue is with the display of the mathematics produced by MathJax, include a screen snapshot that illustrates the problem, when possible.
 * Check your browser console window for any error messages, and include them here.
 * Include the MathJax configuration you are using, and the script tag that loads MathJax itself.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
