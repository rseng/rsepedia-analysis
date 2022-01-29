# Create React App [![Build Status](https://travis-ci.org/facebook/create-react-app.svg?branch=master)](https://travis-ci.org/facebook/create-react-app) [![PRs Welcome](https://img.shields.io/badge/PRs-welcome-green.svg)](https://github.com/facebook/create-react-app/pulls)

Create React apps with no build configuration.

- [Creating an App](#creating-an-app) – How to create a new app.
- [User Guide](https://facebook.github.io/create-react-app/) – How to develop apps bootstrapped with Create React App.

Create React App works on macOS, Windows, and Linux.<br>
If something doesn’t work, please [file an issue](https://github.com/facebook/create-react-app/issues/new).

## Quick Overview

```sh
npx create-react-app my-app
cd my-app
npm start
```

_([npx](https://medium.com/@maybekatz/introducing-npx-an-npm-package-runner-55f7d4bd282b) comes with npm 5.2+ and higher, see [instructions for older npm versions](https://gist.github.com/gaearon/4064d3c23a77c74a3614c498a8bb1c5f))_

Then open [http://localhost:3000/](http://localhost:3000/) to see your app.<br>
When you’re ready to deploy to production, create a minified bundle with `npm run build`.

<p align='center'>
<img src='https://cdn.rawgit.com/facebook/create-react-app/27b42ac/screencast.svg' width='600' alt='npm start'>
</p>

### Get Started Immediately

You **don’t** need to install or configure tools like Webpack or Babel.<br>
They are preconfigured and hidden so that you can focus on the code.

Just create a project, and you’re good to go.

## Creating an App

**You’ll need to have Node 8.10.0 or later on your local development machine** (but it’s not required on the server). You can use [nvm](https://github.com/creationix/nvm#installation) (macOS/Linux) or [nvm-windows](https://github.com/coreybutler/nvm-windows#node-version-manager-nvm-for-windows) to easily switch Node versions between different projects.

To create a new app, you may choose one of the following methods:

### npx

```sh
npx create-react-app my-app
```

_([npx](https://medium.com/@maybekatz/introducing-npx-an-npm-package-runner-55f7d4bd282b) comes with npm 5.2+ and higher, see [instructions for older npm versions](https://gist.github.com/gaearon/4064d3c23a77c74a3614c498a8bb1c5f))_

### npm

```sh
npm init react-app my-app
```

_`npm init <initializer>` is available in npm 6+_

### Yarn

```sh
yarn create react-app my-app
```

_`yarn create` is available in Yarn 0.25+_

It will create a directory called `my-app` inside the current folder.<br>
Inside that directory, it will generate the initial project structure and install the transitive dependencies:

```
my-app
├── README.md
├── node_modules
├── package.json
├── .gitignore
├── public
│   ├── favicon.ico
│   ├── index.html
│   └── manifest.json
└── src
    ├── App.css
    ├── App.js
    ├── App.test.js
    ├── index.css
    ├── index.js
    ├── logo.svg
    └── serviceWorker.js
```

No configuration or complicated folder structures, just the files you need to build your app.<br>
Once the installation is done, you can open your project folder:

```sh
cd my-app
```

Inside the newly created project, you can run some built-in commands:

### `npm start` or `yarn start`

Runs the app in development mode.<br>
Open [http://localhost:3000](http://localhost:3000) to view it in the browser.

The page will automatically reload if you make changes to the code.<br>
You will see the build errors and lint warnings in the console.

<p align='center'>
<img src='https://cdn.rawgit.com/marionebl/create-react-app/9f62826/screencast-error.svg' width='600' alt='Build errors'>
</p>

### `npm test` or `yarn test`

Runs the test watcher in an interactive mode.<br>
By default, runs tests related to files changed since the last commit.

[Read more about testing.](https://facebook.github.io/create-react-app/docs/running-tests)

### `npm run build` or `yarn build`

Builds the app for production to the `build` folder.<br>
It correctly bundles React in production mode and optimizes the build for the best performance.

The build is minified and the filenames include the hashes.<br>

Your app is ready to be deployed.

## User Guide

You can find detailed instructions on using Create React App and many tips in [its documentation](https://facebook.github.io/create-react-app/).

## How to Update to New Versions?

Please refer to the [User Guide](https://facebook.github.io/create-react-app/docs/updating-to-new-releases) for this and other information.

## Philosophy

- **One Dependency:** There is just one build dependency. It uses Webpack, Babel, ESLint, and other amazing projects, but provides a cohesive curated experience on top of them.

- **No Configuration Required:** You don't need to configure anything. A reasonably good configuration of both development and production builds is handled for you so you can focus on writing code.

- **No Lock-In:** You can “eject” to a custom setup at any time. Run a single command, and all the configuration and build dependencies will be moved directly into your project, so you can pick up right where you left off.

## What’s Included?

Your environment will have everything you need to build a modern single-page React app:

- React, JSX, ES6, TypeScript and Flow syntax support.
- Language extras beyond ES6 like the object spread operator.
- Autoprefixed CSS, so you don’t need `-webkit-` or other prefixes.
- A fast interactive unit test runner with built-in support for coverage reporting.
- A live development server that warns about common mistakes.
- A build script to bundle JS, CSS, and images for production, with hashes and sourcemaps.
- An offline-first [service worker](https://developers.google.com/web/fundamentals/getting-started/primers/service-workers) and a [web app manifest](https://developers.google.com/web/fundamentals/engage-and-retain/web-app-manifest/), meeting all the [Progressive Web App](https://facebook.github.io/create-react-app/docs/making-a-progressive-web-app) criteria. (_Note: Using the service worker is opt-in as of `react-scripts@2.0.0` and higher_)
- Hassle-free updates for the above tools with a single dependency.

Check out [this guide](https://github.com/nitishdayal/cra_closer_look) for an overview of how these tools fit together.

The tradeoff is that **these tools are preconfigured to work in a specific way**. If your project needs more customization, you can ["eject"](https://facebook.github.io/create-react-app/docs/available-scripts#npm-run-eject) and customize it, but then you will need to maintain this configuration.

## Popular Alternatives

Create React App is a great fit for:

- **Learning React** in a comfortable and feature-rich development environment.
- **Starting new single-page React applications.**
- **Creating examples** with React for your libraries and components.

Here are few common cases where you might want to try something else:

- If you want to **try React** without hundreds of transitive build tool dependencies, consider [using a single HTML file or an online sandbox instead](https://reactjs.org/docs/try-react.html).

- If you need to **integrate React code with a server-side template framework** like Rails, Django or Symfony, or if you’re **not building a single-page app**, consider using [nwb](https://github.com/insin/nwb), or [Neutrino](https://neutrino.js.org/) which are more flexible. For Rails specifically, you can use [Rails Webpacker](https://github.com/rails/webpacker). For Symfony, try [Symfony's Webpack Encore](https://symfony.com/doc/current/frontend/encore/reactjs.html).

- If you need to **publish a React component**, [nwb](https://github.com/insin/nwb) can [also do this](https://github.com/insin/nwb#react-components-and-libraries), as well as [Neutrino's react-components preset](https://neutrino.js.org/packages/react-components/).

- If you want to do **server rendering** with React and Node.js, check out [Next.js](https://github.com/zeit/next.js/) or [Razzle](https://github.com/jaredpalmer/razzle). Create React App is agnostic of the backend, and just produces static HTML/JS/CSS bundles.

- If your website is **mostly static** (for example, a portfolio or a blog), consider using [Gatsby](https://www.gatsbyjs.org/) instead. Unlike Create React App, it pre-renders the website into HTML at the build time.

- Finally, if you need **more customization**, check out [Neutrino](https://neutrino.js.org/) and its [React preset](https://neutrino.js.org/packages/react/).

All of the above tools can work with little to no configuration.

If you prefer configuring the build yourself, [follow this guide](https://reactjs.org/docs/add-react-to-an-existing-app.html).

## Contributing

We'd love to have your helping hand on `create-react-app`! See [CONTRIBUTING.md](CONTRIBUTING.md) for more information on what we're looking for and how to get started.

## React Native

Looking for something similar, but for React Native?<br>
Check out [Expo CLI](https://github.com/expo/expo-cli).

## Acknowledgements

We are grateful to the authors of existing related projects for their ideas and collaboration:

- [@eanplatter](https://github.com/eanplatter)
- [@insin](https://github.com/insin)
- [@mxstbr](https://github.com/mxstbr)

## License

Create React App is open source software [licensed as MIT](https://github.com/facebook/create-react-app/blob/master/LICENSE).
---
title: 'LiberTEM: Software platform for scalable multidimensional data processing in transmission electron microscopy'
tags:
  - Python
  - transmission electron microscopy
  - distributed
  - big data
  - MapReduce
authors:
  - name: Alexander Clausen
    orcid: 0000-0002-9555-7455
    affiliation: 1
  - name: Dieter Weber
    orcid: 0000-0001-6635-9567
    affiliation: 1
  - name: Karina Ruzaeva
    affiliation: 1
    orcid: 0000-0003-3610-0989
  - name: Vadim Migunov
    affiliation: "1, 3"
    orcid: 0000-0002-6296-4492
  - name: Anand Baburajan
    affiliation: 5
    orcid: 0000-0002-2870-366X
  - name: Abijith Bahuleyan
    affiliation: 5
    orcid: 0000-0001-5045-5650
  - name: Jan Caron
    affiliation: 1
    orcid: 0000-0002-0873-889X
  - name: Rahul Chandra
    affiliation: 2
    orcid: 0000-0003-2079-5368
  - name: Sayandip Halder
    affiliation: 6
    orcid: 0000-0003-0224-3746
  - name: Magnus Nord
    affiliation: 4
    orcid: 0000-0001-7981-5293
  - name: Knut Müller-Caspary
    affiliation: 1
    orcid: 0000-0002-2588-7993
  - name: Rafal E. Dunin-Borkowski
    affiliation: 1
    orcid: 0000-0001-8082-0647
affiliations:
 - name: Forschungszentrum Jülich, Ernst Ruska-Centre for Microscopy and Spectroscopy with Electrons
   index: 1
 - name: Chandigarh University
   index: 2
 - name: Central Facility for Electron Microscopy (GFE), RWTH Aachen University
   index: 3
 - name: University of Antwerp, EMAT
   index: 4
 - name: Government Engineering College Sreekrishnapuram
   index: 5
 - name: Jadavpur University
   index: 6

date: 12 December 2019
bibliography: paper.bib
---

# Summary

Increases in the data rates of detectors for electron microscopy (EM) have
outpaced increases in network, mass storage and memory bandwidth by two orders
of magnitude between 2009 and 2019 [@Weber2018]. The LiberTEM open source
platform [@Clausen2020] is designed to match the growing performance
requirements of EM data processing [@Weber2020].

# Motivation

The data rate of the fastest detectors for electron microscopy that are
available in 2019 exceeds 50 GB/s, which is faster than the memory bandwidth of
typical personal computers (PCs) at this time. Applications from ten years
before that ran smoothly on a typical PC have evolved into numerical analysis of
complex multidimensional datasets [@Ophus2019] that require distributed
processing on high-performance systems. Furthermore, electron microscopy is
interactive and visual, and experiments performed inside electron microscopes
(so-called in situ experiments) often rely on fast on-line data processing as
the experimental parameters need to be adjusted based on the observation
results. As a consequence, modern data processing systems for electron
microscopy should be designed for very high throughput in combination with short
response times for interactive GUI use and closed-loop feedback. That requires
fundamental changes in the architecture and programming model, and consequently
in the implementation of algorithms and user interfaces for electron microscopy
applications.

# Description

The LiberTEM open source platform for high-throughput distributed processing of
large-scale binary data sets is developed to fulfill these demanding
requirements: Very high throughput on distributed systems, in combination with a
responsive, interactive interface. The current focus for application development
is electron microscopy. Nevertheless, LiberTEM is suitable for any kind of
large-scale binary data that has a hyper-rectangular array layout, notably data
from synchrotrons and neutron sources.

LiberTEM uses a simplified MapReduce [@Dean2008] programming model. It is
designed to run and perform well on PCs, single server nodes, clusters and cloud
services. On clusters it can use fast distributed local storage on
high-performance SSDs. That way it achieves very high aggregate IO performance
on a compact and cost-efficient system built from stock components. On a cluster
with eight microblade nodes we could show a mass storage throughput of 46 GB/s
for a virtual detector calculation.

LiberTEM is supported on Linux, Mac OS X and Windows. Other platforms that allow
installation of Python 3 and the required packages will likely work as well. The
GUI is running in a web browser.

Based on its processing architecture, LiberTEM offers implementations for
various applications of electron microscopy data. That includes basic
capabilities such as integrating over ranges of the input data (virtual
detectors and virtual darkfield imaging, for example), and advanced applications
such as data processing for strain mapping, amorphous materials and phase change
materials. More applications will be added as development progresses.

Compared to established MapReduce-like systems like Apache Spark [@Zaharia2016]
or Apache Hadoop [@Patel2012], it offers a data model that is similar to NumPy
[@Walt2011], suitable for typical binary data from area detectors, as opposed to
tabular data in the case of Spark and Hadoop. It includes interfaces to the
established Python-based numerical processing tools, supports a number of
relevant file formats for EM data, and features optimized data paths for
numerical data that eliminate unnecessary copies and allow cache-efficient
processing.

Compared to tools such as Dask arrays for NumPy-based distributed computations
[@Rocklin2015], LiberTEM is developed towards low-latency interactive feedback
for GUI display as well as future applications for high-throughput distributed
live data processing. As a consequence, data reduction operations in LiberTEM
are not defined as top-down operations like in the case of Dask arrays that are
then broken down into a graph of lower-level operations, but as explicitly
defined bottom-up streaming operations that work on small portions of the input
data. That way, LiberTEM can work efficiently on smaller data portions that fit
into the L3 cache of typical CPUs.

When compared to Dask arrays that try to emulate NumPy arrays as closely as
possible, the data and programming model of LiberTEM is more rigid and closely
linked to the way how the data is structured and how reduction operations are
performed in the back-end. That places restrictions on the implementation of
operations, but at the same time it is easier to understand, control and
optimize how a specific operation is executed, both in the back-end and in
user-defined operations. In particular, it is easier to implement complex
reductions in such a way that they are performed efficiently with a single pass
over the data.

The main focus in LiberTEM has been achieving very high throughput, responsive
GUI interaction and scalability for both offline and live data processing. These
requirements resulted in a distinctive way of handling data in combination with
a matching programming model. Ensuring compatibility and interoperability with
other solutions like [Gatan Microscopy Suite
(GMS)](http://www.gatan.com/products/tem-analysis/gatan-microscopy-suite-software),
Dask, HyperSpy [@Pena2019], pyXem [@Johnstone2019],
[pixStem](https://pixstem.org/) and others is work in progress. They use a
similar data model, which makes integration possible, in principle. As an
example, LiberTEM can be run from within development versions of an upcoming GMS
release that includes an embedded Python interpreter, and it can already generate
efficient Dask.distributed arrays from the data formats it supports.

The online documentation with more details on installation, use, architecture,
applications, supported formats, performance benchmarks and more can be found at
<https://libertem.github.io/LiberTEM/>.

Live data processing for interactive microscopy and automation of experiments is
currently under development. The architecture and programming model of LiberTEM
are already developed in such a way that current applications will work on live
data streams without modification as soon as back-end support is implemented.

# Acknowledgements

This project has received funding from the European Research Council (ERC) under
the European Union’s Horizon 2020 research and innovation programme (grant
agreements No 780487 - VIDEO and No 856538 - 3D MAGiC).

This project has received funding from the European Union’s Horizon 2020
research and innovation programme under grant agreements No 686053 - CritCat and
No 823717 – ESTEEM3.

We gratefully acknowledge funding from the Initiative and Networking Fund of the
Helmholtz Association within the Helmholtz Young Investigator Group moreSTEM
under Contract No. VH-NG-1317 at Forschungszentrum Jülich in Germany.

We gratefully acknowledge funding from the Information & Data Science Pilot
Project "Ptychography 4.0" of the Helmholtz Association.

We kindly acknowledge funding from Google Summer of Code 2019 under the umbrella
of the Python software foundation.

STEMx equipment and software for 4D STEM data acquisition with K2 IS camera
courtesy of Gatan Inc.

Forschungszentrum Jülich is supporting LiberTEM with funding for personnel,
access to its infrastructure and administrative support.

Furthermore, we wish to thank a large number of people who contributed in
various ways to LiberTEM. [We maintain a full and continuously updated list of creators and
contributors online.](https://libertem.github.io/LiberTEM/acknowledgments.html)

# References
