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
|gitter|_ |azure|_ |github|_ |codeclimate|_ |precommit|_ |joss|_ |zenodo|_

.. |gitter| image:: https://badges.gitter.im/Join%20Chat.svg
.. _gitter: https://gitter.im/LiberTEM/Lobby

.. |azure| image:: https://dev.azure.com/LiberTEM/LiberTEM/_apis/build/status/LiberTEM.LiberTEM?branchName=master
.. _azure: https://dev.azure.com/LiberTEM/LiberTEM/_build/latest?definitionId=3&branchName=master

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1477847.svg
.. _zenodo: https://doi.org/10.5281/zenodo.1477847

.. |github| image:: https://img.shields.io/badge/GitHub-GPL--3.0-informational
.. _github: https://github.com/LiberTEM/LiberTEM/

.. |codeclimate| image:: https://api.codeclimate.com/v1/badges/dee042f64380f64737e5/maintainability
.. _codeclimate: https://codeclimate.com/github/LiberTEM/LiberTEM

.. |joss| image:: https://joss.theoj.org/papers/10.21105/joss.02006/status.svg
.. _joss: https://doi.org/10.21105/joss.02006

.. |precommit| image:: https://results.pre-commit.ci/badge/github/LiberTEM/LiberTEM/master.svg
.. _precommit: https://results.pre-commit.ci/latest/github/LiberTEM/LiberTEM/master

LiberTEM is an open source platform for high-throughput distributed processing
of large-scale binary data sets and live data streams using a modified
`MapReduce programming model <https://en.wikipedia.org/wiki/MapReduce>`_. The
current focus is `pixelated
<https://en.wikipedia.org/wiki/Scanning_transmission_electron_microscopy#Universal_detectors_(4D_STEM)>`_
scanning transmission electron microscopy (`STEM
<https://en.wikipedia.org/wiki/Scanning_transmission_electron_microscopy>`_)
:cite:`doi:10.1002/9783527808465.EMC2016.6284,Ophus_2019` and scanning electron
beam diffraction data.

MapReduce-like processing allows to specify an algorithm through two functions:
One function that is mapped on portions of the input data, and another function
that merges (reduces) a partial result from this mapping step into the complete
result. A wide range of TEM and 4D STEM processing tasks can be expressed in
this fashion, see `Applications`_.

The UDF interface of LiberTEM offers a standardized, versatile API to decouple
the mathematical core of an algorithm from details of data source, parallelism,
and use of results. Mapping and merging can be performed in any order and with
different subdivisions of the input data, including running parts of the
calculation concurrently. That means the same implementation can be used in a
wide range of modalities, including massive scaling on clusters. Since each
merge step produces an intermediate result, this style of processing is suitable
for displaying live results from a running calculation in a GUI application and
for `processing live data streams <https://github.com/LiberTEM/LiberTEM-live>`_.
A closed-loop feedback between processing and instrument control can be realized
as well. See `User-defined functions
<https://libertem.github.io/LiberTEM/udf.html>`_ for more details on the
LiberTEM UDF interface.

The LiberTEM back-end offers `high throughput and scalability
<https://libertem.github.io/LiberTEM/architecture.html>`_ on PCs, single server
nodes, clusters and cloud services. On clusters it can use fast distributed
local storage on high-performance SSDs. That way it achieves `very high
aggregate IO performance
<https://libertem.github.io/LiberTEM/performance.html>`_ on a compact and
cost-efficient system built from stock components. All CPU cores and CUDA
devices in a system can be used in parallel.

LiberTEM is supported on Linux, Mac OS X and Windows. Other platforms that allow
installation of Python 3.6+ and the required packages will likely work as well. The
GUI is running in a web browser.

Installation
------------

The short version:

.. code-block:: shell

    $ virtualenv -p python3 ~/libertem-venv/
    $ source ~/libertem-venv/bin/activate
    (libertem-venv) $ python -m pip install "libertem[torch]"

    # optional for GPU support
    # See also https://docs.cupy.dev/en/stable/install.html
    (libertem-venv) $ python -m pip install cupy

Please see `our documentation
<https://libertem.github.io/LiberTEM/install.html>`_ for details!

Alternatively, to run the `LiberTEM Docker image
<https://libertem.github.io/LiberTEM/deployment/clustercontainer.html>`_:

.. code-block:: shell

    $ docker run -p localhost:9000:9000 --mount type=bind,source=/path/to/your/data/,dst=/data/,ro libertem/libertem

or

.. code-block:: shell

    $ singularity run docker://libertem/libertem -- /venv/bin/libertem-server

Deployment for offline data processing on a single-node system for a local user
is thoroughly tested and can be considered stable. Deployment on a cluster is
experimental and still requires some additional work, see `Issue #105
<https://github.com/LiberTEM/LiberTEM/issues/105>`_. Back-end support for live data processing
is still experimental as well, see https://github.com/LiberTEM/LiberTEM-live.

Applications
------------

Since LiberTEM is programmable through `user-defined functions (UDFs)
<https://libertem.github.io/LiberTEM/udf.html>`_, it can be used for a wide
range of processing tasks on array-like data and data streams. The following
applications have been implemented already:

- Virtual detectors (virtual bright field, virtual HAADF, center of mass :cite:`Krajnak2016`,
  custom shapes via masks)
- `Analysis of amorphous materials <https://libertem.github.io/LiberTEM/app/amorphous.html>`_
- `Strain mapping <https://libertem.github.io/LiberTEM-blobfinder/>`_
- `Off-axis electron holography reconstruction <https://libertem.github.io/LiberTEM/app/holography.html>`_
- `Single Side Band ptychography <https://ptychography-4-0.github.io/ptychography/>`_

Some of these applications are available through an `interactive web GUI
<https://libertem.github.io/LiberTEM/usage.html#gui-usage>`_. Please see `the
applications section <https://libertem.github.io/LiberTEM/applications.html>`_
of our documentation for details!

The Python API and user-defined functions (UDFs) can be used for complex
operations such as arbitrary linear operations and other features like data
export. Example Jupyter notebooks are available in the `examples directory
<https://github.com/LiberTEM/LiberTEM/tree/master/examples>`_. If you are having
trouble running the examples, please let us know by filing an issue or
by `joining our Gitter chat <https://gitter.im/LiberTEM/Lobby>`_.

LiberTEM is suitable as a high-performance processing backend for other
applications, including live data streams. `Contact us
<https://gitter.im/LiberTEM/Lobby>`_ if you are interested!

LiberTEM is evolving rapidly and prioritizes features following user demand and
contributions. Currently we are working on `live data processing
<https://github.com/LiberTEM/LiberTEM-live>`_, `integration with Dask arrays and
Hyperspy <https://github.com/LiberTEM/LiberTEM/issues/922>`_, support for sparse
data, and implementing analysis methods for various applications of pixelated
STEM and other large-scale detector data. If you like to influence the direction
this project is taking, or if you'd like to `contribute
<https://libertem.github.io/LiberTEM/contributing.html>`_, please join our
`gitter chat <https://gitter.im/LiberTEM/Lobby>`_ and our `general mailing list
<https://groups.google.com/forum/#!forum/libertem>`_.

File formats
------------

LiberTEM currently opens most file formats used for pixelated STEM. See `our
general information on loading data
<https://libertem.github.io/LiberTEM/formats.html>`_ and `format-specific
documentation
<https://libertem.github.io/LiberTEM/reference/dataset.html#formats>`_ for more
information!

- Raw binary files
- Thermo Fisher EMPAD detector :cite:`Tate2016` files
- `Quantum Detectors MIB format <http://quantumdetectors.com/wp-content/uploads/2017/01/1532-Merlin-for-EM-Technical-Datasheet-v2.pdf>`_
- Nanomegas .blo block files
- Direct Electron DE5 files (HDF5-based) and Norpix SEQ files for `DE-Series <http://www.directelectron.com/de-series/>`_ detectors
- `Gatan K2 IS <https://web.archive.org/web/20180809021832/http://www.gatan.com/products/tem-imaging-spectroscopy/k2-camera>`_ raw format
- Stacks of Gatan DM3 and DM4 files (via `openNCEM <https://github.com/ercius/openNCEM>`_)
- FRMS6 from PNDetector pnCCD cameras :cite:`Simson2015` (currently alpha, gain correction still needs UI changes)
- FEI SER files (via `openNCEM <https://github.com/ercius/openNCEM>`_)
- MRC (via `openNCEM <https://github.com/ercius/openNCEM>`_)
- HDF5-based formats such as Hyperspy files, NeXus and EMD
- TVIPS binary files
- Please contact us if you are interested in support for an additional format!

Detectors (experimental)
------------------------

Currently the Quantum Detectors Merlin camera is supported for live processing.
Support for the Gatan K2 IS camera is in a prototype state. Please
`contact us <https://gitter.im/LiberTEM/Lobby>`_ if you are interested in this
feature! See https://github.com/LiberTEM/LiberTEM-live for more details on live
processing.

License
-------

LiberTEM is licensed under GPLv3. The I/O parts are also available under the MIT
license, please see LICENSE files in the subdirectories for details.
DECODING 12 BIT UINTs
=====================

Image data of the Gatan K2 camera is encoded as 12 bit packed unsigned integers, i.e. each
pixel value occupies 1.5 bytes. Since this format cannot be processed directly by the CPU, it has 
to be unpacked to a suitable format first, for example as uint16.

On CPU-native raw data, the LiberTEM processing back-end achieves staggering throughput 
beyond 12 GB/s with suitable numerics packages like a BLAS implementation or pytorch. An 
optimized 12 bit decoder is one critical component to achieve such performance levels with 
the native K2 format.

The simple algorithm in unpack-12-alex.py achieves decent performance of about 1 GB/s with numba.
However, this is an order of magnitude slower than the speed at which LiberTEM processes suitable
data directly. For that reason we have investigated options to increase speed.

Why the simple algorithm is limited
-----------------------------------

Modern CPUs achieve their optimal throughput when they can use their SIMD instruction sets, 
optimize cache use and read data in predictable patterns aligned with 32-bit or 64-bit boundaries, 
among other things. The simple algorithm is problematic in that context: It reads the input in 
three-byte portions so that each loop iteration reads with a different alignment relative to 32-bit or 
64-bit boundaries. Furthermore, the number of iterations is unknown at compile time. Such a pattern 
makes it hard for a compiler to vectorize the loop, i.e. combine blocks of several loop iterations 
into sets of SIMD instructions.

Approach to vectorize bit unpacking
-----------------------------------

The bit pattern of the 12 bit data stream aligns with the 32-bit or 64-bit boundaries in a repeating
fashion. Specifically, every 96/192 bits (3x uint32/uint64; 8/16 x 12 bit), the input stream and the 
CPU-native view on the data have the same alignment. The general case for any bit-length of the input
data can be calculated by dividing the source resp. working bit length by their greatest common divisor
(gcd). This yields the number of data words in a block, seen from working resp. source perspective.

That means each such 96/192 blocks can be processed with the same sequence of CPU instructions, 
i.e. SIMD instructions can be used to process several blocks together in parallel. 
In theory, an advanced compiler could perhaps in the future recognize such a pattern without human 
intervention. In practice, the code has to be written explicitly to process the data in such 
blocks so that the compiler can recognize and vectorize the pattern.

The conversion gets more complicated with the endian-ness of input data and CPU interpretation. Since x86
CPUs are little-endian, i.e. have the lowest order byte at the lowest address in a data word, the bit significance
in a 32 bit uint is 7 6 5 4 3 2 1 0 | 15 14 13 12 11 10 9 8 | 23 22 21 20 19 18 17 16 | 31 30 29 28 27 26 25 24.
The input data of the K2 is little endian as well, following 
a7 a6 a5 a4 a3 a2 a1 a0 | b3 b2 b1 b0 a11 a10 a9 a8 | b11 b10 b9 b8 b7 b6 b5 b4 and so on. 
If we want good performance and use SIMD instructions, the idea is to read the input data aligned as 
uint32 or uint64. On a LITTLE-ENDIAN machine, the 12 bit portions can be packed and unpacked with relatively simple 
masking and shifting operations. On a big endian machine, the various bits of the input stream have to 
be fished out in smaller units using AND with a mask, shifted, and finally ORd their appropriate 
place in the uint16 output word. That is very much possible with SIMD instructions, but 
hard to understand for a human. 

The Jupyter notebook encode_decode.ipynb contains code to calculate the appropriate indices,
shifts and masks to map every bit of an input block to their appropriate place in the output sequence
for all bit lengths and all permutations of big and little endian. For verification and as a sample,
a hand-written example to encode and decode a block of 12 bit little-endian uints on a little-endian
machine is included in the notebook.

unpack-12-7.cpp and unpack-12-4.cpp contain code that can cover the general case for processing
big endian data on a little endian machine, i.e. arbitrary combinations of input bit length and 
working registers for this endian-ness. This code is indeed 
vectorized by gcc. However, this is likely not optimal yet. In many cases several single-bit 
operations can be ORd together when they have same indices and shifts, and only differ by 
their mask bit. In the compiled version of unpack-12-7.cpp and unpack-12-4.cpp, there was no indication that the
optimizer recognized that several operations can be grouped -- but that's not known for sure. 

For that reason versions with hand-written code that uses the optimized sequence for 12 bit big-endian
on little-endian calculated with the Python version were written. It remains to be tested which version is 
actually faster, in the end.

General considerations for optimizer-friendly code
--------------------------------------------------

Any compiler optimization works best when as much as possible is known or precalculated at 
compile time. The C/C++ code contains various measures that should help the compiler with 
that task. The code generally contains intermediate loops with a number of iterations that is 
known at compile time to help with vectorization and test for a potential "sweet spot" for the block size.
Code to cover an unaligned remaining block is omitted for simplicity. Helper functions are 
designed and declared as constexpr to motivate the compiler to the maximum to evaluate 
them at compile time. That is indeed successful: much of the source code is not part of the 
output machine code or is inlined as compact instruction sequences if the function input 
values depend on the input data. Constants are used instead of variables where possible.

It seemed to help the optimizer when input and output data are not read / written directly 
from/to the buffer, but when they are first copied to a local array, processed into a second 
local array and then copied to the output array.

The various implementations in C++
----------------------------------

How to group and order the operations for a block has room for variation.
The inner loop can process an entire block, or several inner loops are executed sequentially 
to process a specific word of the input or output data in each block for a number of 
blocks in a row. These possibilities are explored in the various C++ implementations with the goal to see
what the compiler makes of them. Furthermore, data types and block sizes can be varied
in the code by changing the appropriate values.

None of the implementations are tested for correctness. They likely contain issues with 
index calculation that don't show up on compile time. Furthermore, they are not 
benchmarked yet. The main goal was to explore how to write code that is vectorized and optimized.

Exploring the compiled result
-----------------------------

Have a look at <https://godbolt.org/>. Select language and compiler, copy and paste the code, 
and see what machine code is generated.

Conclusion
----------

The modified algorithm can indeed be vectorized by the compiler in several of the variations.
Which one is faster is far from certain. That might depend on many different factors, including
cache efficiency for code and data. It is even possible that the simple implementation is the 
fastest, in the end. The work presented here can perhaps serve as a starting point to 
investigate this further.
.. _`concepts`:

Concepts
========

LiberTEM is developed to fulfill the requirements for 4D STEM data processing
and can be used for other processing tasks for n-dimensional detector data as
well.

In 4D STEM, an electron beam is scanned over a sample and for each position,
a 2D image is recorded. That means there are two spatial axes, and two "detector" axes.
We also call the spatial axes "scanning axes" or, more generally, navigation axes. This roughly
corresponds to the batch axis in machine learning frameworks.

The "detector" axes are also called "signal axes", and a single 2D image is also called a frame.

Axis order
----------

We generally follow the numpy convention for axis order, so for a 4D data set,
you could have a :code:`(ny, nx, sy, sx)` tuple describing the shape.

MATLAB users should keep one thing in mind. The navigation axes in Python is the transpose of that of MATLAB. 
In Python, the indices increase linearly with the row. A 3x3 Python matrix is represented in the following way:
 
.. code-block:: python

    [[0,1,2],
    [3,4,5],
    [6,7,8]]
	
`The official "NumPy for Matlab users" documentation`_ might be helpful for Matlab users.

Coordinate system
-----------------

LiberTEM works in pixel coordinates corresponding to array indices. That means
(0, 0) is on the top left corner, the x axis points to the right and the y axis
points to the bottom. This follows the usual plotting conventions for pixel
data.

LiberTEM uses a right-handed coordinate system, which means the z axis points away and positive
rotations are therefore clock-wise.

.. note::
    Please note that the data can be saved with different relation of physical coordinates and
    pixel coordinates. Notably, MIB reference files from Quantum Detectors Merlin cameras have their
    y axis inverted when displayed with LiberTEM. LiberTEM generally
    doesn't deal with such transformations in the numerical back-end.

    In :meth:`~libertem.api.Context.create_com_analysis`, a capability to flip the y axis and rotate
    the shift coordinates is added in version 0.6.0 to support processing MIB files and
    calculate results with physical meaning in electron microscopy, such as the curl and divergence.
    See also :issue:`325`.

    Discussion regarding full support for physical units can be found in :issue:`121`.

Multidimensional data
---------------------

While our GUI is currently limited to 2D visualizations, the Python API does not have that
limitation. You can load data of arbitraty dimensionality and specify an application-specific shape
using the GUI or the Python API, provided our I/O routines support the format. Most of our methods
are currently built for 2D image data, so it should be no problem to load and process for
example time series.

If you want to process data with, for example, 1D or 3D samples, you can write
:ref:`UDFs <user-defined functions>`. Note that in that case, a "frame" is no longer 2D!

.. _The official "NumPy for Matlab users" documentation: https://numpy.org/doc/1.18/user/numpy-for-matlab-users.html#numpy-for-matlab-users
.. _contributing:

Contributing
============

LiberTEM is intended and designed as a collaboratively developed platform for
data analysis. That means all our development is coordinated openly, mostly on
our `GitHub repository <https://github.com/LiberTEM/LiberTEM/>`_ where our code
is hosted. Any suggestions, Issues, bug reports, discussions and code
contributions are highly appreciated! Please let us know if you think we can
improve on something, be it code, communication or other aspects.

Development principles
----------------------

We have a `rather extensive and growing list of things to work on
<https://github.com/LiberTEM/LiberTEM/issues>`_ and therefore have to prioritize
our limited resources to work on items with the largest benefit for our user
base and project. Supporting users who contribute code is most important to us.
Please contact us for help! Furthermore, we prioritize features that create
direct benefits for many current users or open significant new applications.
Generally, we follow user demand with our developments.

For design of new features we roughly follow the `lead user method
<https://en.wikipedia.org/wiki/Lead_user>`_, which means that we develop new
features closely along a non-trivial real-world application in order to make
sure the developments are appropriate and easy to use in practice. The interface
for :ref:`user-defined functions`, as an example, follows the requirements
around implementing and running complex algorithms like `strain mapping
<https://libertem.github.io/LiberTEM-blobfinder/examples.html>`_ for distributed
systems.

Furthermore we value a systematic approach to development with requirements
analysis and evaluation of design options as well as iterative design with fast
test and review cycles.

Code contributions
------------------

We are using `pull requests
<https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests>`_
to accept contributions. Each pull request should focus on a single issue, to
keep the number of changes small and reviewable. To keep your changes organized
and to prevent unrelated changes from disturbing your pull request, create a new
branch for each pull request.

All pull requests should come from a user's personal fork since we don't push to
the main repository for development. See the `GitHub documentation on how to
create and manage forks
<https://docs.github.com/en/get-started/quickstart/fork-a-repo>`_
for details.

Before creating a pull request, please make sure all tests still pass. See
`Running the Tests`_ for more information. You should also update the test suite
and add test cases for your contribution. See the section `Code coverage`_ below
on how to check if your new code is covered by tests.

To make sure our code base stays readable and consistent, we follow a `Code Style`_.

Please update ``packaging/creators.json`` with your author information when you
contribute to LiberTEM for the first time. This helps us to keep track of all
contributors and give credit where credit is due! Please let us know if you
wouldn't like to be credited. Our :ref:`authorship` describes in more detail how
we manage authorship of LiberTEM and related material.

If you are changing parts of LiberTEM that are currently not covered by tests,
please consider writing new tests! When changing example code, which is not run
as part of the tests, make sure the example still runs.

When adding or changing a feature, you should also update the corresponding
documentation, or add a new section for your feature. Follow the current
documentation structure, or ask the maintainers where your new documentation
should end up. When introducing a feature, it is okay to start with a draft
documentation in the first PR, if it will be completed later. Changes of APIs
should update the corresponding docstrings.

Please include version information if you add or change a feature in order to
track and document changes. We use a rolling documentation that documents
previous behavior as well, for example *This feature was added in Version
0.3.0.dev0* or *This describes the behavior from Version 0.3.0.dev0 and onwards.
The previous behavior was this and that*. If applicable, use
`versionadded <https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-versionadded>`_
and related directives.

The changelog for the development branch is maintained as a collection of files
in the :code:`docs/source/changelog/*/` folder structure. Each change should get
a separate file to avoid merge conflicts. The files are merged into the
master changelog when creating a release.

The following items might require an
update upon introducing or changing a feature:

* Changelog snippet in :code:`docs/source/changelog/*/`
* Docstrings
* Examples
* Main Documentation

When you have submitted your pull request, someone from the LiberTEM
organization will review your pull request, and may add comments or ask
questions. 

In case your PR touches I/O code, an organization member may run
the I/O tests with access to test data sets on a separate Azure Agent,
using the following comment on the PR:

.. code-block:: text

    /azp run libertem.libertem-data

If everything is good to go, your changes will be merged and you can
delete the branch you created for the pull request.

.. seealso:: `Guide on understanding the GitHub flow <https://guides.github.com/introduction/flow/>`_

.. seealso:: `How to make a simple GitHub PR (video) <https://www.youtube.com/watch?v=cysuuUtbC6E>`_

.. _`running tests`:

Running the tests
-----------------

Our tests are written using pytest. For running them in a repeatable manner, we
are using tox. Tox automatically manages virtualenvs and allows testing on
different Python versions and interpreter implementations.

This makes sure that you can run the tests locally the same way as they are run
in continuous integration.

After `installing tox <https://tox.readthedocs.io/en/latest/install.html>`_, you
can run the tests on all Python versions by simply running tox:

.. code-block:: shell

    $ tox

Or specify a specific environment you want to run:

.. code-block:: shell

    $ tox -e py36

For faster iteration, you can also run only a part of the test suite, without
using tox. To make this work, first install the test requirements into your
virtualenv:

.. code-block:: shell

   (libertem) $ python -m pip install -r test_requirements.txt

Now you can run pytest on a subset of tests, for example:

.. code-block:: shell

   (libertem) $ pytest tests/test_analysis_masks.py

Or you can run tests in parallel, which may make sense if you have a beefy
machine with many cores and a lot of RAM:

.. code-block:: shell

   (libertem) $ pytest -n auto tests/

See the `pytest documentation
<https://docs.pytest.org/en/latest/how-to/usage.html#specifying-which-tests-to-run>`_
for details on how to select which tests to run. Before submitting a pull
request, you should always run the whole test suite.

Some tests are marked with `custom markers
<https://docs.pytest.org/en/latest/example/markers.html>`_, for example we have
some tests that take many seconds to complete. To select tests to run by these
marks, you can use the `-m` switch. For example, to only run the slow tests:

.. code-block:: shell

   $ tox -- -m slow

By default, these slow tests are not run. If you want to run both slow and all
other tests, you can use a boolean expression like this:

.. code-block:: shell

   $ tox -- -m "slow or not slow"

Another example, to exclude both slow and functional tests:

.. code-block:: shell

   $ tox -- -m "not functional and not slow"

In these examples, ``--`` separates the the arguments of tox (left of ``--``)
from the arguments for pytest on the right. List of marks used in our test
suite:

- `slow`: tests that take much longer than 1 second to run
- `functional`: tests that spin up a local dask cluster

Example notebooks
~~~~~~~~~~~~~~~~~

The example notebooks are also run as part of our test suite using `nbval
<https://nbval.readthedocs.io/en/latest/>`_. The
output saved in the notebooks is compared to the output of re-running the
notebook. When writing an example notebook, this sometimes requires some work,
because certain things will change from run to run. Please check `the nbval
documentation
<https://nbviewer.jupyter.org/github/computationalmodelling/nbval/blob/master/docs/source/index.ipynb>`_
to understand how to ignore such values. See also the :code:`nbval_sanitize.cfg`
file for our currently ignored patterns.

If your notebook outputs floating point values, you may get spurious failures
from precision issues. You can set the precision using the :code:`%precision` ipython
magic, which should be used *after* importing numpy.

You can run the notebook tests as follows:

.. code-block:: shell

    $ TESTDATA_BASE_PATH=/path/to/the/test/data tox -e notebooks

You will need access to certain sample data sets; as most of them are not
published yet, please contact us to get access!

CUDA
~~~~

To run tests that require CuPy using tox, you can specify the CUDA version with the test environment:

.. code-block:: shell

    $ tox -e py36-cuda101

Code coverage
~~~~~~~~~~~~~

After running the tests, you can inspect the test coverage by opening
`htmlcov/index.html` in a web browser. When creating a pull request, the change
in coverage is also reported by the codecov bot. Ideally, the test coverage
should go up with each pull request, at least it should stay the same.

.. _`benchmarking`:

Benchmarking
~~~~~~~~~~~~

LiberTEM uses `pytest-benchmark
<https://pytest-benchmark.readthedocs.io/en/latest/usage.html>`_ to benchmark
certain performance-critical parts of the code. You can run the benchmarks
ad-hoc using

.. code-block:: shell

   $ pytest benchmarks/

The benchmarks for Numba compilation time are disabled by default since Numba
caches compilation results, i.e. one has to make sure that benchmarked functions
were not previously run in the same interpreter. To run them, you can use

.. code-block:: shell

   $ pytest -m compilation benchmarks/

In order to record a complete benchmark run for later comparison, you can use

.. note::
   This requires :code:`tox>=3.15` since we are using generative section names
   in :code:`tox.ini`.

.. code-block:: shell

   $ tox -e benchmark
   $ # alternatively
   $ tox -e benchmark-cuda101
   $ tox -e benchmark-cuda102

This saves the benchmark data as a JSON file in a subfolder of
:code:`benchmark_results`. A process to commit such results and report them in a
convenient fashion is to be developed. See :issue:`198`, feedback welcome!

.. versionadded:: 0.6.0
   First benchmarks included to help resolve :issue:`814`, benchmark coverage will grow over time.

Running tests for the client
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run the testsuite for the client, first install the JavaScript/TypeScript dependencies:

.. code-block:: shell

   $ cd client/
   $ npm install

Then, in the same directory, to run the tests execute:

.. code-block:: shell

   $ npm test -- --coverage

This will run all tests and report code coverage. If you want to run the tests
while developing the client, you can run them in watch mode, which is the
default:

.. code-block:: shell

   $ cd client/
   $ npm test

Code style
----------

We try to keep our code `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_
-compliant, with line-length relaxed to 100 chars, and some rules ignored. See
the flake8 section in :code:`setup.cfg` for the current PEP8 settings. As a
general rule, try to keep your changes in a similar style as the surrounding
code.

Before submitting a pull request, please check the code style by running:

.. code-block:: shell

   $ pre-commit run

You may need to install `pre-commit <https://pre-commit.com/>`_ into your
Python environment first. You can use the following to automatically run
the pre-commit hooks before each commit:

.. code-block:: shell

   $ pre-commit install --install-hooks

We recommend using an editor that can check code style on the fly, such as
`Visual Studio Code <https://code.visualstudio.com/docs/python/linting>`__.

Mypy type annotations
~~~~~~~~~~~~~~~~~~~~~

We are starting to introduce type annotations and checking them in CI with
`mypy <http://mypy-lang.org/>`_.
Adding type annotations improves developer experience, especially by improving
auto completion and type information on hover in IDEs. Type checking is
currently quite lax and opt-in. See the file
:code:`.mypy-checked` for the list of Python files that currently opt in.
When adding new code, please consider adding new modules to this list.

The checks are run with pre-commit on changed files that opt in.
You can run mypy locally on all files that opt in with:

.. code-block:: shell

   $ pre-commit run --all-files mypy

Please note that in many cases the type for classes is specified with a string
instead of the class itself. That allows to import classes for typing only if
type checking is performed. See `the section on forward references in PEP484
<https://www.python.org/dev/peps/pep-0484/#forward-references>`_ for more
information.

For general information on type annotations in Python, including best
practices, please also see `Static Typing with Python
<https://typing.readthedocs.io/en/latest/>`_.

Docstrings
~~~~~~~~~~

The `NumPy docstring guide
<https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_ is
our guideline for formatting docstrings. We are testing docstring code examples
in Continuous Integration using `doctest
<https://docs.python.org/3/library/doctest.html>`_. You can test files by hand
by running :code:`pytest --doctest-modules <pathspec>`.

Building the documentation
--------------------------

Documentation building is also done with tox, see above for the basics. It
requires manual `installation of pandoc <https://pandoc.org/installing.html>`_
on the build system since pandoc can't be installed reliably using pip. To start
the live building process:

.. code-block:: shell

    $ tox -e docs

You can then view the documentation at http://localhost:8008, which will
be rebuilt every time a change to the documentation source code is detected.
Note that changes to the Python source code don't trigger a rebuild, so if
you are working on docstrings, you may have to manually trigger a rebuild,
for example by saving one of the :code:`.rst` files.

You can include code samples with the `doctest sphinx extension
<https://www.sphinx-doc.org/en/master/usage/extensions/doctest.html>`_ and test
them with

.. code-block:: shell

    $ tox -e docs-check

.. _`building the client`:

Building the GUI (client)
-------------------------

The LiberTEM GUI (also called the client) is written in TypeScript, using a combination of
React/Redux/Redux-Saga. The client communicates with the Python API server using
both HTTP and websockets. Because browsers can't directly execute TypeScript,
there is a build step involved, which translates the TypeScript code into
JavaScript that is then understood by the browser. This build step is needed
both for development and then again for building the production version.

If you would like to contribute to the client, you first need to set up the
development environment. For this, first install Node.js. On Linux, we recommend
to `install via package manager
<https://nodejs.org/en/download/package-manager/>`_, on Windows `the installer
<https://nodejs.org/en/download/>`_ should be fine. Choose the current LTS
version.

One you have Node.js installed, you should have the :code:`npm` command available
in your path. You can then install the needed build tools and dependencies by
changing to the client directory and running the install command:

.. code-block:: shell

   $ cd client/
   $ npm install

.. note::

   It is always a good idea to start development with installing the current
   dependencies with the above command. Having old versions of dependencies
   installed may cause the build to fail or cause unpredictable failures.

Once this command finished without errors, you can start a development server
(also from the client directory):

.. code-block:: shell

   $ npm run start

This server watches all source files for changes and automatically starts the
build process. The development server, which listens on port 3000, will only be
able to serve requests for JavaScript and other static files. For handling HTTP
API requests you still need to run the Python :code:`libertem-server` process on
the default port (9000) alongside the development server:

.. code-block:: shell

   $ libertem-server --no-browser

This allows proxying the HTTP API requests from the front-end server to the API
server without opening an additional browser window that could interfere with
the development server.

To learn more about the build process, please see `the README in the client
directory <https://github.com/LiberTEM/LiberTEM/blob/master/client/README.md>`_.

You can then use any editor you like to change the client source files, in the
:code:`client/src` directory. We recommend `Visual Studio Code
<https://code.visualstudio.com/>`_ for its excellent TypeScript support.

To simplify development and installing from a git checkout, we currently always
ship a production build of the client in the git repository. Please always open
your pull request for the client as WIP and include a rebuilt production build
after the PR is approved and ready to merge. You can create it using a tox
shortcut:

.. code-block:: shell

   $ tox -e build_client

This will build an optimized production version of the client and copy it into
:code:`src/libertem/web/client`. This version will then be used when you start a
libertem-server without the client development proxy in front.

Advanced
--------

See more:

.. toctree::
   :maxdepth: 2

   releasing
   dist_tests
.. _authorship:

Authorship policy
=================

Since the performance of scientists is evaluated based on the number and quality
of their publications, authorship is an important topic for software with a
focus on scientific application, such as LiberTEM. Our goal is to give
appropriate credit to each contribution so that the work of contributors becomes
visible as a form of scientific output and is considered in the evaluation of
their productivity as scientists. We hope that this can support and encourage
Open Source software development as a form of scientific work.

Following established practices for scientific papers, we distinguish between
creators, who are persons that contributed in a way that establishes
co-authorship of LiberTEM, and contributors, who are persons that contributed in
other ways, for example through discussions or help.

Creators are persons who fulfill at least one of these criteria:

* Contributed at least one commit to the repository.
* Contributed material such as code or documentation in other ways that don't
  appear directly as commit, including prototype code that served as a basis for
  code in LiberTEM.
* Contributed to the management and direction, for example
  active management and support for contributors who dedicate a significant part
  of their working time to LiberTEM.

We maintain two files, ``packaging/creators.json`` and
``packaging/contributors.json``, where we work to keep track of creators and
contributors. For transparency, we include a short statement of the nature of
the contribution. If you feel like you or someone else should be listed there or
should be represented differently, please file a pull request or contact us!
Accidental discrepancies happen, and many eyes and helping hands help us to keep
the information up to date.

Please note that the GitHub breakdown of contributions is not always reflecting
reality. As an example, commits are only linked to GitHub accounts if the Git
client is configured correctly for that.

All people listed in ``creators.json`` are included as authors in publications
of the software itself, for example `uploads on Zenodo
<https://doi.org/10.5281/zenodo.1477847>`_.

Publications about LiberTEM
---------------------------

For other publications about LiberTEM, such as scientific papers, all
persons with at least one contribution within the past 18 months that qualifies
them for being a creator of LiberTEM will be contacted during the drafting
process of the publication, for example through their GitHub handle or via
e-mail, and offered co-authorship. If they actively consent to being a co-author
in a timely manner, they are included as co-authors under the rules of that
particular publication and medium. For scientific papers, that typically
includes an obligation to approve drafts before submission and being accountable
for all aspects of the work. Co-authors who are not responsive may be excluded
if they don't respond to a reminder, for example to approve a draft, in a timely
fashion.

For publications that only cover specific aspects of LiberTEM, for example a
particular feature, only contributions to those particular aspects are
considered to select creators who will be offered co-authorship.

Co-authors who are not creators of LiberTEM can be included in publications
about LiberTEM following established practices for authorship, for example as
contributors to the content of the publication or to other presented material,
such as scientific applications.

Ordering of author lists
------------------------

The creators and contributors are listed alphabetically in
:ref:`acknowledgments` with a short statement about their contribution.

Alphabetical ordering of author lists is uncommon for scientific papers.
Instead, the position on the author list is used to indicate the relative amount
of contribution of an author. The author that contributed most of the content is
listed first and the author that contributed most guidance is listed last.
Casual readers, in particular encountering an abbreviated author list in a
reference within a citing paper, will assign most credit to the first author. An
alphabetical author list would therefore be unfair towards main contributors
with a name in the middle of the alphabet.

In order to resolve this issue and assign prestigious author positions to the
people who deserve them in a transparent fashion, the authors agree among each
other about their authorship positions for each publication individually.

Since Zenodo assigns DOIs and allows to export a citation for reference managers
that are commonly used in scientific publishing, a Zenodo upload is treated as a
scientific publication and the authors are included in the order of
``packaging/creators.json``. This allows to assign prestigious first author
positions to main contributors when LiberTEM is cited in scientific papers.

Authorship questions should be resolved through discussion in Issues and change
proposals in the form of Pull Requests. Our goal is amicable cooperation and
proper credit for all contributions.

If you have questions or would like to suggest changes to this policy, please
contact us! See :pr:`460` for the initial discussion that lead to establishing
this policy.
.. _`applications`:

Applications
============

.. note::
    Starting from 0.4.0, the application-specific code of LiberTEM will be
    spun out step by step as sub-packages that can be installed independent of
    LiberTEM core. See :ref:`packages` for the current overview of sub-packages.

.. toctree::
   :maxdepth: 2

   app/phasecontrast
   app/amorphous
   app/strain
   app/phasechange
   app/holography

* `Ptychography (external) <https://ptychography-4-0.github.io/ptychography/>`_
.. note::
    The features described below are experimental and under development.

.. testsetup:: *

    import numpy as np
    
    from libertem.api import Context
    from libertem.executor.inline import InlineJobExecutor

    ctx = Context(executor=InlineJobExecutor())
    dataset = ctx.load("memory", datashape=(16, 16, 16, 16), sig_dims=2)

.. _`dask`:

Dask integration
================

.. versionadded:: 0.9.0

Since version 0.9, LiberTEM supports seamless integration with workflows that
are based on Dask arrays. That means Dask arrays can serve as input for LiberTEM, and
UDF computations can produce Dask arrays. Additionally,
:meth:`~libertem.contrib.daskadapter.make_dask_array` can create Dask arrays from LiberTEM datasets.

This can be used to combine features from Hyperspy and LiberTEM in a single analysis workflow:

.. toctree::

   hyperspy-integration

Scheduler
---------

LiberTEM uses Dask.distributed as the default method to execute tasks.
Unless instructed otherwise, LiberTEM keeps the default Dask scheduler as-is and only
uses the Dask :code:`Client` internally to make sure existing workflows keep running
as before. However, for a closer integration it can be beneficial to use the same scheduler
for both LiberTEM and other Dask computations. There are several options for that:

:Set LiberTEM Dask cluster as default Dask scheduler:
    * Use :code:`Context.make_with('dask-make-default')`
    * Pass :code:`client_kwargs={'set_as_default': True}` to
      :meth:`~libertem.executor.dask.DaskJobExecutor.connect` or
      :meth:`~libertem.executor.dask.DaskJobExecutor.make_local`
:Use existing Dask scheduler:
    * Use :code:`Context.make_with('dask-integration')` to start an executor
      that is compatible with the current Dask scheduler.
:Use dask.delayed:
    * :class:`libertem.executor.delayed.DelayedJobExecutor` can
      return UDF computations as Dask arrays. The scheduler will only be
      determined when :code:`compute()` is called downstream using the normal
      mechanisms of Dask.

.. _daskarray:

Load datasets as Dask arrays
----------------------------

The :meth:`~libertem.contrib.daskadapter.make_dask_array` function can generate
a `distributed Dask array <https://docs.dask.org/en/latest/array.html>`_ from a
:class:`~libertem.io.dataset.base.DataSet` using its partitions as blocks. The
typical LiberTEM partition size is close to the optimum size for Dask array
blocks under most circumstances. The dask array is accompanied with a map of
optimal workers. This map should be passed to the :meth:`compute` method in
order to construct the blocks on the workers that have them in local storage.

.. testcode:: make_dask_array

    from libertem.contrib.daskadapter import make_dask_array

    # Construct a Dask array from the dataset
    # The second return value contains information
    # on workers that hold parts of a dataset in local
    # storage to ensure optimal data locality
    dask_array, workers = make_dask_array(dataset)

    result = dask_array.sum(axis=(-1, -2)).compute(workers=workers)

Load Dask arrays as datasets
----------------------------

LiberTEM datasets can be created from Dask arrays by using
:meth:`libertem.api.context.load` with filetype :code:`'dask'`. See
:ref:`daskds` for details. The performance can vary and depends on chunking,
executor, I/O method and Dask array creation method. Please `contact us
<https://gitter.im/LiberTEM/Lobby>`_ or `create an Issue
<https://github.com/LiberTEM/LiberTEM/issues/new>` for questions, bug reports
and other feedback!

This basic example shows running a UDF on a Dask array:

.. testcode:: from_dask

    import dask.array as da
    from libertem.udf.sum import SumUDF

    # Create a Context that integrates well with
    # the current Dask scheduler
    ctx = Context.make_with('dask-integration')

    a = da.random.random((23, 42, 17, 19))
    
    ds = ctx.load('dask', a, sig_dims=2)
    res = ctx.run_udf(dataset=ds, udf=SumUDF())

    assert np.allclose(
        a.sum(axis=(0, 1)).compute(),
        res['intensity'].raw_data
    )

See also :ref:`executors` on how to set up compatible schedulers for
both Dask and LiberTEM.

Create Dask arrays with UDFs
----------------------------

Using a :class:`~libertem.executor.delayed.DelayedJobExecutor` with a
:class:`~libertem.api.Context` lets :class:`~libertem.api.Context.run_udf`
return results as Dask arrays.
In addition to the usual :attr:`~libertem.common.buffers.BufferWrapper.data`
and :attr:`~libertem.common.buffers.BufferWrapper.raw_data` properties, which
provide results as numpy arrays eagerly, the results are made available
as dask arrays using the attributes
:attr:`~libertem.executor.utils.dask_buffer.DaskResultBufferWrapper.delayed_data`
and :attr:`~libertem.executor.utils.dask_buffer.DaskResultBufferWrapper.delayed_raw_data`.
The computation is only performed when the
:code:`compute()` method is called on a Dask array result.

.. testcode:: to_dask

    from libertem.api import Context
    from libertem.executor.delayed import DelayedJobExecutor
    
    from libertem.udf.sumsigudf import SumSigUDF

    ctx = Context(executor=DelayedJobExecutor())
    res = ctx.run_udf(dataset=dataset, udf=SumSigUDF())

    print(res['intensity'].delayed_data)

.. testoutput:: to_dask
   
   dask.array<reshape, shape=(16, 16), dtype=float32, chunksize=(..., ...), chunktype=numpy.ndarray>

This method allows to create large intermediate results of :code:`kind='nav'`
efficiently in the form of Dask arrays. The partial results can stay as
ephemeral array chunks on the workers only while they are needed for downstream
calculations instead of transporting them through the cluster and instantiating
the full result on the main node.

For calls to :class:`~libertem.api.Context.run_udf` which run multiple UDFs
or have with UDFs which return multiple results, it is strongly recommended
to compute results with a single call to :code:`compute()`, in order to re-use
computation over the dataset. By default, Dask does not cache intermediate
results from prior runs, so individual calls to :code:`compute()` require
a complete re-run of the UDFs that were passed to :code:`run_udf`.


Merge function for Dask array results
-------------------------------------

LiberTEM already uses an efficient default merging method to create Dask arrays
for UDFs that only use :code:`kind='nav'` buffers and don't specify their own
:code:`merge()`.

For all UDFs that define their own :code:`merge()`,
:class:`~libertem.executor.delayed.DelayedJobExecutor` will use this
existing :code:`udf.merge()` function to assemble the final results, in the same
way it is used to assemble partial results in each partition. This is carried
out by wrapping the :code:`udf.merge()` in :code:`dask.delayed` calls.

However, the user can also specify a
:meth:`~libertem.udf.base.UDFMergeAllMixin.merge_all()` method on their UDF that
combines all partial results from the workers to the complete result in a single
step. This allows the :class:`~libertem.executor.delayed.DelayedJobExecutor` to
produce a streamlined task tree, which then gives Dask greater scope to
parallelise the merge step and reduce data transfers.

This example shows the task tree first with the built-in wrapper for :code:`udf.merge()`
and then the stramlined one after defining :code:`merge_all()`. In tests the streamlined
variant was twice as fast as the built-in wrapper.

.. testsetup:: merge_all

    import numpy as np
    from libertem.udf.base import UDF
    from libertem.executor.delayed import DelayedJobExecutor

.. testcode:: merge_all

    class MySumUDF(UDF):
        def get_result_buffers(self):
            return {
                'intensity': self.buffer(kind='sig', dtype=self.meta.input_dtype)
            }

        def process_tile(self, tile):
            self.results.intensity[:] += np.sum(tile, axis=0)

        def merge(self, dest, src):
            dest.intensity[:] += src.intensity
    
    ctx = Context(executor=DelayedJobExecutor())
    result = ctx.run_udf(udf=MySumUDF(), dataset=dataset)

    result['intensity'].delayed_data.visualize()

.. image:: ./images/tree-default-merge.png
    :width: 300px
    :alt: Task tree with wrapper for standard merge.

.. testcode:: merge_all

    class MySumMergeUDF(MySumUDF):
        # Define the merge_all() method for the UDF above
        def merge_all(self, ordered_results):
            # List and not generator for NumPy dispatch to work
            chunks = [b.intensity for b in ordered_results.values()]
            # NumPy will dispatch the stacking to the appropriate method
            # for the chunks.
            # See also https://numpy.org/doc/stable/user/basics.dispatch.html
            stack = np.stack(chunks)
            # Perform computation on the stacked chunks
            # equivalent to the normal merge()
            intensity = stack.sum(axis=0)
            
            # Return a dictionary mapping buffer name to new content
            return {'intensity': intensity}
    
    result2 = ctx.run_udf(udf=MySumMergeUDF(), dataset=dataset)

    result2['intensity'].delayed_data.visualize()

.. image:: ./images/tree-merge-all.png
    :width: 300px
    :alt: Task tree with :code:`merge_all()`.

.. testcleanup:: merge_all
    
    assert np.allclose(
        result['intensity'].raw_data,
        result2['intensity'].raw_data,
    )

The argument :code:`ordered_results` is an ordered dictionary of all partial
results for that UDF indexed by the slice for the corresponding dataset
partition in the flattened navigation dimension with ROI applied. The order of
the partial results is such that the slices are increasing through the dataset
navigation dimension, so the merge method can safely concatenate the results in
the case of :code:`'nav'`-shaped results.

CUDA and scheduling
-------------------

A native LiberTEM Dask cluster uses resource tags to schedule work on CPUs or
CUDA devices based on an UDF's capability. For Dask integration a fallback was
implemented that allows running computations that can run on a CPU on a native
LiberTEM Dask cluster without requiring resource tags. However, CUDA-only computations
will require passing the appropriate tags as :code:`resources` argument to :code:`compute()`.

:meth:`libertem.executor.delayed.DelayedJobExecutor.get_resources_from_udfs` returns
the appropriate resources for a given set of UDFs based on their capabilities.
Tips and tricks
===============

This is a collection of various helpful tips that don't fit in elsewhere.

.. _`ssh forwarding`:

Using SSH forwarding
--------------------

As there is no built-in authentication yet, LiberTEM should not listen on a network
port where untrusted parties have access. You can use ssh port forwarding from `localhost` instead
to access LiberTEM from a different computer.

For example with conda:

.. code-block:: shell

     $ ssh -L 9000:localhost:9000 <remote-hostname> "conda activate libertem; libertem-server"

Or, with virtualenv:

.. code-block:: shell

     $ ssh -L 9000:localhost:9000 <remote-hostname> "/path/to/virtualenv/bin/libertem-server"

This makes LiberTEM, which is running on `remote-hostname`, available on your
local host via http://localhost:9000/

Alternatively, you can launch and access LiberTEM on remote systems through
:ref:`jupyter integration`.

Activating ipywidgets in Jupyter
--------------------------------

Some examples use :code:`ipywidgets` in notebooks, most notably the fast
:class:`~libertem.viz.bqp.BQLive2DPlot`. In some cases the corresponding Jupyter
Notebook extension `has to be activated manually
<https://ipywidgets.readthedocs.io/en/stable/user_install.html#installing-in-classic-jupyter-notebook>`_:

.. code-block:: shell

    $ jupyter nbextension enable --py widgetsnbextension
    
For :class:`~libertem.viz.bqp.BQLive2DPlot` to function correctly in JupyterHub
and possibly JupyterLab, the packages :code:`bqplot` and :code:`bqplot-image-gl`
have to be installed in both the notebook environment and in the environment for
the notebook server. This might require admin privileges.

Running in a top-level script
-----------------------------

Since LiberTEM uses multiprocessing, the `script entry point may have to be
protected
<https://docs.python.org/3/library/multiprocessing.html#the-spawn-and-forkserver-start-methods>`_,
most notably on Windows:

.. testcode::

    if __name__ == '__main__':
        # Here goes starting a LiberTEM Context etc
        ...

This is not necessary if LiberTEM is used in a Jupyter notebook or IPython.

Gatan Digital Micrograph and other embedded interpreters
--------------------------------------------------------

If LiberTEM is run from within an embedded interpreter, the following steps
should be taken. This is necessary for Python scripting in Gatan Digital
Micrograph (GMS), for example.

The variable :code:`sys.argv` `may not be set in embedded interpreters
<https://bugs.python.org/issue32573>`_, but it is expected by the
:code:`multiprocessing` module when spawning new processes. This workaround
guarantees that :code:`sys.argv` is set `until this is fixed upstream
<https://github.com/python/cpython/pull/12463>`_:

.. testsetup::

    import sys

.. testcode::

    if not hasattr(sys, 'argv'):
        sys.argv  = []

Furthermore, the `correct executable for spawning subprocesses
<https://docs.python.org/3/library/multiprocessing.html#multiprocessing.set_executable>`_
has to be set.

.. testsetup::

    import multiprocessing
    import sys
    import os

.. testcode::

    multiprocessing.set_executable(
        os.path.join(sys.exec_prefix, 'pythonw.exe'))  # Windows only

In GMS the script may have to run in an additional thread since loading SciPy in
a GMS background thread doesn't work. See https://www.gatan.com/python-faq for
more information.

.. testcode::

    import threading

    def main():
        # Here goes the actual script
        ...

    if __name__ == '__main__':
        # Start the workload "main()" in a thread and wait for it to finish
        th = threading.Thread(target=main)
        th.start()
        th.join()

See `our examples folder
<https://github.com/LiberTEM/LiberTEM/tree/master/examples>`_ for a number of
scripts that work in GMS!

.. _`show warnings`:

Show deprecation warnings
-------------------------

Many warning messages via the :code:`warnings` built-in module are suppressed by
default, including in interactive shells such as IPython and Jupyter. If you'd
like to be informed early about upcoming backwards-incompatible changes, you
should activate deprecation warnings. This is recommended since LiberTEM is
under active development.

.. testcode::

    import warnings

    warnings.filterwarnings("default", category=DeprecationWarning)
    warnings.filterwarnings("default", category=PendingDeprecationWarning)

.. _`profiling tests`:

Profiling long-running tests
----------------------------

Since our code base and test coverage is growing continuously, we should make
sure that our test suite remains efficient to finish within reasonable time
frames.

You can find the five slowest tests in the output of Tox, see :ref:`running tests`
for details. If you are using :code:`pytest` directly, you can use the
:code:`--durations` parameter:

.. code-block:: text

    (libertem) $ pytest --durations=10 tests/
    (...)
    ================= slowest 10 test durations =============================
    31.61s call     tests/udf/test_blobfinder.py::test_run_refine_affinematch
    17.08s call     tests/udf/test_blobfinder.py::test_run_refine_sparse
    16.89s call     tests/test_analysis_masks.py::test_numerics_fail
    12.78s call     tests/server/test_job.py::test_run_job_delete_ds
    10.90s call     tests/server/test_cancel.py::test_cancel_udf_job
     8.61s call     tests/test_local_cluster.py::test_start_local
     8.26s call     tests/server/test_job.py::test_run_job_1_sum
     6.76s call     tests/server/test_job.py::test_run_with_all_zeros_roi
     6.50s call     tests/test_analysis_masks.py::test_numerics_succeed
     5.75s call     tests/test_analysis_masks.py::test_avoid_calculating_masks_on_client
    = 288 passed, 66 skipped, 6 deselected, 2 xfailed, 7 warnings in 260.65 seconds =

Please note that functional tests which involve starting a local cluster have
long lead times that are hard to avoid.

In order to gain more information on what slows down a particular test, you can
install the `pytest-profiling extension
<https://github.com/man-group/pytest-plugins/tree/master/pytest-profiling>`_ and
use it to profile individual slow tests that you identified before:

.. code-block:: text

    (libertem) $ pytest --profile tests/udf/test_blobfinder.py::test_run_refine_affinematch
    (...)
    749921 function calls (713493 primitive calls) in 5.346 seconds

    Ordered by: cumulative time
    List reduced from 1031 to 20 due to restriction <20>

    ncalls  tottime  percall  cumtime  percall filename:lineno(function)
         1    0.000    0.000    5.346    5.346 runner.py:76(pytest_runtest_protocol)
     44/11    0.000    0.000    5.344    0.486 hooks.py:270(__call__)
     44/11    0.000    0.000    5.344    0.486 manager.py:65(_hookexec)
     44/11    0.000    0.000    5.344    0.486 manager.py:59(<lambda>)
     44/11    0.001    0.000    5.344    0.486 callers.py:157(_multicall)
         1    0.000    0.000    5.331    5.331 runner.py:83(runtestprotocol)
         3    0.000    0.000    5.331    1.777 runner.py:172(call_and_report)
         3    0.000    0.000    5.330    1.777 runner.py:191(call_runtest_hook)
         3    0.000    0.000    5.329    1.776 runner.py:219(from_call)
         3    0.000    0.000    5.329    1.776 runner.py:198(<lambda>)
         1    0.000    0.000    5.138    5.138 runner.py:119(pytest_runtest_call)
         1    0.000    0.000    5.138    5.138 python.py:1355(runtest)
         1    0.000    0.000    5.138    5.138 python.py:155(pytest_pyfunc_call)
         1    0.004    0.004    5.137    5.137 test_blobfinder.py:149(test_run_refine_affinematch)
         5    0.159    0.032    3.150    0.630 generate.py:6(cbed_frame)
       245    0.001    0.000    2.989    0.012 masks.py:98(circular)
       245    0.046    0.000    2.988    0.012 masks.py:8(_make_circular_mask)
       245    0.490    0.002    2.941    0.012 masks.py:280(radial_bins)
       245    0.152    0.001    2.229    0.009 masks.py:212(polar_map)
        25    0.001    0.000    1.968    0.079 blobfinder.py:741(run_refine)

    =============================== 1 passed, 1 warnings in 7.81 seconds ============================

.. _`os mismatch`:

Platform-dependent code and remote executor
-------------------------------------------

Platform-dependent code in a lambda function or nested function can lead to
incompatibilities when run on an executor with remote workers, such as the
:class:`~libertem.executor.dask.DaskJobExecutor`. Instead, the function should
be defined as part of a module, for example as a stand-alone function or as a
method of a class. That way, the correct remote implementation for
platform-dependent code is used on the remote worker since only a reference to
the function and not the implementation itself is sent over.

Benchmark Numba compilation time
--------------------------------

One has to capture the very first execution of a jitted function and compare it
with subsequent executions to measure its compilation time. By default,
pytest-benchmark performs calibration runs and possibly warmup rounds that don't
report the very first run.

The only way to completely disable this is to use the `pedantic mode
<https://pytest-benchmark.readthedocs.io/en/latest/pedantic.html>`_ specifying
no warmup rounds and two rounds with one iteration each:

.. code-block:: python

   @numba.njit
    def hello():
        return "world"


    @pytest.mark.compilation
    @pytest.mark.benchmark(
        group="compilation"
    )
    def test_numba_compilation(benchmark):
        benchmark.extra_info["mark"] = "compilation"
        benchmark.pedantic(hello, warmup_rounds=0, rounds=2, iterations=1)

That way the maximum is the first run with compilation, and the minimum is the
second one without compilation. Tests are marked as compilation tests in the
extra info as well to aid later data evaluation. Note that the compilation tests
will have poor statistics since it only runs once. If you have an idea on how to
collect better statistics, please `let us know
<https://github.com/LiberTEM/LiberTEM/issues/new>`_!


Simulating slow systems with control groups
-------------------------------------------

Under Linux, it is possible to simulate a slow system using control groups:

.. code-block:: shell

    sudo cgcreate -g cpu:/slow
    sudo cgset -r cpu.cfs_period_us=1000000 slow
    sudo cgset -r cpu.cfs_quota_us=200000 slow
    sudo chown root:<yourgroup> /sys/fs/cgroup/cpu,cpuacct/slow
    sudo chmod 664 /sys/fs/cgroup/cpu,cpuacct/slow

Then, as a user, you can use :code:`cgexec` to run a command in that control group:

.. code-block:: shell

    cgexec -g cpu:slow pytest tests/

This is useful, for example, to debug test failures that only seem to happen in CI
or under heavy load. Note that tools like :code:`cgcreate` only work with cgroups v1,
with newer distributions using cgroups v2 you may have to adapt these instructions.

.. _`jupyter install`:

Jupyter
-------

To use the Python API from within a Jupyter notebook, you can install Jupyter
into your LiberTEM virtual environment.

.. code-block:: shell

    (libertem) $ python -m pip install jupyter

You can then run a local notebook from within the LiberTEM environment, which
should open a browser window with Jupyter that uses your LiberTEM environment.

.. code-block:: shell

    (libertem) $ jupyter notebook

.. _`jupyterhub install`:

JupyterHub
----------

If you'd like to use the Python API from a LiberTEM virtual environment on a
system that manages logins with JupyterHub, you can easily `install a custom
kernel definition
<https://ipython.readthedocs.io/en/stable/install/kernel_install.html>`_ for
your LiberTEM environment.

First, you can launch a terminal on JupyterHub from the "New" drop-down menu in
the file browser. Alternatively you can execute shell commands by prefixing them
with "!" in a Python notebook.

In the terminal you can create and activate virtual environments and perform the
LiberTEM installation as described above. Within the activated LiberTEM
environment you additionally install ipykernel:

.. code-block:: shell

    (libertem) $ python -m pip install ipykernel

Now you can create a custom ipython kernel definition for your environment:

.. code-block:: shell

    (libertem) $ python -m ipykernel install --user --name libertem --display-name "Python (libertem)"

After reloading the file browser window, a new Notebook option "Python
(libertem)" should be available in the "New" drop-down menu. You can test it by
creating a new notebook and running

.. code-block:: python

    In [1]: import libertem

See also :ref:`jupyter integration` for launching the web GUI from JupyterHub or JupyterLab.
Citing
======

To help us spread the word, please credit and cite LiberTEM in publications where it has been significant. 
Resources for citing LiberTEM are linked through the DOI badge below.

|zenodo|_

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1477847.svg
.. _zenodo: https://doi.org/10.5281/zenodo.1477847

Please refer to our :ref:`authorship` for details regarding authorship and credit for contributors... _`acknowledgments`:

Acknowledgments
===============

We are very grateful for your continuing support for LiberTEM!

Please help us keeping these lists up-to-date and complete! If you feel that you
should be listed here, please contact us. We are grateful for every
contribution, and if your contribution is not listed here we'd like to extend
our apologies and update this as soon as possible.

Creators
~~~~~~~~

The following people in alphabetical order contributed to source code,
documentation, design and management following our :ref:`authorship`.

.. include:: autogenerated/creators.rst

Contributions
~~~~~~~~~~~~~

The following people in alphabetical order contributed to the LiberTEM project
in other ways.

.. include:: autogenerated/contributors.rst

Notable upstream projects
~~~~~~~~~~~~~~~~~~~~~~~~~

`Python <https://www.python.org>`_, `PyData universe <https://pydata.org/>`_,
`Dask.distributed <https://distributed.dask.org/>`_, `PyTorch
<https://pytorch.org/>`_, `NumPy <https://numpy.org/>`_, `OpenBLAS
<https://www.openblas.net/>`_, `Click <https://click.palletsprojects.com/>`_,
`Tornado web <https://www.tornadoweb.org/>`_, `Matplotlib
<https://matplotlib.org/>`_, `Pillow <https://pillow.readthedocs.io/>`_, `H5Py
<https://www.h5py.org/>`_, `Numba <https://numba.pydata.org/>`_, `Psutil
<https://psutil.readthedocs.io/>`_,
`Ncempy`_

`TypeScript <https://www.typescriptlang.org/>`_, `React
<https://reactjs.org/>`_, `React Window <https://react-window.vercel.app/>`_, `Redux
<https://redux.js.org/>`_, `Redux-saga <https://redux-saga.js.org/>`_, `Semantic
UI <https://semantic-ui.com/>`_

Not dependencies, but notable related projects or useful tools: `Hyperspy
<https://hyperspy.org/>`_, `NeXus <https://www.nexusformat.org/>`_, `Apache
Spark <https://spark.apache.org/>`_, `Hadoop file system
<http://hadoop.apache.org/docs/stable/hadoop-project-dist/hadoop-hdfs/HdfsDesign.html>`_,
`Godbolt compiler explorer <https://godbolt.org/>`_, `FIO
<https://github.com/axboe/fio>`_, `pyXem
<https://pyxem.github.io/pyxem-website/>`_, `Nion Swift <https://nionswift.readthedocs.io/>`_,
`Alpaka <http://alpaka-group.github.io/alpaka/>`_

Code licensing
~~~~~~~~~~~~~~

Parts of the SEQ reader are based on the PIMS project,
which has the following license:

.. code-block:: text

    Copyright (c) 2013-2014 PIMS contributors
    https://github.com/soft-matter/pims
    All rights reserved

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
          documentation and/or other materials provided with the distribution.
        * Neither the name of the soft-matter organization nor the
          names of its contributors may be used to endorse or promote products
          derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Confirm message on exit using ctrl+c is implemented by referencing Jupyter
notebook project, which has the following license:

.. code-block:: text

    - Copyright (c) 2001-2015, IPython Development Team
    - Copyright (c) 2015-, Jupyter Development Team

    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.

    Redistributions in binary form must reproduce the above copyright notice, this
    list of conditions and the following disclaimer in the documentation and/or
    other materials provided with the distribution.

    Neither the name of the Jupyter Development Team nor the names of its
    contributors may be used to endorse or promote products derived from this
    software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Funding
~~~~~~~

.. raw:: html

    <style>
        .libertem-ack-item { display: flex; }
        .libertem-ack-item img { display: block; }
        .libertem-ack-item > a { flex-shrink: 0; display: block; padding-right: 20px; }
    </style>

.. |itemstart|  raw:: html

    <div class="libertem-ack-item">

.. |itemend|  raw:: html

    </div>

LiberTEM kindly acknowledges funding and support from the following sources:

ERC Proof-of-Concept grant VIDEO
................................

|itemstart|

.. image:: ./images/EU.jpg
    :width: 100px
    :alt: European Union flag
    :target: https://cordis.europa.eu/project/id/780487

This project has received funding from the European Research Council (ERC) under
the European Union’s Horizon 2020 research and innovation programme (`grant
agreement No 780487
<https://cordis.europa.eu/project/id/780487>`_).
|itemend|

CritCat
.......

|itemstart|

.. image:: ./images/EU.jpg
    :width: 100px
    :alt: European Union flag

This project has received funding from the European Union's Horizon 2020
research and innovation programme under `grant agreement No 686053
<http://www.critcat.eu/>`_.
|itemend|

ESTEEM3
.......

|itemstart|

.. image:: ./images/EU.jpg
    :width: 100px
    :alt: European Union flag

This project has received funding from the European Union's Horizon 2020
research and innovation programme under grant agreement No 823717 – `ESTEEM3
<https://cordis.europa.eu/project/rcn/220936/factsheet/en>`_.
|itemend|

ERC Synergy grant 3D MAGiC
..........................

|itemstart|

.. image:: ./images/EU.jpg
    :width: 100px
    :alt: European Union flag

This project has received funding from the European Research Council (ERC) under
the European Union’s Horizon 2020 research and innovation programme (grant
agreement No 856538).
|itemend|

moreSTEM
........

|itemstart|

.. image:: ./images/Helmholtz.png
    :width: 100px
    :alt: Helmholtz Gemeinschaft Deutscher Forschungszentren

We gratefully acknowledge funding from the `Initiative and Networking Fund of
the Helmholtz Association
<https://www.helmholtz.de/en/about-us/the-association/initiating-and-networking/>`_
within the `Helmholtz Young Investigator Group moreSTEM
<https://morestem.fz-juelich.de/>`_ under Contract No. VH-NG-1317 at
Forschungszentrum Jülich in Germany.
|itemend|

Ptychography 4.0
................

|itemstart|

.. image:: ./images/Helmholtz-lower.png
    :width: 100px
    :alt: Helmholtz Gemeinschaft Deutscher Forschungszentren

We gratefully acknowledge funding from the `Information & Data Science Pilot
Project
<https://www.helmholtz.de/en/research/information-data-science/information-data-science-pilot-projects/pilot-projects-2/>`_
"Ptychography 4.0" of the Helmholtz Association.
|itemend|

Google Summer of Code
.....................

|itemstart|

.. image:: images/GSoC-icon-192.png
    :width: 100px
    :alt: Google Summer of Code logo

We kindly acknowledge funding from `Google Summer of Code 2019 and 2020
<https://summerofcode.withgoogle.com/>`_ under the `umbrella of the Python
software foundation <https://python-gsoc.org/>`_.
|itemend|

Gatan Inc.
..........

|itemstart|

.. image:: images/Gatan-logo-vertical.png
    :width: 100px
    :alt: Gatan Inc.

STEMx equipment and software for 4D STEM data acquisition with K2 IS camera
courtesy of `Gatan Inc <https://www.gatan.com/>`_.
|itemend|

AIDAS
.....

|itemstart|

.. image:: images/AIDAS-logo.png
    :width: 100px
    :alt: AIDAS

LiberTEM development is supported by `AIDAS
<https://www.fz-juelich.de/SharedDocs/Meldungen/PORTAL/EN/2021/2021-06-24-AIDAS_en.html>`_.
|itemend|

Forschungszentrum Jülich, Ernst-Ruska Centrum
.............................................

|itemstart|

.. image:: ./images/FZJ.jpg
    :width: 100px
    :alt: Forschungszentrum Jülich GmbH
    :target: https://www.fz-juelich.de/er-c/EN/Home/home_node.html

Forschungszentrum Jülich is supporting LiberTEM with funding for personnel,
access to its infrastructure and administrative support.
|itemend|
Releasing
=========

This document describes release procedures and infrastructure that is relevant
for advanced contributors. See :ref:`contributing` for information on regular
contributions.

Tagging a version
-----------------

Install dependencies from :code:`scripts/requirements.txt`,
which are used by :code:`scripts/release`. Then call the script with
the :code:`bump` command, with the new version as parameter:

.. code-block:: shell

    $ ./scripts/release bump v0.3.0rc0 --tag

If you are bumping to a .dev0 suffix, omit :code:`--tag` and only pass :code:`--commit`:

.. code-block:: shell

    $ ./scripts/release bump v0.4.0.dev0 --commit

.. note::
   In normal development, the version in the master branch will be x.y.z.dev0,
   if the next expected version is x.y.z. When starting the release process, it
   will be bumped up to x.y.zrc0 (note: no dot before rc!) and possibly
   additional release candidates afterwards (rc1, ..., rcN). These release candidates
   are done mostly to assure our release scripts work as expected and for doing
   additional QA. See below for our QA process.

Issue template: release checklist
---------------------------------

When planning a release, create a new issue with the following checklist:

.. code-block:: text

    # Release checklist

    Issues and pull requests to be considered for this release:
    
    * #XXX
    * #YYY
    * #ZZZ

    ## Before (using a release candidate package)

    * [ ] Review open issues and pull requests
    * [ ] Run full CI pipeline, including slow tests, on [Azure DevOps](https://dev.azure.com/LiberTEM/LiberTEM/_build?definitionId=3)
    * [ ] Handle deprecation, search the code base for `DeprecationWarning`
          that are supposed to be removed in that release.
    * [ ] GUI dependency update with `npm install`
    * [ ] Review https://github.com/LiberTEM/LiberTEM/security/dependabot and update dependencies
    * [ ] Full documentation review and update, including link check using
          ``sphinx-build -b linkcheck "docs/source" "docs/build/html"``
    * [ ] Run complete test suite, including slow tests that are deactivated by default
          and tests that require sample files.
    * [ ] Update the expected version in notes on changes, i.e. from `0.3.0.dev0`
          to `0.3.0` when releasing version 0.3.0.
    * [ ] Update and review change log in `docs/source/changelog.rst`, merging
          snippets in `docs/source/changelog/*/` as appropriate.
    * [ ] Update the JSON files in the ``packaging/`` folder with author and project information
    * [ ] Edit `pytest.ini` to exclude flaky tests temporarily from release builds
    * [ ] Create a release candidate using `scripts/release`. See `scripts/release --help` for details.
    * [ ] `Confirm that wheel, tar.gz, and AppImage are built for the release candidate on
          GitHub <https://github.com/LiberTEM/LiberTEM/releases>`_
    * [ ] Confirm that a new version with the most recent release candidate is created in the
          `Zenodo.org sandbox <https://sandbox.zenodo.org/record/367108>`_ that is ready for submission.
    * [ ] Install release candidate packages in a clean environment
          (for example:
          `python -m pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple 'libertem==0.2.0rc11'`)
    * [ ] Test the release candidate docker image
    * [ ] For the GUI-related items, open in an incognito window to start from a clean slate
    * [ ] Correct version info displayed in info dialogue?
    * [ ] Link check in version info dialogue
    * [ ] Make sure you have test files of all supported types available
        * [ ] Include floats, ints, big endian, little endian, complex raw data
    * [ ] Open each test file
        * [ ] Are parameters recognized correctly, as far as implemented?
        * [ ] Any bad default values?
        * [ ] Does the file open correctly?
        * [ ] Have a look at the dataset info dialogue. Reasonable values?
    * [ ] Perform all analyses on each test file.
        * [ ] Does the result change when the input parameters are changed?
        * [ ] All display channels present and looking reasonable?
        * [ ] Reasonable performance?
        * [ ] Use pick mode.
    * [ ] Re-open all the files
        * [ ] Are the files listed in "recent files"?
        * [ ] Are the parameters filled from the cache correctly?
    * [ ] Try opening all file types with wrong parameters
        * [ ] Proper understandable error messages?
    * [ ] Pick one file and confirm keyboard and mouse interaction for all analyses
        * [ ] Correct bounds check for keyboard and mouse?
    * [ ] Check what happens when trying to open non-existent files or directories in the GUI.
        * [ ] Proper understandable error message?
        * [ ] Possible to continue working?
    * [ ] Shut down libertem-server while analysis is running
        * [ ] Shut down within a few seconds?
        * [ ] All workers reaped?
    * [ ] Check what happens when trying to open non-existent files by scripting.
        * [ ] Proper understandable error message? TODO automate?
    * [ ] Check what happens when opening all file types with bad parameters by scripting
        * [ ] Proper understandable error message? TODO automate?
    * [ ] Run libertem-server on Windows, connect to a remote dask cluster running on Linux,
      open all file types and perform an analysis for each file type.
    * [ ] Use the GUI while a long-running analysis is running
        * [ ] Still usable, decent response times?
    * [ ] Confirm that pull requests and issues are handled as intended, i.e. milestoned and merged
      in appropriate branch.
    * [ ] Final version bump: `./scripts/release bump v0.3.0 --tag`, push to github
    * [ ] After pipeline finishes, write minimal release notes for the [release](https://github.com/liberTEM/LiberTEM/releases) and publish the GitHub release

    ## After releasing on GitHub

    * [ ] Confirm that all release packages are built and release notes are up-to-date
    * [ ] Install release package
    * [ ] Confirm correct version info
    * [ ] Confirm package upload to PyPI
    * [ ] Confirm images and tags on https://hub.docker.com/r/libertem/libertem
    * [ ] Publish new version on zenodo.org
    * [ ] Update documentation with new links, if necessary
        * [ ] Add zenodo badge for the new release to Changelog page
    * [ ] Send announcement message on mailing list
    * [ ] Edit `pytest.ini` to include flaky tests again
    * [ ] Bump version in master branch to next .dev0 (`./scripts/release bump v0.X.0.dev0 --commit`)
    * [ ] Add to institutional publication databases
    * [ ] Add the current LiberTEM version to [CVL](https://github.com/Chasdfracterisation-Virtual-Laboratory/CharacterisationVL-Software>) - add both the singularity and the .desktop file!
.. _`usage documentation`:

GUI usage
=========

Starting the LiberTEM server
----------------------------

The LiberTEM GUI is based on a client-server architecture. To use the LiberTEM GUI, you need to
have the server running on the machine where your data is available. For using LiberTEM from
Python scripts, this is not necessary, see :ref:`api documentation`.

After :doc:`installing LiberTEM <install>`, activate the virtualenv or conda environment.

You can then start the LiberTEM server by running:

.. code-block:: shell

    (libertem) $ libertem-server

By default, this starts the server on http://localhost:9000, which you can verify by the
log output::

    [2018-08-08 13:57:58,266] INFO [libertem.web.server.main:886] listening on localhost:9000

It will then try to open your default web browser to this URL.

There are a few command line options available:

.. include:: autogenerated/libertem-server.help
    :literal:

.. versionadded:: 0.4.0
    :code:`--browser` / :code:`--no-browser` option was added, open browser by default.
.. versionadded:: 0.6.0
    :code:`-l, --log-level` was added.
.. versionadded:: 0.8.0
    :code:`-t, --token-path` was added and :code:`-h, --host` was re-enabled.
.. versionadded:: 0.9.0.dev0
    :code:`--preload` and :code:`--insecure` were added.

To access LiberTEM remotely, you can use :ref:`use SSH forwarding <ssh forwarding>`
or our :ref:`jupyter integration`, if you already have JupyterHub or JupyterLab
set up on a server.

Connecting
----------

.. note::

   The GUI is tested to work on Firefox and Chromium-based browsers for now. If you
   cannot use a compatible browser for some reason, please `file an issue <https://github.com/liberTEM/LiberTEM/issues>`_!

After starting the server, you can open the GUI in your browser. If it didn't open
automatically, you can access it by default at http://localhost:9000 . At the beginning,
the GUI shows a prompt to create a local cluster or connect to a running one.
The number of workers is preset with a number that will likely give optimal
performance on the given machine. You can also select which CUDA devices to use, if you have
any (needs to have a working cupy installation).

..  figure:: ./images/use/create.png

Opening data
------------

After connection to a cluster, LiberTEM shows a button to start browsing for
available files. On a local cluster that's simply the local filesystem.

.. note:: See :ref:`sample data` for publicly available datasets.

..  figure:: ./images/use/browse.png

This opens the file browser dialogue. On top it shows the current directory,
below it lists all files and subdirectories in that directory. You select an
entry by clicking once on it. You can move up one directory with the ".." entry
on top of the list. The file browser is still very basic. Possible improvements
are discussed in `Issue #83 <https://github.com/LiberTEM/LiberTEM/issues/83>`_.
Contributions are highly appreciated! This example opens an HDF5 file :cite:`Zeltmann2019`.

..  figure:: ./images/use/open.png

You can also bookmark locations you frequently need to access, using the
star icon. The bookmarks are then found under "Go to...".

..  figure:: ./images/use/star.png

After selecting a file, you set the type in the drop-down menu at the top of the
dialogue above the file name. After that you set the appropriate parameters that
depend on the file type. Clicking on "Load Dataset" will open the file with the
selected parameters. The interface and internal logic to find good presets based
on file type and available metadata, validate the inputs and display helpful
error messages is still work in progress. Contributions are highly appreciated!

See :ref:`Loading using the GUI` for more detailed instructions and
format-specific information.

..  figure:: ./images/use/type.png

Running analyses
----------------

Once a dataset is loaded, you can add analyses to it. As an example we choose a
"Ring" analysis, which implements a ring-shaped virtual detector.

..  figure:: ./images/use/add_analysis.png

..  figure:: ./images/use/adjust.png


This analysis shows two views on your data: the two detector dimensions on
the left, the scanning dimensions on the right, assuming a 4D-STEM dataset.
For the general case, we also call the detector dimensions the *signal
dimensions*, and the scanning dimensions the *navigation dimensions*.
See also :ref:`concepts` for more information on axes and coordinate system.

Directly after
adding the analysis, LiberTEM starts calculating an average of all the detector
frames. The average is overlaid with the mask representing the virtual detector. The view on the right
will later show the result of applying the mask to the data. In the beginning it
is empty. The first processing might take a while depending on file size and I/O
performance. Fast SSDs and enough RAM to keep the working files in the file
system cache are highly recommended for a good user experience.

You can adjust the virtual detector by dragging the handles in the GUI. Below it
shows the parameters in numerical form. This is useful to extract positions, for
example for scripting.

After clicking "Apply", LiberTEM performs the calculation and shows the result
in scan coordinates on the right side.

..  figure:: ./images/use/apply.png

Instead of average, you can select "Standard Deviation". This calculates
standard deviation of all detector frames.

..  figure:: ./images/use/std_dev.png

If you are interested in individual frames rather than the average, you can
switch to "Pick" mode in the "Mode" drop-down menu directly below the detector
window.

..  figure:: ./images/use/pick.png

In "Pick" mode, a selector appears in the result frame on the right. You can
drag it around with the mouse to see the frames live in the left window. The
picked coordinates are displayed along with the virtual detector parameters
below the frame window on the left.

..  figure:: ./images/use/pick_frame.png

If you are interested in a limited region, the ROI dropdown provides the option
to select a rectangular region. For example if you select "Rect", the
average/standard deviation is calculated over all images that lie inside selected
rectangle. You can adjust the rectangle by dragging the handles in the GUI.

..  figure:: ./images/use/rect.png

Some analyses, such as the Center of Mass (COM) analysis, can render the result
in different ways. You can select different result channels in the "Channel" drop-down menu
below the right window.

..  figure:: ./images/use/image.png

.. _`download results`:

Downloading results
-------------------

After an analysis has finished running, you can download the results. Clicking the download button
below the analysis will open a dialog:

..  figure:: ./images/use/download-btn.png

In the download dialog, you can choose between different file formats, and separately
download the available results.

..  figure:: ./images/use/download-modal.png

You can also download a Jupyter notebook corresponding to the analysis and
continue working with the same parameters using scripting.

.. figure:: ./images/use/download-jupyter.png

It's also possible to copy individual cells of Jupyter notebook directly from GUI, with an option
to copy the complete source code.

.. figure:: ./images/use/copy-jupyter.png

Keyboard controls
~~~~~~~~~~~~~~~~~

You can use arrow keys to change the coordinate parameters of any analysis. To
do this, click on the handle you want to modify, and then use the arrow keys to
move the handle. Hold shift to move in larger steps.

Application-specific documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For more applications, like strain mapping and crystallinity analysis, please
see the :doc:`Applications <applications>` section.
.. _`performance`:

Performance
===========

LiberTEM is designed and optimized for very high throughput in order to process large-scale datasets with fast turn-around. 

..  figure:: ./images/performance/IO-bound-madmax.png

    Near linear throughput scaling with the number of nodes for an IO-bound virtual detector calculation on a 480 GiB file using direct IO on Supermicro Microcloud 5038MD-H8TRF, Intel(R) Xeon(R) CPU D-1541 @ 2.10GHz, Intel Ethernet Controller 10G X550T, RAID 0 of 2x Samsung SSD 970 EVO 2TB, 32 GB RAM, cluster head node with Intel(R)Xeon(R) W-2195 CPU @ 2.30GHz.

..  figure:: ./images/performance/CPU-bound-single-node.png

    Benchmark results performing a CPU-bound virtual detector calculation using memory-mapped IO on files in the file system cache.
.. _`architecture`:

Architecture
============

.. image:: ./images/architecture.svg

LiberTEM currently focuses on pixelated STEM and scanning electron beam
diffraction data processing, both interactive and offline. The processing
back-end supports any n-dimensional binary data. As concrete supported
operations, we started with everything that can be expressed as the application
of one or more masks and summation, i.e. virtual detector, center of mass etc.
These operations are embarrassingly parallel and can be scaled to a distributed
system very well. Furthermore, we support :ref:`user-defined functions` that
process well-defined subsets of a dataset following a simplified `MapReduce
programming model <https://en.wikipedia.org/wiki/MapReduce>`_.

For our task, data locality is one of the most important factors for achieving
good performance and scalability. With a traditional distributed storage
solution (like Lustre or NFS), the network will quickly become the bottleneck.

LiberTEM is distributing the data to the local storage of each compute node. One
possible implementation is using the `Hadoop filesystem (HDFS)`_, although we
are `working on a transparent caching layer
<https://github.com/LiberTEM/LiberTEM/issues/136>`_ as an alternative. The
general idea is to split the dataset into (usually disjoint) partitions, which
are assigned to worker nodes.

The execution is structured into Tasks and Jobs. A Job represents the
computation on a whole dataset, and is divided into Tasks for each partition.
The scheduler executes Tasks on the available worker nodes. For fast execution
on each node, the Task reads the data in small Tiles (~1MB).

For distributing the workload, we are using `dask.distributed
<https://distributed.readthedocs.io/>`_. The `Future` API allows us to control
our computation in a flexible way, with little overhead. With dask Futures, we
can assure that computation on a partition of the dataset is scheduled on the
node(s) that hold the partition on their local storage.

.. _Hadoop filesystem (HDFS): https://hadoop.apache.org/docs/r3.1.0/


For ingesting data into the cluster, a `caching layer
<https://github.com/LiberTEM/LiberTEM/issues/136>`_ (WIP) will transparently
read a dataset from a primary source (via a shared network file system, HTTP,
...) and stores it on fast local storage in a format that is best suited for
efficient processing. The cached data can also be pre-processed, for example for
offset correction or applying a gain map.

An important part of the architecture is the API server. Through the API server,
the client gets access to the resources of the cluster, by running analyses. It
uses a protocol based on HTTP and/or websockets. Processing is initiated by HTTP
calls, and results are streamed back to the browser via web sockets.

The API server keeps some state in memory, like information about the currently
opened datasets. Traditionally this would be done with an external in-memory
database, but for ease of deployment, we decided to integrate this into the API
server.

Processing is done in an asynchronous fashion; if you start an analysis using
the HTTP API, the request immediately returns, but you get notifications about
status changes and results on the websocket channel, or you can explicitly query
the API server about a specific analysis. API calls in the Python API are
synchronous for keeping it easy to use. Please `contact us
<https://gitter.im/LiberTEM/Lobby>`_ or `add a comment to Issue #216
<https://github.com/LiberTEM/LiberTEM/issues/216>`_ if you are interested in an
asynchronous Python API for LiberTEM!

As UI, we use a web-based interface. This allows LiberTEM to work
in cloud environment as well as locally on a single node. We can benefit from
existing FOSS frameworks and infrastructure for communication, authentication
etc. of the web.

LiberTEM is also suited for running on your laptop or workstation. In this case, 
all parts can run on a single node. We can also skip the caching step, if the data
is already available locally.

When taking care to avoid needless copying and buffering, we can achieve native
throughput on each node. With NVMe SSDs, this means we can :ref:`process multiple
gigabytes per second per node <performance>`.


Mathematical expression for applying masks
------------------------------------------

The most basic operation with pixelated STEM data is multiplying each frame
element-wise with a mask of the same size and summing up the result for each
frame. This way one can implement virtual detectors, center of mass and
darkfield imaging, to name a few examples.

Since each frame is processed individually and there's a 1:1 match between 
mask element and detector pixel, the processing can be simplified by
flattening the 4D dataset to a 2D matrix, i.e. a list of vectors where each
vector is one flattened detector frame.
Correspondingly, each mask is flattened to a 1D vector as well.
The result of element-wise  multiplication and summation is a 1D vector for 
each mask, where each entry corresponds to the result of one frame. 
To display the result, the data can be re-shaped to the scan dimension again.

With those extra dimensions flattened, let ``c`` be the number of pixels per
frame, ``n`` the number of frames and ``m`` the number of masks. Let (``A``) be
a ``n x c`` matrix with the scanned data, and ``B`` a ``c x m`` matrix with the
masks. Applying the masks in ``B`` to the detector data in ``A`` can be
expressed as a `rectangular matrix product
<https://en.wikipedia.org/wiki/Matrix_multiplication#Definition>`_ of ``A`` and
``B``.

That means we can use the BLAS function `GEMM
<https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms#Level_3>`_. to
process the flattened data and we can choose from  wide array of highly
optimized implementations of this operation. LiberTEM supports both dense and
sparse masks (`sparse.pydata.org <https://sparse.pydata.org>`_ and `scipy.sparse
<https://docs.scipy.org/doc/scipy/reference/sparse.html>`_) for this purpose.
This functionality is available through
:meth:`~libertem.udf.masks.ApplyMasksUDF`,
:meth:`~libertem.api.Context.create_mask_analysis` and a number of more
specialized analysis functions in :class:`~libertem.api.Context`.
.. _`sample data`:

================
Sample Datasets
================

Bullseye and circular probe diffraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Scanning convergent beam electron diffraction data of gold nanoparticles
(:code:`4DSTEM_experiment/data/datacubes/polyAu_4DSTEM/data`) and simulated
strained gold (:code:`4DSTEM_experiment/data/datacubes/simulation_4DSTEM/data`)
with one file using a standard circular aperture and another using a bullseye
aperture.

:Link:
    https://zenodo.org/record/3592520 :cite:`ophus_colin_2019_3592520,Zeltmann2019`
:Format:
    HDF5 (uint16)
:Dimension:
    4D (100, 84, 250, 250)
:Size:
    2.1 GB

Electron Bessel beam diffraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Scanning convergent beam electron diffraction with ring-shaped aperture and
overlapping diffraction orders.

:Link:
    https://zenodo.org/record/2566137 :cite:`giulio_guzzinati_2019_2566137,Guzzinati2019`
:Format:
    Stack of DM3 (currently only scripting)
:Dimension:
    3D
:Size:
    2.6 GB

.. _`hires STO`:

High-resolution 4D STEM dataset of SrTiO3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This dataset can be used to test various analysis methods for high-resolution 4D
STEM, including phase contrast methods such as ptychography.

:Link:
    https://zenodo.org/record/5113449 :cite:`strauch_achim_2021_5113449`
:Format:
    MIB (6 bit)
:Dimension:
    4D (128, 128, 256, 256)
:Size:
    177 MB

.. _`synthetic STO`:

Synthetic 4D STEM dataset based on a SrTiO3 supercell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This dataset allows to investigate phase contrast methods for 4D scanning
transmission electron microscopy, such as ptychography.

:Link:
    https://zenodo.org/record/5113235 :cite:`strauch_achim_2021_5113235`
:Format:
    RAW (float32)
:Dimension:
    4D (100, 100, 596, 596)
:Size:
    14.2 GB

Creating random data
~~~~~~~~~~~~~~~~~~~~~~~

Random data can be generated in the following way. It should be kept in mind
that the data generated in this way can only be used for simple testing as it
has no physical significance.

**Raw file:**

.. testsetup:: sampledataraw

    import os
    import tempfile
    raw_temp = tempfile.TemporaryDirectory()
    os.chdir(raw_temp.name)

.. testcode:: sampledataraw

    # Create sample raw file
    import numpy as np
    sample_data = np.random.randn(16, 16, 16, 16).astype("float32")
    sample_data.tofile("raw_sample.raw")

.. testcode:: sampledataraw

    # Load through Python API
    from libertem.api import Context
    if __name__ == '__main__':
      ctx = Context()
      ds = ctx.load("raw", path="./raw_sample.raw", nav_shape=(16, 16), dtype="float32", sig_shape=(16, 16))

.. testcleanup:: sampledataraw

    os.chdir("..")
    raw_temp.cleanup()

**HDF5 file:**

.. testsetup:: sampledatahdf5

    import os
    import tempfile
    hdf5_temp = tempfile.TemporaryDirectory()
    os.chdir(hdf5_temp.name)

.. testcode:: sampledatahdf5

    # Create sample HDF5 file
    import h5py
    import numpy as np
    file = h5py.File('hdf5_sample.h5','w')
    sample_data = np.random.randn(16,16,16,16).astype("float32")
    dataset = file.create_dataset("dataset",(16,16,16,16), data=sample_data)
    file.close()

.. testcode:: sampledatahdf5

    # Load through Python API
    from libertem.api import Context
    if __name__ == '__main__':
      ctx = Context()
      ds = ctx.load("hdf5", path="./hdf5_sample.h5", ds_path="/dataset")

.. testcleanup:: sampledatahdf5

    os.chdir("..")
    hdf5_temp.cleanup()

Alternatively, you can enter the parameters (scan_size, dtype, detector_size)
directly into the load dialog of the GUI. For more details on loading, please
check :ref:`loading data`.
GSoC 2020 ideas
===============

LiberTEM is participating in the `Google Summer of Code
<https://summerofcode.withgoogle.com/>`_ as a sub-organization of the `Python
Software Foundation <https://python-gsoc.org/>`_. As a student, you can get paid
by Google for three months, have fun working on an interesting open source
software project, gain real-world development experience, and do something that
looks nice on your CV!

* Check out our description and project ideas below
* Contact us if you'd like to work on LiberTEM
* Prepare a `proposal <https://python-gsoc.org/index.html#apply>`_ together with us
* You submit your application at the Google Summer of Code homepage to the Python
  Software Foundation organization, naming LiberTEM as the sub-organization.

Why LiberTEM
--------------

`LiberTEM <.>`_ is an open source platform for high-throughput distributed
processing of pixelated scanning transmission electron microscopy (STEM) data.
It is created to deal with the terabytes of data that modern high-speed
high-resolution detectors for electron microscopy can produce. Our
:doc:`architecture` page describes in more detail how exactly it works.

..  figure:: ./images/Principle.png
    :scale: 50%
    :alt: In pixelated STEM, a full diffraction image is recorded for each scan position.

    *In pixelated STEM, a sample is scanned with a focused electron beam, and a
    full image of the transmitted beam is recorded for each scan position. The
    result is a four-dimensional data hypercube. This application can generate
    tremendous amounts of data from high-resolution scans with a high-speed
    high-resolution detector.*

The project started in the beginning of 2018 and is currently attracting more
and more users because it is orders of magnitude faster and more scalable than
established solutions.

Working on LiberTEM will give you experience in developing distributed systems
for high-performance data processing with Python. You can learn how to profile
an application and optimize performance in a targeted way. See
:ref:`performance` for benchmarking results. LiberTEM has its roots in electron
microscopy, but can be adopted for other tasks that involve high-throughput
data-parallel processing of very large binary data sets.

..  figure:: ./images/Future.png
    :alt: Envisioned future architecture of LiberTEM

    *LiberTEM currently implements distributed offline data processing as shown
    on the right of this figure, and is designed to be extended to
    high-throughput distributed live data processing as illustrated on the
    left.*

If you work on our GUI, you'll learn how a responsive web application for big
data analytics can be built with a front-end based on TypeScript, React and
Redux, and an asynchronous Python back-end based on Tornado and
dask.distributed.

Working on the application side will give you experience in Python-based big
data analytics of large-scale binary data sets with a focus on imaging, physics
and materials science with industry-leading throughput and efficiency.

About us
--------

Alex is an experienced software engineer, systems administrator, Python
programmer, web developer, and expert on profiling and performance optimization.
He focuses on the implementation side of LiberTEM. 

* https://github.com/sk1p


Dieter has an interdisciplinary background in materials science, computer
science, product development, product management and business administration. He
is taking care of the application and business side of LiberTEM. 

* https://github.com/uellue
* https://www.facebook.com/uellue
* https://www.linkedin.com/in/uellue/

We'd be happy to share our experience with you!

How to reach us
---------------

The easiest path is our Gitter channel: https://gitter.im/LiberTEM/Lobby

E-Mail: `Dieter Weber <mailto:d.weber@fz-juelich.de>`_ `Alexander Clausen
<mailto:a.clausen@fz-juelich.de>`_

Just drop a message! We are based in Germany (UTC+1 / UTC+2) and are generally
active during normal working hours.

Getting started
---------------

If you have questions, please ask freely: Supporting users and contributors has
a high priority for us and your questions help us improve our documentation.

Installation
~~~~~~~~~~~~

Please see our :ref:`installation` instructions for details! Forking our
repository, cloning the fork and :ref:`installing from a git clone` are the
recommended setup if you will be contributing significant amounts of code. Our
:ref:`contributing` page has some information that can
help you get started with development.

Currently, we are still working on getting suitable sample files online. Please
contact us to get interesting sample data to work on!

What to work on
~~~~~~~~~~~~~~~

Our `issue tracker can give you a broad overview
<https://github.com/LiberTEM/LiberTEM/issues>`_ of what we have on our plate.
We've marked a number of `Good first issues
<https://github.com/LiberTEM/LiberTEM/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22>`_
that might be relatively easy to solve and can help you get introduced to our
code base. Further below we list a few specific ideas.


Writing your GSoC application
-----------------------------

The PYTHON GSOC website has a good overview for the application process:
https://python-gsoc.org/index.html#apply https://python-gsoc.org/students.html
with links to additional resources. Please remember to name the sub-organization
*LiberTEM* in your proposal!

To get an inspiration on how to write your application, Mailman
(link broken :samp:`https://turnbull.sk.tsukuba.ac.jp/Blog/SPAM.txt`) gives a few general ideas.
For us it is most important to know why you'd like to contribute to LiberTEM,
what your experience level is and what you'd like to learn in order to find a
good match for a project. Generally, we like targeted applications and
applicants who contact us directly early-on. We are happy to work with you for
writing up a project idea!

For Python GSoC it is a general requirement to already contribute a pull request
to a sub-organization before submitting a project idea. Please contact us if
you'd like some help with that! `Improving our examples
<https://github.com/LiberTEM/LiberTEM/tree/master/examples>`_ in code,
description and presentation would be both relatively easy and quite useful for
us. You could perform the update with our help, hunt
down discrepancies and suggest improvements. Please contact us for the
corresponding data to run the examples!

Other small and relatively straightforward issues include :issue:`81,267,506`

Project ideas
-------------

These are somewhat larger work items. Some of them can keep you busy for the
entire time. Please feel free to suggest your own ideas as well! Just working on
a number of smaller features and getting a good cross-sectional experience of
LiberTEM can work as well.

1. **Beginner**: Implement rotation in GUI and back-end for center of mass analysis :issue:`31`.
    Currently, the implementation assumes that the detector edges are parallel
    to the sample X and Y coordinate directions. This is mostly, but not always
    the case. In this project you implement an additional parameter for rotation
    both in the back-end and GUI. That includes a discussion with users of the
    center of mask analysis on how to determine and verify this parameter and
    what the interface should be like.

    *Skills:* Communication, Python, NumPy and TypeScript programming, simple
    web GUI development.

    *Domain knowledge:* simple linear algebra, basic optics.

    *Primary contact:* Dieter (@uellue)

2. **Beginner**: Code snippet for analysis in GUI :issue:`158`.
    Currently, our GUI only offers limited capabilities. Most notably, color
    scale, scale bars and exporting results are not implemented. In this
    project, you implement functions that generate a code snippet in the GUI
    ready for copy & paste. Pasting this snippet in, for example, a Jupyter
    notebook allows to use the analysis with the same parameters as in the GUI
    in a scripting environment that gives users more flexibility.

    *Skills:* Python

    *Domain knowledge:* Learning some LiberTEM basics

    *Primary contact:* Dieter (@uellue), Alex (@sk1p)

3. **Intermediate**: Implement an analysis workflow for `RDF mapping <https://publikationen.bibliothek.kit.edu/1000056485/5249497>`_.
    This can give you experience with the product development, design and
    application side of software engineering, and applied data science. A major
    part of the work is first figuring out *what* to implement together with our
    users and domain experts, and then *how* to implement it. You can decide how
    far you take it: A detailed requirements document, a technical
    specification, a prototype, or a full production-grade implementation? All
    of that is useful for us.

    *Skills:* Communication, software development methodology, Python and NumPy programming.
    
    *Domain knowledge:* Math, statistics, image processing and physics are of advantage.

    *Primary contact:* Dieter (@uellue)

4. **Intermediate**: Allow reshaping datasets into a custom shape :issue:`441`.
    Data in files is not always stored in the native shape, or different
    representations may be possible depending on the application. The dataset
    implementation and GUI should allow specifying a different shape than the
    layout in the dataset.

    *Skills:* Python, NumPy and TypeScript programming.

    *Domain knowledge:* None

    *Primary contact:* Alex (@sk1p)


5. (*Removed since being implemented*: Set number of threads and workers dynamically for UDFs :issue:`546`.)

6. **Beginner/Intermediate/Advanced**: Compression survey :issue:`387`.
    Analyze high-throughput compression techniques, dive into lz4/zstd, blosc
    etc., compare against existing file formats.
    
    *Beginner level*: Test a number of established compression algorithms on typical
    data sets in terms of compression ratio, compression speed and decompression speed.

    *Intermediate level*: Implement the compression in the LiberTEM caching layer.

    *Advanced*: Explore your own ideas regarding compression.

    With this project you can improve your understanding of compression
    techniques for the easier levels, and low-level optimization and programming
    for the advanced level.

    *Skills:* Programming in Python, profiling. C or Numba programming for advanced level.
    
    *Domain knowledge:* Good understanding of computer architecture for the advanced level.

    *Contact:* Dieter (@uellue), Alex (@sk1p)

7. **Intermediate**: Explore automated benchmarks in detail :issue:`198`.
    This will help us to catch performance regressions. In our experience,
    running a benchmark requires a reproducible, undisturbed environment and
    comparison to good reference data. For that reason we see it as more
    challenging than automated tests for functionality and correctness. You
    could run benchmarks in CI and observe variance, and record and present
    benchmark results over time.

    *Skills:* Programming, profiling, visualization.
    
    *Domain knowledge:* Continuous integration and automation tools.

    *Primary contact:* Alex (@sk1p)

8. **Intermediate**: Editor for masks :issue:`47`.
    Currently, the masks in the GUI are limited to a few simple shapes, while
    the back-end allows arbitrary masks. You could implement an online mask
    editor to give users more flexibility on designing masks. Part of the task
    would be a requirements analysis with experts for the scientific
    application, and an analysis if any existing code like
    https://react-designer.github.io/react-designer/
    https://two.js.org/examples/ or http://fabricjs.com/controls-customization
    can possibly be used. This project would be mostly implemented in
    TypeScript.

    *Skills:* Programming in TypeScript, GUI development, basic computer graphics knowledge.
    
    *Domain knowledge:* --

    *Contact:* Dieter (@uellue), Alex (@sk1p)

9. **Intermediate**: Deploy LiberTEM with kubernetes :issue:`105,484`.
    Help us set up a helm chart and documentation to deploy a LiberTEM cluster
    with kubernetes. The subject is fairly new to us and we'd appreciate your
    help, in particular if you already have experience with kubernetes.

    *Skills:* Systems administration and automation.
    
    *Domain knowledge:* kubernetes

    *Primary contact:* Alex (@sk1p)

10. **Intermediate/Advanced**: Proper schemas, validation and automatic form generation for analysis parameters :issue:`316`.
     This feature will make it easier to implement new types of analysis in the
     GUI. This is a cross-section through Python and TypeScript, though we could
     also split off the more react-y part. Does not require NumPy knowledge, or
     domain knowledge. Python/TypeScript required. General WebDev experience
     could help.

     *Skills:* Systematic thinking and abstraction, Python and TypeScript programming, web development.

     *Domain knowledge:* --

     *Primary contact:* Alex (@sk1p)

11. **Intermediate/Advanced**: Custom math kernel for bit masks :issue:`26`.
     Currently, binary masks are first converted to floating point and then used
     in a dot product. NumPy uses GEMM from a third-party BLAS implementation for
     this. This could be accelerated significantly with a Numba-based custom GEMM
     implementation that can work on bit masks directly. Furthermore, such a
     custom Numba-based GEMM kernel has potential other uses in LiberTEM:
     :issue:`555`.

     *Skills:* Python, Numba

     *Domain knowledge:* Optimization, efficient matrix product implementations.

     *Contact:* Dieter (@uellue), Alex (@sk1p)

12. **Advanced**: Live visualization of large binary data :issue:`134`.
     Basically an efficient/zoomable/user-friendly/fully-featured replacement for
     our visualization. Requires a cross-section of different technologies from
     Python/numpy/threading over HTTP/websockets to Canvas/WebGL. Could be spun
     off into its own project if it is successful! This is a larger project that
     can be split into smaller individual parts. If you are interested, we should
     discuss about setting a scope that suits your interests.

     *Skills:* Python and TypeScript programming, web development, asynchronous
     and parallel programming, numerical processing, visualization.
    
     *Domain knowledge:* Experience with similar projects and frameworks like for
     example `GR <https://gr-framework.org/>`_ desirable. Knowledge of `GIS
     <https://en.wikipedia.org/wiki/Geographic_information_system>`_ could
     potentially be useful.

     *Contact:* Dieter (@uellue), Alex (@sk1p)

13. **Advanced**: Enable user-defined functions based on WebAssembly :issue:`199`.
     This would allow users to write user-defined functions in their favorite
     compiled language and is a step towards using LiberTEM independent of
     Python.

     *Skills:* Python and compiled languages.
    
     *Domain knowledge:* Experience with WebAssembly would be useful.

     *Contact:* Dieter (@uellue), Alex (@sk1p)


.. _`loading data`:

Loading data
============

To efficiently handle files larger than main memory, LiberTEM never loads the
whole data set into memory. Calling the :meth:`~libertem.api.Context.load`
function only opens the data set and gives back a handle; running an analysis
with :meth:`~libertem.api.Context.run` or a UDF with
:meth:`~libertem.api.Context.run_udf` then streams the data from mass storage.

See :ref:`sample data` for publicly available datasets for testing.

There are two main ways of opening a data set in LiberTEM: using the GUI, or the
Python API.

Loading through the API
~~~~~~~~~~~~~~~~~~~~~~~

In the API, you can use :meth:`libertem.api.Context.load`. The general
pattern is:

.. code-block:: python

   ctx = Context()
   ctx.load("typename", path="/path/to/some/file", arg1="val1", arg2=42)

So, you need to specify the data set type, the path, and dataset-specific
arguments. These arguments are documented below.

For most file types, it is possible to automatically detect the type and
parameters, which you can trigger by using :code:`"auto"` as file type:

.. code-block:: python

   ctx.load("auto", path="/path/to/some/file")

For the full list of supported file formats with links to their reference
documentation, see :ref:`supported formats` below.

.. _`Loading using the GUI`:

Loading using the GUI
~~~~~~~~~~~~~~~~~~~~~

Using the GUI, mostly the same parameters need to be specified, although some
are only available in the Python API. Tuples (for example for :code:`nav_shape`)
have to be entered as separated values into the fields. You can hit a comma to jump to
the next field. We follow the NumPy convention here and specify the "fast-access" dimension
last, so a value of :code:`42`, :code:`21` would mean the same as specifying
:code:`(42, 21)` in the Python API, setting :code:`y=42` and :code:`x=21`. Note that the GUI
is currently limited to 2D visualizations, while the scripting API can handle more
general cases.

See also :ref:`the concepts section <concepts>`.

Common parameters
~~~~~~~~~~~~~~~~~

There are some common parameters across data set types:

`name`
  The name of the data set, for display purposes. Only used in the GUI.
`nav_shape`
  In the GUI, we generally support visualizing data containing rectangular 2D scans. For
  all the dataset types, you can specify a nav_shape as a tuple `(y, x)`. If the dataset
  isn't 4D, the GUI can reshape it to 4D. When using the Python API, you are free to
  use n-dimensional `nav_shape`, if the data set and chosen analysis supports it.
`sig_shape`
  In the GUI, you can specify shape of the detector as :code:`height`, :code:`width`, but
  when using the Python API, it can be of any dimensionality.
`sync_offset`
  You can specify a `sync_offset` to handle synchronization or acquisition problems.
  If it's positive, `sync_offset` number of frames will be skipped from start.
  If it's negative, `abs(sync_offset)` number of blank frames will be inserted at start.
`io_backend`
  Different methods for I/O are available in LiberTEM, which can influence performance. 
  See :ref:`io backends` for details.
  

.. _`supported formats`:

Supported formats
~~~~~~~~~~~~~~~~~

LiberTEM supports the following file formats out of the box, see links for details:

* :ref:`mib`
* :ref:`raw binary`
* :ref:`dm format`
* :ref:`empad`
* :ref:`k2is`
* :ref:`frms6`
* :ref:`blo`
* :ref:`ser`
* :ref:`hdf5`
* :ref:`seq`
* :ref:`mrc`
* :ref:`tvips`

Furthermore, two alternative mechanisms exist for interfacing LiberTEM with data loaded
elsewhere in Python via other libraries:

- a memory data set can be constructed from a NumPy array for testing
  purposes. See :ref:`memory` for details.
- a Dask data set can be constructed from a Dask array. Depending on the
  method used to construct the source array this can achieve good performance.
  See :ref:`daskds` for details.
  
Running distributed tests
=========================

To test the behavior of the distributed parts of LiberTEM, there are now a
few automated test cases available. Because of the distributed nature of the
tests, they need to be run across multiple network nodes (computers, VMs, or
containers). We ship a :code:`docker-compose` configuration to spin up a
scheduler and two workers as separate docker containers.

Requirements
------------

To run the distributed tests, you need to install (a current version of) :code:`docker`
and :code:`docker-compose`. See
`the official docker documentation <https://docs.docker.com/install/>`_ to get started.
To start the test environment, change to the :code:`packaging/docker/` directory and
run the following:

.. code:: shell

   $ bash dist-build.sh
   $ docker-compose down && docker-compose run tests ; docker-compose down

Note that you need to run the above command as a user that can access the docker daemon,
for example by being in the :code:`docker` group. The containers are running as long as
the :code:`docker-compose` command is running, so you can stop them using :code:`Ctrl-C`.

You can then run the distributed tests using (from the same directory):

.. code:: shell

   $ docker-compose run --rm tests

As we are running the tests in a docker container, the environment of the `tests` container
will automatically match the environment of the workers.

After changing LiberTEM code, you need to rebuild and restart the containers. Just cancel the first
:code:`docker-compose up` run using CTRL-C and start from the top.

The rebuild should be faster than the initial build, which is accomplished by careful
use of the layer caching feature of docker. This also means that you may need to update
the :code:`packaging/docker/requirements.txt` file by running the provided
script :code:`update_reqs.sh` when the dependencies of LiberTEM change, or when new
versions of dependencies are released. In CI, this is done automatically.
.. LiberTEM documentation master file, created by
   sphinx-quickstart on Wed May  2 15:54:09 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

LiberTEM - Open Pixelated STEM platform
=======================================

.. include:: ../../README.rst
..   see also: https://muffinresearch.co.uk/selectively-including-parts-readme-rst-in-your-docs/
..   later, maybe :start-after: inclusion-marker-do-not-remove


Documentation
=============

.. toctree::
   :maxdepth: 2
   :caption: User manual
   
   install
   usage
   applications
   formats
   sample_datasets
   api
   udf
   dask
   changelog
   citing
   acknowledgments

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/index
   tips
   packages
   concepts
   performance
   why_python

.. toctree::
   :maxdepth: 2
   :caption: For developers

   dev/setup.rst
   contributing
   debugging
   architecture
   dev/how-io-works
   dev/executors
   gsoc
   authorship

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Citations
---------

.. bibliography:: references-libertem.bib
Debugging
=========

.. testsetup:: *

    import numpy as np
    from libertem import api
    from libertem.executor.inline import InlineJobExecutor

    ctx = api.Context(executor=InlineJobExecutor())
    data = np.random.random((16, 16, 32, 32)).astype(np.float32)
    dataset = ctx.load("memory", data=data, sig_dims=2)
    roi = np.random.choice([True, False], dataset.shape.nav)

There are different parts of LiberTEM which can be debugged with different tools and methods.

Debugging the Web GUI
---------------------

For debugging the GUI, you can use all standard debugging tools for web development. Most useful
in this context are the `Chrome DevTools <https://developer.chrome.com/docs/devtools/>`_
or `Firefox Developer Tools <https://developer.mozilla.org/en-US/docs/Tools>`_, which can be
accessed by pressing F12. You can extend these with additional panels
`for React <https://reactjs.org/blog/2019/08/15/new-react-devtools.html>`_
and `for Redux <https://github.com/reduxjs/redux-devtools>`_.

These tools can be used for inspecting all frontend-related processes, from network traffic
up to rendering behavior. Take note of the :code:`/api/events/` websocket connection, where all
asynchronous notification and results will be transferred.

Note that you should always debug using the development build of the GUI, using :code:`npm start`,
as described in :ref:`the contributing section <building the client>`. Otherwise the debugging
experience may be painful, for example worse error output from react, minified source and
minified component names, ...

Debugging the API server
------------------------

If the API server returns a server error (500), the detailed exception should be logged
in the output of :code:`libertem-server`. You can also try
`enabling the debug mode of tornado <https://www.tornadoweb.org/en/stable/guide/running.html#debug-mode-and-automatic-reloading>`_
(there is currently no command line flag for this, so you need to change
:py:mod:`libertem.web.server` accordingly.)

If an analysis based on the exception alone is inconclusive,
you can try to reproduce the problem using the Python API and follow the instructions below.

.. _`debugging udfs`:

Debugging UDFs or other Python code
-----------------------------------

If you are trying to write a UDF, or debug other Python parts of LiberTEM, you can
instruct LiberTEM to use simple single-threaded execution using the
:class:`~libertem.executor.inline.InlineJobExecutor`.

.. testsetup::

    from libertem.udf.logsum import LogsumUDF

    udf = LogsumUDF()

.. testcode::

   from libertem.executor.inline import InlineJobExecutor
   from libertem import api as lt

   ctx = lt.Context(executor=InlineJobExecutor())

   ctx.run_udf(dataset=dataset, udf=udf)


You can then use all usual debugging facilities, including
`pdb <https://docs.python.org/3.7/library/pdb.html>`_ and
`the %pdb magic of ipython/Jupyter <https://ipython.org/ipython-doc/3/interactive/magics.html#magic-pdb>`_.

The :class:`libertem.executor.inline.InlineJobExecutor` uses a single CPU core
by default. It can be switched to GPU processing to test CuPy-enabled UDFs by
calling :meth:`libertem.common.backend.set_use_cuda` with the device ID to use.
:code:`libertem.common.backend.set_use_cpu(0)` switches back to CPU processing.

.. testsetup::

    from libertem.udf.masks import ApplyMasksUDF

    udf = ApplyMasksUDF(mask_factories=[lambda:np.ones(dataset.shape.sig)])

.. testcode::

   from libertem.executor.inline import InlineJobExecutor
   from libertem import api as lt
   from libertem.utils.devices import detect
   from libertem.common.backend import set_use_cpu, set_use_cuda

   ctx = lt.Context(executor=InlineJobExecutor())

   d = detect()
   if d['cudas'] and d['has_cupy']:
       set_use_cuda(d['cudas'][0])
   ctx.run_udf(dataset=dataset, udf=udf)
   set_use_cpu(0)

If a problem is only reproducible using the default executor, you will have to follow the
`debugging instructions of dask-distributed <https://docs.dask.org/en/latest/debugging.html>`_.
As the API server can't use the synchronous :class:`~libertem.executor.inline.InlineJobExecutor`,
this is also the case when debugging problems that only occur in context of the API server.

Debugging failing test cases
----------------------------

When a test case fails, there are some options to find the root cause:

The :code:`--pdb` command line switch of pytest can be used to automatically
drop you into a PDB prompt in the failing test case, where you will either land
on the failing :code:`assert` statement, or the place in the code where an
exception was raised.

This does not help if the test case only fails in CI. Here, it may be easier to
use logging. Because we call pytest with the :code:`--log-level=DEBUG`
parameter, the failing test case output will have a section containing the
captured logging output.

You can sprinkle the code with `log.debug(...)` calls that output the relevant
variables. In some cases you may also leave the logging statements in the code
even after the problem is fixed, depending on the overhead.
.. _`user-defined functions`:

User-defined functions (UDFs)
=============================

.. testsetup:: *

   import numpy as np
   from libertem import api
   from libertem.executor.inline import InlineJobExecutor

   from libertem.udf import UDF as RealUDF

   # We override UDF in such a way that it can be used
   # without implementing all methods
   class UDF(RealUDF):
      def get_result_buffers(self):
         return {}

      def process_frame(self, frame):
         pass

   class YourUDF(RealUDF):
      def get_result_buffers(self):
         return {'buf1': self.buffer(kind="nav")}

      def process_frame(self, frame):
         self.results.buf1[:] = 42

   ctx = api.Context(executor=InlineJobExecutor())
   data = np.random.random((16, 16, 32, 32)).astype(np.float32)
   dataset = ctx.load("memory", data=data, sig_dims=2)
   roi = np.random.choice([True, False], dataset.shape.nav)
   udf = YourUDF()

A common case for analyzing big EM data sets is running a reduction operation
on each individual detector frame or other small subsets of a data set and then
combining the results of these reductions to form the complete result. This should
cover a wide range of use cases, from simple mathematical operations, for
example statistics, to complex image processing and analysis, like feature extraction.

The user-defined functions (UDF) interface of LiberTEM allows you to define and run your
own reduction functions easily, without having to worry about parallelizing,
I/O, or the details of buffer management. This corresponds to
a simplified `MapReduce programming model <https://en.wikipedia.org/wiki/MapReduce>`_,
where the intermediate re-keying and shuffling step is omitted.

LiberTEM ships with some :ref:`utility UDFs <utilify udfs>` that implement
general functionality:

* :ref:`Sum <sum udf>`
* :ref:`Logsum <logsum udf>`
* :ref:`StdDev <stddev udf>`
* :ref:`SumSig <sumsig udf>`
* :ref:`Masks and other linear operations <masks udf>`
* :ref:`Pick <pick udf>`

Also, LiberTEM includes :ref:`ready-to-use application-specific UDFs
<applications>`.

It can be helpful to review :ref:`some general concepts <concepts>` before
reading the following sections.

Getting started
---------------

The easiest way of running a function over your data is using the 
:meth:`~libertem.api.Context.map` method of the LiberTEM API. For example,
to calculate the sum over the last signal axis:

.. testcode:: autoudf

   import functools
   import numpy as np

   result = ctx.map(
      dataset=dataset,
      f=functools.partial(np.sum, axis=-1)
   )
   # access the result as NumPy array:
   np.array(result)
   # or, alternatively:
   result.data

The function specified via the :code:`f` parameter is called for each frame / diffraction pattern.
See :ref:`auto UDF` for more details. This is most suited for simple functions; once you have
parameters or want to re-use some data across function calls, you should create a
:class:`~libertem.udf.base.UDF` subclass instead.

:func:`functools.partial` is a higher-order function that allows to create a new
function by wrapping an existing function and passing additional parameters to
it. In this case, the resulting call to :func:`numpy.sum` within
:code:`ctx.map(...)` is :code:`numpy.sum(frame, axis=-1)`. See
https://docs.python.org/3/library/functools.html#functools.partial for more
details.

Example notebook
----------------

See the following notebook for a demonstration of basic UDF functionality. It
can be downloaded `from our example collection on GitHub
<https://github.com/LiberTEM/LiberTEM/blob/master/examples/Introduction%20to%20UDFs.ipynb>`_.

.. toctree::

   udf/introduction

.. _`how UDFs work`:

How UDFs works
--------------

.. image:: ./images/udf-diagram.svg
   :width: 450
   :alt: Diagram with the data flow of a UDF.

To allow for parallel processing, data is first divided into partitions along the navigation axes,
which are worked on by different worker processes. Then, for each frame of a partition, a
user-defined function :meth:`~libertem.udf.base.UDFFrameMixin.process_frame` is called,
which is free to do any imaginable processing.

As a result of splitting the data set into partitions, the results then need to be merged
back together. This is accomplished by calling the :meth:`~libertem.udf.base.UDF.merge` method
after all frames of a partition are processed.

In pseudocode, data is processed in the following way:

.. code-block:: python

   result = empty
   for partition in get_partitions(dataset):
      partition_result = empty
      for frame, frame_slice in get_frames(partition):
         frame_result = process_frame(frame)
         partition_result[frame_slice] = frame_result
      merge(dest=result, src=partition_result)

In reality, the loop over partitions is run in parallel using multiple worker processes,
potentially :ref:`on multiple computers <architecture>`. The loop over individual frames is
run in the worker processes, and the merge function is run in the main process, accumulating the
results, every time the results for a partition are available. 

In addition to :meth:`~libertem.udf.base.UDFFrameMixin.process_frame`, there are two more methods
available for overriding, to work on larger/different units of data at the same time:
:meth:`~libertem.udf.base.UDFTileMixin.process_tile`
and :meth:`~libertem.udf.base.UDFPartitionMixin.process_partition`. They can be used for optimizing
some operations, and are documented in the :ref:`advanced topics <advanced udf>` section.

More about UDFs
---------------

Now would be a good time to :ref:`read more about implementing UDFs <implement
udf>` and :ref:`advanced UDF functionality <advanced udf>`. The :ref:`general
section on debugging <debugging udfs>` helps with resolving issues. Once you
have your UDF working, you can proceed to :ref:`UDF profiling <udf profiling>`
to gain insights into the efficiency of your UDF.

.. toctree::

   udf/basic
   udf/advanced
   udf/profiling

.. seealso::

   :ref:`udf reference`
      API documentation for UDFs
.. _packages:

Package overview
================

LiberTEM is currently being split into several application-specific sub-packages
to keep dependencies manageable and facilitate using LiberTEM code in other
projects and vice-versa. See also :issue:`261`.

This list provides an overview of the current sub-packages and is constantly updated.

:libertem-blobfinder:
    Since 0.4.0: Correlation-based peak finding and strain mapping.
    https://libertem.github.io/LiberTEM-blobfinder/

:libertem-holo:
    Since 0.6.0: Off-axis electron holography (work in progress)
    https://github.com/LiberTEM/LiberTEM-holo/

:ptychography40:
    Ptychography implementations within the Ptychography 4.0 project.
    https://ptychography-4-0.github.io/ptychography/
.. _`installation`:

Installation
============

.. note::
    LiberTEM can currently be used on Python 3.6, 3.7, 3.8 and >=3.9.3.

    If you would like to install the latest development version, please also
    see :ref:`installing from a git clone`.

The short version
-----------------

.. code-block:: shell

    $ virtualenv -p python3 ~/libertem-venv/
    $ source ~/libertem-venv/bin/activate
    (libertem-venv) $ python -m pip install "libertem[torch]"

    # optional for GPU support
    # See also https://docs.cupy.dev/en/stable/install.html
    (libertem-venv) $ python -m pip install cupy

For details, please read on!

Linux and Mac OS X
------------------

AppImage
~~~~~~~~

On Linux, the easiest method is to use the provided AppImage. Just download the
AppImage file from `our releases page on GitHub
<https://github.com/LiberTEM/LiberTEM/releases>`_, mark it executable and run
the AppImage. See also the `official documentation
<https://docs.appimage.org/user-guide/run-appimages.html>`_. Continue by reading
the :ref:`usage documentation`.

Creating an isolated Python environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To provide an isolated environment for LiberTEM and its dependencies, you can
use virtualenvs or conda environments. This is important if you want to use
different Python applications that may have conflicting dependencies, and it
allows to quickly re-create an environment in case things go sideways.

Using virtualenv
################

You can use `virtualenv <https://virtualenv.pypa.io/>`_ or `venv
<https://docs.python.org/3/tutorial/venv.html>`_ if you have a system-wide
compatible Python installation. For Mac OS X, using `conda`_ is recommended.

To create a new virtualenv for LiberTEM, you can use the following command:

.. code-block:: shell

    $ virtualenv -p python3 ~/libertem-venv/

If multiple Python versions are installed, replace :code:`python3` with 
:code:`python3.6` or a later version.

Replace :code:`~/libertem-venv/` with any path where you would like to create
the venv. You can then activate the virtualenv with

.. code-block:: shell

    $ source ~/libertem-venv/bin/activate

Afterwards, your shell prompt should be prefixed with :code:`(libertem-venv)` to
indicate that the environment is active:

.. code-block:: shell

    (libertem-venv) $

For more information about virtualenv, for example if you are using a shell
without :code:`source`, please `refer to the virtualenv documentation
<https://virtualenv.pypa.io/en/stable/user_guide.html>`_. If you are often
working with virtualenvs, using a convenience wrapper like `virtualenvwrapper
<https://virtualenvwrapper.readthedocs.io/en/latest/>`_ is recommended.

Continue by `installing from PyPI`_.

.. _`conda`:

Using conda
###########

If you are already using conda, or if you don't have a system-wide compatible
Python installation, you can create a conda environment for LiberTEM.

This section assumes that you have `installed anaconda or miniconda
<https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation>`_
and that your installation is working.

You can create a new conda environment to install LiberTEM with the following
command:

.. code-block:: shell

    $ conda create -n libertem python=3.9

To install or later run LiberTEM, activate the environment with the following
command (see also :ref:`install on windows` if applicable):

.. code-block:: shell

    $ conda activate libertem

Afterwards, your shell prompt should be prefixed with :code:`(libertem)` to
indicate that the environment is active:

.. code-block:: shell

    (libertem) $

Now the environment is ready to install LiberTEM.

For more information about conda, see their `documentation about creating and
managing environments
<https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.

.. _`installing from PyPI`:

Installing from PyPI
~~~~~~~~~~~~~~~~~~~~

To install the latest release version, you can use pip. Activate the Python
environment (conda or virtualenv) and install using:

.. code-block:: shell

    (libertem) $ python -m pip install libertem

This should install LiberTEM and its dependencies in the environment. Please
continue by reading the :ref:`usage documentation`.

PyTorch
~~~~~~~

LiberTEM can use `PyTorch <https://pytorch.org/>`_ for processing if it is
available. Otherwise it uses NumPy as a fallback. We've experienced up to 2x
speed-ups with PyTorch compared to a default NumPy installation. For that reason
we recommend `installing PyTorch <https://pytorch.org/>`_. We currently use
PyTorch only on the CPU. Contributions to use GPUs as well are very welcome!

You can let pip install PyTorch automatically by using the torch variant, for
example from PyPI:

.. code-block:: shell

    (libertem) $ python -m pip install "libertem[torch]"

CuPy
~~~~

GPU support is based on `CuPy <https://cupy.dev/>`_. See
https://docs.cupy.dev/en/stable/install.html#installing-cupy for installation of
precompiled binary packages (recommended). :code:`python -m pip install
"libertem[cupy]"` installs CuPy from source, which requires a build chain and
can be time-consuming.

.. versionadded:: 0.6.0

Other extra packages
~~~~~~~~~~~~~~~~~~~~

.. versionchanged:: 0.4.0
    A number of LiberTEM applications are being spun out as sub-packages that
    can be installed separately. See :ref:`packages` for an overview.

The full grid matching routines in :py:mod:`libertem.analysis.fullmatch` depend
on `HDBSCAN <https://hdbscan.readthedocs.io/en/latest/>`_. This is an optional
dependency because of installation issues on some platforms.

Updating
~~~~~~~~

When installed from PyPI via pip, you can update like this:

.. code-block:: shell

    (libertem) $ python -m pip install -U libertem

This should install a new version of LiberTEM and update all requirements that
have changed.

After updating the installation, you can run the updated version by restarting
the libertem-server and afterwards reloading all browser windows that are
running the LiberTEM GUI. In other environments, like Jupyter notebooks, you
need to restart the Python interpreter to make sure the new version is used,
for example by restarting the ipython kernel.

.. _`install on windows`:

Windows
-------

The recommended method to install LiberTEM on Windows is based on `Miniconda 64
bit with a compatible Python version <https://conda.io/miniconda.html>`_.
This installs a Python distribution.

The installation and running of LiberTEM on Windows with the
Anaconda Prompt is very similar to `Using conda`_ on Linux or Mac OS X.

Differences:

* You might have to install pip into your local LiberTEM conda environment to
  make sure that ``pip install`` installs packages into your local environment and
  not into the global Anaconda base environment. This helps to avoid permission
  issues and interference between environments.

.. code-block:: shell

    (libertem) > conda install pip

Docker and Singularity
----------------------

.. versionadded:: 0.9.0

A `Docker image with a LiberTEM installation
<https://hub.docker.com/r/libertem/libertem/tags>`_ is available on
Docker hub. See :ref:`containers` for more details.

Troubleshooting
---------------

If you are having trouble with the installation, please let us know by
either `filing an issue  <https://github.com/liberTEM/LiberTEM/issues>`_
or by asking on `our Gitter channel <https://gitter.im/LiberTEM/Lobby>`_.

Integration and deployment
--------------------------

.. toctree::
    :maxdepth: 2

    deployment/jupyter
    deployment/clustercontainer
Why Python?
===========

Python is well established in the scientific community. It is usable both for
developers, but also for many microscopists.

Its high-level construct allow for fast iteration and prototyping. There is an
extensive ecosystem of packages for scientific computation leveraging existing
native libraries. Python has good interoperability with low-level languages
like C, which means it is well suited as a glue language for existing low-level
routines, without introducing inefficiencies such as copies of large buffers on
each interaction.


Isn't Python slow?
------------------

Yes, it can be slow, but we mostly use Python for setting up the computation, creating buffers,
setting parameters, etc. We use it as a glue language for native parts
(libhdfs3, numpy/OpenBLAS, ...) or use `Numba <https://numba.pydata.org/>`_ for critical
code sections.

See for example this profile, visualized as a flamegraph:

.. image:: ./images/read_from_hdfs_profile.png

The workload is similar to what LiberTEM will later do: it reads files from HDFS
and does matrix multiplication on the data.

Most of the time is spent reading the file (block on the left: `sys_read`) or
actually performing the matrix multiplication (center block: anything containing `dgemm`).
The Python parts are mostly in the narrow (= little time) but high (= deep call stacks)
pillar on the right. The dask scheduler is also visible in the profile, but takes up
less than 2% of the total samples.

Note the `swapper` part on the right: this was a full-system profile, so unrelated
things like `swapper` or `intel_idle` are also included. 


But what about (multicore) scaling?
-----------------------------------

``NumPy`` releases the GIL, so multiple threads can work at the same time. Even if
this were not the case, we can still use the multiprocessing workers of ``dask.distributed``
and scale to multiple cores. See also the :doc:`performance` section.
Changelog
=========

.. testsetup:: *

    from libertem import api
    from libertem.executor.inline import InlineJobExecutor

    ctx = api.Context(executor=InlineJobExecutor())
    dataset = ctx.load("memory", datashape=(16, 16, 16, 16), sig_dims=2)

.. _continuous:

0.9.0.dev0
##########

.. toctree::
  :glob:

  changelog/*/*

.. _`v0-8-0`:

.. _latest:

0.8.0 / 2021-10-04
##################

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5547992.svg
   :target: https://doi.org/10.5281/zenodo.5547992

This release mainly contains improvements of center of mass / first moment
analysis and support for starting the web GUI from JupyterHub or JupyterLab.

New features
------------

* Support for center of mass with annular masks in :meth:`~libertem.api.Context.create_com_analysis`,
  :class:`~libertem.analysis.com.COMAnalysis` and the GUI (:issue:`633`, :pr:`1089`).
* Support in the GUI for specifying rotation of scan against detector and
  flipping the detector y axis (:pr:`1087`, :issue:`31`). Previously this was only
  supported in the Python API.
* Tweaks and instructions for JupyterHub and JupyterLab integration in LiberTEM, see :ref:`jupyter integration` (:pr:`1074`).
  New package `LiberTEM/LiberTEM-jupyter-proxy <https://github.com/LiberTEM/LiberTEM-jupyter-proxy>`_
  for interfacing.
* In the web API, support was added to re-run visualization only, without
  re-running UDFs for an analysis. This allows for almost instant feedback
  for some operations, like changing CoM parameters.
* Added token-based authentication. For now, it is only usable via
  integrations like Jupyter. It will be extended to local/manual usage
  later (:pr:`1074`, :issue:`1097`). Please comment on :issue:`1097` if local/manual use
  would be beneficial for you so that it is prioritized accordingly.
* SEQ dataset: Added support for loading excluded pixels from XML (:issue:`805`, :pr:`1077`).
  See :class:`~libertem.io.dataset.seq.SEQDataSet` for more information. Also
  support both :code:`*.seq.seq` and :code:`*.seq` as extension for the main SEQ file
  to find files with matching base name that contain correction data (:issue:`1120`, :pr:`1121`).

Bugfixes
--------

* Assert that the :code:`files` argument to :class:`~libertem.io.dataset.dm.DMDataSet` is actually a list or tuple,
  to prevent iterating over a string path (:pr:`1058`).
* Escape globs to support special characters in file names for multi-file
  datasets (:issue:`1066`, :pr:`1067`).
* Make sure multithreading in the main process still works properly after
  launching a :class:`~libertem.api.Context` (:issue:`1053`, :pr:`1100`).
* Allow custom plots to return RGB as plot data, for example a color
  wheel for vector fields (:issue:`1052`, :pr:`1101`).
* Adjust partition count to match the number of CPU compute workers,
  not total workers to prevent residual partitions (:issue:`1086`, :pr:`1103`).
* Correct partition shape for ROI in :class:`~libertem.udf.base.UDFMeta` (:pr:`1109`).
* Fix memory leak: Don't submit dynamically generated callables directly to the distributed cluster,
  as they are cached in an unbounded cache (:issue:`894,964`, :pr:`1119`).

Documentation
-------------

* Note on handling HDF5 files with non-standard compression
  in :class:`~libertem.io.dataset.hdf5.H5DataSet` (:pr:`1059`).
* Link to two more public datasets: :ref:`hires STO` and :ref:`synthetic STO` (:pr:`1073`).

Misc
----

* Speed up coordinate calculation (:issue:`1108`, :pr:`1109`).
* Make sure tasks are scheduled dynamically on available workers if they have uneven
  run time to benefit more from GPUs (:pr:`1107`).
* Cache loaded libraries to reduce overhead of setting the thread count (:issue:`1117`, :pr:`1118`).

Many thanks to our new contributors Levente Puskás for the excluded pixel loading and to
Matthew Bryan for figuring non-standard compression in HDF5 and improving DM
input validation. Congratulations to Alex for closing the long-standing CoM issue :issue:`31`
and for enabling easy and secure access to the web interface on shared IT infrastructure.

.. _`v0-7-1`:

0.7.1 / 2021-07-08
##################

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5083032.svg
   :target: https://doi.org/10.5281/zenodo.5083032

This is a bugfix release that ensures compatibility with the upcoming numba 0.54
release.

Our custom numba caching makes some assumptions about numba internals, which
have changed in numba 0.54. This fixes compatibility with numba 0.54, and also
makes sure we fail gracefully for future changes (:issue:`1060`, :pr:`1061`).

.. _`v0-7-0`:

0.7.0 / 2021-06-10
##################

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4923277.svg
   :target: https://doi.org/10.5281/zenodo.4923277

This release introduces features that are essential for live data processing,
but can be used for offline processing as well: Live plotting, API for bundled
execution of several UDFs in one run, iteration over partial UDF results, and
asynchronous UDF execution. Features and infrastructure that are specific to
live processing are included in the `LiberTEM-live
<https://github.com/LiberTEM/LiberTEM-live/>`_ package, which will be released
soon.

New features
------------

* Support for postprocessing of results on the main node after merging partial
  results. This adds :meth:`~libertem.udf.base.UDF.get_results` and the
  :code:`use` parameter to :meth:`~libertem.udf.base.UDF.buffer`. See :ref:`udf
  final post processing` for details (:pr:`994`, :pr:`1003`, :issue:`1001`).
* Obtain partial results from each merge step iteratively as a generator
  using :meth:`~libertem.api.Context.run_udf_iter`. See :ref:`partial` and an
  `example
  <https://github.com/LiberTEM/LiberTEM/blob/master/examples/async.ipynb>`_ for
  details (:pr:`1011`)!
* Run multiple UDFs in one pass over a single :code:`DataSet` by passing a
  list of UDFs instead of one UDF in :meth:`~libertem.api.Context.run_udf` and
  :meth:`~libertem.api.Context.run_udf_iter` (:pr:`1011`).
* Allow usage from an asynchronous context with the new :code:`sync=False`
  argument to :meth:`~libertem.api.Context.run_udf` and
  :meth:`~libertem.api.Context.run_udf_iter`. See :ref:`partial` and an `example
  <https://github.com/LiberTEM/LiberTEM/blob/master/examples/async.ipynb>`_ for
  details (:issue:`216`, :pr:`1011`)!
* Live plotting using the new :code:`plots` parameter for
  :meth:`~libertem.api.Context.run_udf` and
  :meth:`~libertem.api.Context.run_udf_iter`, as well as live plotting classes
  documented in :ref:`viz reference`. Pass :code:`plots=True` for simple usage.
  See :ref:`plotting` as well as `an example
  <https://github.com/LiberTEM/LiberTEM/blob/master/examples/live-plotting.ipynb>`_
  for the various possibilities for advanced usage (:issue:`980`, :pr:`1011`).
* Allow some UDF-internal threading. This is mostly
  interesting for ad-hoc parallelization on top of the
  :class:`~libertem.executor.inline.InlineJobExecutor` and live processing that
  currently relies on the :class:`~libertem.executor.inline.InlineJobExecutor`
  for simplicity, but could also be used for hybrid multiprocess/multithreaded
  workloads. Threads for numba, pyfftw, OMP/MKL are automatically
  controlled. The executor makes the number of allowed threads available as
  :attr:`libertem.udf.base.UDFMeta.threads_per_worker` for other threading
  mechanisms that are not controlled automatically (:pr:`993`).
* K2IS: reshaping, sync offset and time series support. Users can now specify a
  :code:`nav_shape`, :code:`sig_shape` and :code:`sync_offset` for a K2IS data
  set, and load time series data (:pr:`1019`, :issue:`911`). Many thanks to
  `@AnandBaburajan <https://github.com/AnandBaburajan>`_ for implementing this
  feature!
* Support for Python >=3.9.3, use Python 3.9 in AppImage (:issue:`914`, :pr:`1037,1039`).

Bugfixes
--------

* UDF: Consistently use attribute access in :code:`UDF.process_*()`, :code:`UDF.merge()`,
  :code:`UDF.get_results()` etc. instead of mixing it with :code:`__getitem__()`
  dict-like access. The previous method still works, but triggers a :class:`UserWarning`
  (:issue:`1000`, :pr:`1003`).
* Also allow non-sliced assignment, for example
  :code:`self.results.res += frame` (:issue:`1000`, :pr:`1003`).
* Better choice of :code:`kind='nav'` buffer fill value outside ROI.

  * String : Was :code:`'n'`, now :code:`''`
  * bool : Was :code:`True`, now :code:`False`
  * integers : Was smallest possible value, now :code:`0`
  * objects : was :code:`np.nan`, now :code:`None` (:pr:`1011`)

* Improve performance for chunked HDF5 files, especially compressed HDF5 files
  which have a chunking in both navigation dimensions. They were causing
  excessive read amplification (:pr:`984`).
* Fix plot range if only zero and one other value are present
  in the result, most notably boolean values (:issue:`944`, :pr:`1011`).
* Fix axes order in COM template: The components in the field are (x, y)
  while the template had them as (y, x) before (:pr:`1023`).

Documentation
-------------

* Update Gatan Digital Micrograph (GMS) examples to work with the current GMS and
  LiberTEM releases and demonstrate the new features. (:issue:`999`,
  :pr:`1002,1004,1011`). Many thanks to Winnie from Gatan for helping to work
  around a number of issues!
* Restructure UDF documentation (:pr:`1034`).
* Document coordinate meta information (:issue:`928`, :pr:`1034`).

Obsolescence
------------

* Removed deprecated blobfinder and :code:`FeatureVecMakerUDF` as
  previously announced. Blobfinder is available as a separate package at
  https://github.com/liberTEM/LiberTEM-blobfinder. Instead of
  :code:`FeatureVecMakerUDF`, you can use a sparse matrix and
  :code:`ApplyMasksUDF` (:pr:`979`).
* Remove deprecated :code:`Job` interface as previously announced.
  The functionality was ported to the more capable UDF interface :pr:`978`.



.. _`v0-6-0`:

0.6.0 / 2021-02-16
##################

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4543704.svg
   :target: https://doi.org/10.5281/zenodo.4543704

We are pleased to announce the latest LiberTEM release, with many improvements
since 0.5. We would like to highlight the contributions of our GSoc 2020
students `@AnandBaburajan <https://github.com/AnandBaburajan>`_ (reshaping and
sync offset correction) and `@twentyse7en <https://github.com/twentyse7en>`_,
(Code generation to replicate GUI analyses in Jupyter notebooks) who implemented
significant improvements in the areas of I/O and the user interface.

Another highlight of this release is experimental support of NVidia GPUs, both
via CuPy and via native libraries. The API is ready to be used, including support
in the GUI. Performance optimization is still to be done (:issue:`946`).
GPU support is activated for all mask-based analyses (virtual detector and
Radial Fourier) for testing purposes, but will not bring a noticeable
improvement of performance yet. GPU-based processing did show significant
benefits for computationally heavy applications like the SSB implementation in
https://github.com/Ptychography-4-0/ptychography.

A lot of work was done to implement tiled reading, resulting in a
new I/O system. This improves performance in many circumstances, especially
when dealing with large detector frames. In addition, a correction module was
integrated into the new I/O system, which can correct gain, subtract a dark
reference, and patch pixel defects on the fly. See below for the full
changelog!

New features
------------

* I/O overhaul

  * Implement tiled reading for most file formats
    (:issue:`27`, :issue:`331`, :issue:`373`, :issue:`435`).
  * Allow UDFs that implement :code:`process_tile` to influence the tile
    shape by overriding :meth:`libertem.udf.base.UDF.get_tiling_preferences`
    and make information about the tiling scheme available to the UDF through
    :attr:`libertem.udf.base.UDFMeta.tiling_scheme`. (:issue:`554`,
    :issue:`247`, :issue:`635`).
  * Update :code:`MemoryDataSet` to allow testing with different
    tile shapes (:issue:`634`).
  * Added I/O backend selection (:pr:`896`), which allows users to select the best-performing
    backend for their circumstance when loading via the new :code:`io_backend`
    parameter of :code:`Context.load`. This fixes a K2IS performance regression
    (:issue:`814`) by disabling any readahead hints by default. Additionaly, this fixes
    a performance regression (:issue:`838`) on slower media (like HDDs), by
    adding a buffered reading backend that tries its best to linearize I/O per-worker.
    GUI integration of backend selection is to be done.
  * For now, direct I/O is no longer supported, please let us know if this is an
    important use-case for you (:issue:`716`)!

* Support for specifying logging level from CLI (:pr:`758`).
* Support for Norpix SEQ files (:issue:`153`, :pr:`767`).
* Support for MRC files, as supported by ncempy (:issue:`152`, :pr:`873`).
* Support for loading stacks of 3D DM files (:pr:`877`). GUI integration still to be done.
* GUI: Filebrowser improvements: users can star directories in the file browser for easy navigation (:pr:`772`).
* Support for running multiple UDFs "at the same time", not yet exposed in public APIs (:pr:`788`).
* GUI: Users can add or remove scan size dimensions according to the dataset's shape (:pr:`779`).
* GUI: Shutdown button to stop server, useful for example for JupyterHub integration (:pr:`786`).
* Infrastructure for consistent coordinate transforms are added in
  :mod:`libertem.corrections.coordinates` and :mod:`libertem.utils`. See also a
  description of coordinate systems in :ref:`concepts`.
* :meth:`~libertem.api.Context.create_com_analysis` now allows to specify a :code:`flipped y axis`
  and a scan rotation angle to deal with MIB files and scan rotation correctly. (:issue:`325`, :pr:`786`).
* Corrections can now be specified by the user when running a UDF (:pr:`778,831,939`).
* Support for loading dark frame and gain map that are sometimes shipped with SEQ data sets.
* GPU support: process data on CPUs, CUDA devices or both (:pr:`760`, :ref:`udf cuda`).
* Spinning out holography to a separate package is in progress: https://github.com/LiberTEM/LiberTEM-holo/
* Implement CuPy support in :class:`~libertem.udf.holography.HoloReconstructUDF`, currently deactivated due to :issue:`815` (:pr:`760`).
* GUI: Allows the user to select the GPUs to use when creating a new local cluster (:pr:`812`).
* GUI: Support to download Jupyter notebook corresponding to an analysis
  made by a user in GUI (:pr:`801`).
* GUI: Copy the Jupyter notebook cells corresponding to the
  analysis directly from GUI, including cluster connection details (:pr:`862`, :pr:`863`)
* Allow reshaping datasets into a custom shape. The :code:`DataSet` implementations (currently except HDF5 and K2IS)
  and GUI now allow specifying :code:`nav_shape` and :code:`sig_shape`
  parameters to set a different shape than the layout in the
  dataset (:issue:`441`, :pr:`793`).
* All :code:`DataSet` implementations handle missing data
  gracefully (:issue:`256`, :pr:`793`).
* The :code:`DataSet` implementations (except HDF5 and K2IS)
  and GUI now allow specifying a :code:`sync_offset` to
  handle synchronization/acquisition problems (:pr:`793`).
* Users can access the coordinates of a tile/partition slice
  through :attr:`~libertem.udf.base.UDFMeta.coordinates` (:issue:`553`, :pr:`793`).
* Cache warmup when opening a data set: Precompiles jit-ed functions on a single process per node, in a controlled manner,
  preventing CPU oversubscription. This improves further through implementing caching for functions which capture other functions
  in their closure (:pr:`886`, :issue:`798`).
* Allow selecting lin and log scaled visualization for sum, stddev, pick and single mask analyses 
  to handle data with large dynamic range. This adds key :code:`intensity_lin` to
  :class:`~libertem.analysis.sum.SumResultSet`, :class:`~libertem.analysis.sum.PickResultSet`
  and the result of :class:`~libertem.analysis.sd.SDAnalysis`.
  It adds key :code:`intensity_log` to :class:`~libertem.analysis.sum.SingleMaskResultSet`.
  The new keys are chosen to not affect existing keys
  (:issue:`925`, :pr:`929`).
* Tuples can be added directly to :code:`Shape` objects. Right
  addition adds to the signal dimensions of the :code:`Shape`
  object while left addition adds to the navigation
  dimensions (:pr:`749`)

Bugfixes
--------

* Fix an off-by-one error in sync offset for K2IS data (drive-by change in :pr:`706`).
* Missing-directory error isn't thrown if it's due to last-recent-directory not being available (:pr:`748`).
* GUI: when cluster connection fails, reopen form with parameters user submitted (:pr:`735`).
* GUI: Fixed the glitch in file opening dialogue by disallowing parallel browsing before loading is concluded (:pr:`752`).
* Handle empty ROI and extra_shape with zero. Empty result buffers of the appropriate shape are returned if the ROI
  is empty or :code:`extra_shape` has a zero (:pr:`765`)
* Improve internals of :mod:`libertem.corrections.detector` and
  :mod:`libertem.corrections.corrset` to better support correction
  of many dead pixels. (:pr:`890`, :issue:`889`)
* Handle single-frame partitions in combination with aux data.
  Instead of squeezing the aux buffer, reshape to the correct shape (:issue:`791`, :pr:`902`).
* Libertem-server can now be started from Bash on Windows (:pr:`731`)
* Fix reading without a copy from multi-file datasets. The start offset of the file was
  not taken account when indexing into the memory maps (:issue:`903`).
* Improve performance and reduce memory consumption of point analysis.
  Custom right hand side matrix product to reduce memory consumption and
  improve performance of sparse masks, such as point analysis. See also
  `scipy/13211 <https://github.com/scipy/scipy/issues/13211>`_ (:issue:`917`, :pr:`920`). 
* Fix stability issue with multiple dask clients. :code:`dd.as_completed` needs
  to specify the :code:`loop` to work with multiple :code:`dask.distributed` clients (:pr:`921`).
* GUI: Snap to pixels in point selection analysis. Consistency between point
  selection and picking (:issue:`926`, :pr:`927`).
* Open datasets with autodetection, positional and keyword arguments.
  Handle keyword and positional arguments to :code:`Context.load('auto', ...)`
  correctly (:issue:`936`, :pr:`938`).

Documentation
-------------

* Switched to the readthedocs sphinx theme, improving the overall
  documentation structure. The developer documentation is now in
  a separate section from the user documentation.

Misc
----

* Command line options can also be accessed with shorter alternatives (:pr:`757`).
* Depend on Numba >= 0.49.1 to support setting Numba thread count (:pr:`783`), bumped to 0.51
  to support caching improvements (:pr:`886`).
* libertem-server: Ask for confirmation if the user press ctrl+c. Can immediately stop using
  another ctrl+c (:pr:`781`).
* Included `pytest-benchmark <https://pytest-benchmark.readthedocs.io/en/latest/usage.html>`_
  to integrate benchmarks in the test infrastructure. See :ref:`benchmarking` for details (:pr:`819`).
* The X and Y components for the color wheel visualization in Center of
  Mass and Radial Fourier Analysis are swapped to match the axis convention in
  empyre. This just changes the color encoding in the visualization and not the
  result (:pr:`851`).

Deprecations
------------

* The :code:`tileshape` parameter of :code:`DataSet` implementations is deprecated in
  favor of tileshape negotiation and will be ignored, if given (:issue:`754`, :pr:`777`).
* Remove color wheel code from :code:`libertem.viz` and replace with imports from empyre.
  Note that these functions expect three vector components instead of two (:pr:`851`).
* The new and consistent :code:`nav_shape` and :code:`sig_shape` parameters should be used
  when loading data. The old :code:`scan_size` and :code:`detector_size` parameters,
  where they existed, are still recognized (:pr:`793`).

.. _`v0-5-1`:

0.5.1 / 2020-08-12
##################

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3982290.svg
   :target: https://doi.org/10.5281/zenodo.3982290

Bugfixes
--------

* Allow installation with latest dask distributed on Python 3.6 and 3.7

.. _`v0-5-0`:

0.5.0 / 2020-04-23
##################

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3763313.svg
   :target: https://doi.org/10.5281/zenodo.3763313

New features
------------

* In addition to tuples, :class:`~libertem.common.shape.Shape` objects can be used as
  :code:`extra_shape` parameter for :meth:`libertem.udf.base.UDF.buffer` and
  :meth:`libertem.udf.base.UDF.aux_data` now. (:pr:`694`)
* Progress bar support based on :code:`tqdm` that can be enabled by passing
  :code:`progress=True` to :meth:`libertem.api.Context.run_udf`,
  :meth:`libertem.api.Context.run` and :meth:`libertem.api.Context.map`: :ref:`progress bar`. (:pr:`613,670,655`)
* Include explicit support for Direct Electron's DE5 format based on HDF5. (:pr:`704`)
* GUI: Downloadable results as HDF5, NPZ, TIFF, and RAW. See
  :ref:`download results` for details. (:pr:`665`)
* :meth:`libertem.api.Context.load` now automatically detects file
  type and parameters if :code:`filetype="auto"` is passed. (:pr:`610,621,734`)
* Relocatable GUI: Allow LiberTEM to run from different URL prefixes, allowing integration into,
  for example, JupyterLab. (:pr:`697`)
* Run :meth:`~libertem.udf.base.UDFPreprocessMixin.preprocess` also before merge on
  the main node to allocate or initialize buffers, in addition to running on the
  workers (:pr:`624`).
* No need to set thread count environment variables anymore since the thread count
  for OpenBLAS, OpenMP, Intel MKL and pyFFTW is now set on the workers at run-time.
  Numba support will be added as soon as Numba 0.49 is released. (:pr:`685`).

Bugfixes
--------

* A large number of usability improvements (:pr:`622,639,641,642,659,666,690,699,700,704`).
  Thanks and credit to many new contributors from GSoC!
* Fixed the buggy "enable Direct I/O" checkbox of the RAW dataset and
  handle unsupported operating systems gracefully. (:pr:`696,659`)


Documentation
-------------

* Added screenshots and description of ROI
  and stddev features in usage docs (:pr:`669`)
* Improved instructions for installing LiberTEM
  (general: :pr:`664`; for development: :pr:`598`)
* Add information for downloading and generating sample
  datasets: :ref:`sample data`. (:pr:`650,670,707`)

Obsolescence
------------

* Parameters :code:`crop_detector_to` and :code:`detector_size_raw` of
  :class:`libertem.io.dataset.raw.RawFileDataSet` are deprecated and will be removed
  after 0.6.0. Please specify :code:`detector_size` instead or use a specialized DataSet, for example for EMPAD.
* :class:`libertem.udf.feature_vector_maker.FeatureVecMakerUDF` is deprecated
  and will be removed in 0.6.0. Use :class:`~libertem.udf.masks.ApplyMasksUDF`
  with a sparse stack of single pixel masks or a stack generated by
  :meth:`libertem_blobfinder.common.patterns.feature_vector` instead.
  (:pr:`618`)

Misc
----

* Clustering analysis

  + Use a connectivity matrix to only cluster neighboring pixels,
    reducing memory footprint while improving speed and quality (:pr:`618`).
  + Use faster :class:`~libertem.udf.masks.ApplyMasksUDF` to generate feature
    vector (:pr:`618`).

* :class:`~libertem.udf.stddev.StdDevUDF`

  + About 10x speed-up for large frames (:pr:`625,640`)
  + Rename result buffers of :class:`~libertem.udf.stddev.StdDevUDF`,
    :meth:`~libertem.udf.stddev.run_stddev` and
    :meth:`~libertem.udf.stddev.consolidate_result` from :code:`'sum_frame'` to
    :code:`'sum'`, :code:`'num_frame'` to :code:`'num_frames'` (:pr:`640`)
  + Resolve ambiguity between variance and sum of variances in result buffer names of
    :class:`~libertem.udf.stddev.StdDevUDF`,
    :meth:`~libertem.udf.stddev.run_stddev` and
    :meth:`~libertem.udf.stddev.consolidate_result`. (:pr:`640`)

* LiberTEM works with Python 3.8 for experimental use. A context using a remote Dask.Distributed cluster
  can lead to lock-ups or errors with Python 3.8. The default local Dask.Distributed context works.
* Improve performance with large tiles. (:pr:`649`)
* :class:`~libertem.udf.sum.SumUDF` moved to the :mod:`libertem.udf` folder (:pr:`613`).
* Make sure the signal dimension of result buffer slices can be
  flattened without creating an implicit copy (:pr:`738`, :issue:`739`)

Many thanks to the contributors to this release: :user:`AnandBaburajan`,
:user:`twentyse7en`, :user:`sayandip18`, :user:`bdalevin`, :user:`saisunku`,
:user:`Iamshankhadeep`, :user:`abiB27`, :user:`sk1p`, :user:`uellue`

.. _`v0-4-1`:

0.4.1 / 2020-02-18
##################

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3674003.svg
   :target: https://doi.org/10.5281/zenodo.3674003

This is a bugfix release, mainly constraining the :code:`msgpack` dependency,
as distributed is not compatible to version 1.0 yet. It also contains
important fixes in the HDF5 dataset.

Bugfixes
--------

* Fix HDF5 with automatic tileshape (:pr:`608`)
* Fix reading from HDF5 with roi beyond the first partition (:pr:`606`)
* Add version constraint on msgpack

.. _`v0-4-0`:

0.4.0 / 2020-02-13
##################

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3666686.svg
   :target: https://doi.org/10.5281/zenodo.3666686

The main points of this release are the :ref:`job deprecation` and restructuring
of our packaging, namely :ref:`extracting the blobfinder module <restructuring-0-4>`.

New features
------------

* :code:`dtype` support for UDFs :ref:`udf dtype` (:issue:`549`, :pr:`550`)
* Dismiss error messages via keyboard: allows pressing the escape key to close all currently open error messages (:issue:`437`)
* ROI doesn't have any effect if in pick mode, so we hide the dropdown in that case (:issue:`511`)
* Make tileshape parameter of HDF5 DataSet optional (:pr:`578`)
* Open browser after starting the server. Enabled by default, can be disabled using --no-browser (:issue:`81`, :pr:`580`)
* Implement :class:`libertem.udf.masks.ApplyMasksUDF` as a replacement of ApplyMasksJob (:issue:`549`, :pr:`550`)
* Implement :class:`libertem.udf.raw.PickUDF` as a replacement of PickFrameJob (:issue:`549`, :pr:`550`)
 
Bug fixes
---------

* Fix FRMS6 in a distributed setting. We now make sure to only do I/O in methods that are running on worker nodes (:pr:`531`).
* Fixed loading of nD HDF5 files. Previously the HDF5 DataSet was hardcoded for
  4D data - now, arbitraty dimensions should be supported (:issue:`574`, :pr:`567`)
* Fix :code:`DaskJobExecutor.run_each_host`. Need to pass :code:`pure=False` to ensure multiple runs of the function (:pr:`528`).

Obsolescence
------------

* Because HDFS support is right now not tested (and to my knowledge also not
  used) and the upstream :code:`hdfs3` project is not actively maintained, remove
  support for HDFS. :code:`ClusterDataSet` or :code:`CachedDataSet` should be used
  instead (:issue:`38`, :pr:`534`).

Misc
----

* Depend on distributed>=2.2.0 because of an API change. (:pr:`577`)
* All analyses ported from Job to UDF back-end. The Job-related code remains for now for comparison purposes (:issue:`549`, :pr:`550`)

.. _`job deprecation`:

Job API deprecation
-------------------

The original Job API of LiberTEM is superseded by the new :ref:`user-defined
functions` API with release 0.4.0. See :issue:`549` for a detailed overview
of the changes. The UDF API brings the following advantages:

* Support for regions of interest (ROIs).
* Easier to implement, extend and re-use UDFs compared to Jobs.
* Clean separation between back-end implementation details and application-specific code.
* Facilities to implement non-trivial operations, see :ref:`advanced udf`.
* Performance is at least on par.

For that reason, the Job API has become obsolete. The existing public
interfaces, namely :meth:`libertem.api.Context.create_mask_job` and
:meth:`libertem.api.Context.create_pick_job`, will be supported in LiberTEM for
two more releases after 0.4.0, i.e. including 0.6.0. Using the Job API will
trigger deprecation warnings starting with this release. The new
:class:`~libertem.udf.masks.ApplyMasksUDF` replaces
:class:`~libertem.job.masks.ApplyMasksJob`, and :class:`~libertem.udf.raw.PickUDF`
replaces :class:`~libertem.job.raw.PickFrameJob`.

The Analysis classes that relied on the Job API as a back-end are already ported
to the corresponding UDF back-end. The new back-end may lead to minor
differences in behavior, such as a change of returned dtype. The legacy code for
using a Job back-end will remain until 0.6.0 and can be activated during the
transition period by setting :code:`analysis.TYPE = 'JOB'` before running.

From :class:`~libertem.job.masks.ApplyMasksJob` to :class:`~libertem.udf.masks.ApplyMasksUDF`
.............................................................................................

Main differences:

* :class:`~libertem.udf.masks.ApplyMasksUDF` returns the result with the first
  axes being the dataset's navigation axes. The last dimension is the mask
  index. :class:`~libertem.job.masks.ApplyMasksJob` used to return transposed
  data with flattened navigation dimension.
* Like all UDFs, running an :class:`~libertem.udf.masks.ApplyMasksUDF` returns a
  dictionary. The result data is accessible with key :code:`'intensity'` as a
  :class:`~libertem.common.buffers.BufferWrapper` object.
* ROIs are supported now, like in all UDFs.

.. testsetup:: jobdeprecation

    import numpy as np
    import libertem
    import matplotlib.pyplot as plt

    def all_ones():
        return np.ones((16, 16))

    def single_pixel():
        buf = np.zeros((16, 16))
        buf[7, 7] = 1
        return buf

Previously with :class:`~libertem.job.masks.ApplyMasksJob`:

.. code-block:: python

    # Deprecated!
    mask_job = ctx.create_mask_job(
      factories=[all_ones, single_pixel],
      dataset=dataset
    )
    mask_job_result = ctx.run(mask_job)

    plt.imshow(mask_job_result[0].reshape(dataset.shape.nav))

Now with :class:`~libertem.udf.masks.ApplyMasksUDF`:

.. testcode:: jobdeprecation

    mask_udf = libertem.udf.masks.ApplyMasksUDF(
      mask_factories=[all_ones, single_pixel]
    )
    mask_udf_result = ctx.run_udf(dataset=dataset, udf=mask_udf)

    plt.imshow(mask_udf_result['intensity'].data[..., 0])

From :class:`~libertem.job.raw.PickFrameJob` to :class:`~libertem.udf.raw.PickUDF`
..................................................................................

:class:`~libertem.job.raw.PickFrameJob` allowed to pick arbitrary contiguous
slices in both navigation and signal dimension. In practice, however, it was
mostly used to extract single complete frames.
:class:`~libertem.udf.raw.PickUDF` allows to pick the *complete* signal
dimension from an arbitrary non-contiguous region of interest in navigation
space by specifying a ROI.

If necessary, more complex subsets of a dataset can be extracted by constructing
a suitable subset of an identity matrix for the signal dimension and using it
with ApplyMasksUDF and the appropriate ROI for the navigation dimension.
Alternatively, it is now easily possible to implement a custom UDF for this
purpose. Performing the complete processing through an UDF on the worker nodes
instead of loading the data to the central node may be a viable alternative as
well.

:class:`~libertem.udf.raw.PickUDF` now returns data in the native :code:`dtype`
of the dataset. Previously, :class:`~libertem.job.raw.PickFrameJob` converted to
floats.

Using :meth:`libertem.api.Context.create_pick_analysis` continues to be the
recommended convenience function to pick single frames.

.. _`restructuring-0-4`:

Restructuring into sub-packages
-------------------------------

We are currently restructuring LiberTEM into packages that can be installed and
used independently, see :issue:`261`. This will be a longer process and changes
the import locations.

* `Blobfinder <https://libertem.github.io/LiberTEM-blobfinder/>`_ is the first
  module separated in 0.4.0.
* See :ref:`packages` for a current overview of sub-packages.

For a transition period, importing from the previous locations is supported but
will trigger a :code:`FutureWarning`. See :ref:`show warnings` on how to
activate deprecation warning messages, which is strongly recommended while the
restructuring is ongoing.

.. _`v0-3-0`:

0.3.0 / 2019-12-12
##################

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3572855.svg
   :target: https://doi.org/10.5281/zenodo.3572855

New features
------------

* Make OOP based composition and subclassing easier for
  :class:`~libertem.udf.blobfinder.correlation.CorrelationUDF` (:pr:`466`)
* Introduce plain circular match pattern :class:`~libertem.udf.blobfinder.patterns.Circular` (:pr:`469`)
* Distributed sharded dataset :class:`~libertem.io.dataset.cluster.ClusterDataSet` (:issue:`136`, :issue:`457`)
* Support for caching data sets :class:`~libertem.io.dataset.cached.CachedDataSet`
  from slower storage (NFS, spinning metal) on fast local storage (:pr:`471`)
* :ref:`Clustering` analysis (:pr:`401,408` by :user:`kruzaeva`).
* :class:`libertem.io.dataset.dm.DMDataSet` implementation based on ncempy (:pr:`497`)

  * Adds a new :meth:`~libertem.executor.base.JobExecutor.map` executor primitive. Used to concurrently
    read the metadata for DM3/DM4 files on initialization.
  * Note: no support for the web GUI yet, as the naming patterns for DM file series varies wildly. Needs
    changes in the file dialog.

* Speed up of up to 150x for correlation-based peak refinement in
  :mod:`libertem.udf.blobfinder.correlation` with a Numba-based pipeline (:pr:`468`)
* Introduce :class:`~libertem.udf.blobfinder.correlation.FullFrameCorrelationUDF` which
  correlates a large number (several hundred) of small peaks (10x10) on small
  frames (256x256) faster than
  :class:`~libertem.udf.blobfinder.correlation.FastCorrelationUDF` and
  :class:`~libertem.udf.blobfinder.correlation.SparseCorrelationUDF` (:pr:`468`)
* Introduce :class:`~libertem.udf.UDFPreprocessMixin` (:pr:`464`)
* Implement iterator over :class:`~libertem.analysis.base.AnalysisResultSet` (:pr:`496`)
* Add hologram simulation
  :func:`libertem.utils.generate.hologram_frame` (:pr:`475`)
* Implement Hologram reconstruction UDF
  :class:`libertem.udf.holography.HoloReconstructUDF` (:pr:`475`)

Bug fixes
---------

* Improved error and validation handling when opening files with GUI (:issue:`433,442`)
* Clean-up and improvements of :class:`libertem.analysis.fullmatch.FullMatcher` (:pr:`463`)
* Ensure that RAW dataset sizes are calculated as int64 to avoid integer overflows (:pr:`495`, :issue:`493`)
* Resolve shape mismatch issue and simplify dominant order calculation in Radial Fourier Analysis (:pr:`502`)
* Actually pass the :code:`enable_direct` parameter from web API to the DataSet

Documentation
-------------

* Created :ref:`authorship` (:pr:`460,483`)
* Change management process (:issue:`443`, :pr:`451,453`)
* Documentation for :ref:`crystallinity map` and :ref:`clustering` analysis (:pr:`408` by :user:`kruzaeva`)
* Instructions for profiling slow tests (:issue:`447`, :pr:`448`)
* Improve API reference on Analysis results (:issue:`494`, :pr:`496`)
* Restructure and update the API reference for a number of UDFs and
  other application-specific code (:issue:`503`, :pr:`507,508`)

Obsolescence
------------

* The Job interface is planned to be replaced with an implementation based on UDFs in one of the upcoming releases.

Misc
----

* Split up the blobfinder code between several files to reduce file size (:pr:`468`)

.. _`v0-2-2`:

0.2.2 / 2019-10-14
##################

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3489385.svg
   :target: https://doi.org/10.5281/zenodo.3489385

Point release to fix a number of minor issues, most notably PR :pr:`439` that
should have been merged for version 0.2.

Bug fixes
---------

* Trigger a timeout when guessing parameters for HDF5 takes too long (:issue:`440` , :pr:`449`)
* Slightly improved error and validation handling when opening files with GUI (:commit:`ec74c1346d93eff58d9e2201a7ead5af7aa7cf44`)
* Recognize BLO file type (:issue:`432`)
* Fixed a glitch where negative peak elevations were possible (:pr:`446`)
* Update examples to match 0.2 release (:pr:`439`)

.. _`v0-2-1`:

0.2.1 / 2019-10-07
##################

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3474968.svg
   :target: https://doi.org/10.5281/zenodo.3474968

Point release to fix a bug in the Zenodo upload for production releases.

.. _`v0-2-0`:

0.2.0 / 2019-10-07
##################

This release constitutes a major update after almost a year of development.
Systematic change management starts with this release.

This is the `release message <https://groups.google.com/d/msg/libertem/p7MVoVqXOs0/vP_tu6K7CwAJ>`_: 

User-defined functions
----------------------

LiberTEM 0.2 offers a new API to define a wide range of user-defined reduction
functions (UDFs) on distributed data. The interface and implementation offers a
number of unique features:

* Reductions are defined as functions that are executed on subsets of the data.
  That means they are equally suitable for distributed computing, for interactive
  display of results from a progressing calculation, and for handling live data¹.
* Interfaces adapted to both simple and complex use cases: From a simple map()
  functionality to complex multi-stage reductions.
* Rich options to define input and output data for the reduction functions, which
  helps to implement non-trivial operations efficiently within a single pass over
  the input data.
* Composition and extension through object oriented programming
* Interfaces that allow highly efficient processing: locality of reference, cache
  efficiency, memory handling

Introduction: https://libertem.github.io/LiberTEM/udf.html

Advanced features: https://libertem.github.io/LiberTEM/udf/advanced.html

A big shoutout to Alex (:user:`sk1p`) who developed it! 🏆

¹User-defined functions will work on live data without modification as soon as
LiberTEM implements back-end support for live data, expected in 2020.

Support for 4D STEM applications
--------------------------------

In parallel to the UDF interface, we have implemented a number of applications
that make use of the new facilities:

* Correlation-based peak finding and refinement for CBED (credit: Karina Ruzaeva :user:`kruzaeva`)
* Strain mapping
* Clustering
* Fluctuation EM
* Radial Fourier Series (advanced Fluctuation EM)

More details and examples: https://libertem.github.io/LiberTEM/applications.html

Extended documentation
----------------------

We have greatly improved the coverage of our documentation:
https://libertem.github.io/LiberTEM/index.html#documentation

Fully automated release pipeline
--------------------------------

Alex (:user:`sk1p`) invested a great deal of effort into fully automating our release
process. From now on, we will be able to release more often, including service
releases. 🚀

Basic dask.distributed array integration
----------------------------------------

LiberTEM can generate efficient dask.distributed arrays from all supported
dataset types with this release. That means it should be possible to use our high-performance file
readers in applications outside of LiberTEM.

File formats
------------

Support for various file formats has improved. More details:
https://libertem.github.io/LiberTEM/formats.html

.. _`v0-1-0`:

0.1.0 / 2018-11-06
##################

Initial release of a minimum viable product and proof of concept.

Support for applying masks with high throughput on distributed systems with
interactive web GUI display and scripting capability.
.. _`api documentation`:

Python API
==========

The Python API is a concise API for using LiberTEM from Python code. It is suitable both
for interactive scripting, for example from Jupyter notebooks, and for usage
from within a Python application or script.

Basic example
-------------

This is a basic example to load the API, create a local cluster, load a file and
run an analysis. For complete examples on how to use the Python API, please see the
Jupyter notebooks in `the example directory
<https://github.com/LiberTEM/LiberTEM/tree/master/examples>`_.

For more details, please see :ref:`loading data`, :ref:`dataset api` and
:ref:`format-specific reference<formats>`. See :ref:`sample data` for publicly available datasets.

.. include:: /../../examples/basic.py
    :code:

Custom processing routines
--------------------------

To go beyond the included capabilities of LiberTEM, you can implement your own
using :ref:`user-defined functions`.

Reference
---------

For a full reference, please see :ref:`reference`.

.. _`executors`:

Executors
---------

.. versionadded:: 0.9.0
    The executor API is internal. Since choice and parameters of executors
    are important for integration with Dask and other frameworks,
    they are now documented. Only the names and creation methods for
    executors are reasonably stable. The rest of the API is subject to
    change without notice. For that reason it is documented in the developer
    section and not in the API reference.

The default executor is :class:`~libertem.executor.dask.DaskJobExecutor` that
uses the dask.distributed scheduler. To support all LiberTEM features and
achieve optimal performance, the methods provided by LiberTEM to start a
dask.distributed cluster should be used. However, LiberTEM can also run on a
"vanilla" dask.distributed cluster. Please note that dask.distributed clusters
that are not created by LiberTEM might use threading or a mixture of threads and
processes.

The :class:`~libertem.executor.inline.InlineJobExecutor` runs all tasks
synchronously in the current thread. This is useful for debugging and for
special applications such as running UDFs that perform their own multithreading
efficiently or for other non-standard use that requires tasks to be executed
sequentially and in order.

See also :ref:`threading` for more information on multithreading in UDFs.

.. versionadded:: 0.9.0

The :class:`~libertem.executor.concurrent.ConcurrentJobExecutor` runs all tasks
using :mod:`python:concurrent.futures`. Using the
:class:`python:concurrent.futures.ThreadPoolExecutor` allows sharing large
amounts of data as well as other resources between main thread and workers
efficiently, but is severely slowed down by the Python `global interpreter lock
<https://wiki.python.org/moin/GlobalInterpreterLock>`_ under many circumstances.
Furthermore, it can create thread safety issues such as
https://github.com/LiberTEM/LiberTEM-blobfinder/issues/35.

For special applications, the
:class:`~libertem.executor.delayed.DelayedJobExecutor` can use `dask.delayed
<https://docs.dask.org/en/stable/delayed.html>`_ to delay the processing. This
is experimental, see :ref:`dask` for more details. It might use threading as
well, depending on the Dask scheduler that is used by :code:`compute()`.

Common executor choices
.......................

.. versionadded:: 0.9.0

:meth:`libertem.api.Context.make_with` provides a convenient shortcut to start a
:class:`~libertem.api.Context` with common executor choices. See the
:meth:`API documentation <libertem.api.Context.make_with>`
for available options!

Connect to a cluster
....................

See :ref:`cluster` on how to start a scheduler and workers.

.. Not run with docs-check since it doesn't play well with launching a multiprocessing
   cluster

.. code-block:: python

    from libertem import api
    from libertem.executor.dask import DaskJobExecutor

    # Connect to a Dask.Distributed scheduler at 'tcp://localhost:8786'
    with DaskJobExecutor.connect('tcp://localhost:8786') as executor:
        ctx = api.Context(executor=executor)
        ...

.. _`cluster spec`:

Customize CPUs and CUDA devices
...............................

To control how many CPUs and which CUDA devices are used, you can specify them as follows:

.. Not run with docs-check since it doesn't play well with launching a multiprocessing
   cluster

.. code-block:: python

    from libertem import api
    from libertem.executor.dask import DaskJobExecutor, cluster_spec
    from libertem.utils.devices import detect

    # Find out what would be used, if you like
    # returns dictionary with keys "cpus" and "cudas", each with a list of device ids
    devices = detect()

    # Example: Deactivate CUDA devices by removing them from the device dictionary
    devices['cudas'] = []

    # Example: Deactivate CuPy integration
    devices['has_cupy'] = False

    # Example: Use 3 CPUs. The IDs are ignored at the moment, i.e. no CPU pinning
    devices['cpus'] = range(3)

    # Generate a spec for a Dask.distributed SpecCluster
    # Relevant kwargs match the dictionary entries
    spec = cluster_spec(**devices)
    # Start a local cluster with the custom spec
    with DaskJobExecutor.make_local(spec=spec) as executor:
        ctx = api.Context(executor=executor)
        ...

Please see :ref:`dask executor` for a reference of the Dask-based executor.
[Docs] Info on multithreading
=============================

* Information on multithreading added to UDF docs in :ref:`threading` (:pr:`1170`).
[Feature] Automated parameter guess for Center of Mass
======================================================

* Experimental helper function :meth:`libertem.analysis.com.guess_corrections`
  to guess parameters for Center of Mass analysis (:pr:`1111`).
* Add a button in the GUI to call :meth:`libertem.analysis.com.guess_corrections`
  and update the GUI parameters from the result (:pr:`1172`).
[Feature] DaskDataSet
=====================
* DataSet class enabling LiberTEM to process a Dask array through the
  UDF interface, using the array chunking for the partition structure.
  :class:`libertem.io.dataset.dask.DaskDataSet` (:pr:`1137`).Add support for some RAW MIB Quad formats
=========================================
* For now, we support :code:`1x1` and :code:`2x2` layouts, with 1bit and 6bit counter depth.
  Support for other layouts and bit depths can be added on demand (:pr:`1169`, :issue:`1135`).
Add :code:`--preload` option to :code:`libertem-server` and :code:`libertem-worker`
===================================================================================

* Make it work as documented in `HDF5 docs
  <https://libertem.github.io/LiberTEM/reference/dataset.html#hdf5>`_, follow
  `Dask worker preloading
  <https://docs.dask.org/en/stable/how-to/customize-initialization.html#preload-scripts>`_
  (:pr:`1151`).[Feature] Executors and Dask integration
========================================

* See :ref:`dask` for a more detailed description of Dask integration.
* Add :class:`~libertem.executor.delayed.DelayedJobExecutor`,
  :class:`~libertem.executor.concurrent.ConcurrentJobExecutor` and
  :meth:`libertem.executor.integration.get_dask_integration_executor`.
* Add :meth:`libertem.api.Context.make_with` to make common executor choices more
  accessible to users.
* Allow using an existing Dask Client as well as setting the LiberTEM Client as default
  Dask scheduler for integration purposes.
* Add :meth:`libertem.udf.base.UDFMergeAllMixin.merge_all` API to combine intermediate
  UDF task results to chunked Dask arrays efficiently.
* Restructure and extend documentation of executors.
* :pr:`1170`, :issue:`1146,922`
Add :code:`UDFMeta.sig_slice` and :code:`UDFMeta.tiling_scheme_idx`
===================================================================

* These attributes can be used for performant access to the current signal
  slice - mostly important for throughput-limited analysis (:pr:`1167`, :issue:`1166`).
Re-add support for direct I/O
=============================

* Direct I/O was previously only supported as a special case for raw files,
  now it is supported for all native dataset formats we support - notable
  exceptions are  HDF5, MRC, and SER (:pr:`1129`, :issue:`716`).
Add support for reading TVIPS binary files
==========================================

* Binary :code:`*_NNN.tvips` files can now be opened with LiberTEM (:pr:`1179`).
Fix :code:`make_dask_array` with :code:`roi`
============================================

* :code:`make_dask_array` no longer crashes when a :code:`roi` is specified
  (:issue:`933`).
Fix running CoM analysis on a linescan dataset
==============================================

* CoM analyses implicitly calculate x- and y- gradients of the centre
  positions through np.gradient. This previously failed for nav
  dimensions of length 1. The new behaviour is to only return div/curl
  results in the COMResultSet when they are defined.
  (:issue:`1138`, :issue:`1139`)Improve performance with large parameters
=========================================

* Previously, parameters were baked into the :code:`UDFTask` objects, so they were
  transferred multiple times for a single UDF run. To allow for larger parameters,
  they are now handled separately from the function that is run (:pr:`1143`).
[Misc] Preloading
=================

* Start using :mod:`libertem.preload` again and import :code:`hdf5plugin` if
  present so that users don't have to specify this common selection of HDF5
  filters as preload themselves. (:pr:`1160`)
  [Misc] Docker image available
=============================

* A `Docker image with a LiberTEM installation
  <https://hub.docker.com/r/libertem/libertem/tags>`_
  is available on DockerHub now (:pr:`1144`, :issue:`484`).
.. _`jupyter integration`:

Jupyter integration
===================

.. versionadded:: 0.8.0

The web-based LiberTEM GUI can be integrated into existing JupyterLab and
JupyterHub installations. That way, the existing authentication and user
management infrastructure of Jupyter can be re-used to offer LiberTEM as
a service. The LiberTEM GUI will have access to the same data that is available
to Jupyter notebooks.

.. note::

    Currently, Jupyter integration is only supported on Unix, not on Windows
    because of outstanding issues in some dependencies. See `simpervisor/6
    <https://github.com/jupyterhub/simpervisor/issues/6>`_ and
    `jupyter-server-proxy/181
    <https://github.com/jupyterhub/jupyter-server-proxy/pull/181>`_ for more
    details.

Installation
------------

.. note::

    Currently, a fork of :code:`jupyter-server-proxy` needs to be used, which
    includes a fix to allow restarting LiberTEM after a graceful shutdown, see
    also
    `jupyter-server-proxy/215
    <https://github.com/jupyterhub/jupyter-server-proxy/pull/215>`_. 
    The fork also includes other unreleased but necessary changes from the
    master branch of :code:`jupyter-server-proxy`.
    Installation can be done with the following commands, after having followed
    the installation instructions below:

    .. code-block:: shell
        
        (jupyter-venv) $ git clone git@github.com:sk1p/jupyter-server-proxy.git
        (jupyter-venv) $ cd jupyter-server-proxy && git checkout websocket-fix-plus-restart
        (jupyter-venv) $ pip install .

    Note that this needs to be done in the environment where 
    :code:`libertem-jupyter-proxy` was installed.


As a first step, you will need to install
`LiberTEM-jupyter-proxy <https://github.com/LiberTEM/LiberTEM-jupyter-proxy>`_
into the Python environment where JupyterHub or JupyterLab is installed. Then,
you will also need an installation of LiberTEM, which can live in either the same
Python environment, or in its own.

In case of a separate Python environment for LiberTEM, you will need to tell
LiberTEM-jupyter-proxy how to start the :code:`libertem-server`. This is done
with a simple JSON configuration file. The configuration file needs to be saved in the
Python environment where LiberTEM-jupyter-proxy is installed. If the virtualenv
that contains the Jupyter installation is created at
:code:`/opt/jupyter/venv/`, the configuration filename should be
:code:`/opt/jupyter/venv/etc/libertem_jupyter_proxy.json`. The contents should
look like this:

.. code-block:: json

   {"libertem_server_path": "/opt/libertem/venv/bin/libertem-server"}

This allows to separate the Jupyter installation from the LiberTEM installation.
The given path can also point to a wrapper script, for example for further environment
preparations, like loading modules on an HPC system. The wrapper script should forward
any arguments it reveices to the :code:`libertem-server` call.

Usage
-----

When the installation is finished, you can find the LiberTEM GUI as an entry
in the :code:`New` menu in JupyterHub:

..  figure:: ../images/use/jupyter-hub.png

Or as a launcher icon in JupyterLab:

..  figure:: ../images/use/jupyter-lab.png
Containers and clusters
=======================

.. versionadded:: 0.9.0
    `LiberTEM repository on Docker Hub
    <https://hub.docker.com/r/libertem/libertem/tags>`_ with images
    for public use.

.. note::
    The LiberTEM server will only bind to localhost by default, unless
    token-based authentication is enabled or the :code:`--insecure` flag is
    provided. It should run behind a reverse proxy that supplies the token and
    adds encryption, authentication and authorization when access from an
    untrusted network is desired. :ref:`jupyter integration` can be used for
    that.
    
    Third-party commands like :code:`dask-scheduler` may bind to public
    interfaces by default, which can expose a machine to remote code execution.
    For that reason the interface should always be specified and limited to
    trusted networks.
    
    Furthermore, the default command of the LiberTEM Docker container starts a
    server that binds to all interfaces to facilitate integration. Docker runs
    containers in an isolated environment and requires the user to expose the
    port explicitly. Singularity, however, will run the container like a regular
    user program. That means it exposes an insecure LiberTEM server to all
    interfaces when running the default container command. For that reason the
    command to run should always be specified explicitly when using Singularity.

.. _`containers`:

Containers
----------

Docker images can be found in `the LiberTEM repository on Docker Hub
<https://hub.docker.com/r/libertem/libertem/tags>`_. LiberTEM is
installed in a virtual environment in :code:`/venv/` in the Docker image. The
executables :code:`libertem-server`, :code:`dask-scheduler` and
:code:`libertem-worker` can be found in :code:`/venv/bin/`, consequently. The
container runs :code:`/venv/bin/libertem-server --host 0.0.0.0 --insecure --port 9000`
by default.

When using Docker, you can run and expose the LiberTEM server to
:code:`localhost` while accessing local data like this:

.. code-block:: shell

    $ docker run -p 127.0.0.1:9000:9000 \
      --mount type=bind,source=/path/to/your/data/,dst=/data/,ro libertem/libertem

To use the Docker image and Singularity to start :code:`libertem-server`:

.. code-block:: shell

    $ singularity run docker://libertem/libertem -- /venv/bin/libertem-server

Available versions
..................

The tag "latest" (default) points to the stable release with the highest version
number, and "continuous" points to the current master. Version tags for all
stable releases are available as well. See `the LiberTEM repository on Docker
Hub <https://hub.docker.com/r/libertem/libertem/tags>`_ for details.

Updating
........

You can update to the latest release like this:

.. code-block:: shell

    $ docker pull libertem/libertem

or

.. code-block:: shell

    $ singularity pull docker://libertem/libertem

.. _`cluster`:

Starting a custom cluster
-------------------------

LiberTEM can connect to a running Dask cluster. To start a cluster on
:code:`localhost`, first run a scheduler:

.. code-block:: shell

    (libertem) $ dask-scheduler --host localhost

GPU support in LiberTEM requires specific resource tags and environment settings
on the dask workers. The easiest way to start workers with the appropriate
settings is

.. code-block:: shell

    (libertem) $ libertem-worker tcp://localhost:8786

There are a few command line options available:

.. include:: ../autogenerated/libertem-worker.help
    :literal:

.. versionadded:: 0.6.0
.. versionadded:: 0.9.0.dev0
    :code:`--preload` was added.

For a cluster setup, you can run the scheduler on the appropriate network interface and
run workers on all cluster nodes to connect to the scheduler.

You can then connect to the cluster's scheduler URL in the LiberTEM web GUI.

For easier deployment of in container-based environments, you can also use the
Docker image.

Example: Start a scheduler and workers in an isolated environment with Docker.

.. code-block:: shell

    $ docker run --mount type=bind,source=/path/to/your/data/,dst=/data/,ro \
      libertem/libertem /venv/bin/dask-scheduler
    $ docker run --mount type=bind,source=/path/to/your/data/,dst=/data/,ro \
      libertem/libertem /venv/bin/libertem-worker tcp://<scheduler-addr>:8786

Example: Start a scheduler and workers in the context
of the local user with Singularity.

.. code-block:: shell

    $ singularity run docker://libertem/libertem -- /venv/bin/dask-scheduler --host localhost
    $ singularity run docker://libertem/libertem -- /venv/bin/libertem-worker tcp://localhost:8786
.. _`holography app`:

Off-axis electron holography
============================

LiberTEM has implementations for both :ref:`hologram simulation <holo-sim>` and 
:ref:`hologram reconstruction <holo-reconstruct>` for off-axis electron holography.

.. versionadded:: 0.3.0

.. _holo-sim:

Hologram simulation
-------------------
Holograms can be simulated using the method described by Lichte et al. :cite:`Lichte2008`
The simulator includes simulation of holograms with Gaussian and Poisson noise, without effect of
Fresnel fringes of the biprism. The simulator requires amplitude and phase images being provided. Those can be
calculated as in example below in which for amplitude a sphere is assumed, the same sphere is used
for the mean inner potential (MIP) contribution to the phase and in addition to the quadratic long-range
phase shift originating from the centre of the sphere:

.. testsetup:: *

    from libertem import api
    from libertem.executor.inline import InlineJobExecutor

    ctx = api.Context(executor=InlineJobExecutor())

.. testcode::

   import numpy as np
   import matplotlib.pyplot as plt
   from libertem.utils.generate import hologram_frame

   # Define grid
   sx, sy = (256, 256)
   mx, my = np.meshgrid(np.arange(sx), np.arange(sy))
   # Define sphere region
   sphere = (mx - 33.)**2 + (my - 103.)**2 < 20.**2
   # Calculate long-range contribution to the phase
   phase = ((mx - 33.)**2 + (my - 103.)**2) / sx / 40.
   # Add mean inner potential contribution to the phase
   phase[sphere] += (-((mx[sphere] - 33.)**2 \
                      + (my[sphere] - 103.)**2) / sx / 3 + 0.5) * 2.
   # Calculate amplitude of the phase
   amp = np.ones_like(phase)
   amp[sphere] = ((mx[sphere] - 33.)**2 \
                  + (my[sphere] - 103.)**2) / sx / 3 + 0.5

   # Plot
   f, ax = plt.subplots(1, 2)
   ax[0].imshow(amp, cmap='gray')
   ax[0].title.set_text('Amplitude')
   ax[0].set_axis_off()
   ax[1].imshow(phase, cmap='viridis')
   ax[1].title.set_text('Phase')
   ax[1].set_axis_off()

.. image:: ./images/holography/amplitude_phase.png

To generate the object hologram, :code:`amp` and :code:`phase` should be passed to the :code:`holo_frame`
function as follows:

.. testcode::

   holo = hologram_frame(amp, phase)

To generate the vacuum reference hologram, use an array of ones for amplitude and zero for phase:

.. testcode::

   ref = hologram_frame(np.ones_like(phase), np.zeros_like(phase))

   # Plot
   f, ax = plt.subplots(1, 2)
   ax[0].imshow(holo, cmap='gray')
   ax[0].title.set_text('Object hologram')
   ax[0].set_axis_off()
   ax[1].imshow(ref, cmap='gray')
   ax[1].title.set_text('Reference hologram')
   ax[1].set_axis_off()

.. image:: ./images/holography/holograms.png

.. _holo-reconstruct:

Hologram reconstruction
-----------------------

LiberTEM can be used to reconstruct off-axis electron holograms using the Fourier space method.
The processing involves the following steps:

* Fast Fourier transform
* Filtering of the sideband in Fourier space and cropping (if applicable)
* Centering of the sideband
* Inverse Fourier transform.

The reconstruction can be accessed through the :class:`~libertem.udf.holography.HoloReconstructUDF` class.
To demonstrate the reconstruction capability, two datasets can be created from the holograms
simulated above as follows:

.. testcode::

   from libertem.io.dataset.memory import MemoryDataSet
   from libertem.udf.holography import HoloReconstructUDF

   dataset_holo = MemoryDataSet(data=holo.reshape((1, sx, sy)),
                                tileshape=(1, sx, sy),
                                num_partitions=1, sig_dims=2)
   dataset_ref = MemoryDataSet(data=ref.reshape((1, sx, sy)),
                               tileshape=(1, sx, sy),
                               num_partitions=1, sig_dims=2)

The reconstruction requires knowledge about the position of the sideband and the size of the
sideband filter which will be used in the reconstruction. The position of the sideband can be
estimated from the Fourier transform of the vacuum reference hologram:

.. testcode::

   # Plot FFT and the sideband position
   plt.imshow(np.log(np.abs(np.fft.fft2(ref))))
   plt.plot(26., 44., '+r')
   plt.axis('off')
   plt.title('FFT of the reference hologram')

   # Define position
   sb_position = [44, 26]

.. image:: ./images/holography/FFT_reference.png

The radius of sideband filter is typically chosen as either half of the distance between the sideband and
autocorrelation for strong phase objects or as one third of the distance for weak phase objects. Assuming
a strong phase object, one can proceed as follows:

.. testcode::

   sb_size = np.hypot(sb_position[0], sb_position[1]) / 2.

Since in off-axis electron holography, the spatial resolution is determined by the interference
fringe spacing rather than by the sampling of the original images, the reconstruction would typically
involve changing the shape of the data.
For medium magnification holography the size of the reconstructed images can be typically set to the size
(diameter) of the sideband filter. (For high-resolution holography reconstruction, typically binning factors of
1-4 are used.) Therefore, the output shape can be defined as follows:

.. testcode::

   output_shape = (int(sb_size * 2), int(sb_size * 2))

Finally the :class:`~libertem.udf.holography.HoloReconstructUDF` class can be used to reconstruct the object and
reference holograms:

.. testcode::

   # Create reconstruction UDF:
   holo_udf = HoloReconstructUDF(out_shape=output_shape,
                                 sb_position=sb_position,
                                 sb_size=sb_size)

   # Reconstruct holograms, access data directly
   w_holo = ctx.run_udf(dataset=dataset_holo,
                        udf=holo_udf)['wave'].data
   w_ref = ctx.run_udf(dataset=dataset_ref,
                       udf=holo_udf)['wave'].data

   # Correct object wave using reference wave
   w = w_holo / w_ref

   # Calculate plot phase shift and amplitude
   amp_r = np.abs(w)
   phase_r = np.angle(w)

   # Plot amplitude
   f, ax = plt.subplots(1, 2)
   ax[0].imshow(amp)
   ax[0].title.set_text('Input amplitude')
   ax[0].set_axis_off()
   ax[1].imshow(amp_r[0])
   ax[1].title.set_text('Reconstructed amplitude')
   ax[1].set_axis_off()

.. image:: ./images/holography/amp_comparison.png

One sees that the reconstructed amplitude has artifacts due to digital Fourier processing. Those are typical for
synthetic data. One of the ways to get synthetic data closer to the experimental would be adding noise.
Comparing phase images, one should keep in mind that phase is typically wrapped in an interval :math:`[0; 2\pi)`.
To unwrap phase one can do the following:

.. testcode::

   from skimage.restoration import unwrap_phase

   # Unwrap phase:
   phase_unwrapped = unwrap_phase(phase_r[0])

   # Plot
   f, ax = plt.subplots(1, 3)
   ax[0].imshow(phase, cmap='viridis')
   ax[0].title.set_text('Input phase')
   ax[0].set_axis_off()
   ax[1].imshow(phase_r[0])
   ax[1].title.set_text('Reconstructed phase')
   ax[1].set_axis_off()
   ax[2].imshow(phase_unwrapped, cmap='viridis')
   ax[2].title.set_text('Reconstructed phase (unwrapped)')
   ax[2].set_axis_off()

.. image:: ./images/holography/phase_comparison.png

In addition to the capabilities demonstrated above, the :class:`~libertem.udf.holography.HoloReconstructUDF`
class can take smoothness of sideband (SB) filter as fraction of the SB size (:code:`sb_smoothness=0.05` is default).
Also, the :code:`precision` argument can be used (:code:`precision=False`) to reduce the calculation precision
to :code:`float32` and :code:`complex64` for the output. Note that depending of NumPy backend, even with reduced
precision the FFT function used in the reconstruction may internally calculate results with double
precision. In this case reducing precision will only affect the size of the output rather than the
speed of processing.
.. `phasecontrast app`:

Phase contrast
==============

Center of mass / first moment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can get a first moment :cite:`ROSE1976251` visualization of your data set by
selecting "Center of mass" in the "Add analysis" dropdown:

..  figure:: ../images/app/com-example.png

Take note of the "Channel" drop down under the right image, where you can select
different visualizations derived from the vector field.

Center of mass is also available in the :class:`~libertem.api.Context` API.
Please see `this example
<https://github.com/LiberTEM/LiberTEM/blob/master/examples/center_of_mass.ipynb>`_
and the reference documentation for
:meth:`~libertem.api.Context.create_com_analysis`.

Ptychographic reconstruction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please see `the ptychography 4.0 project
<https://ptychography-4-0.github.io/ptychography/algorithms.html>`_ for
ptychography algorithms implemented as LiberTEM UDFs.
.. _`phasechange app`:

Phase-change materials
======================

Phase-change materials are some of the most promising materials for data-storage applications. In this application section, a pixelated STEM dataset from a AgInSbTe phase-change material specimen will be processed.
Distinguishing between amorphous and crystalline regions can be hard in microstructure analysis with 4-D-STEM because of low contrast between them during the two most common method of visualization: bright and dark field imaging.
The main difference of pixelated diffraction patterns for each of the beam positions the presence of additional (non zero-order) diffraction peaks in crystalline frames, while the amorphous frames have only the zero-order peak.
Detection of the positions of additional peaks for all frames of crystalline regions will allow gaining information on the crystal structure.

Next, methods to distinguish crystalline and amorphous regions in phase-change materials will be described.

.. note::
    See :ref:`phasechange` for the API reference.

.. _`crystallinity map`:

Crystallinity map
~~~~~~~~~~~~~~~~~
First, two different types of the regions, crystalline (low resistance) and amorphous (high resistance), should be distinguished.
Crystalline regions are characterized by the presence of additional non-zero diffraction peaks on diffraction image, while the diffraction image (the frame) of amorphous regions contains only the zero-order diffraction peak on a diffuse scattering background.
The presence of additional periodic diffraction peaks leads to differences in the Fourier spectrum of the crystalline frames compared to amorphous frames in the intermediate frequency range.
Integration of all frames Fourier spectra over a ring (in intermediate frequency range) and using the result of the integration as intensity for each position for sample visualization (similar to virtual dark field image, but in Fourier space) allows distinguishing amorphous and crystalline regions of the sample.

GUI use:
--------

You can select "FFT Fourier" from the "Add Analysis" menu in the GUI to obtain crystallinity map with high contrast between crystalline and amorphous regions.

..  figure:: ./images/phasechange/mode.PNG 

Use the checkbox in the middle if you would like to enable/disable of masking of zero-order peak out.
In case of enabling zero-order diffraction peak removal: use the controls on the middle to position the disk-shaped selector over the region of average frame you'd like to mask out.
Then, adjust the radii on the left to position the ring-shaped selector over the region of the average of frames spectra you'd like to integrate over and then click "Apply". 

.. figure:: ./images/phasechange/newdatawithmasking.PNG

In the right side, the resulting image will be generated. To check the quality of the detection you can use pick mode which is located under the middle image.
In the case of correct parameters settings, crystalline regions will be brighter than amorphous.

Crystalline frame:

.. figure:: ./images/phasechange/newdatapickcryst.PNG

Amorphous frame:

.. figure:: ./images/phasechange/newdatapickam.PNG

Tips for enabling/disabling of zero-order diffraction peak removal:

In case of the sharp and bright zero-order peak you will need to mask it out, as it described above (spectrum of such a peak will be a 2d sinc function which will spread all over the frames spectra)

.. figure:: ./images/phasechange/newdatawithoutmasking.PNG

So, the crystallinity map will not provide any contrast.

But, in the case of sinc shaped zero-order peak, you do not need to remove it. The influence of it will be presented in low-frequency range and will be masked out during the integration region choosing

.. figure:: ./images/phasechange/interf.PNG

Crystalline frame:

.. figure:: ./images/phasechange/pickcrystalline.PNG

Amorphous frame:

.. figure:: ./images/phasechange/pickam.PNG

.. _clustering:

Clustering
~~~~~~~~~~

.. versionadded:: 0.3.0

To further categorize the crystalline regions according to their lattice
orientation, clustering, based on non-zero diffraction peaks positions, can be
used. The scripting interface allows to cluster membrane, amorphous and
crystalline regions with different lattice orientation in a more efficient way,
taking into account regions of interests. To look at full analysis `follow this
link to a Jupyter notebook <pcmclustering.ipynb>`_.

.. toctree::

   pcmclustering

GUI use:
--------

To make preliminary analysis for parameters choice and a brief look at the
result, you can select "Clustering" from the "Add Analysis" menu in the GUI.

Then you can choose the region in the navigation space in the middle to select
the standard deviation (SD) calculation region (recommendation: try to select as
much of the specimen you can, avoiding zones without usable data, such as
membrane or vacuum). The frame regions which will have higher intensity in the
standard deviation image will be assumed as possible positions of non-zero order
diffraction peaks. Adjust the position of the ring and its radii in the left
side to select the region for detecting peaks and click "Apply". Then, for each
of coordinates of the peaks on SD image, the decision about the presence or
absence of a peak at this position for each frame will be made. Next, as a
result, the feature vector will be generated for each frame to further use in
clustering.

.. figure:: ./images/phasechange/first.PNG

After the first result will be shown, you can readjust the region of peak
detection (in the left side), SD calculation region (in the middle) and
parameters, which are hidden in the parameters section. For radii readjustment,
you can use SD over ROI mode. To see the correct ROI, use the "SD over ROI (rect) mode".
To choose ROI for feature vector calculation and clustering, use the middle section.

.. figure:: ./images/phasechange/second.PNG

You can check the clustering result by using Pick mode.

Example of cluster with one lattice orientation:

.. figure:: ./images/phasechange/1.PNG

Example of cluster with another lattice orientation:

.. figure:: ./images/phasechange/2.PNG

.. _`amorphous app`:

Amorphous materials
===================

.. note::

    See :ref:`fem api` and :ref:`radial fourier api` for API references


LiberTEM offers methods to determine the local order or crystallinity of
amorphous and nanocrystalline materials.

Fluctuation EM
~~~~~~~~~~~~~~

The diffraction pattern of amorphous materials show a characteristic ring of
intensity around the zero-order diffraction peak that is a result of near-range
ordering of the atoms in the material. Any local ordering within the interaction
volume of beam and sample will lead to increased "speckle", i.e. intensity
fluctuations within this ring, since regions with local ordering will diffract
the beam to preferred directions similar to a crystal. For Fluctuation EM
:cite:`Gibson1997`, the standard deviation within this ring is calculated to
obtain a measure of this local ordering.

GUI use
-------

You can select `FEM` from the "Add Analysis" menu in the GUI to add a
Fluctuation EM Analysis.

..  figure:: ./images/fem/fem-select.png

Use the controls on the left to position the ring-shaped selector over the
region you'd like to analyze and then click "Apply".

..  figure:: ./images/fem/fem.png

Pick mode or average over ROI comes in handy to inspect individual frames or
regions of interest:

..  figure:: ./images/fem/fem-pick.png

Sample data: Metallic glass courtesy of Shuai Wei
<shuai.wei@physik.rwth-aachen.de>, RWTH Aachen; Alexander Kuball, Universität
des Saarlandes; Hongchu Du <h.du@fz-juelich.de>, ER-C, Forschungszentrum Jülich.

Scripting interface
-------------------

LiberTEM supports this with the :class:`~libertem.udf.FEM.FEMUDF` and
:meth:`~libertem.udf.FEM.run_fem` method.

.. _`radialfourier app`:

Radial Fourier Series
~~~~~~~~~~~~~~~~~~~~~

.. toctree::

   radialfourier

Fluctuation EM doesn't evaluate the spatial distribution of intensity. It only
works if enough intensity accumulates in each pixel so that local ordering leads
to larger intensity variations than just random noise, i.e. if statistical
variations from shot noise average out and variations introduced by the sample
dominate.

If most of the detector pixels are only hit with one electron or none, the
standard deviation between detector positions in a region of interest is the
same, even if pixels that received electrons are spatially close, while other
regions received no intensity. That means Fluctuation EM doesn't produce
contrast between amorphous and crystalline regions if the detector has a high
resolution, and/or if only low intensity is scattered into the region of
interest.

The Radial Fourier Series Analysis solves this problem by calculating a Fourier
series in a ring-shaped region instead of just evaluating the standard
deviation. The reference :cite:`6980942` describes a previous application of
this method as a descriptor for feature extraction from images. The angle of a
pixel relative to the user-defined center point of the diffraction pattern is
used as a phase angle for the Fourier series.

Since `diffraction patterns usually show characteristic
<http://xrayweb.chem.ou.edu/notes/symmetry.html>`_ `symmetries
<https://en.wikipedia.org/wiki/Friedel%27s_law>`_, the strength of the Fourier
coefficients of orders 2, 4 and 6 highlight regions with crystalline order for
even the lowest intensities. With the `relationship between variance in real
space and power spectral density in frequency space
<https://en.wikipedia.org/wiki/Parseval%27s_theorem>`_, the sum of all
coefficients that are larger than zero is equivalent to the standard deviation,
i.e. Fluctuation EM. Only summing coefficients from lower orders corresponds
Fluctuation EM on a smoothened dataset.

The phase angle of selected coefficients, for example first or second order, can
be used for in-plane orientation mapping similar to :cite:`Panova2019`.

Please note that this method is new and experimental and still needs to be
validated further. If you are interested in helping us with applications and
validations, that would be highly appreciated!

GUI use
-------

You can select "Radial Fourier" from the "Add Analysis" menu in the GUI to add a
Radial Fourier Series Analysis.

..  figure:: ./images/radialfourier/select.png 

Use the controls on the left to position the ring-shaped selector over the
region you'd like to analyze and then click "Apply".

..  figure:: ./images/radialfourier/radialfourier.png

Under the right-hand side image you can select the channel to display. Available
are predominant order, absolute value, phase angle and a cubehelix vector field
visualization of the complex number. The orders larger than zero are all plotted
on the same range and are normalized by the zero-order component for each scan
position to make the components comparable and eliminate the influence of
absolute intensity variations in the visual display.

You can select entries with the arrow keys on the keyboard in case the menu is
outside the browser window. Your help with a more user-friendly GUI
implementation would be highly appreciated!

..  figure:: ./images/radialfourier/radialfourier-channel.png

Pick mode comes in handy to inspect individual frames:

..  figure:: ./images/radialfourier/radialfourier-pick.png

Sample data: Metallic glass courtesy of Shuai Wei
<shuai.wei@physik.rwth-aachen.de>, RWTH Aachen; Alexander Kuball, Universität
des Saarlandes; Hongchu Du <h.du@fz-juelich.de>, ER-C, Forschungszentrum Jülich.

Scripting interface
-------------------

The scripting interface
:meth:`~libertem.api.Context.create_radial_fourier_analysis` and
:class:`~libertem.analysis.radialfourier.RadialFourierAnalysis` allows to
calculate more than one bin at a time and influence the number of orders that
are calculated. It relies on the sparse matrix back-end for :code:`ApplyMasksUDF` and allows
to calculate many orders for many bins at once with acceptable efficiency.

Having a fine-grained analysis with many orders calculated as a function of
radius allows for a number of additional processing and analysis steps. `Follow
this link to a Jupyter notebook. <radialfourier.ipynb>`_

.. rubric:: Reference

See :ref:`fem api` and :ref:`radial fourier api` for API references!
.. _`strain mapping`:

Strain mapping
==============

.. versionchanged:: 0.4.0
    Blobfinder has been spun out as a separate package. Documentation is
    available at https://libertem.github.io/LiberTEM-blobfinder/
Setting up a development environment
====================================

When you want to develop LiberTEM itself, or run the latest git version, the installation works a
bit differently as opposed to installing from PyPI.
As a first step, please follow the instructions for creating a Python virtual environment from
the :ref:`installation` section. Then, instead of installing from PyPI, follow the instructions below.

Prerequisites
~~~~~~~~~~~~~

In addition to a Python installation, there are a few system dependencies needed when contributing:

* `pandoc <https://pandoc.org/installing.html>`_ and `graphviz <https://graphviz.org/download/>`_
  are needed to build and test the documentation.
* Node.js is needed for rebuilding the LiberTEM GUI. On Linux, we recommend
  to `install via package manager
  <https://nodejs.org/en/download/package-manager/>`_, on Windows `the installer
  <https://nodejs.org/en/download/>`_ should be fine. Choose the current LTS
  version.

.. _`installing from a git clone`:

Installing from a git clone
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
    Distinguish between installing a released version and installing the latest
    development version. Both :ref:`installing from PyPI` and :ref:`installing from a git
    clone` use pip, but they do fundamentally different things. :code:`python -m pip
    install libertem` downloads the latest release from PyPI.

    Changing directory to a git clone and running :code:`python -m pip install -e .`
    installs from the local directory in editable mode. "Editable mode" means
    that the source directory is "linked" into the current Python environment
    rather than copied. That means changes in the source directory are
    immediately active in the Python environment.

    Installing from a git clone in editable mode is the correct setup for
    development work and using :ref:`the latest features in the development
    branch <continuous>`. Installing from PyPI is easier and preferred for
    production use and for new users.

If you want to follow the latest development, you should install LiberTEM from
a git clone. As a prerequisite, you need to have git installed on your system. On Linux,
we suggest using the git package that comes with your package manager. On Windows, you can use one
of the many available clients, like  `git for windows <https://gitforwindows.org/>`_, 
`GitHub Desktop <https://desktop.github.com/>`_, `TortoiseGit <https://tortoisegit.org/>`_,
or the git integration of your development environment.

.. code-block:: shell

    $ git clone https://github.com/LiberTEM/LiberTEM

Or, if you wish to contribute to LiberTEM, create a fork of the LiberTEM repository
by following these steps:

#. Log into your `GitHub <https://github.com/>`_ account.

#. Go to the `LiberTEM GitHub <https://github.com/liberteM/LiberTEM/>`_ home page.

#. Click on the *fork* button:

    ..  figure:: ../images/forking_button.png

#. Clone your fork of LiberTEM from GitHub to your computer

.. code-block:: shell

    $ git clone https://github.com/your-user-name/LiberTEM

More information about `forking a repository
<https://docs.github.com/en/get-started/quickstart/fork-a-repo>`_.
For a beginner-friendly introduction to git and GitHub, consider going through
the following resources:

* This `free course <https://www.udacity.com/course/version-control-with-git--ud123>`_
  covers the essentials of using Git.
* Practice `pull request <https://github.com/firstcontributions/first-contributions>`_
  in a safe sandbox environment.
* Sample `workflow <https://docs.astropy.org/en/latest/development/workflow/development_workflow.html>`_
  for contributing code.

Activate the Python environment (conda or virtualenv) and change to the newly
created directory with the clone of the LiberTEM repository. Now you can start
the LiberTEM installation. Please note the dot at the end, which indicates the
current directory!

.. code-block:: shell

    (libertem) $ python -m pip install -e .

This should download the dependencies and install LiberTEM in the environment.
Please continue by reading the :ref:`usage documentation`.

Installing extra dependencies works just like when installing LiberTEM from PyPI:

.. code-block:: shell

    (libertem) $ python -m pip install -e .[torch,hdbscan,cupy]

Updating
~~~~~~~~

If you have installed from a git clone, you can easily update it to the current
status. Open a command line in the base directory of the LiberTEM clone and
update the source code with this command:

.. code-block:: shell

    $ git pull

The installation with :code:`python -m pip install -e` has installed LiberTEM in `"editable"
mode <https://pip.pypa.io/en/stable/cli/pip_install/#editable-installs>`_.
That means the changes pulled from git are active immediately. Only if the
requirements for installed third-party packages have changed, you should re-run
:code:`python -m pip install -e .` in order to install any missing packages.

Setting up tox on Windows
~~~~~~~~~~~~~~~~~~~~~~~~~

We are using tox to run our tests.
On Windows with Anaconda, you have to create named aliases for the Python
interpreter before you can run :literal:`tox` so that tox finds the python
interpreter where it is expected. Assuming that you run LiberTEM with Python
3.6, place the following file as :literal:`python3.6.bat` in your LiberTEM conda
environment base folder, typically
:literal:`%LOCALAPPDATA%\\conda\\conda\\envs\\libertem\\`, where the
:literal:`python.exe` of that environment is located.

.. code-block:: bat

    @echo off
    REM @echo off is vital so that the file doesn't clutter the output
    REM execute python.exe with the same command line
    @python.exe %*

To execute tests with Python 3.7, you create a new environment with Python 3.7:

.. code-block:: shell

    > conda create -n libertem-3.7 python=3.7

Now you can create :literal:`python3.7.bat` in your normal LiberTEM environment
alongside :literal:`python3.6.bat` and make it execute the Python interpreter of
your new libertem-3.7 environment:

.. code-block:: bat

    @echo off
    REM @echo off is vital so that the file doesn't clutter the output
    REM execute python.exe in a different environment
    REM with the same command line
    @%LOCALAPPDATA%\conda\conda\envs\libertem-3.7\python.exe %*

See also:
https://tox.readthedocs.io/en/latest/developers.html#multiple-python-versions-on-windows
.. _`executor api`:

LiberTEM executors
==================

All access to data and processing is done by an executor that implements the
:class:`~libertem.executor.base.JobExecutor` interface to run functions and
tasks. That allows to modify where and how processing is done, including running
on a cluster or in a single thread, without changes in other parts of LiberTEM.
See :ref:`executors` for an overview from a user's perspective.

.. versionadded:: 0.9.0
    The executor API is internal. Since choice and parameters of executors
    are important for integration with Dask and other frameworks,
    they are now documented. Only the names and creation methods for
    executors are reasonably stable. The rest of the API is subject to
    change without notice.

.. _`dask executor`:

Dask.Distributed
................

.. automodule:: libertem.executor.dask
    :members:

.. automodule:: libertem.executor.integration
    :members:

Inline
......

.. automodule:: libertem.executor.inline
    :members:

Concurrent
..........

.. automodule:: libertem.executor.concurrent
    :members:

Dask.delayed
............

.. automodule:: libertem.executor.delayed
    :members:
How does I/O work in LiberTEM?
==============================

Many algorithms benefit from :ref:`tiled` where the same slice of the signal
dimension is processed for several frames in a row. In many cases, algorithms
have specific minimum and maximum sizes in signal dimension, navigation
dimension or total size where they operate efficiently. Smaller sizes might
increase overheads, while larger sizes might reduce cache efficiency.

At the same time, file formats might operate well within specific size and
shape limits. The :ref:`k2is` raw format is a prime example where data is saved
in tiled form and can be processed efficiently in specific tile sizes and
shapes that follow the native layout. Furthermore, some formats require
decoding or corrections by the CPU, such as :ref:`frms6`, where tiles that fit
the L3 cache can speed up subsequent processing steps. Requirements from the
I/O method such as alignment and efficient block sizes are taken into account
as well.

The LiberTEM I/O back-end negotiates a tiling scheme between UDF and dataset
that fulfills requirements from both UDF and dataset side as far as possible.
However, it is not always guaranteed that the supplied data will fall within
the requested limits.

.. versionadded:: 0.6.0
  This guide is written for version 0.6.0

High-level overview
~~~~~~~~~~~~~~~~~~~

- Once per :code:`DataSet`, the :code:`UDFRunner` negotiates a :code:`TilingScheme` using the
  :code:`Negotiator` class
- Data is split into partitions which are read from independently. Usually
  they split the navigation axis.
- The :code:`TilingScheme` is then passed on to :code:`Partition.get_tiles`,
  which then yields :code:`DataTiles` that match the given
  :code:`TilingScheme`. Note that the tiles MUST match the slicing in the signal dimensions, but
  the tile depth may differ from the scheme (for example, if partition depth isn't evenly divisible by
  tile depth, or if there are other technical reasons the :code:`DataSet` can't return the given depth)
- Under the hood, the :code:`Partition`...
   - instantiates an :code:`IOBackend`, which has a reference to a :code:`Decoder`
   - generates read ranges, which are passed on to the :code:`IOBackend`
   - delegates :code:`get_tiles` to the :code:`IOBackend`
- The I/O process can be influenced by passing a subclass
  of :code:`FileSet` to the :code:`Partition` and overriding :code:`FileSet.get_read_ranges`,
  implementing a :code:`Decoder`, or even completely overriding
  the :code:`Partition.get_tiles` functionality.
- There are currently three I/O backends implemented: :code:`MMapBackend`,
  :code:`BufferedBackend`, and :code:`DirectBackend`
  which are useful for different storage media and purposes
- :code:`MMapBackend.get_tiles` has two modes of operation: either it returns a reference to the
  tiles "straight" from the file, without copying or decoding, or it
  uses the read ranges and copies/decodes the tiles in smaller units.
- When reading the tiles "straight", the read ranges are not used, instead
  only the slice information for each tile is used. That also means that this
  mode only works for very simple formats, when reading without a :code:`roi`
  and when not doing any :code:`dtype` conversion or decoding.
- For most formats, a :code:`sync_offset` can be specified, which can be used to
  correct for synchronization errors by inserting blank frames,
  or ignoring one or more frames, at the beginning or at the end of the data set.

Read ranges
-----------

In :code:`FileSet.get_read_ranges`, the reading parameters (:code:`TilingScheme`, :code:`roi` etc.)
are translated into one or more byte ranges (offset, length) for each tile.
You can imagine it as translating pixel/element positions into byte offsets.

Each range corresponds to a read operation on a single file, and with multiple
read ranges per tile, ranges for a single tile can correspond to reads from multiple files.
This is important when reading from a data set with many small files - we can
still generate deep tiles for efficient processing.

There are some built-in common parameters in :code:`FileSet`, like
:code:`frame_header_bytes`, :code:`frame_footer_bytes`, which can be used to easily
implement formats where the reading just needs to skip a few bytes for each
frame header/footer.

If you need more influence over how data is read, you can override
:code:`FileSet.get_read_ranges` and return your own read ranges. You can use
the :code:`make_get_read_ranges` function to re-use much of the tiling logic,
or implement this yourself. Using :code:`make_get_read_ranges` you can either
override just the :code:`px_to_bytes` part, or :code:`read_ranges_tile_block` for whole
tile blocks. This is done by passing in njit-ed functions to :code:`make_get_read_ranges`.
:code:`make_get_read_ranges` should only be called on module-level to enable
caching of the numba compilation.

Read ranges are generated as an array with the following shape::

    :code:`(number_of_tiles, rr_per_tile, rr_num_entries)`

:code:`rr_per_tile` here is the maximum number of read ranges per tile - there
can be tiles that are smaller than this, for example at the end of a partition.
:code:`rr_num_entries` is at least 3 and contains at least the values
:code:`(file_idx, start, stop)`. This means to read :code:`stop - start`
bytes, beginning at offset :code:`start`, from the file :code:`file_idx`
in the corresponding :code:`FileSet`.

Overriding :code:`DataSet`\ s are free to add additional fields to the end, for
example if the decoding functions need additional information.

As an example when you would generate custom read ranges, have a look at the
implementations for MIB, K2IS, and FRMS6 - they may not have a direct 1:1 mapping
to a numpy :code:`dtype`, or the pixels may need to be re-ordered after decoding.

Notes for implementing a :class:`~libertem.io.dataset.base.DataSet`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Read file header(s) in :meth:`~libertem.io.dataset.base.DataSet.initialize` -
  make sure to do the actual I/O in a function dispatched via the
  :code:`JobExecutor` that is passed to :code:`initialize`.
  See also :ref:`os mismatch` regarding platform-dependent code.
- Implement :meth:`~libertem.io.dataset.base.DataSet.check_valid` - this will
  be run on a worker node
- Implement :meth:`~libertem.io.dataset.base.DataSet.get_msg_converter` - the
  :class:`~libertem.web.messages.MessageConverter` class returned is responsible
  for parsing parameters passed to the Web API and converting them to a Python
  representation that can be passed to the
  :class:`~libertem.io.dataset.base.DataSet` constructor.
- Implement :meth:`~libertem.io.dataset.base.DataSet.get_cache_key` - the cache
  key must be different for :code:`DataSet`\ s that return different data.
- Implement :meth:`~libertem.io.dataset.base.DataSet.get_partitions`. You may
  want to use the helper function
  :meth:`~libertem.io.dataset.base.DataSet.get_slices` to generate
  slices for a specified number of partitions.
  :meth:`~libertem.io.dataset.base.DataSet.get_partitions` should yield either
  :class:`~libertem.io.dataset.base.BasePartition` instances or instances of
  your own subclass (see below). The same is true for the
  :class:`~libertem.io.dataset.base.FileSet` that is passed to each
  partition - you possibly have to implement your own subclass.
- Implement :meth:`~libertem.io.dataset.base.DataSet.get_decoder` to return
  an instance of :class:`~libertem.io.dataset.base.Decoder`. Only needed if
  the data is saved in a data type that is not directly understood by numpy
  or numba. See below for details.
- Implement :meth:`~libertem.io.dataset.base.DataSet.get_base_shape`. This
  is only needed if the data format imposes any constraints on how the data can be
  read in an efficient manner, for example if data is saved in blocks. The tileshape
  that is negotiated before reading will be a multiple of the base shape in
  all dimensions.
- Implement :meth:`~libertem.io.dataset.base.DataSet.adjust_tileshape`. This
  is needed if you need to "veto" the generated tileshape, for example if your dataset
  has constraints that can't be expressed by the base shape.

Subclass :class:`~libertem.io.dataset.base.BasePartition`
---------------------------------------------------------

- Override :meth:`~libertem.io.dataset.base.BasePartition.get_tiles` if you need to
  use completely custom I/O logic.

Implementing a :class:`~libertem.io.dataset.base.Decoder`
---------------------------------------------------------

This may be needed if the raw data is not directly supported
by numpy or numba. Mostly your decoder will return a different
:code:`decode` function in :meth:`~libertem.io.dataset.base.Decoder.get_decode`.
You can also return different decode functions, depending on the
concrete data set you are currently reading. For example, this may be needed if there
are different data representations generated by different detector modes.
You can also instruct the :code:`IOBackend` to clear the read
buffer before calling :code:`decode` by returning :code:`True` from 
:meth:`~libertem.io.dataset.base.Decoder.do_clear`. This can be needed
if different read ranges contribute to the same part of the output buffer
and the :code:`decode` function accumulates into the buffer instead of slice-assigning.

The :code:`decode` function will be called for each read range that was
generated by the :code:`get_read_ranges` method described above... _`advanced udf`:

User-defined functions: advanced topics
---------------------------------------

.. testsetup:: *

    import numpy as np
    from libertem import api
    from libertem.executor.inline import InlineJobExecutor
    from libertem.viz.base import Dummy2DPlot

    ctx = api.Context(executor=InlineJobExecutor(), plot_class=Dummy2DPlot)
    data = np.random.random((16, 16, 32, 32)).astype(np.float32)
    dataset = ctx.load("memory", data=data, sig_dims=2)
    roi = np.random.choice([True, False], dataset.shape.nav)

The UDF interface offers a wide range of features to help implement advanced
functionality and to optimize the performance of an UDF. These features are
optional in order to keep UDFs that don't need them simple.

See :ref:`user-defined functions` for an introduction to basic topics.

.. _tiled:

Tiled processing
----------------

Many operations can be significantly optimized by working on stacks of frames.
You can often perform `loop nest optimization
<https://en.wikipedia.org/wiki/Loop_nest_optimization>`_ to improve the
`locality of reference <https://en.wikipedia.org/wiki/Locality_of_reference>`_,
for example using `numba <https://numba.pydata.org/>`_, or using an optimized
NumPy function.

As an example, applying a gain map and subtracting dark frames can be up to an
order of magnitude faster when properly optimized compared to a naive NumPy
implementation. These optimizations are only possible if you have access to data
from more than one frame.

For very large frames, another problem arises: a stack of frames would be too
large to efficiently handle, as it would no longer fit into even the L3 cache,
which is the largest cache in most CPUs. For these cases, we support a tiled
reading and processing strategy. Tiled means we slice the frame into disjoint
rectangular regions. A tile then is the data from a single rectangular region
for multiple frames.

For example, in case of K2IS data, frames have a shape of :code:`(1860, 2048)`.
When reading them with the tiled strategy, a single tile will contain data from
16 subsequent frames, and each rectangle has a shape of :code:`(930, 16)`, which
is the natural block size for K2IS data. That means the tiles will have a shape
of :code:`(16, 930, 16)`, and processing 16 frames from the data set means
reading 256 individual tiles.

Loading a tile of this size as float32 data still fits comfortably into usual L3
CPU caches (~1MB per core), and thus enables efficient processing. As a comparison, a
whole :code:`(1860, 2048)` frame is about 15MB large, and accessing it
repeatedly means having to load data from the slower main memory.

.. note::
    You may have noticed that we talk about block sizes of 1MB as efficient in
    the L3 cache, but many CPUs have larger L3 caches. As the L3 cache is shared
    between cores, and LiberTEM tries to use multiple cores, the effectively
    available L3 cache has to be divided by number of cores.

.. _`slice example`:

Real-world example
~~~~~~~~~~~~~~~~~~

The :class:`libertem_blobfinder.udf.correlation.SparseCorrelationUDF` uses
:meth:`~libertem.udf.base.UDFTileMixin.process_tile` to implement a custom version of
a :class:`~libertem.udf.masks.ApplyMasksUDF` that works on log-scaled data. The
mask stack is stored in a :class:`libertem.common.container.MaskContainer` as part of
the task data. Note how the :code:`self.meta.slice` property of type
:class:`~libertem.common.slice.Slice` is used to extract the region from the mask
stack that matches the tile using the facilities of a
:class:`~libertem.common.container.MaskContainer`. After reshaping, transposing and log
scaling the tile data into the right memory layout, the mask stack is applied to
the data with a dot product. The result is *added* to the buffer in order to
merge it with the results of the other tiles because addition is the correct
merge function for a dot product. Other operations would require a different
merge function here, for example :meth:`numpy.max()` if a global maximum is to
be calculated.

.. testsetup::

    from libertem_blobfinder.base.correlation import log_scale

.. testcode::

    def process_tile(self, tile):
        tile_slice = self.meta.slice
        c = self.task_data['mask_container']
        tile_t = np.zeros(
            (np.prod(tile.shape[1:]), tile.shape[0]),
            dtype=tile.dtype
        )
        log_scale(tile.reshape((tile.shape[0], -1)).T, out=tile_t)

        sl = c.get(key=tile_slice, transpose=False)
        self.results.corr[:] += sl.dot(tile_t).T

.. _`udf post processing`:

Partition processing
--------------------

Some algorithms can benefit from processing entire partitions, for example if
they require several passes over the data. In most cases, :ref:`tiled
processing<tiled>` will be faster because it uses the L3 cache more efficiently.
For that reason, per-partition processing should only be used if there are clear
indications for it. Implementing
:meth:`~libertem.udf.base.UDFPartitionMixin.process_partition` activates
per-partition processing for an UDF.

Precedence
----------

The UDF interface looks for methods in the order
:meth:`~libertem.udf.base.UDFTileMixin.process_tile`,
:meth:`~libertem.udf.base.UDFFrameMixin.process_frame`,
:meth:`~libertem.udf.base.UDFPartitionMixin.process_partition`. For now, the first in
that order is executed. In the future, composition of UDFs may allow to use
different methods depending on the circumstances.
:meth:`~libertem.udf.base.UDFTileMixin.process_tile` is the most general method and
allows by-frame and by-partition processing as well.

Post-processing of partition results
------------------------------------

Post-processing allows to perform additional processing steps once the data of a
partition is completely processed with
:meth:`~libertem.udf.base.UDFFrameMixin.process_frame`,
:meth:`~libertem.udf.base.UDFTileMixin.process_tile` or
:meth:`~libertem.udf.base.UDFPartitionMixin.process_partition`. Post-processing is
particularly relevant for tiled processing since that allows to combine the
performance benefits of tiled processing for a first reduction step with
subsequent steps that require reduced data from complete frames or even a
complete partition.

Real-world example from
:class:`libertem_blobfinder.udf.correlation.SparseCorrelationUDF` which
evaluates the correlation maps that have been generated with the dot product in
the previous processing step and places the results in additional result
buffers:

.. testsetup::

    from libertem_blobfinder.base.correlation import evaluate_correlations

.. testcode::

    def postprocess(self):
        steps = 2 * self.params.steps + 1
        corrmaps = self.results.corr.reshape((
            -1,  # frames
            len(self.params.peaks),  # peaks
            steps,  # Y steps
            steps,  # X steps
        ))
        peaks = self.params.peaks
        (centers, refineds, peak_values, peak_elevations) = self.output_buffers()
        for f in range(corrmaps.shape[0]):
            evaluate_correlations(
                corrs=corrmaps[f], peaks=peaks, crop_size=self.params.steps,
                out_centers=centers[f], out_refineds=refineds[f],
                out_heights=peak_values[f], out_elevations=peak_elevations[f]
            )

The :meth:`libertem.udf.base.UDFPostprocessMixin.postprocess` method is called
for each partition on the worker process, before the results from different
partitions have been merged.

.. _`udf final post processing`:

Post-processing after merging
-----------------------------

If you want to implement a post-processing step that is run on the main node
after merging result buffers, you can override
:meth:`libertem.udf.base.UDF.get_results`:

.. testsetup::

    from libertem.udf import UDF

.. testcode::

    class AverageUDF(UDF):
        """
        Like SumUDF, but also computes the average
        """
        def get_result_buffers(self):
            return {
                'sum': self.buffer(kind='sig', dtype=np.float32),
                'num_frames': self.buffer(kind='single', dtype=np.uint64),
                'average': self.buffer(kind='sig', dtype=np.float32, use='result_only'),
            }

        def process_frame(self, frame):
            self.results.sum[:] += frame
            self.results.num_frames[:] += 1

        def merge(self, dest, src):
            dest.sum[:] += src.sum
            dest.num_frames[:] += src.num_frames

        def get_results(self):
            return {
                # NOTE: 'sum' omitted here, will be returned unchanged
                'average': self.results.sum / self.results.num_frames,
            }

    ctx.run_udf(dataset=dataset, udf=AverageUDF())

Note that :meth:`UDF.get_result_buffers` returns a placeholder entry for the
:code:`average` result using :code:`use='result_only'`, which is then filled in
:code:`get_results`.  We don't need to repeat those buffers that should be
returned unchanged; if you want to omit a buffer from the results completely,
you can declare it as private with :code:`self.buffer(..., use='private')` in
:code:`get_result_buffers`.

:meth:`UDF.get_results` should return the results as a dictionary of numpy
arrays, with the keys matching those returned by
:meth:`UDF.get_result_buffers`.

When returned from :meth:`Context.run_udf`, all results are wrapped into
:code:`BufferWrapper` instances. This is done primarily to get convenient
access to a version of the result that is suitable for visualization, even if
a :code:`roi` was used, but still allow access to the raw result using
:attr:`BufferWrapper.raw_data` attribute.

The detailed rules for buffer declarations, :code:`get_result_buffers` and :code:`get_results` are:

1) All buffers are declared in :code:`get_result_buffers`
2) If a buffer is only computed in :code:`get_results`, it should be marked via
   :code:`use='result_only'` so it isn't allocated on workers
3) If a buffer is only used as intermediary result, it should be marked via :code:`use='private'`
4) Not including a buffer in :code:`get_results` means it will either be passed on
   unchanged, or dropped if :code:`use='private'`
5) It's an error to omit an :code:`use='result_only'` buffer in :code:`get_results`
6) It's an error to include a :code:`use='private'` buffer in :code:`get_results`
7) All results are returned from :code:`Context.run_udf` as :code:`BufferWrapper` instances
8) By default, if :code:`get_results` is not implemented, :code:`use='private'` buffers are dropped,
   and others are passed through unchanged

.. versionadded:: 0.7.0
   :meth:`UDF.get_results` and the :code:`use` argument for :meth:`UDF.buffer` were added.

Pre-processing
---------------

Pre-processing allows to initialize result buffers before processing or merging.
This is particularly useful to set up :code:`dtype=object` buffers, for example
ragged arrays, or to initialize buffers for operations where the neutral element
is not 0. :meth:`libertem.udf.base.UDFPreprocessMixin.preprocess` is executed after
all buffers are allocated, but before the data is processed. On the worker nodes
it is executed with views set for the whole partition masked by the current ROI.
On the central node it is executed with views set for the whole dataset masked
by the ROI. 

.. versionadded:: 0.3.0

.. versionchanged:: 0.5.0
    :meth:`libertem.udf.base.UDFPreprocessMixin.preprocess` is executed on the main
    node, too. Views for aux data are set correctly on the main node. Previously,
    it was only executed on the worker nodes.

AUX data
--------

If a parameter is an instance of :class:`~libertem.common.buffers.BufferWrapper`
that was created using the :meth:`~libertem.udf.base.UDF.aux_data` class method, the
UDF interface will interpret it as auxiliary data. It will set the views for
each tile/frame/partition accordingly so that accessing the parameter returns a
view of the auxiliary data matching the data portion that is currently being
processed. That way, it is possible to pass parameters individually for each
frame or to mask the signal dimension.

Note that the :class:`~libertem.common.buffers.BufferWrapper` instance for AUX
data should always be created using the :meth:`~libertem.udf.base.UDF.aux_data` class
method and not directly by instantiating a
:class:`~libertem.common.buffers.BufferWrapper` since
:meth:`~libertem.udf.base.UDF.aux_data` ensures that it is set up correctly.

For masks in the signal dimension that are used for dot products in combination
with per-tile processing, a :class:`~libertem.common.container.MaskContainer` allows
to use more advanced slicing and transformation methods targeted at preparing
mask stacks for optimal dot product performance.

Task data
---------

A UDF can generate task-specific intermediate data on the worker nodes by
defining a :meth:`~libertem.udf.base.UDF.get_task_data` method. The result is
available as an instance of :class:`~libertem.udf.base.UDFData` in
:code:`self.task_data`. Depending on the circumstances, this can be more
efficient than making the data available as a parameter since it avoids
pickling, network transport and unpickling.

This non-trivial example from
:class:`libertem_blobfinder.udf.correlation.SparseCorrelationUDF` creates
a :class:`~libertem.common.container.MaskContainer` based on the parameters in
:code:`self.params`. This :class:`~libertem.common.container.MaskContainer` is then
available as :code:`self.task_data['mask_container']` within the processing
functions.

.. testsetup::

    from libertem.common.container import MaskContainer
    import libertem.masks as masks

.. testcode::

    def get_task_data(self):
        match_pattern = self.params.match_pattern
        crop_size = match_pattern.get_crop_size()
        size = (2 * crop_size + 1, 2 * crop_size + 1)
        template = match_pattern.get_mask(sig_shape=size)
        steps = self.params.steps
        peak_offsetY, peak_offsetX = np.mgrid[-steps:steps + 1, -steps:steps + 1]

        offsetY = self.params.peaks[:, 0, np.newaxis, np.newaxis] + peak_offsetY - crop_size
        offsetX = self.params.peaks[:, 1, np.newaxis, np.newaxis] + peak_offsetX - crop_size

        offsetY = offsetY.flatten()
        offsetX = offsetX.flatten()

        stack = functools.partial(
            masks.sparse_template_multi_stack,
            mask_index=range(len(offsetY)),
            offsetX=offsetX,
            offsetY=offsetY,
            template=template,
            imageSizeX=self.meta.dataset_shape.sig[1],
            imageSizeY=self.meta.dataset_shape.sig[0]
        )
        # CSC matrices in combination with transposed data are fastest
        container = MaskContainer(mask_factories=stack, dtype=np.float32,
            use_sparse='scipy.sparse.csc')

        kwargs = {
            'mask_container': container,
            'crop_size': crop_size,
        }
        return kwargs

.. testcleanup::

    from libertem_blobfinder.udf.correlation import SparseCorrelationUDF
    from libertem_blobfinder.common.patterns import RadialGradient

    class TestUDF(SparseCorrelationUDF):
        pass

    # Override methods with functions that are defined above

    TestUDF.process_tile = process_tile
    TestUDF.postprocess = postprocess
    TestUDF.get_task_data = get_task_data

    u = TestUDF(
        peaks=np.array([(8, 8)]),
        match_pattern=RadialGradient(2),
        steps=3
    )
    ctx.run_udf(dataset=dataset, udf=u)

Meta information
----------------

Advanced processing routines may require context information about the processed
data set, ROI and current data portion being processed. This information is
available as properties of the :attr:`libertem.udf.base.UDF.meta` attribute of type
:class:`~libertem.udf.base.UDFMeta`.

Input data shapes and types
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Common applications include allocating buffers with a :code:`dtype` or shape
that matches the dataset or partition via
:attr:`~libertem.udf.base.UDFMeta.dataset_dtype`,
:attr:`~libertem.udf.base.UDFMeta.input_dtype`,
:attr:`~libertem.udf.base.UDFMeta.dataset_shape` and
:attr:`~libertem.udf.base.UDFMeta.partition_shape`.

Device class
~~~~~~~~~~~~

.. versionadded:: 0.6.0

The currently used compute device class can be accessed through
:attr:`libertem.udf.base.UDFMeta.device_class`. It defaults to 'cpu' and can be 'cuda'
for UDFs that make use of :ref:`udf cuda` support.

ROI and current slice
~~~~~~~~~~~~~~~~~~~~~

For more advanced applications, the ROI and currently processed data portion are
available as :attr:`libertem.udf.base.UDFMeta.roi` and
:attr:`libertem.udf.base.UDFMeta.slice`. This allows to replace the built-in masking
behavior of :class:`~libertem.common.buffers.BufferWrapper` for result buffers
and aux data with a custom implementation. The :ref:`mask container for tiled
processing example<slice example>` makes use of these attributes to employ a
:class:`libertem.common.container.MaskContainer` instead of a :code:`shape="sig"`
buffer in order to optimize dot product performance and support sparse masks.

The slice is in the reference frame of the dataset, masked by the current ROI,
with flattened navigation dimension. This example illustrates the behavior by
implementing a custom version of the :ref:`simple "sum over sig" example
<sumsig>`. It allocates a custom result buffer that matches the navigation
dimension as it appears in processing:

.. testcode::

    import numpy as np

    from libertem.udf import UDF

    class PixelsumUDF(UDF):
        def get_result_buffers(self):
            if self.meta.roi is not None:
                navsize = np.count_nonzero(self.meta.roi)
            else:
                navsize = np.prod(self.meta.dataset_shape.nav)
            return {
                'pixelsum_nav_raw': self.buffer(
                    kind="single",
                    dtype=self.meta.dataset_dtype,
                    extra_shape=(navsize, ),
                )
            }

        def merge(self, dest, src):
            dest.pixelsum_nav_raw[:] += src.pixelsum_nav_raw

        def process_frame(self, frame):
            np_slice = self.meta.slice.get(nav_only=True)
            self.results.pixelsum_nav_raw[np_slice] = np.sum(frame)

.. testcleanup::

    pixelsum = PixelsumUDF()
    res = ctx.run_udf(dataset=dataset, udf=pixelsum, roi=roi)

    assert np.allclose(res['pixelsum_nav_raw'].data, dataset.data[roi].sum(axis=(1, 2)))

Coordinates
~~~~~~~~~~~

.. versionadded:: 0.6.0

The coordinates of the current frame, tile or partition within the true dataset
navigation dimension, as opposed to the current slice that is given in flattened
nav dimensions with applied ROI, is available through
:attr:`~libertem.udf.base.UDFMeta.coordinates`. The following UDF simply
collects the coordinate info for demonstration purposes. A real-world example
that uses the coordinates is `the UDF implementation of single side band
ptychography
<https://github.com/Ptychography-4-0/ptychography/blob/master/src/ptychography40/reconstruction/ssb/udf.py>`_.

.. testcode::

    import numpy as np

    from libertem.udf import UDF

    class CoordUDF(UDF):
        def get_result_buffers(self):
            # Declare a buffer that fits the coordinates,
            # i.e. one int per nav axis for each nav position
            nav_dims = len(self.meta.dataset_shape.nav)
            return {
                'coords': self.buffer(
                    kind="nav",
                    dtype=int,
                    extra_shape=(nav_dims, ),
                )
            }

        def process_tile(self, tile):
            # Simply copy the coordinates into
            # the result buffer
            self.results.coords[:] = self.meta.coordinates

    my_roi = np.zeros(dataset.shape.nav, dtype=bool)
    my_roi[7, 13] = True
    my_roi[11, 3] = True

    res = ctx.run_udf(
        dataset=dataset,
        udf=CoordUDF(),
        roi=my_roi
    )

    assert np.all(
        res['coords'].raw_data == np.array([(7, 13), (11, 3)])
    )

.. _`udf dtype`:

Preferred input dtype
---------------------

.. versionadded:: 0.4.0

UDFs can override :meth:`~libertem.udf.base.UDF.get_preferred_input_dtype` to
indicate a "lowest common denominator" compatible dtype. The actual input dtype
is determined by combining the indicated preferred dtype with the input
dataset's native dtype using :func:`numpy.result_type`. The default preferred
dtype is :attr:`numpy.float32`. Returning :attr:`UDF.USE_NATIVE_DTYPE`, which is
currently identical to :code:`numpy.bool`, will switch to the dataset's native
dtype since :code:`numpy.bool` behaves as a neutral element in
:func:`numpy.result_type`.

If an UDF requires a specific dtype rather than only preferring it, it should
override this method and additionally check the actual input type, throw an
error when used incorrectly and/or implement a meaningful conversion in its
processing routine since indicating a preferred dtype doesn't enforce it. That
way, unsafe conversions are performed explicitly in the UDF rather than
indirectly in the back-end.

.. _`udf cuda`:

CuPy support
------------

.. versionadded:: 0.6.0

LiberTEM can use CUDA devices through `CuPy <https://cupy.dev/>`_. Since
CuPy largely replicates the NumPy array interface, any UDF that uses NumPy for
its main processing can likely be ported to use both CPUs and CUDA devices in
parallel. Some adjustments are often necessary to account for minor differences
between NumPy and CuPy. CuPy is most beneficial for compute-heavy tasks with
good CUDA math library support such as large Fourier transforms or matrix
products.

In order to activate CuPy processing, a UDF can overwrite the
:meth:`~libertem.udf.base.UDF.get_backends` method. By default this returns
:code:`('numpy',)`, indicating only NumPy support. By returning :code:`('numpy',
'cupy')` or :code:`('cupy',)`, a UDF activates being run on both CUDA and CPU
workers, or exclusively on CUDA workers. Using :code:`cuda` instead of
:code:`cupy` schedules on CUDA workers, but without using the CuPy library. This
is useful for running code that uses CUDA in a different way, for example
integration of C++ CUDA code, and allows to skip installation of CuPy in this
situation.

The :attr:`libertem.udf.base.UDF.xp` property points to the :code:`numpy` or
:code:`cupy` module, depending which back-end is currently used. By using
:code:`self.xp` instead of the usual :code:`np` for NumPy, one can write UDFs
that use the same code for CUDA and CPU processing.

Result buffers can be declared as device arrays by setting
:code:`self.buffer(..., where='device')` in
:meth:`~libertem.udf.base.UDF.get_result_buffers`. That allows to keep data in
the device until a partition is completely processed and the result is exported
to the leader node.

The input argument for :code:`process_*()` functions is already provided as a
CuPy array instead of NumPy array if CuPy is used.

A UDF should only use one GPU at a time. If :code:`cupy` is used, the correct
device to use is set within CuPy in the back-end and should not be modified in
the UDF itself. If :code:`cuda` is used, it is the responsibility of the user to
set the device ID to the value returned by
:meth:`libertem.common.backend.get_use_cuda`. The environment variable
:code:`CUDA_VISIBLE_DEVICES` can be set `before` any CUDA library is loaded to
control which devices are visible.

The :meth:`~libertem.api.Context.run_udf` method allows setting the
:code:`backends` attribute to :code:`('numpy',)` :code:`('cupy',)` or :code:`('cuda',)` to
restrict execution to CPU-only or CUDA-only on a hybrid cluster. This is mostly
useful for testing.

.. _`threading`:

Threading
---------

By default, LiberTEM uses multiprocessing with one process per CPU core for
offline processing. In that scenario, UDFs should only use a single thread to
avoid oversubscription.

However, when running with a single-process single-thread executor like
:class:`~libertem.executor.inline.InlineJobExecutor`, multiple threads can be
used. This is particularly relevant to optimize performance for live processing
in `LiberTEM-live <https://libertem.github.io/LiberTEM-live/>`_ since that may
only be using a single process. The thread count for many common numerics
libraries is set automatically, see
:attr:`~libertem.udf.base.UDFMeta.threads_per_worker`. For other cases the
thread count on a worker should be set by the user according to
:code:`self.meta.threads_per_worker`.

Multithreaded executors are introduced with release 0.9, see :ref:`executors`.
They run UDF functions in parallel threads within the same process. This can
introduce issues with thread safety, for example shared objects being changed
concurrently by multiple threads. The LiberTEM internals and core UDFs are
tested to work with these executors, but user code may break unexpectedly.
`PyFFTW interface caching is a known issue of this category
<https://github.com/LiberTEM/LiberTEM-blobfinder/issues/35>`_. For that reason,
the threaded executors should be considererd experimental for the time being.
Furthermore, setting and re-setting any global variable, for example the thread
count of an external library, should be protected with a reentrant locking
mechanism.

The pyFFTW cache is disabled with threaded executors because of this known bug.
That can have a negative impact on performance. For performance optimization
with pyFFTW, users could use the `builder interface of PyFFTW
<https://pyfftw.readthedocs.io/en/latest/source/pyfftw/builders/builders.html>`_
or `use the native FFTW object interface
<https://pyfftw.readthedocs.io/en/latest/source/pyfftw/pyfftw.html#pyfftw.FFTW>`_.

Running multiple LiberTEM :class:`~libertem.api.Context` objects resp. executors
in parallel threads is not tested and can lead to unexpected interactions.

.. _auto UDF:

Auto UDF
--------

The :class:`~libertem.udf.AutoUDF` class and :meth:`~libertem.api.Context.map`
method allow to run simple functions that accept a frame as the only parameter
with an auto-generated :code:`kind="nav"` result buffer over a dataset ad-hoc
without defining an UDF class. For more advanced processing, such as custom
merge functions, post-processing or performance optimization through tiled
processing, defining an UDF class is required.

As an alternative to Auto UDF, you can use the
:meth:`~libertem.contrib.daskadapter.make_dask_array` method to create
a `dask.array <https://docs.dask.org/en/latest/array.html>`_ from
a :class:`~libertem.io.dataset.base.DataSet` to perform calculations. See
:ref:`Integration with Dask arrays<daskarray>` for more details.

The :class:`~libertem.udf.AutoUDF` class determines the output shape and type
by calling the function with a mock-up frame of the same type and shape as
a real detector frame and converting the return value to a NumPy array. The
:code:`extra_shape` and :code:`dtype` parameters for the result buffer are
derived automatically from this NumPy array.

Additional constant parameters can be passed to the function via
:meth:`functools.partial`, for example. The return value should be much smaller
than the input size for this to work efficiently.

Example: Calculate sum over the last signal axis.

.. testcode::

    import functools

    result = ctx.map(
        dataset=dataset,
        f=functools.partial(np.sum, axis=-1)
    )

    # or alternatively:
    from libertem.udf import AutoUDF

    udf = AutoUDF(f=functools.partial(np.sum, axis=-1))
    result = ctx.run_udf(dataset=dataset, udf=udf)
.. _`udf profiling`:

.. testsetup:: *

    from libertem.api import Context


Profiling UDFs
==============

To bring your UDF to production, it may be necessary to look into the performance and
efficiency of your function. There may be some low-hanging fruit to optimize away, but
they only really become visible once you profile your code. Just looking at code and guessing
which parts are expensive can be misleading and cause you to spend time optimizing the wrong part
of your code. Always measure!

Profiling means instrumenting the execution of a program and finding out which parts
of your code use how much of a given resource (for example, CPU time, or memory). 

Prerequisite: :code:`InlineJobExecutor`
---------------------------------------

By default, your UDF will be run in a multi-process or multi-threaded
environment. This makes profiling a challenge, since the usual profiling tools
will only capture information on the main process that performs the final
reduction and not the processes that do most of the work.

So, in order to use the default Python profiling mechanisms, all parts of the UDF
should be executed in the same single thread. LiberTEM provides what it calls
:class:`~libertem.executor.inline.InlineJobExecutor` for this purpose. Improving
the execution time in a single-threaded executor will almost always improve the
multi-process execution time. Nevertheless, you should confirm the impact of
changes in your production environment.

To use the :class:`~libertem.executor.inline.InlineJobExecutor`, pass it to the
:class:`~libertem.api.Context` on construction:

.. testcode::
   
   from libertem.executor.inline import InlineJobExecutor
   ctx = Context(executor=InlineJobExecutor())

Then, you can continue as usual, loading data, executing your UDF, etc.

Using progress bar
---------------------

A progress bar is a graphical tool that shows a far a process has progressed.
It can be used to assess the partitioning of data, since the progress is
indicated on a per partition basis in the following way:

.. code-block:: text
     
    res = ctx.run_udf(udf=fit_udf, dataset=ds, roi=roi, progress=True)
    > 100%|██████████| 18/18 [01:17<00:00,  4.29s/it]

This run is using 18 partitions, and took 4.29s per partition. 

Line profiling using `line_profiler`
------------------------------------------

Profiling comes in different forms, some will show you how much time each function
of a whole program took, others can give you a line-by-line overview for functions
you are interested in, like `line_profiler`. For most UDFs this will be sufficient
for performance analysis.

First, you need to install the `line_profiler` package via pip, having your
conda environment or virtualenv activated:

.. code-block:: shell

   (libertem) $ python -m pip install line_profiler

Then the easiest way to get started is to use the IPython/Jupyter integration of
`line_profiler`. Put :code:`%load_ext line_profiler` somewhere in your notebook,
then you can use the :code:`%lprun` magic command to run your UDF while profiling:

.. code-block:: text

   %lprun -f YourUDF.get_task_data -f YourUDF.process_frame ctx.run_udf(dataset=dataset, YourUDF())

Note the repeated :code:`-f` arguments - you can list any methods of your UDF you are
interested in, or even other Python functions you are directly or indirectly using. If you
had some complicated code in :code:`YourUDF.merge`, you would include :code:`-f YourUDF.merge`
in the :code:`%lprun` call.

This is how the output could look like:

.. code-block:: text

   Timer unit: 1e-06 s

   Total time: 0.002283 s
   File: /home/clausen/source/libertem/src/libertem/udf/holography.py
   Function: get_task_data at line 145

   Line #      Hits         Time  Per Hit   % Time  Line Contents
   ==============================================================
      145                                               def get_task_data(self):
      146                                                   """
      147                                                   Updates `task_data`
      148                                           
      149                                                   Returns
      150                                                   -------
      151                                                   kwargs : dict
      152                                                   A dictionary with the following keys:
      153                                                       kwargs['aperture'] : array-like
      154                                                       Side band filter aperture (mask)
      155                                                       kwargs['slice'] : slice
      156                                                       Slice for slicing FFT of the hologram
      157                                                   """
      158                                           
      159         2         48.0     24.0      2.1          out_shape = self.params.out_shape
      160         2         51.0     25.5      2.2          sy, sx = self.meta.partition_shape.sig
      161         2          5.0      2.5      0.2          oy, ox = out_shape
      162         2          7.0      3.5      0.3          f_sampling = (1. / oy, 1. / ox)
      163         2        292.0    146.0     12.8          sb_size = self.params.sb_size * np.mean(f_sampling)
      164         2        261.0    130.5     11.4          sb_smoothness = sb_size * self.params.sb_smoothness * np.mean(f_sampling)
      165                                           
      166         2       1172.0    586.0     51.3          f_freq = freq_array(out_shape)
      167         2        263.0    131.5     11.5          aperture = aperture_function(f_freq, sb_size, sb_smoothness)
      168                                           
      169         2         64.0     32.0      2.8          y_min = int(sy / 2 - oy / 2)
      170         2         37.0     18.5      1.6          y_max = int(sy / 2 + oy / 2)
      171         2         32.0     16.0      1.4          x_min = int(sx / 2 - ox / 2)
      172         2         30.0     15.0      1.3          x_max = int(sx / 2 + oy / 2)
      173         2          8.0      4.0      0.4          slice_fft = (slice(y_min, y_max), slice(x_min, x_max))
      174                                           
      175                                                   kwargs = {
      176         2          4.0      2.0      0.2              'aperture': aperture,
      177         2          6.0      3.0      0.3              'slice': slice_fft
      178                                                   }
      179         2          3.0      1.5      0.1          return kwargs

   Total time: 63.748 s
   File: /home/clausen/source/libertem/src/libertem/udf/holography.py
   Function: process_frame at line 181

   Line #      Hits         Time  Per Hit   % Time  Line Contents
   ==============================================================
      181                                               def process_frame(self, frame):
      182                                                   """
      183                                                   Reconstructs holograms outputting results into 'wave'
      184                                           
      185                                                   Parameters
      186                                                   ----------
      187                                                   frame
      188                                                      single frame (hologram) of the data
      189                                                   """
      190        16        154.0      9.6      0.0          if not self.params.precision:
      191                                                       frame = frame.astype(np.float32)
      192                                                   # size_x, size_y = self.params.out_shape
      193        16         81.0      5.1      0.0          frame_size = self.meta.partition_shape.sig
      194        16         58.0      3.6      0.0          sb_pos = self.params.sb_position
      195        16         66.0      4.1      0.0          aperture = self.task_data.aperture
      196        16         52.0      3.2      0.0          slice_fft = self.task_data.slice
      197                                           
      198        16   59291808.0 3705738.0     93.0          fft_frame = fft2(frame) / np.prod(frame_size)
      199        16    2189960.0 136872.5      3.4          fft_frame = np.roll(fft_frame, sb_pos, axis=(0, 1))
      200                                           
      201        16    2258700.0 141168.8      3.5          fft_frame = fftshift(fftshift(fft_frame)[slice_fft])
      202                                           
      203        16        816.0     51.0      0.0          fft_frame = fft_frame * aperture
      204                                           
      205        16       5957.0    372.3      0.0          wav = ifft2(fft_frame) * np.prod(frame_size)
      206        16        364.0     22.8      0.0          self.results.wave[:] = wav

Things to note:

 * :code:`get_task_data` takes a very small amount of time, compared to :code:`process_frame`. It does
   not make sense to concentrate on optimizing :code:`get_task_data` at all, in this case!
 * In :code:`process_frame`, the :code:`fft2` call takes up most time, so that is where
   we should direct our efforts. Improving, for example, the calls to :code:`fftshift` would give us
   a max speed-up of a few percent - and only, if we manage to dramatically improve their execution time!
 * `line_profiler` doesn't give information about individual expressions - sometimes you have to
   put expressions on their own line to see their individual contributions to the execution time. See
   the :code:`fft2` and :code:`np.prod` calls on the hottest line in the profile!
 * After successfully improving on the profiled times, always re-run with profiling disabled and without
   :class:`~libertem.executor.inline.InlineJobExecutor` and measure the total time, for example using
   :code:`%%time`. This makes sure that your optimizations actually work in a production environment!
 * The usual benchmarking rules apply - for example, try to run the profiling on an otherwise idle system,
   otherwise you can get noisy results.
 * Single-threaded execution can be quite slow compared to using LiberTEM in production - if it is too slow
   for your taste, you can run your UDF on a subset of your data using a :ref:`region of interest <udf roi>`.

.. seealso::

   `Python Data Science Handbook <https://jakevdp.github.io/PythonDataScienceHandbook/01.07-timing-and-profiling.html#Line-By-Line-Profiling-with-%lprun>`_
      The Python Data Science Handbook has a section on profiling and timing, including `line_profiler`.

   `Official documentation for line_profiler <https://github.com/rkern/line_profiler>`_
      All information on how to use `line_profiler`, including using it from different contexts.

   :ref:`Profiling long-running tests <profiling tests>`
      Information on how to profile the execution time of test cases.

   :ref:`Debugging UDFs`
      Using the :code:`InlineJobExecutor` to debug problems in your UDF.
.. testsetup:: *

   import numpy as np
   from libertem import api
   from libertem.executor.inline import InlineJobExecutor

   from libertem.udf import UDF as RealUDF

   # We override UDF in such a way that it can be used
   # without implementing all methods
   class UDF(RealUDF):
      def get_result_buffers(self):
         return {}

      def process_frame(self, frame):
         pass

   class YourUDF(RealUDF):
      def get_result_buffers(self):
         return {'buf1': self.buffer(kind="nav")}

      def process_frame(self, frame):
         self.results.buf1[:] = 42

   ctx = api.Context(executor=InlineJobExecutor())
   data = np.random.random((16, 16, 32, 32)).astype(np.float32)
   dataset = ctx.load("memory", data=data, sig_dims=2)
   roi = np.random.choice([True, False], dataset.shape.nav)
   udf = YourUDF()


.. _`implement udf`:

Implementing a UDF
------------------

The workflow for implementing a UDF starts with subclassing
:class:`~libertem.udf.base.UDF`. In the simplest case, you need to implement the
:meth:`~libertem.udf.base.UDF.get_result_buffers` method and 
:meth:`~libertem.udf.base.UDFFrameMixin.process_frame`.

There are two very common patterns for reductions, reducing over the navigation axes
into a common accumulator for all frames, keeping the shape of a single frame,
or reducing over the signal axes and keeping the navigation axes.

Declaring buffers
~~~~~~~~~~~~~~~~~

A UDF can implement one of these reductions or combinations. To handle indexing
for you, LiberTEM needs to know about the structure of your reduction. You can
build this structure in the :meth:`~libertem.udf.base.UDF.get_result_buffers`
method, by declaring one or more buffers.

These buffers can have a :code:`kind` declared, which corresponds to the two
reduction patterns above: :code:`kind="sig"` for reducing over the navigation
axes (and keeping the signal axes), and :code:`kind="nav"` for reducing over the
signal axes and keeping the navigation axes. There is a third,
:code:`kind="single"`, which allows to declare buffers with custom shapes that
don't correspond directly to the data set's shape.

It is also possible to append additional axes to the buffer's shape using the
:code:`extra_shape` parameter.

:meth:`~libertem.udf.base.UDF.get_result_buffers` should return a :code:`dict` which maps
buffer names to buffer declarations. You can create a buffer declaration by calling
the :meth:`~libertem.udf.base.UDF.buffer` method.

The buffer name is later used to access the buffer via :code:`self.results.<buffername>`,
which returns a view into a NumPy array. For this to work, the name has to be a valid Python
identifier.

Examples of buffer declarations:

.. testcode:: getresultbuffers

   def get_result_buffers(self):
      # Suppose our dataset has the shape (14, 14, 32, 32),
      # where the first two dimensions represent the navigation
      # dimension and the last two dimensions represent the signal
      # dimension.

      buffers = {
         # same shape as navigation dimensions of dataset, plus two
         # extra dimensions of shape (3, 2). The full shape is
         # (14, 14, 3, 2) in this example. This means this buffer can
         # store an array of shape (3, 2) for each frame in the dataset.
         "nav_buffer": self.buffer(
               kind="nav",
               extra_shape=(3, 2),
               dtype="float32",
         ),

         # same shape as signal dimensions of dataset, plus an extra
         # dimension of shape (2,). Consequently, the full shape is
         # (32, 32, 2) in this example. That means we can store two
         # float32 values for each pixel of the signal dimensions.
         "sig_buffer": self.buffer(
               kind="sig",
               extra_shape=(2,),
               dtype="float32",
         ),

         # buffer of shape (16, 16); shape is unrelated to dataset shape
         "single_buffer": self.buffer(
               kind="single",
               extra_shape=(16, 16),
               dtype="float32",
         ),

      }

      return buffers

.. testcleanup:: getresultbuffers

   class TestUDF(UDF):
      def merge(self, dest, src):
         pass

   TestUDF.get_result_buffers = get_result_buffers

   u = TestUDF()
   ctx.run_udf(dataset=dataset, udf=u)

See below for some more real-world examples.

All NumPy dtypes are supported for buffers. That includes the :code:`object`
dtype for arbitrary Python variables. The item just has to be pickleable with
:code:`cloudpickle`.

Note that buffers are only designed to pass lightweight intermediate results
and thus, it is important that the size of the buffer remains small. Having too
large buffers can lead to significant decline in performance.

Implementing the processing function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now to the actual core of the processing: implementing
:meth:`~libertem.udf.UDFFrameMixin.process_frame`. The method signature looks like this:

.. testcode:: processing

   class ExampleUDF(UDF):
      def process_frame(self, frame):
         ...

.. testcleanup:: processing

   u = ExampleUDF()
   ctx.run_udf(dataset=dataset, udf=u)

The general idea is that you get a single frame from the data set, do your processing,
and write the results to one of the previously declared buffers, via :code:`self.results.<buffername>`.
When accessing a :code:`kind="nav"` buffer this way, you automatically get a view into the buffer
that corresponds to the current frame that is being processed. In case of :code:`kind="sig"`
or :code:`kind="single"`, you get the whole buffer.

Intuitively, with :code:`kind="sig"` (and :code:`kind="single"`), you are most
likely implementing an operation like :code:`buf = f(buf, frame)`. That is, you
are computing a new result based on a single (new) frame and the results from all
previous frames, and overwrite the results with the new value(s).

With :code:`kind="nav"`, you compute independent results for each frame,
which are written to different positions in the result buffer. Because of the independence
between frames, you don't need to merge with a previous value; the result is simply written
to the correct index in the result buffer (via the aforementioned view).

As an easy example, let's have a look at a function that simply sums up each frame
to a single value. This is a :code:`kind="nav"` reduction, as we sum over all values
in the signal dimensions:

.. _`sumsig`:
   
.. testcode:: sumsig

   import numpy as np   
   from libertem.udf import UDF


   class SumOverSig(UDF):
      def get_result_buffers(self):
         """
         Describe the buffers we need to store our results:
         kind="nav" means we want to have a value for each coordinate
         in the navigation dimensions. We name our buffer 'pixelsum'.
         """
         return {
            'pixelsum': self.buffer(
               kind="nav", dtype="float32"
            )
         }

      def process_frame(self, frame):
         """
         Sum up all pixels in this frame and store the result in the
         `pixelsum` buffer. `self.results.pixelsum` is a view into the
         result buffer we defined above, and corresponds to the entry
         for the current frame we work on. We don't have to take care
         of finding the correct index for the frame we are processing
         ourselves.
         """
         self.results.pixelsum[:] = np.sum(frame)

   res = ctx.run_udf(
      udf=SumOverSig(),
      dataset=dataset,
   )

   # to access the named buffer as a NumPy array:
   res['pixelsum'].data

On a 4D data set, this operation is roughly equivalent to :code:`np.sum(arr, axis=(2, 3))`.

Merging partial results
~~~~~~~~~~~~~~~~~~~~~~~

As :ref:`described above <how UDFs work>`, data from multiple partitions is
processed in parallel. That also means that we need a way of merging partial
results into the final result. In the example above, we didn't need to do anything:
we only have a :code:`kind="nav"` buffer, where merging just means assigning the
result of one partition to the right slice in the final result. This is done by
the default implementation of :meth:`~libertem.udf.base.UDF.merge`. 

In case of :code:`kind="sig"` buffers and the corresponding reduction, assignment would
just overwrite the result from the previous partition with the one from the current partition,
and is not the correct operation. So let's have a look at the merge method:

.. testcode:: merge

   class ExampleUDF(UDF):
      def merge(self, dest, src):
         pass

.. testcleanup:: merge

   u = ExampleUDF()
   ctx.run_udf(dataset=dataset, udf=u)

:code:`dest` is the result of all previous merge calls, and :code:`src` is the
result from a single new partition. Your :code:`merge` implementation should read from both
:code:`dest` and :code:`src` and write the result back to :code:`dest`.

Here is an example demonstrating :code:`kind="sig"` buffers and the :code:`merge` function:
   
.. testcode:: realmerge

   import numpy as np
   from libertem.udf import UDF


   class MaxUDF(UDF):
      def get_result_buffers(self):
         """
         Describe the buffers we need to store our results:
         kind="sig" means we want to have a value for each coordinate
         in the signal dimensions (i.e. a value for each pixel of the
         diffraction patterns). We name our buffer 'maxbuf'.
         """
         return {
            'maxbuf': self.buffer(
               kind="sig", dtype=self.meta.input_dtype
            )
         }

      def preprocess(self):
         """
         Initialize buffer with neutral element for maximum.
         """
         self.results.maxbuf[:] = -np.inf

      def process_frame(self, frame):
         """
         In this function, we have a frame and the buffer `maxbuf`
         available, which we declared above. This function is called
         for all frames / diffraction patterns in the data set. The
         maxbuf is a partial result, and all partial results will
         later be merged (see below).

         In this case, we determine the maximum from the current
         maximum and the current frame, for each pixel in the
         diffraction pattern.

         Notes:

         - You cannot rely on any particular order of frames this function
           is called in.
         - Your function should be pure, that is, it should not have side
           effects beyond modifying the content of result buffers or task data,
           and should only depend on it's input parameters, including
           the UDF object :code:`self`.
         """
         self.results.maxbuf[:] = np.maximum(frame, self.results.maxbuf)

      def merge(self, dest, src):
         """
         merge two partial results, from src into dest
         """
         dest.maxbuf[:] = np.maximum(dest.maxbuf, src.maxbuf)

   res = ctx.run_udf(
      udf=MaxUDF(),
      dataset=dataset,
   )

   # to access the named buffer as a NumPy array:
   res['maxbuf'].data


For more complete examples, you can also have a look at the functions
implemented in the sub-modules of :code:`libertem.udf` and at our :ref:`packages`.

Passing parameters
~~~~~~~~~~~~~~~~~~

By default, keyword arguments that are passed to the constructor of a UDF are
available as properties of :code:`self.params`. This allows clean handling and passing
of parameters in distributed execution scenarios, see below.

.. testsetup:: params

   def correlate_peaks(frame, peaks):
      pass

   peaks = None

.. testcode:: params

    class MyUDF(UDF):

        def process_frame(self, frame):
            result = correlate_peaks(frame, self.params.peaks)
            ...

    udf = MyUDF(peaks=peaks, other=...)

.. testcleanup:: params

   def get_result_buffers():
      return {}

   udf.get_result_buffers = get_result_buffers

   ctx.run_udf(dataset=dataset, udf=udf)

Declaring a constructor
~~~~~~~~~~~~~~~~~~~~~~~

In order to document parameters of an UDF and avoid typos with parameter names,
it is good practice to define a constructor for any UDF that will be re-used and
shared. Since functions of an UDF are executed both on the main node and on
worker processes, UDFs will be pickled and sent over the network in the process.
To avoid transferring unwanted state or unnecessary or unpicklable member
variables, a clean copy of an UDF is created before transfer. This clean copy is
created by re-instantiating the UDF class with the parameters that were passed
to the constructor of the UDF base class.

That means a user-defined constructor has to fulfill two conditions:

1. It has to pass any parameters to the superclass constructor.
2. It has to accept exactly the parameters that it passed to the superclass
   constructor whenever a clean copy is created and behave the same way as in the
   original call.

That means modification of parameters other than assigning values for default parameters
can lead to complications and should be avoided, if possible.

.. testcode:: constructor

    class MyParamUDF(UDF):
        '''
        This UDF demonstrates how to define a constructor for a UDF.

        Parameters
        ----------
        my_param
            A parameter that is passed into the UDF for demonstration purposes.
            It is mirrored in the `demo` result buffer.
        '''
        def __init__(self, my_param=None):
            # Assigning a default parameter works,
            # provided the UDF accepts it in subsequent calls
            # exactly as it passed it to the superclass:
            if my_param is None:
                my_param = "Eins, zwei drei!"
            
            # !!!DON'T!!! self.my_param = my_param
            # !!!DON'T!!! super().__init__()
            # We have to pass it to the superclass constructor instead.
            # This makes sure it will be available via self.params and
            # on clean copies.

            # !!!DON'T!!! super().__init__(other_param=my_param)
            # This would trigger a TypeError when
            # MyParamUDF(other_param=my_param) is called for a clean copy.

            # DO pass all parameters to the superclass exactly as this
            # class would accept them
            super().__init__(my_param=my_param)

        # The rest of this UDF just mirrors the parameter
        # back in a result buffer for demonstration.

        def get_result_buffers(self):
            '''
            Declare a result buffer for the parameter
            '''
            return {
                'demo': self.buffer(
                    kind='single',
                    dtype=object
                ),
            }

        def process_frame(self, frame):
            '''
            We assign the parameter to the result buffer
            '''
            self.results.demo[:] = self.params.my_param
            
        def merge(self, dest, src):
            '''
            We pass through the result buffer
            '''
            dest.demo[:] = src.demo


    res = ctx.run_udf(dataset=dataset, udf=MyParamUDF())
    assert res['demo'].data[0] == "Eins, zwei drei!"

    res2 = ctx.run_udf(dataset=dataset, udf=MyParamUDF(my_param=42))
    assert res2['demo'].data[0] == 42


Initializing result buffers
~~~~~~~~~~~~~~~~~~~~~~~~~~~

To allow a UDF to initialize a result buffer to the correct values,
the method :meth:`~libertem.udf.base.UDFPreprocessMixin.preprocess`
can be implemented. It is run once per partition and assigning to
:code:`kind="nav"` result buffers will assign to the results of the
whole partition. See :code:`MaxUDF` above for an example.

.. _`progress bar`:

Running UDFs
------------

As shown in the examples above, the :meth:`~libertem.api.Context.run_udf` method
of :class:`~libertem.api.Context` is used to run UDFs. Usually, you only need to
pass an instance of your UDF and the dataset you want to run on:
   
.. testcode:: run

    udf = YourUDF(param1="value1")
    res = ctx.run_udf(udf=udf, dataset=dataset)

You can also enable a progress bar:

.. testcode:: run

    udf = YourUDF(param1="value1")
    res = ctx.run_udf(udf=udf, dataset=dataset, progress=True)

:meth:`~libertem.api.Context.run_udf` returns a :code:`dict`, having the buffer
names as keys (as defined in :meth:`~libertem.udf.UDF.get_result_buffers`) and
:class:`~libertem.common.buffers.BufferWrapper` instances as values. You
can use these in any place you would use a NumPy array, for example as an argument to
NumPy functions, or you can explicitly convert them to NumPy arrays by accessing
the :code:`.data` attribute, or by calling :meth:`numpy.array`:
   
.. testcode:: run

   import numpy as np

   res = ctx.run_udf(udf=udf, dataset=dataset)
   # convert to NumPy array, assuming we declared a buffer
   # with name `buf1`:
   arr = res['buf1'].data
   arr = np.array(res['buf1'])

   # or directly treat as array:
   np.sum(res['buf1'])

If you want to run multiple independent UDFs on a single :code:`DataSet`,
you can pass in a list of UDFs to :meth:`~libertem.api.Context.run_udf`. This can be faster
than making two passes over the whole :code:`DataSet`:

.. testcode:: run

   from libertem.udf.sum import SumUDF

   # results are returned as a tuple, so we can unpack them here into
   # `res` and `res_sum`:
   res, res_sum = ctx.run_udf(udf=[udf, SumUDF()], dataset=dataset)


.. _`udf roi`:

Regions of interest
~~~~~~~~~~~~~~~~~~~

In addition, you can pass the :code:`roi` (region of interest) parameter, to
run your UDF on a selected subset of data. :code:`roi` should be a NumPy array
containing a bool mask, having the shape of the navigation axes of the dataset.
For example, to process a random subset of a 4D-STEM dataset:
   
.. testcode:: run

   import numpy as np

   # If your dataset has shape `(14, 14, 32, 32)` with two signal
   # and two navigation dimensions, `dataset.shape.nav`
   # translates to `(14, 14)`.
   roi = np.random.choice(a=[False, True], size=dataset.shape.nav)
   ctx.run_udf(udf=udf, dataset=dataset, roi=roi)

Note that the result array only contains values for the selected indices, all
other indices are set to :code:`nan` (or, if the dtype doesn't support nan,
some other, not further defined value). It is best to limit further processing
to the same :code:`roi`.

You can also access a flat array that is not filled up with :code:`nan` using
:code:`.raw_data`:

.. testcode:: run

   res = ctx.run_udf(udf=udf, dataset=dataset, roi=roi)
   res['buf1'].raw_data

.. _plotting:

Live Plotting
-------------

.. versionadded:: 0.7.0

LiberTEM can display a live plot of the UDF results. In the most simple case,
this can be done by setting :code:`plots=True` in
:meth:`~libertem.api.Context.run_udf`. 

.. testsetup:: live

    from libertem.udf.sum import SumUDF
    udf = SumUDF()

.. testcode:: live

    ctx.run_udf(dataset=dataset, udf=udf, plots=True)

See the following items for a full demonstration, including setting up fully
customized plots. The API reference can be found in :ref:`viz reference`.

.. toctree::

    liveplotting

.. _partial:

Partial results
---------------

.. versionadded:: 0.7.0

Instead of only getting the whole result after the UDF has finished running, you
can also use :meth:`~libertem.api.Context.run_udf_iter` to get a generator for
partial results:

.. testsetup:: partial

    from libertem.udf.sum import SumUDF


.. testcode:: partial

    for udf_results in ctx.run_udf_iter(dataset=dataset, udf=SumUDF()):
        # ... do something interesting with `udf_results`:
        a = np.sum(udf_results.buffers[0]['intensity'])

    # after the loop, `udf_results` contains the final results as usual

While the UDF execution is running, the UDF object should not be modified since
that leads to undefined behavior. In particular, nested or concurrent execution
of the same UDF objects must be avoided since it modifies the buffers that are
allocated internally while a UDF is running.

Asynchronous execution
----------------------

It is also possible to integrate LiberTEM into an async script or application by
passing :code:`sync=False` to :meth:`~libertem.api.Context.run_udf_iter` or
:meth:`~libertem.api.Context.run_udf`:

.. Not run with docs-check since we can't easily test async code there...

.. code-block:: python

    async for udf_results in ctx.run_udf_iter(dataset=dataset, udf=udf, sync=False):
        # ... do something interesting with `udf_results`:
        a = np.sum(udf_results[0]['intensity'])

    # or the version without intermediate results:
    udf_results = await ctx.run_udf(dataset=dataset, udf=udf, sync=False)

See the items below for a more comprehensive demonstration and documentation:

.. toctree::

    asyncDask API
--------

LiberTEM can generate an efficient Dask array from data sets. See :ref:`daskarray` for details!

.. automodule:: libertem.contrib.daskadapter
   :members:
   :undoc-members:

Dask result buffer
------------------

.. autoclass:: libertem.executor.utils.dask_buffer.DaskResultBufferWrapper
   :members:
   :undoc-members:
.. _`viz reference`:

Visualization
-------------

Matplotlib
~~~~~~~~~~

.. automodule:: libertem.viz.mpl
   :members:

bqplot
~~~~~~

.. automodule:: libertem.viz.bqp
   :members:

GMS
~~~~~~

.. automodule:: libertem.viz.gms
   :members:

Base classes and functions
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: libertem.viz.base
   :members:
.. _`corrections`:

Corrections
===========

LiberTEM includes correction facilities to substract a dark frame, multiply with a gain map
and patch out defect pixels (dead/hot). These corrections are applied on-the-fly when running
a UDF, both in the GUI and via the Python API.

The following data set formats ship with some of
this correction data:

* FRMS6: dark frame is loaded from the first part of the data set
* SEQ: dark frame and gain map are loaded from MRC sidecar files
  :code:`<basename>.dark.mrc` and :code:`<basename>.gain.mrc`, and
  bad pixels are loaded from :code:`<basename>.Config.Metadata.xml`.

In the GUI, all corrections that are supplied by the data set will be applied. In the Python API,
the user can decide to pass their own corrections to apply, via the :code:`corrections` parameter
of :code:`Context.run()` and :code:`Context.run_udf()`. It expects a :class:`~libertem.corrections.CorrectionSet` object,
for example:

.. testsetup:: *

    import numpy as np
    from libertem.executor.inline import InlineJobExecutor
    from libertem.udf.sum import SumUDF
    from libertem import api

    ctx = api.Context(executor=InlineJobExecutor())
    data = np.random.random((16, 16, 32, 32)).astype(np.float32)
    dataset = ctx.load("memory", data=data, sig_dims=2)

.. testcode::

    from libertem.corrections import CorrectionSet
    import sparse

    # excluded pixels are passed as a sparse COO matrix, which can be built
    # in different ways, here is one way:
    excluded = np.zeros((32, 32), dtype=bool)
    excluded[5, 16] = 1
    excluded = sparse.COO(excluded)

    ctx.run_udf(udf=SumUDF(), dataset=dataset, corrections=CorrectionSet(
        dark=np.zeros((32, 32)),
        gain=np.ones((32, 32)),
        excluded_pixels=excluded,
    ))

It can also be empty to disable corrections:

.. testcode::

    ctx.run_udf(udf=SumUDF(), dataset=dataset, corrections=CorrectionSet())

.. autoclass:: libertem.corrections.CorrectionSet
Distributed caching layer
=========================

:code:`CachedDataSet` can be used to temporarily store your data on a
local fast storage device (i.e. SSD), if your source data is stored on a slower
device (i.e. NFS or a slower spinning disk).

.. automodule:: libertem.io.dataset.cached
   :members:
   :undoc-members:
Mask creation and manipulation
------------------------------

.. automodule:: libertem.masks
   :members:
   :undoc-members:
.. _`dataset api`:

Data Set API
============

This API allows to load and handle data on a distributed system efficiently. Note that you should
not directly use most dataset methods, but rather use the more high-level tools available, for
example user-defined functions.

See :ref:`our documentation on loading data <loading data>` for a high-level introduction.

.. _`formats`:

Formats
-------

.. _`mib`:

Merlin Medipix (MIB)
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: libertem.io.dataset.mib.MIBDataSet

.. _`raw binary`:

Raw binary files
~~~~~~~~~~~~~~~~

.. autoclass:: libertem.io.dataset.raw.RawFileDataSet

.. _`dm format`:

Digital Micrograph (DM3, DM4) files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: libertem.io.dataset.dm.DMDataSet

.. _`empad`:

EMPAD
~~~~~

.. autoclass:: libertem.io.dataset.empad.EMPADDataSet

.. _`k2is`:

K2IS
~~~~

.. autoclass:: libertem.io.dataset.k2is.K2ISDataSet

.. _`frms6`:

FRMS6
~~~~~

.. autoclass:: libertem.io.dataset.frms6.FRMS6DataSet

.. _`blo`:

BLO
~~~

.. autoclass:: libertem.io.dataset.blo.BloDataSet

.. _`ser`:

SER
~~~

.. autoclass:: libertem.io.dataset.ser.SERDataSet

.. _`hdf5`:

HDF5
~~~~

.. autoclass:: libertem.io.dataset.hdf5.H5DataSet

.. _`seq`:

Norpix SEQ
~~~~~~~~~~

.. autoclass:: libertem.io.dataset.seq.SEQDataSet

.. _`mrc`:

MRC
~~~

.. autoclass:: libertem.io.dataset.mrc.MRCDataSet

.. _`tvips`:

TVIPS
~~~~~

.. autoclass:: libertem.io.dataset.tvips.TVIPSDataSet

.. _`memory`:

Memory data set
~~~~~~~~~~~~~~~

.. autoclass:: libertem.io.dataset.memory.MemoryDataSet

.. _`daskds`:

Dask data set
~~~~~~~~~~~~~~~

.. autoclass:: libertem.io.dataset.dask.DaskDataSet

Internal DataSet API
--------------------

.. automodule:: libertem.io.dataset.base
   :members:
   :undoc-members:
Internal API
------------

MaskContainer
~~~~~~~~~~~~~

:class:`libertem.common.container.MaskContainer` helps to implement highly efficient
mask application operations, such as virtual detector, center of mass or feature
vector calculations.

.. versionchanged:: 0.4.0
    Moved from :mod:`libertem.job.masks` to :mod:`libertem.common.container` to
    use it in UDFs and prepare deprecation of the Job interface.

.. automodule:: libertem.common.container
   :members: MaskContainer

Shapes and slices
~~~~~~~~~~~~~~~~~

These classes help to manipulate shapes and slices of n-dimensional binary data to
facilitate the MapReduce-like processing of LiberTEM. See :ref:`concepts` for
a high-level introduction.

.. automodule:: libertem.common.shape
   :members:
   :undoc-members:

.. automodule:: libertem.common.slice
   :members:
   :undoc-members:

CPU and CUDA devices
~~~~~~~~~~~~~~~~~~~~

These methods get and set information that controls on which devices a
computation runs.

.. automodule:: libertem.common.backend
   :members:
   :undoc-members:

.. automodule:: libertem.utils.devices
   :members:
   :undoc-members:
.. _reference:

Python API reference
====================

.. toctree::
   :maxdepth: 2

   api
   udf
   dataset
   io_backends
   application
   masks
   internals
   dask
   cluster-ds
   caching
   corrections
   viz

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
.. _`application api`:

Application-specific APIs
=========================

.. note::
    Starting from 0.4.0, the application-specific code of LiberTEM will be
    spun out step by step as sub-packages that can be installed independent of
    LiberTEM core. See :ref:`packages` for the current overview of sub-packages.

.. toctree::
    :maxdepth: 2

    app/amorphous
    app/strain
    app/phasechange
    app/holography
.. _`udf reference`:

UDF API reference
-----------------

.. _`define udf ref`:

Defining UDFs
~~~~~~~~~~~~~

See :ref:`user-defined functions` for an introduction and in-depth explanation.

Mixins for processing methods
#############################

.. autoclass:: libertem.udf.base.UDFFrameMixin
    :members:

.. autoclass:: libertem.udf.base.UDFTileMixin
    :members:

.. autoclass:: libertem.udf.base.UDFPartitionMixin
    :members:

Base UDF class
##############

.. autoclass:: libertem.udf.base.UDF
    :members:

Meta information
################

.. autoclass:: libertem.udf.base.UDFMeta
    :members:

Pre- and postprocessing
#######################

.. autoclass:: libertem.udf.base.UDFPreprocessMixin
    :members:

.. autoclass:: libertem.udf.base.UDFPostprocessMixin
    :members:

.. _`run udf ref`:

Running UDFs
~~~~~~~~~~~~

Three methods of :class:`libertem.api.Context` are relevant for running user-defined functions:

.. autoclass:: libertem.api.Context
   :members: run_udf,run_udf_iter,map
   :noindex:

Result type for UDF result iterators:

.. autoclass:: libertem.udf.base.UDFResults
    :members:

.. _`buffer udf ref`:

Buffers
~~~~~~~

:class:`~libertem.common.buffers.BufferWrapper` objects are used to manage data in the context of user-defined functions.

.. automodule:: libertem.common.buffers
   :members:
   :undoc-members:

.. _`utilify udfs`:

Included utility UDFs
~~~~~~~~~~~~~~~~~~~~~

Some generally useful UDFs are included with LiberTEM:

.. note::
    See :ref:`application api` for application-specific UDFs and analyses.

.. _`sum udf`:

Sum of frames
#############

.. autoclass:: libertem.udf.sum.SumUDF
    :members:

.. _`logsum udf`:

Sum of log-scaled frames
########################

.. autoclass:: libertem.udf.logsum.LogsumUDF
    :members:

.. _`stddev udf`:

Standard deviation
##################

.. autoclass:: libertem.udf.stddev.StdDevUDF
    :members:

.. autofunction:: libertem.udf.stddev.run_stddev

.. autofunction:: libertem.udf.stddev.consolidate_result

.. _`sumsig udf`:

Sum per frame
#############

.. autoclass:: libertem.udf.sumsigudf.SumSigUDF
    :members:

.. _`masks udf`:

Apply masks
###########

.. autoclass:: libertem.udf.masks.ApplyMasksUDF
    :members:

.. _`pick udf`:

Load data
#########

.. autoclass:: libertem.udf.raw.PickUDF
    :members:

NoOp
####

.. autoclass:: libertem.udf.base.NoOpUDF
    :members:

.. _`io backends`:

I/O Backends
============

By default, on Windows, a buffered strategy is chosen, whereas Linux and Mac
OS use unbuffered memory-mapped I/O. You can pass an instance of
:class:`~libertem.io.dataset.base.backend.IOBackend` as the :code:`io_backend`
parameter of :meth:`~libertem.api.Context.load` to use a different backend.
This allows you to override the default backend choice, or set parameters
for the backend.

Note that some file formats can't support different I/O backends, such as HDF5
or SER, because they are implemented using third-party reading libraries
which perform their own I/O.

Available I/O backends
----------------------

BufferedBackend
~~~~~~~~~~~~~~~

.. autoclass:: libertem.io.dataset.base.BufferedBackend
    :noindex:

MmapBackend
~~~~~~~~~~~

.. autoclass:: libertem.io.dataset.base.MMapBackend
    :noindex:
 
DirectBackend
~~~~~~~~~~~~~

.. autoclass:: libertem.io.dataset.base.DirectBackend
    :noindex:
Cluster data set
================

:code:`ClusterDataSet` is a distributed data set, storing individual partitions
using raw binary files.

.. automodule:: libertem.io.dataset.cluster
   :members:
   :undoc-members:
LiberTEM Context API
--------------------

See :ref:`api documentation` for introduction and complete examples.

Context
~~~~~~~

.. automodule:: libertem.api
   :members:

Analysis API
~~~~~~~~~~~~

.. automodule:: libertem.analysis.base
   :members: Analysis, AnalysisResult, AnalysisResultSet

.. automodule:: libertem.analysis.masks
   :members: MasksResultSet, SingleMaskResultSet

.. automodule:: libertem.analysis.com
   :members: COMResultSet

.. automodule:: libertem.analysis.sum
   :members: SumResultSet

.. automodule:: libertem.analysis.raw
   :members: PickResultSet
Off-axis electron holography
============================

.. versionadded:: 0.3.0

The off-axis holography applications (see :ref:`holography app` for the application examples) are realized in
two modules: UDF for off axis electron holography reconstruction and utility function for hologram simulations.

Hologram reconstruction
-----------------------

The reconstruction module contains class for reconstruction of off-axis holograms using Fourier-space method
which implies following processing steps:

* Fast Fourier transform
* Filtering of the sideband in Fourier space and cropping (if applicable)
* Centering of the sideband
* Inverse Fourier transform.

.. automodule:: libertem.udf.holography
   :members:
   :undoc-members:

Hologram simulation
-------------------

.. automodule:: libertem.utils.generate
   :members: hologram_frame
.. _`phasechange`:

Phase-change materials reference
================================

.. note::

    See :ref:`phasechange app` for an overview and description of these UDFs.

Crystallinity map
~~~~~~~~~~~~~~~~~

.. automodule:: libertem.udf.crystallinity
   :members:
   :exclude-members: get_result_buffers, get_task_data
Amorphous materials reference
=============================

.. note::

    See :ref:`amorphous app` for an overview and description of the amorphous applications.

.. _`fem api`:

Fluctuation EM
--------------

This module contains the UDF for applying FEM to a single ring (mostly useful for interactive use).

.. automodule:: libertem.udf.FEM
   :members:
   :exclude-members: get_result_buffers, get_task_data

.. _`radial fourier api`:

Radial Fourier Analysis
-----------------------

This module contains the radial fourier series analysis, for analysing frequencies and
symmetries of diffraction patterns.

.. automodule:: libertem.analysis.radialfourier
   :members: RadialFourierResultSet
Correlation-based peak finding and strain mapping reference
===========================================================

.. _`blobfinder api`:

Blobfinder
----------

.. deprecated:: 0.4.0
    Blobfinder has moved to its own package LiberTEM-blobfinder with a new
    structure. Please see
    https://libertem.github.io/LiberTEM-blobfinder/index.html for installation
    instructions and documentation of the new structure. Imports from
    :code:`libertem.udf.blobfinder` trigger a :code:`FutureWarning` starting from
    0.4.0 and are supported until LiberTEM release 0.6.0.

.. _`matching api`:

Matching
--------

These modules contain classes and helper functions that extract and manipulate lattices from correlation results.

.. automodule:: libertem.analysis.gridmatching
   :members:
   :show-inheritance:

.. automodule:: libertem.analysis.fullmatch
   :members:
   :show-inheritance:
LiberTEM
========

**LiberTEM creators**: Alexander Clausen¹, Dieter Weber¹, Karina Ruzaeva¹, Vadim
Migunov¹², Jan Caron¹, Rahul Chandra³, Magnus Nord⁴, Colin Ophus⁵, Simon Peter,
Jay van Schyndel⁶, Jaeweon Shin⁷, Knut Müller-Caspary¹, Rafal Dunin-Borkowski¹

**Presentation**: Alexander Clausen¹, Dieter Weber¹

¹Jülich Research Centre, ²RWTH Aachen University, ³Chandigarh University,
⁴University of Antwerp, ⁵Lawrence Livermore National Lab, ⁶Monash University
eResearch Centre, ⁷ETH Zürich

Starting from 2009, increases in the data rates of TEM cameras have outpaced
increases in network, mass storage and memory bandwidth by two orders of
magnitude. This immense increase in performance opens new doors for applications
and, at the same time, requires adequate IT systems to handle acquisition,
storage and processing that go well beyond solutions based on personal
computers.

The LiberTEM open source platform is developed to address these challenges. It
allows high-throughput distributed processing of large-scale binary data sets
using a simplified MapReduce programming model. This programming model allows to
process live data streams and offline data with the same algorithms while
generating a stream of intermediate results for interactive GUI display. At the
same time, it scales to very high throughput. In a benchmark we have already
shown 46 GB/s on a microcloud cluster with eight nodes.

The current focus for application development is pixelated scanning transmission
electron microscopy (STEM) and scanning electron beam diffraction data.
Nevertheless, LiberTEM is suitable for processing any type of n-dimensional
numerical data.

LiberTEM can be used with a web-based GUI, via Python scripting, or embedded in
other applications as a data processing back-end. It is supported on Windows,
Linux and Mac OS X. Several typical TEM file formats are supported.

This contribution will first introduce the general architecture and interfaces,
demonstrate how it can be applied to typical tasks in electron microscopy, and
show examples for embedding in an application as a processing back-end.

More information on LiberTEM is available at
https://libertem.github.io/LiberTEM/index.html

Full continuously updated list of acknowledgments:
https://libertem.github.io/LiberTEM/acknowledgments.html
