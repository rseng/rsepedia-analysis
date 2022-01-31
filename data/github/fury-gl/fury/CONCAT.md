---
title: 'FURY: advanced scientific visualization'
tags:
  - Python
  - scientific visualization
  - 3D rendering
  - Shaders
  - GLSL
  - Vulkan
authors:
  - name: Eleftherios Garyfallidis
    orcid: 0000-0002-3991-2774
    affiliation: 1
  - name: Serge Koudoro
    orcid: 0000-0002-9819-9884
    affiliation: 1
  - name: Javier Guaje
    orcid: 0000-0003-3534-3794
    affiliation: 1
  - name: Marc-Alexandre Côté
    orcid: 0000-0002-5147-7859
    affiliation: 3
  - name: Soham Biswas
    orcid: 0000-0002-8449-2107
    affiliation: 4
  - name: David Reagan
    orcid: 0000-0002-8359-7580
    affiliation: 5
  - name: Nasim Anousheh
    orcid: 0000-0002-5931-7753
    affiliation: 1
  - name: Filipi Silva
    orcid: 0000-0002-9151-6517
    affiliation: 2
  - name: Geoffrey Fox
    orcid: 0000-0003-1017-1391
    affiliation: 1
  - name: FURY Contributors
    affiliation: 6
affiliations:
 - name: Department of Intelligent Systems Engineering, Luddy School of Informatics, Computing and Engineering, Indiana University, Bloomington, IN, USA
   index: 1
 - name: Network Science Institute, Indiana University, Bloomington, IN, USA
   index: 2
 - name: Microsoft Research, Montreal, Canada
   index: 3
 - name: Department of Computer Science and Engineering, Institute of Engineering and Management, Kolkata, India
   index: 4
 - name: Advanced Visualization Lab, University Information Technology Services, Indiana University, Bloomington, IN, USA
   index: 5
 - name: Anywhere in the Universe
   index: 6

date: 8 April 2021
bibliography: paper.bib
---



# Summary

Free Unified Rendering in pYthon (FURY), is a community-driven, open-source, and high-performance scientific visualization library that harnesses the graphics processing unit (GPU) for improved speed, precise interactivity, and visual clarity. FURY provides an integrated API in Python that allows UI elements and 3D graphics to be programmed together. FURY is designed to be fully interoperable with most projects of the Pythonic ecosystem that use NumPy [@harris2020array] for processing numerical arrays. In addition, FURY uses core parts of VTK [@schroeder1996visualization] and enhances them using customized shaders. FURY provides access to the latest technologies such as raytracing, signed distance functionality, physically based rendering, and collision detection for direct use in research. More importantly, FURY enables students and researchers to script their own 3D animations in Python and simulate dynamic environments.


# Statement of need

The massive amount of data collected and analyzed by scientists in several disciplines requires powerful tools and techniques able to handle these whilst still managing efficiently the computational resources available. In some particular disciplines, these datasets not only are large but also encapsulate the dynamics of their environment, increasing the demand for resources. Although 3D visualization technologies are advancing quickly [@sellers2016vulkan], their sophistication and focus on non-scientific domains makes it hard for researchers to use them.  In other words, most of the existing 3D visualization and computing APIs are low-level (close to the hardware) and made for professional specialist developers.  Because of these issues, there is a significant barrier to many scientists and these powerful technologies are rarely deployed to everyday research practices.

Therefore, FURY is created to address this necessity of high-performance 3D scientific visualization in an easy-to-use API fully compatible with the Pythonic ecosystem.

# FURY Architecture

FURY is built to be modular, scalable, and to respect software engineering principles including a well-documented codebase and unit integration testing. The framework runs in all major operating systems including multiple Linux distributions, Windows, and macOS. Also, it can be used on the desktop and the web. The framework contains multiple interconnected engines, modules, API managers as illustrated in \autoref{fig:architecture}.

![The FURY framework contains multiple interconnected engines to bring forward advanced visualization capabilities. Additionally, it contains an integrated user interface module and an extendable I/O module. One of the most important classes is the Scene Manager that connects the actors to the shaders, animations, and interactors for picking 3D objects. The actors are directly connected to NumPy arrays with vertices, triangles, and connectivity information that is provided by the core engine. These are then connected to the physics and networks  engines.\label{fig:architecture}](https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/fury_paper/architecture.png)


**Rendering Engine**: This engine includes managers like scene, animation, shader, and picking manager. The scene manager allows the visual objects to appear on a canvas. The picking manager allows selecting specific objects in the scene. The animation manager allows users to script their own 3D animations and videos with timelines allowing objects to act in specific times. Lastly, the shader manager provides several interfaces to different elements in the OpenGL rendering pipeline. This manager allows developers to add customized shaders snippets to the existing shaders included in the API.

**Core Engine**: This engine contains utilities for building actors from primitives and transforming them. A primitive is an object that describes its shape and connectivity with elements such as vertices and triangles.

**Physics Engine**: This engine allows us to either build collision mechanisms as used in molecular dynamics or integrate well-established engines such as Bullet [@coumans2013bullet] and NVIDIA PhysX [@harris2009cuda].

**Networks Engine**:  This engine allows for the creation and use of graph systems and layouts.

**Integrated User Interfaces Module**: FURY contains its own user interfaces. This module provides a range of UI 2D / 3D elements such as buttons, combo boxes, and file dialogues. Nevertheless, users can easily connect to other known GUIs such as Qt or IMGUI if necessary.

**I/O module**: FURY supports a range of file formats from the classic OBJ format to the more advanced GLTF format that can be used to describe a complete scene with many actors and animations.

**Interoperability**: FURY can be used together with projects such as SciPy [@virtanen2020scipy], Matplotlib [@hunter2007matplotlib], pandas [@mckinney2010data], scikit-learn [@pedregosa2011scikit], NetworkX [@hagberg2008exploring], PyTorch [@paszke2019pytorch] and TensorFlow [@abadi2016tensorflow].

FURY’s visualization API can be compared with VisPy [@campagnola2015vispy], glumpy [@rougier2015glumpy], Mayavi [@ramachandran2011mayavi], and others. VisPy and glumpy directly connect to OpenGL. FURY uses OpenGL through Python VTK, which can be advantageous because it can use the large stack of visualization algorithms available in VTK. This is similar to Mayavi, however, FURY provides an easy and efficient way to ease interaction with 3D scientific data via integrated user interface elements and allows to reprogram the low-level shaders for the creation of stunning effects (see \autoref{fig:features}) not available in VTK. Historically, FURY had also a different path than these libraries as it was originally created for heavy-duty medical visualization purposes for DIPY [@garyfallidis2014dipy]. As the project grew it spinned off as an independent project with applications across the domains of science and engineering including visualization of nanomaterials and robotics simulations.




![**Top**. Dynamic changes are shown as diffused waves on the surface of the horse visualized with FURY. Showing here 4 frames at 4 different time points (t1−t4). A vertex and fragment shader are used to calculate in real-time the mirroring texture and blend its colors with the blue-yellow wave. **Bottom**. In FURY we create actors that contain multiple visual objects controlled by NumPy arrays.  Here an actor is generating 5 superquadrics with different properties (e.g. colors, directions, metallicity) by injecting the information as NumPy arrays in a single call.  This is one of the important design choices that make FURY easier to use but also faster to render. Actors in FURY can contain many objects. The user can select any of the objects in the actor. Here the user selected the first object (spherical superquadric).\label{fig:features}](https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/fury_paper/features.png)



# Acknowledgements
FURY is partly funded through NSF #1720625 Network for Computational Nanotechnology - Engineered nanoBIO Node [@klimeck2008nanohub].


# References
<h1 align="center">
  <br>
  <a href="https://www.fury.gl"><img src="https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/fury-logo.png" alt="FURY" width="200"></a>
  <br>Free Unified Rendering in Python<br>

</h1>

<h4 align="center">A software library for scientific visualization in Python.
</h4>

<p align="center">
<a href="https://dev.azure.com/fury-gl/fury/_build/latest?definitionId=1&branchName=master"><img src="https://dev.azure.com/fury-gl/fury/_apis/build/status/fury-gl.fury?branchName=master">
</a>
<a href="https://pypi.python.org/pypi/fury"><img src="https://img.shields.io/pypi/v/fury.svg"></a>
<a href="https://anaconda.org/conda-forge/fury"><img src="https://anaconda.org/conda-forge/fury/badges/version.svg"></a>
<a href="https://codecov.io/gh/fury-gl/fury"><img src="https://codecov.io/gh/fury-gl/fury/branch/master/graph/badge.svg"></a>
<a href="https://app.codacy.com/app/fury-gl/fury?utm_source=github.com&utm_medium=referral&utm_content=fury-gl/fury&utm_campaign=Badge_Grade_Dashboard"><img src="https://api.codacy.com/project/badge/Grade/922600af9f94445ead5a12423b813576"></a>
<a href="https://doi.org/10.21105/joss.03384"><img src="https://joss.theoj.org/papers/10.21105/joss.03384/status.svg"></a>

</p>

<p align="center">
  <a href="#general-information">General Information</a> •
  <a href="#key-features">Key Features</a> •
  <a href="#installation">Installation</a> •
  <a href="#how-to-use">How to use</a> •
  <a href="#credits">Credits</a> •
  <a href="#contribute">Contribute</a> •
  <a href="#credits">Citing</a>
</p>

|         |         |         |
|:--------|:--------|:--------|
| <a href="#"><img src="https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/ws_smaller.gif" alt="FURY Networks" width="400px"></a> | <a href="#"><img src="https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/swarming_simulation.gif" alt="swarming simulation" width="400px"></a> | <a href="#"><img src="https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/shaders_horse.gif" alt="shaders horse" width="400px"></a> |
| *Network Visualization*          | *Swarming/flocking simulation based on simple boids rules*  |  *Easy shader effect integration.*  |
| <a href="#"><img src="https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/shaders_sdf.gif" alt="sdf" width="400px"></a>  | <a href="#"><img src="https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/collides_simulation.gif" alt="Collides simulation" width="400px"></a> | <a href="#"><img src="https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/physics_bricks_fast.gif" alt="Physics bricks" width="400px"></a> |
| *Ray Marching and Signed Distance Functions* | *Particle collisions* | *Interoperability with the [pyBullet](https://pybullet.org/wordpress/) library.*  |
| <a href="#"><img src="https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/ui_tab.gif" alt="UI Tabs" width="400px"></a>  | <a href="#"><img src="https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/shaders_dragon_skybox.gif" alt="Shaders dragon skybox" width="400px"></a>  | <a href="#"><img src="https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/picking_engine.gif" alt="Picking object" width="400px"></a> |
| *Custom User Interfaces* |  *Shaders and SkyBox integration*  | *Easy picking manager* |


# General Information

- **Website and Documentation:** https://fury.gl
- **Tutorials:** https://fury.gl/latest/auto_tutorials/index.html
- **Demos:** https://fury.gl/latest/auto_examples/index.html
- **Blog:**  https://fury.gl/latest/blog.html
- **Mailing list:** https://mail.python.org/mailman3/lists/fury.python.org
- **Official source code repo:** https://github.com/fury-gl/fury.git
- **Download releases:** https://pypi.org/project/fury/
- **Issue tracker:** https://github.com/fury-gl/fury/issues
- **Free software:** 3-clause BSD license
- **Community:** Come to chat on [Discord](https://discord.gg/6btFPPj)

# Key Features

- Custom User Interfaces
- Physics Engines API
- Custom Shaders
- Interactive local and Remote rendering in Jupyter Notebooks
- Large amount of Tutorials and Examples

# Installation

## Dependencies

FURY requires:

- Numpy (>=1.7.1)
- Vtk (>=8.1.2)
- Scipy (>=1.2.0)
- Pillow>=5.4.1

## Releases

`pip install fury` or `conda install -c conda-forge fury`

## Development

### Installation from source

**Step 1.** Get the latest source by cloning this repo:

    git clone https://github.com/fury-gl/fury.git

**Step 2.** Install requirements:

    pip install -r requirements/default.txt

**Step 3.** Install fury

As a [local project installation](https://pip.pypa.io/en/stable/reference/pip_install/#id44) using:

    pip install .

Or as an ["editable" installation](https://pip.pypa.io/en/stable/reference/pip_install/#id44) using:

    pip install -e .

**If you are developing fury you should go with editable installation.**

**Step 4:** Enjoy!

For more information, see also [installation page on fury.gl](https://fury.gl/latest/installation.html)

## Testing

After installation, you can install test suite requirements:

    pip install -r requirements/test.txt

And to launch test suite:

    pytest -svv fury


# How to use

There are many ways to start using FURY:

- Go to [Getting Started](https://fury.gl/latest/getting_started.html)
- Explore our [Tutorials](https://fury.gl/latest/auto_tutorials/index.html) or [Demos](https://fury.gl/latest/auto_examples/index.html).


# Credits

Please, go to [contributors page](https://github.com/fury-gl/fury/graphs/contributors) to see who have been involved in the development of FURY.


# Contribute

We love contributions!

You've discovered a bug or something else you want to change - excellent! Create an [issue](https://github.com/fury-gl/fury/issues/new)!

# Citing

If you are using FURY in your work then do cite [this paper](https://doi.org/10.21105/joss.03384). By citing FURY, you are helping sustain the FURY ecosystem.

    Garyfallidis, Eleftherios, Serge Koudoro, Javier Guaje, Marc-Alexandre Côté, Soham Biswas,
    David Reagan, Nasim Anousheh, Filipi Silva, Geoffrey Fox, and Fury Contributors.
    "FURY: advanced scientific visualization." Journal of Open Source Software 6, no. 64 (2021): 3384.
    https://doi.org/10.21105/joss.03384


```css
    @article{Garyfallidis2021,
        doi = {10.21105/joss.03384},
        url = {https://doi.org/10.21105/joss.03384},
        year = {2021},
        publisher = {The Open Journal},
        volume = {6},
        number = {64},
        pages = {3384},
        author = {Eleftherios Garyfallidis and Serge Koudoro and Javier Guaje and Marc-Alexandre Côté and Soham Biswas and David Reagan and Nasim Anousheh and Filipi Silva and Geoffrey Fox and Fury Contributors},
        title = {FURY: advanced scientific visualization},
        journal = {Journal of Open Source Software}
    }
```
# pip requirements files

## Index

- [default.txt](default.txt) Default requirements
- [docs.txt](docs.txt) Documentation requirements
- [optional.txt](optional.txt) Optional requirements
- [test.txt](test.txt) Requirements for running test suite

## Examples

### Installing requirements

```bash
$ pip install -U -r requirements/default.txt
$ pip install -U -r requirements/optional.txt
```

or 

```bash
conda install --yes --file=requirements/default.txt --file=requirements/optional.txt
```

### Running the tests

```bash
$ pip install -U -r requirements/default.txt
$ pip install -U -r requirements/test.txt
```

or 

```bash
conda install --yes --file=requirements/default.txt --file=requirements/test.txt
```

### Running the Docs

```bash
$ pip install -U -r requirements/default.txt
$ pip install -U -r requirements/optional.txt
$ pip install -U -r requirements/docs.txt
```
---
name: Generic issue template
about: Describe this issue template's purpose here.
title: ''
labels: ''
assignees: ''

---

## Description
[Please provide a general introduction to the issue/proposal.]

[If reporting a bug, attach the entire traceback from Python and follow the way to reproduce below]
[If proposing an enhancement/new feature, provide links to related articles, reference examples, etc.]


## Way to reproduce
[If reporting a bug, please include the following important information:]
- [ ] Code example
- [ ] Relevant images (if any)
- [ ] Operating system and versions (run `python -c "from fury import get_info; print(get_info())"`)
---
name: GSOC request
about: Suggest an idea for a GSOC project
title: "[GSOC]"
labels: ":rocket: :snake: GSOC 2020"
assignees: ''

---

**Please describe your idea.**
A clear and concise description of what your idea is.  [...]

**Describe the level needed**
A clear and concise description of what you want to happen.

**Requested Mentor**
Cite or tag the desired mentor for this project.

**Technologies**
List the technologies needed for this idea.
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
# Documentation Generation

## Index

-   ``source``: Contains main *.rst files
-   ``tutorials``: python script describe how to use the api
-   ``examples``: Fury app showcases
-   ``build``: Contains the generated documentation

## Doc generation steps:

### Installing requirements

```bash
$ pip install -U -r requirements/default.txt
$ pip install -U -r requirements/optional.txt
$ pip install -U -r requirements/docs.txt
```

### Generate all the Documentation

Go to the `docs` folder and run the following commands to generate it.

#### Under Linux and OSX

```bash
$ make -C . clean && make -C . html
```

To generate the documentation without running the examples:

```bash
$ make -C . clean && make -C . html-no-examples
```
#### Under Windows

```bash
$ make clean
$ make html
```

To generate the documentation without running the examples:

```bash
$ make clean
$ make html-no-examples
```============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/fury-gl/fury/issues.

If you are reporting a bug, please include:

* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "feature"
is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

FURY could always use more documentation, whether
as part of the official FURY docs, in docstrings,
or even on the web in blog posts, articles, and such.
FURY uses [Sphinx](http://www.sphinx-doc.org/en/stable/index.html) to generate documentation.
Please follow the [numpy coding style](https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard) - and of course - [PEP8](https://www.python.org/dev/peps/pep-0008/)
for docstring documentation.



Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/fury-gl/fury/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `FURY` for local development.

1. Fork the `FURY` repo on GitHub.
2. Clone your fork locally::

    $ git clone https://github.com/your_name_here/fury.git

3. Add a tracking branch which can always have the last version of `FURY`::

    $ git remote add fury-gl https://github.com/fury-gl/fury.git
    $ git fetch fury-gl
    $ git branch fury-gl-master --track fury-gl/master
    $ git checkout fury-gl-master
    $ git pull

4. Create a branch from the last dev version of your tracking branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

5. Install it locally::

    $ pip install --user -e .

6. Now you can make your changes locally::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Install the required packages for running the unittests::

    $ pip install -r requirements/optional.txt
    $ pip install -r requirements/test.txt

8. When you're done making changes, check that your changes pass flake8 and pytest::

    $ flake8 fury
    $ pytest -svv fury

   To get flake8 and pytest, just pip install them into your virtualenv.

9. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.md.
3. The pull request should work for Python 3.6, 3.7, 3.8, 3.9 and for PyPy. Check
   https://github.com/fury-gl/fury/actions
   and make sure that the tests pass for all supported Python versions.

Publishing Releases
--------------------

Checklist before Releasing
~~~~~~~~~~~~~~~~~~~~~~~~~~

* Review the open list of `FURY issues <https://github.com/fury-gl/fury/issues>`_.  Check whether there are
  outstanding issues that can be closed, and whether there are any issues that
  should delay the release.  Label them !

* Check whether there are no build failing on `Travis`.

* Review and update the release notes.  Review and update the :file:`Changelog`
  file.  Get a partial list of contributors with something like::

      git shortlog -nse v0.1.0..

  where ``v0.1.0`` was the last release tag name.

  Then manually go over ``git shortlog v0.1.0..`` to make sure the release notes
  are as complete as possible and that every contributor was recognized.

* Use the opportunity to update the ``.mailmap`` file if there are any duplicate
  authors listed from ``git shortlog -ns``.

* Add any new authors to the ``AUTHORS`` file.

* Check the copyright years in ``docs/source/conf.py`` and ``LICENSE``

* Generate release notes. Go to ``docs/source/ext`` and run ``github_tools.py`` script the following way::

    $ python github_tools.py --tag=v0.1.0 --save --version=0.2.0

  This command will generate a new file named ``release0.2.0.rst`` in ``release_notes`` folder.

* Check the examples and tutorial - we really need an automated check here.

* Make sure all tests pass on your local machine (from the ``<fury root>`` directory)::

    cd ..
    pytest -s --verbose --doctest-modules fury
    cd fury # back to the root directory

* Check the documentation doctests::

    cd docs
    make -C . html
    cd ..

* The release should now be ready.

Doing the release
~~~~~~~~~~~~~~~~~

* Update release-history.rst in the documentation if you have not done so already.
  You may also highlight any additions, improvements, and bug fixes.

* Type git status and check that you are on the master branch with no uncommitted code.

* Now it's time for the source release. Mark the release with an empty commit, just to leave a marker.
  It makes it easier to find the release when skimming through the git history::

    git commit --allow-empty -m "REL: vX.Y.Z"

* Tag the commit::

    git tag -am 'Second public release' vX.Y.Z  # Don't forget the leading v

  This will create a tag named vX.Y.Z. The -a flag (strongly recommended) opens up a text editor where
  you should enter a brief description of the release.

* Verify that the __version__ attribute is correctly updated::

    import fury
    fury.__version__  # should be 'X.Y.Z'

  Incidentally, once you resume development and add the first commit after this tag, __version__ will take
  on a value like X.Y.Z+1.g58ad5f7, where +1 means “1 commit past version X.Y.Z” and 58ad5f7 is the
  first 7 characters of the hash of the current commit. The letter g stands for “git”. This is all managed
  automatically by versioneer and in accordance with the specification in PEP 440.

* Push the new commit and the tag to master::

    git push origin master
    git push origin vX.Y.Z

* Register for a PyPI account and Install twine, a tool for uploading packages to PyPI::

    python3 -m pip install --upgrade twine

* Remove any extraneous files::

    git clean -dfx

  If you happen to have any important files in your project directory that are not committed to git,
  move them first; this will delete them!

* Publish a release on PyPI::

    python setup.py sdist
    python setup.py bdist_wheel
    twine upload dist/*


* Check how everything looks on pypi - the description, the packages.  If
  necessary delete the release and try again if it doesn't look right.

* Set up maintenance / development branches

  If this is this is a full release you need to set up two branches, one for
  further substantial development (often called 'trunk') and another for
  maintenance releases.

  * Branch to maintenance::

      git co -b maint/X.Y.Z


    Push with something like ``git push upstream-rw maint/0.6.x --set-upstream``

  * Start next development series::

      git co main-master


    Next merge the maintenace branch with the "ours" strategy.  This just labels
    the maintenance branch `info.py` edits as seen but discarded, so we can
    merge from maintenance in future without getting spurious merge conflicts::

       git merge -s ours maint/0.6.x

    Push with something like ``git push upstream-rw main-master:master``

  If this is just a maintenance release from ``maint/0.6.x`` or similar, just
  tag and set the version number to - say - ``0.6.2.dev``.

* Push the tag with ``git push upstream-rw 0.6.0``

Other stuff that needs doing for the release
============================================

* Checkout the tagged release, build the html docs and upload them to
  the github pages website::

    make upload

* Announce to the mailing lists.  With fear and trembling.========
Credits
========

Core Development Team
---------------------

* Eleftherios Garyfallidis, Indiana University, IN, USA
* Marc-Alexandre Côté, Microsoft Research, Montreal, QC, CA
* David Reagan, Indiana University, IN, USA
* Ranveer Aggarwal, Microsoft, Hyderabad, Telangana, India
* Serge Koudoro, Indiana University, IN, USA
* Karandeep Singh Juneja, India

Contributors
------------

* Javier Guaje
* Soham Biswas
* Lenix Lobo
* Nasim Anousheh
* Ariel Rokem, University of Washington, WA, USA
* Matthew Brett, Birmingham University, Birmingham, UK
* Charles Poirier
* Yaroslav Halchenko
* Kesshi jordan
* Prashil
* Matthew Brett
* Kevin Sitek
* Bago Amirbekian
* Omar Ocegueda
* Gregory R. Lee
* Stefan van der Walt
* Enes Albay
* Shahnawaz Ahmed
* Jon Haitz Legarreta Gorroño
* Guillaume Favelier
* Etienne St-Onge
* Bishakh Ghosh
* Bago Amirbekian
* Christopher Nguyen
* Alexandre Gauvin
* Jiri Borovec
* Scott Trinkle
* Jakob Wasserthal
* Etienne St-Onge
* Gauvin Alexandre
* Ian Nimmo-Smith
* Yaroslav Halchenko
* Vivek Choudhary
* Saransh Jain
* Shreyas Bhujbal
* Ibrahim Anis
* Filipi Nascimento Silva
* Melina Raglin
* Tushar
* Naman Bansal
* Gottipati Gautam
* Liam Donohue
* Devanshu Modi
* Sanjay Marreddi
* Sajag Swami
* aju tamang
* Aman Soni
* Pietro Astolfi
* Haran.. _fury_pybullet:

FURY - pyBullet Integration Guide
=================================

.. contents::
    :local:
    :depth: 3

**Official docs:**
  * `FURY <https://fury.gl/latest/reference/index.html>`__

  * `pyBullet <https://docs.google.com/document/d/10sXEhzFRSnvFcl3XxNGhnD4N2SedqwdAvK3dsihxVUA/edit#>`__

.. note:: All elements are in SI units.

Simple Rigid Body Dynamics
**************************

Necessary Imports
-----------------
The following imports are necessary for physics simulations:

+-----------------------+---------------------------------------------------------------+
|        Imports        |         Usage                                                 |
+=======================+===============================================================+
|         Numpy         |  Creation of arrays and conversion of radians to degrees.     |
+-----------------------+---------------------------------------------------------------+
|         FURY          |  Window and Actor API is used to visualize the simulation.    |
+-----------------------+---------------------------------------------------------------+
|         pyBullet      |  Physics simulation.                                          |
+-----------------------+---------------------------------------------------------------+
|         Itertools     |  The Counter iterator for keeping track of simulation steps.  |
+-----------------------+---------------------------------------------------------------+


.. code-block:: python

  import numpy as np
  from fury import window, actor
  import itertools
  import pybullet as p

Connection Mode
---------------

*"After importing the PyBullet module, the first thing to do is 'connecting' to the physics simulation. PyBullet is designed around a client-server driven API, with a client sending commands and a physics server returning the status. PyBullet has some built-in physics servers: DIRECT and GUI."*

In our case we use **DIRECT** connection as the visualization will be handled by FURY.

.. code-block:: python

  client = p.connect(p.DIRECT)

.. note:: Keeping track of physics client ID is optional unless multiple physics clients are used. In order to observe the same simulation in pybullet, replace ``p.DIRECT`` with ``p.GUI``.

Disconnection
-------------

PyBullet Physics client can be shutdown by the following command:

.. code-block:: python

  p.disconnect()

Setting Gravity
---------------

Global :py:class:`.Scene()` gravity can be set using the following command:

.. code-block:: python

  # Gravity vector.
  gravity_x = 0
  gravity_y = 0
  gravity_z = -10
  p.setGravity(gravity_x, gravity_y, gravity_z)

Creating Objects
----------------

The following criterion must be fulfilled in order to create an object which is in sync with both FURY and pyBullet:


+-----------------------+----------------------------------------------------------------------+
|       Object Actor    |         The actor which will be rendered by FURY                     |
+-----------------------+----------------------------------------------------------------------+
|      Collision Shape  |  The shape used by pybullet for collision simulations.               |
|                       |  **Optional** if collision simulation is not required.               |
+-----------------------+----------------------------------------------------------------------+
|       Multi-Body      |  The object that will be tracked by pybullet for general simulations.|
+-----------------------+----------------------------------------------------------------------+

The following is a snippet for creating a spherical ball of radius = 0.3

.. code-block:: python

  ###### Creating BALL
  # Ball actor
  ball_actor = actor.sphere(centers = np.array([[0, 0, 0]]),
                            colors=np.array([1,0,0]),
                            radii=0.3)

  # Collision shape for the ball.
  ball_coll = p.createCollisionShape(p.GEOM_SPHERE,
                                     radius=0.3)

  # Creating a Multibody which will be tracked by pybullet.
  ball = p.createMultiBody(baseMass=3,
                           baseCollisionShapeIndex=ball_coll,
                           basePosition=[2, 0, 1.5],
                           baseOrientation=[ 0, 0, 0, 1 ])

.. warning:: Centers for the actor must be set to ``(0, 0, 0)`` or else the simulation will be offset by that particular value.

Changing Object Dynamics
------------------------

Object dynamics such as mass, lateral_friction, damping, inertial_pos, inertial_orn, restitution, rolling friction etc can be changed. The following snippet shows how to change the lateral_friction and coeff of restitution of the same ball:

.. code-block:: python

  p.changeDynamics(ball, -1, lateralFriction=0.3, restitution=0.5)

.. note:: The second parameter is ``linkIndex`` which is for bodies having multiple links or joints. Passing -1 means applying changes to the base object.

Adding objects to the scene
---------------------------

Objects can be added simply by adding their respective actors to the scene.

.. code-block:: python

  scene = window.Scene()
  scene.add(ball_actor)

Application of Force/Torque
---------------------------

External force or torque to a body can be applied using applyExternalForce and applyExternalTorque. For e.g

.. code-block:: python

  p.applyExternalForce(ball, -1,
                       forceObj=[-2000, 0, 0],
                       posObj=ball_pos,
                       flags=p.WORLD_FRAME)

Here, the first argument refers to the object, the second one refers to the link, ``forceObj`` = force vector, ``posObj`` = Position Vector of the application of force. [Not applicable for ``applyExternalTorque``].

.. code-block:: python

  p.applyExternalTorque(ball, -1,
                       forceObj=[-2000, 0, 0],
                       flags=p.WORLD_FRAME)

Enabling collision
------------------

By default, collision detection is enabled between different dynamic moving bodies. The following snippet can be used to enable/disable collision explicitly between a pair of objects.

.. code-block:: python

  enableCol = 1
  p.setCollisionFilterPair(ball, brick, -1, -1, enableCol)

Here, we enable the collision between a ball and a brick object.

Creation of Show Manager
------------------------

A ``window.ShowManager`` and ``itertools.count`` instance must be created before defining the timer callback function and setting it to initialize.

.. code-block:: python

  # Create a show manager.
  showm = window.ShowManager(scene,
                          size=(900, 768), reset_camera=False,
                          order_transparent=True)
  showm.initialize()
  # Counter iterator for tracking simulation steps.
  counter = itertools.count()

Syncing properties of actors
----------------------------

The position and orientation of the actors in FURY can be updated by the values generated in pybullet during simulation. The following snippet updates all required parameters.

.. code-block:: python

  # Get the position and orientation of the ball.
  ball_pos, ball_orn = p.getBasePositionAndOrientation(ball)

  # Set position and orientation of the ball.
  ball_actor.SetPosition(*ball_pos)
  orn_deg = np.degrees(p.getEulerFromQuaternion(ball_orn))
  ball_actor.SetOrientation(*orn_deg)

``ball`` and ``ball_actor`` can be replaced by the appropriate object and actor.

Creation of Timer Callback
--------------------------

To simulate physics we need to call ``p.stepSimulation()`` in order to simulate a single step of physics simulation. Therefore, in order to update actors and simulate steps at each interval, we need to create a timer callback. At this point one can perform any operation that they feel like during each step of the simulation. This is also the appropriate section for the user to define all syncing activities required by the actors and render the scene accordingly. The following can be an example snippet:

.. code-block:: python

  # Counter iterator for tracking simulation steps.
  counter = itertools.count()

  # Variable for tracking applied force.
  apply_force = True

  # Create a timer callback which will execute at each step of simulation.
  def timer_callback(_obj, _event):
      global apply_force
      cnt = next(counter)
      showm.render()
      # Get the position and orientation of the ball.
      ball_pos, ball_orn = p.getBasePositionAndOrientation(ball)

      # Apply force for 5 times for the first step of simulation.
      if apply_force:
          # Apply the force.
          p.applyExternalForce(ball, -1,
                                forceObj=[-2000, 0, 0],
                                posObj=ball_pos,
                                flags=p.WORLD_FRAME)
          apply_force = False

      # Set position and orientation of the ball.
      ball_actor.SetPosition(*ball_pos)
      orn_deg = np.degrees(p.getEulerFromQuaternion(ball_orn))
      ball_actor.SetOrientation(*orn_deg)
      ball_actor.RotateWXYZ(*ball_orn)

      # Simulate a step.
      p.stepSimulation()

      # Exit after 2000 steps of simulation.
      if cnt == 2000:
          showm.exit()

  # Add the timer callback to showmanager.
  # Increasing the duration value will slow down the simulation.
  showm.add_timer_callback(True, 10, timer_callback)

Initiating the simulation
-------------------------

Once everything is set up, one can execute ``showm.start()`` to start the simulation.

Rendering multiple objects by a single actor
--------------------------------------------

Rendering multiple similar objects by a single actor is possible by manually updating the vertices of the individual objects. The said procedure will be demonstrated with the help of the brick wall simulation example where each brick is rendered by a single actor.
Firstly, we need to define the following parameters:

+-------------------------+-----------------------+-------------------------------------------------------------------------+
|         Variable        |        Shape          |                             Description                                 |
+=========================+=======================+=========================================================================+
|    nb_objects           |     1, 1              |   Number of objects to be rendered                                      |
+-------------------------+-----------------------+-------------------------------------------------------------------------+
|    object_centers       |     nb_objects, 3     |   To keep track of the centers in the xyz coordinate system. [x, y, z]  |
+-------------------------+-----------------------+-------------------------------------------------------------------------+
|    object_directions    |     nb_objects, 3     |   Array to track directions.                                            |
+-------------------------+-----------------------+-------------------------------------------------------------------------+
|    object_orientations  |     nb_objects, 4     |   Array to track orientations in quaternions. [x, y, z, w]              |
+-------------------------+-----------------------+-------------------------------------------------------------------------+
|    object_colors        |     nb_bricks, 3      |   Array to track colors.                                                |
+-------------------------+-----------------------+-------------------------------------------------------------------------+
|    object_collision     |     1, 1              |   Collision shape of the objects.                                       |
+-------------------------+-----------------------+-------------------------------------------------------------------------+

.. warning:: ``object_directions`` & ``object_orientations`` must be updated together or else orientation of objects in both the worlds may not be in sync.

Once we are ready with the above variables and array, we can proceed further to render the objects both in the FURY and pybullet world:

Rendering objects in FURY
~~~~~~~~~~~~~~~~~~~~~~~~~

To render objects in the FURY world we simply call the respective actors. For this example we call ``actor.box`` for rendering the bricks:

.. code-block:: python

  brick_actor_single = actor.box(centers=brick_centers,
                              directions=brick_directions,
                              scales=brick_sizes,
                              colors=brick_colors)

  scene.add(brick_actor_single)

Render Pybullet Objects
~~~~~~~~~~~~~~~~~~~~~~~

Now to render pybullet objects we simply create a list of multibodies:

.. code-block:: python

  bricks[i] = p.createMultiBody(baseMass=0.5,
                                baseCollisionShapeIndex=brick_coll,
                                basePosition=center_pos,
                                baseOrientation=brick_orn)

Syncing objects
~~~~~~~~~~~~~~~

Now in order to calculate and the vertices we execute the following snippet:

.. code-block:: python

  vertices = utils.vertices_from_actor(brick_actor_single)
  num_vertices = vertices.shape[0]
  num_objects = brick_centers.shape[0]
  sec = int(num_vertices / num_objects)

+-------------------+---------------------------------------------------------+
|      Vertices     |      Array storing vertices of all the objects.         |
+===================+=========================================================+
|    num_vertices   |  Number of vertices required to render the objects.     |
+-------------------+---------------------------------------------------------+
|    num_objects    |  Number of objects rendered                             |
+-------------------+---------------------------------------------------------+
|       sec         |  Number of vertices required to render a single object. |
+-------------------+---------------------------------------------------------+


Now the pybullet and FURY objects can be synced together by the following snippet:

.. code-block:: python

  def sync_brick(object_index, multibody):
    pos, orn = p.getBasePositionAndOrientation(multibody)

    rot_mat = np.reshape(
        p.getMatrixFromQuaternion(
            p.getDifferenceQuaternion(orn, brick_orns[object_index])),
        (3, 3))

    vertices[object_index * sec: object_index * sec + sec] = \
        (vertices[object_index * sec: object_index * sec + sec] -
        brick_centers[object_index])@rot_mat + pos

    brick_centers[object_index] = pos
    brick_orns[object_index] = orn


In order to Sync correctly, we do the following:

#. First we get the current position and orientation of the objects in the pybullet world with the help of ``p.getBasePositionAndOrientation``.
#. Then we calculate the difference between two quaternions using ``p.getDifferenceFromQuarternion``.
#. The said difference is then passed to ``p.getMatrixFromQuaternion`` to calculate the rotation matrix.
#. Now the method returns a tuple of size 9. Therefore we finally need to reshape the said tuple into a 3x3 matrix with the help of ``np.reshape``.
#. Next, we slice the necessary part of the vertices which render our desired object.
#. Then we bring it back to the origin by subtracting their centers.
#. After that we perform matrix multiplication of the rotation matrix and the vertices to orient the object.
#. After orientation we bring the object to its new position.
#. Finally we update the centers and the orientation of the object.

Lastly, we call this function in our timer callback to sync the objects correctly.

.. note:: VTK has an in-built method to handle gimbal locks therefore using ``actor.SetOrientation`` may lead to unwanted spinning simulations each time a gimbal lock is experienced. Hence, it is always advisable to use vertices and its corresponding rotation matrix to set the orientation.

Rendering Joints
----------------

.. image:: https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/physics_joints.png
    :align: center

A simulated robot as described in a URDF file has a base, and optionally links connected by joints. Each joint connects one parent link to a child link. At the root of the hierarchy there is a single root parent that we call base. The base can be either fully fixed, 0 degrees of freedom, or fully free, with 6 degrees of freedom. Since each link is connected to a parent with a single joint, the number of joints is equal to the number of links. Regular links have link indices in the range ``[0..getNumJoints()]`` Since the base is not a regular 'link', we use the convention of -1 as its link index. We use the convention that joint frames are expressed relative to the parent center of mass inertial frame, which is aligned with the principal axis of inertia. To know more how joints are implemented in pybullet refer the official docs.

We can create and sync joints in pybullet and FURY by following a few simple steps:

Firstly, in order to create objects with multiple joints we need to keep track of the following parameters:


+-----------------------------+--------------------+------------------------------------------+
|     Vertices                |      Shape         |             Description                  |
+=============================+====================+==========================================+
|     nb_links                |       1,1          |  Number of links to be rendered.         |
+-----------------------------+--------------------+------------------------------------------+
|    link_masses              |     nb_links       |  Masses of the links.                    |
+-----------------------------+--------------------+------------------------------------------+
|  linkCollisionShapeIndices  |     nb_links       |  Array tracking the collision shape IDs. |
+-----------------------------+--------------------+------------------------------------------+
|  linkVisualShapeIndices     |     nb_links       |  Optional as we won't be using           |
|                             |                    |  pybullet’s GUI render.                  |
+-----------------------------+--------------------+------------------------------------------+
|  linkPositions              |     nb_links, 3    |  Position of the links in [x, y, z].     |
+-----------------------------+--------------------+------------------------------------------+
|  linkOrientations           |     nb_links, 4    |  Orientation of the links in             |
|                             |                    |  [x, y, z, w].                           |
+-----------------------------+--------------------+------------------------------------------+
|  linkInertialFramePositions |     nb_links, 3    |  Position of the inertial frame of the   |
|                             |                    |  links.                                  |
+-----------------------------+--------------------+------------------------------------------+
|  linkInertialFrameOrns      |     nb_links, 4    |  Orientation of the inertial frame of    |
|                             |                    |  the links.                              |
+-----------------------------+--------------------+------------------------------------------+
|  indices                    |     nb_link        |  Link ID each corresponding link is      |
|                             |                    |  supposed to attach at.                  |
+-----------------------------+--------------------+------------------------------------------+
|  jointTypes                 |     nb_link        |  The type of joint between the links.    |
|                             |                    |  Multiple joint types are available.     |
+-----------------------------+--------------------+------------------------------------------+
|  axis                       |     nb_links, 3    |  The axis at which each link is supposed |
|                             |                    |  to rotate.                              |
+-----------------------------+--------------------+------------------------------------------+
|  linkDirections             |     nb_links, 3    |  Direction vector required to render     |
|                             |                    |  links in FURY.                          |
+-----------------------------+--------------------+------------------------------------------+

Extra Arrays such as ``linkHeights``, ``linkRadii`` etc may be required based on the link shape.
**Base link** is rendered separately, hence the above parameters should not contain information about the base link.

Now separately create definitions for the base link using the following parameters. Once we are ready with the required link parameters and definition, we can create a multibody to be rendered in the pybullet world. We can do so using ``p.createMultiBody``. Here’s a snippet:

.. code-block:: python

  rope = p.createMultiBody(base_mass,
                     	   base_shape,
                     	   visualShapeId,
                     	   basePosition,
                     	   baseOrientation,
                     	   linkMasses=link_Masses,
                          linkCollisionShapeIndices=linkCollisionShapeIndices,
                     	   linkVisualShapeIndices=linkVisualShapeIndices,
                     	   linkPositions=linkPositions,
                     	   linkOrientations=linkOrientations,
              	          linkInertialFramePositions=linkInertialFramePositions,
                 	    linkInertialFrameOrientations=linkInertialFrameOrns,
                     	   linkParentIndices=indices,
                     	   linkJointTypes=jointTypes,
                     	   linkJointAxis=axis)

Once we are done with the multibody we can create the actor to render the links:

.. code-block:: python

  rope_actor = actor.cylinder(centers=linkPositions,
                        directions=linkDirections,
                        colors=np.random.rand(n_links, 3),
                        radius=radii,
                        heights=link_heights, capped=True)

We can sync the joints using the following code snippet:

.. code-block:: python

  def sync_joints(actor_list, multibody):
    for joint in range(p.getNumJoints(multibody)):
        pos, orn = p.getLinkState(multibody, joint)[4:6]

        rot_mat = np.reshape(
        	p.getMatrixFromQuaternion(
            	p.getDifferenceQuaternion(orn, linkOrientations[joint])),
        	(3, 3))

    	vertices[joint * sec: joint * sec + sec] =\
        	(vertices[joint * sec: joint * sec + sec] -
         	linkPositions[joint])@rot_mat + pos

    	linkPositions[joint] = pos
    	linkOrientations[joint] = orn

Here, we determine the total number of joints using ``p.getNumJoints`` and run a loop to iterate through all the joints present within the object. Once we get access to a particular joint we use the ``p.getLinkState`` to get various information about a particular joint. Within the list of information we have access to positions and orientation of the joints at index 4 and 5. So we perform the query to get the position and orientation of the joints. After that the process of translation and rotation are the same as shown here.

---------------------

Examples
********

Brick Wall Simulation
---------------------

.. image:: https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/physics_bricks_multi_actor.gif
    :align: center

The code for the above simulation can be found `here <https://github.com/fury-gl/fury/blob/master/docs/examples/physics_using_pybullet/viz_brick_wall.py>`__.

Ball Collision Simulation
-------------------------

.. image:: https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/physics_collision.gif
    :align: center

The code for the above simulation can be found `here <https://github.com/fury-gl/fury/blob/master/docs/examples/physics_using_pybullet/viz_ball_collide.py>`__.

Brick Wall Simulation(Single Actor)
-----------------------------------

.. image:: https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/physics_bricks_fast.gif
    :align: center

The code for the above simulation can be found `here <https://github.com/fury-gl/fury/blob/master/docs/examples/physics_using_pybullet/viz_brick_wall.py>`__.

Chain Simulation
----------------

.. image:: https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/physics_chain.gif
    :align: center

The code for the above simulation can be found `here <https://github.com/fury-gl/fury/blob/master/docs/examples/physics_using_pybullet/viz_chain.py>`__.

Wrecking Ball Simulation
------------------------

.. image:: https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/physics_wrecking_ball.gif
    :align: center

The code for the above simulation can be found `here <https://github.com/fury-gl/fury/blob/master/docs/examples/physics_using_pybullet/viz_wrecking_ball.py>`__.

Domino Simulation
-----------------

.. image:: https://raw.githubusercontent.com/fury-gl/fury-communication-assets/main/physics_domino.gif
    :align: center

The code for the above simulation can be found `here <https://github.com/fury-gl/fury/blob/master/docs/examples/physics_using_pybullet/viz_domino.py>`__.
.. _community:

=========
Community
=========

Join Us!
--------

.. raw:: html

    <ul style="list-style-type:none;">
        <li style="display: block"><a href='https://discord.gg/6btFPPj'><i class="fab fa-discord fa-fw"></i> Discord</a></li>
        <li style="display: block"><a href='https://mail.python.org/mailman3/lists/fury.python.org'><i class="fa fa-envelope fa-fw"></i> Mailing list</a></li>
        <li style="display: block"><a href='https://github.com/fury-gl/fury'><i class="fab fa-github fa-fw"></i> Github</a></li>
    <ul>

Contributors
------------

.. raw:: html

    <div id="github_visualization_main_container">
        <div class="github_visualization_visualization_container">
            <div class="github_visualization_basic_stats_container">
                <div class="github_visualization_basic_stats" id="github_visualization_repo_stars">
                    <span class="stat-value banner-start-link">{{ basic_stats["stargazers_count"] }}</span> Stars
                    <img class="basic_stat_icon" src="_static/images/stars.png">
                </div>
                <div class="github_visualization_basic_stats" id="github_visualization_repo_forks">
                    <span class="stat-value">{{ basic_stats["forks_count"] }}</span> Forks
                    <img class="basic_stat_icon" src="_static/images/forks.png">
                </div>
                <div class="github_visualization_basic_stats" id="github_visualization_repo_contributors_count">
                    <span class="stat-value">{{ contributors["total_contributors"] }}</span> Contributors
                    <img class="basic_stat_icon" src="_static/images/contributors.png">
                </div>
                <div class="github_visualization_basic_stats" id="github_visualization_repo_commits_count">
                    <span class="stat-value">{{ contributors["total_commits"] }}</span> Commits
                    <img class="basic_stat_icon" src="_static/images/commits.png">
                </div>
            </div>
        <div id="github_visualization_contributors_wrapper">
        {% for contributor in contributors["contributors"] %}
        <a href="{{ contributor.html_url }}" target="_blank">
        <div class="github_visualization_contributor_info">
            <img class="github_visualization_contributor_img" src="{{ contributor.avatar_url }}">
            {% if contributor.fullname %}
            <span class="github_visualization_contributor_name">{{ contributor.fullname }}</span>
            {% else %}
            <span class="github_visualization_contributor_name">{{ contributor.username }}</span>
            {% endif %}
            <span class="github_visualization_contributor_commits">Commits: {{ contributor.nb_commits }}</span>
            <span class="github_visualization_contributor_additions"> ++{{ contributor.total_additions }}</span>
            <span class="github_visualization_contributor_deletions"> --{{contributor.total_deletions }}</span>
        </div>
        </a>
        {% endfor %}
            </div>
        </div>
    </div>============
Installation
============

FURY supports Python 3.5+. You can currently still use Python 2.7 but it will soon stop being supported as the Python 2.7 end of life is on December 31st 2019.

Dependencies
------------

The mandatory dependencies are:

- numpy >= 1.7.1
- vtk >= 8.1.0
- scipy >= 0.9

The optional dependencies are:

- matplotlib >= 2.0.0
- dipy >= 0.16.0


Installation with PyPi
----------------------

In a terminal, issue the following command

.. code-block:: shell

    pip install fury

Installation with Conda
-----------------------

Our conda package is on the Conda-Forge channel. You will need to run the following command

.. code-block:: shell

    conda install -c conda-forge fury

Installation via Source
-----------------------

**Step 1.** Get the latest source by cloning this repo

.. code-block:: shell

    git clone https://github.com/fury-gl/fury.git

**Step 2.** Install requirements

.. code-block:: shell

    pip install -r requirements/default.txt

**Step 3.** Install fury via

.. code-block:: shell

    pip install .

or

.. code-block:: shell

    pip install -e .

**Step 4:** Enjoy!

Test the Installation
---------------------

You can check your installation via this command

.. code-block:: shell

    python -c "from fury import get_info; print(get_info())"

This command will give you important information about FURY's installation. The next step will be to run a :doc:`tutorial <auto_tutorials/index>`.

Running the Tests
-----------------

Let's install all required packages for the running the test

.. code-block:: shell

    pip install -r requirements/default.txt
    pip install -r requirements/test.txt

There are two ways to run FURY tests:

- From the command line. You need to be on the FURY package folder

.. code-block:: shell

    pytest -svv fury

- To run a specific test file

.. code-block:: shell

    pytest -svv fury/tests/test_actor.py

- To run a specific test directory

.. code-block:: shell

    pytest -svv fury/tests

- To run a specific test function

.. code-block:: shell

    pytest -svv -k "test_my_function_name"

Running the Tests Offscreen
---------------------------

FURY is based on VTK which uses OpenGL for all its rendering. For a headless rendering, we recommend to install and use Xvfb software on linux or OSX.
Since Xvfb will require an X server (we also recommend to install XQuartz package on OSX). After Xvfb is installed you have 2 options to run FURY tests:

- First option

.. code-block:: shell

    export DISPLAY=:0
    Xvfb :0 -screen 1920x1080x24 > /dev/null 2>1 &
    pytest -svv fury

- Second option

.. code-block:: shell

    export DISPLAY=:0
    xvfb-run --server-args="-screen 0 1920x1080x24" pytest -svv fury


Populating our Documentation
----------------------------

Folder Structure
~~~~~~~~~~~~~~~~

Let’s start by showcasing the ``docs`` folder structure:

| fury
| ├── docs
| │   ├── build
| │   ├── make.bat
| │   ├── Makefile
| │   ├── Readme.md
| │   ├── upload_to_gh-pages.py
| │   ├── demos
| │   ├── tutorials
| │   ├── experimental
| │   └── source
| ├── requirements.txt
| ├── fury
| │   ├── actor.py
| │   ├── ...
| │
| │── ...
|
|

In our ``docs`` folder structure above:

- ``source`` is the folder that contains all ``*.rst`` files.
- ``tutorials`` is the directory where we have all python scripts that describe how to use the api.
- ``demos`` being the FURY app showcases.
- ``experimental`` directory contains experimental Python scripts. The goal is to keep a trace of expermiental work.

Building the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Step 1.** Install all required packages for the documentation generation

.. code-block:: shell

    pip install -U -r requirements/default.txt
    pip install -U -r requirements/optional.txt
    pip install -U -r requirements/docs.txt

**Step 2.** Go to the ``docs`` folder and run the following command to generate it (Linux and macOS)

.. code-block:: shell

    make -C . clean && make -C . html

To generate the documentation without running the examples

.. code-block:: shell

    make -C . clean && make -C . html-no-examples

or under Windows

.. code-block:: shell

    make clean
    make html

To generate the documentation without running the examples under Windows

.. code-block:: shell

    make clean
    make html-no-examples


**Step 3.** Congratulations! the ``build`` folder has been generated! Go to ``build/html`` and open with browser ``index.html`` to see your generated documentation.===============
Getting Started
===============

Start by importing FURY.

.. code-block:: python

    import numpy as np
    from fury import window, actor, ui, io, utils

To import a model, use :py:func:`.io.load_polydata`. Currently supported formats include OBJ, VTK, FIB, PLY, STL and XML.
Let us include the ``suzanne`` model used by Blender

.. code-block:: python

    suzanne = io.load_polydata('suzanne.obj')
    suzanne = utils.get_polymapper_from_polydata(suzanne)
    suzanne = utils.get_actor_from_polymapper(suzanne)

Set the opacity of the model::

    suzanne.GetProperty().SetOpacity(0.5)

Let's create some random variables for the cylinder parameters

.. code-block:: python

    centers = np.random.rand(2, 3)
    directions = np.random.rand(2, 3)
    heights = np.random.rand(2)
    colors = np.random.rand(2, 3)

Now, we create a cylinder::

    cylinders = actor.cylinder(centers, directions, colors, heights=heights)

Anything that has to be rendered needs to be added to the scene so let's create a :py:class:`.Scene()`::

    scene = window.Scene()

We set the window scene variables e.g. (width, height)::

    showm = window.ShowManager(scene, size=(1024,720), reset_camera=False)
    showm.initialize()

We add a text block to add some information::

    tb = ui.TextBlock2D()
    tb.message = "Hello Fury"

The function :py:meth:`.Scene.add()` is used to add the created objects to the scene to be rendered::

    scene.add(suzanne)
    scene.add(cylinders)
    scene.add(tb)

Start the rendering of the scene::

    showm.start()
test

.. toctree::
   :maxdepth: 2
   :caption: About Fury
   :hidden:

   introduction
   License <symlink/license.rst>
   blog

.. toctree::
   :maxdepth: 2
   :caption: Programming with Fury
   :hidden:

   installation
   getting_started
   auto_tutorials/index
   auto_examples/index
   fury-pybullet

.. toctree::
   :maxdepth: 2
   :caption: Developers Guide
   :hidden:

   reference/index
   symlink/contributing.rst
   release-history

.. toctree::
   :maxdepth: 1
   :caption: Join our community!
   :hidden:

   Contributors <community.rst>
   Discord Community <https://discord.gg/6btFPPj>
   Mailing list <https://mail.python.org/mailman3/lists/fury.python.org>
   Source Code <https://github.com/fury-gl/fury>
   File a bug <https://github.com/fury-gl/fury/issues/new/choose>
   Feature Request <https://github.com/fury-gl/fury/issues/new/choose>=====
About
=====

Overview
--------

Free Unified Rendering in pYthon

What is FURY?
-------------

It is a community-driven, open-source, and high-performance scientific visualization library that harnesses the graphics processing unit (GPU) for improved speed, precise interactivity, and visual clarity.

Statement of Need
-----------------

FURY was created to address the growing necessity of high-performance 3D scientific visualization in an easy-to-use API fully compatible with the Pythonic ecosystem. To achieve this FURY takes ideas from CGI (Computer-Generated Imagery) and game engines to then be deployed for usage in everyday research practices.

Mission Statement
-----------------

The purpose of FURY is to make it easier to perform advanced visualizations and animations without compromising performance. We aim to build software that is:

* clearly written
* clearly explained
* a good fit for the underlying ideas
* a natural home for collaboration

We hope that, if we fail to do this, you will let us know and we will try and make it better.

Features
--------

- **Efficient**
- **Robust**
- **Multiplatform**


License
-------

FURY is distributed under the BSD 3 License

Credits
-------

Go to :doc:`Community page <community>` to see who have been involved in the development of FURY.

Bug reports and support
-----------------------

Please report any issues via https://github.com/fury-gl/fury/issues. All types of issues are welcome including bug reports, documentation typos, feature requests and so on.
====
Blog
===================
Release History
===============

For a full list of the features implemented in the most recent release cycle, check out the release notes.

.. toctree::
   :maxdepth: 1

   release_notes/releasev0.7.1
   release_notes/releasev0.7.0
   release_notes/releasev0.6.1
   release_notes/releasev0.6.0
   release_notes/releasev0.5.1
   release_notes/releasev0.4.0
   release_notes/releasev0.3.0
   release_notes/releasev0.2.0
   release_notes/releasev0.1.4
   release_notes/releasev0.1.3
   release_notes/releasev0.1.1
   release_notes/releasev0.1.0
.. include:: ../../../CONTRIBUTING.rstFURY 0.4.0 Released
===================

.. post:: October 29 2019
   :author: skoudoro
   :tags: fury
   :category: release


The FURY project is happy to announce the release of FURY 0.4.0!
FURY is a free and open source software library for scientific visualization and 3D animations.

You can show your support by `adding a star <https://github.com/fury-gl/fury/stargazers>`_ on FURY github project.

The **major highlights** of this release are:

.. include:: ../../release_notes/releasev0.4.0.rst
    :start-after: --------------
    :end-before: Details

.. note:: The complete release notes are :ref:`available here <releasev0.4.0>`

**To upgrade or install FURY**

Run the following command in your terminal::

    pip install --upgrade fury

or::

    conda install -c conda-forge fury


**Questions or suggestions?**

For any questions go to http://fury.gl, or send an e-mail to fury@python.org
We can also join our `discord community <https://discord.gg/6btFPPj>`_

We would like to thanks to :ref:`all contributors <community>` for this release:

.. include:: ../../release_notes/releasev0.4.0.rst
    :start-after: commits.
    :end-before: We closed


On behalf of the :ref:`FURY developers <community>`,

Serge K.

Success on Brain Art Competition using FURY
===========================================

.. post:: June 19 2019
   :author: skoudoro
   :tags: shader
   :category: news


Congratulations to `David <https://github.com/dmreagan>`_ who won the `OHBM BrainArt (MELPOMENE category) <https://www.neurobureau.org/galleries/brain-art-competition-2019-2/>`_ by using DIPY and FURY!!!

As you can see below, really beautiful streamlines effect created during shader experiment/fail!

.. image:: Melpomene_Brain-Streamlines-shader-experiment_by_David-Reagan.jpg
   :alt: Brain Streamlines Shader Experiment
   :align: center
FURY 0.3.0 Released
===================

.. post:: August 2 2019
   :author: skoudoro
   :tags: fury
   :category: release


The FURY project is happy to announce the release of FURY 0.3.0!
FURY is a free and open source software library for scientific visualization and 3D animations.

You can show your support by `adding a star <https://github.com/fury-gl/fury/stargazers>`_ on FURY github project.

The **major highlights** of this release are:

.. include:: ../../release_notes/releasev0.3.0.rst
    :start-after: --------------
    :end-before: Details

.. note:: The complete release notes are :ref:`available here <releasev0.3.0>`

**To upgrade or install FURY**

Run the following command in your terminal::

    pip install --upgrade fury

or::

    conda install -c conda-forge fury


**Questions or suggestions?**

For any questions go to http://fury.gl, or send an e-mail to fury@python.org
We can also join our `discord community <https://discord.gg/6btFPPj>`_

We would like to thanks to :ref:`all contributors <community>` for this release:

.. include:: ../../release_notes/releasev0.3.0.rst
    :start-after: commits.
    :end-before: We closed


On behalf of the :ref:`FURY developers <community>`,

Serge K.
FURY 0.2.0 Released
===================

.. post:: March 8 2019
   :author: skoudoro
   :tags: fury
   :category: release


The FURY project is happy to announce the release of FURY 0.2.0!
FURY is a free and open source software library for scientific visualization and 3D animations.

You can show your support by `adding a star <https://github.com/fury-gl/fury/stargazers>`_ on FURY github project.

The **major highlights** of this release are:

.. include:: ../../release_notes/releasev0.2.0.rst
    :start-after: --------------
    :end-before: Details

.. note:: The complete release notes are :ref:`available here <releasev0.2.0>`

**To upgrade or install FURY**

Run the following command in your terminal::

    pip install --upgrade fury

or::

    conda install -c conda-forge fury


**Questions or suggestions?**

For any questions go to http://fury.gl, or send an e-mail to fury@python.org
We can also join our `discord community <https://discord.gg/6btFPPj>`_

We would like to thanks to :ref:`all contributors <community>` for this release:

.. include:: ../../release_notes/releasev0.2.0.rst
    :start-after: commits.
    :end-before: We closed


On behalf of the :ref:`FURY developers <community>`,

Serge K.
FURY 0.1.3 Released
===================

.. post:: October 31 2018
   :author: skoudoro
   :tags: fury
   :category: release


The FURY project is happy to announce the release of FURY 0.1.3!
FURY is a free and open source software library for scientific visualization and 3D animations.

You can show your support by `adding a star <https://github.com/fury-gl/fury/stargazers>`_ on FURY github project.

The **major highlights** of this release are:

.. include:: ../../release_notes/releasev0.1.3.rst
    :start-after: --------------

.. note:: The complete release notes are :ref:`available here <releasev0.1.3>`

**To upgrade or install FURY**

Run the following command in your terminal::

    pip install --upgrade fury

or::

    conda install -c conda-forge fury


**Questions or suggestions?**

For any questions go to http://fury.gl, or send an e-mail to fury@python.org
We can also join our `discord community <https://discord.gg/6btFPPj>`_


On behalf of the :ref:`FURY developers <community>`,

Serge K.
FURY 0.1.0 Released
===================

.. post:: September 21 2018
   :author: skoudoro
   :tags: fury
   :category: release


The FURY project is happy to announce the release of FURY 0.1.0!
FURY is a free and open source software library for scientific visualization and 3D animations.

You can show your support by `adding a star <https://github.com/fury-gl/fury/stargazers>`_ on FURY github project.

The **major highlights** of this release are:

.. include:: ../../release_notes/releasev0.1.0.rst
    :start-after: --------------

.. note:: The complete release notes are :ref:`available here <releasev0.1.0>`

**To upgrade or install FURY**

Run the following command in your terminal::

    pip install --upgrade fury

or::

    conda install -c conda-forge fury


**Questions or suggestions?**

For any questions go to http://fury.gl, or send an e-mail to fury@python.org
We can also join our `discord community <https://discord.gg/6btFPPj>`_


On behalf of the :ref:`FURY developers <community>`,

Serge K.
FURY 0.1.4 Released
===================

.. post:: November 26 2018
   :author: skoudoro
   :tags: fury
   :category: release


The FURY project is happy to announce the release of FURY 0.1.4!
FURY is a free and open source software library for scientific visualization and 3D animations.

You can show your support by `adding a star <https://github.com/fury-gl/fury/stargazers>`_ on FURY github project.

The **major highlights** of this release are:

.. include:: ../../release_notes/releasev0.1.4.rst
    :start-after: --------------

.. note:: The complete release notes are :ref:`available here <releasev0.1.4>`

**To upgrade or install FURY**

Run the following command in your terminal::

    pip install --upgrade fury

or::

    conda install -c conda-forge fury


**Questions or suggestions?**

For any questions go to http://fury.gl, or send an e-mail to fury@python.org
We can also join our `discord community <https://discord.gg/6btFPPj>`_


On behalf of the :ref:`FURY developers <community>`,

Serge K.
A Stadia-like system for data visualization
===========================================

.. post:: June 13 2021
   :author: Bruno Messias
   :tags: google
   :category: gsoc

Hi all! In this post I'll talk about the PR
`#437 <https://github.com/fury-gl/fury/pull/437>`__.

There are several reasons to have a streaming system for data
visualization. Because I’m doing a PhD in a developing country I always
need to think of the cheapest way to use the computational resources
available. For example, with the GPUs prices increasing, it’s necessary
to share a machine with a GPU with different users in different
locations. Therefore, to convince my Brazilian friends to use FURY I
need to code thinking inside of the (a) low-budget scenario.

To construct the streaming system for my project I’m thinking about the
following properties and behaviors:

#. I want to avoid blocking the code execution in the main thread (where
   the vtk/fury instance resides).
#. The streaming should work inside of a low bandwidth environment.
#. I need an easy way to share the rendering result. For example, using
   the free version of ngrok.

To achieve the property **1.** we need to circumvent the GIL problem.
Using the threading module alone it’s not good enough because we can’t
use the python-threading for parallel CPU computation. In addition, to
achieve a better organization it’s better to define the server system as
an uncoupled module. Therefore, I believe that multiprocessing-lib in
python will fit very well for our proposes.

For the streaming system to work smoothly in a low-bandwidth scenario we
need to choose the protocol wisely. In the recent years the WebRTC
protocol has been used in a myriad of applications like google hangouts
and Google Stadia aiming low latency behavior. Therefore, I choose the
webrtc as my first protocol to be available in the streaming system
proposal.

To achieve the third property, we must be economical in adding
requirements and dependencies.

Currently, the system has some issues, but it's already working. You can
see some tutorials about how to use this streaming system
`here <https://github.com/devmessias/fury/tree/feature_fury_stream_client/docs/tutorials/04_stream>`__.
After running one of these examples you can easily share the results and
interact with other users. For example, using the ngrok For example,
using the ngrok

::

     ./ngrok http 8000  
    

| 

How does it works?
------------------

The image below it's a simple representation of the streaming system.

|image1|

As you can see, the streaming system is made up of different processes
that share some memory blocks with each other. One of the hardest part
of this PR was to code this sharing between different objects like VTK,
numpy and the webserver. I'll discuss next some of technical issues that
I had to learn/circumvent.

Sharing data between process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We want to avoid any kind of unnecessary duplication of data or
expensive copy/write actions. We can achieve this economy of
computational resources using the multiprocessing module from python.

multiprocessing RawArray
^^^^^^^^^^^^^^^^^^^^^^^^

| The
  `RawArray <https://docs.python.org/3/library/multiprocessing.html#multiprocessing.sharedctypes.RawArray>`__
  from multiprocessing allows to share resources between different
  processes. However, there are some tricks to get a better performance
  when we are dealing with RawArray's. For example, `take a look at my
  PR in a older
  stage. <https://github.com/devmessias/fury/tree/6ae82fd239dbde6a577f9cccaa001275dcb58229>`__
  In this older stage my streaming system was working well. However, one
  of my mentors (Filipi Nascimento) saw a huge latency for
  high-resolutions examples. My first thought was that latency was
  caused by the GPU-CPU copy from the opengl context. However, I
  discovered that I've been using RawArray's wrong in my entire life!
| See for example this line of code
  `fury/stream/client.py#L101 <https://github.com/devmessias/fury/blob/6ae82fd239dbde6a577f9cccaa001275dcb58229/fury/stream/client.py#L101>`__
  The code below shows how I've been updating the raw arrays

::

   raw_arr_buffer[:] = new_data

This works fine for small and medium sized arrays, but for large ones it
takes a large amount of time, more than GPU-CPU copy. The explanation
for this bad performance is available here : `Demystifying sharedctypes
performance. <https://stackoverflow.com/questions/33853543/demystifying-sharedctypes-performance>`__
The solution which gives a stupendous performance improvement is quite
simple. RawArrays implements the buffer protocol. Therefore, we just
need to use the memoryview:

::

   memview(arr_buffer)[:] = new_data

The memview is really good, but there it's a litle issue when we are
dealing with uint8 RawArrays. The following code will cause an exception:

::

   memview(arr_buffer_uint8)[:] = new_data_uint8

There is a solution for uint8 rawarrays using just memview and cast
methods. However, numpy comes to rescue and offers a simple and a 
generic solution. You just need to convert the rawarray to a np
representation in the following way:

::

   arr_uint8_repr = np.ctypeslib.as_array(arr_buffer_uint8)
   arr_uint8_repr[:] = new_data_uint8

You can navigate to my repository in this specific `commit
position <https://github.com/devmessias/fury/commit/b1b0caf30db762cc018fc99dd4e77ba0390b2f9e>`__
and test the streaming examples to see how this little modification
improves the performance.

Multiprocessing inside of different Operating Systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Serge Koudoro, who is one of my mentors, has pointed out an issue of the
streaming system running in MacOs. I don't know many things about MacOs,
and as pointed out by Filipi the way that MacOs deals with
multiprocessing is very different than the Linux approach. Although we
solved the issue discovered by Serge, I need to be more careful to
assume that different operating systems will behave in the same way. If
you want to know more,I recommend that you read this post `Python:
Forking vs
Spawm <https://britishgeologicalsurvey.github.io/science/python-forking-vs-spawn/>`__.
And it's also important to read the official documentation from python.
It can save you a lot of time. Take a look what the
official python documentation says about the multiprocessing method

|image2| Source:\ https://docs.python.org/3/library/multiprocessing.html

.. |image1| image:: https://user-images.githubusercontent.com/6979335/121934889-33ff1480-cd1e-11eb-89a4-562fbb953ba4.png
.. |image2| image:: https://user-images.githubusercontent.com/6979335/121958121-b0ebb780-cd39-11eb-862a-37244f7f635b.png
Third week of coding!
=====================

.. post:: June 28 2021
   :author: Sajag Swami
   :tags: google
   :category: gsoc

Welcome to the fourth weekly check-in. I'll be sharing my progress for the third week of coding.

What did you do this week?
--------------------------

I made a document with code snippets and visuals to show how one can use
some vtk classes in python for molecular visualization. Classes of
interest:

-  vtkMolecule (store atomic information about the molecule).
-  vtkSimpleBondPerceiver (calculate bonding info for a vtkMolecule).
-  vtkMoleculeMapper (mapper to draw vtkMolecule object).
-  vtkPeriodicTable (stores chemical data sourced from the Blue Obelisk
   Data).

Link to the document: `Molecular_viz_vtk`_. In addition to the
document, I read some research papers recommended by my mentors to
understand some other (and potentially better) methods of ribbon
visualization. Tried to implement vtkProteinRibbonFilter usage without
using vtkPDBReader but was unsuccessful in this endeavor.

What is coming up next week?
----------------------------

Three goals for next week:

#. Implement vtkProteinRibbonFilter usage without using vtkPDBReader.
#. Make a class for vtkMolecule which can store molecular data and pass
   it on to different function for rendering purposes.
#. Read papers on surface model.

Did you get stuck anywhere?
---------------------------

Implementing vtkProteinRibbonFilter usage via vtkPolyData without using
vtkPDBReader has confounded me for some time now.

.. _Molecular_viz_vtk: https://docs.google.com/document/d/1LC2MgT9mUQK0Yo9hsI4lWqaTXHWAkSNxyBKWGAqHqe8/edit

``Au Revoir!``
Week#10: Accordion UI, Support for sprite sheet animations
======================

.. post:: August 09 2021
   :author: Antriksh Misri
   :tags: google
   :category: gsoc

What did I do this week?
------------------------
Below are the tasks that I worked on:

* `Added Accordion2D to UI sub-module <https://github.com/fury-gl/fury/pull/487>`_ : This PR adds the Accordion UI to the UI sub-module. This UI inherits from the Tree2D UI and can only be merged once the Tree2D UI is in. Here's a screenshot for reference:

    .. image:: https://i.imgur.com/klI4Tb5.png
        :width: 200
        :height: 200

* `Adding X, Y, Z Layouts <https://github.com/fury-gl/fury/pull/486>`_ :  It was pointed out in last week's meeting that in 3D space horizontal/vertical means nothing. Instead X, Y, Z are used, so, these three layouts were added on top of horizontal/vertical layouts. They also have functionality of changing the direction i.e. reverse the stacking order.
* `Added support of sprite sheet animation in Card2D <https://github.com/fury-gl/fury/pull/398>`_ : The image in Card2D was static in nature and wasn't very interesting. So, to make things a bit interesting support for animated images were added. These animations are played from a sprite sheet or a texture atlas. A buffer of processed sprite chunks is maintained and with the help of a timer callback the image in the card is updated after a certain delay which is dependent of the frame rate. Below is the demonstration:

    .. image:: https://i.imgur.com/DliSpf0.gif
        :width: 200
        :height: 200

* **Researching more about Freetype/Freetype-GL**: Apart from coding stuff, i did some more research on custom font using freetype and freetype-gl. I found some examples that used the python bindings of the c++ library and displayed custom fonts that were transformable i.e. can be rotated by some angle. Hopefully I can create a working example by this weeks meeting.

Did I get stuck anywhere?
-------------------------
No, I did not get stuck anywhere.

What is coming up next week?
----------------------------
Next week I will finish up my remaining work. Which includes addressing all PR reviews and adding some more features.

**See you guys next week!**Second week of coding!
======================

.. post:: June 21 2021
   :author: Sajag Swami
   :tags: google
   :category: gsoc

Welcome to the third weekly check-in. I'll be sharing my progress for the second week of coding.

What did you do this week?
--------------------------
I created an example to demonstrate how one can render multiple bonds (double and triple). This required me to write an algorithm to detect bonding.
I used `this blog <https://www.kaggle.com/aekoch95/bonds-from-structure-data>`_ as a reference and made a few tweaks of my own to detect the presence of double/triple bonds from interatomic distances.
The math involved in generating the coordinates of bonds was quite intriguing. Preview:
  
  .. figure:: https://user-images.githubusercontent.com/65067354/122672109-7d040c80-d1e7-11eb-815d-1d07fe47bbc4.png
    :width: 300
    :height: 300

    molecules rendered: Ethane, Ethene, Ethyne (from left to right)

In addition to this, I tried understanding the codebase of vtkMolecule, vtkSimpleBondPerceiver, vtkMoleculeMapper, vtkPeriodicTable and was able to render bond-stick models and stick models using it.
This will be of great help although it's rather slow in rendering large molecules (using shaders to improve its speed will be crucial if it's to be utilised).


  .. figure:: https://github.com/SunTzunami/gsoc2021_blog_data/blob/master/visuals/week2_wire_rep.png?raw=true
    :width: 300
    :height: 300

    Stick representation using vtkMoleculeMapper



  .. figure:: https://raw.githubusercontent.com/SunTzunami/gsoc2021_blog_data/master/visuals/week2_bs_rep.png
    :width: 300
    :height: 300

    Ball and Stick representation using vtkMoleculeMapper

What is coming up next week?
----------------------------
Try to see if the above models can be implemented using shaders. Try implementing the ribbon model using the vtkProteinRibbonFilter. The rest will be decided in the meeting with the mentors.

Did you get stuck anywhere?
---------------------------
Predicting bonds had been a problem since the past few weeks, it was resolved to a large extent by vtkSimpleBondPerceiver (the only limitation of vtkSimpleBondPerceiver being its inability to predict multiple bonds).

``Au Revoir!``
Eighth coding week!
=======================

.. post:: August 02 2021
   :author: Sajag Swami
   :tags: google
   :category: gsoc

Welcome to the ninth weekly check-in. I'll be sharing my progress for the eighth week of coding.

What did you do this week?
--------------------------

#. Updated `PR #452`_: Had an extra meeting with the mentors in which we fine-tuned the `molecular module and optimised the code to make it more pythonic.

#. I was able to generate vertices and triangles for Solvent Excluded Surfaces (SES) by using a bioconda package called `msms`_. It's based on `this paper`_ by Michel F. Sanner, Arthur J. Olson & Jean-Claude Spehner. The vertices and triangles were then sent to surface actor to generate a surface.

	 .. figure:: https://user-images.githubusercontent.com/65067354/128756004-553d1880-b6e1-4a43-99fa-5bd6a2ee70d4.png
	    :width: 300
	    :height: 300

	    SES surface generated via msms and surface actor

#. Added my GSoC blogs to the FURY blogs directory. (`PR #475`_)

Other goals will be decided in the meeting with mentors.

What is coming up next week?
----------------------------

#. Research about recent papers having good (fast) algorithms to create the molecular surfaces.
#. Create tutorials to explain how to use molecular module.

Did you get stuck anywhere?
---------------------------

No.

.. _PR #452: https://github.com/fury-gl/fury/pull/452
.. _msms: https://anaconda.org/bioconda/msms
.. _this paper: https://onlinelibrary.wiley.com/doi/10.1002/%28SICI%291097-0282%28199603%2938%3A3%3C305%3A%3AAID-BIP4%3E3.0.CO%3B2-Y
.. _PR #475: https://github.com/fury-gl/fury/pull/475

``Au Revoir!``
Weekly Check-In #7
===================

.. post:: July 19 2021
   :author: Bruno Messias
   :tags: google
   :category: gsoc


What did I do this week?
------------------------

-  `PR fury-gl/helios#16
   (merged): <https://github.com/fury-gl/helios/pull/16>`__ Helios IPC
   network layout support for MacOs

-  `PR fury-gl/helios#17
   (merged): <https://github.com/fury-gl/helios/pull/17>`__ Smooth
   animations for IPC network layout algorithms

   Before this commit was not possible to record the positions to have a
   smooth animations with IPCLayout approach. See the animation below

   |image1|

   After this PR now it's possible to tell Helios to store the evolution
   of the network positions using the record_positions parameter. This
   parameter should be passed on the start method. Notice in the image
   below how this gives to us a better visualization

   |image2|

-  `PR fury-gl/helios#13
   (merged) <https://github.com/fury-gl/helios/pull/13>`__ Merged the
   forceatlas2 cugraph layout algorithm

Did I get stuck anywhere?
-------------------------

I did not get stuck this week.

What is coming up next?
-----------------------

Probably, I'll work more on Helios. Specifically I want to improve the
memory management system. It seems that some shared memory resources are
not been released when using the IPCLayout approach.

.. |image1| image:: https://user-images.githubusercontent.com/6979335/126175596-e6e2b415-bd79-4d99-82e7-53e10548be8c.gif
.. |image2| image:: https://user-images.githubusercontent.com/6979335/126175583-c7d85f0a-3d0c-400e-bbdd-4cbcd2a36fed.gif
Fourth week of coding!
======================

.. post:: July 5 2021
   :author: Sajag Swami
   :tags: google
   :category: gsoc

Welcome to the fifth weekly check-in. I'll be sharing my progress for the fourth week of coding.

What did you do this week?
--------------------------

Created a `PR`_ for the molecular module. Enables the ability to create
three types of molecular representations:

#. Space-filling model aka calotte model and CPK model.
#. Stick model.
#. Ball and Stick model.

What is coming up next week?
----------------------------

Mentors suggested changes to be made to the molecular module which I
shall make. Other goals to be decided post mid-week meeting.

Did you get stuck anywhere?
---------------------------

Sending protein data to ProteinRibbonFilter via a vtkPolyData has been
unsuccessful so far.

.. _PR: https://github.com/fury-gl/fury/pull/452

``Au Revoir!``
Welcome to my GSoC Blog!
========================

.. post:: June 8 2021
   :author: Sajag Swami
   :tags: google
   :category: gsoc

Hi all! 
I'm Sajag Swami, a sophomore at Indian Institute of Technology, Roorkee. This summer, I will be working on adding a new functionality to **FURY** which shall enable users to 
visualise various types of proteins via different representations like Richardson aka Ribbon diagrams and molecular surface diagrams. 
As a part of my stretch goals, I’ll try to expand protein diagrams via other representations including:

1. Stick
2. Ball and stick
3. Wire
4. Pipes and Planks
5. Sphere

What did you do during the Community Bonding Period?
----------------------------------------------------
I had weekly meetings with my mentors and other core team members. In the first meeting I got acquainted with the team members and learnt about the organisation and its goal/vision.
In the later meetings we discussed about various representations of proteins and how to go about implementing them in FURY.
We discussed about various libraries which can be used to parse PDB and PDBx files.
I made a `document <https://docs.google.com/document/d/1mSoAWyXlLNrCa3hN-hiP35Lj7rURYMk5jFnWZbZp70s>`_ for the same to list pros and cons of using each library. 
I worked upon my `previous PR <https://github.com/fury-gl/fury/pull/404>`_ too during the community bonding period and fixed its docstring syntax.

As my college ended early courtesy covid, I had extra time during which I experimented and learnt more about PDB and PDBx files - the details they contain and how to parse them. 
A small backbone visualisation  of 1mb0 protein made on FURY by extracting coordinate data of its alpha carbons:

.. figure:: https://github.com/SunTzunami/gsoc2021_blog_data/blob/master/visuals/week1_backbone.png?raw=true
  :align: center

What is coming up next week?
----------------------------
I have two major goals for the next week:

1. Make an actor for the space filling model of the proteins and make PR for the same which will also include the unit tests and a small tutorial for the users.
2. Try to understand the documentation of vtkProteinRibbonFilter which will prove beneficial in generating Ribbon diagrams.

``Au Revoir!``
FURY 0.7.0 Released
===================

.. post:: March 13 2021
   :author: skoudoro
   :tags: fury
   :category: release


The FURY project is happy to announce the release of FURY 0.7.0!
FURY is a free and open source software library for scientific visualization and 3D animations.

You can show your support by `adding a star <https://github.com/fury-gl/fury/stargazers>`_ on FURY github project.

This Release is mainly a maintenance release. The **major highlights** of this release are:

.. include:: ../../release_notes/releasev0.7.0.rst
    :start-after: --------------
    :end-before: Details

.. note:: The complete release notes are available :ref:`here <releasev0.7.0>`

**To upgrade or install FURY**

Run the following command in your terminal::

    pip install --upgrade fury

or::

    conda install -c conda-forge fury


**Questions or suggestions?**

For any questions go to http://fury.gl, or send an e-mail to fury@python.org
We can also join our `discord community <https://discord.gg/6btFPPj>`_

We would like to thanks to :ref:`all contributors <community>` for this release:

.. include:: ../../release_notes/releasev0.7.0.rst
    :start-after: commits.
    :end-before: We closed


On behalf of the :ref:`FURY developers <community>`,

Serge K.
Weekly Check-In #1
==================

.. post:: June 08 2021
   :author: Bruno Messias
   :tags: google
   :category: gsoc

Hi everyone! My name is Bruno Messias currently I'm a Ph.D student at
USP/Brazil. In this summer I'll develop new tools and features for
FURY-GL. Specifically, I'll focus into developing a system for
collaborative visualization of large network layouts using FURY and VTK.

What did I do this week?
------------------------

In my first meeting the mentors explained the rules and the code of
conduct inside the FURY organization. We also made some modifications in
the timeline and discussed the next steps of my project. I started
coding during the community bonding period. The next paragraph shows my
contributions in the past weeks

-  `A FURY/VTK webrtc stream system proposal:`_ to the second part of my
   GSoC project I need to have a efficiently and easy to use streaming
   system to send the graph visualizations across the Internet. In
   addition, I also need this to my Ph.D. Therefore, I’ve been working a
   lot in this PR. This PR it’s also help me to achieve the first part
   of my project. Because I don’t have a computer with good specs in my
   house and I need to access a external computer to test the examples
   for large graphs.
-  Minor improvements into the `shader markers PR`_ and `fine tunning
   open-gl state PR`_.

Did I get stuck anywhere?
-------------------------

I’ve been stuck into a performance issue (copying the opengl framebuffer
to a python rawarray) which caused a lot of lag in the webrtc streamer.
Fortunately, I discovered that I’ve been using rawarrays in the wrong
way. My `commit`_ solved this performance issue.

What is coming up next?
-----------------------

In this week I'll focus on finish the #432 and #422 pull-requests.

.. _`A FURY/VTK webrtc stream system proposal:`: https://github.com/fury-gl/fury/pull/437
.. _shader markers PR: https://github.com/fury-gl/fury/pull/422
.. _fine tunning open-gl state PR: https://github.com/fury-gl/fury/pull/432/
.. _commit: https://github.com/fury-gl/fury/pull/437/commits/b1b0caf30db762cc018fc99dd4e77ba0390b2f9e%20
Sixth week of coding!
=====================

.. post:: July 19 2021
   :author: Sajag Swami
   :tags: google
   :category: gsoc

Welcome to the seventh weekly check-in. I'll be sharing my progress for the sixth week of coding.

What did you do this week?
--------------------------

#. Updated `Molecular module`_: made it more pythonic, implemented
   ribbon actor, added support to pass numpy arrays (earlier, atomic
   data could only be added by using the add_atom).
#. Created `PR #462`_ to:

   -  Update the helical motion animation to use a single line actor,
      added textblocks to display velocity of the particle.

      |image1|

   -  For brownian motion animation, I removed rotation(azimuth) and box
      actor, added textblock to display the number of particles and to
      show the simulation steps.

      |image2|

#. Updated surface animation (used gridUI, added multiple animations).

   |image3|

#. Created a `topic`_ on vtk discourse forum to query about gaps in
   bonds (tried resolving it by manipulating vtkProperties:
   BackfaceCulling, FrontfaceCulling but was unsuccessful).
#. Read about molecular surface (theory behind it).

What is coming up next week?
----------------------------

#. Update molecular module by adding tests, ribbon actor.
#. Try to implement molecular surface representation.
#. Interactivity of the molecules.

Did you get stuck anywhere?
---------------------------

I didn't get stuck anywhere this week.

.. _Molecular module: https://github.com/fury-gl/fury/pull/452
.. _PR #462: https://github.com/fury-gl/fury/pull/462
.. _topic: https://discourse.vtk.org/t/vtkmoleculemapper-gaps-in-bonds-on-zooming-in/6183

.. |image1| image:: https://user-images.githubusercontent.com/65067354/126033284-882ed6fd-fcc3-4a1c-8dfd-3220908859b1.png
   :width: 400px
   :height: 300px
.. |image2| image:: https://user-images.githubusercontent.com/65067354/126033291-da68cb0d-b856-48ad-9aa4-c46621052267.png
   :width: 400px
   :height: 400px
.. |image3| image:: https://user-images.githubusercontent.com/65067354/126061012-b183a47d-ed5e-4026-938b-4124da291524.png
   :width: 400px
   :height: 400px

``Au Revoir!``
First week of coding!
=====================

.. post:: June 14 2021
   :author: Sajag Swami
   :tags: google
   :category: gsoc

Welcome to the second weekly check-in. I'll be sharing my progress for the first week of coding.

What did you do this week?
--------------------------
I implemented the space filling model for proteins and created a PR for the same. Preview:

.. figure:: https://user-images.githubusercontent.com/65067354/121518963-b92cb580-ca0e-11eb-8232-3512edc04670.png

   (protein rendered: `3pgk <https://www.rcsb.org/structure/3pgk>`_)
   
The PR has: 

1. Actor for space_filling_model.

2. Two examples where I show how to visualize the proteins:

   a. In `example 1 <https://github.com/fury-gl/fury/pull/439/files#diff-2c9d065c4d4873b6ce534137cfd990cea495faffd249ff35cf51e36749883534>`_, I parse a PDBx file myself and extract the atomic info essential for constructing the model which is then used by the actor to visualize it.
   
   b. In `example 2 <https://github.com/fury-gl/fury/pull/439/files#diff-68e69b9f24731ed981cd91763f3dd078aa2bf9a4da638d561352a9cf37cfd29c>`_, I parse a PDB file by using `Biopython module <http://biopython.org/>`_ and extract the atomic info essential for constructing the model which is then used by the actor to visualize it.

I created a basic test for the actor which needs to be improved. I'll discuss how to improve the test with the mentors.

What is coming up next week?
----------------------------
I have two major goals for the next week:

1. Make an actor for the space filling model of the proteins and make PR for the same which will also include the unit tests and a small tutorial for the users.
2. Try to understand the documentation of vtkProteinRibbonFilter which will prove beneficial in generating Ribbon diagrams.

Did you get stuck anywhere?
---------------------------
I tried to create a class in python which inherits from a vtk class called vtkMoleculeReaderBase but was unsucessful in this endeavour. I'll try to find a workaround.

``Au Revoir!``
Google Summer of Code 2021 - Final Report - Bruno Messias  
=========================================================

.. post:: August 23 2021
   :author: Bruno Messias
   :tags: google
   :category: gsoc


Abstract
--------

We have changed some points of my project in the first meeting.
Specifically, we focused the efforts into developing a streaming system
using the WebRTC protocol that could be used in more generic scenarios
than just the network visualization. In addition to that, we have opted
to develop the network visualization for fury as a separated repository
and package available `here <https://github.com/fury-gl/helios>`__. The
name Helios was selected for this new network visualization system based
on the Fury rendering pipeline.

Proposed Objectives
-------------------

-  Create a streaming system (stadia-like) for FURY

   -  Should work in a low-bandwidth scenario
   -  Should allow user interactions and collaboration across the
      Internet using a web-browser

-  Helios Network System objectives:

   -  Implement the Force-Directed Algorithm with examples
   -  Implement the ForceAtlas2 algorithm using cugraph with examples
   -  Implement Minimum-Distortion Embeddings algorithm (PyMDE) and
      examples
   -  Non-blocking network algorithms computation avoiding the GIL using
      the Shared Memory approach
   -  Create the documentation and the actions for the CI

-  Stretch Goals:

   -  Create an actor in FURY to draw text efficiently using shaders
   -  Add support to draw millions of nodes using FURY
   -  Add support to control the opengl state on FURY

Objectives Completed
--------------------

-  .. rubric:: Create a streaming system (stadia-like) for FURY
      :name: create-a-streaming-system-stadia-like-for-fury

   To construct the streaming system for my project we have opted to
   follow three main properties and behaviors:

   1. avoid blocking the code execution in the main thread (where the
      vtk/fury instance resides)
   2. work inside of a low bandwidth environment
   3. make it easy and cheap to share the rendering result. For example,
      using the free version of ``ngrok``

   To achieve the first property we need to circumvent the GIL and allow
   python code to execute in parallel. Using the threading module alone
   is not good enough to reach real parallelism as Python calls in the
   same process can not execute concurrently. In addition to that, to
   achieve better organization it is desirable to define the server
   system as an uncoupled module from the rendering pipeline. Therefore,
   I have chosen to employ the multiprocessing approach for that. The
   second and third property can be only achieved choosing a suitable
   protocol for transfering the rendered results to the client. We have
   opted to implement two streaming protocols: the MJPEG and the WebRTC.
   The latter is more suitable for low-bandwidth scenarios [1].

   The image below shows a simple representation of the streaming
   system.

.. raw:: html

      <center>
      <img alt="..." height="400"
         src="https://user-images.githubusercontent.com/6979335/121934889-33ff1480-cd1e-11eb-89a4-562fbb953ba4.png"/>
      </center>

   The video below shows how our streaming system works smothly and can
   be easily integrated inside of a Jupyter notebook.

`Video: WebRTC Streaming +
Ngrok <https://user-images.githubusercontent.com/6979335/130284952-2ffbf117-7119-4048-b7aa-428e0162fb7a.mp4>`__

`Video: WebRTC Streaming +
Jupyter <https://user-images.githubusercontent.com/6979335/130284261-20e84622-427e-4a59-a46f-6a33f5473025.mp4>`__

*Pull Requests:* \* https://github.com/fury-gl/fury/pull/480

-  .. rubric:: 2D and 3D marker actor
      :name: d-and-3d-marker-actor

   This feature gave FURY the ability to efficiently draw millions of
   markers and impostor 3D spheres. This feature was essential for the
   development of Helios. This feature work with signed distance fields
   (SDFs) you can get more information about how SDFs works here [4] .

   The image below shows 1 million of markers rendered using an Intel
   HD graphics 3000.

.. raw:: html

      <center>
         <img src="https://user-images.githubusercontent.com/6979335/116004971-70927780-a5db-11eb-8363-8c0757574eb4.png"/>
      </center>

-  .. rubric:: Fine-Tunning the OpenGl State
      :name: fine-tunning-the-opengl-state

   Sometimes users may need to have finer control on how OpenGL will
   render the actors. This can be useful when they need to create
   specialized visualization effects or to improve the performance.

   In this PR I have worked in a feature that allows FURY to control the
   OpenGL context created by VTK

   *Pull Request:*

   -  https://github.com/fury-gl/fury/pull/432

-  .. rubric:: Helios Network Visualization Lib: Network Layout
      Algorithms
      :name: helios-network-visualization-lib-network-layout-algorithms

   **Case 1:** Suppose that you need to monitor a hashtag and build a
   social graph. You want to interact with the graph and at the same
   time get insights about the structure of the user interactions. To
   get those insights you can perform a node embedding using any kind of
   network layout algorithm, such as force-directed or minimum
   distortion embeddings.

   **Case 2:** Suppose that you are modelling a network dynamic such as
   an epidemic spreading or a Kuramoto model. In some of those network
   dynamics a node can change the state and the edges related to the
   node must be deleted. For example, in an epidemic model a node can
   represent a person who died due to a disease. Consequently, the
   layout of the network must be recomputed to give better insights.

   In the described cases, if we want a better (UX) and at the same time
   a more practical and insightful application of Helios, the employed
   layout algorithms should not block any kind of computation in the
   main thread.

   In Helios we already have a lib written in C (with a python wrapper)
   which performs the force-directed layout algorithm using separated
   threads avoiding the GIL problem and consequently avoiding blocking
   the main thread. But what about the other open-source network layout
   libs available on the internet? Unfortunately, most of those libs
   have not been implemented like Helios force-directed methods and
   consequently, if we want to update the network layout the Python
   interpreter will block the computation and user interaction in your
   network visualization.

   My solution for having PyMDE and CuGraph-ForceAtlas not blocking the
   main thread was to break the network layout method into two different
   types of processes: A and B and communicate both process using the
   Shared Memory approach. You can more information about this PR
   through my following posts [2], [3].

The image below show an example that I made and is available at
https://github.com/fury-gl/helios/blob/main/docs/examples/viz_mde.py

|image2| *Pull Requests:*

-  **MDE Layout:** https://github.com/fury-gl/helios/pull/6

-  **CuGraph ForceAtlas2** https://github.com/fury-gl/helios/pull/13

-  **Force-Directed and MDE improvements**
   https://github.com/fury-gl/helios/pull/14

-  .. rubric:: Helios Network Visualization Lib: Visual Aspects
      :name: helios-network-visualization-lib-visual-aspects

I’ve made several stuffs to give Helios a better visual aspects. One of
them was to give a smooth real-time network layout animations. Because
the layout computations happens into a different process that the
process responsible to render the network was necessary to record the
positions and communicate the state of layout between both process.

The GIF below shows how the network layout through IPC behaved before
these modification

.. raw:: html
   
   <center>
   <img src="https://user-images.githubusercontent.com/6979335/125310065-a3a9f480-e308-11eb-98d9-0ff5406a0e96.gif"/>
   </center>

below, you can see how after those modifications the visual aspect is
better.

.. raw:: html
   
   <center>
   <img alt="..." height="300" 
   src="https://user-images.githubusercontent.com/6979335/126175583-c7d85f0a-3d0c-400e-bbdd-4cbcd2a36fed.gif"/>
   </center>

*Pull Requests:*

-  **OpenGL SuperActors:** https://github.com/fury-gl/helios/pull/1

-  **Fixed the flickering effect**
   https://github.com/fury-gl/helios/pull/10

-  **Improvements in the network node visual aspects**
   https://github.com/fury-gl/helios/pull/15

-  **Smooth animations when using IPC layouts**
   https://github.com/fury-gl/helios/pull/17

-  .. rubric:: Helios Network Visualization Lib: CI and Documentation
      :name: helios-network-visualization-lib-ci-and-documentation

Because Helios was an project that begins in my GSoC project It was
necessary to create the documentation, hosting and more. Now we have a
online documentation available at https://heliosnetwork.io/ altough the
documentation still need some improvements.

The Helios Logo which was developed by
Filipi Nascimento.

.. raw:: html
   <img alt="Helios Network Logo" height="100" src="https://fury-gl.github.io/helios-website/_images/logo.png"/>

*Pull Requests:*

-  **CI and pytests:** https://github.com/fury-gl/helios/pull/5,
   https://github.com/fury-gl/helios/pull/20

-  **Helios Logo, Sphinx Gallery and API documentation**
   https://github.com/fury-gl/helios/pull/18

-  **Documentation improvements:**
   https://github.com/fury-gl/helios/pull/8

-  .. rubric:: Objectives in Progress
      :name: objectives-in-progress

-  .. rubric:: Draw texts on FURY and Helios
      :name: draw-texts-on-fury-and-helios

   This two PRs allows FURY and Helios to draw millions of characters in
   VTK windows instance with low computational resources consumptions. I
   still working on that, finishing the SDF font rendering which the
   theory behinds was developed here [5].

   *Pull Requests:*

   -  https://github.com/fury-gl/helios/pull/24

   -  https://github.com/fury-gl/fury/pull/489

      .. raw:: html

         <center>
         <img alt="..." height="400" src="https://user-images.githubusercontent.com/6979335/129643743-6cb12c06-3415-4a02-ba43-ccc97003b02d.png"/>
         </center>

-  .. rubric:: GSoC weekly Blogs
      :name: gsoc-weekly-blogs

   Weekly blogs were added to the FURY Website.

   *Pull Requests:*

   -  **First Evaluation:** https://github.com/fury-gl/fury/pull/476
   -  **Second Evaluation:** TBD

Timeline
--------

+----------+-----------------------------+-----------------------------+
| Date     | Description                 | Blog Link                   |
+==========+=============================+=============================+
| Week     | Welcome to my weekly Blogs! | `Weekly Check-in            |
| 1(08-    |                             | #1 <https://blogs.python-   |
| 06-2021) |                             | gsoc.org/en/demvessiass-blo |
|          |                             | g/weekly-check-in-1-21/>`__ |
+----------+-----------------------------+-----------------------------+
| Week     | Post #1: A Stadia-like      | `Weekly Check-in            |
| 2(14-    | system for data             | #                           |
| 06-2021) | visualization               | 2 <https://blogs.python-gso |
|          |                             | c.org/en/demvessiass-blog/p |
|          |                             | ost-1-a-stadia-like-system- |
|          |                             | for-data-visualization/>`__ |
+----------+-----------------------------+-----------------------------+
| Week     | 2d and 3d fake impostors    | `Weekly Check-in            |
| 3(21-    | marker; fine-tunning        | #3 <https://blogs.python-   |
| 06-2021) | open-gl state; Shared       | gsoc.org/en/demvessiass-blo |
|          | Memory support for the      | g/weekly-check-in-3-15/>`__ |
|          | streaming system;           |                             |
|          | first-version of helios:    |                             |
|          | the network visualization   |                             |
|          | lib for helios              |                             |
+----------+-----------------------------+-----------------------------+
| Week     | Post #2: SOLID, monkey      | `Weekly Check-in            |
| 4(28-    | patching a python issue and | #4                          |
| 06-2020) | network layouts through     |  <https://blogs.python-gsoc |
|          | WebRTC                      | .org/en/demvessiass-blog/po |
|          |                             | st-2-solid-monkey-patching- |
|          |                             | a-python-issue-and-network- |
|          |                             | layouts-through-webrtc/>`__ |
+----------+-----------------------------+-----------------------------+
| Week     | Code refactoring; 2d        | `Weekly Check-in            |
| 5(05-    | network layouts for Helios; | #5 <https://blogs.python-   |
| 07-2021) | Implemented the Minimum     | gsoc.org/en/demvessiass-blo |
|          | distortion embedding        | g/weekly-check-in-5-14/>`__ |
|          | algorithm using the IPC     |                             |
|          | approach                    |                             |
+----------+-----------------------------+-----------------------------+
| Week     | Post #3: Network layout     | `Weekly Check-in            |
| 6(12-    | algorithms using IPC        | #6 <https://blogs.py        |
| 07-2020) |                             | thon-gsoc.org/en/demvessias |
|          |                             | s-blog/post-3-network-layou |
|          |                             | t-algorithms-using-ipc/>`__ |
+----------+-----------------------------+-----------------------------+
| Week     | Helios IPC network layout   | `eekly Check-in             |
| 7(19-    | algorithms support for      | #7 <https://blogs.python-   |
| 07-2020) | MacOs; Smooth animations    | gsoc.org/en/demvessiass-blo |
|          | for IPC layouts;            | g/weekly-check-in-7-14/>`__ |
|          | ForceAtlas2 network layout  |                             |
|          | using cugraph/cuda          |                             |
+----------+-----------------------------+-----------------------------+
| Week     | Helios CI, Helios           | `Weekly Check-in            |
| 8(26-    | documentation               | #8 <https://blogs.python    |
| 07-2020) |                             | -gsoc.org/en/demvessiass-bl |
|          |                             | og/weekly-check-in-8-9/>`__ |
+----------+-----------------------------+-----------------------------+
| Week     | Helios documentation;       | `Weekly Check-in            |
| 9(02-    | improved the examples and   | #9 <https://blogs.python-   |
| 08-2020) | documentation of the WebRTC | gsoc.org/en/demvessiass-blo |
|          | streaming system and made   | g/weekly-check-in-9-16/>`__ |
|          | some improvements in the    |                             |
|          | compatibility removing some |                             |
|          | dependencies                |                             |
+----------+-----------------------------+-----------------------------+
| Week     | Helios documentation        | `Weekly Check-in            |
| 10(09-   | improvements; found and     | #10 <https://blogs.python-g |
| 08-2020) | fixed a bug in fury w.r.t.  | soc.org/en/demvessiass-blog |
|          | the time management system; | /weekly-check-in-10-12/>`__ |
|          | improved the memory         |                             |
|          | management system for the   |                             |
|          | network layout algorithms   |                             |
|          | using IPC                   |                             |
+----------+-----------------------------+-----------------------------+
| Week     | Created a PR that allows    | `Weekly Check-in            |
| 11(16-   | FURY to draw hundred of     | #11 <https://blogs.python-g |
| 08-2020) | thousands of characters     | soc.org/en/demvessiass-blog |
|          | without any expensive GPU;  | /weekly-check-in-11-13/>`__ |
|          | fixed the flickering effect |                             |
|          | on the streaming system;    |                             |
|          | helios node labels feature; |                             |
|          | finalizing remaining PRs    |                             |
+----------+-----------------------------+-----------------------------+

Detailed weekly tasks, progress and work done can be found
`here <https://blogs.python-gsoc.org/en/demvessiass-blog/>`__.

References
~~~~~~~~~~

[1] ( Python GSoC - Post #1 - A Stadia-like system for data
visualization - demvessias s Blog, n.d.;
https://blogs.python-gsoc.org/en/demvessiass-blog/post-1-a-stadia-like-system-for-data-visualization/

[2] Python GSoC - Post #2: SOLID, monkey patching a python issue and
network layouts through WebRTC - demvessias s Blog, n.d.;
https://blogs.python-gsoc.org/en/demvessiass-blog/post-2-solid-monkey-patching-a-python-issue-and-network-layouts-through-webrtc/

[3] Python GSoC - Post #3: Network layout algorithms using IPC -
demvessias s Blog,
n.d.)https://blogs.python-gsoc.org/en/demvessiass-blog/post-3-network-layout-algorithms-using-ipc/

[4] Rougier, N.P., 2018. An open access book on Python, OpenGL and
Scientific Visualization [WWW Document]. An open access book on Python,
OpenGL and Scientific Visualization. URL
https://github.com/rougier/python-opengl (accessed 8.21.21).

[5] Green, C., 2007. Improved alpha-tested magnification for vector
textures and special effects, in: ACM SIGGRAPH 2007 Courses on -
SIGGRAPH ’07. Presented at the ACM SIGGRAPH 2007 courses, ACM Press, San
Diego, California, p. 9. https://doi.org/10.1145/1281500.1281665

.. |image2| image:: https://user-images.githubusercontent.com/6979335/125310065-a3a9f480-e308-11eb-98d9-0ff5406a0e96.gif
Week #10: SDF Fonts 
===================

.. post:: August 09 2021
   :author: Bruno Messias
   :tags: google
   :category: gsoc


What did I do this week?
------------------------

FURY/Helios
^^^^^^^^^^^

-  `PR fury-gl/helios#22
   : <https://github.com/fury-gl/helios/pull/22>`__ Helios Documentation
   Improvements.
-  `PR fury-gl/helios#23: <https://github.com/fury-gl/helios/pull/23>`__
   A PR that makes helios IPCLayout system compatible with Windows.

FURY
^^^^

-  `PR fury-gl/fury#484: I've found and fixed a bug in FURY time
   managment system <https://github.com/fury-gl/fury/pull/484>`__
-  `PR fury-gl/fury#437: <https://github.com/fury-gl/fury/pull/437>`__

   -  Fixed the tests on Windows
   -  Improve the streaming memory managment system for IPC
      communication

-  I've developing a feature that will allows FURY to draw hundreds
   thousands of labels using texture maps and signed distance functions.
   Until now I've a sketch that at least is able to draw the labels
   using the markers billboards and bitmap fonts |image1|
-  `PR fury-gl/fury#432: <https://github.com/fury-gl/fury/pull/432>`__
   minor improvements
-  `PR #474 <https://github.com/fury-gl/fury/pull/474>`__ Helped to
   review this PR

Did I get stuck anywhere?
-------------------------

I did not get stuck this week.

What is coming up next?
-----------------------

I’ll discuss that with my mentors tomorrow.

.. |image1| image:: https://user-images.githubusercontent.com/6979335/128761833-53f53e2c-5bc0-4ff3-93c4-0ad01dc7d8eb.png
Week #9: More Layouts!
======================

.. post:: August 02 2021
   :author: Antriksh Misri
   :tags: google
   :category: gsoc

What did I do this week?
------------------------
Below are the tasks that I worked on:

* `Added Horizontal/Vertical Layout to the layout module <https://github.com/fury-gl/fury/pull/480>`_ : These PRs add support for Horizontal/Vertical layouts. These layouts allow the actors to be placed in a horizontal/vertical stack.

    .. image:: https://user-images.githubusercontent.com/54466356/127688192-8412b604-d6c7-4da9-87c4-dfae044a136e.png
        :width: 200
        :height: 200
    
    .. image:: https://user-images.githubusercontent.com/54466356/127620054-7e14f86e-1579-46ac-b4a6-ac98c109094a.png
        :width: 200
        :height: 200

* `Finalizing Card2D UI element <https://github.com/fury-gl/fury/pull/398>`_ : As panel border, wrap overflow PRs were merged this week I updated the Card2D UI to take advantage of these features.
* `Added GSoC blog posts <https://github.com/fury-gl/fury/pull/477>`_ : Added GSoC blog posts in .rst format for the FURY's blog website. Also reviewed the blog posts of other members.
* `Added support for dragging by label text/icon in Tree2D UI <https://github.com/fury-gl/fury/pull/460>`_ : Added support for dragging TreeNode2D by the label text/icon. This will help making the Tree2D as well as TreeNode2D UIs more mobile.

Did I get stuck anywhere?
-------------------------
For now I am not stuck anywhere but I have yet to start my work on freetype this could pose some trouble.

What is coming up next week?
----------------------------
Next week I will finish the remaining UI elements which includes Accordion2D, SpinBox2D.

**See you guys next week!**Week #1: Welcome to my weekly Blogs!
====================================

.. post:: June 08 2021
   :author: Antriksh Misri
   :tags: google
   :category: gsoc

Hi everyone! I am **Antriksh Misri**. I am a Pre-Final year student at MIT Pune. This summer, I will be working on **Layout Management** under FURY's `UI <https://fury.gl/latest/reference/fury.ui.html>`_ module as my primary goal. This includes addition of different classes under Layout Management to provide different layouts/arrangements in which the UI elements can be arranged. As my stretch goals, I will be working on a couple of UI components for FURY’s UI module. These 2D and 3D components will be sci-fi like as seen in the movie “**Guardians of The Galaxy**”. My objective for the stretch goals would be to develop these UI components with their respective test and tutorials such that it adds on to the UI module of FURY and doesn’t hinder existing functionalities/performance.

What did I do this week?
------------------------
During the community bonding period I got to know the mentors as well as other participants. We had an introductory meeting, in which the rules and code of conduct was explained. Also, my proposal was reviewed and modified slightly. Initially, I had to develop UI elements as my primary goal and I had to work on layout management as my stretch goals but the tasks were switched. Now I have to work on Layout Management as my primary task and develop UI in the stretch goals period. I also started coding before hand to actually make use of this free period. I worked on different PR's which are described below:-

* `Added tests for Layout module <https://github.com/fury-gl/fury/pull/434>`_ : The layout module of FURY didn't had any tests implemented, so I made this PR to add tests for **Layout** & **GridLayout** class.
* `Complied available classes for Layout Management in different libraries <https://docs.google.com/document/d/1zo981_cyXZUgMDA9QdkVQKAHTuMmKaixDRudkQi4zlc/edit>`_ : In order to decide the behavior and functionality of Layout Management in FURY, I made a document that has all classes available in different libraries to manage layout of UI elements. This document also contains code snippets for these classes.
* `Resize Panel2D UI on WindowResizeEvent <https://github.com/antrikshmisri/fury/tree/panel-resize>`_ : Currently, the **Panel2D** UI is not responsive to window resizing which means its size is static. In this branch I implemented this feature.

Did I get stuck anywhere?
-------------------------
I got stuck at Panel resizing feature. I couldn't figure out how to propagate the window invoked events to a specific actor. Fortunately, the mentors helped me to solve this problem by using **partial** from **functools**.

What is coming up next?
-----------------------
The next tasks will be decided in this week's open meeting with the mentors.

**See you guys next week!**Tenth coding week!
=======================

.. post:: August 16 2021
   :author: Sajag Swami
   :tags: google
   :category: gsoc

Welcome to the eleventh weekly check-in. I'll be sharing my progress for the tenth week of coding.

What did you do this week?
--------------------------

#. Implemented `this paper`_ to generate Van der Waals surface and solvent-accessible surface (PR created: `PR #492`_). It was a good learning experience because the first time I read the paper, I didn't understand the underlying math, it all seemed alien to me. I had to read it many times, read about the algorithms used and understand the terminologies. I had a meeting with the mentors to understand a bit of the theory which proved to be quite fruitful as I understood how to go about making the space-filling model. `This`_ blog was helpful in understanding how to use vtkMarchingCubes with numpy arrays. One of the earliest SAS rendering looked like this (this implementation was not strictly according to the paper):

	.. figure:: https://user-images.githubusercontent.com/65067354/129559593-baf201bf-720c-45f7-9269-3b31954efd5e.png
	    :width: 300
	    :height: 300
	    
	    Notice that it's rather rough

   Current implementation (this implementation was according to the paper):

	.. figure:: https://user-images.githubusercontent.com/65067354/129560374-14180b22-14b2-449b-88a6-b3140226418d.png
	    :width: 300
	    :height: 300

	    grid dimenstions = 256 × 256 × 256, used smoothing algorithms recommended by vtk

I also understood how to go about rendering volumes. I think that the ability to render volumes with FURY will be a cool capability and I'll discuss my implementation and request the mentors for feedback and ideas in the weekly meeting. Example of volume rendering:

	.. figure:: https://user-images.githubusercontent.com/65067354/129562606-50a9f0cf-e16d-4501-b0fa-a0038fda406b.png
	    :width: 300
	    :height: 300

	    grid dimenstions = 256 × 256 × 256

What is coming up next week?
----------------------------

I'll try to get `PR #452`_ merged. Documentation work to be done as GSoC coding period has come to an end.

Did you get stuck anywhere?
---------------------------

No.

.. _this paper: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0008140
.. _PR #492: https://github.com/fury-gl/fury/pull/492
.. _This: ttps://pyscience.wordpress.com/2014/09/11/surface-extraction-creating-a-mesh-from-pixel-data-using-python-and-vtk/
.. _PR #452: https://github.com/fury-gl/fury/pull/452

``Au Revoir!``
Weekly Check-In #3
==================

.. post:: June 21 2021
   :author: Bruno Messias
   :tags: google
   :category: gsoc



What did you do this week?
--------------------------

-  `PR fury-gl/fury#422
   (merged): <https://github.com/fury-gl/fury/pull/422/commits/8a0012b66b95987bafdb71367a64897b25c89368>`__
   Integrated the 3d impostor spheres with the marker actor.
-  `PR fury-gl/fury#422
   (merged): <https://github.com/fury-gl/fury/pull/422>`__ Fixed some
   issues with my maker PR which now it's merged on fury.
-  `PR fury-gl/fury#432 <https://github.com/fury-gl/fury/pull/432>`__
   I've made some improvements in my PR which can be used to fine tune 
   the opengl state on VTK.
-  `PR fury-gl/fury#437 <https://github.com/fury-gl/fury/pull/437>`__
   I've made several improvements in my streamer proposal for FURY related to memory management.


-  `PR fury-gl/helios#1 <https://github.com/fury-gl/helios/pull/1>`__
   First version of async network layout using force-directed.

Did I get stuck anywhere?
-------------------------

A python-core issue
~~~~~~~~~~~~~~~~~~~

I've spent some hours trying to discover this issue. But now it's solved
through the commit
`devmessias/fury/commit/071dab85 <https://github.com/devmessias/fury/commit/071dab85a86ec4f97eba36721b247ca9233fd59e>`__

The `SharedMemory <https://docs.python.org/3/library/multiprocessing.shared_memory.html>`__
from python>=3.8 offers a new a way to share memory resources between
unrelated process. One of the advantages of using the SharedMemory
instead of the RawArray from multiprocessing is that the SharedMemory
allows to share memory blocks without those processes be related with a
fork or spawm method. The SharedMemory behavior allowed to achieve our
jupyter integration and `simplifies the use of the streaming
system <https://github.com/fury-gl/fury/pull/437/files#diff-7680a28c3a88a93b8dae7b777c5db5805e1157365805eeaf2e58fd12a00df046>`__.
However, I saw a issue in the shared memory implementation.

Let’s see the following scenario:

::

   1-Process A creates a shared memory X
   2-Process A creates a subprocess B using popen (shell=False)
   3-Process B reads X
   4-Process B closes X
   5-Process A kills B
   4-Process A closes  X
   5-Process A unlink() the shared memory resource X

The above scenario should work flawless. Calling unlink() in X is the right way as
discussed in the python official documentation. However, there is a open
issue  related the unlink method

-  `Issue:
   https://bugs.python.org/issue38119 <https://bugs.python.org/issue38119>`__
-  `PR
   python/cpython/pull/21516 <https://github.com/python/cpython/pull/21516>`__

Fortunately, I could use a
`monkey-patching <https://bugs.python.org/msg388287>`__ solution to fix
that meanwhile we wait to the python-core team to fix the
resource_tracker (38119) issue.

What is coming up next?
-----------------------

I'm planning to work in the
`fury-gl/fury#432 <https://github.com/fury-gl/fury/pull/432>`__ and
`fury-gl/helios#1 <https://github.com/fury-gl/helios/pull/1>`__.
Week #11: Removing the flickering effect
========================================

.. post:: August 16 2021
   :author: Bruno Messias
   :tags: google
   :category: gsoc

What did I do this week?
------------------------

FURY
^^^^

-  `PR fury-gl/fury#489: <https://github.com/fury-gl/fury/pull/489>`__

  This PR give to FURY three
  pre-built texture maps using different fonts. However, is quite easy
  to create new fonts to be used in a visualization.
| It’s was quite hard to develop the shader code and find the correct
  positions of the texture maps to be used in the shader. Because we
  used the freetype-py to generate the texture and packing the glyps.
  However, the lib has some examples with bugs. But fortunelly, now
  everthing is woking on FURY. I’ve also created two different examples
  to show how this PR works.

   The first example, viz_huge_amount_of_labels.py, shows that the user can 
   draw hundreds of thounsands of characters.


   |image2|

   The second example, viz_billboad_labels.py, shows the different behaviors of the label actor. In addition, presents 
   to the user how to create a new texture atlas font to be used across different visualizations.

-  `PR fury-gl/fury#437: <https://github.com/fury-gl/fury/pull/437>`__

   -  Fix: avoid multiple OpenGl context on windows using asyncio
         The streaming system must be generic, but opengl and vtk behaves in uniques ways in each Operating System. Thus, can be tricky 
         to have the same behavior acrros different OS. One hard stuff that we founded is that was not possible to use my 
         TimeIntervals objects (implemented with threading module) with vtk. The reason for this impossibility is because we can't use 
         vtk in windows in different threads. But fortunely, moving from the threading (multithreading) to the asyncio approcach (concurrency) 
         have fixed this issue and now the streaming system is ready to be used anywhere.

   -  Flickering:
  
         Finally, I could found the cause of the flickering effect on the streaming system. 
         This flickering was appearing only when the streaming was created using the Widget object. 
         The cause seems to be a bug or a strange behavior from vtk. 
         Calling   iren.MouseWheelForwardEvent() or iren.MouseWheelBackwardEvent() 
         inside of a thread without invoking the
         Start method from a vtk instance produces a memory corruption.
         Fortunately, I could fix this behavior and now the streaming system is
         working without this glitch effect.


FURY/Helios
^^^^^^^^^^^

-  `PR fury-gl/helios#24
   : <https://github.com/fury-gl/helios/pull/24>`__

This uses the
`PRfury-gl/fury#489: <https://github.com/fury-gl/fury/pull/489>`__ to
give the network label feature to helios. Is possible to draw node
labels, update the colors, change the positions at runtime. In addition,
when a network layout algorithm is running this will automatically
update the node labels positions to follow the nodes across the screen.

|image1|

-  `PR fury-gl/helios#23:
   Merged. <https://github.com/fury-gl/helios/pull/23>`__

This PR granted compatibility between IPC Layouts and Windows. Besides
that , now is quite easier to create new network layouts using inter
process communication

Did I get stuck anywhere?
-------------------------

I did not get stuck this week.

.. |image1| image:: https://user-images.githubusercontent.com/6979335/129642582-fc6785d8-0e4f-4fdd-81f4-b2552e1ff7c7.png
.. |image2| image:: https://user-images.githubusercontent.com/6979335/129643743-6cb12c06-3415-4a02-ba43-ccc97003b02d.png
Week #4: Adding Tree UI to the UI module
========================================

.. post:: June 28 2021
   :author: Antriksh Misri
   :tags: google
   :category: gsoc

What did I do this week?
------------------------
This week I had very few tasks to do, almost all of them revolved around UI elements and adding stuff to the UI module. Earlier, it was pointed out that due to some design issues, importing certain modules into others caused circular imports which led to importing the specific modules inside a class/method which is not the best approach. This will be resolved as soon as the PR that fixes this issue is reviewed/merged in the codebase. In the meantime, all the PR's related to UI will be on hold, which is why this I had few tasks this week. The tasks are described below in detail:

* `Addition of watcher class in UI <https://github.com/fury-gl/fury/pull/448>`_ :This is finally done, as described in the previous blogs this was something that was on hold for a long time. Primarily, due to other tasks I couldn't work on this but this week due to less tasks I was able to complete the watcher class and create a PR. This PR adds support for a watcher class in the UI elements. The purpose of this class is to monitor a particular attribute from the UI element after it has been added to the scene. If the attribute changes in the real time, a user defined callback is triggered and the scene is force rendered. Currently, if any attribute of the UI element changes after it is added to the scene it does not get updated accordingly. The only way to update the UI element would be to add a custom user hook that will be triggered when a particular event that can change the attribute is invoked. This is highly ambiguous as some unmonitored event can easily change many attributes of the UI element. Also it would be really hard to add user hooks for so many events. The watcher class does this automatically, it monitors the attribute for changes and if the attribute changes, a user defined callback is triggered. If this is something that is required in the UI module, then in the future a good addition would be to monitor the UI element instance as a whole instead of a single attribute .
* `Addition of Tree UI in the UI module <https://github.com/antrikshmisri/fury/blob/bb45d1c5b6fc0495dfe4d7814fab9aefbf9f7727/fury/ui.py#L5249>`_ : Another task for this week was to work on either Tree UI or the Accordion UI. I chose to work on Tree UI as it is very interesting to implement and the logic for Tree is almost similar to that of an Accordion. So far, I have managed to implement TreeNode2D. The Tree UI contains several nodes and each node can have its own sub-nodes/children. Also, each node has an expand/collapse button which can be used to chow/hide the underlying children. The Tree UI would take some sort of data structure that contains nodes/sub-nodes and convert each node to TreeNode2D and add all the processed node to the main Panel. So far this the result I have achieved: 

    .. image:: https://i.imgur.com/WIMWsrp.png
        :width: 200
        :height: 200

    .. image:: https://i.imgur.com/u33D7Qi.png
        :width: 200
        :height: 200
* `Resize Panel2D on window resizing <https://github.com/fury-gl/fury/pull/446>`_ : This PR adds support for resizing Panel2D on WindowResizeEvent. This means that the Panle2D resizes itself with respect to the changed window size. It also retains its maximum possible size and does not overflow. Also, this PR adds support for resizing the Panel2D for the bottom right corner. A placeholder button is placed at the bottom right corner of the Panel2D and when it is dragged by the placeholder the Panel2D resize accordingly. Below is an example:

    .. image:: https://i.imgur.com/87PN7TQ.gif
        :width: 200
        :height: 200
* Also, I did some testing of GridLayout when placed inside a resizable Panel2D. This will need to be worked on before advancing any further. Currently the elements present in the Panel2D do not resize properly w.r.t the changed panel size. Hopefully, this will be fixed in the future PRs.

Did I get stuck anywhere?
-------------------------
Fortunately, I did not get stuck this week.

What is coming up next?
-----------------------
The tasks for the next week will be decided in this weeks meeting.

**See you guys next week!**Google Summer of Code
=====================

.. post:: March 9 2021
   :author: skoudoro
   :tags: google
   :category: gsoc


FURY is participating in the `Google Summer of Code 2021 <https://summerofcode.withgoogle.com/>`_ under the umbrella of the `Python Software Foundation <https://python-gsoc.org/>`_.

FURY is a free and open source software library for scientific visualization and 3D animations. FURY contains many tools for visualizing a series of scientific data including graph and imaging data.

A list of project ideas and application info is on our `GitHub Wiki <https://github.com/fury-gl/fury/wiki/Google-Summer-of-Code-2021>`_.

If you are interested in talking to us about projects, applications join us to our `discord community <https://discord.gg/aXRZmmM>`_ or drop us a line on our `mailing list <https://mail.python.org/mailman3/lists/fury.python.org>`_.

Be part of our community and Enjoy your summer of code!

Serge K.Fifth week of coding!
=====================

.. post:: July 12 2021
   :author: Sajag Swami
   :tags: google
   :category: gsoc

Welcome to the sixth weekly check-in. I'll be sharing my progress for the fifth week of coding.

What did you do this week?
--------------------------

1. Generalised the vtkProteinRibbonFilter implementation.
2. Updated the molecular module based on suggestions of team members
   and mentors (`PR #452`_).
3. Updated wave function animation (`PR #362`_).

 .. figure:: https://user-images.githubusercontent.com/65067354/125155195-d4105800-e17b-11eb-9e6d-2b66ba7a8f6e.gif
    :width: 300
    :height: 300

    an animation


What is coming up next week?
----------------------------

1. Update molecular module based on team members' suggestions and add
   tests for the same.
2. Add protein ribbon implementation to the molecular module.
3. Begin working on molecular surface model.

Did you get stuck anywhere?
---------------------------

No! **I was finally able to generalise the vtkProteinRibbonFilter implementation!!** I'm
grateful to the mentors for keeping a meeting and for helping me debug
the code. I figured out most of the stuff courtesy the meeting.

.. _PR #452: https://github.com/fury-gl/fury/pull/452
.. _PR #362: https://github.com/fury-gl/fury/pull/362
   
``Au Revoir!``
Week #3: Adapting GridLayout to work with UI
============================================

.. post:: June 21 2021
   :author: Antriksh Misri
   :tags: google
   :category: gsoc

What did I do this week?
------------------------
This week my tasks revolved around layout and UI elements. The primary goal for this week was to adapt the GridLayout to work with different UI elements. Currently, GridLayout just supports vtk actors and not UI elements, my task was to modify the class to support UI elements. The other tasks for this week are described below in detail:

* `Adapt GridLayout to support UI elements <https://github.com/fury-gl/fury/pull/443>`_ : This was the main task for the week and the aim for this was to actually modify GridLayout to support UI elements. This was not possible before because GridLayout only supported vtk actors (because of certain methods only being provided by vtk actors). I modified the main class itself along with some utility functions. The problem that I faced during this was circular imports. Currently, the structure of FURY doesn't allow certain modules to be imported into other modules because of circular imports. A way to get around this was to actually import the modules inside the methods but this is not ideal always. This will be fixed in the future PRs where the UI module will be redesigned. I also added support for grid position offsetting, which basically means that the position of the UI elements that are placed in the Grid can be offset by a global offset passed in the constructor of GridLayout class. Below is an example showing the current state of GridLayout with different UI elements. I also created a brief example to demonstrate how to use GridLayout of different cellshapes with UI elements link to which is `here <https://github.com/fury-gl/fury/pull/443/files#diff-853d17c3134e7d22de88523bb787dc05d52ec798dc2111aa0419dfd5d634350a>`_.

    .. image:: https://i.imgur.com/EX2cN1i.png
        :width: 200
        :height: 200
* `Reviewed the FileDialog2D PR <https://github.com/fury-gl/fury/pull/294>`_ : This PR added support for FileDialog2D in the UI module. The PR needed to be reviewed in order to merge it as soon as other required PRs were merged. One of the mentors already reviewed the PR briefly my job was to review the PR for any remaining bugs.
* `Study #422 PR to understand contours around the drawn markers <https://github.com/fury-gl/fury/pull/422>`_ : In my previous week's tasks I created a PR to add support for borders in Panel2D. The borders were individually customizable just like in CSS which meant 4 Rectangle2D objects were needed to represent border in each direction. This is not ideal for a scenario where a lot of Panel2D are present in the scene as it can be performance taxing. A possible solution for this was to actually look how this was implemented in the #422. This PR allowed drawing millions of markers in one call that too from the GPU. Interestingly, each marker had a contour surrounding it which is exactly what we needed for Panel2D. This is something that can be considered in the future for border implementation in other complex UI elements.
* I also continued my work on the watcher class that I mentioned in the previous week's blog. The work for this is almost done and just needs some tests implemented, which should be done soon.

Did I get stuck anywhere?
-------------------------
Fortunately, I did not get stuck this week.

What is coming up next?
-----------------------
Next week I would probably continue to work on GridLayout and possibly other layouts as well, other tasks will be decided in the next meeting.

**See you guys next week!**

.. image:: https://developers.google.com/open-source/gsoc/resources/downloads/GSoC-logo-horizontal.svg
   :height: 50
   :align: center
   :target: https://summerofcode.withgoogle.com/projects/#6653942668197888

.. image:: https://www.python.org/static/community_logos/python-logo.png
   :width: 40%
   :target: https://blogs.python-gsoc.org/en/nibba2018s-blog/

.. image:: https://python-gsoc.org/logos/FURY.png
   :width: 25%
   :target: https://fury.gl/latest/community.html

Google Summer of Code Final Work Product
========================================

.. post:: August 23 2021
   :author: Antriksh Misri
   :tags: google
   :category: gsoc

-  **Name:** Antriksh Misri
-  **Organisation:** Python Software Foundation
-  **Sub-Organisation:** FURY
-  **Project:** `FURY: Create new user interface widget <https://github.com/fury-gl/fury/wiki/Google-Summer-of-Code-2021#project-3-create-new-user-interface-widget>`_

Proposed Objectives
-------------------

* Add support for Layouts in UI elements
* Add support for Horizontal Layout
* Add support for Vertical Layout
* Add support for Layout along X, Y, Z axes.
* Stretch Goals:

  * Add Tree2D UI element to the UI sub-module
  * Add Accordion2D UI element to the UI sub-module
  * Add SpinBox2D UI element to the UI sub-module

Objectives Completed
--------------------


-  **Add support for Horizontal Layout**

   Added support for Horizontal Layout in the layout module. This layout allows the user to stack actors in a horizontal fashion. Primarily, should be used for laying out UI elements as there is no meaning of horizontal/vertical in 3D space.

   *Pull Requests:*

   -  **Horizontal Layout:** https://github.com/fury-gl/fury/pull/480
   -  **Ribbon Representation demo:** https://github.com/fury-gl/fury/pull/480

- **Add support for Vertical Layout**

  Added support for Vertical Layout in the layout module. This layout allows the user to stack actors in a vertical fashion. Primarily, should be used for laying out UI elements as there is no meaning of horizontal/vertical in 3D space.

  *Pull Requests:*

  - **Vertical Layout:** https://github.com/fury-gl/fury/pull/479
  - **Vertical Layout demo:** https://github.com/fury-gl/fury/pull/479

- **Add support for Layout along X, Y, Z axes**

  Added support for Layout along x, y, z axes. Allows user to layout different actors along any given axes. Also it allows users to switch the stacking order by passing a axis+ or axis- to the constructor.

  *Pull Requests:*

  - **X, Y, Z axes Layout:** https://github.com/fury-gl/fury/pull/486
  - **X, Y, Z axes Layout demo:** https://github.com/fury-gl/fury/pull/486

- **Add Tree2D UI element to the UI sub-module**

  Added Tree2D UI element to the UI sub-module. This allows user to visualize some data in a hierarchical fashion. Each node inside the tree can have N child nodes and the depth can be infinite. Each node can be clicked to trigger a user defined callback to perform some action. Tests and two demos were added for this UI element. Below is a screenshot for reference:

  .. image:: https://camo.githubusercontent.com/dd23b7c8503e4d01c80f2d9e84ee173e06c61eeb7c348c35aeadc75f722647ca/68747470733a2f2f692e696d6775722e636f6d2f4e49334873746c2e706e67
        :width: 200
        :height: 200

  *Pull Requests:*

  - **Tree2D UI element:** https://github.com/fury-gl/fury/pull/460
  - **Tree2D UI element demo:** https://github.com/fury-gl/fury/pull/460

- **Add Accordion2D UI element to the UI sub-module**

  Added Accordion2D to the UI sub-module. This Ui element allows users to visulize data in a tree with depth of one. Each node has a title and a content panel. The children for each node can be N if and only if the children are not nodes themselves. The child UIs can be placed inside the content panel by passing some coordinates, which can be absolute or normalized w.r.t the node content panel size. Tests and two demos were added for this UI element. Below is a screenshot for reference

  .. image:: https://camo.githubusercontent.com/9395d0ea572d7f253a051823f02496450c9f79d19ff0baf32841ec648b6f2860/68747470733a2f2f692e696d6775722e636f6d2f7854754f645a742e706e67
        :width: 200
        :height: 200

  *Pull Requests:*

  - **Accordion2D UI element:** https://github.com/fury-gl/fury/pull/487
  - **Accordion2D UI element demo:** https://github.com/fury-gl/fury/pull/487

Objectives in Progress
----------------------

-  **Add support for Layout in UI elements**

   Currently all the available layouts are only available for actors i.e. of type vtkActor2D. In order to add support for the layouts in UI elements there needs to be some tweaking in the base Layout class. Currently, the PR that adds these functionalities in stalling because of some circular imports. These will hopefully be fixed soon and as soon as the circular imports are fixed, the PR will be merged.

   *Pull Requests:*

   - **Add support for Layout in UI elements:** https://github.com/fury-gl/fury/pull/443

-  **Method to process and load sprite sheets**

   This method adds support for loading and processing a sprite sheet. This will be very useful in playing animations from a n*m sprite sheet. This also has a flag to convert the processed chunks into vtkimageData which can be directly used to update the texture in some UI elements. The primary use of this method will in a tutorial for Card2D, wherein, the image panel of the card will play the animation directly from the sprite sheet.

   *Pull Requests:*

   - **Method to process and load sprite sheets:** https://github.com/fury-gl/fury/pull/491

Other Objectives
----------------

-  **Add Card2D UI element to UI sub-module**

   Added Card2D UI element to the UI sub-module. A Card2D is generally divided into two parts i.e. the image content and the text content. This version of card has an image which can be fetched from a URL and the text content which is yet again divided into two parts i.e. the title and the body. The space distribution between the image and the text content is decided by a float between 0 and 1. A value of 0 means the image takes up no space and a value of 1 means the image consumes the whole space. Below is a demonstration:

   .. image:: https://camo.githubusercontent.com/a2e461352799b6490088de15ac041162d7bf8adf9c07485ea921b525fecd0a8e/68747470733a2f2f692e696d6775722e636f6d2f446c69537066302e676966
        :width: 200
        :height: 200
 
   *Pull Requests:*

   - **Add Card2D UI element to UI sub-module:**  https://github.com/fury-gl/fury/pull/398

-  **Resize Panel2D with WindowResizeEvent or from corner placeholder**

   Currently, the size of the Panel2D is static and cannot be changed dynamically. The size is passed in during the initialization and cannot be changed easily at runtime. This PR adds support for resizing the Panel2D dynamically by adding a placeholder icon at the bottom right corner of the panel. This icon can be click and dragged on to change the size accordingly. Other than this, the panel also retains a specific size ratio when the window is resized. This means if the window is resized in any direction the panel adapts itself w.r.t the updated size. This is done by adding relevant observers for the WindowResizeEvent and binding the relevant callback to it. Below is a quick demonstration:

    .. image:: https://camo.githubusercontent.com/3b1bf6a1b6522a6079055ff196551362fcf89a41b35ac4b32315ce02333e496d/68747470733a2f2f692e696d6775722e636f6d2f3837504e3754512e676966
        :width: 200
        :height: 200

   *Pull Requests:*

   - **Resize Panel2D with WindowResizeEvent or from corner placeholder:**  https://github.com/fury-gl/fury/pull/446

-  **Added the watcher class to UI**

   This PR adds support for a watcher class in the UI elements. The purpose of this class is to monitor a particular attribute from the UI element after it has been added to the scene. If the attribute changes in the real time, a user defined callback is triggered and the scene is force rendered.

   *Pull Requests:*

   - **Added wathcer class to the UI sub-module:**  https://github.com/fury-gl/fury/pull/448

-  **Added support for borders in Panel2D**

   The Panel2D previously, didn't support any sort of effect, the main reason behind this is that, all UI elements are individual entities that are comprised of different actors. These are not the widgets provided by vtk and in order to have some effects provided by vtk shaders must be involved. This obviously makes the whole system very complicated. The border on the other hand uses 4 Rectangle2Ds to draw the 4 borders. This makes the whole process easier but makes the Panel2D very performance heavy as we are adding 5 actors to the scene. Future iterations will replace these rectangles by textures, that way we don't compromise performance and we can have different patterns in the border. Below is a demonstration:

   .. image:: https://user-images.githubusercontent.com/54466356/121709989-bd340280-caf6-11eb-9b8a-81c65260d277.png
        :width: 200
        :height: 200
 
   *Pull Requests:*

   - **Added support for borders in Panel2D:**  https://github.com/fury-gl/fury/pull/441

-  **GSoC weekly Blogs**

    Weekly blogs were added for FURY's Website.

    *Pull Requests:*

    - **First Evaluation:** https://github.com/fury-gl/fury/pull/477

    - **Second Evaluation:** https://github.com/fury-gl/fury/pull/494

Timeline
--------

+-----------------------+------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------+
| Date                  | Description                                                      | Blog Link                                                                                                                                             |
+=======================+==================================================================+=======================================================================================================================================================+
| Week 1(08-06-2021)    | Welcome to my weekly Blogs!                                      | `Weekly Check-in #1 <https://blogs.python-gsoc.org/en/antrikshmisris-blog/week-1-welcome-to-my-weekly-blogs/>`__                                      |
+-----------------------+------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------+
| Week 2(14-06-2021)    | Feature additions in UI and IO modules                           | `Weekly Check-in #2 <https://blogs.python-gsoc.org/en/antrikshmisris-blog/week-2-feature-additions-in-ui-and-io-modules/>`__                          |
+-----------------------+------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------+
| Week 3(21-06-2021)    | Adapting GridLayout to work with UI                              | `Weekly Check-in #3 <https://blogs.python-gsoc.org/en/antrikshmisris-blog/week-3-adapting-gridlayout-to-work-with-ui/>`__                             |
+-----------------------+------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------+
| Week 4(28-06-2021)    | Adding Tree UI to the UI module                                  | `Weekly Check-in #4 <https://blogs.python-gsoc.org/en/antrikshmisris-blog/week-4-adding-tree-ui-to-the-ui-module/>`__                                 |
+-----------------------+------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------+
| Week 5(05-07-2021)    | Rebasing all PR's w.r.t the UI restructuring, Tree2D, Bug Fixes  | `Weekly Check-in #5 <https://blogs.python-gsoc.org/en/antrikshmisris-blog/week-5-rebasing-all-pr-s-w-r-t-the-ui-restructuring-tree2d-bug-fixes/>`__   |
+-----------------------+------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------+
| Week 6(12-07-2021)    | Bug fixes, Working on Tree2D UI                                  | `Weekly Check-in #6 <https://blogs.python-gsoc.org/en/antrikshmisris-blog/week-6-bug-fixes-working-on-tree2d-ui/>`__                                  |
+-----------------------+------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------+
| Week 7(19-07-2021)    | Finalizing the stalling PR's, finishing up Tree2D UI.            | `Weekly Check-in #7 <https://blogs.python-gsoc.org/en/antrikshmisris-blog/week-7-finalizing-the-stalling-pr-s-finishing-up-tree2d-ui/>`__             |
+-----------------------+------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------+
| Week 8(26-07-2020)    | Code Cleanup, Finishing up open PR's, Continuing work on Tree2D. | `Weekly Check-in #8 <https://blogs.python-gsoc.org/en/antrikshmisris-blog/week-8-code-cleanup-finishing-up-open-pr-s-continuing-work-on-tree2d/>`__   |
+-----------------------+------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------+
| Week 9(02-08-2021)    | More Layouts!                                                    | `Weekly Check-in #9 <https://blogs.python-gsoc.org/en/antrikshmisris-blog/week-9-more-layouts/>`__                                                    |
+-----------------------+------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------+
| Week 10(09-08-2021)   | Accordion UI, Support for sprite sheet animations.               | `Weekly Check-in #10 <https://blogs.python-gsoc.org/en/antrikshmisris-blog/week-10-accordion-ui-support-for-sprite-sheet-animations/>`__              |
+-----------------------+------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------+
| Week 11(16-08-2021)   | More tutorials for Accordion2D, Finalizing remaining PRs.        | `Weekly Check-in #11 <https://blogs.python-gsoc.org/en/antrikshmisris-blog/week-11-2/>`__                                                             |
+-----------------------+------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------+



Detailed weekly tasks and work done can be found
`here <https://blogs.python-gsoc.org/en/antrikshmisris-blog/>`_.Week #8: Code Cleanup, Finishing up open PRs, Continuing work on Tree2D
========================================================================

.. post:: July 26 2021
   :author: Antriksh Misri
   :tags: google
   :category: gsoc

What did I do this week?
------------------------
This week I had to work on the open PRs specifically work on the bugs that were pointed out in last week's meeting. Along side the bugs I had to continue the work on Tree2D UI element. Below is the detailed description of what I worked on this week:

* `Added new tutorial, code clean-up, bug fixes <https://github.com/fury-gl/fury/pull/460>`_ : The Tree2D had some issues with its resizing of child nodes. The size for the nodes was calculated automatically based on the vertical size occupied by its children but this could be problematic when working with sliders or UI elements that take up a lot of vertical size. To avoid this the children sizes are calculated relative to each other and the vertical size is calculated such that all children fit in perfectly. Besides this, a multiselection flag has been added that decides whether multiple child nodes can be selected or not.
* `Adding tests for corner resize callback <https://github.com/fury-gl/fury/pull/446>`_ : This PR is almost done as it was decided that WindowsResizeEvent will be ignored for now. Which leaves us with corner resizing, the callback for corner resizing didn't have any tests so the recording was redone and tests for the corner resize callback was added.
* `Fixing the failing CI's for #443 <https://github.com/fury-gl/fury/pull/443>`_ : The solution that ended up working was creating custom objects for testing of is_ui method. With this update there will be no circular dependencies and no failing CI's.
* `Addressing all comments regarding #442 <https://github.com/fury-gl/fury/pull/442>`_ : In the last meeting, a bug was pointed out wherein the text wouldn't wrap as expected. The reason for this was some minor calculation mistakes. The whole wrap_overflow method was redone and now everything works as expected. Hopefully, no bugs pop up during this week's meeting.
* `Addressing some minor comments <https://github.com/fury-gl/fury/pull/441>`_ : This PR is almost done too, there were some minor changes that were required to be addressed before this could be merged. So, these comments were fixed and hopefully this will be merged soon.
* Using different fonts using FreeType python API: A major thing that FURY needs right now is using different fonts on the fly. This is more complicated than it seems, in case of browser environment this is not a problem as browsers can support and render all fonts using various techniques. In case of a desktop environment, we need to generate the bitmap for the fonts and then use them in form of textures. For now I have created a small example that generates these bitmaps from a python API called freetype-py, the fonts are fetched from google fonts and then they are displayed as textures.
* **Starting working on Vertical Layout**: As majority of PRs are almost done, I started working on Vertical Layout. This will be hihgly inspired from the Grid Layout with obvious differences. The same techniques are being used in the Tree2D so this shouldn't be difficult to implement.

Did I get stuck anywhere?
-------------------------
The failing CI's for Grid Layout Pr needed some custom objects to remove circular dependencies. I couldn't figure out where should these custom objects go but fortunaltely the mentors showed me a quick example of where it should go.

What is coming up next week?
----------------------------
Next week I will continue my work on some other UI's and the remaining Layouts.

**See you guys next week!**SOLID, monkey patching  a python issue and  network vizualization through WebRTC
================================================================================

.. post:: July 05 2021
   :author: Bruno Messias
   :tags: google
   :category: gsoc



These past two weeks I’ve spent most of my time in the `Streaming System
PR <https://github.com/fury-gl/fury/pull/437>`__ and the `Network Layout
PR <https://github.com/fury-gl/helios/pull/1/>`__ . In this post I’ll
focus on the most relevant things I’ve made for those PRs.

Streaming System
----------------

**Pull
request** : \ `fury-gl/fury/pull/437 <https://github.com/fury-gl/fury/pull/437/>`__.

Code Refactoring
~~~~~~~~~~~~~~~~

Abstract class and SOLID
^^^^^^^^^^^^^^^^^^^^^^^^

The past weeks I've spent some time refactoring the code to see what
I’ve done let’ s take a look into this
`fury/blob/b1e985.../fury/stream/client.py#L20 <https://github.com/devmessias/fury/blob/b1e985bd6a0088acb4a116684577c4733395c9b3/fury/stream/client.py#L20>`__,
the FuryStreamClient Object before the refactoring.

The code is a mess. To see why this code is not good according to SOLID
principles let’s just list all the responsibilities of FuryStreamClient:

-  Creates a RawArray or SharedMemory to store the n-buffers
-  Creates a RawArray or SharedMemory to store the information about
   each buffer
-  Cleanup the shared memory resources if the SharedMemory was used
-  Write the vtk buffer into the shared memory resource
-  Creates the vtk callbacks to update the vtk-buffer

That’s a lot and those responsibilities are not even related to each
other. How can we be more SOLID[1]? An obvious solution is to create a
specific object to deal with the shared memory resources. But it's not
good enough because we still have a poor generalization since this new
object still needs to deal with different memory management systems:
rawarray or shared memory (maybe sockets in the future). Fortunately, we
can use the python Abstract Classes[2] to organize the code.

To use the ABC from python I first listed all the behaviors that should
be mandatory in the new abstract class. If we are using SharedMemory or
RawArrays we need first to create the memory resource in a proper way.
Therefore, the GenericImageBufferManager must have a abstract method
create_mem_resource. Now take a look into the ImageBufferManager inside
of
`stream/server/server.py <https://github.com/devmessias/fury/blob/c196cf43c0135dada4e2c5d59d68bcc009542a6c/fury/stream/server/server.py#L40>`__,
sometimes it is necessary to load the memory resource in a proper way.
Because of that, the GenericImageBufferManager needs to have a
load_mem_resource abstract method. Finally, each type of
ImageBufferManager should have a different cleanup method. The code
below presents the sketch of the abstract class


.. code-block:: python
   from abc import ABC, abstractmethod

   GenericImageBufferManager(ABC):
       def __init__(
               self, max_window_size=None, num_buffers=2, use_shared_mem=False):
            ...
       @abstractmethod
       def load_mem_resource(self):
           pass
       @abstractmethod
       def create_mem_resource(self):
           pass
       @abstractmethod
       def cleanup(self):
           pass

Now we can look for those behaviors inside of FuryStreamClient.py and
ImageBufferManger.py that does not depend if we are using the
SharedMemory or RawArrays. These behaviors should be methods inside of
the new GenericImageBufferManager.



.. code-block:: python

   # code at: https://github.com/devmessias/fury/blob/440a39d427822096679ba384c7d1d9a362dab061/fury/stream/tools.py#L491

   class GenericImageBufferManager(ABC):
       def __init__(
               self, max_window_size=None, num_buffers=2, use_shared_mem=False)
           self.max_window_size = max_window_size
           self.num_buffers = num_buffers
           self.info_buffer_size = num_buffers*2 + 2
           self._use_shared_mem = use_shared_mem
            # omitted code
       @property
       def next_buffer_index(self):
           index = int((self.info_buffer_repr[1]+1) % self.num_buffers)
           return index
       @property
       def buffer_index(self):
           index = int(self.info_buffer_repr[1])
           return index
       def write_into(self, w, h, np_arr):
           buffer_size = buffer_size = int(h*w)
           next_buffer_index = self.next_buffer_index
            # omitted code

       def get_current_frame(self):
           if not self._use_shared_mem:
           # omitted code
           return self.width, self.height, self.image_buffer_repr

       def get_jpeg(self):
           width, height, image = self.get_current_frame()
           if self._use_shared_mem:
           # omitted code
           return image_encoded.tobytes()

       async def async_get_jpeg(self, ms=33):
          # omitted code
       @abstractmethod
       def load_mem_resource(self):
           pass

       @abstractmethod
       def create_mem_resource(self):
           pass

       @abstractmethod
       def cleanup(self):
           Pass

With the
`GenericImageBufferManager <https://github.com/devmessias/fury/blob/440a39d427822096679ba384c7d1d9a362dab061/fury/stream/tools.py#L491>`__
the
`RawArrayImageBufferManager <https://github.com/devmessias/fury/blob/440a39d427822096679ba384c7d1d9a362dab061/fury/stream/tools.py#L609>`__
and
`SharedMemImageBufferManager <https://github.com/devmessias/fury/blob/440a39d427822096679ba384c7d1d9a362dab061/fury/stream/tools.py#L681>`__
is now implemented with less duplication of code (DRY principle). This
makes the code more readable and easier to find bugs. In addition, later
we can implement other memory management systems in the streaming system
without modifying the behavior of FuryStreamClient or the code inside of
server.py.

I’ve also applied the same SOLID principles to improve the CircularQueue
object. Although the CircularQueue and FuryStreamInteraction were not
violating the S from SOLID, the head-tail buffer from the CircularQueue
must have a way to lock the write/read if the memory resource is busy.
Meanwhile the
`multiprocessing.Arrays <https://docs.python.org/3/library/multiprocessing.html#multiprocessing.Array>`__
already has a context which allows lock (.get_lock()) SharedMemory
dosen’t[2]. The use of abstract class allowed me to deal with those
peculiarities. `commit
358402e <https://github.com/fury-gl/fury/pull/437/commits/358402ea2f06833f66f45f3818ccc3448b2da9cd>`__

Using namedtuples to grant immutability and to avoid silent bugs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The circular queue and the user interaction are implemented in the
streaming system using numbers to identify the type of event (mouse
click, mouse weel, ...) and where to store the specific values
associated with the event , for example if the ctrl key is pressed or
not. Therefore, those numbers appear in different files and locations:
tests/test_stream.py, stream/client.py, steam/server/app_async.py. This
can be problematic because a typo can create a silent bug. One
possibility to mitigate this is to use a python dictionary to store the
constant values, for example

.. code-block:: python

   EVENT_IDS = {
        "mouse_move" : 2, "mouse_weel": 1, #...
   }

But this solution has another issue, anywhere in the code we can change
the values of EVENT_IDS and this will produce a new silent bug. To avoid
this I chose to use
`namedtuples <https://docs.python.org/3/library/collections.html#collections.namedtuple>`__
to create an immutable object which holds all the constant values
associated with the user interactions.
`stream/constants.py <https://github.com/devmessias/fury/blob/b1e985bd6a0088acb4a116684577c4733395c9b3/fury/stream/constants.py#L59>`__

The namedtuple has several advantages when compared to dictionaries for
this specific situation. In addition, it has a better performance. A
good tutorial about namedtuples it’s available here
https://realpython.com/python-namedtuple/

Testing
~~~~~~~

My mentors asked me to write tests for this PR. Therefore, this past
week I’ve implemented the most important tests for the streaming system:
`/fury/tests/test_stream.py <https://github.com/devmessias/fury/blob/440a39d427822096679ba384c7d1d9a362dab061/fury/tests/test_stream.py>`__

Most relevant bugs
~~~~~~~~~~~~~~~~~~

As I discussed in my `third
week <https://blogs.python-gsoc.org/en/demvessiass-blog/weekly-check-in-3-15/>`__
check-in there is an open issue related to SharedMemory in python.
This"bug" happens in the streaming system through the following scenario

.. code-block:: bash 

   1-Process A creates a shared memory X
   2-Process A creates a subprocess B using popen (shell=False)
   3-Process B reads X
   4-Process B closes X
   5-Process A kills B
   4-Process A closes  X
   5-Process A unlink() the shared memory resource 

In python, this scenario translates to

.. code-block:: python

   from multiprocessing import shared_memory as sh
   import time
   import subprocess
   import sys

   shm_a = sh.SharedMemory(create=True, size=10000)
   command_string = f"from multiprocessing import shared_memory as sh;import time;shm_b = sh.SharedMemory('{shm_a.name}');shm_b.close();"
   time.sleep(2)
   p = subprocess.Popen(
       [sys.executable, '-c', command_string],
       stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
   p.wait()
   print("\nSTDOUT")
   print("=======\n")
   print(p.stdout.read())
   print("\nSTDERR")
   print("=======\n")
   print(p.stderr.read())
   print("========\n")
   time.sleep(2)
   shm_a.close()
   shm_a.unlink()

Fortunately, I could use a monkey-patching[3] solution to fix that;
meanwhile we're waiting for the python-core team to fix the
resource_tracker (38119) issue [4].

Network Layout (Helios-FURY)
----------------------------

**Pull
request**\ `fury-gl/helios/pull/1 <https://github.com/fury-gl/helios/pull/1/>`__

Finally, the first version of FURY network layout is working as you can 
see in the video below.

In addition, this already can be used with the streaming system allowing
user interactions across the internet with WebRTC protocol.

One of the issues that I had to solve to achieve the result presented in
the video above was to find a way to update the positions of the vtk
objects without blocking the main thread and at the same time allowing
the vtk events calls. My solution was to define an interval timer using
the python threading module:
`/fury/stream/tools.py#L776 <https://github.com/devmessias/fury/blob/440a39d427822096679ba384c7d1d9a362dab061/fury/stream/tools.py#L776>`__,
`/fury/stream/client.py#L112 <https://github.com/devmessias/fury/blob/440a39d427822096679ba384c7d1d9a362dab061/fury/stream/client.py#L112>`__
`/fury/stream/client.py#L296 <https://github.com/devmessias/fury/blob/440a39d427822096679ba384c7d1d9a362dab061/fury/stream/client.py#L296>`__

Refs:
-----

-  [1] A. Souly,"5 Principles to write SOLID Code (examples in Python),"
   Medium, Apr. 26, 2021.
   https://towardsdatascience.com/5-principles-to-write-solid-code-examples-in-python-9062272e6bdc
   (accessed Jun. 28, 2021).
-  [2]"[Python-ideas] Re: How to prevent shared memory from being
   corrupted ?"
   https://www.mail-archive.com/python-ideas@python.org/msg22935.html
   (accessed Jun. 28, 2021).
-  [3]“Message 388287 - Python tracker."
   https://bugs.python.org/msg388287 (accessed Jun. 28, 2021).
-  [4]“bpo-38119: Fix shmem resource tracking by vinay0410 · Pull
   Request #21516 · python/cpython," GitHub.
   https://github.com/python/cpython/pull/21516 (accessed Jun. 28,
   2021).

  Network layout algorithms using IPC
===================================

.. post:: July 12 2021
   :author: Bruno Messias
   :tags: google
   :category: gsoc

Hi all. In the past weeks, I’ve been focusing on developing Helios; the
network visualization library for FURY. I improved the visual aspects of
the network rendering as well as implemented the most relevant network
layout methods.

In this post I will discuss the most challenging task that I faced to
implement those new network layout methods and how I solved it.

The problem: network layout algorithm implementations with a blocking behavior
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Case 1:** Suppose that you need to monitor a hashtag and build a
social graph. You want to interact with the graph and at the same time
get insights about the structure of the user interactions. To get those
insights you can perform a node embedding using any kind of network
layout algorithm, such as force-directed or minimum distortion
embeddings.

**Case 2:** Suppose that you are modelling a network dynamic such as an
epidemic spreading or a Kuramoto model. In some of those network
dynamics a node can change the state and the edges related to the node
must be deleted. For example, in an epidemic model a node can represent
a person who died due to a disease. Consequently, the layout of the
network must be recomputed to give better insights.

In described cases if we want a better (UX) and at the same time a more
practical and insightful application of Helios layouts algorithms
shouldn’t block any kind of computation in the main thread.

In Helios we already have a lib written in C (with a python wrapper)
which performs the force-directed layout algorithm using separated
threads avoiding the GIL problem and consequently avoiding the blocking.
But and the other open-source network layout libs available on the
internet? Unfortunately, most of those libs have not been implemented
like Helios force-directed methods and consequently, if we want to
update the network layout the python interpreter will block the
computation and user interaction in your network visualization. How to
solve this problem?

Why is using the python threading is not a good solution?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One solution to remove the blocking behavior of the network layout libs
like PyMDE is to use the threading module from python. However, remember
the GIL problem: only one thread can execute python code at once.
Therefore, this solution will be unfeasible for networks with more than
some hundreds of nodes or even less! Ok, then how to solve it well?

IPC using python
~~~~~~~~~~~~~~~~

As I said in my previous posts I’ve created a streaming system for data
visualization for FURY using webrtc. The streaming system is already
working and an important piece in this system was implemented using the
python SharedMemory from multiprocessing. We can get the same ideas from
the streaming system to remove the blocking behavior of the network
layout libs.

My solution to have PyMDE and CuGraph-ForceAtlas without blocking was to
break the network layout method into two different types of processes: A
and B. The list below describes the most important behaviors and
responsibilities for each process

**Process A:**

-  Where the visualization (NetworkDraw) will happen
-  Create the shared memory resources: edges, weights, positions, info..
-  Check if the process B has updated the shared memory resource which
   stores the positions using the timestamp stored in the info_buffer
-  Update the positions inside of NetworkDraw instance

**Process B:**

-  Read the network information stored in the shared memory resources:
   edges , weights, positions
-  Execute the network layout algorithm
-  Update the positions values inside of the shared memory resource
-  Update the timestamp inside of the shared memory resource

I used the timestamp information to avoid unnecessary updates in the
FURY/VTK window instance, which can consume a lot of computational
resources.

How have I implemented the code for A and B?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because we need to deal with a lot of different data and share them
between different processes I’ve created a set of tools to deal with
that, take a look for example in the `ShmManagerMultiArrays
Object <https://github.com/fury-gl/helios/blob/main/helios/layouts/ipc_tools.py#L111>`__
, which makes the memory management less painful.

I'm breaking the layout method into two different processes. Thus I’ve
created two abstract objects to deal with any kind of network layout
algorithm which must be performed using inter-process-communication
(IPC). Those objects are:
`NetworkLayoutIPCServerCalc <https://github.com/devmessias/helios/blob/a0a24525697ec932a398db6413899495fb5633dd/helios/layouts/base.py#L65>`__
; used by processes of type B and
`NetworkLayoutIPCRender <https://github.com/devmessias/helios/blob/a0a24525697ec932a398db6413899495fb5633dd/helios/layouts/base.py#L135>`__
; which should be used by processes of type A.

I’ll not bore you with the details of the implementation. But let’s take
a look into some important points. As I’ve said saving the timestamp
after each step of the network layout algorithm. Take a look into the
method \_check_and_sync from NetworkLayoutIPCRender
`here <https://github.com/fury-gl/helios/blob/a0a24525697ec932a398db6413899495fb5633dd/helios/layouts/base.py#L266>`__.
Notice that the update happens only if the stored timestamp has been
changed. Also, look at this line
`helios/layouts/mde.py#L180 <https://github.com/fury-gl/helios/blob/a0a24525697ec932a398db6413899495fb5633dd/helios/layouts/mde.py#L180>`__,
the IPC-PyMDE implementation This line writes a value 1 into the second
element of the info_buffer. This value is used to inform the process A
that everything worked well. I used that info for example in the tests
for the network layout method, see the link
`helios/tests/test_mde_layouts.py#L43 <https://github.com/fury-gl/helios/blob/a0a24525697ec932a398db6413899495fb5633dd/helios/tests/test_mde_layouts.py#L43>`__

Results
~~~~~~~

Until now Helios has three network layout methods implemented: Force
Directed , Minimum Distortion Embeddings and Force Atlas 2. Here
`docs/examples/viz_helios_mde.ipynb <https://github.com/fury-gl/helios/blob/a0a24525697ec932a398db6413899495fb5633dd/docs/examples/viz_helios_mde.ipynb>`__
you can get a jupyter notebook that I’ve a created showing how to use
MDE with IPC in Helios.

In the animation below we can see the result of the Helios-MDE
application into a network with a set of anchored nodes.

|image1|

Next steps
~~~~~~~~~~

I’ll probably focus on the Helios network visualization system.
Improving the documentation and testing the ForceAtlas2 in a computer
with cuda installed. See the list of opened
`issues <https://github.com/fury-gl/helios/issues>`__

Summary of most important pull-requests:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  IPC tools for network layout methods (helios issue #7)
   `fury-gl/helios/pull/6 <https://github.com/fury-gl/helios/pull/6>`__
-  New network layout methods for fury (helios issue #7)
   `fury-gl/helios/pull/9 <https://github.com/fury-gl/helios/pull/9>`__
   `fury-gl/helios/pull/14 <https://github.com/fury-gl/helios/pull/14>`__
   `fury-gl/helios/pull/13 <https://github.com/fury-gl/helios/pull/13>`__
-  Improved the visual aspects and configurations of the network
   rendering(helios issue #12)
   https://github.com/devmessias/helios/tree/fury_network_actors_improvements
-  Tests, examples and documentation for Helios (helios issues #3 and
   #4)
   `fury-gl/helios/pull/5 <https://github.com/fury-gl/helios/pull/5>`__
-  Reduced the flickering effect on the FURY/Helios streaming system
   `fury-gl/helios/pull/10 <https://github.com/fury-gl/helios/pull/10>`__
   `fury-gl/fury/pull/437/commits/a94e22dbc2854ec87b8c934f6cabdf48931dc279 <https://github.com/fury-gl/fury/pull/437/commits/a94e22dbc2854ec87b8c934f6cabdf48931dc279>`__

.. |image1| image:: https://user-images.githubusercontent.com/6979335/125310065-a3a9f480-e308-11eb-98d9-0ff5406a0e96.gif


Weekly Check-In #8
==================

.. post:: July 26 2021
   :author: Bruno Messias
   :tags: google
   :category: gsoc


What did I do this week?
------------------------

-  `PR fury-gl/helios#18 (merged):`_ Helios Documentation

   I’ve been working on Helios documentation. Now it’s available
   online at https://fury-gl.github.io/helios-website |image1|

-   `PR fury-gl/helios#17 (merged):`_ Helios CI for tests and code
   coverage

Did I get stuck anywhere?
-------------------------

I did not get stuck this week.

What is coming up next?
-----------------------

I’ll discuss that with my mentors tomorrow.

.. _`PR fury-gl/helios#18 (merged):`: https://github.com/fury-gl/helios/pull/18
.. _`PR fury-gl/helios#17 (merged):`: https://github.com/fury-gl/helios/pull/17

.. |image1| image:: https://fury-gl.github.io/helios-website/_images/logo.pngWeek #09: Sphinx custom summary
===============================

.. post:: August 02 2021
   :author: Bruno Messias
   :tags: google
   :category: gsoc


What did I do this week?
------------------------

FURY/Helios
^^^^^^^^^^^


-  `PR fury-gl/helios#22
   : <https://github.com/fury-gl/helios/pull/22>`__ Helios Documentation
   Improvements.
   I’ve spent some time studying sphinx in order to discover how I could create a 
   custom summary inside of a template module.
   
FURY
^^^^
Added my GSoC blogs to the FURY blogs as requested by my mentors.
-  `PR fury-gl/fury#437: <https://github.com/fury-gl/fury/pull/437>`__

   - Docstrings improvements
   - Covered more tests
   - Covered tests using optional dependencies.
   - Aiortc now it’s not a mandatory dependency
   - improvements in memory management

- PR #432 Fixed some typos, improved the tests and docstrings
- `PR fury-gl/fury#474: <https://github.com/fury-gl/fury/pull/474>`__
- Helped to review and made some suggestions to the PR #474 made by @mehabhalodiya.


Did I get stuck anywhere?
-------------------------

I did not get stuck this week.

What is coming up next?
-----------------------

I’ll discuss that with my mentors tomorrow.FURY 0.7.0 Released
===================

.. post:: August 03 2021
   :author: skoudoro
   :tags: fury
   :category: release


The FURY project is happy to announce the release of FURY 0.7.1!
FURY is a free and open source software library for scientific visualization and 3D animations.

You can show your support by `adding a star <https://github.com/fury-gl/fury/stargazers>`_ on FURY github project.

This Release is mainly a maintenance release. The **major highlights** of this release are:

.. include:: ../../release_notes/releasev0.7.1.rst
    :start-after: --------------
    :end-before: Details

.. note:: The complete release notes are available :ref:`here <releasev0.7.1>`

**To upgrade or install FURY**

Run the following command in your terminal::

    pip install --upgrade fury

or::

    conda install -c conda-forge fury


**Questions or suggestions?**

For any questions go to http://fury.gl, or send an e-mail to fury@python.org
We can also join our `discord community <https://discord.gg/6btFPPj>`_

We would like to thanks to :ref:`all contributors <community>` for this release:

.. include:: ../../release_notes/releasev0.7.1.rst
    :start-after: commits.
    :end-before: We closed


On behalf of the :ref:`FURY developers <community>`,

Serge K.
Week #11: Finalizing open Pull Requests
=======================================

.. post:: August 16 2021
   :author: Antriksh Misri
   :tags: google
   :category: gsoc

What did I do this week?
------------------------
Below are the tasks that I worked on:

* `Created PR for sprite sheet animation <https://github.com/fury-gl/fury/pull/491>`_ : This PR adds support for playing animations from a sprite sheet. This feature will be used in Card2D to create a tutorial in which the card will show the animation in the image box. Previously, the utility functions for this were added directly inside the tutorial but now they are refactored to go in their respective modules.
* `Finalized the x, y, z layouts <https://github.com/fury-gl/fury/pull/486>`_ : The PR that adds these layouts needed some updates for it to work as intended. These changes were added and this PR is ready to go.
* `Resolved all conflicts in the GridLayout PR <https://github.com/fury-gl/fury/pull/443>`_ : As the Horizontal and Vertical layouts were merged this week the GridLayout PR had got some conflicts. These conflicts were resolved and the PR is almost ready.
* **Continuing the work on custom font rendering** : In the last meeting, a few points were brought up. Firstly, to position each glyph to their respective location in the atlas a seperate module is used which is freetype-gl. The python bindings for this module are not available which means either we have to write the bindings ourselves or the freetype team will be emailed about this and they will add bindings for that. On the other hand, I looked how latex is rendered in matplotlib. `This <https://github.com/matplotlib/matplotlib/blob/3a4fdea8d23207d67431973fe5df1811605c4132/lib/matplotlib/text.py#L106>`_ is the Text class that is used to represent the string that is to be drawn and `This <https://github.com/matplotlib/matplotlib/blob/3a4fdea8d23207d67431973fe5df1811605c4132/lib/matplotlib/artist.py#L94>`_ is the class that it inherits from. Everything is handled internally in matplotlib, to draw the rasterized text `this <https://github.com/matplotlib/matplotlib/blob/3a4fdea8d23207d67431973fe5df1811605c4132/lib/matplotlib/text.py#L672>`_ function is used. The text can be rendered in two ways, the first one is by using the default renderer and the second way is by using PathEffectRenderer that is used to add effects like outlines, anti-aliasing etc. It is a very rigid way of rendering text and is designed to be used internally.

Did I get stuck anywhere?
-------------------------
No, I did not get stuck anywhere.

What is coming up next week?
----------------------------
Hopefully everything is resolved by the end of this week and next week I will hopefully submit my final code in a gist format.

**See you guys next week!**Week #5: Rebasing all PRs w.r.t the UI restructuring, Tree2D, Bug Fixes
========================================================================

.. post:: July 05 2021
   :author: Antriksh Misri
   :tags: google
   :category: gsoc

What did I do this week?
------------------------
The UI restructuring was finally merged this week. This means UI is now a submodule in itself and provides different UI elements based on their types which are, core, elements, containers and some helper methods/classes. So, this week my main tasks were to rebase and fix merge conflicts of my open PR's. Other than that, I had to continue my work on Tree2D UI element and finish the remaining aspects of it. Also, I had to create an example demonstrating how to use the newly added UI element. Many use cases were discussed in the open meeting like, an image gallery which displays preview image on the title and when expanded reveals the full scale image. I am thinking of adding multiple examples for the Tree2D to brainstorm on its certain features. Also, I had this annoying bug in Panel2D which didn't allow it to be resized from the bottom right corner. It was resizing from the top right corner. I had to address this bug as well. Below are the tasks in detail:

* `Rebasing all PRs w.r.t the UI restructuring <https://github.com/fury-gl/fury/pulls/antrikshmisri>`_: As discussed in the earlier blogs, due to circular imports and large size of the UI module, a bit of restructuring was required. This week the PR that converts the UI into a sub module was finally merged. This meant I had to fix all the merge conflicts and rebase all UI related PR's. So, I rebased all the UI related PR's and fixed the merge conflicts. Currently, there are still some PR's that need some attention as still some of the imports are circular in nature. This means if the issue is correct then some more restructuring is required, which will be hopefully done in the near future.
* `Continuing the work on Tree2D <https://github.com/antrikshmisri/fury/blob/86b16ba3f74c3bdcf9aab58f546b37b919254cd1/fury/ui/elements.py#L3278>`_ : This week I continued my work on Tree2D, TreeNode2D. I had to fix/add multiple features on both the classes but my priority was to fix the looks of the UI element as well as make it easier for the user to manipulate the UI element. The first thing that I fixed was node offsetting, when a node is collapsed and expanded the nodes below the current node should also offset. Previously, only the child nodes within the same parents were offset and not the nodes/parent beyond that. With some minor adjusting, now the nodes are offset recursively and all child/parent nodes below the current nodes are offset. Secondly, before only a node could be added to a node which meant it wasn't any easy way to add any other UI element to a node but with some updates/fixes any UI element can be added to a node. Also, the Tree2D lacked some properties/methods to easily manipulate it. So, i added some properties/methods that allow to easily/efficiently manipulate individual node inside the Tree2D. Below is the current state of the Tree2D. In the below tree, two panels are added to a child node after the tree has been initialized. Also, the coordinated of the child elements are totally fluid i.e they can be placed anywhere inside the content panel by passing normalized or absolute coordinates.

    .. image:: https://i.imgur.com/rQQvLqs.png
        :width: 200
        :height: 200

* Fixed Panel2D bottom corner resizing: Previously, the panel would not resize from the bottom left corner but it would resize from top right corner. I didn't understand what was going wrong and was stuck on this for a long time. But I finally figured out the problem, I was calculating the Y-offset wrong as well as the panel resized from the top side instead of bottom. With some minor tweaking the bug was gone and the panel resizes correctly from the bottom right corner.

Did I get stuck anywhere?
-------------------------
I got stuck on recording events for the updated panel UI element. The panel updates w.r.t the window size but I couldn't figure out how to record the events invoked by the window. Unfortunately, I still haven't figured out how this will be done. My guess is that I have to propagate the event first to the interactor and then to the UI element.

What is coming up next?
-----------------------
I will probably finish up the GridLayout, Tree2D UI along side some other UI's. This will be decided in the next meeting.

**See you guys next week!**Weekly Check-In #5
===================

.. post:: July 05 2021
   :author: Bruno Messias
   :tags: google
   :category: gsoc

What did you do this week?
--------------------------

`fury-gl/fury PR#437: WebRTC streaming system for FURY`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Before the `8c670c2`_ commit, for some versions of MacOs the
   streaming system was falling in a silent bug. I’ve spent a lot of
   time researching to found a cause for this. Fortunately, I could found
   the cause and the solution. This troublesome MacOs was falling in a
   silent bug because the SharedMemory Object was creating a memory
   resource with at least 4086 bytes indepedent if I've requested less
   than that. If we look into the MultiDimensionalBuffer Object
   (stream/tools.py) before the 8c670c2 commit we can see that Object
   has max_size parameter which needs to be updated if the SharedMemory
   was created with a "wrong" size.

`fury-gl/helios PR 1: Network Layout and SuperActors`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the past week I've made a lot of improvements in this PR, from
performance improvements to visual effects. Bellow are the list of the
tasks related with this PR:

-  - Code refactoring.
-  - Visual improvements: Using the UniformTools from my pull request
   `#424`_ now is possible to control all the visual characteristics at
   runtime.
-  - 2D Layout: Meanwhile 3d network representations are very usefully
   for exploring a dataset is hard to convice a group of network
   scientists to use a visualization system which dosen't allow 2d
   representations. Because of that I started to coding the 2d behavior
   in the network visualization system.
-  - Minimum Distortion Embeddings examples: I've created some examples
   which shows how integrate pymde (Python Minimum Distortion
   Embeddings) with fury/helios. The image below shows the result of
   this integration: a "perfect" graph embedding

.. image:: https://user-images.githubusercontent.com/6979335/124524052-da937e00-ddcf-11eb-83ca-9b58ca692c2e.png

What is coming up next week?
----------------------------

I'll probably focus on the `heliosPR#1`_. Specifically, writing tests
and improving the minimum distortion embedding layout.

Did you get stuck anywhere?
---------------------------

I did not get stuck this week.

.. _`fury-gl/fury PR#437: WebRTC streaming system for FURY`: https://github.com/fury-gl/fury/pull/427
.. _8c670c2: https://github.com/fury-gl/fury/pull/437/commits/8c670c284368029cdb5b54c178a792ec615e4d4d
.. _`fury-gl/helios PR 1: Network Layout and SuperActors`: https://github.com/fury-gl/helios/pull/1
.. _#424: https://github.com/fury-gl/fury/pull/424
.. _heliosPR#1: Ninth coding week!
=======================

.. post:: August 09 2021
   :author: Sajag Swami
   :tags: google
   :category: gsoc

Welcome to the tenth weekly check-in. I'll be sharing my progress for the ninth week of coding.

What did you do this week?
--------------------------

#. Updated `PR #452`_ :

	- Made ribbon representation faster.
	- Added an actor to display bounding box around the molecule.

	 .. figure:: https://user-images.githubusercontent.com/65067354/128624529-03c026be-7f80-4792-b57e-eceeb1767ec2.png
	    :width: 300
	    :height: 300

	    Bounding Box

#. Made a tutorial which showcases the abilities of molecular module (will create a PR after molecular module is merged).

#. I'm trying to implement a native implementation of molecular surfaces in FURY. Currently searching for recent research papers to find good algorithms to generate the molecular surfaces (the ones I'd collected in the research period were archaic and rather time consuming). The papers that I've read so far seem a tad bit intimidating as I've never done math related to this domain yet. Implementing them will be a good learning experience I reckon.

What is coming up next week?
----------------------------

#. Try to create a native implementation of molecular surface.
#. Small fixes to `PR #362`_, `PR #462`_.

Did you get stuck anywhere?
---------------------------

No.

.. _PR #452: https://github.com/fury-gl/fury/pull/452
.. _PR #362: https://github.com/fury-gl/fury/pull/362
.. _PR #462: https://github.com/fury-gl/fury/pull/462

``Au Revoir!``
Seventh week of coding!
=======================

.. post:: July 26 2021
   :author: Sajag Swami
   :tags: google
   :category: gsoc

Welcome to the eighth weekly check-in. I'll be sharing my progress for the seventh week of coding.

What did you do this week?
--------------------------

#. Updated `PR #452`_:

   -  Added ribbon actor to the molecular module.
   -  Added tests for all functions in the molecular module.

#. Updated `PR #462`_: Addressed the reviews of team members and
   mentors, added new features.

   |image1|

#. Updated `PR #362`_: Addressed the feedbacks of team members and
   mentors.

What is coming up next week?
----------------------------

#. Work more on molecular module, meeting with mentors and core team on
   Thursday to optimize the module and merge `PR #452`_.
#. Start working upon molecular surface model.

Did you get stuck anywhere?
---------------------------

No.

.. _PR #452: https://github.com/fury-gl/fury/pull/452
.. _PR #462: https://github.com/fury-gl/fury/pull/462
.. _PR #362: https://github.com/fury-gl/fury/pull/362

.. |image1| image:: https://user-images.githubusercontent.com/65067354/126382288-b755c01d-8010-43ab-87db-2f1a4fb5b015.png
   :width: 300px
   :height: 300px

``Au Revoir!``
Week #6: Bug fixes, Working on Tree2D UI
========================================

.. post:: July 12 2021
   :author: Antriksh Misri
   :tags: google
   :category: gsoc

What did I do this week?
------------------------
This week I had relatively few tasks and most of them were to fix some bugs/design flaws that were discussed in last week's meeting. Other than that, I had to implement a few things in the Tree2D UI element that will be discussed in detail below. I also had to update some existing PRs in order to make things work well. Below are the list of things I worked on this week:

* `Extracted Button2D class from elements to core <https://github.com/fury-gl/fury/pull/459>`_ : Button2D was placed in elements during the UI restructuring. The problem with that was that Button2D was needed in other UI elements outside UI elements present in elements in Panel2D. So, it was decided to create a rule that only the UI elements that do not depend on any other UI element are allowed to be placed in core UI elements. Button2D does not depend on any other UI element so it was extracted from elements to core.

* `Adapting GridLayout to work with UI elements <https://github.com/fury-gl/fury/pull/443>`_ : This was a PR that aimed to add support for UI elements to be placed in a grid fashion. the problem was that there still some circular imports even after UI restructuring, primarily because UI was being imported in the layout module that in turn indirectly imported some UI elements making them circularly dependent. To remove the circular imports, it was decided to determine the UI element by checking for a add_to_scene method attribute in the instance. I updated the existing PR to implement the same.

* `Continuing my work on Tree2D <https://github.com/fury-gl/fury/pull/460>`_: The Tree2D lacked some important things related to design and visual aspect of it. Before, if the children of any node exceeded its height they would just overflow. To prevent this I came up with a few solutions two of which were to either add a scrollbar on the overflowing node or to simply auto resize the parent node. Currently, there is no global API for the scrollbar and it has to be manually setup in a UI element, this will be hopefully implemented in the near future probably using layout management. Till then the auto resizing has been implemented for the nodes. In future, an option for scrollbar will be added.

Did I get stuck anywhere?
-------------------------
I am still stuck with adding tests for panel resizing PR. As it needs windows events to be recorded as well. I tried to propagate the event to the interactor first but that just lead to that particular event being registered multiple times. I will try to find a workaround for it.

What is coming up next?
-----------------------
If the Tree2D gets merged by this week then I'll probably work on other UI elements. Other tasks will be decided in the next meeting.

**See you guys next week!**Week #7: Finalizing the stalling PRs, finishing up Tree2D UI.
==============================================================

.. post:: July 19 2021
   :author: Antriksh Misri
   :tags: google
   :category: gsoc

What did I do this week?
------------------------
This week I had limited tasks to do, mostly tasks related to existing PRs. Other than some minor fixes I had to implement some more things in Tree2D which included some minor UI fixes, some changes in tutorial, adding tests. Below is the detailed description of what I worked on this week:

* `Tests, tutorial changes, UI fixes for Tree2D <https://github.com/fury-gl/fury/pull/460>`_ : The Tree2D lacked some things like proper UI resizing, relative indentation, tests for the UI class. These were added with this PR. Currently, the indentation, resizing needs some improvement, which will be fixed after feedback from this week's meeting. Also, tests for Tree2D, TreeNode2D were added as well.
* `Updating Panel2D tests, re-recording the events <https://github.com/fury-gl/fury/pull/446>`_ : This PR is almost done with just some tests blocking the PR. The tests were added this week, but tests for some callbacks that are associated with window event are still not added. This is because there is no way to count the WindowResizeEvent without actually using the API of the window provided by the OS. This can become very complicated very soon so, these tests may be added in the future.
* `Fixing the failing CI's for #443 <https://github.com/fury-gl/fury/pull/443>`_ : The CI was failing on this PR and needed some fixing which was done this week. This PR still needs some refactoring before the all CI's pass. This will hopefully be fixed before this week's meeting.
* `Addressing all comments regarding #442 <https://github.com/fury-gl/fury/pull/442>`_ : Previously, it was pointed out that the some code can be extracted into a function and can be reused in other methods. So, this week the extracted method was updated to reuse even more code and now almost no code is repeated.
* `Adding has_border flag in Panel2D <https://github.com/fury-gl/fury/pull/441>`_ : Adding a has_border flag in Panel2D: Previously, to create the borders 4 Rectangle2D's were used and they were created everytime even when border_width was set to 0. This would take a lot of wasted system resources. To fix this, a flag is added in the the constructor which is by default set to False. If false, the borders are not initialized and the resources are saved.

Did I get stuck anywhere?
-------------------------
Fortunately, this week I didn't get stuck anywhere.

**See you guys next week!**Week #2: Feature additions in UI and IO modules
===============================================

.. post:: June 13 2021
   :author: Antriksh Misri
   :tags: google
   :category: gsoc

What did I do this week?
------------------------
This week I had to work on 3 PRs as well as some documentation. I really enjoyed this week's work as the tasks were really interesting. The aim for these PRs were to actually add a couple of features in the UI as well as the IO module, which includes, adding support for border in Panel2D, adding support for network/URL images in load_image method in IO module, adding resizing Panel2D from bottom right corner, completing the document with layout solutions provided by Unity/Unreal engine. Below are the PRs that I worked on:

* `Added support for URL image in load_image <https://github.com/fury-gl/fury/pull/440>`_ : The load_image of IO module didn't support network /URL images, so I made this PR to add support for the same.
* `Added support for border in Panel2D <https://github.com/fury-gl/fury/pull/441>`_ : This PR was made in association with the Card2D PR. This PR adds support for border in Panel2D. The borders are individually customizable just like in CSS. This PR needs a little tweaking in terms of getters/setters. The same support needs to be added in Rectangle2D.
* `Complete the document with layout solutions provided by Unity/Unreal engine <https://docs.google.com/document/d/1zo981_cyXZUgMDA9QdkVQKAHTuMmKaixDRudkQi4zlc/edit?usp=sharing>`_ : Completed the document with layout solutions provided by Unity/Unreal Engine.
* Behind the scenes I also worked on a Watcher class for the UI elements. The purpose of the watcher would be to monitor the UI elements for any changes after they have been added to the scene. A PR should be up by 2-3 days.

Did I get stuck anywhere?
-------------------------
I had a minor issue with the tests for the **IO** module. When running the tests for IO module using **pytest 5.0.0** resulted in Window fatal error, this was a sideeffect of pytest 5.0.0 wherein support for **faulthandler** was added. This error was suppressed by using certain flags while running the tests.

What is coming up next?
-----------------------
Next week I would probably work on adapting the **GridLayout** with UI elements, some other tasks that will be decided in the next meeting.

**See you guys next week!**.. image:: https://developers.google.com/open-source/gsoc/resources/downloads/GSoC-logo-horizontal.svg
   :height: 50
   :align: center
   :target: https://summerofcode.withgoogle.com/projects/#6653942668197888

.. image:: https://www.python.org/static/community_logos/python-logo.png
   :width: 40%
   :target: https://blogs.python-gsoc.org/en/nibba2018s-blog/

.. image:: https://python-gsoc.org/logos/FURY.png
   :width: 25%
   :target: https://fury.gl/latest/community.html

Google Summer of Code Final Work Product
========================================

.. post:: August 23 2021
   :author: Sajag Swami
   :tags: google
   :category: gsoc

-  **Name:** Sajag Swami
-  **Organisation:** Python Software Foundation
-  **Sub-Organisation:** FURY
-  **Project:** `FURY: Ribbon and Molecular Surface Representations for 
   Proteins <https://github.com/fury-gl/fury/wiki/Google-Summer-of-Code-2021>`_

Proposed Objectives
-------------------

* Ribbon Representation
* Molecular Surface Representation
* Stretch Goals:

  * Stick Representation
  * Ball and stick Representation
  * Wire Representation
  * Pipes and Planks Representation
  * Sphere Representation

Objectives Completed
--------------------


-  **Ribbon Representation**

   Ribbon diagrams, also known as Richardson diagrams,
   are 3D schematic representations of protein structure. Ribbon diagrams are
   generated by interpolating a smooth curve through the polypeptide backbone.
   α-helices are shown as coiled ribbons. β-strands as sheets, and non-repetitive
   coils or loops as lines or thin strips. It was implemented by using
   `vtkProteinRibbonFilter`. Generating a `vtkPolyData` of appropriate format
   required by `vtkProteinRibbonFilter` was initially unclear due to lack of
   examples. I was able to understand what kind of output the filter required 
   after a meeting with mentors. Tests were added and a demo was created.

   *Pull Requests:*

   -  **Ribbon representation:** https://github.com/fury-gl/fury/pull/452
   -  **Ribbon Representation demo:** https://github.com/fury-gl/fury/pull/452


- **Ball and Stick Representation**

  The ball-and-stick model is a molecular model of a chemical substance which
  displays both the three-dimensional position of the atoms and the bonds between
  them. The atoms are typically represented by spheres, connected by tubes which
  represent the bonds. It was created by using  `vtkOpenGLMoleculeMapper`. 
  Added `vtkSimpleBondPerceiver` for detection of bonds. Tests were added and a 
  demo was created.

  *Pull Requests:*

  - **Ball and Stick Representation:** https://github.com/fury-gl/fury/pull/452
  - **Ball and Stick Representation demo:** https://github.com/fury-gl/fury/pull/452

- **Stick Representation**

  Stick model is a special case of Ball and Stick model where atomic radius of all
  molecules is set equal to the radius of tubes used to create bonds. It was created
  by using  `vtkOpenGLMoleculeMapper`. Tests were added and a demo was created.

  *Pull Requests:*

  - **Stick Representation:** https://github.com/fury-gl/fury/pull/452
  - **Stick Representation demo:** https://github.com/fury-gl/fury/pull/452

- **Sphere Representation**

  In chemistry, a space-filling model, also known as a calotte or sphere model, is a
  type of three-dimensional (3D) molecular model where the atoms are represented by
  spheres whose radii are proportional to the radii of the atoms. It was created by
  using `vtkOpenGLMoleculeMapper`. Tests were added and a demo was created.

  *Pull Requests:*

  - **Sphere Representation:** https://github.com/fury-gl/fury/pull/452
  - **Sphere Representation demo:** https://github.com/fury-gl/fury/pull/452



Objectives in Progress
----------------------

-  **Molecular Surfaces**

   There are three types of molecular surfaces:

   - Van der Waals
   - Solvent Accessible
   - Solvent Excluded

   Currently the first two molecular surfaces i.e. Van der Waals and Solvent
   Accessible are implemented. The code is based on the paper "Generating
   Triangulated Macromolecular Surfaces by Euclidean Distance Transform" by
   Dong Xu and Yang Zhang.

   *Pull Requests:*

   - **Molecular Surfaces Implementation:** https://github.com/fury-gl/fury/pull/492


Other Objectives
----------------

-  **2D Animated Surfaces**

   This was a simple demonstration that animated Two-Dimensional (2D) functions using FURY. 
   Created a grid of x-y coordinates and mapped the heights (z-values) to the corresponding x, y 
   coordinates to generate the surfaces. Used colormaps to color the surfaces.

   *Pull Requests:*

   - **Animated Surfaces:**  https://github.com/fury-gl/fury/pull/362

-  **Updated miscellaneous animations**

   -  Updated the demo of helical motion to stop using multiple line actors as discussed in the meeting.
   -  Updated the demo of brownian motion to make it more scientifically useful (removed unnecessary rotation of camera 
      during animation and box actor).
   -  Display simulation data for brownian motion and helical motion animations (number of simulated steps for brownian 
      motion and velocity of the particle for helical motion). 
   -  Created utility functions to make the code understandable and used these in emwave, helical and brownian 
      animations.

   *Pull Requests:*

   - **Updated helical, brownian, emwave animations:**  https://github.com/fury-gl/fury/pull/462

-  **GSoC weekly Blogs**

    Weekly blogs were added for FURY's Website.

    *Pull Requests:*

    - **First Evaluation:** https://github.com/fury-gl/fury/pull/475

    - **Second Evaluation:** https://github.com/fury-gl/fury/pull/493

Timeline
--------

+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Date                  | Description                                                      | Blog Link                                                                                          |
+=======================+==================================================================+====================================================================================================+
| Week 1(08-06-2021)    | Welcome to my GSoC Blog!                                         | `Weekly Check-in #1 <https://blogs.python-gsoc.org/en/suntzunamis-blog/weekly-check-in-1-11/>`__   |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 2(14-06-2021)    | First Week of coding: sphere model.                              | `Weekly Check-in #2 <https://blogs.python-gsoc.org/en/suntzunamis-blog/weekly-check-in-2-11/>`__   |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 3(21-06-2021)    | Bonding algorithms, Ball and Stick model progress.               | `Weekly Check-in #3 <https://blogs.python-gsoc.org/en/suntzunamis-blog/weekly-check-in-3-13/>`__   |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 4(28-06-2021)    | VTK molecular visualization classes.                             | `Weekly Check-in #4 <https://blogs.python-gsoc.org/en/suntzunamis-blog/weekly-check-in-4-14/>`__   |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 5(05-07-2021)    | Genesis of `molecular` module.                                   | `Weekly Check-in #5 <https://blogs.python-gsoc.org/en/suntzunamis-blog/weekly-check-in-5-13/>`__   |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 6(12-07-2021)    | Ribbon representation, updated `molecular` module (more pythonic)| `Weekly Check-in #6 <https://blogs.python-gsoc.org/en/suntzunamis-blog/weekly-check-in-6-18/>`__   |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 7(19-07-2021)    | More features to `molecular`, updated misc. animations.          | `Weekly Check-in #7 <https://blogs.python-gsoc.org/en/suntzunamis-blog/weekly-check-in-7-16/>`__   |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 8(26-07-2020)    | Ribbon to `molecular`, tests for `molecular`, animated surfaces. | `Weekly Check-in #8 <https://blogs.python-gsoc.org/en/suntzunamis-blog/weekly-check-in-8-11/>`__   |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 9(02-08-2021)    | Optimized `molecular` with mentors, GSoC blogs to FURY docs.     | `Weekly Check-in #9 <https://blogs.python-gsoc.org/en/suntzunamis-blog/weekly-check-in-9-11/>`__   |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 10(09-08-2021)   | Bounding box, `molecular` tutorial, molecular surfaces progress. | `Weekly Check-in #10 <https://blogs.python-gsoc.org/en/suntzunamis-blog/weekly-check-in-10-11/>`__ |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 11(16-08-2021)   | Molecular Surfaces (VDW, SAS) implementation.                    | `Weekly Check-in #11 <https://blogs.python-gsoc.org/en/suntzunamis-blog/weekly-check-in-11-9/>`__  |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+



Detailed weekly tasks and work done can be found
`here <https://blogs.python-gsoc.org/en/suntzunamis-blog/>`_.
ComboBox2D, TextBlock2D, Clipping Overflow.
===========================================

.. post:: July 19 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

Hello and welcome to my 8th weekly check-in. I will be sharing my progress with ``ComboBox2D`` and ``TextBlock2D`` UI components. After solving ``TextBlock2D`` sizing issue, I was finally able to complete the implementation of combo box. I also faced some problems while refactoring ``TextBlock2D``. The official repository of my sub-org, FURY can always be found `here <https://github.com/fury-gl/fury/>`_.

What did you do this week?
--------------------------
This week I finished the implementation of ``TextBlock2D`` and ``ComboBox2D``. I also added necessary tests and tutorials to aid the implementation. It still needs some color design, which will be decided later on. Here's an overview of the tutorial:

.. image:: https://user-images.githubusercontent.com/29832615/87884567-a8f19d80-ca2c-11ea-8fd2-e7e37b30602f.gif

Its a simple tutorial where I change the color of the label based on the option selected on the combo box.

Now while refactoring TextBlock2D, I noticed that the Text Overflow workaround on ListBox2D broke and started acting weird. This was mainly because the size implementation was changed and the code that I wrote a few months back wasn't relevant anymore. So I removed the previous code and had to re-implement the same. I later found out that I would be needing such an implementation for most of my UI components so I decided to create a separate method for clipping overflowing texts.

I had to figure out a specific algorithm to achieve this as the height and width of each characters were different and the combined width of multiple characters were not equal to the sum of their widths. I decided to go for a binary search implementation as it was faster compared to a simple linear checking algorithm. The said algorithm works as expected and is more effective compared to its previous implementation. I tweaked the previous combo box example to showcase this algorithm.

.. image:: https://user-images.githubusercontent.com/29832615/87884568-aabb6100-ca2c-11ea-9ab8-b05bdb8b0631.gif

The extremely long text is now clipped correctly.

What is coming up next week?
----------------------------
Next week, I have a couple of things to work on. Firstly, the single actor wall brick simulation is yet to be completed. Once I am done with that I will continue working on Tab UI and try to finish its skeletal implementation.

Did you get stuck anywhere?
---------------------------
I did get stuck regarding the previous implementation of Text Overflow in ListBox2D. Before the current implementation it looked something like this:

.. image:: https://user-images.githubusercontent.com/29832615/87430848-6b8fa900-c603-11ea-87b8-327f6e7f2ee0.png

Apart from that, I did not face any major issues.

``Thank you for reading, see you next week!!``
Multiple SDF primitives.
============================

.. post:: July 13 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

Make sure to check out Project `FURY <https://github.com/fury-gl/fury>`_

Hey Everyone! 
This week, multiple SDF primitives.

What did you do this week?
--------------------------
The objective of this week was to understand the capabilities of the multiple SDF primitives within the same cube project. To get a clearer understanding of what we can achieve with this approach, i added support for real time rotating primitives. The output of 2 cubes with rotating SDFs is shown below.

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-7.gif


The code for the above render is available at the `branch <https://github.com/lenixlobo/fury/tree/SDF-Experiments>`_

What is coming up next week?
----------------------------
The task for next week is to implement features in the SDF actor and also to fix minor bugs within the shaders. Once everything is done as planned, the PR will be merged for users to access the features.

Did you get stuck anywhere?
---------------------------
No issues this weekComboBox2D Progress!!
=====================

.. post:: June 14 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

Hello and welcome to my third weekly check-in, I will be sharing my progress with the project so far. In case you wondered, the sub-org that I am affiliated to is FURY. Make sure to check out the official repository `here <https://github.com/fury-gl/fury/>`_.

What did you do this week?
--------------------------
This week my objective was to work on the sizing and positioning issue regarding the sub-components of the ComboBox2D UI element. After countless debugging sessions and my mentor's support, I was successfully able to fix the issue. I also developed helpful methods and callbacks for this element to allow users to interact with them in a user friendly manner. I also wrote tests for the same. So far the said progress can be summarized via this gif:

.. image:: https://user-images.githubusercontent.com/29832615/84592637-cc8d5b00-ae64-11ea-9ff3-c1ce2095f7f2.gif

What is coming up next week?
----------------------------
Unfortunately, one of the sub-components ``TextBlock2D``, didn't have a resizing feature that I could make use of for this new UI component. Thus, I need to add that feature on a different PR. This feature will also be required while building other UI elements therefore adding this feature is currently my top priority. There's also a bug regarding the scrollbar that needs to be fixed. The scrollbar overshoots the specified combobox area after new items are appended to the combobox's list of items during runtime. Hopefully I will be able to get them done by next week.

Did you get stuck anywhere?
---------------------------
I was really confused with the coordinate system of Panel2D that was the main reason why components were misplaced. I also faced some issues regarding the scrollbar callback as it wasn't scrolling the elements properly, the items were stuck in place. So far I was able to fix them. Apart from these, I didn't face any major issues.

``Thank you for reading, see you next week!!``
Part of the Journey is the end unless its Open Source!
======================================================

.. post:: August 23 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

Hello and welcome to my final weekly check-in. Today officially marks the end of  the coding period for GSoC 2020. I enjoyed every bit of it. This was a life-changing experience for me and now I observe and interpret everything from a different perspective altogether. I have grown into a better developer and a person since GSoC. I would like to thank all my mentors and especially Serge for his immense support and mentorship. I would love to contribute to fury even after GSoC is over but unfortunately my semester break is over so I wont be as active as I was during the summer.

Now, regarding work I will be sharing my progress with the File Dialog UI component. The official repository of my sub-org can be found `here <https://github.com/fury-gl/fury/>`_.

What did you do this week?
--------------------------
This week I worked on the File Dialog UI component. Fury previously had a FileMenu component which could browse through the file system but we needed a dialog like implementation for it so that its easier for the users to read and write files during runtime. I tried implementing a simple design for it. It specifically has two modes, one for saving files and the other for writing files. The implementation can be demonstrated as follows:

Open Dialog:
^^^^^^^^^^^^

.. image :: https://user-images.githubusercontent.com/29832615/90978632-df12c780-e56c-11ea-8517-6243ea06bdd2.gif

Save Dialog:
^^^^^^^^^^^^

.. image :: https://user-images.githubusercontent.com/29832615/90978638-eafe8980-e56c-11ea-835a-3a82ccee2973.gif

What is coming up next week?
----------------------------
Next week I will start with my final GSoC documentation and code submission. I will also try to implement the tests and tutorials for File Dialog or any further changes requested by my mentors. If I am not able to finish it within the next week, I will get it done after GSoC.

Did you get stuck anywhere?
---------------------------
I did not face any major issues this week.

``Thank you all for your love and support. ❤️😄``
More Shaders!!
=====================

.. post:: August 02 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

Make sure to check out Project `FURY <https://github.com/fury-gl/fury>`_

Hey ! 
This week, More Shaders!

What did you do this week?
--------------------------
After merging the SDF actor last week , the mentors and I decided to work on another project idea which was discussed prior to GSoC . So the next step assigned to me was to look into the current shader framework and create immersive visualizations using shaders. Being interested in volumetric rendering I looked into volumetric cloud rendering algorithms and created a basic shader-based render in FURY.


The output for the same is given below:

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-10a.gif


.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-10b.gif

The shader demos are available `here <https://github.com/lenixlobo/fury/tree/shader-demos>`_

What is coming up next week?
----------------------------
The main focus now is to discuss and work on shader demos which an demonstrate the capabilities of FURY 

Did you get stuck anywhere?
---------------------------
No issues were faced this weekTab UI, TabPanel2D, Tab UI Tutorial.
====================================

.. post:: July 26 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

Hello and welcome to my 9th weekly check-in. I will be sharing my progress with ``TabUI`` and its corresponding tutorial. The official repository of my sub-org, FURY can always be found `here <https://github.com/fury-gl/fury/>`_.

What did you do this week?
--------------------------
This week I finished the basic implementation of ``TabUI``. Apart from that I was also able to finish the tutorial which showcases different features of this said UI. With the help of this UI one can have multiple panels containing different UI elements within them. The basic implementation can be demonstrated as follows:

.. image:: https://user-images.githubusercontent.com/29832615/88484746-87456880-cf8e-11ea-9e96-9cba111b90d3.gif

After finishing with the basic implementation I moved on to create its tutorial. For that, I decided to combine 3 of the existing UI tutorials to be rendered with the help of Tab UI. I implemented the following in individual tabs:

- Controlling a cube with the help of ``LineSlider2D`` and ``RingSlider2D``.
- Using a ``CheckBox`` to render a cylinder or a sphere or both.
- Using ``ComboBox`` to set the color of a label.

The above tutorial can be demonstrated as follows:

.. image:: https://user-images.githubusercontent.com/29832615/88481324-6a9e3600-cf78-11ea-8c5b-e26bf158388a.gif

What is coming up next week?
----------------------------
Next week I will continue working on the Physics engine integration. Previously we were facing a problem regarding the uncontrollable spinning of objects rendered by a single actor. There must be an issue with mismatching axes alignment of FURY and pyBullet. The spinning problem is as following:

.. image:: https://user-images.githubusercontent.com/29832615/88485303-87476780-cf92-11ea-850c-cc63e1376ef8.gif

Did you get stuck anywhere?
---------------------------
I did get stuck with the collapsing functionality of Tab UI and the uncontrollable spinning of the bricks in the Physics simulation. Apart from that I did not face any major issues.

``Thank you for reading, see you next week!!``
Spherical harmonics
===========================

.. post:: June 28 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

Make sure to check out Project `FURY <https://github.com/fury-gl/fury>`_

Hey ! 
This week, Spherical harmonics!

What did you do this week?
--------------------------
The main task for the week was to include an implementation of spherical harmonics (upto the order of 4) as a FURY actor. This was the initial milestone to be achieved to work towards the support of using spherical harmonics as an visualization technique. I have added the GIFs for both the renders below. I also worked on a occlusion based lighting model.

Spherical harmonics for different values of order and degree:

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-5a.gif

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-5b.gif


What is coming up next week?
----------------------------
The next task to add support for the user to be able to render different spherical harmonics by passing arguments

Did you get stuck anywhere?
---------------------------
Spherical harmonics involve a lot of complicated math behind the hood. So the initial days were spent understanding the math .I was confused while working on the implementation but eventually got it working.
Translation, Reposition, Rotation.
==================================

.. post:: July 5 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

Hello and welcome to my 6th weekly check-in. The first evaluation period officially ends and I am very excited to move on to the second coding period. I will be sharing my progress with handling specific object's properties among various multiple objects rendered by a single actor. I am mainly focusing on making it easier to translate, rotate and reposition a particular object, so that I can use them to render physics simulations more efficiently. The official repository of my sub-org, FURY can always be found `here <https://github.com/fury-gl/fury/>`_.

What did you do this week?
--------------------------
Last week I worked on physics simulations rendered in FURY with the help of pyBullet. Now the simulations were highly un-optimized, specially the brick wall simulation as each brick was rendered by its own actor. In other words, 1 brick = 1 actor. Now my objective was to render all the bricks using a single actor, but before jumping into the simulation I had to figure out how to modify specific properties of an individual object. Thanks to my mentor's `PR <https://github.com/fury-gl/fury/pull/233>`_, I was able to experiment my implementations quickly.

Translation:
^^^^^^^^^^^^

.. image:: https://user-images.githubusercontent.com/29832615/86536066-5085b080-bf02-11ea-9bcd-9e555adc2ca1.gif

The algorithm behind translation is to first identify the vertices of the object, then bring the vertices to the origin by subtracting their centers and then adding the displacement vector. The said operation can be achieved by the following snippet:

.. code-block::

    # Update vertices positions
    vertices[object_index * sec: object_index * sec + sec] = \
        (vertices[object_index * sec: object_index * sec + sec] - centers[object_index]) + transln_vector​


Rotation:
^^^^^^^^^

.. image:: https://user-images.githubusercontent.com/29832615/86536065-4fed1a00-bf02-11ea-815d-f7f297165c53.gif

The algorithm behind rotation is to first calculate the difference between the vertices and the center of the object. Once we get the resultant matrix, we matrix multiply it with the rotation matrix and then we further add the centers back to it so that we preserve the position of the object. Rotation matrix can be defined as:

.. image:: https://wikimedia.org/api/rest_v1/media/math/render/svg/242deb7010fd504134a6cacab3d0ef4ce02e7613

where gamma, beta and alpha corresponds to the angle of rotation along Z-axis, Y-axis and X-axis.

.. code-block:: python

    def get_R(gamma, beta, alpha):
        """ Returns rotational matrix.
        """
        r = [
            [np.cos(alpha)*np.cos(beta), np.cos(alpha)*np.sin(beta)*np.sin(gamma) - np.sin(alpha)*np.cos(gamma),
            np.cos(alpha)*np.sin(beta)*np.cos(gamma) + np.sin(alpha)*np.sin(gamma)],
            [np.sin(alpha)*np.cos(beta), np.sin(alpha)*np.sin(beta)*np.sin(gamma) + np.cos(alpha)*np.cos(gamma),
            np.sin(alpha)*np.sin(beta)*np.cos(gamma) - np.cos(alpha)*np.sin(gamma)],
            [-np.sin(beta), np.cos(beta)*np.sin(gamma), np.cos(beta)*np.cos(gamma)]
        ]
        r = np.array(r)
        return r

    vertices[object_index * sec: object_index * sec + sec] = \
        (vertices[object_index * sec: object_index * sec + sec] -
        centers[object_index])@get_R(0, np.pi/4, np.pi/4) + centers[object_index]


Reposition:
^^^^^^^^^^^

.. image:: https://user-images.githubusercontent.com/29832615/86536063-4ebbed00-bf02-11ea-8592-a695d7b91426.gif

Repositioning is similar to that of translation, except in this case, while repositioning we update centers with the new position value.

.. code-block:: python

    new_pos = np.array([1, 2, 3])

    # Update vertices positions
    vertices[object_index * sec: object_index * sec + sec] = \
        (vertices[object_index * sec: object_index * sec + sec] -
        centers[object_index]) + new_pos

    centers[object_index] = new_pos

What is coming up next week?
----------------------------
Currently, I am yet to figure out the orientation problem. Once I figure that out I will be ready to implement simulations without any major issues. I am also tasked with creating a wrecking ball simulation and a quadruped robot simulation.

Did you get stuck anywhere?
---------------------------
I did face some problems while rotating objects. My mentors suggested me to implement it via rotation matrix. I still haven't figured out the orientation problem, which I plan to work on next. Apart from these I did not face any major issues.

``Thank you for reading, see you next week!!``
First week of coding!!
======================

.. post:: June 7 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

Hey everyone!
This week : Geometry Shaders!

What did you do this week?
--------------------------
To get a better understanding of the working of the shader pipeline, the mentors assigned me a challenging task of implementing a Dynamic Texture. The basic idea is to create a 'volumetric' texture by stacking layer of textures. Such an example is an ideal use case for a geometry shader. Since i had not much prior experience with Geometry shaders before, i spent the initial days going through existing implementations of similar ideas in OpenGL/DirectX.
After working on the code, the final image rendered is given below.

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-2.png

I created a PR for the fur texture which is available `here <https://github.com/lenixlobo/fury/blob/Dynamic-Texture/docs/experimental/viz_dynamictex.py>`_


What is coming up next week?
----------------------------
The current PR has some minor bugs which need to be fixed. The next step would be to review the code and find the solutions for the bugs. Also we are looking into ideas on optimization for faster rendering time.

The next week will be spent looking into ray marching algorithms and adding them to the current code base as possible alternatives for FURY Actor primitives.

Did you get stuck anywhere?
---------------------------
Nothing major.
May the Force be with you!!
===========================

.. post:: June 28 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

Hello and welcome to my 5th weekly check-in, I will be sharing my progress of pyBullet Physics Engine Integration with FURY. The official repository of my sub-org, FURY can always be found `here <https://github.com/fury-gl/fury/>`_.

What did you do this week?
--------------------------
Last week due to the ``vtkTextActor`` sizing issue I was not able to continue my work with ``ComboBox2D`` UI element. Thus, I decided to work on the "Physics Engine Integration" part of my project. It took me quite a while to understand the terminology of various methods and algorithms used for detection and physics simulation. Nevertheless, I was able to create a couple of rigid body examples to showcase the integration procedure. For physics calculations we used pyBullet and for rendering we used FURY. In other words, pyBullet will handle the backend part of the simulation and FURY will handle the frontend. I have documented the entire integration process `here <https://docs.google.com/document/d/1XJcG1TL5ZRJZDyi8V76leYZt_maxGp0kOB7OZIxKsTA/edit?usp=sharing>`__ in detail.

Ball Collision Simulation:
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: https://user-images.githubusercontent.com/29832615/85949638-988e5b80-b975-11ea-85ba-16c1f78dec89.gif

For the first example, I created a simple collision between two balls in which two spheres were created both in FURY and pyBullet world and then both the worlds are connected by syncing the position and orientation of the said bodies. Timer callback is created to update the positions and to simulate Physics for each step.

Brick Wall Simulation:
^^^^^^^^^^^^^^^^^^^^^^

.. image:: https://raw.githubusercontent.com/Nibba2018/testing-fury/master/2020-06-26_21-15-47.gif

For the second example I tried to increase the complexity of the simulations by increasing the number of dynamic objects. Here a brick-wall is created using 100 bricks and then a ball is thrown at it. The rest of it is simulated. The same concepts from the first example is used to render the second one. Timer callback is created to update position of all the objects and to simulate Physics for each step.

What is coming up next week?
----------------------------
In the above examples I used a separate actor for each object which is highly un-optimized. Thus, my mentor suggested me to render all the bricks using a single actor, so that the simulation is more optimized. I am not very confident in changing the position and orientation of different objects rendered by a single actor. Therefore, I will have to research a bit more about it. Apart from that I will also try to work on Constraint based simulation examples if possible.

Did you get stuck anywhere?
---------------------------
The pyBullet documentation isn't well written for cross rendering, hence it was a bit of a challenge for me to implement the integration. I also faced a problem regarding the offset of actors between the FURY world and pyBullet world. Both use different coordinate systems for rendering and simulation because of which I had a constant offset between the objects during simulations. I was able to fix it by converting one coordinate system to the other. Apart from this I did not have any major issues.

``Thank you for reading, see you next week!!``
TextBlock2D Progress!!
======================

.. post:: June 21 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

Hello and welcome to my 4th weekly check-in, I will be sharing my progress with TextBlock2D UI component. The official repository of my sub-org, FURY can always be found `here <https://github.com/fury-gl/fury/>`_.

What did you do this week?
--------------------------
This week my objective was to work on the ``TextBlock2D`` resizing feature and also work on the scrollbar bug in ``ComboBox2D`` UI component. The ``TextBlock2D`` component makes use of ``vtkTextActor`` and its methods available in VTK for handling all Text related properties such as font size, color, background color etc. Now, going through the official `docs <https://vtk.org/doc/nightly/html/classvtkTextActor.html>`_ I figured out that the there were three different scaling modes for the said actor:

- TEXT_SCALE_MODE_NONE
- TEXT_SCALE_MODE_PROP
- TEXT_SCALE_MODE_VIEWPORT

The first scaling mode was currently implemented in FURY's ``TextBlock2D`` code base. This was done mainly because the only way for the user to create a text block was by providing the font size as parameter. We couldn't use the third option as the FURY API tries to maintain abstraction by not exposing ``vtkViewport`` parameter. Now in order to allow size as a parameter I had to use the second scaling option which is ``TEXT_SCALE_MODE_PROP``. With  this mode one can set the ``Position2`` parameter to the desired value of height and width. But the problem is I cannot use this for the combobox element as the background size will be inconsistent with respect to font size and text length.

.. image:: https://user-images.githubusercontent.com/29832615/85226809-12af6500-b3f7-11ea-8a75-7e86d40701d1.png

.. image:: https://user-images.githubusercontent.com/29832615/85226812-14792880-b3f7-11ea-8edd-1df25382a48f.png

Therefore, as a solution we agreed to add a separate Rectangle2D UI component as the background for Text Block along with vtkTextActor. With this implementation one can easily manipulate the background irrespective of the text properties. But this had another issue, the user now had 2 different ways of initializing a TextBlock2D. Therefore, when the user uses size a constructor parameter, then everything works in sync, but the same is not true for the font size parameter. This is mainly because we do not have a specific way of determining the size of the actor based on font size. My mentors have agreed to look into this problem and advised me to focus on something else instead.

What is coming up next week?
----------------------------
Next week I am planning to start with Physics engine integration of FURY with pyBullet. My main focus would be to create simple simulations first before initiating the integration process. If the size issue is solved before I move into Physics engine then I would complete the ComboBox UI and start with Tab UI instead. I have also fixed the scrollbar bug.

Did you get stuck anywhere?
---------------------------
The main problem that I faced so far was regarding the variable background size of the Text Block component. Also the ``GetPosition2`` method of ``vtkTextActor`` would always return the same value irrespective of the font size parameter passed as a constructor argument. Therefore, it was not possible for me to determine the actual size or bounding box of the said text actor.

``Thank you for reading, see you next week!!``
First week of coding!!
======================

.. post:: June 7 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

Hello and welcome to my second weekly check-in, I will be sharing my progress for the first week of coding.

What did you do this week?
--------------------------
This week my objective was to create the skeleton model for ComboBox2D UI element along with its required sub-components and their respective default callbacks. So far, I have successfully been able to implement them. Apart from some issues that I was facing regarding the sizing and position of the sub-components.

What is coming up next week?
----------------------------
Next week, my main focus will be on fixing the sizing issue and I will also be working on respective default and specialized class methods required by the component and its sub-components for proper functioning. I will also create tests for the same.

Did you get stuck anywhere?
---------------------------
The main problem that I faced while coding the UI element was regarding the sizing and positioning of the sub-components. The sub-components appeared flipped and the TextBlock2D element had no option to resize based on a particular width & height definition.

.. image:: https://user-images.githubusercontent.com/29832615/82089691-b713fc80-9711-11ea-8d72-efaeaac1044c.png

As you can see, the ComboBox UI element doesn't look the way its supposed to. Along with this issue the scrollbars callback doesn't work either. Hopefully, I will be able to get these issues fixed by next week.

``Thank you for reading, see you next week!!``
FURY 0.5.1 Released
===================

.. post:: April 9 2020
   :author: skoudoro
   :tags: fury
   :category: release


The FURY project is happy to announce the release of FURY 0.5.1!
FURY is a free and open source software library for scientific visualization and 3D animations.

You can show your support by `adding a star <https://github.com/fury-gl/fury/stargazers>`_ on FURY github project.

The **major highlights** of this release are:

.. include:: ../../release_notes/releasev0.5.1.rst
    :start-after: --------------
    :end-before: Details

.. note:: The complete release notes are available :ref:`here <releasev0.5.1>`

**To upgrade or install FURY**

Run the following command in your terminal::

    pip install --upgrade fury

or::

    conda install -c conda-forge fury


**Questions or suggestions?**

For any questions go to http://fury.gl, or send an e-mail to fury@python.org
We can also join our `discord community <https://discord.gg/6btFPPj>`_

We would like to thanks to :ref:`all contributors <community>` for this release:

.. include:: ../../release_notes/releasev0.5.1.rst
    :start-after: commits.
    :end-before: We closed


On behalf of the :ref:`FURY developers <community>`,

Serge K.
Weekly Check-in #1
==================

.. post:: May 30 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

Welcome to my GSoC Blog!!!
--------------------------
Hey everyone!
This is my blog for this summer’s GSoC @ PSF
I am Lenix Lobo, an undergraduate student from India, and this summer I will be working with project Fury under the umbrella of the Python Software foundation.  I will be working on improving the current shader framework.

What did you do during the Community Bonding Period?
----------------------------------------------------
Since most of the places including my university are closed due to the pandemic outbreak, I decided to get a head start and start with the project early. During the community bonding period, I had video conference meetings with my mentors scheduled every week on Wednesday. During these meetings i interacted with the mentors to have a coherent understanding of how the project design and implementation will be managed over the course of the entire period. 

Since my project involves a lot of theoretical understanding of concepts such as ray marching, I spent the period going through the theory of each topic.This week also involved going through the documentation for shaders used in VTK. 

What is coming up next week?
----------------------------
The next task assigned to me is to go through the theory of geometry shaders and to create a example using the same. 
Spherical harmonics, Continued.
==================================

.. post:: July 5 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

Make sure to check out Project `FURY <https://github.com/fury-gl/fury>`_

Hey ! 
Spherical harmonics, Continued!

What did you do this week?
--------------------------
Last week I added a basic implementation of Spherical harmonics based actors. However, the implementation was quite restricted and we needed to add support for more accurate generation of spherical harmonics. So the task assigned this week was to implement the spherical harmonics function within the shader rather than passing variables as uniforms. This was quite an challenging task as it involved understanding of mathematical formulae and implementing them using existing GLSL functions.
The output of the implementation is shown below :

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-6.gif

While , i was able to complete the task the frame rate for the generated output was quite lower than expected. 
The code for the above render is available at the `branch <https://github.com/lenixlobo/fury/tree/Spherical-Harmonics>`_

What is coming up next week?
----------------------------
The next task is to discuss possible performance improvements with the mentors and also look into alternative ideas to add spherical harmonics as actors in FURY.

Did you get stuck anywhere?
---------------------------
Spherical harmonics involve a lot of complicated math behind the hood as a result the generated output has a very poor frame rate. Currently, we are looking into improving this.
Improvements in SDF primitives.
===========================================

.. post:: July 20 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

Make sure to check out Project `FURY <https://github.com/fury-gl/fury>`_

Hey Everyone! 
This week, Improvements in SDF primitives.

What did you do this week?
--------------------------
Over the past few weeks i have been working on adding SDF based actors in FURY. The task for this week was to complete the remaining features such as support for scaling and rotation based on user direction for the actors. The main objective is to complete the SDF actor and make it ready to merge into the FURY codebase. Below are the outputs after added the improvements to the SDF primitives.

Scaling:

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-8b.gif

Rotation based on passed directions:

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-8a.gif

The code for the above render is available at the `branch <https://github.com/lenixlobo/fury/tree/SDF-Experiments>`_

What is coming up next week?
----------------------------
Since the features are almost done the task for next week is clean the code and test for bugs. And then to eventually merge the  sdf_actor into the fury codebase.

Did you get stuck anywhere?
---------------------------
No major issues this week.
Outline Picker
=====================

.. post:: August 17 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

Make sure to check out Project `FURY <https://github.com/fury-gl/fury>`_

Hey ! 
This week, Picking Outline!

What did you do this week?
--------------------------
We needed a Outline feature in FURY to indicate which model we choose in the scene. So the task assigned was to find options to achieve this. There were two ways to do this, 1. Using shader and 2. Using Vtk PolyData Silhouette. Despite trying multiple implementation methods the shader approach was not working . I also tried using VTKs inbuilt function , but there is a bug when i use some models. When i choose a model, it renders outline for every polygon , which is not what we want to achieve. The bug is shown below:


Below are the outputs of the techniques i worked on :

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-12.gif


The shader demos are available `here <https://github.com/lenixlobo/fury/tree/shader-demos>`_

What is coming up next week?
----------------------------
With the end of GSoC approaching soon, the next task is to create a PR which can help new users to test different shaders using UI to get started. 

Did you get stuck anywhere?
---------------------------
I still was not able to figure out how we can achieve the outline effect. Am currently looking into other approaches we could useWrecking Ball Simulation, Scrollbars Update, Physics Tutorials.
===============================================================

.. post:: August 16 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

Hello and welcome to my 12th weekly check-in. In this blog I will be discussing my progress with the wrecking ball simulation and my scrollbar separation work. Apart from this I have also finalized the physics simulation tutorials and have created a Pull Request to finally get it merged with the official repository. The official repository of my sub-org can be found `here <https://github.com/fury-gl/fury/>`_.

What did you do this week?
--------------------------
This week I was mainly focusing on the wrecking ball simulation. This simulation is basically the combination of chain simulation and brick wall simulation. A sphere attached to a chain smashes a "NxNxN" brick wall. The simulation is as follows:

.. image :: https://user-images.githubusercontent.com/29832615/90336291-84232280-dff8-11ea-869b-21a99b203c31.gif

There's a rendering bug with the cylinders because of which the chain segments look weird. My mentors confirmed that this bug originated from VTK's `cylinder source` method and they are currently working on it to fix it. The simulation will render correctly once that bug is fixed.

Regarding the scrollbar separation task, I was able to fix those callback issues that I was facing. The mouse callbacks on the scrollbar now work correctly:

.. image :: https://user-images.githubusercontent.com/29832615/90337280-1af2dd80-dfff-11ea-94c4-508121307583.gif

I have also created a pull request to add the following physics simulations with proper documentation to the main repository:

- Brick Wall Simulation
- Ball Collision Simulation
- Chain Simulation
- Wrecking Ball Simulation

What is coming up next week?
----------------------------
Currently I am a bit confused with the implementation of scrollbars with UI components. I had a meeting with my mentor and he decided to help me out with this. So most probably I will be working with the scrollbar component and its implementation. Next week will also be the final week for GSoC 2020 before the evaluations start so I would work on getting the final documentation and formalities ready.

Did you get stuck anywhere?
---------------------------
Apart from the scrollbar implementation idea, I did not face any major issues.

``Thank you for reading, see you next week!!``
Merging SDF primitives.
===========================================

.. post:: July 27 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

Make sure to check out Project `FURY <https://github.com/fury-gl/fury>`_

Hey Everyone! 
This week, Merging SDF primitives.

What did you do this week?
--------------------------
Since GSoC started I have been working on adding support for raymarching based SDF actors as primitives in the FURY codebase. This week with the release of FURY 0.6.0 , the task assigned to me was to complete the remaining parts of the SDF actor including tests and tutorial. THe SDF actor is now part of the FURY actor and can be accessed using sdf_actor.
Currently we support , ellipsoids, spheres and torus as primitive options. As expected, SDF based actors have shown tremendous performance improvements over traditional polygon based actor.

Despite using 100,000 torus the FPS is higher than 60 :

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-9.gif

10,000 actors :

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-9b.gif

I also made a tutorial for new users to get started `here <https://fury.gl/latest/auto_tutorials/04_shaders/viz_sdfactor.html#sphx-glr-auto-tutorials-04-shaders-viz-sdfactor-py>`_

What is coming up next week?
----------------------------
Now that the SDF actor is merged , the next step is to focus on spherical harmonics and i will also be working on creating shader visualization to showcase the features of FURY

Did you get stuck anywhere?
---------------------------
This week involved alot of work , including making tests, tutorial and looking for bugs but everything went smoothly .
.. image:: https://developers.google.com/open-source/gsoc/resources/downloads/GSoC-logo-horizontal.svg
   :height: 50
   :align: center

.. image:: https://www.python.org/static/community_logos/python-logo.png
   :width: 40%

.. image:: https://python-gsoc.org/logos/FURY.png
   :height: 30



Google Summer of Code 2020 Final Work Product
=============================================

.. post:: August 24 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

-  **Name:** Lenix Lobo
-  **Organisation:** Python Software Foundation
-  **Sub-Organisation:** FURY
-  **Project:** `FURY - Improve Shader Framework <https://github.com/fury-gl/fury/wiki/Google-Summer-of-Code-2020>`__

Introduction
------------
The current shader framework for FURY is based on VTK and lacks documentation to get started which can be overwhelming for new users. The objective of this project is to enable users to be easily able to understand and use the shader framework to render stunning visual representations of data. The project involves programming vertex and fragment shaders to generate effects for more immersive visualization.

Proposed Objectives
-------------------
**Adding SDF actor to the API**

This actor uses raymarching to model actors using SDF. The method provides several actors including `ellipsoid`, `sphere` and `torus`.
**Shader demos**

 Use the FURY shader system to create and visualize different shading algorithms. Implementations include `SphereMap`, `Toon`, `Gooch` and `Vertex Noise`

Unsubmitted Functionalities
---------------------------
**Spherical Harmonics using Shaders.**

The spherical harmonics algorithm is used to generate spherical surfaces using biases and coefficients computed. The general approach to achieve this is computationally expensive. The idea proposed was to leverage the GPU hardware using shaders to provide a faster more efficient alternative to the current implementations. The second month of the coding period was devoted to the same task but unfortunately, the observed performance was quite unsatisfactory than the expected performance. Moreover, the output shape of the geometry was distorted. It was then decided to continue the work after the GSoC period and prioritize the task at hand.

The Work in Progress can be accessed here. https://github.com/lenixlobo/fury/tree/Spherical-Harmonics

**Dynamic Texture using Geometry Shader**

Geometry Shaders provide a lot of flexibility to users to create custom geometry behaviors such as instancing. The idea was to create a dynamic Fur/Hair effect on top of a FURY actor. Unfortunately due to missing documentation on VTK geometry shaders and lack of resources, the project was not completed during the GSoC period. However, I will continue to try to solve the issue.

The source code for the current progress can be accessed here. https://github.com/lenixlobo/fury/tree/Dynamic-Texture


Objectives Completed
--------------------
**SDF based Actor**

  The objective here was to provide an alternative approach to users to use SDF modeled actors in the scene. This actor is modeled using the raymarching algorithm which provides much better performance than conventional polygon-based actors. Currently, the shapes supported include ellipsoid, sphere and torus

  *Pull Requests:*
  **SDF Actor method:** https://github.com/fury-gl/fury/pull/250

**Multiple SDF Actor**

  The objective was to create a method through which multiple SDF primitives are rendered within a single cube. This task helped us explore the limitations of the shader system and also benchmarking the performance.

  *Pull Requests:*
  **MultiSDF Shader:** https://github.com/fury-gl/fury/blob/master/docs/experimental/viz_multisdf.py

**Shader Demos**

  The task here was to create a pull request showcasing the capabilities of the FURY shader system and to also provide examples or new users to get started with integrating custom shaders into the scenes.

  *Pull Requests:*
  **Shader Demos:** https://github.com/fury-gl/fury/pull/296



Other Objectives
----------------
- **Tutorials**

   Create Tutorials for new users to get familiar with the Shader System

   *Pull Requests:*
   - **Shader UI Tutorial**

   https://github.com/fury-gl/fury/pull/296

   -**SDF Actor Tutorial**

   https://github.com/fury-gl/fury/pull/267

- **GSoC weekly Blogs**

  Weekly blogs were added for FURY's Website.

  *Pull Requests:*
  - **First & Second Evaluation:**

  https://github.com/fury-gl/fury/pull/250
  https://github.com/fury-gl/fury/pull/267

  - **Third Evaluation:**

  https://github.com/fury-gl/fury/pull/296


Timeline
--------

====================  ============================================================  ===========================================================================================
Date                  Description                                                   Blog Link
====================  ============================================================  ===========================================================================================
Week 1(30-05-2020)    Welcome to my GSoC Blog!                                      `Weekly Check-in #1 <https://blogs.python-gsoc.org/en/lenixlobos-blog/gsoc-blog-week-1/>`__
Week 2(07-06-2020)    Geometry Shaders!                                             `Weekly Check-in #2 <https://blogs.python-gsoc.org/en/lenixlobos-blog/weekly-check-in-week-2/>`__
Week 3(14-06-2020)    Ray Marching!                                                 `Weekly Check-in #3 <https://blogs.python-gsoc.org/en/lenixlobos-blog/weekly-check-in-week-3/>`__
Week 4(21-06-2020)    RayMarching Continued                                         `Weekly Check-in #4 <https://blogs.python-gsoc.org/en/lenixlobos-blog/weekly-check-in-week-4/>`__
Week 5(28-06-2020)    Spherical Harmonics                                           `Weekly Check-in #5 <https://blogs.python-gsoc.org/en/lenixlobos-blog/weekly-check-in-week-5/>`__
Week 6(05-07-2020)    Spherical Harmonics Continued                                 `Weekly Check-in #6 <https://blogs.python-gsoc.org/en/lenixlobos-blog/weekly-check-in-week-6/>`__
Week 7(12-07-2020)    Multiple SDF Primitives                                       `Weekly Check-in #7 <https://blogs.python-gsoc.org/en/lenixlobos-blog/weekly-check-in-week-7/>`__
Week 8(19-07-2020)    Improvements in SDF primitives                                `Weekly Check-in #8 <https://blogs.python-gsoc.org/en/lenixlobos-blog/weekly-check-in-week-8/>`__
Week 9(26-07-2020)    Merging SDF Actor and Benchmarks!                             `Weekly Check-in #9 <https://blogs.python-gsoc.org/en/lenixlobos-blog/weekly-check-in-week-9/>`__
Week 10(02-08-2020)   More Shaders                                                  `Weekly Check-in #10 <https://blogs.python-gsoc.org/en/lenixlobos-blog/weekly-check-in-week-10/>`__
Week 11(08-08-2020)   Even More Shaders                                             `Weekly Check-in #11 <https://blogs.python-gsoc.org/en/lenixlobos-blog/weekly-check-in-week-11/>`__
Week 12(16-08-2020)   Picking Outline                                               `Weekly Check-in #12 <https://blogs.python-gsoc.org/en/lenixlobos-blog/weekly-check-in-week-12/>`__
Week 13(23-08-2020)   Final Week                                                    `Weekly Check-in #13 <https://blogs.python-gsoc.org/en/lenixlobos-blog/weekly-check-in-week-13/>`__
====================  ============================================================  ===========================================================================================


Detailed weekly tasks and work done can be found
`here <https://blogs.python-gsoc.org/en/lenixlobos-blog/>`__.
FURY 0.6.1 Released
===================

.. post:: August 18 2020
   :author: skoudoro
   :tags: fury
   :category: release


The FURY project is happy to announce the release of FURY 0.6.1!
FURY is a free and open source software library for scientific visualization and 3D animations.

You can show your support by `adding a star <https://github.com/fury-gl/fury/stargazers>`_ on FURY github project.

This Release is mainly a maintenance release. The **major highlights** of this release are:

.. include:: ../../release_notes/releasev0.6.1.rst
    :start-after: --------------
    :end-before: Details

.. note:: The complete release notes are available :ref:`here <releasev0.6.1>`

**To upgrade or install FURY**

Run the following command in your terminal::

    pip install --upgrade fury

or::

    conda install -c conda-forge fury


**Questions or suggestions?**

For any questions go to http://fury.gl, or send an e-mail to fury@python.org
We can also join our `discord community <https://discord.gg/6btFPPj>`_

We would like to thanks to :ref:`all contributors <community>` for this release:

.. include:: ../../release_notes/releasev0.6.1.rst
    :start-after: commits.
    :end-before: We closed


On behalf of the :ref:`FURY developers <community>`,

Serge K.
FURY 0.6.0 Released
===================

.. post:: July 20 2020
   :author: skoudoro
   :tags: fury
   :category: release


The FURY project is happy to announce the release of FURY 0.6.0!
FURY is a free and open source software library for scientific visualization and 3D animations.

You can show your support by `adding a star <https://github.com/fury-gl/fury/stargazers>`_ on FURY github project.

The **major highlights** of this release are:

.. include:: ../../release_notes/releasev0.6.0.rst
    :start-after: --------------
    :end-before: Details

.. note:: The complete release notes are available :ref:`here <releasev0.6.0>`

**To upgrade or install FURY**

Run the following command in your terminal::

    pip install --upgrade fury

or::

    conda install -c conda-forge fury


**Questions or suggestions?**

For any questions go to http://fury.gl, or send an e-mail to fury@python.org
We can also join our `discord community <https://discord.gg/6btFPPj>`_

We would like to thanks to :ref:`all contributors <community>` for this release:

.. include:: ../../release_notes/releasev0.6.0.rst
    :start-after: commits.
    :end-before: We closed


On behalf of the :ref:`FURY developers <community>`,

Serge K.
.. image:: https://developers.google.com/open-source/gsoc/resources/downloads/GSoC-logo-horizontal.svg
   :height: 50
   :align: center
   :target: https://summerofcode.withgoogle.com/projects/#6653942668197888

.. image:: https://www.python.org/static/community_logos/python-logo.png
   :width: 40%
   :target: https://blogs.python-gsoc.org/en/nibba2018s-blog/

.. image:: https://python-gsoc.org/logos/FURY.png
   :width: 25%
   :target: https://fury.gl/latest/community.html

Google Summer of Code Final Work Product
========================================

.. post:: August 24 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

-  **Name:** Soham Biswas
-  **Organisation:** Python Software Foundation
-  **Sub-Organisation:** FURY
-  **Project:** `FURY - Create new UI Widgets & Physics Engine
   Integration <https://github.com/fury-gl/fury/wiki/Google-Summer-of-Code-2020>`_

Proposed Objectives
-------------------

-  ComboBox
-  Tab UI
-  File Dialog Improvements

Modified Objectives
-------------------

-  Combobox
-  Tab UI
-  File Dialog Improvements
-  Double Click Callback
-  TextBlock2D Improvements
-  Scrollbars as a Standalone Component
-  Physics Engine Integration

Objectives Completed
--------------------

-  **ComboBox2D UI Component**

   A combobox is a commonly used graphical user interface widget.
   Traditionally, it is a combination of a drop-down list or list box and a
   single-line textbox, allowing the user to select a value from the list.
   The term "combo box" is sometimes used to mean "drop-down list".
   Respective components, tests and tutorials were created.

   *Pull Requests:*

   -  **Combobox UI component:** https://github.com/fury-gl/fury/pull/240
   -  **Combobox UI Tutorial:** https://github.com/fury-gl/fury/pull/246

-  **Tab UI Component**

   In interface design, a tabbed document interface or Tab is a graphical
   control element that allows multiple documents or panels to be contained
   within a single window, using tabs as a navigational widget for
   switching between sets of documents. Respective components, tests and
   tutorials were created.

   *Pull Requests:*

   -  **Tab UI component:** https://github.com/fury-gl/fury/pull/252
   -  **Tab UI tutorial:** https://github.com/fury-gl/fury/pull/275

-  **Double Click Callback**

   Double click callbacks aren't implemented in VTK by default so they need
   to be implemented manually. With my mentor's help I was able to
   implement double click callbacks for all the three mouse buttons
   successfully.

   *Pull Requests:*

   -  **Adding Double Click Callback:**
      https://github.com/fury-gl/fury/pull/231

-  **TextBlock2D Improvements**

   The previous implementation of ``TextBlock2D`` was lacking a few
   features such as size arguments and text overflow. There was no specific
   way to create Texts occupying a said height or width area. Apart from
   that UI components like ``ListBoxItem2D``, ``FileMenu2D`` etc had an
   issue where text would overflow from their specified width. In order to
   tackle these problems, a modification was done to ``TextBlock2D`` to
   accept size as an argument and a new method was added to clip
   overflowing text based on a specified width and to replace the
   overflowing characters with ``...``.

   *Pull Requests:*

   -  **Setting Bounding Box for TextBlock2D:**
      https://github.com/fury-gl/fury/pull/248
   -  **Clip Text Overflow:** https://github.com/fury-gl/fury/pull/268

-  **Physics Engine Integration**

   Optional support for Physics engine integration of Pybullet was added to
   Fury. Pybullet's engine was used for the simulations and FURY was used
   for rendering the said simulations. Exhaustive examples were added to
   demonstrate various types of physics simulations possible using pybullet
   and fury. The said examples are as follows:

   -  Brick Wall Simulation

      -  Explains how to render and simulate external forces, objects and
         gravity.

   -  Ball Collision Simulation

      -  Explains how collisions work and how to detect said collisions.

   -  Chain Simulation

      -  Explains how to render and simulate joints.

   -  Wrecking Ball Simulation

      -  A more complicated simulation that combines concepts explained by
         the other examples.

   Apart from that, a document was created to explain the integration
   process between pybullet and fury in detail.

   *Pull Requests:*

   -  **Physics Simulation Examples:**
      https://github.com/fury-gl/fury/pull/287
   -  **Fury-Pybullet Integration Docs:**
      https://docs.google.com/document/d/1XJcG1TL5ZRJZDyi8V76leYZt_maxGp0kOB7OZIxKsTA/edit?usp=sharing

Objectives in Progress
----------------------

-  **Scrollbars as a standalone component**

   The previous implementation of scrollbars were hard coded into
   ``ListBox2D``. Therefore, it was not possible to use scrollbars with any
   other UI component. Apart from that, scrollbars in terms of design were
   limited. Creating a horizontal scrollbar was not possible. The objective
   of this PR is to make scrollbars separate so that other UI elements can
   also make use of it.

   Currently, the skeletal and design aspects of the scrollbars are
   implemented but the combination of scrollbars with other UI components
   are still in progress.

   *Pull Requests:*

   -  **Scrollbars as a Standalone API:**
      https://github.com/fury-gl/fury/pull/285

-  **File Dialog Improvements**

   Currently, we have access to ``FileMenu2D`` which allows us to navigate
   through the filesystem but it does not provide a user friendly Dialog to
   read and write files in Fury. Hence the idea is to create a file dialog
   which can easily open or save file at runtime. As of now, ``Open`` and
   ``Save`` operations are implemented. Corresponding tests and tutorials
   are in progress.

   *Pull Requests:*

   -  **File Dialog UI component:**
      https://github.com/fury-gl/fury/pull/294

Other Objectives
----------------

-  **Radio Checkbox Tutorial using FURY API**

   The objects for Radio button and Checkbox tutorial were rendered using
   VTK's method by a fellow contributor so I decided to replace them with
   native FURY API. The methods were rewritten keeping the previous commits
   intact.

   *Pull Requests:*

   -  **Radio Checkbox tutorial using FURY API:**
      https://github.com/fury-gl/fury/pull/281

-  **GSoC weekly Blogs**

   Weekly blogs were added for FURY's Website.

   *Pull Requests:*

   -  **First & Second Evaluation:**
      https://github.com/fury-gl/fury/pull/272
   -  **Third Evaluation:** https://github.com/fury-gl/fury/pull/286

Timeline
--------

+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Date                  | Description                                                      | Blog Link                                                                                          |
+=======================+==================================================================+====================================================================================================+
| Week 1(30-05-2020)    | Welcome to my GSoC Blog!!                                        | `Weekly Check-in #1 <https://blogs.python-gsoc.org/en/nibba2018s-blog/weekly-check-in-1-5/>`__     |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 2(07-06-2020)    | First Week of Coding!!                                           | `Weekly Check-in #2 <https://blogs.python-gsoc.org/en/nibba2018s-blog/weekly-check-in-2-3/>`__     |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 3(14-06-2020)    | ComboBox2D Progress!!                                            | `Weekly Check-in #3 <https://blogs.python-gsoc.org/en/nibba2018s-blog/weekly-check-in-3-4/>`__     |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 4(21-06-2020)    | TextBlock2D Progress!!                                           | `Weekly Check-in #4 <https://blogs.python-gsoc.org/en/nibba2018s-blog/weekly-check-in-4-4/>`__     |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 5(28-06-2020)    | May the Force be with you!!                                      | `Weekly Check-in #5 <https://blogs.python-gsoc.org/en/nibba2018s-blog/weekly-check-in-5-4/>`__     |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 6(05-07-2020)    | Translation, Reposition, Rotation.                               | `Weekly Check-in #6 <https://blogs.python-gsoc.org/en/nibba2018s-blog/weekly-check-in-6-7/>`__     |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 7(12-07-2020)    | Orientation, Sizing, Tab UI.                                     | `Weekly Check-in #7 <https://blogs.python-gsoc.org/en/nibba2018s-blog/weekly-check-in-7-4/>`__     |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 8(19-07-2020)    | ComboBox2D, TextBlock2D, ClippingOverflow.                       | `Weekly Check-in #8 <https://blogs.python-gsoc.org/en/nibba2018s-blog/weekly-check-in-8-2/>`__     |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 9(26-07-2020)    | Tab UI, TabPanel2D, Tab UI Tutorial.                             | `Weekly Check-in #9 <https://blogs.python-gsoc.org/en/nibba2018s-blog/weekly-check-in-9-4/>`__     |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 10(02-08-2020)   | Single Actor, Physics, Scrollbars.                               | `Weekly Check-in #10 <https://blogs.python-gsoc.org/en/nibba2018s-blog/weekly-check-in-10-2/>`__   |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 11(09-08-2020)   | Chain Simulation, Scrollbar Refactor,Tutorial Update.            | `Weekly Check-in #11 <https://blogs.python-gsoc.org/en/nibba2018s-blog/weekly-check-in-11-1/>`__   |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 12(16-08-2020)   | Wrecking Ball Simulation, ScrollbarsUpdate, Physics Tutorials.   | `Weekly Check-in #12 <https://blogs.python-gsoc.org/en/nibba2018s-blog/weekly-check-in-12/>`__     |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+
| Week 13(23-08-2020)   | Part of the Journey is the end unless itsOpen Source!            | `Weekly Check-in #13 <https://blogs.python-gsoc.org/en/nibba2018s-blog/weekly-check-in-13/>`__     |
+-----------------------+------------------------------------------------------------------+----------------------------------------------------------------------------------------------------+

Detailed weekly tasks and work done can be found
`here <https://blogs.python-gsoc.org/en/nibba2018s-blog/>`__.
Shader Showcase
=====================

.. post:: August 24 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

Make sure to check out Project `FURY <https://github.com/fury-gl/fury>`_

Hey Everyone!
Today marked the official end of the coding period for Google Summer of Code 2020. On this day I would like to take the opportunity to thank all my mentors and Soham who have shown immense support during this time and helped me grow not only to be a better programmer but also to be a better team member. While the GSoC period ends, I will try my best to be active and contribute to the project and help it grow.
Cheers!

What did you do this week?
--------------------------
This being the final week of GSoC , my task was to create a PR which showcases not only the shader capabilities of the project but also to create a example which integrates both the UI and shader system of project FURY. This example can help new users to get familiar with both the UI and shaders.
Apart from this i also worked on a Toon Shader.

The output for the above task is given below :


.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-13.gif


The shader demos are available `here <https://github.com/lenixlobo/fury/tree/shader-demos>`_

What is coming up next week?
----------------------------
The next week I will work on the final GSoC documentation which explains what I worked on throughout the GSoC period. In case of any changes are requested by the mentors I will also try to implement them.

Did you get stuck anywhere?
---------------------------
With the help of Soham and the mentors this week went smoothly.Orientation, Sizing, Tab UI.
============================

.. post:: July 12 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

Hello and welcome to my 7th weekly check-in. I will be sharing my progress with single actor physics simulation and ``TextBlock2D`` sizing issue which was pending for quite a while now. I will also be sharing my thoughts regarding the ``TAB UI`` component. The official repository of my sub-org, FURY can always be found `here <https://github.com/fury-gl/fury/>`_.

What did you do this week?
--------------------------
This week I was working with the orientation problem I mentioned in my previous check-in. I did have some problems regarding it, but thanks to my mentor I was able to figure it out and implement it with the help of scipy's Quaternion to Rotation Matrix `method <https://github.com/scipy/scipy/blob/v1.5.1/scipy/spatial/transform/rotation.py#L174-L1978>`_. After understanding the orientation implementation, I spent some time regarding the design of the TAB UI component. I expect the design to be something similar as follows:

.. image:: https://user-images.githubusercontent.com/29832615/87254337-906b0b80-c49f-11ea-93f3-3af0f2d8de10.png

Fortunately, we have good news this week. The ``TextBlock2D`` sizing issue that we were facing for quite a while has now been fixed. We simply have to keep track of the ``scene`` parameter in the ``_add_to_scene`` method of each UI component. This scene parameter inherits ``vtkViewport``, hence it allows us to determine the necessary details about the environment of the 3D objects.

What is coming up next week?
----------------------------
As the sizing issue regarding ``TextBlock2D`` has been fixed, I am very much eager to finish the remaining work left. ``ComboBox2D`` and ``TextBlock2D`` were left incomplete because of the issue. Therefore, my main objective next week would be to finish them first. Once I am done with it, I will move on to ``Tab UI``. If not, I will continue with the physics simulation.

Did you get stuck anywhere?
---------------------------
Apart from the orientation and sizing problem, I did not face any major issues.

``Thank you for reading, see you next week!!``
Google Summer of Code
=====================

.. post:: January 5 2020
   :author: skoudoro
   :tags: google
   :category: gsoc


FURY is participating in the `Google Summer of Code 2020 <https://summerofcode.withgoogle.com/>`_ under the umbrella of the `Python Software Foundation <https://python-gsoc.org/>`_.

FURY is a free and open source software library for scientific visualization and 3D animations. FURY contains many tools for visualizing a series of scientific data including graph and imaging data.

A list of project ideas and application info is on our `GitHub Wiki <https://github.com/fury-gl/fury/wiki/Google-Summer-of-Code-2020>`_.

If you are interested in talking to us about projects, applications join us to our `discord community <https://discord.gg/aXRZmmM>`_ or drop us a line on our `mailing list <https://mail.python.org/mailman3/lists/fury.python.org>`_.

Be part of our community and Enjoy your summer of code!

Serge K.Chain Simulation, Scrollbar Refactor, Tutorial Update.
======================================================

.. post:: August 09 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

Hello and welcome to my 11th weekly check-in. In this blog I will be discussing my progress with multiple topics related to physics and ui components. I was actively working on a couple of things, specifically Joint simulations in pyBullet and scrollbar UI component. I also took up the responsibility to complete an incomplete Pull Request which was pending for quite a while. The official repository of my sub-org can be found `here <https://github.com/fury-gl/fury/>`_.

What did you do this week?
--------------------------
The first thing that I did this week was to figure out joint simulations in pybullet. Due to lack of proper documentation I was not aware that Joints are kept stiff by default, hence I had no idea what was wrong with my simulations. Thankfully, while I was browsing pybullet forums, I found this `post <https://pybullet.org/Bullet/phpBB3/viewtopic.php?f=24&t=13035>`_ regarding rope simulations when I realized that I had to explicitly set the friction force to prevent stiffness among the Joints. Keeping this point in mind I was able to simulate the following Chain of hexagonal prisms:

.. image:: https://user-images.githubusercontent.com/29832615/89737601-b7613100-da8f-11ea-947f-a96c66caefae.gif

This week I was mainly supposed to  work on refactoring scrollbars as a standalone component. I have made some progress for now. I am able to render the scrollbars properly, with proper size, orientation and color but I am having some issues regarding its scrolling callbacks. I need to look further into it. Here's a brief glimpse:

.. image:: https://user-images.githubusercontent.com/29832615/89738159-28a2e300-da94-11ea-9167-e825f82edf98.png

This particular `PR <https://github.com/fury-gl/fury/pull/208>`_ by a fellow contributor was pending for quite a while, so I decided to look into it and complete it. The objective of the PR was to add examples for the ``CheckBox`` and ``RadioButton`` UI components, but the problem was that the objects were not rendered using FURY API in the tutorial, so I decided to complete that. It was already a well made tutorial. I only had to replace the appropriate functions with FURY's API calls.

The ``CheckBox`` tutorial:

.. image:: https://user-images.githubusercontent.com/29832615/89438967-20326b80-d767-11ea-8f47-e7711e900c9f.gif

There's still some internal issues while updating the colors of the cube which is currently being worked on by my mentors.

The ``RadioButton`` tutorial:

.. image:: https://user-images.githubusercontent.com/29832615/89438999-2e808780-d767-11ea-8b08-2a36a05294bc.gif

What is coming up next week?
----------------------------
Next week I will continue working on the scrollbar component and try to fix the issues that I am having with its callbacks. I will also try to work on the wrecking ball simulation.

Did you get stuck anywhere?
---------------------------
Apart from the scrollbar callbacks and stiff joints, I did not face any major issues.

``Thank you for reading, see you next week!!``
Raymarching continued
======================

.. post:: June 21 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

Make sure to check out Project `FURY <https://github.com/fury-gl/fury>`_

Hey ! 
Raymarching continued

What did you do this week?
--------------------------
As you read in my last blog post, while the SDF primitives were working , there was slight deformation in the render. So the main focus for this week was working on solving the deformation bug. Initially the suspect was possibly the coordinate space in which the ray marching algorithm was being computed, however after testing multiple combination of transformations the issue wasn't solved. To avoid getting stuck too long on a single bug, I decided to simultaneously work on any alternatives to the current approach. So i started working on the 2nd approach. The idea was to render multiple primitives in a single cube rather than one SDF per cube. This turned out to be highly beneficial as while implementing this , i realized what was causing the deformation .

I have added the GIFs for both the renders below. I also worked on a lambertian lighting model to create a more realistic render.



Multiple Primitives within a single cube:

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-4a.gif

Solved deformation with added lambertian Lighting: 

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-4b.gif

The code for the above render is available at the `branch <https://github.com/lenixlobo/fury/tree/SDF-Experiments>`_

What is coming up next week?
----------------------------
The next task assigned is to add support for spherical harmonics as primitives.

Did you get stuck anywhere?
---------------------------
I was stuck on the deformation issue for most of the week, but was eventually able to solve that.
Welcome to my GSoC Blog!!!
==========================

.. post:: May 30 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

Hello Everyone, this is **Soham Biswas** currently in 2nd year pursuing my Bachelor’s(B.Tech) degree in Computer Science & Engineering from Institute of Engineering & Management, Kolkata. I have been selected for GSoC' 20 at sub-org **FURY** under the umbrella organization of **Python Software Foundation**. I will be working on building sci-fi-like 2D and 3D interfaces and provide physics engine integration under project titled "Create new UI widgets & Physics Engine Integration".

What did you do during the Community Bonding Period?
----------------------------------------------------
Due to the pandemic outbreak and the country wide lockdown in India, many places including my university were closed and therefore I decided to get a head start and start with the project early. During the community bonding period, we had video conference meetings with our mentors and the project's core team. We interacted with each other and discussed the implementation details and their respective deadlines for the entire event period. We will be having such meetings every week on Wednesday in order to update ourselves about the progress of our respective tasks.

I completed the remaining Pull Requests that I had pending before the GSoC students announcement. I also reviewed other issues and pull requests to make sure everything remains up-to-date.

What is coming up next week?
----------------------------
Currently, I am focusing on building the ComboBox2D UI element. I will try to get the skeleton model, required sub-components and their respective default callbacks done by next week.

While working with my previous PR related to *Double Click* callbacks, I faced an issue where I was unable to create and implement User Events properly in VTK. Thankfully, my mentors helped me out I was able to implement double click callbacks for all three mouse buttons successfully.

``See you next week, cheers!!``
Raymarching!!
=====================

.. post:: June 14 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

Make sure to check out Project `FURY <https://github.com/fury-gl/fury>`_

Hey ! 
This week, Raymarching!

What did you do this week?
--------------------------
This was an exciting week as i got to learn and implement the ray marching algorithm in the FURY repo. In the weekly meeting, the mentors suggested adding support for SDF modelled actors as an alternative to the existing FURY actors. After going through a few implementations of ray marching in GLSL, i proceeded with the implementation in VTK. After being able to render a torus , the next logical step was to add support for multiple actors in the same window. The below render shows support for multiple SDF actors :

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-3.gif

The code for the above render is available at the `branch <https://github.com/lenixlobo/fury/tree/SDF-Experiments>`_

What is coming up next week?
----------------------------
In the above output, there is some deformation in some of the cubes, The next step is to get rid of this deformation .
Also i will be working on adding lighting within the shaders for a slightly more realistic experience.

Did you get stuck anywhere?
---------------------------
Going through and understanding theVTK documentation was quite a challenging task, however whenever i was confused the doubts were immediately cleared by the mentors
Single Actor, Physics, Scrollbars.
====================================

.. post:: August 02 2020
   :author: Soham Biswas
   :tags: google
   :category: gsoc

Hello and welcome to my 10th weekly check-in. Second evaluation ended this week and now we move on to our 3rd and final coding period. In today's check-in I will be sharing my progress with the single actor physics simulation that I facing some issues with and I will also be discussing my future plans regarding UI components. The official repository of my sub-org, FURY can always be found `here <https://github.com/fury-gl/fury/>`__.

What did you do this week?
--------------------------
This week I was able to figure out the uncontrollable spinning problem that I was facing while rendering physics simulations. Specifically the simulation where a brick wall was rendered by a single actor. The spinning problem was as follows:

.. image:: https://user-images.githubusercontent.com/29832615/88485303-87476780-cf92-11ea-850c-cc63e1376ef8.gif

Here's how the fixed simulation looks like:

.. image:: https://user-images.githubusercontent.com/29832615/89126963-946ed400-d507-11ea-93cd-aad3a9f59ab0.gif

I was facing this particular issue because I was directly syncing the orientation of the objects in pyBullet world to the objects in the Fury world. So I decided to apply the change in orientation instead and it worked. In order to achieve this I had to keep track of the bricks' orientation at each step of the simulation, sync the change and then update the tracked orientation. Thankfully, pybullet had convenient tools to achieve this. Here's a snippet on how to update individual objects rendered by a single actor:

.. code-block:: python

    def sync_brick(object_index, multibody):
        pos, orn = p.getBasePositionAndOrientation(multibody)

        rot_mat = np.reshape(
            p.getMatrixFromQuaternion(
                p.getDifferenceQuaternion(orn, brick_orns[object_index])),
            (3, 3))

        vertices[object_index * sec: object_index * sec + sec] = \
            (vertices[object_index * sec: object_index * sec + sec] -
            brick_centers[object_index])@rot_mat + pos

        brick_centers[object_index] = pos
        brick_orns[object_index] = orn

All the necessary information is updated `here <https://docs.google.com/document/d/1XJcG1TL5ZRJZDyi8V76leYZt_maxGp0kOB7OZIxKsTA/edit?usp=sharing>`_.

What is coming up next week?
----------------------------
Currently, the scrollbars are native to ``ListBox2D`` only. We are planning to separate scrollbars from ``ListBox2D`` to create a standalone UI component. This was in progress previously but was later discontinued, so I was given the responsibility to complete it. After this we plan to improve File Dialog capabilities later on.

Did you get stuck anywhere?
---------------------------
I did not face any major issues but it took me some time to understand and evaluate the existing discontinued `PR <https://github.com/fury-gl/fury/pull/222>`_ regarding scrollbar separation.

``Thank you for reading, see you next week!!``
More Shaders!!
=====================

.. post:: August 09 2020
   :author: Lenix Lobo
   :tags: google
   :category: gsoc

Make sure to check out Project `FURY <https://github.com/fury-gl/fury>`_

Hey ! 
This week, More Shaders!

What did you do this week?
--------------------------
The task assigned for this week was to explore more shader techniques which could be implemented using FURY and which would demonstrate the capability of FURY shader system. So i decided to work on implementing shading examples such as Gooch shading and reflection shader using textures.


Below are the outputs of the techniques i worked on :

.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-11a.gif


.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-11b.gif


.. image:: https://raw.githubusercontent.com/lenixlobo/fury-outputs/master/blog-week-11c.gif

The shader demos are available `here <https://github.com/lenixlobo/fury/tree/shader-demos>`_

What is coming up next week?
----------------------------
The next week will involve working on more such demos which can demonstrate the capabilities of FURY 

Did you get stuck anywhere?
---------------------------
No issues were faced this week.. _releasev0.2.0:

==================================
 Release notes v0.2.0 (2019-03-08)
==================================

Quick Overview
--------------
* Replace ``fury.window.Renderer`` by ``fury.window.Scene``
* Add stereo support
* Add GridUI object
* Increase tests coverage and code quality

Details
-------

GitHub stats for 2018/11/26 - 2019/03/08 (tag: v0.1.4)

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 5 authors contributed 186 commits.

* David Reagan
* Eleftherios Garyfallidis
* Jon Haitz Legarreta Gorroño
* Marc-Alexandre Côté
* Serge Koudoro


We closed a total of 60 issues, 20 pull requests and 40 regular issues;
this is the full list (generated with the script 
:file:`ext/github_tools.py`):

Pull Requests (20):

* :ghpull:`61`: [Fix] fetch github api from documentation
* :ghpull:`55`: Add stereo rendering support
* :ghpull:`50`: [NF] Doc version
* :ghpull:`60`: backward compatibilities
* :ghpull:`59`: [NF] Add get_info function
* :ghpull:`58`: Tests addition
* :ghpull:`57`: Add tests to utils module 
* :ghpull:`52`: [Fix] add filterwarnings
* :ghpull:`56`: Increase Codacy rank
* :ghpull:`32`: The grid - add actors in a grid and interact with each one of them independently
* :ghpull:`47`: BUG: Fix `TensorSlicerActor` actor opacity property not being enforced.
* :ghpull:`48`: ENH: Exercise the `PeakSlicerActor` `opacity` explicitly.
* :ghpull:`43`: BUG: Fix `peaks_slicer` actor properties not being enforced.
* :ghpull:`44`: BUG: Fix elementwise comparison deprecation warning.
* :ghpull:`42`: [Fix] viz_surface
* :ghpull:`40`: Re-enable transparency test, change colors
* :ghpull:`39`: removing widget module
* :ghpull:`21`: Add depth_cue and fake_tube to simulate tubes with lines
* :ghpull:`30`: Add doc generation on Travis
* :ghpull:`28`: Renaming Renderer to Scene

Issues (40):

* :ghissue:`61`: [Fix] fetch github api from documentation
* :ghissue:`55`: Add stereo rendering support
* :ghissue:`50`: [NF] Doc version
* :ghissue:`60`: backward compatibilities
* :ghissue:`59`: [NF] Add get_info function
* :ghissue:`58`: Tests addition
* :ghissue:`8`: dipy.viz.colormap crash on single fibers
* :ghissue:`57`: Add tests to utils module 
* :ghissue:`52`: [Fix] add filterwarnings
* :ghissue:`46`: Hide/Ignore numpy_vtk support warning
* :ghissue:`56`: Increase Codacy rank
* :ghissue:`32`: The grid - add actors in a grid and interact with each one of them independently
* :ghissue:`49`: Add a Codacy badge to README.rst
* :ghissue:`47`: BUG: Fix `TensorSlicerActor` actor opacity property not being enforced.
* :ghissue:`48`: ENH: Exercise the `PeakSlicerActor` `opacity` explicitly.
* :ghissue:`43`: BUG: Fix `peaks_slicer` actor properties not being enforced.
* :ghissue:`22`: Peak slicer doesn't honor linewidth parameter
* :ghissue:`37`: Fix DeprecationWarning
* :ghissue:`44`: BUG: Fix elementwise comparison deprecation warning.
* :ghissue:`45`: Change Miniconda version
* :ghissue:`42`: [Fix] viz_surface
* :ghissue:`41`: module 'fury.window' has no attribute 'Scene'
* :ghissue:`6`: VTK and Python 3 support in fvtk
* :ghissue:`40`: Re-enable transparency test, change colors
* :ghissue:`2`: Dipy visualization (fvtk) crash when saving series of images
* :ghissue:`4`: fvtk contour function ignores voxsz parameter
* :ghissue:`1`: fvtk.label won't show up if called twice 
* :ghissue:`39`: removing widget module
* :ghissue:`21`: Add depth_cue and fake_tube to simulate tubes with lines
* :ghissue:`3`: Dipy visualization with missing (?) affine parameter
* :ghissue:`5`: How to resolve python-vtk6 link issues in Ubuntu 
* :ghissue:`29`: Added surface function
* :ghissue:`30`: Add doc generation on Travis
* :ghissue:`23`: DOC: sphinx_gallery master branch is required
* :ghissue:`28`: Renaming Renderer to Scene
* :ghissue:`26`: Rename Renderer to Scene
* :ghissue:`24`: VTK dependency on installation
* :ghissue:`11`: Reorienting peak_slicer and ODF_slicer
* :ghissue:`14`: dipy test failed on mac osx sierra with ananoda python.
* :ghissue:`17`: dipy test failed on mac osx sierra with ananoda python.
.. _releasev0.7.1:

==============================
 Release notes v0.7.1
 ==============================

Quick Overview
--------------

* FURY paper added.
* Fast selection of multiple objects added.
* UI refactored.
* Tests coverage increased.
* New actor (Marker) added.
* New primitive (Triangular Prism) added.
* Demos added and updated.
* Large Documentation Update.


Details
-------

GitHub stats for 2021/03/13 - 2021/08/03 (tag: v0.7.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 15 authors contributed 211 commits.

* Amit Chaudhari
* Antriksh Misri
* Bruno Messias
* Daniel S. Katz
* Eleftherios Garyfallidis
* Gurdit Siyan
* Javier Guaje
* Jhalak Gupta
* LoopThrough-i-j
* MIHIR
* Praneeth Shetty
* Sajag Swami
* Serge Koudoro
* Hariharan Ayappane


We closed a total of 89 issues, 35 pull requests and 54 regular issues;
this is the full list (generated with the script
:file:`tools/github_stats.py`):

Pull Requests (35):

* :ghpull:`475`: Gsoc blog 2021
* :ghpull:`476`: Google Summer of Code blog posts
* :ghpull:`477`: added blog posts for GSoC'21
* :ghpull:`442`: added method to wrap overflowing text
* :ghpull:`441`: Added border support in Panel2D
* :ghpull:`466`: two more small bib changes
* :ghpull:`464`: Paper Dan's Comments
* :ghpull:`459`: extracted Button2D class from `elements` to `core`
* :ghpull:`430`: Surface actor colormap fix
* :ghpull:`456`: Updated about documentation
* :ghpull:`455`: Fixed bibtex
* :ghpull:`454`: Added missing DOIs and URLs
* :ghpull:`451`: Typo fix
* :ghpull:`447`: UI refactoring
* :ghpull:`438`: Fast selection of multiple objects in 3D using GPU acceleration
* :ghpull:`420`: added an example about graph-tool  and nested stochastic block model
* :ghpull:`422`: This allow to draw markers using shaders
* :ghpull:`444`: Remove deprecated functions
* :ghpull:`440`: added support for URL image in ImageContainer2D
* :ghpull:`356`: Render a video on an actor.
* :ghpull:`436`: [Fix] Update Azure pipeline for windows
* :ghpull:`434`: WIP: added tests for layout module
* :ghpull:`426`: Allows to define the priority of a shader_callback and obtain the vtkEventId
* :ghpull:`394`: Fixed warnings in test_utils.py
* :ghpull:`415`: update sk orcid
* :ghpull:`413`: add nanohub doi
* :ghpull:`412`: fix paper doi
* :ghpull:`386`: FURY paper for Journal of Open Source Software (JOSS)
* :ghpull:`371`: Textbox2d  special character support
* :ghpull:`408`: Updating the Missing Parenthesis
* :ghpull:`406`: Removed unused library in FURY tutorial
* :ghpull:`405`: Updating Redirecting Issues in Readme
* :ghpull:`399`: Resolve warnings #317 & and Fix Issue: #355
* :ghpull:`393`: added primitive and actor for triangular prism, added tests too
* :ghpull:`396`: #317 Fixing Warnings during test : test_actors.py

Issues (54):

* :ghissue:`407`: UI Textbox background doesn't resize according to text in it.
* :ghissue:`421`: Implementing the Resizing of Listbox
* :ghissue:`416`: Fixing the Resizing Background issue of TextBox2D UI.
* :ghissue:`475`: Gsoc blog 2021
* :ghissue:`476`: Google Summer of Code blog posts
* :ghissue:`477`: added blog posts for GSoC'21
* :ghissue:`442`: added method to wrap overflowing text
* :ghissue:`441`: Added border support in Panel2D
* :ghissue:`466`: two more small bib changes
* :ghissue:`464`: Paper Dan's Comments
* :ghissue:`445`: [WIP] Example to show how to render multiple bonds
* :ghissue:`410`: added BulletList to UI
* :ghissue:`459`: extracted Button2D class from `elements` to `core`
* :ghissue:`429`: Colormap not working as intended with surface actor
* :ghissue:`430`: Surface actor colormap fix
* :ghissue:`450`: Issue with references related to JOSS review
* :ghissue:`456`: Updated about documentation
* :ghissue:`455`: Fixed bibtex
* :ghissue:`454`: Added missing DOIs and URLs
* :ghissue:`453`: Add missing DOIs and URLs
* :ghissue:`451`: Typo fix
* :ghissue:`439`: [WIP] Space filling model
* :ghissue:`447`: UI refactoring
* :ghissue:`438`: Fast selection of multiple objects in 3D using GPU acceleration
* :ghissue:`420`: added an example about graph-tool  and nested stochastic block model
* :ghissue:`422`: This allow to draw markers using shaders
* :ghissue:`444`: Remove deprecated functions
* :ghissue:`440`: added support for URL image in ImageContainer2D
* :ghissue:`356`: Render a video on an actor.
* :ghissue:`436`: [Fix] Update Azure pipeline for windows
* :ghissue:`434`: WIP: added tests for layout module
* :ghissue:`403`: Creating test for layout module
* :ghissue:`411`: Added Layout test file
* :ghissue:`426`: Allows to define the priority of a shader_callback and obtain the vtkEventId
* :ghissue:`417`: Fixing pep issues
* :ghissue:`394`: Fixed warnings in test_utils.py
* :ghissue:`415`: update sk orcid
* :ghissue:`414`: Duplicate ORCIDs in the JOSS paper
* :ghissue:`413`: add nanohub doi
* :ghissue:`412`: fix paper doi
* :ghissue:`386`: FURY paper for Journal of Open Source Software (JOSS)
* :ghissue:`371`: Textbox2d  special character support
* :ghissue:`409`: Segmentation Fault When Running Fury Tests
* :ghissue:`408`: Updating the Missing Parenthesis
* :ghissue:`406`: Removed unused library in FURY tutorial
* :ghissue:`405`: Updating Redirecting Issues in Readme
* :ghissue:`375`: Visuals for some parametric 2D functions
* :ghissue:`317`: Track and fix warnings during tests.
* :ghissue:`355`: [Vulnerability Bug] Used blacklisted dangerous function call that can lead to RCE
* :ghissue:`399`: Resolve warnings #317 & and Fix Issue: #355
* :ghissue:`393`: added primitive and actor for triangular prism, added tests too
* :ghissue:`395`: FURY installation conflict
* :ghissue:`396`: #317 Fixing Warnings during test : test_actors.py
* :ghissue:`358`: Updated io.py
.. _releasev0.1.4:

==================================
 Release notes v0.1.4 (2018-11-26)
==================================

Quick Overview
--------------

This is a maintenance release

* Add vtk.utils.color on window package
* Update Readme, codecov, travis..... _releasev0.4.0:

=========================================
 Release notes v0.4.0 (2019-10-29)
=========================================

Quick Overview
--------------

* Enable Anti aliasing and frame rate features
* Add multiples actors (arrow, box, ...)
* Glyph extentions
* Remove Nose dependency
* Replace Appveyor by Azure pipeline for Windows
* Update Documentation, examples and tutorials
* Increase tests coverage and code quality

Details
-------
GitHub stats for 2019/08/02 - 2019/10/29 (tag: v0.3.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 4 authors contributed 169 commits.

* Eleftherios Garyfallidis
* Etienne St-Onge
* Javier Guaje
* Serge Koudoro


We closed a total of 32 issues, 15 pull requests and 17 regular issues;
this is the full list (generated with the script
:file:`tools/github_stats.py`):

Pull Requests (15):

* :ghpull:`109`: Show frame rate and enable anti-aliasing
* :ghpull:`112`: Use glyph for several actors
* :ghpull:`110`: Experimental work
* :ghpull:`107`: Increase coverage for test_actors
* :ghpull:`106`: [Fix] doc generation
* :ghpull:`104`: Made it simpler to change color in ListBox
* :ghpull:`105`: Replacing appveyor by Azure
* :ghpull:`103`: Remove nose dependency
* :ghpull:`101`: Update travis
* :ghpull:`102`: from assert to npt_assert_equal
* :ghpull:`98`: [Fix] layout with small cells
* :ghpull:`97`: Bug fix for double slider
* :ghpull:`100`:  fix snapshot when size is not square
* :ghpull:`92`: [Fix] update travis to manage pip
* :ghpull:`94`: [miniconda] move to https

Issues (17):

* :ghissue:`109`: Show frame rate and enable anti-aliasing
* :ghissue:`112`: Use glyph for several actors
* :ghissue:`66`: Directed Arrows
* :ghissue:`110`: Experimental work
* :ghissue:`107`: Increase coverage for test_actors
* :ghissue:`106`: [Fix] doc generation
* :ghissue:`104`: Made it simpler to change color in ListBox
* :ghissue:`105`: Replacing appveyor by Azure
* :ghissue:`103`: Remove nose dependency
* :ghissue:`101`: Update travis
* :ghissue:`102`: from assert to npt_assert_equal
* :ghissue:`98`: [Fix] layout with small cells
* :ghissue:`97`: Bug fix for double slider
* :ghissue:`96`: Double slider handles not at right position when window starts
* :ghissue:`100`:  fix snapshot when size is not square
* :ghissue:`92`: [Fix] update travis to manage pip
* :ghissue:`94`: [miniconda] move to https
.. _releasev0.5.1:

=========================================
 Release notes v0.5.1 (2020-04-01)
=========================================

Quick Overview
--------------

* Remove python 2 compatibility
* Added texture management
* Added multiples primitives.
* Added multiples actors (contour_from_label, billboard...)
* Huge improvement of multiple UI (RangeSlider, ...)
* Improved security (from md5 to sha256)
* Large documentation update, examples and tutorials
* Increased tests coverage and code quality

Details
-------
GitHub stats for 2019/10/29 - 2020/04/02 (tag: v0.4.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 20 authors contributed 407 commits.

* ChenCheng0630
* Devanshu Modi
* Eleftherios Garyfallidis
* Etienne St-Onge
* Filipi Nascimento Silva
* Gottipati Gautam
* Javier Guaje
* Jon Haitz Legarreta Gorroño
* Liam Donohue
* Marc-Alexandre Côté
* Marssis
* Naman Bansal
* Nasim
* Saransh Jain
* Serge Koudoro
* Shreyas Bhujbal
* Soham Biswas
* Vivek Choudhary
* ibrahimAnis
* lenixlobo


We closed a total of 153 issues, 49 pull requests and 104 regular issues;
this is the full list (generated with the script
:file:`tools/github_stats.py`):

Pull Requests (49):

* :ghpull:`227`: [Fix] update streamlines default color
* :ghpull:`210`: Added contour_from_label method
* :ghpull:`225`: update tutorial folder structure
* :ghpull:`223`: [Fix] sphere winding issue
* :ghpull:`218`: Changed options attribute from list to dict and updated respective tests
* :ghpull:`220`: bumping scipy version to 1.2.0
* :ghpull:`213`: Utils vtk
* :ghpull:`215`: Remove more than one actors at once
* :ghpull:`207`: updated fetcher
* :ghpull:`206`: [FIX] avoid in-place replacements
* :ghpull:`203`: Namanb009 windowtitlefix
* :ghpull:`204`: Vertical Layout for RangeSlider
* :ghpull:`190`: Add initial state to checkbox
* :ghpull:`201`: [FIX] icons flipping
* :ghpull:`181`: Vertical Layout for LineDoubleSlider2D
* :ghpull:`198`: Utils test and winding order algorithm
* :ghpull:`192`: Tetrahedron, Icosahedron primitives
* :ghpull:`189`: Added dynamic text positioning
* :ghpull:`194`: [FIX] Update superquadrics test
* :ghpull:`182`: [Doc] Reshape the documentation
* :ghpull:`177`: [Fix] Flipping during save
* :ghpull:`191`: DOC: Fix `actor.line` parameter type and add `optional` keyword
* :ghpull:`173`: Fixing Text Overflow of ListBox2D
* :ghpull:`167`: Animated Network Visualization Example
* :ghpull:`165`: Vertical Layout for LineSlider2D
* :ghpull:`154`: Added Shader tutorial
* :ghpull:`153`: Sep viz ui
* :ghpull:`132`: Add Billboard actor
* :ghpull:`164`: Documentation
* :ghpull:`163`: Spelling error
* :ghpull:`157`: Corrected Disk2D comments
* :ghpull:`148`: Replace md5 by sha 256
* :ghpull:`145`: DOC: Fix `io:load_image` and `io:save_image` docstrings
* :ghpull:`144`: STYLE: Change examples `README` file extension to reStructuredText
* :ghpull:`143`: STYLE: Improve the requirements' files' style.
* :ghpull:`139`: [Fix] some docstring for doc generation
* :ghpull:`140`: [DOC] Add demo for showing an network
* :ghpull:`136`: Started new tutorial about using normals to make spiky spheres
* :ghpull:`134`: Add event parameter on add_window_callback method in ShowManager class.
* :ghpull:`129`: update loading and saving IO for polydata
* :ghpull:`131`: Add Superquadric primitives and actors
* :ghpull:`130`: Adding Sphere primitives
* :ghpull:`128`: Update Deprecated function
* :ghpull:`126`: Add basic primitives
* :ghpull:`125`: Add Deprecated decorator
* :ghpull:`124`: Texture utilities and actors
* :ghpull:`118`: Remove python2 compatibility
* :ghpull:`120`: Replace pickle with JSON for "events_counts" dict serialization
* :ghpull:`115`: Release 0.4.0 preparation

Issues (104):

* :ghissue:`150`: Re-compute Bounds in Slicer
* :ghissue:`227`: [Fix] update streamlines default color
* :ghissue:`135`: Backward compatibilities problem with streamtube
* :ghissue:`77`: contour_from_label
* :ghissue:`210`: Added contour_from_label method
* :ghissue:`225`: update tutorial folder structure
* :ghissue:`223`: [Fix] sphere winding issue
* :ghissue:`137`: Issues with provided spheres
* :ghissue:`152`: Improve checkbox options cases
* :ghissue:`218`: Changed options attribute from list to dict and updated respective tests
* :ghissue:`76`: Improve Checkbox options access
* :ghissue:`219`: Issue occur when I Start testing the project
* :ghissue:`220`: bumping scipy version to 1.2.0
* :ghissue:`217`: Transformed options attribute from list to dict and updated respective tests
* :ghissue:`213`: Utils vtk
* :ghissue:`179`: Utility functions are needed for getting numpy arrays from actors
* :ghissue:`212`: Namanb009 issue 133 fix
* :ghissue:`214`: Namanb009 Remove mulitple actors
* :ghissue:`215`: Remove more than one actors at once
* :ghissue:`211`: Namanb009 hexadecimal color support
* :ghissue:`187`: New utility functions are added in utils.py and tests are added in te…
* :ghissue:`209`: Namanb009 viz_ui.py does not show render window when run
* :ghissue:`207`: updated fetcher
* :ghissue:`206`: [FIX] avoid in-place replacements
* :ghissue:`203`: Namanb009 windowtitlefix
* :ghissue:`202`: Window Title name does not change
* :ghissue:`204`: Vertical Layout for RangeSlider
* :ghissue:`190`: Add initial state to checkbox
* :ghissue:`75`: Improve Checkbox initialisation
* :ghissue:`201`: [FIX] icons flipping
* :ghissue:`199`: Loading of Inverted icons using read_viz_icons
* :ghissue:`181`: Vertical Layout for LineDoubleSlider2D
* :ghissue:`175`: LineDoubleSlider2D vertical layout
* :ghissue:`198`: Utils test and winding order algorithm
* :ghissue:`192`: Tetrahedron, Icosahedron primitives
* :ghissue:`189`: Added dynamic text positioning
* :ghissue:`176`: Allowing to change text position on Sliders
* :ghissue:`185`: NF: winding order in utils
* :ghissue:`170`: NF: adding primitive stars, 3D stars, rhombi.
* :ghissue:`195`: Added dynamic text position on sliders
* :ghissue:`194`: [FIX] Update superquadrics test
* :ghissue:`171`: bug-in-image 0.1
* :ghissue:`182`: [Doc] Reshape the documentation
* :ghissue:`156`: Test Case File Updated
* :ghissue:`155`: There are libraries we have to install not mentioned in the requirement.txt file to run the test case.
* :ghissue:`122`: Documentation not being rendered correctly
* :ghissue:`177`: [Fix] Flipping during save
* :ghissue:`160`: Saved Images are vertically Inverted
* :ghissue:`193`: Merge pull request #2 from fury-gl/master
* :ghissue:`191`: DOC: Fix `actor.line` parameter type and add `optional` keyword
* :ghissue:`178`: changed text position
* :ghissue:`188`: Added dynamic text positioning
* :ghissue:`173`: Fixing Text Overflow of ListBox2D
* :ghissue:`15`: viz.ui.ListBoxItem2D text overflow
* :ghissue:`166`: Build Native File Dialogs
* :ghissue:`180`: Native File Dialog Text Overflow Issue
* :ghissue:`186`: add name
* :ghissue:`184`: Added winding order algorithm to utils
* :ghissue:`183`: Added star2D and 3D, rhombicuboctahedron to tests_primitive
* :ghissue:`54`: generating directed arrows
* :ghissue:`174`: List box text overflow
* :ghissue:`167`: Animated Network Visualization Example
* :ghissue:`165`: Vertical Layout for LineSlider2D
* :ghissue:`108`: Slider vertical layout
* :ghissue:`172`: window.show() is giving Attribute error.
* :ghissue:`154`: Added Shader tutorial
* :ghissue:`151`: Prim shapes
* :ghissue:`162`: Winding order 2
* :ghissue:`168`: Prim test
* :ghissue:`158`: nose is missing
* :ghissue:`71`: viz_ui.py example needs expansion
* :ghissue:`153`: Sep viz ui
* :ghissue:`132`: Add Billboard actor
* :ghissue:`164`: Documentation
* :ghissue:`163`: Spelling error
* :ghissue:`161`: Merge pull request #1 from fury-gl/master
* :ghissue:`157`: Corrected Disk2D comments
* :ghissue:`121`: Replace md5 by sha2 or sha3 for security issue
* :ghissue:`148`: Replace md5 by sha 256
* :ghissue:`147`: update md5 to sha256
* :ghissue:`146`: Shapes
* :ghissue:`145`: DOC: Fix `io:load_image` and `io:save_image` docstrings
* :ghissue:`144`: STYLE: Change examples `README` file extension to reStructuredText
* :ghissue:`142`: STYLE: Change examples `README` file extension to markdown
* :ghissue:`143`: STYLE: Improve the requirements' files' style.
* :ghissue:`139`: [Fix] some docstring for doc generation
* :ghissue:`140`: [DOC] Add demo for showing an network
* :ghissue:`136`: Started new tutorial about using normals to make spiky spheres
* :ghissue:`134`: Add event parameter on add_window_callback method in ShowManager class.
* :ghissue:`81`: Add superquadric function in actor.py
* :ghissue:`129`: update loading and saving IO for polydata
* :ghissue:`131`: Add Superquadric primitives and actors
* :ghissue:`130`: Adding Sphere primitives
* :ghissue:`128`: Update Deprecated function
* :ghissue:`126`: Add basic primitives
* :ghissue:`125`: Add Deprecated decorator
* :ghissue:`124`: Texture utilities and actors
* :ghissue:`99`: [WIP] Adding util to get Numpy 3D array of RGBA values
* :ghissue:`118`: Remove python2 compatibility
* :ghissue:`117`: Remove compatibility with python 2
* :ghissue:`123`: WIP: Texture support
* :ghissue:`119`: Improve data Serialization
* :ghissue:`120`: Replace pickle with JSON for "events_counts" dict serialization
* :ghissue:`115`: Release 0.4.0 preparation
.. _releasev0.6.0:

=========================================
 Release notes v0.6.0 (2020-07-20)
=========================================

Quick Overview
--------------

* Added new features: Picking and DoubleClick callback.
* Added Signed Distance Field actor.
* Added a new UI ComboBox.
* Added multiples primitives (Rhombocuboctahedron, ...).
* Huge improvement of multiple UIs and actors.
* Fixed Compatibility with VTK9.
* Large documentation update, examples and tutorials (5 new).
* Added a blog system.
* Increased tests coverage and code quality.

Details
-------

GitHub stats for 2020/04/09 - 2020/07/20 (tag: v0.5.1)

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 9 authors contributed 266 commits.

* Eleftherios Garyfallidis
* Liam Donohue
* Marc-Alexandre Côté
* Melina Raglin
* Naman Bansal
* Serge Koudoro
* Soham Biswas
* Tushar
* Lenix Lobo


We closed a total of 60 issues, 25 pull requests and 35 regular issues;
this is the full list (generated with the script
:file:`tools/github_stats.py`):

Pull Requests (25):

* :ghpull:`241`: Adding  a Blog system to fury
* :ghpull:`268`: Clip Text Overflow
* :ghpull:`267`: SDF tutorial
* :ghpull:`250`: SDF based FURY actors
* :ghpull:`246`: ComboBox UI tutorial
* :ghpull:`265`: Adding new Earth coordinates tutorial
* :ghpull:`240`: ComboBox2D UI element
* :ghpull:`169`: Use primitive for box, cube, square, rectangle actors
* :ghpull:`262`: Solar System Tutorial
* :ghpull:`248`: Setting Bounding Box for TextBlock2D
* :ghpull:`263`: [NF] Deprecated parameters
* :ghpull:`127`: [CI] add Python 3.8
* :ghpull:`255`: Let VTK delete timers on exit.
* :ghpull:`233`: Picking API and picking tutorial
* :ghpull:`261`: Adding Earth Animation Tutorial
* :ghpull:`249`: Octagonal Prism and Frustum Square Pyramid
* :ghpull:`258`: Updating test of order_transparency for compatibility with vtk 9
* :ghpull:`259`: Updated Step 3 in README.rst
* :ghpull:`231`: Adding Double Click Callback
* :ghpull:`256`: Install ssl certificate for azure pipeline windows
* :ghpull:`245`: [FIX]  Compatibility with VTK 9
* :ghpull:`244`: Added new texture tutorial
* :ghpull:`235`: Function to use Hexadecimal color code in Colormap
* :ghpull:`238`: Added Rhombocuboctahedron, 2D and 3D star to primitive
* :ghpull:`237`: update copyright years

Issues (35):

* :ghissue:`241`: Adding  a Blog system to fury
* :ghissue:`268`: Clip Text Overflow
* :ghissue:`264`: Re-implementation of Text Overflow in ListBox2D
* :ghissue:`267`: SDF tutorial
* :ghissue:`247`: PR idea: create SDF alternatives to FURY primitive actors
* :ghissue:`250`: SDF based FURY actors
* :ghissue:`246`: ComboBox UI tutorial
* :ghissue:`265`: Adding new Earth coordinates tutorial
* :ghissue:`240`: ComboBox2D UI element
* :ghissue:`169`: Use primitive for box, cube, square, rectangle actors
* :ghissue:`138`: Box, cone etc. to work similarly to superquadric
* :ghissue:`262`: Solar System Tutorial
* :ghissue:`248`: Setting Bounding Box for TextBlock2D
* :ghissue:`263`: [NF] Deprecated parameters
* :ghissue:`127`: [CI] add Python 3.8
* :ghissue:`51`: Improvements from VTK 8.2.0?
* :ghissue:`255`: Let VTK delete timers on exit.
* :ghissue:`253`: Programs with timers hang on exit. [VTK9] [Linux]
* :ghissue:`233`: Picking API and picking tutorial
* :ghissue:`261`: Adding Earth Animation Tutorial
* :ghissue:`249`: Octagonal Prism and Frustum Square Pyramid
* :ghissue:`258`: Updating test of order_transparency for compatibility with vtk 9
* :ghissue:`254`: unexpected order transparent behavior [VTK9] [Ubuntu 18.04]
* :ghissue:`259`: Updated Step 3 in README.rst
* :ghissue:`251`: Developer installation instructions should describe -e option
* :ghissue:`226`: Adding DoubleClick event
* :ghissue:`231`: Adding Double Click Callback
* :ghissue:`256`: Install ssl certificate for azure pipeline windows
* :ghissue:`245`: [FIX]  Compatibility with VTK 9
* :ghissue:`244`: Added new texture tutorial
* :ghissue:`235`: Function to use Hexadecimal color code in Colormap
* :ghissue:`238`: Added Rhombocuboctahedron, 2D and 3D star to primitive
* :ghissue:`197`: Added Rhombocuboctahedron, 2D and 3D star to primitive
* :ghissue:`237`: update copyright years
* :ghissue:`216`: Utiltiy function to use Hexadecimal color code
.. _releasev0.1.1:

==================================
 Release notes v0.1.1 (2018-10-29)
==================================

Quick Overview
--------------

This is a maintenance release

* Fix error on python 2.7
* Travis integration
* Documentation integration.. _releasev0.7.0:

===================================
 Release notes v0.7.0 (2021/03/13)
===================================

Quick Overview
--------------

* New SDF actors added.
* Materials module added.
* ODF slicer actor performance improved.
* New primitive (Cylinder) added.
* Compatibility with VTK 9 added.
* Five new demos added.
* Large Documentation Update.
* Migration from Travis to Github Action.


Details
-------

GitHub stats for 2020/08/20 - 2021/03/13 (tag: v0.6.1)

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 14 authors contributed 195 commits.

* Eleftherios Garyfallidis
* Serge Koudoro
* Charles Poirier
* Javier Guaje
* Soham Biswas
* Sajag Swami
* Lenix Lobo
* Pietro Astolfi
* Sanjay Marreddi
* Tushar
* ganimtron-10
* haran2001
* Aju100
* Aman Soni


We closed a total of 98 issues, 37 pull requests and 61 regular issues;
this is the full list (generated with the script
:file:`tools/github_stats.py`):

Pull Requests (37):

* :ghpull:`388`: added simulation for brownian motion
* :ghpull:`389`: ENH: peaks_slicer option for asymmetric peaks visualization
* :ghpull:`370`: Materials module including Physically Based Rendering (PBR)
* :ghpull:`385`: fixed the example for superquadric function
* :ghpull:`387`: [fix]   Propagate update_actor
* :ghpull:`382`: Added an sdf for rendering a Capsule actor
* :ghpull:`383`: Minor documentation fix
* :ghpull:`376`: Added animations for some electromagnetic phenomena
* :ghpull:`374`: ENH: Refactor actor.odf_slicer for increased performances
* :ghpull:`373`: Updated actor.py
* :ghpull:`368`: Solving The Minor Documentation Error
* :ghpull:`343`: Adding physics engine integration docs
* :ghpull:`353`: fix: Minor docs changes
* :ghpull:`346`: Fix the sdf bug by checking the arguments passed
* :ghpull:`351`: Opacity bug fix for point and sphere actors
* :ghpull:`350`: modelsuzanne to suzanne
* :ghpull:`348`: Added center forwarding in billboard shaders.
* :ghpull:`341`: Add Option to generate the documentation without examples
* :ghpull:`342`: From Travis to Github Actions
* :ghpull:`339`: Update Readme informations
* :ghpull:`340`: Pass OAuth token through header
* :ghpull:`337`: Add support for clipping side in clip_overflow_text
* :ghpull:`336`: Update UI tutorials.
* :ghpull:`334`: Added Domino-Simulation-file for Review
* :ghpull:`332`: Fixing UI warnings
* :ghpull:`328`: Added cylinder primitive
* :ghpull:`329`: [FIX] Force LUT to be RGB
* :ghpull:`286`: GSoC blogs for Third Evaluation.
* :ghpull:`319`: fixed discord icon bug in documentation
* :ghpull:`311`: Remove python35 from Travis
* :ghpull:`307`: Fixed translating and scaling issues on billboard and SDF actors
* :ghpull:`304`: Blogs for the final review
* :ghpull:`306`: merged basic UI and advanced UI tutorials into one
* :ghpull:`302`: moved physics tutorials to examples under the heading 'Integrate physics using pybullet'
* :ghpull:`303`: FIX vtp reader
* :ghpull:`300`: BF: Out should be varying and alpha is not passed to shader
* :ghpull:`295`: Update fetcher

Issues (61):

* :ghissue:`388`: added simulation for brownian motion
* :ghissue:`389`: ENH: peaks_slicer option for asymmetric peaks visualization
* :ghissue:`370`: Materials module including Physically Based Rendering (PBR)
* :ghissue:`385`: fixed the example for superquadric function
* :ghissue:`387`: [fix]   Propagate update_actor
* :ghissue:`382`: Added an sdf for rendering a Capsule actor
* :ghissue:`383`: Minor documentation fix
* :ghissue:`376`: Added animations for some electromagnetic phenomena
* :ghissue:`374`: ENH: Refactor actor.odf_slicer for increased performances
* :ghissue:`364`: New Animated Network Demo/Example
* :ghissue:`379`: Merge pull request #2 from fury-gl/master
* :ghissue:`361`: Closes #352
* :ghissue:`373`: Updated actor.py
* :ghissue:`372`: Ellipsoid primitive needs to be added in the comment section of sdf actor.
* :ghissue:`369`: Added Special Character Support
* :ghissue:`363`: Minor error in documentation of create_colormap function
* :ghissue:`368`: Solving The Minor Documentation Error
* :ghissue:`366`: added special character support for TextBox2D
* :ghissue:`357`: Patches: vulnerable code that can lead to RCE
* :ghissue:`359`: unwanted objects rendering randomly
* :ghissue:`343`: Adding physics engine integration docs
* :ghissue:`312`: Adding Physics Integration Docs to FURY's Website
* :ghissue:`353`: fix: Minor docs changes
* :ghissue:`346`: Fix the sdf bug by checking the arguments passed
* :ghissue:`310`: Rendering bug in SDF actor when not all primitives are defined
* :ghissue:`351`: Opacity bug fix for point and sphere actors
* :ghissue:`335`: _opacity argument for point doesn't seem to work
* :ghissue:`345`: Fixes the opacity bug for sphere and point actors (unit tests are included)
* :ghissue:`350`: modelsuzanne to suzanne
* :ghissue:`348`: Added center forwarding in billboard shaders.
* :ghissue:`341`: Add Option to generate the documentation without examples
* :ghissue:`342`: From Travis to Github Actions
* :ghissue:`338`: From travis (pricing model changed) to github Actions ?
* :ghissue:`339`: Update Readme informations
* :ghissue:`340`: Pass OAuth token through header
* :ghissue:`315`: Deprecation notice for authentication via URL query parameters
* :ghissue:`337`: Add support for clipping side in clip_overflow_text
* :ghissue:`308`: Clipping overflowing text from the left.
* :ghissue:`336`: Update UI tutorials.
* :ghissue:`334`: Added Domino-Simulation-file for Review
* :ghissue:`309`: Domino Physics Simulation
* :ghissue:`333`: Unable to set up the project locally for python 32bit system
* :ghissue:`332`: Fixing UI warnings
* :ghissue:`239`: Superquadric Slicer
* :ghissue:`328`: Added cylinder primitive
* :ghissue:`318`: Cylinder primitive generation
* :ghissue:`329`: [FIX] Force LUT to be RGB
* :ghissue:`286`: GSoC blogs for Third Evaluation.
* :ghissue:`319`: fixed discord icon bug in documentation
* :ghissue:`313`: Discord icon should appear in doc too
* :ghissue:`311`: Remove python35 from Travis
* :ghissue:`307`: Fixed translating and scaling issues on billboard and SDF actors
* :ghissue:`274`: SDF rendering bug for low scale values
* :ghissue:`304`: Blogs for the final review
* :ghissue:`306`: merged basic UI and advanced UI tutorials into one
* :ghissue:`297`: Update Demos/tutorial
* :ghissue:`302`: moved physics tutorials to examples under the heading 'Integrate physics using pybullet'
* :ghissue:`298`: Wrecking Ball Simulation
* :ghissue:`303`: FIX vtp reader
* :ghissue:`300`: BF: Out should be varying and alpha is not passed to shader
* :ghissue:`295`: Update fetcher
.. _releasev0.6.1:

===================================
 Release notes v0.6.1 (2020-08-20)
===================================

Quick Overview
--------------

* Added Shaders Manager.
* Standardized colors across API.
* Added a new UI Tab.
* Added Physics Engine Tutorial.
* Large documentation update, examples and tutorials (4 new).
* Increased tests coverage and code quality.

Details
-------

GitHub stats for 2020/07/21 - 2020/08/20 (tag: v0.6.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 8 authors contributed 164 commits.

* Eleftherios Garyfallidis
* Javier Guaje
* Lenix Lobo
* Melina Raglin
* Nasim Anousheh
* Serge Koudoro
* Soham Biswas
* Vivek Choudhary


We closed a total of 42 issues, 15 pull requests and 27 regular issues;
this is the full list (generated with the script
:file:`tools/github_stats.py`):

Pull Requests (15):

* :ghpull:`288`: Shader manager
* :ghpull:`292`: Disable VTK Warnings on release version
* :ghpull:`289`: Rename parameter scale to scales
* :ghpull:`287`: Add Physics simulation examples.
* :ghpull:`284`: Colliding particles in the box
* :ghpull:`275`: Tab UI tutorial
* :ghpull:`208`: Added tutorial for RadioButton and CheckBox UI
* :ghpull:`281`: Radio checkbox tutorial using Fury API
* :ghpull:`283`: Standardize colors across API
* :ghpull:`282`: Standardize Colors array name
* :ghpull:`252`: Tab ui
* :ghpull:`279`: Decreasing the size of the sun in solarsystem tutorial
* :ghpull:`273`: Python GSoC Weekly blogs
* :ghpull:`276`: Update Deprecated test
* :ghpull:`272`: Python GSoC Blogs upto 19th July 2020

Issues (27):

* :ghissue:`260`: Changes to shader API in VTK 9
* :ghissue:`116`: Update Shader system
* :ghissue:`288`: Shader manager
* :ghissue:`292`: Disable VTK Warnings on release version
* :ghissue:`270`: Disable VTK Warnings on release version
* :ghissue:`289`: Rename parameter scale to scales
* :ghissue:`236`: Pybullet examples
* :ghissue:`287`: Add Physics simulation examples.
* :ghissue:`205`: Create a tutorial for checkbox/radiobutton UI
* :ghissue:`284`: Colliding particles in the box
* :ghissue:`275`: Tab UI tutorial
* :ghissue:`208`: Added tutorial for RadioButton and CheckBox UI
* :ghissue:`281`: Radio checkbox tutorial using Fury API
* :ghissue:`283`: Standardize colors across API
* :ghissue:`269`: Fixed bug in tutorial with accessing colors for latest release
* :ghissue:`242`: Standardize colors across API
* :ghissue:`243`: Single Color in Primitive Does Not Work
* :ghissue:`271`: Some issues with actor.box after release 0.6.0
* :ghissue:`282`: Standardize Colors array name
* :ghissue:`280`: Unable to extract colors from a Box
* :ghissue:`252`: Tab ui
* :ghissue:`279`: Decreasing the size of the sun in solarsystem tutorial
* :ghissue:`278`: Changing Size of Sun in viz_solar_system.py Tutorial
* :ghissue:`273`: Python GSoC Weekly blogs
* :ghissue:`277`: Sun
* :ghissue:`276`: Update Deprecated test
* :ghissue:`272`: Python GSoC Blogs upto 19th July 2020
.. _releasev0.3.0:

=========================================
 Release notes v0.3.0 (2019-08-02)
=========================================

Quick Overview
--------------
* Add cone actor and update odf actor
* Add Appveyor CI and update MacOS CI
* Update Documentation, examples and tutorials
* Increase tests coverage and code quality

Details
-------

GitHub stats for 2019/03/08 - 2019/08/02 (tag: v0.2.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 7 authors contributed 164 commits.

* Ariel Rokem
* Eleftherios Garyfallidis
* Guillaume Favelier
* Kevin Sitek
* Prashil
* Scott Trinkle
* Serge Koudoro


We closed a total of 39 issues, 15 pull requests and 24 regular issues;
this is the full list (generated with the script
:file:`tools/github_stats.py`):

Pull Requests (15):

* :ghpull:`91`: [Doc] add installation instruction
* :ghpull:`53`: [Fix] resolution limit
* :ghpull:`89`: [Fix] cmap when not have matplotlib
* :ghpull:`88`: Updates and corrections in slicer and ui
* :ghpull:`69`: Add middle button callback
* :ghpull:`86`: Fixes a typo in the viz_ui tutorial.
* :ghpull:`83`: Remove issue template
* :ghpull:`87`: TEST: updated order transparency issue for vtk 8.2.0
* :ghpull:`80`: Add cones as glyphs
* :ghpull:`73`: Add appveyor
* :ghpull:`72`: Update OSX bots on Travis
* :ghpull:`68`: Allow Str for Grid Caption
* :ghpull:`67`: [Fix] Update doc management
* :ghpull:`62`: Directional color odfs
* :ghpull:`31`: new surface function

Issues (24):

* :ghissue:`91`: [Doc] add installation instruction
* :ghissue:`36`: Tests Documentation
* :ghissue:`53`: [Fix] resolution limit
* :ghissue:`13`: window.record() resolution limit
* :ghissue:`89`: [Fix] cmap when not have matplotlib
* :ghissue:`90`: [Fix] dtype problem for x64 machine
* :ghissue:`88`: Updates and corrections in slicer and ui
* :ghissue:`69`: Add middle button callback
* :ghissue:`86`: Fixes a typo in the viz_ui tutorial.
* :ghissue:`84`: Test_order_transparent failed with VTK 8.2.0
* :ghissue:`83`: Remove issue template
* :ghissue:`87`: TEST: updated order transparency issue for vtk 8.2.0
* :ghissue:`85`: Save from active window?
* :ghissue:`79`: add link to fury example gallery in sphinx-gallery readme
* :ghissue:`80`: Add cones as glyphs
* :ghissue:`73`: Add appveyor
* :ghissue:`72`: Update OSX bots on Travis
* :ghissue:`18`: Improve unit tests
* :ghissue:`63`: Improve doc generation
* :ghissue:`68`: Allow Str for Grid Caption
* :ghissue:`67`: [Fix] Update doc management
* :ghissue:`62`: Directional color odfs
* :ghissue:`65`: Directed Arrows
* :ghissue:`31`: new surface function
.. _releasev0.1.3:

=============================================
 Release notes v0.1.2 and v0.1.3 (2018-10-31)
=============================================

Quick Overview
--------------

This is a maintenance release

* Update setup.py
* Remove dependence on requirements.txt.. _releasev0.1.0:

==================================
 Release notes v0.1.0 (2018-09-21)
==================================

Quick Overview
--------------

This initial release is a split from DIPY. It contains:

* from ``dipy.viz.actors`` to ``fury.actors``
* from ``dipy.viz.window`` to ``fury.window``
* from ``dipy.viz.ui`` to ``fury.ui``
* from ``dipy.viz.widget`` to ``fury.widget``
* from ``dipy.viz.interactor`` to ``fury.interactor``
* from ``dipy.viz.colormap`` to ``fury.colormap``
* from ``dipy.viz.utils`` to ``fury.utils``Demos
=====

Below is a gallery of Demos. A bunch of apps powered by FURY.Integrate Physics using pybullet
--------------------------------

These demos show how to use connect FURY with a physics engine library.Introductory
------------

These tutorials show:

- How to combine a timer with an actor
- How to slice data with the slicer actor
- How to use the normals of your data.User Interface Elements
-----------------------

These tutorials show how to create user interfaces elements.Shaders
-------

These tutorials show:

- How to use shaders in FURY actors.
- How to create new user shaders and internal conventions.