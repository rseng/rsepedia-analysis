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
```